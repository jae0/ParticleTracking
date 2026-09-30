using NCDatasets
using Test
using ParticleTracking
const C = ParticleTracking

const TMP = mktempdir(; cleanup = true)

# ---------------------------------------------------------------------------
# Build a synthetic TPXO9-Atlas-shaped file: known amplitudes and phases, so
# the reader can be checked against values computed by hand.
# ---------------------------------------------------------------------------
const LON = collect(-71.0:0.5:-53.0)
const LAT = collect(48.5:-0.5:40.0)          # atlas files use descending latitude
const AMPS = Dict(:M2 => 0.25, :S2 => 0.11)  # u amplitude by constituent
const PHASE = Dict(:M2 => 0.0, :S2 => 0.4)   # u phase (radians), S2 offset to catch sign errors

function write_fake_tpxo(path)
    NCDatasets.NCDataset(path, "c") do ds
        defDim(ds, "con", 2)
        defDim(ds, "lon_u", length(LON))
        defDim(ds, "lat_u", length(LAT))
        defDim(ds, "lon_v", length(LON))
        defDim(ds, "lat_v", length(LAT))

        v = defVar(ds, "cnc", String, ("con",))
        v[:] = ["M2", "S2"]
        v = defVar(ds, "lon_u", Float64, ("lon_u",)); v[:] = LON
        v = defVar(ds, "lat_u", Float64, ("lat_u",)); v[:] = LAT
        v = defVar(ds, "lon_v", Float64, ("lon_v",)); v[:] = LON
        v = defVar(ds, "lat_v", Float64, ("lat_v",)); v[:] = LAT
        v = defVar(ds, "omega", Float64, ("con",))
        v[:] = [C.get_tidal_frequency(:M2), C.get_tidal_frequency(:S2)]

        for (name, amps, phases) in (("uRe", AMPS, PHASE), ("uIm", AMPS, PHASE),
                                      ("vRe", AMPS, PHASE), ("vIm", AMPS, PHASE))
            v = defVar(ds, name, Float64, ("con", "lat_u", "lon_u"))
            for (k, sym) in enumerate((:M2, :S2)), j in eachindex(LAT), i in eachindex(LON)
                # vary with position so a transposed read would be caught
                v[k, j, i] = amps[sym] * cos(phases[sym] + 0.01 * j) +
                             0.1 * amps[sym] * sin(phases[sym])
            end
        end
    end
end

const FILE = joinpath(TMP, "tpxo9_atlas.nc")
write_fake_tpxo(FILE)

@testset "TPXO harmonic reader" begin
    @testset "missing file fails loudly, with no silent fallback" begin
        err = try
            C.read_tidal_harmonics(joinpath(TMP, "absent.nc"))
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("deliberately does not fall back", err.msg)
    end

    @testset "reads the requested constituents and box" begin
        h = C.read_tidal_harmonics(FILE; constituents = (:M2, :S2),
                                   lon_range = (-70.0, -54.0), lat_range = (41.0, 47.0))
        @test h.constituents == [:M2, :S2]
        @test issorted(h.lon) && issorted(h.lat)
        @test all(-70.0 .<= h.lon .<= -54.0)
        @test all(41.0 .<= h.lat .<= 47.0)
        @test size(h.uRe) == (2, length(h.lat), length(h.lon))
        @test size(h.vIm) == size(h.uRe)
        @test h.omega[1] ≈ C.get_tidal_frequency(:M2)
        @test h.omega[2] ≈ C.get_tidal_frequency(:S2)
    end

    @testset "omega is taken from the file, not hardcoded" begin
        h = C.read_tidal_harmonics(FILE; constituents = (:S2,))
        @test length(h.omega) == 1
        @test h.omega[1] ≈ C.get_tidal_frequency(:S2)
    end

    @testset "absent constituent is an error, not a silent zero" begin
        err = try
            C.read_tidal_harmonics(FILE; constituents = (:K1,))
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("K1", err.msg)
    end

    @testset "reconstruction matches the analytic harmonic sum" begin
        h = C.read_tidal_harmonics(FILE; constituents = (:M2, :S2))
        i, j = 3, 4
        t = 12345.0
        u_ref = v_ref = 0.0
        for (k, sym) in enumerate((:M2, :S2))
            ω = C.get_tidal_frequency(sym)
            a = AMPS[sym] * cos(PHASE[sym] + 0.01 * (length(LAT) - j + 1)) +
                0.1 * AMPS[sym] * sin(PHASE[sym])
            # the file stored descending latitude; the reader must have flipped it
            u_ref += a * cos(ω * t) + a * sin(ω * t)
            v_ref += a * cos(ω * t) + a * sin(ω * t)
        end
        u, v = C.tidal_velocity(h, i, j, t)
        @test u ≈ u_ref atol = 1e-12
        @test v ≈ v_ref atol = 1e-12
    end

    @testset "phase matters: dropping it would change the answer" begin
        h = C.read_tidal_harmonics(FILE; constituents = (:M2, :S2))
        t = 98765.0
        us = [C.tidal_velocity(h, 2, 2, t) for _ in 1:1]
        @test us[1][1] != 0.0
        # a quarter period later the M2 part must have changed phase appreciably
        later = C.tidal_velocity(h, 2, 2, t + C.get_tidal_frequency(:M2)^(-1) * pi / 2)
        @test abs(later[1] - us[1][1]) > 1e-3
    end
end
