using Test
using ParticleTracking

const DIR = normpath(joinpath(@__DIR__, "..", "inputs", "tpxo9_atlas_v2"))
const isdir_here = isdir(DIR)

@testset "TPXO elevation harmonics (real archive)" begin
    if !isdir_here
        @warn "TPXO archive not extracted at $DIR; skipping real-data tests"
        @test true
    else
        LONR = (-71.0, -53.0)
        LATR = (40.0, 48.5)

        # The openly mirrored archive in the registry is a *tropical-band subset*: it spans
        # lat -35..15 only. The Scotian Shelf at 40-48.5N is outside it, so the reader must
        # refuse rather than quietly return an empty or wrong box. This test asserts that
        # refusal, which is the behaviour that keeps a bad archive from being used.
        err = try
            read_tidal_elevation_harmonics(DIR; constituents = (:M2,), lon_range = LONR,
                                           lat_range = LATR)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test err !== nothing && occursin("outside the coverage", err.msg)
        @test err !== nothing && occursin("lat", err.msg)

        @testset "a covered box is read correctly (equatorial band)" begin
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2, :S2),
                                               lon_range = (-30.0, 10.0),
                                               lat_range = (-10.0, 10.0))
            @test h.constituents == [:M2, :S2]
            @test issorted(h.lon) && issorted(h.lat)
            @test all(-30.0 .<= h.lon .<= 10.0)
            @test all(-10.0 .<= h.lat .<= 10.0)
            @test size(h.hRe) == (2, length(h.lat), length(h.lon))
            @test h.omega[1] ≈ get_tidal_frequency(:M2)
            @test h.omega[2] ≈ get_tidal_frequency(:S2)
        end

        @testset "longitude wrapped from 0-360 into -180..180" begin
            # A negative-longitude request must match, which it cannot on a raw 0-360 axis.
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2,),
                                               lon_range = (-30.0, -10.0),
                                               lat_range = (-5.0, 5.0))
            @test !isempty(h.lon)
            @test all(-30.0 .<= h.lon .<= -10.0)
        end

        @testset "physical sanity: equatorial amplitude is O(1 m)" begin
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2, :S2),
                                               lon_range = (-30.0, 10.0),
                                               lat_range = (-10.0, 10.0))
            i, j = length(h.lon) ÷ 2, length(h.lat) ÷ 2
            amp = tidal_elevation_amplitude(h, i, j)
            # A file stored in millimetres but read as metres gives ~447 here; a file read in
            # metres but genuinely in metres would give hundreds of metres. Anything outside
            # this band means the unit scaling or the sign convention is wrong.
            @test amp > 0.05
            @test amp < 10.0
        end

        @testset "amplitudes are stored in millimetres, read as metres" begin
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2,),
                                               lon_range = (-30.0, 10.0),
                                               lat_range = (-10.0, 10.0))
            i, j = length(h.lon) ÷ 2, length(h.lat) ÷ 2
            a = tidal_elevation_amplitude(h, i, j)
            @test a < 5.0
            @test a > 0.01
        end

        @testset "reconstruction is a bounded, oscillating signal" begin
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2, :S2),
                                               lon_range = (-30.0, 10.0),
                                               lat_range = (-10.0, 10.0))
            i, j = length(h.lon) ÷ 2, length(h.lat) ÷ 2
            samples = [tidal_elevation(h, i, j, k * 600.0) for k in 0:200]
            amp = tidal_elevation_amplitude(h, i, j)
            @test maximum(abs, samples) <= amp + 1e-9
            @test maximum(samples) - minimum(samples) > 0.5 * amp
        end

        @testset "phase is not discarded" begin
            h = read_tidal_elevation_harmonics(DIR; constituents = (:M2,),
                                               lon_range = (-30.0, 10.0),
                                               lat_range = (-10.0, 10.0))
            i, j = length(h.lon) ÷ 2, length(h.lat) ÷ 2
            quarter = get_tidal_frequency(:M2)^(-1) * pi / 2
            @test abs(tidal_elevation(h, i, j, 0.0) - tidal_elevation(h, i, j, quarter)) > 1e-6
        end

        @testset "a constituent with no file is an error" begin
            err = try
                read_tidal_elevation_harmonics(DIR; constituents = (:NOPE,))
                nothing
            catch e
                e
            end
            @test err isa ErrorException
        end

        @testset "missing directory points at the keyless source" begin
            err = try
                read_tidal_elevation_harmonics("inputs/does_not_exist")
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin("zenodo", err.msg)
        end
    end
end
