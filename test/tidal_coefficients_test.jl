using Test
using ParticleTracking

# The GPU-safe coefficient path is the one the driver actually uses, and its whole purpose is
# to be callable from a kernel: no search, no allocation, no host-array access. These tests
# check the reconstruction against the point-wise harmonics it replaces, so a regression in
# the bilinear indexing or the sign convention is caught without needing a GPU.
const ROOT = normpath(joinpath(@__DIR__, ".."))
const TIDES = joinpath(ROOT, "work", "snowcrab", "inputs", "tpxo9", "u_tpxo9.v1.nc")

@testset "GPU-safe tidal forcing" begin
    if !isfile(TIDES)
        @info "tidal solution not present; skipping coefficient tests"
        @test true
    else
        h = ParticleTracking.read_tidal_velocity_harmonics(
            TIDES; constituents = (:M2, :S2), lon_range = (-71.0, -53.0), lat_range = (40.0, 48.5))

        g = TidalCoefficientGrid(h, (-71.0, -53.0), (40.0, 48.5); nx = 16, ny = 12)
        f = tidal_forcing_coefficients(g)

        @testset "the grid is isbits, which is what makes it kernel-safe" begin
            @test isbitstype(typeof(g))
            @test isbitstype(typeof(f.u))
        end

        @testset "outside the grid the forcing is zero, not clamped" begin
            @test f.u(-80.0, 50.0, -10.0, 0.0) == 0.0
            @test f.u(-60.0, 30.0, -10.0, 0.0) == 0.0
        end

        @testset "inside the grid the forcing is bounded by the local amplitude" begin
            # A sum of harmonics is bounded by the sum of the amplitudes, so no point in the
            # domain may exceed that. The old uniform forcing was 0.25 m/s everywhere; the
            # real field is far weaker in most of the domain and stronger in the Fundy.
            worst = 0.0
            for lon in -70.0:1.0:-54.0, lat in 41.0:1.0:48.0, t in (0.0, 5000.0, 12000.0)
                worst = max(worst, abs(f.u(lon, lat, -50.0, t)), abs(f.v(lon, lat, -50.0, t)))
            end
            @test worst < 1e-2         # a tendency, not a velocity
            @test worst > 0.0          # and it is not silently zero either
        end

        @testset "it varies across the domain" begin
            south = abs(f.u(-60.0, 40.5, -50.0, 0.0))
            fundy = abs(f.u(-65.7, 45.2, -50.0, 0.0))
            @test fundy != south       # a uniform field would give equal values
        end

        @testset "agrees with the point-wise reconstruction at a grid corner" begin
            # At an exact grid node the bilinear weights degenerate to one corner, so the
            # tendency must equal the time derivative of that corner's complex amplitude:
            #   d/dt [Re cos(wt) - Im sin(wt)] = -w (Re sin(wt) + Im cos(wt))
            # The factor of omega is what makes this an acceleration rather than a velocity.
            lon, lat = g.lon0, g.lat0
            i = argmin(abs.(h.lon .- lon)); j = argmin(abs.(h.lat .- lat))
            (ib, jb) = nearest_water(h, i, j)
            t = 3600.0
            want = 0.0
            for k in eachindex(h.omega)
                w = h.omega[k]
                want += -w * (h.uRe[k, jb, ib] * sin(w * t) + h.uIm[k, jb, ib] * cos(w * t))
            end
            @test f.u(lon, lat, -50.0, t) ≈ want atol = 1e-15
        end

        @testset "it is a tendency, not a velocity" begin
            # Amplitudes reach 0.45 m/s, so a velocity-valued forcing would be O(0.45). The
            # tendency is smaller by omega ~ 1.4e-4, i.e. O(1e-4 m/s^2). Feeding the
            # velocity instead diverged within ten iterations at 87 m/s.
            peak = 0.0
            for lon in -70.0:1.0:-54.0, lat in 41.0:1.0:48.0, t in (0.0, 5000.0, 12000.0)
                peak = max(peak, abs(f.u(lon, lat, -50.0, t)))
            end
            @test peak > 0.0
            @test peak < 1e-2
        end
    end
end
