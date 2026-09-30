using Test
using TOML
using SHA
using ParticleTracking
const C = ParticleTracking

const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "inputs are containerised under output_dir" begin
    for f in ("default.toml", "snowcrab.toml", "opendata.toml")
        cfg = TOML.parsefile(joinpath(ROOT, "configs", f))
        o = C.configuration_to_options(cfg)
        @testset "$f" begin
            @test normpath(o.input_dir) == normpath(joinpath(o.output_dir, "inputs"))
            @test !startswith(o.input_dir, "inputs")
            @test normpath(o.input_dir) ==
                  normpath(joinpath("work", "snowcrab", "inputs")) ||
                  normpath(o.input_dir) == normpath(joinpath("outputs", "default", "inputs")) ||
                  normpath(o.input_dir) == normpath(joinpath("work", "opendata", "inputs"))
        end
    end

    @testset "the specific case the requirement names" begin
        cfg = TOML.parsefile(joinpath(ROOT, "configs", "snowcrab.toml"))
        o = C.configuration_to_options(cfg)
        @test endswith(replace(o.input_dir, '\\' => '/'), "work/snowcrab/inputs")
    end

    @testset "provenance resolves the placeholder against output_dir" begin
        dir = mktempdir(; cleanup = true)
        rec0 = read(C.data_provenance(dir), String)
        @test occursin("input_dir", rec0)
        @test !occursin("<output_dir>", rec0)   # no unresolved placeholders survive
        @test occursin("\"present\":false", rec0)   # nothing cached yet

        # A file placed at a *declared* cache path must be detected and fingerprinted.
        ind = joinpath(dir, "inputs")
        mkpath(ind)
        write(joinpath(ind, "bathymetry_active.nc"), UInt8[1, 2, 3, 4])
        rec1 = read(C.data_provenance(dir), String)
        @test occursin("\"present\":true", rec1)
        # 4 bytes of [1,2,3,4] has a known SHA-256, so the fingerprint is checked, not just
        # the presence flag.
        expected = bytes2hex(SHA.sha256(UInt8[1, 2, 3, 4]))
        @test occursin(expected, rec1)
    end
end
