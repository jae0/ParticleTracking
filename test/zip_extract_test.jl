using Test
using ParticleTracking

# The ZIP reader in src/data/manifest.jl has to cope with real archives, which are full of
# fields that disagree with each other. Every quirk asserted here was found by comparing
# against a known-good extraction, not by reading the specification.
const ROOT = normpath(joinpath(@__DIR__, ".."))
const ZIP = joinpath(ROOT, "inputs", "tpxo9_atlas_v2.zip")
const GOOD = joinpath(ROOT, "inputs", "tpxo9_atlas_v2", "TPXO9_atlas_v2",
                      "h_m2_tpxo9_atlas_30_v2.nc")

@testset "ZIP member extraction" begin
    if !isfile(ZIP) || !isfile(GOOD)
        @info "local TPXO archive not present; skipping extraction tests"
        @test true
    else
        member = "TPXO9_atlas_v2/h_m2_tpxo9_atlas_30_v2.nc"
        dest = joinpath(mktempdir(), "out.nc")
        ParticleTracking._extract_zip_member(ZIP, member, dest)

        @test isfile(dest)
        @test filesize(dest) == filesize(GOOD)

        a = read(dest); b = read(GOOD)
        @test length(a) == length(b)
        @test a == b                      # byte-for-byte, not just the same size

        @testset "a member that is not in the archive is a clear error" begin
            err = try
                ParticleTracking._extract_zip_member(ZIP, "nope/nothing.nc", dest)
                nothing
            catch e
                e
            end
            @test err isa ErrorException
        end
    end
end
