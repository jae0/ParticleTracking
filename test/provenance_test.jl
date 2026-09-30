using Test
using JSON3
using ParticleTracking

@testset "Data provenance registry" begin
    @testset "every declared input is either keyless or says what it needs" begin
        # This used to assert `d.keyless` for every entry. That made the field
        # untestable-by-construction: a flag that must always hold can only ever record
        # success, so the one input that genuinely needs an account -- the lateral boundary,
        # which is credentialed precisely because it is the only source carrying currents --
        # had to be marked `keyless = true` to pass, or the suite sat permanently red.
        # Marking it false is the honest record, and the property worth asserting is that a
        # false flag is never used to hide a requirement.
        for (k, d) in DATA_SOURCES
            @test !isempty(d.primary)
            @test !isempty(d.cache)
            @test d.kind in (:gridded, :station, :analytic)
            if !d.keyless
                # Whatever the flag says, the record a reader actually sees must name the
                # constraint, so nobody discovers the requirement from a 401 mid-run. The
                # full record is searched rather than `notes` alone, because that is where a
                # careful writer puts the detail -- but it must be findable somewhere.
                record = join(vcat([d.title, d.primary, d.license, d.notes], d.alternatives), " ")
                @test occursin(r"(?i)credential|account|\.cdsapirc|401", record) ||
                      error("$(k) is marked keyless = false but its record never says what " *
                            "credential or account it requires.")
            end
        end
        # The lateral boundary is the known non-keyless input, and it must stay declared as
        # such: flipping it to true to make a suite green is exactly the overclaim this
        # registry exists to prevent.
        @test data_source(:lateral_boundary).keyless == false
    end

    @testset "the seven physical inputs are all declared" begin
        # A new input added to the model without a registry entry would make a run
        # unreproducible, so the set is asserted explicitly rather than by counting.
        for k in (:bathymetry, :surface_winds, :hydrography, :oxygen,
                  :atmosphere, :lateral_boundary, :tides)
            @test haskey(DATA_SOURCES, k)
        end
    end

    @testset "tides names the verified global source and rejects the regional ones" begin
        d = data_source(:tides)
        @test occursin("zenodo", d.primary)
        @test occursin("CC-BY-4.0", d.license)
        # The correct source must be marked as measured, and the two regional archives that
        # superficially look equivalent must be marked as unusable for this study area.
        joined = join(d.alternatives, " ")
        @test occursin("MEASURED AND CORRECT", joined)
        @test occursin("lat_u spans -90..90", joined)
        @test occursin("DO NOT USE", joined)
        @test occursin("-34.97..14.97", joined)
    end

    @testset "inputs are containerised under output_dir, not a shared top-level dir" begin
        # A shared inputs/ directory let one scenario silently reuse another's bathymetry.
        for (k, d) in DATA_SOURCES
            @test startswith(d.cache, "<output_dir>/inputs/")
            @test !startswith(d.cache, "inputs/")
        end
        @test normpath(input_dir("work/snowcrab")) == joinpath("work", "snowcrab", "inputs")
    end

    @testset "unknown key errors with the valid list" begin
        err = try
            data_source(:nope)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("bathymetry", err.msg)
    end

    @testset "atmospheric forcing declares its actual limitation" begin
        # The atmosphere entry must keep saying the forcing is a single representative point
        # rather than a resolved field, so the registry cannot quietly start claiming spatial
        # structure the fetcher does not provide.
        notes = data_source(:atmosphere).notes
        @test occursin("do NOT vary in space", notes)
        @test occursin("representative point", notes)
    end

    @testset "provenance records what was read, not what was asked for" begin
        dir = mktempdir(; cleanup = true)
        p = data_provenance(dir; config_path = "configs/x.toml",
                            extra = Dict{String, Any}("seed" => 42))
        @test isfile(p)
        rec = JSON3.read(read(p, String))
        @test rec["config"] == "configs/x.toml"
        @test rec["seed"] == 42
        @test haskey(rec["inputs"], "tides")
        @test haskey(rec["inputs"], "bathymetry")
        for (_, v) in rec["inputs"]
            @test haskey(v, "source")
            @test haskey(v, "cache")
            @test haskey(v, "sha256")
        end
    end

    @testset "digest of a real file is stable" begin
        dir = mktempdir(; cleanup = true)
        f = joinpath(dir, "x.bin")
        write(f, UInt8[1, 2, 3])
        d1 = file_digest(f)
        @test length(d1) == 64          # SHA-256 hex
        @test d1 == file_digest(f)
        @test file_digest(joinpath(dir, "absent")) == ""
    end
end
