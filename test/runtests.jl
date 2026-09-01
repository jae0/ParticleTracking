"""
    test/runtests.jl

Automated unit test suite for the ParticleTracking package covering:
1. Synthetic and real data processing & drag laws
2. Grid and immersed boundary construction
3. Hydrodynamic model configuration & stratification
4. Climate change scenarios, PLD models, and thermal mortality
5. DVM swimming velocities, degree-day molting, and settlement filters
6. Eulerian/Lagrangian particle tracking and tidal current superposition
7. Visualization generators
"""

using Test
using Random

# Import necessary symbols from external packages
  # Required external imports for tests
using Oceananigans
using Oceananigans.Grids: LatitudeLongitudeGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Architectures
using DataFrames: DataFrame, nrow
using DataAPI
using DBInterface  # ← ADD THIS
using NCDatasets
using JLD2
 

# Load the ParticleTracking module directly
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."), io = devnull)
using ParticleTracking


@testset "ParticleTracking.jl Test Suite" begin

    @testset "1. Drag Law & Open Data Utilities" begin
        # Test Large & Pond / Wu wind stress conversion
        tau_x, tau_y = wind_speed_to_kinematic_stress(10.0, 0.0)
        @test tau_x > 0.0
        @test tau_y == 0.0

        # High wind regime (>11 m/s)
        tau_x_high, _ = wind_speed_to_kinematic_stress(20.0, 0.0)
        @test tau_x_high > tau_x

        # 2D bilinear regridding test
        src_lon = range(-65.0, -60.0, length = 10)
        src_lat = range(42.0, 46.0, length = 10)
        src_data = [x + y for x in src_lon, y in src_lat]

        tgt_lon = range(-64.5, -60.5, length = 5)
        tgt_lat = range(42.5, 45.5, length = 5)
        regridded = regrid_2d_field(src_lon, src_lat, src_data, tgt_lon, tgt_lat)
        @test size(regridded) == (5, 5)
        @test !any(isnan, regridded)
    end

    @testset "2. Synthetic Data & Inspection" begin
        synth_bathy = generate_synthetic_bathymetry(
            "inputs/test_bathy.nc",
            lon_range = (-65.0, -60.0),
            lat_range = (43.0, 46.0),
            n_lon = 20,
            n_lat = 20
        )
        @test isfile(synth_bathy)

        info = inspect_netcdf(synth_bathy, verbose = false)
        @test haskey(info[:dimensions], "lon")
        @test haskey(info[:dimensions], "lat")
    end

    @testset "3. Grid & Immersed Topography" begin
        base_grid = build_shelf_grid(
            lon_range = (-65.0, -60.0),
            lat_range = (43.0, 46.0),
            z_range = (-500.0, 0.0),
            grid_size = (20, 20, 5)
        )
        @test base_grid isa LatitudeLongitudeGrid

        immersed = build_immersed_grid_from_real_data(base_grid, "inputs/test_bathy.nc")
        @test immersed isa ImmersedBoundaryGrid
    end

    @testset "4. Hydrodynamic Model & Climate Scenarios" begin
        base_grid = build_shelf_grid(
            lon_range = (-65.0, -60.0),
            lat_range = (43.0, 46.0),
            z_range = (-500.0, 0.0),
            grid_size = (10, 10, 5)
        )
        immersed = build_immersed_grid_from_real_data(base_grid, "inputs/test_bathy.nc")
        model = build_hydrodynamic_model(immersed, coriolis_latitude = 44.0)
        @test haskey(model.tracers, :T)
        @test haskey(model.tracers, :S)

        # Apply climate change scenario
        deltas = apply_climate_scenario!(model, scenario = :ssp245, year = 2050)
        @test deltas.ΔT_surface > 0.0
        @test deltas.ΔS_surface < 0.0

        # Test PLD and mortality
        pld_cold = temperature_dependent_pld(2.0)
        pld_warm = temperature_dependent_pld(8.0)
        @test pld_cold > pld_warm # Warmer water accelerates development

        mort_base = larval_thermal_mortality_rate(6.0)
        mort_stressed = larval_thermal_mortality_rate(14.0)
        @test mort_stressed > mort_base # Elevated mortality past 10°C
    end

    @testset "5. Tidal Dynamics & Simpson-Hunter Parameter" begin
        # Test astronomical constituent frequencies
        omega_m2 = get_tidal_frequency(:M2)
        @test omega_m2 ≈ 1.405189e-4

        # Test tidal body forcing generation
        tf = build_tidal_body_forcing(constituents = [:M2], u_amplitudes = Dict(:M2 => 0.3))
        @test haskey(tf, :u)
        @test haskey(tf, :v)
        @test tf.u(0.0, 0.0, 0.0, 0.0) ≈ 0.3 * omega_m2

        # Test harmonic velocity vector
        u_t, v_t = tidal_velocity_vector(0.0, constituents = [:M2], u_amplitudes = Dict(:M2 => 0.25))
        @test u_t ≈ 0.25
        @test v_t ≈ 0.0

        # Test Simpson-Hunter tidal mixing front index
        chi_mixed = simpson_hunter_parameter(30.0, 1.2)   # shallow bank / strong tide
        chi_strat = simpson_hunter_parameter(150.0, 0.2)  # deep shelf / weak tide
        @test chi_mixed < 1.5 # well-mixed
        @test chi_strat > 2.0 # stratified
    end

    @testset "6. DVM, Tidal Superposition & Settlement Filters" begin
        # DVM swimming velocity: day should swim down (negative), night should swim up (positive)
        w_day = diel_vertical_migration_velocity(-10.0, 43200.0; stage = :zoea1) # midday at shallow depth
        w_night = diel_vertical_migration_velocity(-40.0, 0.0; stage = :zoea1)   # midnight at deep depth
        @test w_day < 0.0   # descending
        @test w_night > 0.0 # ascending

        # Tidal current superposition
        u_tide, v_tide = superpose_tidal_velocity(0.1, 0.0, 0.0; u_amp = 0.2, v_amp = 0.1)
        @test u_tide ≈ 0.3

        # Degree-day molting progression
        @test update_larval_stage(:zoea1, 50.0)  == :zoea1          # below zoea1→zoea2 (65 DD)
        @test update_larval_stage(:zoea1, 70.0)  == :zoea2          # above 65 DD
        @test update_larval_stage(:zoea2, 140.0) == :megalopa        # above 130 DD
        @test update_larval_stage(:megalopa, 210.0) == :instar1_settled  # above 200 DD

        # Benthic nursery settlement suitability
        suit_good = evaluate_settlement_suitability(-120.0, 3.5)
        @test suit_good.suitable == true

        suit_abyssal = evaluate_settlement_suitability(-800.0, 2.0)
        @test suit_abyssal.suitable == false # too deep

        suit_warm = evaluate_settlement_suitability(-100.0, 8.5)
        @test suit_warm.suitable == false # too warm
    end

    @testset "7. Lagrangian Particle Tracking & Marine Land Exclusion" begin
        # Test 2D bathymetry interpolator
        bathy_interp = get_bathymetry_interpolator("inputs/test_bathy.nc")
        @test bathy_interp(-62.5, 44.5) isa Float64
        @test bathy_interp(-62.5, 44.5) < 0.0 # Marine elevation

        # Test marine water detection
        @test is_marine_water(-62.5, 44.5, bathymetry = "inputs/test_bathy.nc") == true
        @test is_marine_water(-62.5, 44.5, bathymetry = (lon, lat) -> 10.0) == false # Land

        # Test discrete marine cell extraction
        m_cells = extract_marine_cells("inputs/test_bathy.nc", min_seabed_depth = 50.0)
        @test length(m_cells.lons) > 0
        @test all(m_cells.depths .<= -50.0)

        # Test direct deterministic marine coordinate sampling (0% land probability)
        rng = MersenneTwister(42)
        s_lons, s_lats, s_z = sample_marine_coordinates(
            20,
            "inputs/test_bathy.nc",
            min_seabed_depth = 80.0,
            rng = rng
        )
        @test length(s_lons) == 20
        @test all(s_z .<= -80.0)

        # Initialize particles enforcing minimum seabed depth (100m) over bathymetry
        larvae = initialize_larval_particles(
            15,
            lon_range = (-64.5, -60.5),
            lat_range = (43.5, 45.5),
            depth_range = (-50.0, -15.0),
            min_seabed_depth = 100.0,
            bathymetry = "inputs/test_bathy.nc",
            rng = rng
        )
        @test length(larvae.lon) == 15

        # Verify all particles are placed in marine water >= 100m deep
        for p in 1:15
            z_bed = bathy_interp(larvae.lon[p], larvae.lat[p])
            @test z_bed <= -100.0 # Must be at least 100m deep
            @test larvae.depth[p] >= z_bed # Must be above seabed
            @test larvae.depth[p] <= -1.0  # Must be below surface
        end

        # Test rejection error when requesting depth impossible in land domain
        land_fn(lon, lat) = 50.0 # Positive elevation = land
        @test_throws ErrorException initialize_larval_particles(
            5,
            min_seabed_depth = 100.0,
            bathymetry = land_fn,
            max_attempts = 100,
            rng = rng
        )

        # Test shoreline reflection during transport step
        # Create a flow field strongly directing particles eastward toward a land boundary
        flow_toward_land(lon, lat, z, t) = (1.5, 0.0, 0.0)
        coastal_bathy(lon, lat) = lon > -62.0 ? 50.0 : -150.0 # Land wall at lon > -62.0

        # Initialize particles strictly on the marine side of the coastline (lon in -64.5 to -62.5)
        larvae_coast = initialize_larval_particles(
            10,
            lon_range = (-64.5, -62.5),
            lat_range = (43.5, 45.5),
            depth_range = (-50.0, -15.0),
            bathymetry = coastal_bathy,
            rng = rng
        )

        trajs = track_larval_cohort(
            larvae_coast,
            velocity_fn = flow_toward_land,
            temperature_fn = (x, y, z, t) -> 4.0,
            bathymetry_fn = coastal_bathy,
            total_duration = 3600.0, # 1 hour
            dt = 300.0,
            enable_tides = true,
            enable_molting = true,
            rng = rng
        )

        @test size(trajs.lons) == (10, 13)
        # Verify particles NEVER penetrated land despite eastward flow
        for p in 1:10, s in 1:13
            @test coastal_bathy(trajs.lons[p, s], trajs.lats[p, s]) < 0.0
            @test trajs.lons[p, s] <= -62.0
        end
    end

    @testset "8. Empirical Movement & Dispersion Analysis" begin
        rng = MersenneTwister(42)
        larvae = initialize_larval_particles(20, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.08, 0.04, 0.0),
            total_duration = 3600.0 * 6,
            dt = 600.0,
            rng = rng
        )

        emp = estimate_empirical_movement(
            trajs,
            lon_bins = range(-68.0, -57.0, length = 15),
            lat_bins = range(42.0, 47.0, length = 15)
        )

        @test length(emp.lon_centers) == 14
        @test length(emp.lat_centers) == 14
        @test any(!isnan, emp.u_mean)
        @test any(!isnan, emp.v_mean)
        @test any(!isnan, emp.diffusivity)
    end

    @testset "9. Gridded Recruitment & Thermal Exposure Metrics" begin
        rng = MersenneTwister(42)
        larvae = initialize_larval_particles(20, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.05, 0.02, 0.0),
            temperature_fn = (x, y, z, t) -> 5.0,
            bathymetry_fn = (x, y) -> -150.0,
            total_duration = 3600.0 * 12,
            dt = 600.0,
            enable_molting = true,
            rng = rng
        )

        # ✓ FIX: Infer bins from actual trajectory bounds
        lon_min, lon_max = extrema(trajs.lons)
        lat_min, lat_max = extrema(trajs.lats)
        # Expand slightly to ensure all particles fall within bins
        lon_bins = range(lon_min - 0.5, lon_max + 0.5, length = 30)
        lat_bins = range(lat_min - 0.5, lat_max + 0.5, length = 30)

        rec = compute_gridded_recruitment_metrics(trajs; lon_bins = lon_bins, lat_bins = lat_bins)
        @test sum(rec.release_density) == 20
        @test sum(rec.settlement_density) == 20

        therm = compute_gridded_thermal_metrics(trajs; lon_bins = lon_bins, lat_bins = lat_bins)
        @test any(!isnan, therm.mean_degree_days)
        @test any(!isnan, therm.mean_exposure_temperature)
    end

    @testset "10. Empirical Connectivity & NetCDF / JLD2 Export" begin
        rng = MersenneTwister(42)
        larvae = initialize_larval_particles(15, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.05, 0.02, 0.0),
            total_duration = 3600.0 * 6,
            dt = 600.0,
            rng = rng
        )

        conn = compute_empirical_connectivity(trajs)
        @test size(conn.matrix, 1) == length(conn.strata_names)
        @test all(sum(conn.matrix, dims = 2) .≈ 1.0) # Row normalized probabilities

        # Test multi-layer NetCDF export
        nc_file = "outputs/test_larval_dispersal.nc"
        export_larval_dispersal_netcdf(nc_file, trajectories = trajs)
        @test isfile(nc_file)

        # Inspect exported NetCDF variables
        ds = NCDataset(nc_file)
        @test haskey(ds, "emp_u")
        @test haskey(ds, "emp_diffusivity")
        @test haskey(ds, "settlement_density")
        @test haskey(ds, "mean_degree_days")
        @test haskey(ds, "connectivity_matrix")
        @test haskey(ds, "particle_lon")
        close(ds)

        # Test JLD2 export
        jld_file = "outputs/test_larval_dispersal.jld2"
        export_larval_dispersal_jld2(jld_file, trajectories = trajs)
        @test isfile(jld_file)
    end

    @testset "11. Visualizations" begin
        rng = MersenneTwister(123)
        larvae = initialize_larval_particles(10, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.05, 0.02, 0.0),
            total_duration = 1800.0,
            dt = 300.0,
            rng = rng
        )

        fig_track = plot_particle_trajectories(trajs, output_path = "outputs/test_tracks.png")
        @test isfile("outputs/test_tracks.png")

        fig_dvm = plot_vertical_migration_profiles(trajs, output_path = "outputs/test_dvm.png")
        @test isfile("outputs/test_dvm.png")

        fig_density = plot_larval_dispersal_density(trajs, output_path = "outputs/test_density.png")
        @test isfile("outputs/test_density.png")

        emp = estimate_empirical_movement(trajs)
        fig_emp = plot_empirical_movement_field(emp, output_path = "outputs/test_emp_mov.png")
        @test isfile("outputs/test_emp_mov.png")

        conn = compute_empirical_connectivity(trajs)
        fig_conn = plot_connectivity_matrix(conn, output_path = "outputs/test_conn.png")
        @test isfile("outputs/test_conn.png")

        therm = compute_gridded_thermal_metrics(trajs)
        fig_therm = plot_thermal_exposure_map(therm, output_path = "outputs/test_therm.png")
        @test isfile("outputs/test_therm.png")

        rec = compute_gridded_recruitment_metrics(trajs)
        fig_rec = plot_recruitment_summary(rec, output_path = "outputs/test_rec.png")
        @test isfile("outputs/test_rec.png")

        # Test hydrodynamic model CairoMakie figures
        fig_h_adv = plot_hydrodynamic_advection(nothing, output_path = "outputs/test_hydro_adv.png")
        @test isfile("outputs/test_hydro_adv.png")

        fig_h_trc = plot_hydrodynamic_tracers(nothing, output_path = "outputs/test_hydro_trc.png")
        @test isfile("outputs/test_hydro_trc.png")

        # Test hydrodynamic dataset extraction
        hydro_ds = extract_hydrodynamic_dataset(nothing)
        @test length(hydro_ds.lons) > 0
        @test length(hydro_ds.lats) > 0
        @test size(hydro_ds.u, 1) == length(hydro_ds.lons)
        @test size(hydro_ds.temperature, 2) == length(hydro_ds.lats)
        @test size(hydro_ds.salinity, 1) == length(hydro_ds.lons)
    end

    @testset "12. Architecture & Device Resolution" begin
        # CPU resolution
        @test resolve_architecture(:cpu) isa Oceananigans.Architectures.CPU
        @test resolve_architecture("cpu") isa Oceananigans.Architectures.CPU
        @test resolve_architecture(false) isa Oceananigans.Architectures.CPU
        @test resolve_architecture(Oceananigans.Architectures.CPU()) isa Oceananigans.Architectures.CPU

        # Fallback resolution on non-CUDA systems
        @test resolve_architecture(:gpu, fallback_to_cpu = true) isa Oceananigans.Architectures.AbstractArchitecture
        @test resolve_architecture("cuda", fallback_to_cpu = true) isa Oceananigans.Architectures.AbstractArchitecture

        # Error without fallback on systems without functional CUDA
        if !isdefined(Main, :CUDA) || !Main.CUDA.functional()
            @test_throws ErrorException resolve_architecture(:gpu, fallback_to_cpu = false)
        end
    end

    @testset "13. Interactive Map Dashboard Generation" begin
        rng = MersenneTwister(42)
        larvae = initialize_larval_particles(12, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.04, 0.01, 0.0),
            temperature_fn = (x, y, z, t) -> 4.5,
            bathymetry_fn = (x, y) -> -120.0,
            total_duration = 3600.0 * 8,
            dt = 600.0,
            rng = rng
        )

        # Test CFA polygon loading and point-in-polygon algorithms
        cfa_polys = load_cfa_polygons("inputs")
        @test length(cfa_polys) >= 3
        @test any(p -> p.code == :cfa4x, cfa_polys)
        @test any(p -> p.code == :cfanorth, cfa_polys)
        @test any(p -> p.code == :cfasouth, cfa_polys)

        # Test point in polygon for simple square and CFA polygons
        sq_lons = [0.0, 1.0, 1.0, 0.0, 0.0]
        sq_lats = [0.0, 0.0, 1.0, 1.0, 0.0]
        @test point_in_polygon(0.5, 0.5, sq_lons, sq_lats) == true
        @test point_in_polygon(1.5, 0.5, sq_lons, sq_lats) == false
        @test point_in_polygon(-0.1, 0.5, sq_lons, sq_lats) == false

        # Test CFA 4X point (e.g. Halifax / SW NS waters: -65.5, 43.5)
        cfa4x_poly = first(filter(p -> p.code == :cfa4x, cfa_polys))
        @test length(cfa4x_poly.coordinates) >= 10
        @test point_in_polygon(-65.5, 43.5, cfa4x_poly.lons, cfa4x_poly.lats) == true

        # Test Leaflet HTML export with true CFA polygons and hydrodynamic layers
        html_file = "outputs/test_interactive_tracks.html"
        res_file = export_interactive_tracks_html(html_file, trajectories = trajs, strata_definitions = cfa_polys)
        @test isfile(html_file)
        @test res_file == html_file

        # Verify content contains required Leaflet structure, polygons, and data payloads
        content = read(html_file, String)
        @test occursin("leaflet", lowercase(content))
        @test occursin("const PARTICLES =", content)
        @test occursin("const STRATA =", content)
        @test occursin("\"is_polygon\": true", content)
        @test occursin("L.polygon(s.polygon", content)
        @test occursin("playback-dock", content)
        @test occursin("telemetry-hud", content)
        @test occursin("layer-hydro-advection", content)
        @test occursin("layer-hydro-temp", content)
        @test occursin("layer-hydro-sal", content)
        @test occursin("layer-hydro-speed", content)
        @test occursin("layer-hydro-eta", content)
        @test occursin("layer-hydro-w", content)
        @test occursin("layer-hydro-bathy", content)
        @test occursin("legend-dock", content)
        @test occursin("hydro-opacity-slider", content)
    end

    @testset "14. DuckDB Storage, Multi-Scenario Querying & Ensemble Model Averaging" begin
        test_db_path = "outputs/test_particle_tracking.duckdb"
        if isfile(test_db_path)
            rm(test_db_path, force = true)
        end

        # 1. Initialize DuckDB storage and verify schema
        db = open_duckdb_storage(test_db_path)
        @test isfile(test_db_path)

        # 2. Generate and save synthetic runs across 3 climate scenarios
        rng = MersenneTwister(42)
        cfa_defs = [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]

        scenarios_to_test = [
            (name = "historical", year = 2015, speed = 0.04, temp = 3.5),
            (name = "ssp245",     year = 2050, speed = 0.06, temp = 5.2),
            (name = "ssp585",     year = 2050, speed = 0.08, temp = 7.1)
        ]

        for sc in scenarios_to_test
            larvae = initialize_larval_particles(10, rng = rng)
            trajs = track_larval_cohort(
                larvae,
                velocity_fn = (x, y, z, t) -> (sc.speed, 0.01, 0.0),
                temperature_fn = (x, y, z, t) -> sc.temp,
                bathymetry_fn = (x, y) -> -140.0,
                total_duration = 3600.0 * 6,
                dt = 600.0,
                rng = rng
            )
            conn = compute_empirical_connectivity(trajs, strata_definitions = cfa_defs)
            emp = estimate_empirical_movement(trajs)

            run_id = "test_run_$(sc.name)_$(sc.year)"
            opts_mock = (
                scenario = sc.name,
                projection_year = sc.year,
                n_particles = 10,
                track_duration = 21600.0,
                track_dt = 600.0,
                enable_tides = true,
                enable_dvm = true,
                enable_molting = true,
                diffusivity_h = 10.0,
                diffusivity_v = 1e-4,
                min_seabed_depth = 100.0,
                seed = 42
            )

            saved_id = save_simulation_run!(
                db,
                run_id,
                opts_mock;
                trajectories = trajs,
                connectivity = conn,
                gridded_dispersal = emp,
                notes = "Test run for $(sc.name)"
            )
            @test saved_id == run_id
        end

        # 3. Query archived simulation runs list
        runs_df = list_simulation_runs(db)
        @test nrow(runs_df) == 3
        @test "historical" in runs_df.scenario
        @test "ssp245" in runs_df.scenario
        @test "ssp585" in runs_df.scenario

        # Test scenario filtering
        ssp245_df = list_simulation_runs(db, scenario = "ssp245")
        @test nrow(ssp245_df) == 1
        @test ssp245_df.projection_year[1] == 2050

        # 4. Load trajectories DataFrame and test filtering
        traj_df = load_trajectories_df(db, "test_run_ssp245_2050", max_particles = 3)
        @test nrow(traj_df) == 3 * 37 # 3 particles x 37 time steps
        @test "particle_id" in names(traj_df)
        @test "lon" in names(traj_df)
        @test "stage" in names(traj_df)

        # 5. Load connectivity matrix
        conn_loaded = load_connectivity_matrix(db, "test_run_ssp245_2050")
        @test size(conn_loaded.matrix) == (4, 4)
        @test length(conn_loaded.strata_names) == 4
        @test all(sum(conn_loaded.matrix, dims = 2) .≈ 1.0)

        # 6. Multi-scenario comparative analytics
        comp_df = compare_scenarios(db)
        @test nrow(comp_df) == 3
        @test "mean_temperature_celsius" in names(comp_df)
        @test "mean_settlement_success" in names(comp_df)

        # 7. Ensemble Model Averaging
        ens = compute_ensemble_model_average(db, ["historical", "ssp245", "ssp585"])
        @test size(ens.mean_connectivity) == (4, 4)
        @test size(ens.std_connectivity) == (4, 4)
        @test all(sum(ens.mean_connectivity, dims = 2) .≈ 1.0)
        @test length(ens.weights) == 3
        @test sum(ens.weights) ≈ 1.0
        @test ens.mean_thermal_exposure isa Float64

        # Test custom weighted model averaging
        ens_weighted = compute_ensemble_model_average(
            db,
            ["historical", "ssp245", "ssp585"],
            weights = [0.2, 0.5, 0.3]
        )
        @test sum(ens_weighted.weights) ≈ 1.0
        @test ens_weighted.weights[2] ≈ 0.5

        # 8. Test NamedTuple Trajectory Reconstruction
        trajs_rec = load_trajectories_namedtuple(db, "test_run_ssp245_2050")
        @test size(trajs_rec.lons) == (10, 37)
        @test size(trajs_rec.lats) == (10, 37)
        @test size(trajs_rec.depths) == (10, 37)
        @test length(trajs_rec.alive) == 10
        @test length(trajs_rec.times) == 37
        @test length(trajs_rec.ids) == 10

        # Test loading all scenarios dictionary
        all_scens = load_all_scenario_trajectories(db)
        @test length(all_scens) == 3
        @test :ssp245 in keys(all_scens)

        # 9. Test Hydrodynamic 3D Field Persistence and Retrieval
        glons = [-64.0, -63.0, -62.0]
        glats = [43.0, 44.0]
        gdepths = [-10.0, -50.0]
        u_dummy = fill(0.08, 3, 2, 2)
        v_dummy = fill(-0.04, 3, 2, 2)
        w_dummy = fill(0.001, 3, 2, 2)
        t_dummy = fill(5.5, 3, 2, 2)
        s_dummy = fill(32.8, 3, 2, 2)

        save_hydrodynamic_field!(
            db,
            "test_run_ssp245_2050",
            (scenario = "ssp245", projection_year = 2050);
            grid_lons = glons,
            grid_lats = glats,
            grid_depths = gdepths,
            u = u_dummy,
            v = v_dummy,
            w = w_dummy,
            temperature = t_dummy,
            salinity = s_dummy,
            time_seconds = 3600.0
        )

        hydro_loaded = load_hydrodynamic_field(db, "test_run_ssp245_2050", time_seconds = 3600.0)
        @test size(hydro_loaded.u) == (3, 2, 2)
        @test hydro_loaded.u[1, 1, 1] ≈ 0.08
        @test hydro_loaded.temperature[1, 1, 1] ≈ 5.5

        # 10. Test Gridded Dispersal Loader
        disp_loaded = load_gridded_dispersal(db, "test_run_ssp245_2050")
        @test length(disp_loaded.lon_centers) > 0
        @test length(disp_loaded.lat_centers) > 0

        # 11. Test Multi-Scenario Multi-Layer Interactive Map Export
        multi_html_path = "outputs/test_interactive_multiscenario.html"
        scenarios_mock = Dict(
            "Historical (2015)" => (trajectories = all_scens[:historical], gridded_dispersal = disp_loaded, connectivity = conn_loaded),
            "SSP2-4.5 (2050)" => (trajectories = all_scens[:ssp245], gridded_dispersal = disp_loaded, connectivity = conn_loaded)
        )
        export_interactive_tracks_html(
            multi_html_path;
            scenarios_data = scenarios_mock,
            strata_definitions = cfa_defs,
            title = "Test Multi-Scenario Dashboard"
        )
        @test isfile(multi_html_path)
        html_str = read(multi_html_path, String)
        @test occursin("Historical (2015)", html_str)
        @test occursin("SSP2-4.5 (2050)", html_str)
        @test occursin("layer-density", html_str)
        @test occursin("layer-hydro-advection", html_str)

        # 12. Export DuckDB tables to Apache Parquet
        parquet_dir = "outputs/test_parquet"
        p_files = export_duckdb_to_parquet(db, parquet_dir)
        @test length(p_files) >= 5
        for pf in p_files
            @test isfile(pf)
        end

        # Clean up database connection
        close_duckdb_storage(db)
    end

    @testset "15. Centralized Configuration File Management & Scenario Metadata" begin
        config_path = "inputs/ParticleTracking.config"
        @test isfile(config_path)

        # 1. Load centralized configuration and verify sections
        cfg = load_configuration(config_path)
        @test haskey(cfg, "domain")
        @test haskey(cfg, "grid")
        @test haskey(cfg, "data")
        @test haskey(cfg, "tides")
        @test haskey(cfg, "climate")
        @test haskey(cfg, "hydrodynamics")
        @test haskey(cfg, "biology")
        @test haskey(cfg, "dvm")
        @test haskey(cfg, "molting_and_settlement")
        @test haskey(cfg, "storage")
        @test haskey(cfg, "hardware")
        @test haskey(cfg, "visualization")
        @test haskey(cfg, "paths")

        def_cfg = get_default_configuration()
        @test def_cfg["domain"]["lon_min"] == -71.0
        @test def_cfg["grid"]["nx"] == 50
        @test def_cfg["biology"]["min_seabed_depth"] == 100.0
        @test def_cfg["climate"]["scenario"] == "ssp245"

        @test haskey(cfg["domain"], "lon_min")
        @test haskey(cfg["grid"], "nx")
        @test haskey(cfg["biology"], "min_seabed_depth")
        @test haskey(cfg["climate"], "scenario")

        # 2. Test configuration to options conversion
        opts = configuration_to_options(cfg, n_particles = 42, scenario = :ssp585, min_seabed_depth = 100.0)
        @test opts.n_particles == 42
        @test opts.scenario == :ssp585
        @test opts.min_seabed_depth == 100.0

        # 3. Test options to configuration conversion
        cfg_out = options_to_configuration(opts)
        @test cfg_out["biology"]["n_particles"] == 42
        @test cfg_out["climate"]["scenario"] == "ssp585"

        # 4. Test save and re-load configuration
        test_save_path = "outputs/test_saved.config"
        save_configuration(cfg_out, test_save_path)
        @test isfile(test_save_path)
        cfg_reloaded = load_configuration(test_save_path)
        @test cfg_reloaded["biology"]["n_particles"] == 42

        # 5. Test archiving configuration into DuckDB and reloading it
        test_db_path = "outputs/test_config_meta.duckdb"
        if isfile(test_db_path)
            rm(test_db_path, force = true)
        end
        db = open_duckdb_storage(test_db_path)
        rng = MersenneTwister(42)
        larvae = initialize_larval_particles(5, rng = rng)
        trajs = track_larval_cohort(
            larvae,
            velocity_fn = (x, y, z, t) -> (0.05, 0.02, 0.0),
            total_duration = 3600.0,
            dt = 600.0,
            rng = rng
        )

        run_id = "test_config_run"
        save_simulation_run!(db, run_id, opts; trajectories = trajs, config = cfg_out)

        # Query metadata table directly
        df_run = DataFrame(DBInterface.execute(db, "SELECT config_toml FROM simulation_runs WHERE run_id = '$(run_id)';"))
        @test nrow(df_run) == 1
        @test !ismissing(df_run.config_toml[1])
        @test occursin("ssp585", df_run.config_toml[1])

        # Load run configuration helper
        cfg_loaded = load_run_configuration(db, run_id)
        @test cfg_loaded["biology"]["n_particles"] == 42
        @test cfg_loaded["climate"]["scenario"] == "ssp585"

        close_duckdb_storage(db)

        # 6. Test NetCDF and JLD2 configuration attribute persistence
        nc_test = "outputs/test_config_meta.nc"
        export_larval_dispersal_netcdf(nc_test, trajectories = trajs, config = cfg_out)
        @test isfile(nc_test)
        NCDataset(nc_test, "r") do ds
            @test haskey(ds.attrib, "configuration")
            @test occursin("ssp585", ds.attrib["configuration"])
        end

        jld_test = "outputs/test_config_meta.jld2"
        export_larval_dispersal_jld2(jld_test, trajectories = trajs, config = cfg_out)
        @test isfile(jld_test)
        jld_data = JLD2.load(jld_test)
        @test haskey(jld_data, "configuration")
        @test jld_data["configuration"]["biology"]["n_particles"] == 42
    end

    @testset "16. Coastline Geometry, 0% Land Seeding & CFA Intersection" begin
        # 1. Terrestrial reference points
        @test is_point_on_land(-63.58, 44.65) == true  # Halifax NS
        @test is_point_on_land(-63.28, 45.36) == true  # Truro NS
        @test is_point_on_land(-65.20, 44.40) == true  # Kejimkujik NS
        @test is_point_on_land(-66.12, 43.83) == true  # Yarmouth NS
        @test is_point_on_land(-60.19, 46.13) == true  # Sydney CB
        @test is_point_on_land(-60.43, 47.00) == true  # Cape North CB
        @test is_point_on_land(-63.13, 46.24) == true  # Charlottetown PEI
        @test is_point_on_land(-66.06, 45.27) == true  # Saint John NB
        @test is_point_on_land(-66.64, 45.96) == true  # Fredericton NB
        @test is_point_on_land(-59.90, 43.93) == true  # Sable Island
        @test is_point_on_land(-57.60, 47.65) == true  # Burgeo NL

        # 2. Marine reference points
        @test is_point_on_land(-62.80, 43.90) == false # Emerald Basin
        @test is_point_on_land(-61.50, 43.50) == false # Western Bank
        @test is_point_on_land(-65.50, 42.80) == false # Browns Bank
        @test is_point_on_land(-58.50, 44.50) == false # Banquereau Bank
        @test is_point_on_land(-59.50, 46.80) == false # Sydney Bight
        @test is_point_on_land(-65.80, 45.00) == false # Bay of Fundy
        @test is_point_on_land(-58.00, 46.00) == false # Laurentian Channel

        # 3. Coastline I/O
        coast = load_coastline_polygons("inputs/coastline.dat")
        @test length(coast) >= 6
        @test any(p -> p.code == :nova_scotia_mainland, coast)
        @test any(p -> p.code == :cape_breton, coast)
        @test any(p -> p.code == :prince_edward_island, coast)

        # 4. Larval Particle Seeding - 0% Land Guarantee
        rng = MersenneTwister(42)
        bathy_path = "inputs/bathymetry_active.nc"
        @test isfile(bathy_path)

        n_test = 500
        larvae = initialize_larval_particles(
            n_test,
            lon_range = (-66.0, -60.0),
            lat_range = (43.0, 46.0),
            depth_range = (-50.0, -10.0),
            min_seabed_depth = 50.0,
            bathymetry = bathy_path,
            rng = rng
        )
        @test length(larvae.lon) == n_test
        land_particles = count(p -> is_point_on_land(larvae.lon[p], larvae.lat[p]), 1:n_test)
        @test land_particles == 0

        # 5. CFA Stratum Seeding
        for cfa_sym in [:cfa4x, :cfasouth, :cfanorth]
            larvae_cfa = initialize_larval_particles(
                50,
                stratum = cfa_sym,
                min_seabed_depth = 40.0,
                bathymetry = bathy_path,
                rng = rng
            )
            @test length(larvae_cfa.lon) == 50
            @test count(p -> is_point_on_land(larvae_cfa.lon[p], larvae_cfa.lat[p]), 1:50) == 0
        end

        # 6. CFA Coastline Intersection & Marine Demarcation
        cfas = load_cfa_polygons("inputs", intersect_coastline = true)
        @test length(cfas) == 3
        cfa4x = cfas[findfirst(p -> p.code == :cfa4x, cfas)]
        cfasouth = cfas[findfirst(p -> p.code == :cfasouth, cfas)]
        cfanorth = cfas[findfirst(p -> p.code == :cfanorth, cfas)]

        # Terrestrial points NOT in any CFA
        @test point_in_polygon(-63.58, 44.65, cfa4x.lons, cfa4x.lats) == false
        @test point_in_polygon(-63.58, 44.65, cfasouth.lons, cfasouth.lats) == false
        @test point_in_polygon(-63.28, 45.36, cfa4x.lons, cfa4x.lats) == false
        @test point_in_polygon(-60.19, 46.13, cfanorth.lons, cfanorth.lats) == false

        # Marine points ARE in respective CFAs
        @test point_in_polygon(-65.50, 42.80, cfa4x.lons, cfa4x.lats) == true
        @test point_in_polygon(-62.80, 43.90, cfasouth.lons, cfasouth.lats) == true
        @test point_in_polygon(-59.50, 46.80, cfanorth.lons, cfanorth.lats) == true

        # 7. Spatial Domain Buffering & Envelope Expansion
        dlon_100, dlat_100 = buffer_distance_to_degrees(100.0, 44.5)
        @test isapprox(dlat_100, 100.0 / 111.195, atol = 0.02)
        @test isapprox(dlon_100, dlat_100 / cos(deg2rad(44.5)), atol = 0.02)

        buf_lon, buf_lat = expand_domain_with_buffer((-64.0, -60.0), (43.0, 45.0), buffer_km = 100.0)
        @test buf_lon[1] < -64.0 && buf_lon[2] > -60.0
        @test buf_lat[1] < 43.0 && buf_lat[2] > 45.0

        cfa_env = get_strata_buffered_envelope(cfas, buffer_km = 100.0)
        @test cfa_env.buffer_km == 100.0
        @test cfa_env.lon_range[1] < cfa_env.raw_lon_range[1]
        @test cfa_env.lon_range[2] > cfa_env.raw_lon_range[2]

        # 8. Seeding with 100 km Buffer (0% Land Placement)
        larvae_buf = initialize_larval_particles(
            300,
            stratum = :cfasouth,
            buffer_km = 100.0,
            min_seabed_depth = 50.0,
            bathymetry = bathy_path,
            rng = rng
        )
        @test length(larvae_buf.lon) == 300
        @test count(p -> is_point_on_land(larvae_buf.lon[p], larvae_buf.lat[p]), 1:300) == 0
    end

end
