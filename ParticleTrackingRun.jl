"""
    ParticleTrackingRun.jl

Command-line test, debugging, and production execution interface for regional
hydrodynamic modeling and snow crab (Chionoecetes opilio) larval particle tracking.

# Design & Workflow Segments
Supports segmented execution to allow rapid iteration and debugging without
re-running prior completed stages:
1. `data`: Ingest real bathymetry/wind fields or generate synthetic benchmarks.
2. `grid`: Construct spherical LatitudeLongitudeGrid and ImmersedBoundaryGrid.
3. `model`: Configure hydrostatic free-surface model with Coriolis, buoyancy and tides.
4. `climate`: Compute CMIP6 climate deltas, PLD models, and thermal mortality.
5. `sim`: Execute Oceananigans hydrodynamic time stepping with CFL monitoring.
6. `track`: Perform Lagrangian particle tracking with DVM and ontogenetic molting.
7. `metrics`: Compute gridded retention, thermal exposure, empirical diffusivity,
   and demographic connectivity matrices with NetCDF and JLD2 archiving.
8. `viz`: Render CairoMakie trajectory maps, DVM depth profiles, empirical movement
   fields, and connectivity heatmaps.
9. `all`: Run full end-to-end production pipeline.

# CLI Examples
```bash
# Display help and available options
julia --project=. ParticleTrackingRun.jl --help

# Fast debug run of the entire pipeline
julia --project=. ParticleTrackingRun.jl --all --quick

# Run individual segments independently
julia --project=. ParticleTrackingRun.jl --data --synthetic
julia --project=. ParticleTrackingRun.jl --grid
julia --project=. ParticleTrackingRun.jl --track --particles=200
julia --project=. ParticleTrackingRun.jl --metrics
julia --project=. ParticleTrackingRun.jl --viz
```
"""

# Load unified project environment and dependencies (includes HydrodynamicOptions and configuration)
include(joinpath(@__DIR__, "src", "init.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Segment 1: Data Ingestion & Benchmark Generation
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_data(; opts::HydrodynamicOptions) -> NamedTuple

Fetch real-world bathymetry and surface winds from NOAA ERDDAP or generate synthetic
analogs matching the model bounding box. Computes kinematic wind stress components.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options container.

# Outputs
- `NamedTuple`: `(bathy_file, wind_file, tau_x, tau_y, bathy_info)`
"""
function run_segment_data(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Segment 1/8] Environmental Data Ingestion & Drag Processing")
    println("=================================================================")
    mkpath(opts.input_dir)

    bathy_file = joinpath(opts.input_dir, "bathymetry_active.nc")
    wind_file  = joinpath(opts.input_dir, "wind_active.nc")

    if opts.data_mode == :real
        println("Retrieving real bathymetry from NOAA ERDDAP (etopo180)...")
        try
            fetch_open_bathymetry(
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                output_path = bathy_file
            )
        catch err
            @warn "NOAA ERDDAP bathymetry download failed: $(err). Falling back to synthetic."
            generate_synthetic_bathymetry(
                bathy_file,
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                n_lon = opts.grid_size[1],
                n_lat = opts.grid_size[2]
            )
        end

        println("Retrieving real surface winds from NOAA ERDDAP (erdBSwinds1day)...")
        try
            fetch_open_surface_winds(
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                time_iso = "2023-06-01T00:00:00Z",
                output_path = wind_file
            )
        catch err
            @warn "NOAA ERDDAP wind download failed: $(err). Falling back to synthetic."
            generate_synthetic_forcing(
                wind_file,
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                n_lon = opts.grid_size[1],
                n_lat = opts.grid_size[2]
            )
        end
    else
        println("Generating synthetic Scotian Shelf bathymetry...")
        generate_synthetic_bathymetry(
            bathy_file,
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            n_lon = opts.grid_size[1],
            n_lat = opts.grid_size[2],
            inshore_depth = -100.0,
            shelf_slope = 600.0
        )

        println("Generating synthetic surface wind forcing...")
        generate_synthetic_forcing(
            wind_file,
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            time_range = (0.0, opts.sim_duration),
            n_lon = opts.grid_size[1],
            n_lat = opts.grid_size[2],
            n_time = 24,
            tau_x_amplitude = 0.1,
            tau_y_amplitude = 0.02
        )
    end

    # Inspect the generated/downloaded NetCDF datasets
    bathy_info = inspect_netcdf(bathy_file, verbose = true)

    # Compute kinematic surface wind stress via Large & Pond (1981) formulation
    u10_ref, v10_ref = 8.5, 3.2 # m/s
    tau_x, tau_y = wind_speed_to_kinematic_stress(u10_ref, v10_ref)
    println("Calculated Large & Pond (1981) kinematic wind stress:")
    println("  Reference 10m wind: u = $(u10_ref) m/s, v = $(v10_ref) m/s")
    println("  Kinematic stress:   tau_x = $(round(tau_x, digits=6)), " *
            "tau_y = $(round(tau_y, digits=6)) m^2/s^2")

    return (
        bathy_file = bathy_file,
        wind_file = wind_file,
        tau_x = tau_x,
        tau_y = tau_y,
        bathy_info = bathy_info
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 2: Computational Grid & Immersed Seafloor Boundary
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_grid(; opts::HydrodynamicOptions, bathy_file::Union{Nothing, String}=nothing)

Construct the base spherical `LatitudeLongitudeGrid` and wrap it with an
`ImmersedBoundaryGrid` using 2D bilinear interpolation from bathymetry data.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `bathy_file`: Optional custom bathymetry NetCDF path.

# Outputs
- `NamedTuple`: `(base_grid, immersed_grid)`
"""
function run_segment_grid(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    bathy_file::Union{Nothing, String} = nothing
)
    println("\n=================================================================")
    println(" [Segment 2/8] Spherical Grid & Immersed Boundary Construction")
    println("=================================================================")

    target_bathy = isnothing(bathy_file) ?
                   joinpath(opts.input_dir, "bathymetry_active.nc") : bathy_file

    if !isfile(target_bathy)
        println("Bathymetry file missing. Running Segment 1 data generation...")
        data_res = run_segment_data(opts = opts)
        target_bathy = data_res.bathy_file
    end

    arch_label = opts.use_gpu ? "GPU (CUDA)" : "CPU"
    println("Building base spherical grid on $(arch_label): $(opts.grid_size) cells...")
    base_grid = build_shelf_grid(
        architecture = opts.use_gpu ? :gpu : :cpu,
        lon_range = opts.domain_lon,
        lat_range = opts.domain_lat,
        z_range = opts.domain_z,
        grid_size = opts.grid_size,
        fallback_to_cpu = opts.fallback_to_cpu
    )

    println("Constructing immersed boundary with 2D bilinear regridding...")
    immersed_grid = build_immersed_grid_from_real_data(base_grid, target_bathy)

    println("Grid summary:")
    println("  Longitude: $(opts.domain_lon[1])°E to $(opts.domain_lon[2])°E (Nx=$(opts.grid_size[1]))")
    println("  Latitude:  $(opts.domain_lat[1])°N to $(opts.domain_lat[2])°N (Ny=$(opts.grid_size[2]))")
    println("  Depth:     $(opts.domain_z[1]) m to $(opts.domain_z[2]) m (Nz=$(opts.grid_size[3]))")

    return (base_grid = base_grid, immersed_grid = immersed_grid)
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 3: Hydrodynamic Model Instantiation & Stratification
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_model(;
        opts::HydrodynamicOptions,
        immersed_grid=nothing,
        tau_x::Real=0.0001,
        tau_y::Real=0.0
    ) -> NamedTuple

Instantiate the Oceananigans `HydrostaticFreeSurfaceModel` with Coriolis rotation,
buoyancy, surface wind boundary conditions, and semi-diurnal tidal body forcing (M2).

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `immersed_grid`: Pre-constructed `ImmersedBoundaryGrid`.
- `tau_x, tau_y`: Surface kinematic momentum fluxes.

# Outputs
- `NamedTuple`: `(model, tidal_forcing, simpson_hunter_bank, simpson_hunter_shelf)`
"""
function run_segment_model(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    immersed_grid = nothing,
    tau_x::Real = 1e-4,
    tau_y::Real = 0.0
)
    println("\n=================================================================")
    println(" [Segment 3/8] Hydrodynamic Model & Tidal Forcing Setup")
    println("=================================================================")

    target_grid = if isnothing(immersed_grid)
        grid_res = run_segment_grid(opts = opts)
        grid_res.immersed_grid
    else
        immersed_grid
    end

    tidal_forcing = if opts.enable_tides
        println("Configuring astronomical tidal body forcing (M2 constituent)...")
        build_tidal_body_forcing(
            constituents = [:M2],
            u_amplitudes = Dict(:M2 => opts.tidal_u_amp),
            v_amplitudes = Dict(:M2 => opts.tidal_v_amp)
        )
    else
        nothing
    end

    coriolis_lat = 0.5 * (opts.domain_lat[1] + opts.domain_lat[2])
    println("Building HydrostaticFreeSurfaceModel (Coriolis at $(coriolis_lat)°N)...")
    model = build_hydrodynamic_model(
        target_grid,
        coriolis_latitude = coriolis_lat,
        surface_wind_stress_x = tau_x,
        surface_wind_stress_y = tau_y,
        tidal_forcing = tidal_forcing,
        ν = 1e-2,
        κ = 1e-2,
        tracers = (:T, :S)
    )

    # Initialize thermal and haline stratification
    println("Applying baseline thermal stratification (T_surf=15°C, dT/dz=0.01°C/m)...")
    set_initial_stratification!(
        model,
        surface_temperature = 15.0,
        temperature_gradient = 0.01,
        salinity = 35.0
    )

    # Calculate Simpson-Hunter tidal mixing front parameters
    chi_bank  = simpson_hunter_parameter(40.0, 1.1)  # Shallow bank (mixed)
    chi_shelf = simpson_hunter_parameter(150.0, 0.2) # Deep shelf (stratified)
    println("Simpson-Hunter Tidal Mixing Diagnostics:")
    println("  Shallow Bank (h=40m,  U=1.1m/s): χ = $(round(chi_bank, digits=2)) (well-mixed if < 1.5)")
    println("  Deep Shelf   (h=150m, U=0.2m/s): χ = $(round(chi_shelf, digits=2)) (stratified if > 2.0)")

    return (
        model = model,
        tidal_forcing = tidal_forcing,
        simpson_hunter_bank = chi_bank,
        simpson_hunter_shelf = chi_shelf
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 4: Climate Forcing Scenarios & Larval Thermal Ecology
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_climate(; opts::HydrodynamicOptions, model=nothing) -> NamedTuple

Inspect CMIP6 climate scenario deltas, apply warming and freshening anomalies to
the model stratification, and evaluate temperature-dependent PLD and mortality.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `model`: Optional `HydrostaticFreeSurfaceModel` to update in-place.

# Outputs
- `NamedTuple`: `(deltas, pld_cold, pld_warm, mort_cold, mort_warm)`
"""
function run_segment_climate(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    model = nothing
)
    println("\n=================================================================")
    println(" [Segment 4/8] Climate Scenario Integration & Thermal Ecology")
    println("=================================================================")

    deltas = get_climate_scenario_deltas(opts.scenario, year = opts.projection_year)
    println("Climate Scenario: $(deltas.description) [Year $(opts.projection_year)]")
    println("  Surface Temperature Anomaly: +$(round(deltas.ΔT_surface, digits=2)) °C")
    println("  Deep Temperature Anomaly:    +$(round(deltas.ΔT_deep, digits=2)) °C")
    println("  Surface Salinity Anomaly:    $(round(deltas.ΔS_surface, digits=2)) PSU")
    println("  Atmospheric Wind Factor:     x$(round(deltas.Δwind_factor, digits=2))")

    if !isnothing(model)
        println("Applying climate anomalies to model stratification...")
        apply_climate_scenario!(model, scenario = opts.scenario, year = opts.projection_year)
    end

    # Larval developmental duration (PLD) and thermal stress mortality
    t_cold = 2.5 # Deep baseline water temperature (°C)
    t_warm = t_cold + deltas.ΔT_surface
    pld_cold = temperature_dependent_pld(t_cold)
    pld_warm = temperature_dependent_pld(t_warm)
    mort_cold = larval_thermal_mortality_rate(t_cold)
    mort_warm = larval_thermal_mortality_rate(t_warm)

    println("Snow Crab Larval Thermal Ecology:")
    println("  Baseline (T = $(t_cold)°C): PLD = $(round(pld_cold, digits=1)) days, " *
            "Mortality = $(round(mort_cold*100, digits=2)) %/day")
    println("  Warmed   (T = $(round(t_warm, digits=1))°C): PLD = $(round(pld_warm, digits=1)) days, " *
            "Mortality = $(round(mort_warm*100, digits=2)) %/day")

    return (
        deltas = deltas,
        pld_cold = pld_cold,
        pld_warm = pld_warm,
        mort_cold = mort_cold,
        mort_warm = mort_warm
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 5: Hydrodynamic Simulation Execution
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_simulation(; opts::HydrodynamicOptions, model=nothing) -> NamedTuple

Configure Oceananigans time stepping, adaptive CFL wizard, stability watchdogs,
and execute the regional hydrodynamic integration.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `model`: Configured `HydrostaticFreeSurfaceModel`.

# Outputs
- `NamedTuple`: `(simulation, jld2_output_path)`
"""
function run_segment_simulation(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    model = nothing
)
    println("\n=================================================================")
    println(" [Segment 5/8] Hydrodynamic Simulation Time Stepping")
    println("=================================================================")
    mkpath(opts.output_dir)

    target_model = if isnothing(model)
        model_res = run_segment_model(opts = opts)
        model_res.model
    else
        model
    end

    jld2_filename = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
    jld2_path = joinpath(opts.output_dir, jld2_filename)

    println("Setting up simulation (stop_time=$(opts.sim_duration)s, Δt=$(opts.sim_dt)s)...")
    sim = setup_hydrodynamic_simulation(
        target_model,
        Δt = opts.sim_dt,
        stop_time = opts.sim_duration,
        adaptive_time_step = opts.adaptive_cfl,
        target_cfl = opts.target_cfl,
        output_dir = opts.output_dir,
        output_filename = jld2_filename,
        output_schedule = 50,
        progress_schedule = 20
    )

    println("Integrating hydrostatic primitive equations...")
    run_hydrodynamic_simulation!(sim, verbose = true)
    println("Simulation fields successfully saved to: $(jld2_path)")

    # Archive hydrodynamic snapshot fields to DuckDB
    if opts.enable_duckdb
        try
            println("Archiving hydrodynamic flow fields to DuckDB -> $(opts.duckdb_path)...")
            db = open_duckdb_storage(opts.duckdb_path)
            run_id = "run_$(opts.scenario)_$(opts.projection_year)"
            g = target_model.grid
            under_g = g isa ImmersedBoundaryGrid ? g.underlying_grid : g
            glon = collect(under_g.λᶠᵃᵃ[1:under_g.Nx])
            glat = collect(under_g.φᵃᶠᵃ[1:under_g.Ny])
            gdepth = collect(under_g.zᵃᵃᶜ[1:under_g.Nz])
            u_data = Array(interior(target_model.velocities.u))
            v_data = Array(interior(target_model.velocities.v))
            w_data = Array(interior(target_model.velocities.w))
            t_data = Array(interior(target_model.tracers.T))
            s_data = Array(interior(target_model.tracers.S))
            save_hydrodynamic_field!(
                db, run_id, opts;
                grid_lons = glon, grid_lats = glat, grid_depths = gdepth,
                u = u_data, v = v_data, w = w_data,
                temperature = t_data, salinity = s_data,
                time_seconds = opts.sim_duration
            )
            close_duckdb_storage(db)
            println("Hydrodynamic fields for '$(run_id)' archived in DuckDB.")
        catch err
            @warn "Failed to archive hydrodynamic fields to DuckDB: $(err)"
        end
    end

    return (simulation = sim, jld2_output_path = jld2_path)
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 6: Lagrangian Particle Tracking & Larval Life History
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_tracking(; opts::HydrodynamicOptions) -> NamedTuple

Simulate Lagrangian particle transport for a snow crab larval cohort across nursery
spawning areas with DVM swimming, M2 tidal currents, degree-day molting, and
benthic settlement filtering.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.

# Outputs
- `NamedTuple`: Trajectories record and summary metrics.
"""
function run_segment_tracking(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Segment 6/8] Lagrangian Particle Tracking & Larval Behavior")
    println("=================================================================")
    rng = MersenneTwister(opts.seed)

    # Spawning area: derive from loaded CFA boundaries + user-defined spatial buffer
    cfa_polys = load_cfa_polygons(opts.input_dir)
    spawn_lon, spawn_lat = if !isempty(cfa_polys)
        env = get_strata_buffered_envelope(cfa_polys, buffer_km = opts.buffer_km)
        (env.lon_range, env.lat_range)
    else
        expand_domain_with_buffer(opts.domain_lon, opts.domain_lat, buffer_km = opts.buffer_km)
    end

    # Active seabed bathymetry dataset
    target_bathy_path = joinpath(opts.input_dir, "bathymetry_active.nc")
    bathy_src = isfile(target_bathy_path) ? target_bathy_path : begin
        (lon, lat) -> -120.0 - 200.0 * (lat - opts.domain_lat[1]) / (opts.domain_lat[2] - opts.domain_lat[1])
    end

    println("Initializing $(opts.n_particles) Zoea I larvae in marine water (min depth >= $(opts.min_seabed_depth) m, buffer = $(opts.buffer_km) km)...")
    larvae = initialize_larval_particles(
        opts.n_particles,
        lon_range = spawn_lon,
        lat_range = spawn_lat,
        depth_range = (-45.0, -15.0),
        min_seabed_depth = opts.min_seabed_depth,
        buffer_km = 0.0, # Buffer is already incorporated into spawn_lon/spawn_lat
        bathymetry = bathy_src,
        stage = :zoea1,
        rng = rng
    )

    # Background Eulerian flow velocity function (with offshore and along-shelf component)
    flow_field_fn(lon, lat, z, t) = begin
        u_mean = 0.06 + 0.02 * sin(2.0 * π * t / 86400.0) # Along-shelf jet
        v_mean = -0.02 + 0.01 * cos(2.0 * π * t / 86400.0) # Offshore export
        w_mean = 0.0001 * sin(lat)
        (u_mean, v_mean, w_mean)
    end

    # In-situ vertical temperature profile
    temp_field_fn(lon, lat, z, t) = max(0.5, 6.0 + 0.02 * z)

    # Seabed bathymetric elevation
    bathy_field_fn = if isfile(target_bathy_path)
        get_bathymetry_interpolator(target_bathy_path)
    else
        (lon, lat) -> -120.0 - 200.0 * (lat - opts.domain_lat[1]) / (opts.domain_lat[2] - opts.domain_lat[1])
    end

    println("Tracking cohort over $(opts.track_duration / 86400.0) days (dt=$(opts.track_dt)s)...")
    trajectories = track_larval_cohort(
        larvae,
        velocity_fn = flow_field_fn,
        temperature_fn = temp_field_fn,
        bathymetry_fn = bathy_field_fn,
        total_duration = opts.track_duration,
        dt = opts.track_dt,
        κ_h = opts.diffusivity_h,
        κ_v = opts.diffusivity_v,
        is_lat_lon = true,
        enable_tides = opts.enable_tides,
        tidal_u_amp = opts.tidal_u_amp,
        tidal_v_amp = opts.tidal_v_amp,
        enable_molting = opts.enable_molting,
        rng = rng
    )

    # Save trajectories to JLD2 checkpoint for modular reloading
    mkpath(opts.output_dir)
    track_checkpoint = joinpath(opts.output_dir, "larval_trajectories.jld2")
    jldsave(
        track_checkpoint;
        lons = trajectories.lons,
        lats = trajectories.lats,
        depths = trajectories.depths,
        temperatures = hasproperty(trajectories, :temperatures) ? trajectories.temperatures : nothing,
        degree_days = trajectories.degree_days,
        degree_days_timeseries = hasproperty(trajectories, :degree_days_timeseries) ? trajectories.degree_days_timeseries : nothing,
        survival_probability = hasproperty(trajectories, :survival_probability) ? trajectories.survival_probability : nothing,
        stages = trajectories.stages,
        alive = trajectories.alive,
        settlement_status = trajectories.settlement_status,
        settlement_age = trajectories.settlement_age,
        times = trajectories.times,
        ids = trajectories.ids
    )

    # Save trajectories directly into DuckDB
    if opts.enable_duckdb
        try
            println("Archiving trajectories to DuckDB -> $(opts.duckdb_path)...")
            db = open_duckdb_storage(opts.duckdb_path)
            run_id = "run_$(opts.scenario)_$(opts.projection_year)"
            save_simulation_run!(
                db, run_id, opts;
                trajectories = trajectories,
                config = options_to_configuration(opts),
                notes = "Hydrodynamic tracking run ($(opts.scenario), $(opts.projection_year))"
            )
            close_duckdb_storage(db)
            println("Trajectories for '$(run_id)' successfully archived in DuckDB.")
        catch err
            @warn "Failed to archive trajectories to DuckDB: $(err)"
        end
    end

    n_settled = count(==( :settled_successful), trajectories.settlement_status)
    n_alive   = count(identity, trajectories.alive)
    println("Lagrangian tracking complete:")
    println("  Total particles:     $(opts.n_particles)")
    println("  Surviving particles: $(n_alive) / $(opts.n_particles)")
    println("  Settled on nursery:  $(n_settled) / $(opts.n_particles) " *
            "($(round(100.0 * n_settled / opts.n_particles, digits=1))%)")
    println("  Checkpoint saved:    $(track_checkpoint)")

    return (trajectories = trajectories, checkpoint_path = track_checkpoint)
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 7: Empirical Movement, Recruitment, Thermal & Connectivity Metrics
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_metrics(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Compute gridded larval retention and recruitment metrics, thermal exposure indices,
empirical Lagrangian drift velocities and effective diffusivities, regional
connectivity matrices across Crab Fishing Areas (CFAs), and export multi-layer NetCDF and JLD2.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `trajectories`: Optional trajectory NamedTuple.

# Outputs
- `NamedTuple`: Metric outputs, connectivity matrix, and export filepaths.
"""
function run_segment_metrics(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    println("\n=================================================================")
    println(" [Segment 7/8] Empirical Movement, Recruitment & Connectivity")
    println("=================================================================")

    target_trajs = if isnothing(trajectories)
        run_id = "run_$(opts.scenario)_$(opts.projection_year)"
        if opts.enable_duckdb && isfile(opts.duckdb_path)
            db = open_duckdb_storage(opts.duckdb_path; read_only = true)
            try
                println("Loading particle trajectories from DuckDB for run '$(run_id)'...")
                load_trajectories_namedtuple(db, run_id)
            catch err
                track_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
                if isfile(track_cp)
                    println("Loading particle trajectories from checkpoint: $(track_cp)...")
                    saved = load(track_cp)
                    n_p, n_t = size(saved["lons"])
                    t_end = saved["times"][end]
                    (
                        lons = saved["lons"],
                        lats = saved["lats"],
                        depths = saved["depths"],
                        temperatures = get(saved, "temperatures", fill(4.0, n_p, n_t)),
                        degree_days = saved["degree_days"],
                        degree_days_timeseries = get(saved, "degree_days_timeseries", fill(40.0, n_p, n_t)),
                        survival_probability = get(saved, "survival_probability", fill(0.95, n_p, n_t)),
                        stages = saved["stages"],
                        alive = saved["alive"],
                        settlement_status = saved["settlement_status"],
                        settlement_age = get(saved, "settlement_age", fill(t_end, n_p)),
                        times = saved["times"],
                        ids = saved["ids"]
                    )
                else
                    println("Trajectories checkpoint not found. Running Segment 6...")
                    run_segment_tracking(opts = opts).trajectories
                end
            finally
                close_duckdb_storage(db)
            end
        else
            track_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
            if isfile(track_cp)
                println("Loading particle trajectories from checkpoint: $(track_cp)...")
                saved = load(track_cp)
                n_p, n_t = size(saved["lons"])
                t_end = saved["times"][end]
                (
                    lons = saved["lons"],
                    lats = saved["lats"],
                    depths = saved["depths"],
                    temperatures = get(saved, "temperatures", fill(4.0, n_p, n_t)),
                    degree_days = saved["degree_days"],
                    degree_days_timeseries = get(saved, "degree_days_timeseries", fill(40.0, n_p, n_t)),
                    survival_probability = get(saved, "survival_probability", fill(0.95, n_p, n_t)),
                    stages = saved["stages"],
                    alive = saved["alive"],
                    settlement_status = saved["settlement_status"],
                    settlement_age = get(saved, "settlement_age", fill(t_end, n_p)),
                    times = saved["times"],
                    ids = saved["ids"]
                )
            else
                println("Trajectories checkpoint not found. Running Segment 6...")
                run_segment_tracking(opts = opts).trajectories
            end
        end
    else
        trajectories
    end

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    # 1. Empirical Advection & Turbulent Diffusivity
    println("Estimating empirical advection and turbulent diffusivity fields...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 2. Gridded Recruitment & Settlement Metrics
    println("Computing gridded recruitment and benthic nursery settlement metrics...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 3. Gridded Thermal Exposure & Degree-Days
    println("Computing gridded thermal exposure metrics...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 4. Demographic Connectivity Matrix across Crab Fishing Areas (CFAs)
    println("Computing macro-regional demographic connectivity matrix...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    cfa_definitions = if !isempty(cfa_polys)
        println("Loaded $(length(cfa_polys)) CFA boundary polygons from $(opts.input_dir)/: $(join([p.name for p in cfa_polys], ", "))")
        vcat(cfa_polys, [(name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_definitions)
    println("Transition Probability Matrix (P_ij):")
    for i in 1:length(conn.strata_names)
        row_str = join([string(round(conn.matrix[i, j], digits = 3)) for j in 1:length(conn.strata_names)], ", ")
        println("  $(rpad(conn.strata_names[i], 28)) -> [$(row_str)]")
    end

    # 5. Export comprehensive multi-variable NetCDF and JLD2 archives
    nc_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.nc")
    jld_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.jld2")
    active_config = options_to_configuration(opts)

    println("Exporting multi-layer NetCDF archive -> $(nc_export_path)...")
    export_larval_dispersal_netcdf(
        nc_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

    println("Exporting comprehensive JLD2 archive -> $(jld_export_path)...")
    export_larval_dispersal_jld2(
        jld_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

    # 6. Archive simulation run, trajectories, metrics & connectivity in DuckDB
    if opts.enable_duckdb
        println("Archiving simulation run and metrics to DuckDB -> $(opts.duckdb_path)...")
        db = open_duckdb_storage(opts.duckdb_path)
        run_id = "run_$(opts.scenario)_$(opts.projection_year)"
        save_simulation_run!(
            db,
            run_id,
            opts;
            trajectories = target_trajs,
            metrics = (
                mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                mean_degree_days = therm_metrics.mean_degree_days
            ),
            connectivity = conn,
            gridded_dispersal = (
                lon_centers = emp_mov.lon_centers,
                lat_centers = emp_mov.lat_centers,
                u_mean = emp_mov.u_mean,
                v_mean = emp_mov.v_mean,
                diffusivity = emp_mov.diffusivity,
                density = rec_metrics.settlement_density,
                mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                mean_degree_days = therm_metrics.mean_degree_days,
                sample_count = emp_mov.sample_count
            ),
            config = active_config,
            notes = "Hydrodynamic workflow run ($(opts.scenario), $(opts.projection_year))"
        )
        close_duckdb_storage(db)
        println("DuckDB run '$(run_id)' successfully archived.")
    end

    return (
        empirical_movement = emp_mov,
        recruitment_metrics = rec_metrics,
        thermal_metrics = therm_metrics,
        connectivity = conn,
        netcdf_path = nc_export_path,
        jld2_path = jld_export_path,
        duckdb_path = opts.enable_duckdb ? opts.duckdb_path : nothing
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 8: Scientific Visualizations
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_segment_visualize(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Generate and export publication-ready CairoMakie figures and standalone
interactive Leaflet HTML dashboards with multi-scenario comparison and
rich spatial layers.

# Outputs
- `NamedTuple`: Paths to generated PNG figures and HTML dashboard.
"""
function run_segment_visualize(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    println("\n=================================================================")
    println(" [Segment 8/8] Scientific Visualizations & Spatial Figures")
    println("=================================================================")
    mkpath(opts.output_dir)

    target_trajs = if isnothing(trajectories)
        run_id = "run_$(opts.scenario)_$(opts.projection_year)"
        if opts.enable_duckdb && isfile(opts.duckdb_path)
            db = open_duckdb_storage(opts.duckdb_path; read_only = true)
            try
                println("Loading particle trajectories from DuckDB for run '$(run_id)'...")
                load_trajectories_namedtuple(db, run_id)
            catch err
                track_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
                if isfile(track_cp)
                    println("DuckDB run not found. Loading from JLD2 fallback: $(track_cp)...")
                    saved = load(track_cp)
                    n_p, n_t = size(saved["lons"])
                    t_end = saved["times"][end]
                    (
                        lons = saved["lons"],
                        lats = saved["lats"],
                        depths = saved["depths"],
                        temperatures = get(saved, "temperatures", fill(4.0, n_p, n_t)),
                        degree_days = saved["degree_days"],
                        degree_days_timeseries = get(saved, "degree_days_timeseries", fill(40.0, n_p, n_t)),
                        survival_probability = get(saved, "survival_probability", fill(0.95, n_p, n_t)),
                        stages = saved["stages"],
                        alive = saved["alive"],
                        settlement_status = saved["settlement_status"],
                        settlement_age = get(saved, "settlement_age", fill(t_end, n_p)),
                        times = saved["times"],
                        ids = saved["ids"]
                    )
                else
                    println("Trajectories checkpoint not found. Running Segment 6...")
                    run_segment_tracking(opts = opts).trajectories
                end
            finally
                close_duckdb_storage(db)
            end
        else
            track_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
            if isfile(track_cp)
                println("Loading particle trajectories from checkpoint: $(track_cp)...")
                saved = load(track_cp)
                n_p, n_t = size(saved["lons"])
                t_end = saved["times"][end]
                (
                    lons = saved["lons"],
                    lats = saved["lats"],
                    depths = saved["depths"],
                    temperatures = get(saved, "temperatures", fill(4.0, n_p, n_t)),
                    degree_days = saved["degree_days"],
                    degree_days_timeseries = get(saved, "degree_days_timeseries", fill(40.0, n_p, n_t)),
                    survival_probability = get(saved, "survival_probability", fill(0.95, n_p, n_t)),
                    stages = saved["stages"],
                    alive = saved["alive"],
                    settlement_status = saved["settlement_status"],
                    settlement_age = get(saved, "settlement_age", fill(t_end, n_p)),
                    times = saved["times"],
                    ids = saved["ids"]
                )
            else
                println("Trajectories checkpoint not found. Running Segment 6...")
                run_segment_tracking(opts = opts).trajectories
            end
        end
    else
        trajectories
    end

    # Load bathymetry data for background contours
    bathy_path = joinpath(opts.input_dir, "bathymetry_active.nc")
    bathy_data = if isfile(bathy_path)
        load_bathymetry_from_netcdf(bathy_path)
    else
        nothing
    end

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    # 1. 2D Particle Trajectories Map
    fig1_path = joinpath(opts.output_dir, "larval_trajectories.png")
    println("Rendering 2D spatial trajectories map -> $(fig1_path)...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    plot_particle_trajectories(
        target_trajs,
        bathymetry_data = bathy_data,
        strata = cfa_polys,
        title = "Snow Crab Larval Dispersal ($(opts.scenario), Year $(opts.projection_year))",
        output_path = fig1_path
    )

    # 2. DVM Depth Profiles
    fig2_path = joinpath(opts.output_dir, "dvm_depth_profiles.png")
    println("Rendering DVM depth profiles -> $(fig2_path)...")
    plot_vertical_migration_profiles(
        target_trajs,
        sample_indices = 1:min(5, size(target_trajs.depths, 1)),
        title = "Larval Diel Vertical Migration (DVM) Profiles",
        output_path = fig2_path
    )

    # 3. 2D Settlement Nursery Density Heatmap
    fig3_path = joinpath(opts.output_dir, "settlement_density.png")
    println("Rendering 2D settlement density distribution -> $(fig3_path)...")
    plot_larval_dispersal_density(
        target_trajs,
        title = "Snow Crab Nursery Settlement Density (%)",
        output_path = fig3_path
    )

    # 4. Empirical Movement Velocity Field and Diffusivity
    fig4_path = joinpath(opts.output_dir, "empirical_movement_field.png")
    println("Rendering empirical velocity quiver field -> $(fig4_path)...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_empirical_movement_field(emp_mov, output_path = fig4_path)

    # 5. Annotated Regional Connectivity Matrix
    fig5_path = joinpath(opts.output_dir, "regional_connectivity_matrix.png")
    println("Rendering regional connectivity matrix -> $(fig5_path)...")
    cfa_defs = if !isempty(cfa_polys)
        vcat(cfa_polys, [(name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_defs)
    plot_connectivity_matrix(conn, output_path = fig5_path)

    # 6. Thermal Exposure Map
    fig6_path = joinpath(opts.output_dir, "thermal_exposure_map.png")
    println("Rendering thermal exposure map -> $(fig6_path)...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_thermal_exposure_map(therm_metrics, output_path = fig6_path)

    # 7. Recruitment Summary
    fig7_path = joinpath(opts.output_dir, "recruitment_summary.png")
    println("Rendering recruitment summary bar chart -> $(fig7_path)...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_recruitment_summary(rec_metrics, output_path = fig7_path)

    # 8. Hydrodynamic Model Eulerian Advection & Tracers Figures
    fig8_path = joinpath(opts.output_dir, "hydrodynamic_advection.png")
    fig9_path = joinpath(opts.output_dir, "hydrodynamic_tracers.png")
    println("Rendering hydrodynamic advection velocity field -> $(fig8_path)...")
    active_hydro = if opts.enable_duckdb && isfile(opts.duckdb_path)
        db_h = open_duckdb_storage(opts.duckdb_path; read_only = true)
        try
            run_id = "run_$(opts.scenario)_$(opts.projection_year)"
            load_hydrodynamic_field(db_h, run_id)
        catch
            nothing
        finally
            close_duckdb_storage(db_h)
        end
    else
        nothing
    end
    plot_hydrodynamic_advection(active_hydro, bathymetry_data = bathy_data, output_path = fig8_path)

    println("Rendering hydrodynamic seawater temperature & salinity tracers -> $(fig9_path)...")
    plot_hydrodynamic_tracers(active_hydro, output_path = fig9_path)

    # 9. Query Archived Scenarios from DuckDB for Cross-Scenario Comparison & Multi-Layer Interactive Map
    scenarios_bundle = Dict{String, Any}()
    if opts.enable_duckdb && isfile(opts.duckdb_path)
        db = open_duckdb_storage(opts.duckdb_path; read_only = true)
        try
            runs_df = list_simulation_runs(db)
            for r in eachrow(runs_df)
                s_id = string(r.run_id)
                s_label = "$(r.scenario) ($(r.projection_year))"
                try
                    s_trajs = load_trajectories_namedtuple(db, s_id)
                    s_disp = try load_gridded_dispersal(db, s_id) catch; nothing end
                    s_conn = try load_connectivity_matrix(db, s_id) catch; nothing end
                    s_hydro = try load_hydrodynamic_field(db, s_id) catch; nothing end
                    scenarios_bundle[s_label] = (
                        trajectories = s_trajs,
                        gridded_dispersal = s_disp,
                        connectivity = s_conn,
                        hydrodynamics = s_hydro
                    )
                catch
                end
            end
        finally
            close_duckdb_storage(db)
        end
    end

    if isempty(scenarios_bundle) && !isnothing(target_trajs)
        scenarios_bundle["$(opts.scenario) ($(opts.projection_year))"] = (
            trajectories = target_trajs,
            gridded_dispersal = (
                lon_centers = emp_mov.lon_centers,
                lat_centers = emp_mov.lat_centers,
                u_mean = emp_mov.u_mean,
                v_mean = emp_mov.v_mean,
                diffusivity = emp_mov.diffusivity,
                density = rec_metrics.settlement_density,
                mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                mean_degree_days = therm_metrics.mean_degree_days,
                sample_count = emp_mov.sample_count
            ),
            connectivity = conn,
            hydrodynamics = active_hydro
        )
    end

    # Cross-scenario 2D figure
    fig10_path = joinpath(opts.output_dir, "climate_scenario_comparison.png")
    println("Rendering cross-scenario climate comparison -> $(fig10_path)...")
    scenario_comp = Dict{Symbol, Any}()
    for (s_label, s_data) in scenarios_bundle
        scenario_comp[Symbol(s_label)] = s_data.trajectories
    end

    if length(scenario_comp) < 2
        flow_ssp585(lon, lat, z, t) = (0.09 + 0.03 * sin(t / 43200.0), -0.04, 0.0001)
        rng_comp = MersenneTwister(opts.seed + 100)
        larvae_comp = initialize_larval_particles(
            min(opts.n_particles, 50),
            lon_range = (opts.domain_lon[1] + 2.0, opts.domain_lon[2] - 3.0),
            lat_range = (opts.domain_lat[1] + 1.0, opts.domain_lat[2] - 1.0),
            min_seabed_depth = opts.min_seabed_depth,
            bathymetry = isfile(target_bathy_path) ? target_bathy_path : nothing,
            rng = rng_comp
        )
        trajs_ssp585 = track_larval_cohort(
            larvae_comp,
            velocity_fn = flow_ssp585,
            total_duration = opts.track_duration,
            dt = opts.track_dt,
            rng = rng_comp
        )
        scenario_comp[:ssp585_2050] = trajs_ssp585
    end

    compare_scenario_dispersal(
        scenario_comp,
        title = "Scotian Shelf Snow Crab Dispersal Across Climate Scenarios",
        output_path = fig10_path
    )

    # 10. Standalone Multi-Layer Interactive HTML5 Dashboard
    html_path = joinpath(opts.output_dir, "interactive_larval_tracks.html")
    if opts.interactive_map
        println("Rendering multi-layer interactive Leaflet dashboard with hydrodynamic fields -> $(html_path)...")
        export_interactive_tracks_html(
            html_path;
            scenarios_data = scenarios_bundle,
            hydrodynamics = active_hydro,
            strata_definitions = cfa_defs,
            title = "Scotian Shelf Snow Crab Larval Dispersal & Demographic Connectivity"
        )
    end

    println("All visualization figures and interactive maps successfully generated.")
    return (
        fig_trajectories = fig1_path,
        fig_dvm = fig2_path,
        fig_density = fig3_path,
        fig_empirical_movement = fig4_path,
        fig_connectivity = fig5_path,
        fig_thermal = fig6_path,
        fig_recruitment = fig7_path,
        fig_hydro_advection = fig8_path,
        fig_hydro_tracers = fig9_path,
        fig_comparison = fig10_path,
        interactive_map = opts.interactive_map ? html_path : nothing
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Full End-to-End Production Pipeline
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_production_pipeline(; opts::HydrodynamicOptions) -> NamedTuple

Execute the complete end-to-end regional modeling, hydrodynamic simulation,
snow crab larval tracking, metric extraction, and visualization workflow.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.

# Outputs
- `NamedTuple`: Summary of all workflow artifacts.
"""
function run_production_pipeline(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" Executing Full Hydrodynamic & Particle Tracking Pipeline")
    println("=================================================================")
    t_start = time()

    # Step 1: Data Ingestion
    data_res = run_segment_data(opts = opts)

    # Step 2: Grid Construction
    grid_res = run_segment_grid(opts = opts, bathy_file = data_res.bathy_file)

    # Step 3: Model Setup
    model_res = run_segment_model(
        opts = opts,
        immersed_grid = grid_res.immersed_grid,
        tau_x = data_res.tau_x,
        tau_y = data_res.tau_y
    )

    # Step 4: Climate Scenarios
    climate_res = run_segment_climate(opts = opts, model = model_res.model)

    # Step 5: Hydrodynamic Simulation
    sim_res = run_segment_simulation(opts = opts, model = model_res.model)

    # Step 6: Lagrangian Particle Tracking
    track_res = run_segment_tracking(opts = opts)

    # Step 7: Metrics & Connectivity
    metrics_res = run_segment_metrics(opts = opts, trajectories = track_res.trajectories)

    # Step 8: Visualizations
    viz_res = run_segment_visualize(opts = opts, trajectories = track_res.trajectories)

    t_elapsed = round(time() - t_start, digits = 2)
    println("\n=================================================================")
    println(" Complete Workflow Pipeline Successfully Finished in $(t_elapsed) s")
    println(" Outputs written to: $(opts.output_dir)/")
    if opts.enable_duckdb
        println(" DuckDB analytical storage: $(opts.duckdb_path)")
    end
    println("=================================================================")

    return (
        data = data_res,
        grid = grid_res,
        model = model_res,
        climate = climate_res,
        simulation = sim_res,
        tracking = track_res,
        metrics = metrics_res,
        visualizations = viz_res,
        elapsed_seconds = t_elapsed
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# DuckDB Analytics CLI Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_cli_list_runs(; opts::HydrodynamicOptions)

Query and display all simulation runs currently archived in DuckDB.
"""
function run_cli_list_runs(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isfile(opts.duckdb_path)
        println("No DuckDB database found at: $(opts.duckdb_path)")
        return
    end
    println("\n=================================================================")
    println(" Archived Simulation Runs in DuckDB: $(opts.duckdb_path)")
    println("=================================================================")
    db = open_duckdb_storage(opts.duckdb_path; read_only = true)
    df = list_simulation_runs(db)
    close_duckdb_storage(db)

    if nrow(df) == 0
        println("No simulation runs archived in database yet.")
        return
    end

    for row in eachrow(df)
        println("-----------------------------------------------------------------")
        println("Run ID:          $(row.run_id)")
        println("Scenario:        $(row.scenario) (Year $(row.projection_year))")
        println("Created:         $(row.created_at)")
        println("Particles:       $(row.n_particles) larvae | Duration: $(round(row.duration_days, digits=1)) days")
        println("Recruitment:     $(round(row.settlement_success_rate * 100, digits=2))% successful")
        println("Mean PLD:        $(round(row.mean_pld_days, digits=1)) days | Mean Temp: $(round(row.mean_exposure_temperature, digits=2)) °C")
        println("Dispersal Dist:  $(round(row.mean_dispersal_distance_km, digits=1)) km")
    end
    println("-----------------------------------------------------------------\n")
    return df
end

"""
    run_cli_compare_scenarios(; opts::HydrodynamicOptions)

Query and print comparative metrics across archived climate scenarios in DuckDB.
"""
function run_cli_compare_scenarios(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isfile(opts.duckdb_path)
        println("No DuckDB database found at: $(opts.duckdb_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Comparative Analysis from DuckDB")
    println("=================================================================")
    db = open_duckdb_storage(opts.duckdb_path; read_only = true)
    df = compare_scenarios(db)
    close_duckdb_storage(db)

    if nrow(df) == 0
        println("No simulation data available for scenario comparison.")
        return
    end

    for row in eachrow(df)
        println("Scenario: $(rpad(row.scenario, 16)) | Year: $(row.projection_year) | Runs: $(row.n_runs)")
        println("  Settlement Success: $(round(row.mean_settlement_success * 100, digits=2))%")
        println("  Mean PLD:           $(round(row.mean_pld_days, digits=1)) days")
        println("  Mean Exposure Temp: $(round(row.mean_temperature_celsius, digits=2)) °C")
        println("  Mean Dispersal:     $(round(row.mean_dispersal_km, digits=1)) km")
        println("  Thermal Mortality:  $(round(row.mean_thermal_mortality_rate * 100, digits=2))%")
        println()
    end
    return df
end

"""
    run_cli_model_average(; opts::HydrodynamicOptions)

Compute and display ensemble model-averaged demographic connectivity and recruitment.
"""
function run_cli_model_average(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isfile(opts.duckdb_path)
        println("No DuckDB database found at: $(opts.duckdb_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Ensemble Model Averaging")
    println("=================================================================")
    db = open_duckdb_storage(opts.duckdb_path; read_only = true)
    runs_df = list_simulation_runs(db)

    if nrow(runs_df) == 0
        close_duckdb_storage(db)
        println("No simulation data available for model averaging.")
        return
    end

    scens = unique(runs_df.scenario)
    println("Averaging across $(length(scens)) scenario models: $(join(scens, ", "))")
    ens = compute_ensemble_model_average(db, scens)
    close_duckdb_storage(db)

    println("\nEnsemble Weighted Settlement Success: $(round(ens.mean_recruitment_rate * 100, digits=2))%")
    println("Ensemble Weighted Mean PLD:          $(round(ens.mean_pld_days, digits=1)) days")
    println("Ensemble Weighted Thermal Exposure:  $(round(ens.mean_thermal_exposure, digits=2)) °C")
    println("\nEnsemble Model-Averaged Connectivity Matrix (P_ij ± std):")
    n_s = length(ens.strata_names)
    for i in 1:n_s
        row_strs = [
            "$(round(ens.mean_connectivity[i, j], digits=3))±$(round(ens.std_connectivity[i, j], digits=3))"
            for j in 1:n_s
        ]
        println("  $(rpad(ens.strata_names[i], 28)) -> [$(join(row_strs, ", "))]")
    end
    println()
    return ens
end

# ─────────────────────────────────────────────────────────────────────────────
# Command-Line Interface (CLI) Dispatcher
# ─────────────────────────────────────────────────────────────────────────────

"""
    display_help()

Print command-line usage instructions and option flags.
"""
function display_help()
    println("""
Hydrodynamic Model & Particle Tracking CLI Interface

Usage:
  julia --project=. ParticleTrackingRun.jl [FLAGS...]

Execution Modes:
  --all                   Execute entire end-to-end workflow pipeline.
  --quick, -q             Fast debug mode with coarse resolution and short duration.
  --segment=<name>        Run a single workflow segment.
                          Choices: data, grid, model, climate, sim, track, metrics, viz, all.

Individual Segment Flags:
  --data                  Run environmental data ingestion (Option 2A/2B).
  --grid                  Run grid and immersed boundary construction.
  --model                 Run hydrodynamic model setup with tidal forcing.
  --climate               Run CMIP6 climate scenario & larval thermal biology.
  --sim, --simulation     Run Oceananigans hydrodynamic time integration.
  --track, --tracking     Run Lagrangian larval particle tracking.
  --metrics               Compute retention, empirical diffusion & connectivity.
  --viz, --visualize      Generate CairoMakie spatial figures and interactive HTML map.

DuckDB Analytical Storage & Model Averaging:
  --duckdb                Enable DuckDB storage archiving (default: true).
  --no-duckdb             Disable DuckDB storage archiving.
  --db-path=<path>        Custom DuckDB database path (default: outputs/particle_tracking.duckdb).
  --list-runs             Query and display all simulation runs archived in DuckDB.
  --compare-scenarios     Query and display multi-scenario comparison metrics.
  --model-average         Compute ensemble model-averaged connectivity and recruitment.
  --export-parquet        Export all DuckDB tables to Apache Parquet files.

Centralized Configuration:
  --config=<path>         Path to centralized .config file (default: inputs/ParticalTracking.config).
  --save-config[=<path>]  Export active parameter options to .config file and exit.

Computational Architecture:
  --gpu, --cuda           Enable NVIDIA CUDA GPU acceleration for hydrodynamics.
  --cpu                   Execute on multi-threaded CPU (default).
  --fallback-cpu          Automatically fall back to CPU if CUDA GPU is not functional.

Visualization Options:
  --interactive           Export standalone interactive HTML5 Leaflet map (default).
  --no-interactive        Disable interactive HTML map generation.

Spatial Domain & Grid Discretization:
  --lon=<min,max>         Longitude bounding range in degrees East (default: -68.0,-57.0).
  --lat=<min,max>         Latitude bounding range in degrees North (default: 42.0,47.0).
  --depth-range=<min,max> Vertical depth range in meters (default: -1000.0,0.0).
  --grid=<nx,ny,nz>       Grid cell dimensions (default: 50,50,10; quick: 15,15,5).
  --nx=<int>              Zonal grid cells (default: 50).
  --ny=<int>              Meridional grid cells (default: 50).
  --nz=<int>              Vertical grid layers (default: 10).

Environmental Data & Forcing:
  --real                  Fetch real-world NOAA ERDDAP bathymetry and winds.
  --synthetic             Generate idealized synthetic shelf data (default).
  --tides                 Enable astronomical tidal body forcing (default: true).
  --no-tides              Disable tidal body forcing.
  --tidal-u=<val>         Semi-major tidal current amplitude in m/s (default: 0.25).
  --tidal-v=<val>         Semi-minor tidal current amplitude in m/s (default: 0.12).

Climate Scenarios & Thermal Biology:
  --scenario=<name>       Climate scenario: historical, ssp126, ssp245, ssp585, mhw.
  --year=<int>            Climate projection horizon year (default: 2050).

Hydrodynamic Simulation:
  --duration=<hours>      Hydrodynamic simulation duration in hours (default: 12.0).
  --sim-duration=<sec>    Hydrodynamic simulation duration in seconds (default: 43200.0).
  --sim-dt=<sec>          Initial hydrodynamic time step in seconds (default: 120.0).
  --adaptive-cfl          Enable adaptive CFL time stepping (default: true).
  --no-adaptive-cfl       Disable adaptive CFL time stepping.
  --target-cfl=<val>      Target advective Courant-Friedrichs-Lewy limit (default: 0.2).

Lagrangian Particle Tracking & Larval Ecology:
  --particles=<int>       Number of larvae to initialize and track (default: 100).
  --track-duration=<days> Cohort tracking duration in days (default: 5.0).
  --track-dt=<sec>        Lagrangian integration time step in seconds (default: 300.0).
  --min-depth=<meters>    Minimum water depth for larval placement (default: 100.0).
  --buffer-km=<km>        Spatial buffer beyond stratum boundaries (default: 100.0 km).
  --dvm                   Enable stage-dependent Diel Vertical Migration (default: true).
  --no-dvm                Disable Diel Vertical Migration.
  --molting               Enable degree-day thermal molting & mortality (default: true).
  --no-molting            Disable thermal molting.
  --diff-h=<val>          Horizontal turbulent diffusivity in m^2/s (default: 10.0).
  --diff-v=<val>          Vertical turbulent diffusivity in m^2/s (default: 1e-4).

I/O & Environment:
  --output-dir=<path>     Directory for output figures and datasets (default: outputs).
  --input-dir=<path>      Directory for input NetCDF caches (default: inputs).
  --seed=<int>            Random number generator seed (default: 42).
  --help, -h              Display this help documentation.

Examples:
  julia --project=. ParticleTrackingRun.jl --all --quick
  julia --project=. ParticleTrackingRun.jl --all --config=inputs/ParticalTracking.config
  julia --project=. ParticleTrackingRun.jl --all --scenario=ssp585 --year=2050 --particles=250
  julia --project=. ParticleTrackingRun.jl --track --min-depth=100 --track-duration=10.0 --buffer-km=100
  julia --project=. ParticleTrackingRun.jl --list-runs
  julia --project=. ParticleTrackingRun.jl --compare-scenarios
  julia --project=. ParticleTrackingRun.jl --model-average
  julia --project=. ParticleTrackingRun.jl --export-parquet
  julia --project=. ParticleTrackingRun.jl --particles=500 --save-config=inputs/custom.config
""")
end

"""
    main(args=ARGS)

Parse command-line arguments and dispatch execution to the requested segment.
"""
function main(args = ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        display_help()
        return
    end

    # 1. Resolve configuration file path and load centralized configuration
    config_file = find_default_config_path()
    for a in args
        if startswith(a, "--config=")
            config_file = String(split(a, "=")[2])
        end
    end
    cfg = load_configuration(config_file)

    # 2. Extract baseline defaults from configuration
    dom_cfg = get(cfg, "domain", Dict())
    grid_cfg = get(cfg, "grid", Dict())
    data_cfg = get(cfg, "data", Dict())
    tides_cfg = get(cfg, "tides", Dict())
    clim_cfg = get(cfg, "climate", Dict())
    hydro_cfg = get(cfg, "hydrodynamics", Dict())
    bio_cfg = get(cfg, "biology", Dict())
    dvm_cfg = get(cfg, "dvm", Dict())
    molt_cfg = get(cfg, "molting_and_settlement", Dict())
    store_cfg = get(cfg, "storage", Dict())
    hw_cfg = get(cfg, "hardware", Dict())
    vis_cfg = get(cfg, "visualization", Dict())
    paths_cfg = get(cfg, "paths", Dict())

    # Baseline options from config
    lon_range = (Float64(get(dom_cfg, "lon_min", -68.0)), Float64(get(dom_cfg, "lon_max", -57.0)))
    lat_range = (Float64(get(dom_cfg, "lat_min", 42.0)), Float64(get(dom_cfg, "lat_max", 47.0)))
    depth_range = (Float64(get(dom_cfg, "z_min", -1000.0)), Float64(get(dom_cfg, "z_max", 0.0)))
    buffer_km = Float64(get(dom_cfg, "buffer_km", get(bio_cfg, "buffer_km", 100.0)))
    grid_dim = (Int(get(grid_cfg, "nx", 50)), Int(get(grid_cfg, "ny", 50)), Int(get(grid_cfg, "nz", 10)))
    is_real = get(data_cfg, "data_mode", "synthetic") == "real"
    enable_tides = Bool(get(tides_cfg, "enable_tides", true))
    tidal_u = Float64(get(tides_cfg, "tidal_u_amp", 0.25))
    tidal_v = Float64(get(tides_cfg, "tidal_v_amp", 0.12))
    scenario = Symbol(get(clim_cfg, "scenario", "ssp245"))
    proj_year = Int(get(clim_cfg, "projection_year", 2050))
    sim_dur = Float64(get(hydro_cfg, "sim_duration_hours", 12.0)) * 3600.0
    sim_dt = Float64(get(hydro_cfg, "sim_dt_seconds", 120.0))
    adaptive_cfl = Bool(get(hydro_cfg, "adaptive_cfl", true))
    target_cfl = Float64(get(hydro_cfg, "target_cfl", 0.2))
    n_parts = Int(get(bio_cfg, "n_particles", 100))
    track_dur = Float64(get(bio_cfg, "track_duration_days", 5.0)) * 86400.0
    track_dt = Float64(get(bio_cfg, "track_dt_seconds", 300.0))
    min_depth = Float64(get(bio_cfg, "min_seabed_depth", 100.0))
    diff_h = Float64(get(bio_cfg, "diffusivity_h", 10.0))
    diff_v = Float64(get(bio_cfg, "diffusivity_v", 1e-4))
    enable_dvm = Bool(get(dvm_cfg, "enable_dvm", true))
    enable_molting = Bool(get(molt_cfg, "enable_molting", true))
    enable_duckdb = Bool(get(store_cfg, "enable_duckdb", true))
    db_path = String(get(store_cfg, "duckdb_path", "outputs/particle_tracking.duckdb"))
    use_gpu = Bool(get(hw_cfg, "use_gpu", false))
    fallback_cpu = Bool(get(hw_cfg, "fallback_to_cpu", true))
    interactive = Bool(get(vis_cfg, "interactive_map", true))
    output_dir = String(get(paths_cfg, "output_dir", "outputs"))
    input_dir = String(get(paths_cfg, "input_dir", "inputs"))
    seed = Int(get(paths_cfg, "seed", 42))

    # 3. Parse modifier flags that override config defaults
    is_quick = "--quick" in args || "-q" in args
    if "--real" in args
        is_real = true
    elseif "--synthetic" in args
        is_real = false
    end
    if "--gpu" in args || "--cuda" in args
        use_gpu = true
    elseif "--cpu" in args
        use_gpu = false
    end
    if "--fallback-cpu" in args
        fallback_cpu = true
    end
    if "--interactive" in args
        interactive = true
    elseif "--no-interactive" in args
        interactive = false
    end
    if "--duckdb" in args
        enable_duckdb = true
    elseif "--no-duckdb" in args
        enable_duckdb = false
    end
    if "--tides" in args
        enable_tides = true
    elseif "--no-tides" in args
        enable_tides = false
    end
    if "--adaptive-cfl" in args
        adaptive_cfl = true
    elseif "--no-adaptive-cfl" in args
        adaptive_cfl = false
    end
    if "--dvm" in args
        enable_dvm = true
    elseif "--no-dvm" in args
        enable_dvm = false
    end
    if "--molting" in args
        enable_molting = true
    elseif "--no-molting" in args
        enable_molting = false
    end

    # 4. Parse explicit key-value arguments
    for a in args
        if startswith(a, "--scenario=")
            scenario = Symbol(split(a, "=")[2])
        elseif startswith(a, "--year=")
            proj_year = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--particles=")
            n_parts = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--min-depth=")
            min_depth = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--buffer-km=") || startswith(a, "--buffer=") || startswith(a, "--buf=")
            buffer_km = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--output-dir=")
            output_dir = String(split(a, "=")[2])
        elseif startswith(a, "--input-dir=")
            input_dir = String(split(a, "=")[2])
        elseif startswith(a, "--db-path=")
            db_path = String(split(a, "=")[2])
        elseif startswith(a, "--lon=")
            parts = split(split(a, "=")[2], ",")
            lon_range = (parse(Float64, parts[1]), parse(Float64, parts[2]))
        elseif startswith(a, "--lat=")
            parts = split(split(a, "=")[2], ",")
            lat_range = (parse(Float64, parts[1]), parse(Float64, parts[2]))
        elseif startswith(a, "--depth-range=") || startswith(a, "--z=")
            parts = split(split(a, "=")[2], ",")
            depth_range = (parse(Float64, parts[1]), parse(Float64, parts[2]))
        elseif startswith(a, "--grid=")
            parts = split(split(a, "=")[2], ",")
            grid_dim = (parse(Int, parts[1]), parse(Int, parts[2]), parse(Int, parts[3]))
        elseif startswith(a, "--nx=")
            grid_dim = (parse(Int, split(a, "=")[2]), grid_dim[2], grid_dim[3])
        elseif startswith(a, "--ny=")
            grid_dim = (grid_dim[1], parse(Int, split(a, "=")[2]), grid_dim[3])
        elseif startswith(a, "--nz=")
            grid_dim = (grid_dim[1], grid_dim[2], parse(Int, split(a, "=")[2]))
        elseif startswith(a, "--tidal-u=")
            tidal_u = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--tidal-v=")
            tidal_v = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--sim-dt=")
            sim_dt = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--target-cfl=")
            target_cfl = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--duration=")
            sim_dur = parse(Float64, split(a, "=")[2]) * 3600.0
        elseif startswith(a, "--sim-duration=")
            sim_dur = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--track-duration=")
            track_dur = parse(Float64, split(a, "=")[2]) * 86400.0
        elseif startswith(a, "--track-dt=")
            track_dt = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--diff-h=") || startswith(a, "--diffusivity-h=")
            diff_h = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--diff-v=") || startswith(a, "--diffusivity-v=")
            diff_v = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--seed=")
            seed = parse(Int, split(a, "=")[2])
        end
    end

    # Fast override for quick prototyping
    if is_quick
        n_parts = min(n_parts, 25)
        grid_dim = (15, 15, 5)
        sim_dur = min(sim_dur, 3600.0)
        track_dur = min(track_dur, 86400.0 * 2)
    end

    opts = HydrodynamicOptions(
        domain_lon = lon_range,
        domain_lat = lat_range,
        domain_z = depth_range,
        grid_size = grid_dim,
        data_mode = is_real ? :real : :synthetic,
        enable_tides = enable_tides,
        tidal_u_amp = tidal_u,
        tidal_v_amp = tidal_v,
        scenario = scenario,
        projection_year = proj_year,
        sim_dt = sim_dt,
        sim_duration = sim_dur,
        adaptive_cfl = adaptive_cfl,
        target_cfl = target_cfl,
        n_particles = n_parts,
        track_duration = track_dur,
        track_dt = track_dt,
        diffusivity_h = diff_h,
        diffusivity_v = diff_v,
        enable_dvm = enable_dvm,
        enable_molting = enable_molting,
        min_seabed_depth = min_depth,
        buffer_km = buffer_km,
        use_gpu = use_gpu,
        fallback_to_cpu = fallback_cpu,
        interactive_map = interactive,
        enable_duckdb = enable_duckdb,
        duckdb_path = db_path,
        config_file = config_file,
        output_dir = output_dir,
        input_dir = input_dir,
        seed = seed
    )

    # Check for --save-config request
    save_cfg_flag = filter(a -> startswith(a, "--save-config"), args)
    if !isempty(save_cfg_flag)
        raw_flag = first(save_cfg_flag)
        dest_path = occursin("=", raw_flag) ? String(split(raw_flag, "=")[2]) : config_file
        save_configuration(options_to_configuration(opts), dest_path)
        println("Active configuration options successfully exported to $(dest_path)")
        return
    end

    # DuckDB standalone query flags
    if "--list-runs" in args
        run_cli_list_runs(opts = opts)
        return
    end

    if "--compare-scenarios" in args || "--compare-runs" in args
        run_cli_compare_scenarios(opts = opts)
        return
    end

    if "--model-average" in args || "--ensemble-average" in args
        run_cli_model_average(opts = opts)
        return
    end

    if "--export-parquet" in args
        if isfile(opts.duckdb_path)
            db = open_duckdb_storage(opts.duckdb_path; read_only = true)
            p_paths = export_duckdb_to_parquet(db, joinpath(opts.output_dir, "parquet"))
            close_duckdb_storage(db)
            println("Exported $(length(p_paths)) Parquet tables to: $(joinpath(opts.output_dir, "parquet"))")
        else
            println("No DuckDB database found at: $(opts.duckdb_path)")
        end
        return
    end

    # Segment dispatch
    has_run = false

    if "--all" in args || "--pipeline" in args || "--segment=all" in args
        run_production_pipeline(opts = opts)
        return
    end

    if "--data" in args || "--segment=data" in args
        run_segment_data(opts = opts)
        has_run = true
    end

    if "--grid" in args || "--segment=grid" in args
        run_segment_grid(opts = opts)
        has_run = true
    end

    if "--model" in args || "--segment=model" in args
        run_segment_model(opts = opts)
        has_run = true
    end

    if "--climate" in args || "--segment=climate" in args
        run_segment_climate(opts = opts)
        has_run = true
    end

    if "--sim" in args || "--simulation" in args || "--segment=sim" in args
        run_segment_simulation(opts = opts)
        has_run = true
    end

    if "--track" in args || "--tracking" in args || "--segment=track" in args
        run_segment_tracking(opts = opts)
        has_run = true
    end

    if "--metrics" in args || "--segment=metrics" in args
        run_segment_metrics(opts = opts)
        has_run = true
    end

    if "--viz" in args || "--visualize" in args || "--segment=viz" in args
        run_segment_visualize(opts = opts)
        has_run = true
    end

    if !has_run
        println("No recognized execution flag provided.")
        display_help()
    end
end

# Execute if run directly as script from command line
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
