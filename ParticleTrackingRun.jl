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

# Load the ParticleTracking module
import Pkg
Pkg.activate(@__DIR__, io = devnull)

# ─────────────────────────────────────────────────────────────────────────────
# Fast CLI Help Handler (Dispatched before loading Oceananigans / CairoMakie)
# ─────────────────────────────────────────────────────────────────────────────

"""
    display_help()

Print command-line usage instructions and option flags for `ParticleTrackingRun.jl`.
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

Decoupled Hydrodynamics & Multi-Cohort Tracking:
  --hydro-model=<path>    Target hydrodynamic JLD2 model file (output or input).
                          Default: outputs/hydrodynamics_<scenario>_<year>.jld2.
  --hydro-only            Run hydrodynamics only (Segments 1-5) and save to --hydro-model.
  --track-only            Run larval tracking only (Segments 6-8) using --hydro-model.
  --reuse-hydro           Reuse existing --hydro-model checkpoint if present; else simulate.
  --run-id=<string>       Unique cohort run identifier for DuckDB persistence and figures.

Individual Segment Flags:
  --data                  Run environmental data ingestion (Segment 1).
  --grid                  Run grid and immersed boundary construction (Segment 2).
  --model                 Run hydrodynamic model setup with tidal forcing (Segment 3).
  --climate               Run CMIP6 climate scenario & larval thermal biology (Segment 4).
  --sim, --simulation     Run Oceananigans hydrodynamic time integration (Segment 5).
  --track, --tracking     Run Lagrangian larval particle tracking (Segment 6).
  --metrics               Compute retention, empirical diffusion & connectivity (Segment 7).
  --viz, --visualize      Generate CairoMakie spatial figures and interactive HTML map (Segment 8).

DuckDB Analytical Storage & Model Averaging:
  --duckdb                Enable DuckDB storage archiving (default: true).
  --no-duckdb             Disable DuckDB storage archiving.
  --db-path=<path>        Custom DuckDB database path (default: outputs/particle_tracking.duckdb).
  --list-runs             Query and display all simulation runs archived in DuckDB.
  --compare-scenarios     Query and display multi-scenario comparison metrics.
  --model-average         Compute ensemble model-averaged connectivity and recruitment.
  --export-parquet        Export all DuckDB tables to Apache Parquet files.

Centralized Configuration:
  --config=<path>         Path to centralized .config file (default: inputs/ParticleTracking.config).
  --save-config[=<path>]  Export active parameter options to .config file and exit.

Ecosystem & Species Configurations:
  --snowcrab-settings     Load calibrated Snow Crab (Chionoecetes opilio) parameters:
                          500 larvae, 60-day PLD, bottom release (0.5-3.0m off bed),
                          active ascent (10 mm/s to -10m), DVM, molting (T_base = -1.5°C;
                          65, 130, 200 DD), 100x100x20 grid, Scotian Shelf/slope domain,
                          and persistence to outputs/snowcrab_tracking.duckdb.
                          Additional CLI arguments override these defaults.
  --snowcrab              Alias for --snowcrab-settings.
  --snowcrab-mode         Alias for --snowcrab-settings.
  --real-5yr              Execute 5-Year physical hydrodynamic cycle scenario.
  --climatology-2yr       Execute 2-Year climatological average cycle scenario.
  --climatology-1.5yr     Execute 1.5-Year (18-Month) climatological cycle scenario.
  --compare               Query DuckDB and display comparative scenario analytics.
  --heat-flux=<val>       Summer atmospheric heat flux in W/m² (default: 50.0).

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
  --release-mode=<mode>   Release depth mode: bottom, range, surface (default: bottom).
  --ascent                Enable post-hatch vertical ascent toward surface (default: true).
  --no-ascent             Disable initial vertical ascent.
  --ascent-speed=<val>    Vertical ascent swimming speed in m/s (default: 0.010).
  --ascent-target=<val>   Target depth in meters for ascent completion (default: -10.0).

I/O & Environment:
  --output-dir=<path>     Directory for output figures and datasets (default: outputs).
  --input-dir=<path>      Directory for input NetCDF caches (default: inputs).
  --seed=<int>            Random number generator seed (default: 42).
  --help, -h              Display this help documentation.

Examples:
  # 1. Run hydrodynamics only and save checkpoint:
  julia --project=. ParticleTrackingRun.jl --hydro-only --hydro-model=hydrodynamics1.jld2

  # 2. Track larval cohort reusing pre-computed hydrodynamics:
  julia --project=. ParticleTrackingRun.jl --track-only --hydro-model=hydrodynamics1.jld2 \\
      --run-id=cohort_spring_2020 --particles=500 --ascent

  # 3. Track second cohort with alternate vertical ascent speed:
  julia --project=. ParticleTrackingRun.jl --track-only --hydro-model=hydrodynamics1.jld2 \\
      --run-id=cohort_summer_fast --particles=500 --ascent-speed=0.015

  # 4. Fast end-to-end debug pipeline:
  julia --project=. ParticleTrackingRun.jl --all --quick
""")
end

is_cli_invocation = isempty(PROGRAM_FILE) ||
    lowercase(normpath(abspath(PROGRAM_FILE))) == lowercase(normpath(abspath(@__FILE__))) ||
    endswith(lowercase(PROGRAM_FILE), "particletrackingrun.jl")

if is_cli_invocation && (isempty(ARGS) || "--help" in ARGS || "-h" in ARGS)
    display_help()
    exit(0)
end

using
    Random,
    CairoMakie,
    NCDatasets,
    Downloads,
    DuckDB,
    DataFrames,
    DBInterface,
    Dates,
    Statistics,
    LinearAlgebra,
    TOML,
    JLD2,
    Oceananigans,
    Oceananigans.Units,
    Oceananigans.Utils,
    ParticleTracking

"""
    resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String) -> Tuple{String, String}

Resolve the target hydrodynamic JLD2 output/input path and filename. If `opts.hydro_model_file`
is specified, it is used directly (or resolved relative to `opts.output_dir` if a bare filename
is provided). If empty, returns `(joinpath(opts.output_dir, default_filename), default_filename)`.

# Inputs
- `opts::HydrodynamicOptions`: Configuration parameters.
- `default_filename::String`: Default fallback filename when none is specified.

# Outputs
- `Tuple{String, String}`: `(full_path, file_basename)`
"""
function resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String)
    if isempty(opts.hydro_model_file)
        full_path = joinpath(opts.output_dir, default_filename)
        return (full_path, default_filename)
    else
        raw = opts.hydro_model_file
        full_path = isabspath(raw) || dirname(raw) != "" ? raw : joinpath(opts.output_dir, raw)
        return (full_path, basename(raw))
    end
end

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
            @warn "Primary NOAA ERDDAP bathymetry download failed: $(err). Trying alternate mirror..."
            try
                # Secondary backup endpoint for topography
                backup_bathy_url = "https://www.ncei.noaa.gov/erddap/griddap/etopo190.nc?altitude[($(opts.domain_lat[1])):(($(opts.domain_lat[2]))][($(opts.domain_lon[1])):(($(opts.domain_lon[2]))]"
                Downloads.download(backup_bathy_url, bathy_file)
            catch backup_err
                @warn "All live bathymetry sources failed. Falling back to synthetic topography."
                generate_synthetic_bathymetry(bathy_file, lon_range = opts.domain_lon, lat_range = opts.domain_lat, n_lon = opts.grid_size[1], n_lat = opts.grid_size[2])
            end
        end

        println("Retrieving real surface winds from NOAA ERDDAP / NCEP Reanalysis...")
        try
            fetch_open_surface_winds(
                lon_range = opts.domain_lon, 
                lat_range = opts.domain_lat, 
                time_iso = "2023-06-01T00:00:00Z", 
                output_path = wind_file
            )
        catch err
            @warn "NOAA CoastWatch wind server timed out: $(err). Trying NCEP/NCAR reanalysis fallback..."
            try
                # Alternative reanalysis wind endpoint
                backup_wind_url = "https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis/surface/uwnd.10m.gauss.2023.nc"
                Downloads.download(backup_wind_url, wind_file)
            catch backup_err
                @warn "Live wind servers unreachable. Falling back to synthetic wind forcing."
                generate_synthetic_forcing(wind_file, lon_range = opts.domain_lon, lat_range = opts.domain_lat, n_lon = opts.grid_size[1], n_lat = opts.grid_size[2])
            end
        end
    else
        println("Generating synthetic Scotian Shelf bathymetry...")
        generate_synthetic_bathymetry(bathy_file, lon_range = opts.domain_lon, lat_range = opts.domain_lat, n_lon = opts.grid_size[1], n_lat = opts.grid_size[2], inshore_depth = -100.0, shelf_slope = 600.0 )
        println("Generating synthetic surface wind forcing...")
        generate_synthetic_forcing(
            wind_file,
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            time_range = (0.0, opts.sim_duration),
            n_lon = opts.grid_size[1],
            n_lat = opts.grid_size[2],
            n_time = 24,
            tau_x_amplitude = 1e-4,
            tau_y_amplitude = 2e-5
        )
    end

    bathy_info = inspect_netcdf(bathy_file, verbose = true)
    u10_ref, v10_ref = 8.5, 3.2 
    tau_x, tau_y = wind_speed_to_kinematic_stress(u10_ref, v10_ref)
    println("Calculated Large & Pond (1981) kinematic wind stress:")
    println("  Reference 10m wind: u = $(u10_ref) m/s, v = $(v10_ref) m/s")
    println("  Kinematic stress:   tau_x = $(round(tau_x, digits=6)), tau_y = $(round(tau_y, digits=6)) m^2/s^2")
    return (bathy_file = bathy_file, wind_file = wind_file, tau_x = tau_x, tau_y = tau_y, bathy_info = bathy_info)
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
        println("Configuring astronomical tidal body forcing (M2 + S2 spring-neap envelope)...")
        u_amps = Dict(:M2 => opts.tidal_u_amp, :S2 => 0.44 * opts.tidal_u_amp)
        v_amps = Dict(:M2 => opts.tidal_v_amp, :S2 => 0.42 * opts.tidal_v_amp)
        build_tidal_body_forcing(
            constituents = [:M2, :S2],
            u_amplitudes = u_amps,
            v_amplitudes = v_amps
        )
    else
        nothing
    end

    coriolis_lat = 0.5 * (opts.domain_lat[1] + opts.domain_lat[2])
    println("Building HydrostaticFreeSurfaceModel (Coriolis at $(coriolis_lat)°N, summer surface heat flux)...")
    model = build_hydrodynamic_model(
        target_grid,
        coriolis_latitude = coriolis_lat,
        surface_wind_stress_x = tau_x,
        surface_wind_stress_y = tau_y,
        surface_heat_flux = opts.surface_heat_flux,
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

    default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
    jld2_path, jld2_filename = resolve_hydro_model_path(opts, default_jld2)

    # If track-only or reuse-hydro with an existing checkpoint, skip hydrodynamics integration
    if opts.track_only || (opts.reuse_hydro && isfile(jld2_path))
        if isfile(jld2_path)
            println("Reusing existing hydrodynamic flow solution from: $(jld2_path)")
            return (simulation = nothing, jld2_output_path = jld2_path)
        else
            error("Cannot run in --track-only mode: hydrodynamic model file does not exist: $(jld2_path)\n" *
                  "Please run with --hydro-only first or specify an existing file with --hydro-model=<path>.")
        end
    end

    out_dir_target = dirname(jld2_path)
    mkpath(out_dir_target)

    println("Setting up simulation (stop_time=$(opts.sim_duration)s, Δt=$(opts.sim_dt)s)...")
    sim = setup_hydrodynamic_simulation(
        target_model,
        Δt = opts.sim_dt,
        stop_time = opts.sim_duration,
        adaptive_time_step = opts.adaptive_cfl,
        target_cfl = opts.target_cfl,
        output_dir = out_dir_target,
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
            coords = extract_grid_coordinates(target_model.grid)
            glon, glat, gdepth = coords.lons, coords.lats, coords.depths
            u_data = Array(interior(target_model.velocities.u))
            v_data = Array(interior(target_model.velocities.v))
            w_data = Array(interior(target_model.velocities.w))
            t_data = Array(interior(target_model.tracers.T))
            s_data = Array(interior(target_model.tracers.S))
            eta_data = if hasproperty(target_model, :free_surface) &&
                          hasproperty(target_model.free_surface, :η)
                Array(interior(target_model.free_surface.η))
            else
                nothing
            end

            db = open_duckdb_storage(opts.duckdb_path)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                save_hydrodynamic_field!(
                    db, target_run_id, opts;
                    grid_lons = glon, grid_lats = glat, grid_depths = gdepth,
                    u = u_data, v = v_data, w = w_data,
                    temperature = t_data, salinity = s_data,
                    elevation = eta_data,
                    time_seconds = opts.sim_duration
                )
                println("Hydrodynamic fields for '$(target_run_id)' archived in DuckDB.")
            finally
                close_duckdb_storage(db)
            end
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

    println("Initializing $(opts.n_particles) Zoea I larvae in marine water (mode = $(opts.release_depth_mode), min depth >= $(opts.min_seabed_depth) m)...")
    larvae = initialize_larval_particles(
        opts.n_particles,
        lon_range = spawn_lon,
        lat_range = spawn_lat,
        release_depth_mode = opts.release_depth_mode,
        bottom_offset = opts.bottom_release_offset,
        ascent_target_depth = opts.ascent_target_depth,
        min_seabed_depth = opts.min_seabed_depth,
        buffer_km = 0.0, # Buffer is already incorporated into spawn_lon/spawn_lat
        bathymetry = bathy_src,
        stage = :zoea1,
        rng = rng
    )

    # 4D hydrodynamic and temperature fields: use create_flow_interpolator_from_jld2
    # if simulation JLD2 archive exists; otherwise use analytical Scotian Shelf background.
    default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
    jld2_target_path, _ = resolve_hydro_model_path(opts, default_jld2)

    if opts.track_only && !isfile(jld2_target_path)
        error("Cannot run in --track-only mode: hydrodynamic model file does not exist: $(jld2_target_path)\n" *
              "Please run with --hydro-only first or specify an existing file with --hydro-model=<path>.")
    end

    sim_jld2_candidates = unique([
        jld2_target_path,
        joinpath(opts.output_dir, default_jld2),
        joinpath(opts.output_dir, "simulation_flow.jld2"),
        joinpath(opts.output_dir, "simulation_$(opts.scenario)_$(opts.projection_year).jld2"),
        joinpath("outputs", "simulation_flow.jld2")
    ])
    sim_jld2_path = findfirst(isfile, sim_jld2_candidates)
    flow_interpolator = if !isnothing(sim_jld2_path)
        actual_path = sim_jld2_candidates[sim_jld2_path]
        println("Loading 4D simulated hydrodynamic fields from $(actual_path)...")
        try
            create_flow_interpolator_from_jld2(actual_path)
        catch err
            @warn "Failed to create JLD2 flow interpolator from $(actual_path): $(err). Using analytical jet."
            nothing
        end
    else
        nothing
    end

    flow_field_fn = if !isnothing(flow_interpolator)
        (lon, lat, z, t) -> begin
            f = flow_interpolator(lon, lat, z, t)
            (f.u, f.v, f.w)
        end
    else
        # Scotian Shelf alongshore current: flows southwestward along the shelf edge.
        # u ≈ -0.08 to -0.12 m/s (westward), v ≈ -0.03 to -0.06 m/s (southward/offshore).
        # Tidal oscillation superimposed with ~6 h period at quarter-amplitude.
        # Reference: Loder, J. W. & Petrie, B. (1991), CJFAS.
        (lon, lat, z, t) -> begin
            tidal_phase = 2.0 * π * t / 44712.0  # M2 period ≈ 12.42 h
            depth_decay  = exp(z / 120.0)          # velocity decays with depth
            lon_norm = clamp((lon - opts.domain_lon[1]) /
                             (opts.domain_lon[2] - opts.domain_lon[1]), 0.0, 1.0)
            jet_factor = 0.5 + 0.8 * exp(-((lon_norm - 0.3)^2) / 0.04)
            u_mean = (-0.10 * jet_factor + 0.02 * sin(tidal_phase)) * depth_decay
            v_mean = (-0.04 * jet_factor + 0.01 * cos(tidal_phase)) * depth_decay
            w_mean = 0.00020 * sin(2.0 * π * lon_norm) * depth_decay
            (u_mean, v_mean, w_mean)
        end
    end

    temp_field_fn = if !isnothing(flow_interpolator)
        (lon, lat, z, t) -> begin
            f = flow_interpolator(lon, lat, z, t)
            f.T
        end
    else
        # 3D Scotian Shelf thermal structure: surface mixed layer + CIL + slope water.
        # Consistent with set_initial_stratification! in hydrodynamic_model.jl.
        (lon, lat, z, t) -> begin
            lon_min, lon_max = opts.domain_lon
            lat_min, lat_max = opts.domain_lat
            x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
            y_norm = clamp((lat - lat_min) / (lat_max - lat_min), 0.0, 1.0)
            T_surf   = 14.0 + 8.0 * x_norm - 5.0 * y_norm   # surface mixed layer
            T_cil_min = 1.5
            T_slope   = 8.5
            if z > -20.0
                frac = (z + 20.0) / 20.0
                clamp(T_cil_min + frac * (T_surf - T_cil_min), -1.5, 22.0)
            elseif z > -80.0
                centre_frac = (z + 50.0) / 30.0
                clamp(T_cil_min + 2.5 * centre_frac^2, -1.5, T_surf)
            else
                depth_factor = (abs(z) - 80.0) / 100.0
                clamp(T_cil_min + depth_factor * (T_slope - T_cil_min), T_cil_min, T_slope)
            end
        end
    end

    # Seabed bathymetric elevation
    bathy_field_fn = if isfile(target_bathy_path)
        get_bathymetry_interpolator(target_bathy_path)
    else
        (lon, lat) -> -120.0 - 200.0 * (lat - opts.domain_lat[1]) / (opts.domain_lat[2] - opts.domain_lat[1])
    end

    coastline_path = joinpath(opts.input_dir, "coastline.dat")
    coast_polys = isfile(coastline_path) ? load_coastline_polygons(coastline_path) : nothing

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
        tidal_constituents = [:M2, :S2],
        tidal_u_amplitudes = Dict(:M2 => opts.tidal_u_amp, :S2 => 0.44 * opts.tidal_u_amp),
        tidal_v_amplitudes = Dict(:M2 => opts.tidal_v_amp, :S2 => 0.42 * opts.tidal_v_amp),
        enable_molting = opts.enable_molting,
        enable_bbl = true,
        enable_sinking = true,
        enable_initial_ascent = opts.enable_initial_ascent,
        ascent_speed = opts.ascent_speed,
        ascent_target_depth = opts.ascent_target_depth,
        coastline = coast_polys,
        rng = rng
    )

    # Save trajectories to JLD2 checkpoint for modular reloading
    mkpath(opts.output_dir)
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    track_checkpoint = joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")
    jldsave(
        track_checkpoint;
        lons = trajectories.lons,
        lats = trajectories.lats,
        depths = trajectories.depths,
        temperatures = hasproperty(trajectories, :temperatures) ? trajectories.temperatures : nothing,
        degree_days = trajectories.degree_days,
        degree_days_timeseries = hasproperty(trajectories, :degree_days_timeseries) ? trajectories.degree_days_timeseries : nothing,
        survival_probability = hasproperty(trajectories, :survival_probability) ? trajectories.survival_probability : nothing,
        stage_survival = hasproperty(trajectories, :stage_survival) ? trajectories.stage_survival : nothing,
        stages = trajectories.stages,
        alive = trajectories.alive,
        settlement_status = trajectories.settlement_status,
        settlement_age = trajectories.settlement_age,
        ascent_duration = hasproperty(trajectories, :ascent_duration) ? trajectories.ascent_duration : nothing,
        times = trajectories.times,
        ids = trajectories.ids
    )

    # Maintain default checkpoint as fallback when custom run_id is supplied
    if !isempty(opts.run_id)
        default_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
        try
            cp(track_checkpoint, default_cp, force = true)
        catch
        end
    end

    # Save trajectories directly into DuckDB
    if opts.enable_duckdb
        try
            println("Archiving trajectories to DuckDB -> $(opts.duckdb_path)...")
            db = open_duckdb_storage(opts.duckdb_path)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                save_simulation_run!(
                    db, target_run_id, opts;
                    trajectories = trajectories,
                    config = options_to_configuration(opts),
                    notes = "Hydrodynamic tracking run ($(opts.scenario), $(opts.projection_year))" *
                            (!isempty(opts.run_id) ? " [$(opts.run_id)]" : "")
                )
                println("Trajectories for '$(target_run_id)' successfully archived in DuckDB.")
            finally
                close_duckdb_storage(db)
            end
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
        target_run_id = !isempty(opts.run_id) ? opts.run_id :
            "run_$(opts.scenario)_$(opts.projection_year)"
        tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
        track_cp = isfile(joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")) ?
            joinpath(opts.output_dir, "larval_trajectories$(tag).jld2") :
            joinpath(opts.output_dir, "larval_trajectories.jld2")

        loaded = nothing
        if opts.enable_duckdb && isfile(opts.duckdb_path)
            db = open_duckdb_storage(opts.duckdb_path; read_only = true)
            try
                println("Loading particle trajectories from DuckDB for run '$(target_run_id)'...")
                loaded = load_trajectories_namedtuple(db, target_run_id)
            catch err
                # fallback to JLD2
            finally
                close_duckdb_storage(db)
            end
        end

        if !isnothing(loaded)
            loaded
        elseif isfile(track_cp)
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
        try
            println("Archiving simulation run and metrics to DuckDB -> $(opts.duckdb_path)...")
            db = open_duckdb_storage(opts.duckdb_path)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                save_simulation_run!(
                    db,
                    target_run_id,
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
                    notes = "Hydrodynamic workflow run ($(opts.scenario), $(opts.projection_year))" *
                            (!isempty(opts.run_id) ? " [$(opts.run_id)]" : "")
                )
                println("DuckDB run '$(target_run_id)' successfully archived.")
            finally
                close_duckdb_storage(db)
            end
        catch err
            @warn "Failed to archive simulation run and metrics to DuckDB: $(err)"
        end
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
        target_run_id = !isempty(opts.run_id) ? opts.run_id :
            "run_$(opts.scenario)_$(opts.projection_year)"
        tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
        track_cp = isfile(joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")) ?
            joinpath(opts.output_dir, "larval_trajectories$(tag).jld2") :
            joinpath(opts.output_dir, "larval_trajectories.jld2")

        loaded = nothing
        if opts.enable_duckdb && isfile(opts.duckdb_path)
            db = open_duckdb_storage(opts.duckdb_path; read_only = true)
            try
                println("Loading particle trajectories from DuckDB for run '$(target_run_id)'...")
                loaded = load_trajectories_namedtuple(db, target_run_id)
            catch err
                # fallback to JLD2
            finally
                close_duckdb_storage(db)
            end
        end

        if !isnothing(loaded)
            loaded
        elseif isfile(track_cp)
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
    else
        trajectories
    end

 

    # Load bathymetry data for background contours
    target_bathy_path = joinpath(opts.input_dir, "bathymetry_active.nc")
    bathy_data = if isfile(target_bathy_path)
        load_bathymetry_from_netcdf(target_bathy_path)
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
        try
            db_h = open_duckdb_storage(opts.duckdb_path; read_only = true)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                load_hydrodynamic_field(db_h, target_run_id)
            finally
                close_duckdb_storage(db_h)
            end
        catch
            nothing
        end
    else
        nothing
    end
    plot_hydrodynamic_advection(active_hydro, bathymetry_data = bathy_data, output_path = fig8_path)

    println("Rendering hydrodynamic seawater temperature & salinity tracers -> $(fig9_path)...")
    plot_hydrodynamic_tracers(active_hydro, output_path = fig9_path)

    fig_strat_path = joinpath(opts.output_dir, "hydrodynamic_stratification.png")
    println("Rendering hydrodynamic stratification diagnostics (N², S) -> $(fig_strat_path)...")
    plot_hydrodynamic_stratification(active_hydro, output_path = fig_strat_path)

    fig_diff_path = joinpath(opts.output_dir, "hydrodynamic_diffusion.png")
    println("Rendering hydrodynamic turbulent diffusion & eddy viscosity -> $(fig_diff_path)...")
    plot_hydrodynamic_diffusion(active_hydro, output_path = fig_diff_path)

    fig_sec_path = joinpath(opts.output_dir, "hydrodynamic_section.png")
    println("Rendering hydrodynamic vertical cross-section -> $(fig_sec_path)...")
    plot_hydrodynamic_section(
        active_hydro,
        variable = :temperature,
        coordinate = 44.0,
        section_type = :lat,
        output_path = fig_sec_path
    )

    # 9. Query Archived Scenarios from DuckDB for Cross-Scenario Comparison & Multi-Layer Interactive Map
    scenarios_bundle = Dict{String, Any}()
    if opts.enable_duckdb && isfile(opts.duckdb_path)
        try
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
        catch
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

    # Decoupled hydrodynamics: skip model and simulation if --track-only is active
    model_res = nothing
    climate_res = nothing
    sim_res = nothing

    if !opts.track_only
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
    else
        println("[Pipeline Notice] --track-only active: Skipping hydrodynamic simulation (Segments 3-5).")
    end

    # If --hydro-only is active, exit early after saving hydrodynamic solution
    if opts.hydro_only
        t_elapsed = round(time() - t_start, digits = 2)
        println("\n=================================================================")
        println(" Hydrodynamics-only execution complete in $(t_elapsed) s (--hydro-only specified).")
        println(" Flow solution saved to: $(isnothing(sim_res) ? opts.hydro_model_file : sim_res.jld2_output_path)")
        println(" Skipping larval particle tracking, metrics, and visualization.")
        if opts.enable_duckdb
            close_all_duckdb_storage!()
        end
        println("=================================================================")
        return (
            data = data_res,
            grid = grid_res,
            model = model_res,
            climate = climate_res,
            simulation = sim_res,
            tracking = nothing,
            metrics = nothing,
            visualizations = nothing,
            elapsed_seconds = t_elapsed
        )
    end

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
        close_all_duckdb_storage!()
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
    main(args=ARGS)

Parse command-line arguments and dispatch execution to the requested segment.
"""
function main(args = ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        display_help()
        return
    end

    # 1. Resolve configuration file path and load centralized configuration
    is_snowcrab = "--snowcrab-settings" in args || "--snowcrab" in args || "--snowcrab-mode" in args
    config_file = is_snowcrab ? joinpath("inputs", "snowcrab.config") : find_default_config_path()
    for a in args
        if startswith(a, "--config=")
            config_file = String(split(a, "=")[2])
        end
    end
    cfg = if is_snowcrab && !isfile(config_file)
        get_snowcrab_configuration()
    else
        load_configuration(config_file)
    end

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
    surface_heat_flux = Float64(get(hydro_cfg, "surface_heat_flux", 50.0))
    hydro_model_file = String(get(hydro_cfg, "hydro_model_file", ""))
    hydro_only = Bool(get(hydro_cfg, "hydro_only", false))
    track_only = Bool(get(hydro_cfg, "track_only", false))
    reuse_hydro = Bool(get(hydro_cfg, "reuse_hydro", false))
    run_id_val = String(get(store_cfg, "run_id", ""))
    n_parts = Int(get(bio_cfg, "n_particles", 100))
    track_dur = Float64(get(bio_cfg, "track_duration_days", 5.0)) * 86400.0
    track_dt = Float64(get(bio_cfg, "track_dt_seconds", 300.0))
    min_depth = Float64(get(bio_cfg, "min_seabed_depth", 100.0))
    diff_h = Float64(get(bio_cfg, "diffusivity_h", 10.0))
    diff_v = Float64(get(bio_cfg, "diffusivity_v", 1e-4))
    rel_mode = Symbol(get(bio_cfg, "release_depth_mode", "bottom"))
    bot_off_raw = get(bio_cfg, "bottom_release_offset", [0.5, 3.0])
    bot_off = (Float64(bot_off_raw[1]), Float64(bot_off_raw[2]))
    init_ascent = Bool(get(bio_cfg, "enable_initial_ascent", true))
    asc_spd = Float64(get(bio_cfg, "ascent_speed", 0.010))
    asc_target = Float64(get(bio_cfg, "ascent_target_depth", -10.0))
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
    if "--ascent" in args
        init_ascent = true
    elseif "--no-ascent" in args
        init_ascent = false
    end
    if "--hydro-only" in args
        hydro_only = true
    end
    if "--track-only" in args
        track_only = true
    end
    if "--reuse-hydro" in args
        reuse_hydro = true
    end
    if "--real-5yr" in args
        scenario = :historical
        proj_year = 2020
        is_real = true
        sim_dur = is_quick ? 432000.0 : 157788000.0
        if isempty(run_id_val)
            run_id_val = "snowcrab_real_5yr"
        end
        if isempty(hydro_model_file)
            hydro_model_file = "hydrodynamics_real_5yr.jld2"
        end
    elseif "--climatology-2yr" in args
        scenario = :ssp245
        proj_year = 2022
        is_real = false
        sim_dur = is_quick ? 172800.0 : 63115200.0
        if isempty(run_id_val)
            run_id_val = "snowcrab_climatology_2yr"
        end
        if isempty(hydro_model_file)
            hydro_model_file = "hydrodynamics_climatology_2yr.jld2"
        end
    elseif "--climatology-1.5yr" in args || "--climatology-18mo" in args
        scenario = :ssp245
        proj_year = 2022
        is_real = false
        sim_dur = is_quick ? 172800.0 : 47336400.0
        if isempty(run_id_val)
            run_id_val = "snowcrab_climatology_1.5yr"
        end
        if isempty(hydro_model_file)
            hydro_model_file = "hydrodynamics_climatology_1.5yr.jld2"
        end
    end

    # 4. Parse explicit key-value arguments
    for a in args
        if startswith(a, "--hydro-model=")
            hydro_model_file = String(split(a, "=")[2])
        elseif startswith(a, "--run-id=")
            run_id_val = String(split(a, "=")[2])
        elseif startswith(a, "--scenario=")
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
        elseif startswith(a, "--heat-flux=")
            surface_heat_flux = parse(Float64, split(a, "=")[2])
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
        elseif startswith(a, "--release-mode=")
            rel_mode = Symbol(split(a, "=")[2])
        elseif startswith(a, "--ascent-speed=")
            asc_spd = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--ascent-target=")
            asc_target = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--seed=")
            seed = parse(Int, split(a, "=")[2])
        end
    end

    # Check for mutually exclusive flags
    if hydro_only && track_only
        error("Flags --hydro-only and --track-only are mutually exclusive. Choose one.")
    end

    # Fast override for quick prototyping
    if is_quick
        n_parts = min(n_parts, is_snowcrab ? 50 : 25)
        grid_dim = is_snowcrab ? (40, 40, 10) : (15, 15, 5)
        sim_dur = min(sim_dur, is_snowcrab ? 432000.0 : 3600.0)
        track_dur = min(track_dur, 86400.0 * 5)
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
        surface_heat_flux = surface_heat_flux,
        hydro_model_file = hydro_model_file,
        hydro_only = hydro_only,
        track_only = track_only,
        reuse_hydro = reuse_hydro,
        run_id = run_id_val,
        n_particles = n_parts,
        track_duration = track_dur,
        track_dt = track_dt,
        diffusivity_h = diff_h,
        diffusivity_v = diff_v,
        enable_dvm = enable_dvm,
        enable_molting = enable_molting,
        release_depth_mode = rel_mode,
        bottom_release_offset = bot_off,
        enable_initial_ascent = init_ascent,
        ascent_speed = asc_spd,
        ascent_target_depth = asc_target,
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

    if "--compare" in args || "--compare-scenarios" in args || "--compare-runs" in args
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

    # Standalone decoupled execution modes (when run without explicit individual segment flags)
    if opts.hydro_only && !("--sim" in args || "--simulation" in args || "--segment=sim" in args || "--all" in args)
        println("=================================================================")
        println(" Execution Mode: Hydrodynamic Simulation ONLY (--hydro-only)")
        if !isempty(opts.hydro_model_file)
            println(" Target Model File: $(opts.hydro_model_file)")
        end
        println("=================================================================")
        run_segment_data(opts = opts)
        run_segment_grid(opts = opts)
        m_res = run_segment_model(opts = opts)
        run_segment_climate(opts = opts, model = m_res.model)
        s_res = run_segment_simulation(opts = opts, model = m_res.model)
        println("\nHydrodynamic simulation completed successfully.")
        println("Output flow field saved to: $(s_res.jld2_output_path)")
        if opts.enable_duckdb
            close_all_duckdb_storage!()
        end
        return
    end

    if opts.track_only && !("--track" in args || "--tracking" in args || "--segment=track" in args || "--all" in args)
        println("=================================================================")
        println(" Execution Mode: Larval Tracking ONLY (--track-only)")
        if !isempty(opts.hydro_model_file)
            println(" Reusing Model File: $(opts.hydro_model_file)")
        end
        if !isempty(opts.run_id)
            println(" Active Run ID: $(opts.run_id)")
        end
        println("=================================================================")
        t_res = run_segment_tracking(opts = opts)
        run_segment_metrics(opts = opts, trajectories = t_res.trajectories)
        run_segment_visualize(opts = opts, trajectories = t_res.trajectories)
        println("\nLarval tracking, metrics, and visualizations completed successfully.")
        if opts.enable_duckdb
            close_all_duckdb_storage!()
        end
        return
    end

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
if is_cli_invocation
    main(ARGS)
end
