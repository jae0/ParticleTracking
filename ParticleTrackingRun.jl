# Auto-discover and activate the project environment if not already loaded
import Pkg

if Base.find_package("ParticleTracking") === nothing
    let
        curr_dir = @__DIR__
        repo_root = nothing
        for _ in 1:6
            proj = joinpath(curr_dir, "Project.toml")
            if isfile(proj)
                txt = read(proj, String)
                if occursin("name = \"ParticleTracking\"", txt)  
                    repo_root = curr_dir
                    break
                end
            end
            parent = dirname(curr_dir)
            parent == curr_dir && break
            curr_dir = parent
        end
        if repo_root !== nothing
            Pkg.activate(repo_root)
        else
            Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))
        end
    end
end


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
Pkg.activate(@__DIR__, io = devnull)

Pkg.instantiate()

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
  --data=manifest        Print the provenance registry for every physical input the run
                          consumes: product, keyless source, licence, cache path, and the
                          other acceptable sources with their access requirements. Exits
                          without running. Use this to regenerate the same inputs elsewhere.
  --hydro-only            Run hydrodynamics only (Segments 1-5) and save to --hydro-model.
  --track-only            Run larval tracking only (Segments 6-8) using --hydro-model.
  --reuse-hydro           Reuse existing --hydro-model flow file if completed; else simulate.
  --restart               Gracefully resume hydrodynamic simulation from latest checkpoint (default: true).
  --no-restart, --force-new
                          Do not resume from checkpoint; overwrite and start fresh from t=0.
  --checkpoint            Enable prognostic state checkpoints (default: true).
  --no-checkpoint         Disable prognostic state checkpoints.
  --checkpoint-interval=<val>
                          State checkpoint interval in seconds, hours (e.g. 6h), or days (e.g. 1d).
  --checkpoints-dir=<dir> Custom directory for storing restart checkpoints.
  --run-id=<string>       Unique cohort run identifier for Zarr persistence and figures.

Tidal Forcing & External Datasets:
  --tides                 Enable astronomical tidal velocity forcing (requires TPXO9 solution).
  --no-tides              Disable tidal forcing (default when solution is absent).
  --tides-source=<src>    Tidal forcing source ("tpxo9_atlas", "tpxo9", "analytic").
  --fetch-tides           Download and extract the TPXO9 global tidal solution.

Individual Segment Flags:
  --data                  Run environmental data ingestion (Segment 1).
  --grid                  Run grid and immersed boundary construction (Segment 2).
  --model                 Run hydrodynamic model setup with tidal forcing (Segment 3).
  --climate               Run CMIP6 climate scenario & larval thermal biology (Segment 4).
  --sim, --simulation     Run Oceananigans hydrodynamic time integration (Segment 5).
  --track, --tracking     Run Lagrangian larval particle tracking (Segment 6).
  --metrics               Compute retention, empirical diffusion & connectivity (Segment 7).
 --figures, --regenerate-figures
                         Redraw EVERY figure and visual (larval, biology diagnostics, Eulerian
                         hydrodynamic set, interactive HTML, animation) from stored run results.
                         Reads the trajectory checkpoint and hydrodynamics archive; does not
                         re-simulate. Use after changing a plotting function.
  --viz, --visualize      Generate CairoMakie spatial figures and interactive HTML map (Segment 8).

Zarr Analytical Storage & Model Averaging:
  --zarr                Enable Zarr storage archiving (default: true).
  --no-zarr             Disable Zarr storage archiving.
  --db-path=<path>      Custom Zarr storage path (default: outputs/particle_tracking.zarr).
  --list-runs           Query and display all simulation runs archived in Zarr.
  --compare-scenarios   Query and display multi-scenario comparison metrics.
  --model-average       Compute ensemble model-averaged connectivity and recruitment.

Centralized Configuration:
  --config=<path>         Path to centralized .toml file (default: inputs/ParticleTracking.toml).
  --save-config[=<path>]  Export active parameter options to .toml file and exit.

Ecosystem & Species Configurations:
  --snowcrab-settings     Load calibrated Snow Crab (Chionoecetes opilio) parameters:
                          500 larvae, 60-day PLD, bottom release (0.5-3.0m off bed),
                          active ascent (10 mm/s to -10m), DVM, molting (T_base = -1.5°C;
                          65, 130, 200 DD), 100x100x20 grid, Scotian Shelf/slope domain,
                          and persistence to outputs/snowcrab_tracking.zarr.
                          Additional CLI arguments override these defaults.
  --snowcrab              Alias for --snowcrab-settings.
  --snowcrab-mode         Alias for --snowcrab-settings.
  --tesselated            Snow crab configuration with depth-stratified Voronoi tessellation:
                          N=5000 units, 80% core (50-350m, 1.5km min res), 10% shallow (0-50m,
                          5km min res), 10% deep (>350m, 10km min res) and Zarr archiving.
  --real-5yr              Execute 5-Year physical hydrodynamic cycle scenario.
  --climatology-2yr       Execute 2-Year climatological average cycle scenario.
  --climatology-1.5yr     Execute 1.5-Year (18-Month) climatological cycle scenario.
  --compare               Query Zarr and display comparative scenario analytics.
  --heat-flux=<val>       Summer atmospheric heat flux in W/m² (default: 50.0).

Computational Architecture:
  --gpu, --cuda           Enable NVIDIA CUDA GPU acceleration for hydrodynamics.
  --cpu                   Execute on multi-threaded CPU (default).
  --fallback-cpu          Allow a GPU request to degrade to CPU if CUDA is not functional.
                          Off by default: without it, such a run aborts rather than quietly
                          executing somewhere other than --gpu said.

Visualization & Hydrodynamic Animation Options:
  --interactive           Export standalone interactive HTML5 Leaflet map (default).
  --no-interactive        Disable interactive HTML map generation.
  --animate-hydro, --anim-hydro
                          Render animated MP4/GIF video of hydrodynamic fields.
  --no-animate-hydro      Disable hydrodynamic animation generation (default).
  --anim-variable=<name>  Field variable: dashboard, temperature (T), advection (speed),
                          diffusion (kappa), salinity (S), viscosity (nu),
                          stratification (N2), density (rho), vorticity (zeta),
                          elevation (eta), w, richardson (Ri) (default: dashboard).
  --anim-fps=<int>        Animation playback framerate in frames/sec (default: 10).
  --anim-format=<mp4|gif> Video container format: mp4, gif (default: mp4).
  --anim-depth=<meters>   Depth slice in meters for 2D fields (default: -2.5).
  --anim-output=<path>, --anim-file=<path>
                          Designate custom output video path (e.g. outputs/flow.mp4).
  --anim-overlay-particles, --anim-particles
                          Synchronously overlay Lagrangian larvae drifting with currents.

Spatial Domain & Grid Discretization:
  --lon=<min,max>         Longitude bounding range in degrees East (default: -68.0,-57.0).
  --lat=<min,max>         Latitude bounding range in degrees North (default: 42.0,47.0).
  --depth-range=<min,max> Vertical depth range in meters (default: -1000.0,0.0).
  --grid=<nx,ny,nz>       Grid cell dimensions (default: 50,50,10; quick: 15,15,5).
  --nx=<int>              Zonal grid cells (default: 50).
  --ny=<int>              Meridional grid cells (default: 50).
  --nz=<int>              Vertical grid layers (default: 10).
  --res-scale=<val>       Resolution scaling factor (e.g. 2.5 for ~2.5km calibrated grid).
  --stretched-z           Enable hyperbolic tangent vertical layer stretching (default).
  --uniform-z             Use uniform vertical layer thicknesses.
  --z-file=<path>         Path to vertical layer coordinate CSV file.

Environmental Data & Forcing:
  --real                  Fetch real-world NOAA ERDDAP bathymetry and winds.
  --synthetic             Generate idealized synthetic shelf data (default).
  --tides                 Enable astronomical tidal body forcing (default: true).
  --no-tides              Disable tidal body forcing.
  --tidal-u=<val>         Semi-major tidal current amplitude in m/s (default: 0.25).
  --tidal-v=<val>         Semi-minor tidal current amplitude in m/s (default: 0.12).
  --wind-source=<name>    Atmospheric forcing source, matching [atmosphere] source in the TOML:
                          era5, era5_climatology, mhw. An empty value, or false/none/off,
                          disables atmospheric surface forcing.
  --obc                   Enable NumericalEarth GLORYS12V1 open boundary conditions.
  --no-obc                Disable open boundary conditions.

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
    TaylorSeries,
    CUDA,
    Oceananigans,
    Oceananigans.Units,
    Oceananigans.Utils,
    ParticleTracking,
NumericalEarth
using CSV, DataFrames, Interpolations

# NOTE: `stretched_tanh_z_faces` used to be redefined here, which shadowed the library's
# version (`ParticleTracking.stretched_tanh_z_faces`, src/data/vertical_grid.jl) with a
# duplicate that parsed the vertical-grid CSV incorrectly -- it read `CSV.read(...;
# header=false)[:, 1]`, i.e. the *name* column of a `name,z_bottom,z_top` file, so
# `Vector{Float64}` threw, the error was swallowed by a `catch`, and the run continued on an
# auto-generated tanh grid instead of the requested levels. The library version delegates to
# `load_vertical_grid_csv` and hard-errors on a length mismatch, so it is the one to use.

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

    # The coastline is a fetched input like any other, and is listed in the provenance
    # registry, so it is acquired rather than seeded from a file in this repository. A stale
    # repository copy would let two scenarios mask against different land.
    coastline_file = joinpath(opts.input_dir, "coastline.dat")
    if !isfile(coastline_file)
        src = opts.coastline_source
        println("Retrieving coastline ($(src))...")
        fetch_natural_earth_coastline(
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            resolution = replace(src, "natural_earth_" => ""),
            output_path = coastline_file
        )
    end

    bathy_file = joinpath(opts.input_dir, "bathymetry_active.nc")
    wind_file  = joinpath(opts.input_dir, "wind_active.nc")

    # Data acquisition always reads observed products. There is no "synthetic mode": the product
    # is named directly by `[data] bathy_source` / `wind_source` in the TOML, and a failed
    # download aborts the run rather than substituting an analytic field, because a run must
    # never claim a dataset-forced state it did not read.
    if !isfile(bathy_file)
        println("Retrieving bathymetry ($(opts.bathy_source))...")
        # One named source, no mirror swap, no analytic substitute. `fetch_bathymetry` stops
        # with a message listing the sources this project can actually retrieve.
        fetch_bathymetry(opts.bathy_source; lon_range = opts.domain_lon,
                         lat_range = opts.domain_lat, output_path = bathy_file)
    else
        println("Using existing real bathymetry file: $bathy_file")
    end

    if !isfile(wind_file)
        println("Retrieving surface winds ($(opts.wind_source)) for $(opts.wind_time_iso)...")
        # One named source, no mirror swap, no synthetic substitute. `fetch_surface_winds`
        # stops with a message listing the sources this project can actually retrieve.
        fetch_surface_winds(opts.wind_source; lon_range = opts.domain_lon,
                            lat_range = opts.domain_lat, time_iso = opts.wind_time_iso,
                            output_path = wind_file)
    else
        println("Using existing real surface winds file: $wind_file")
    end

    if opts.enable_tides
        src_str = lowercase(strip(String(opts.tides_source)))
        if !(src_str in ("analytic", "synthetic", "none", "off", "false"))
            tides_file = replace(data_source(:tides).cache, "<output_dir>" => opts.output_dir)
            if !isfile(tides_file)
                println("Retrieving tidal solution ($(opts.tides_source))...")
                fetch_input(:tides, opts.output_dir)
            else
                println("Using existing real tidal solution file: $tides_file")
            end
        end
    end

    # The stress is read back out of the file that was actually fetched. Deriving it from
    # hard-coded reference winds would make the fetch decorative: the numbers driving the
    # model would not come from the ingested product at all.
    w = read_wind_stress(wind_file)
    tau_x, tau_y = w.tau_x, w.tau_y
    bathy_info = (file = bathy_file, exists = true)
    println("Large & Pond (1981) kinematic wind stress, read from $(basename(wind_file)):")
    println("  Domain-mean 10m wind speed: $(round(w.speed10_rms, digits=2)) m/s")
    println("  Domain-mean stress:   tau_x = $(round(tau_x, digits=8)), " *
            "tau_y = $(round(tau_y, digits=8)) m^2/s^2")
    println("  Peak stress magnitude: $(round(w.tau_max, digits=8)) m^2/s^2")
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
    target_bathy = isnothing(bathy_file) ?
                   joinpath(opts.input_dir, "bathymetry_active.nc") : bathy_file

    if !isfile(target_bathy)
        println("Bathymetry file missing. Running Segment 1 data generation...")
        data_res = run_segment_data(opts = opts)
        target_bathy = data_res.bathy_file
    end

    println("\n=================================================================")
    println(" [Segment 2/8] Spherical Grid & Immersed Boundary Construction")
    println("=================================================================")

    arch_label = opts.use_gpu ? "GPU (CUDA)" : "CPU"
    println("Building base spherical grid on $(arch_label): $(opts.grid_size) cells...")

    z_faces = if !isempty(opts.vertical_depths)
        # `[grid] depths` is written positive-downward from the surface (0, 10, 20, ...);
        # Oceananigans wants z faces negative and ascending, so negate and reverse. The config
        # layer has already checked monotonicity, the 0 m surface, and the nz+1 face count.
        f = sort(-Float64.(opts.vertical_depths))
        # Precedence: an explicit depth list DEFINES the column, so nz must agree with it
        # rather than being applied on top. Two numbers for one axis is a configuration
        # disagreement, and picking either silently is how a run ends up with a column nobody
        # asked for.
        nz_from_depths = length(f) - 1
        nz_from_depths != opts.grid_size[3] && error(
            "[grid] nz = $(opts.grid_size[3]) conflicts with [grid] depths, which define " *
            "$(length(f)) faces and therefore Nz = $(nz_from_depths). Set nz = $(nz_from_depths), " *
            "or remove the depth list to let nz alone define the column. Refusing to pick one.")
        println("Using $(length(f)) vertical faces from [grid] depths " *
                "(shallowest dz=$(round(f[end] - f[end-1], digits=1))m, " *
                "deepest dz=$(round(f[2] - f[1], digits=1))m)...")
        f
    elseif !isempty(opts.vertical_grid_file)
        # A named `vertical_grid_file` is the operator's own file. `configuration_to_options`
        # already resolved it to a copy inside this run's input directory and refused to
        # proceed if it was missing, so reaching here means the file is present and staged.
        Lz = abs(opts.domain_z[1] - opts.domain_z[2])
        println("Using vertical layer faces from $(basename(opts.vertical_grid_file)) (Lz=$(Lz)m)...")
        stretched_tanh_z_faces(opts.grid_size[3], Lz; csv_path = opts.vertical_grid_file, strict = true)
    elseif opts.vertical_stretching_mode === :csv
        error(
            "[grid] vertical_stretching_mode = \"csv\" but neither vertical_grid_file nor " *
            "depths is set. Name a file, give a depth list, or use a mode that derives faces " *
            "from the depth range."
        )
    elseif opts.vertical_stretching_mode in (:tanh, :stretched)
        # Two-segment grid: a surface-refined upper segment that actually resolves the shallow
        # active layer, plus a coarse lower segment so the deep basin is still represented.
        # A single tanh over 0 to -5000 m puts its *coarsest* cell at the surface, which left the
        # shallowest cell centre at -258 m and no cell centre at all in the 0-50 m nursery band.
        Lz = abs(opts.domain_z[1] - opts.domain_z[2])
        nz_tot = opts.grid_size[3]
        # The segment split is a second, independent way of saying how many layers there are,
        # so it has to agree with nz rather than quietly dividing it differently.
        if opts.nz_above > 0
            nz_above = opts.nz_above
            nz_below_implied = nz_tot - nz_above
            nz_below_implied < 1 && error(
                "[grid] nz_above = $nz_above leaves only $nz_below_implied layers below it " *
                "given nz = $nz_tot. Need nz_above < nz.")
        else
            nz_above = round(Int, 0.6 * nz_tot)
        end
        nz_below = nz_tot - nz_above
        nz_below < 1 && error(
            "vertical grid split leaves nz_below = $nz_below; need nz_above < nz = $nz_tot.")
        z_break = opts.vertical_break_depth < 0.0 ? opts.vertical_break_depth : -400.0
        println("Applying two-segment vertical coordinates (break=$(z_break)m, " *
                "nz_above=$nz_above, nz_below=$nz_below, Lz=$(Lz)m)...")
        two_segment_z_faces(z_break = z_break, z_min = -Lz,
                            nz_above = nz_above, nz_below = nz_below,
                            scaling_above = 2.0, scaling_below = 0.0)
    else
        nothing
    end

    base_grid = build_shelf_grid(
        architecture = opts.use_gpu ? :gpu : :cpu,
        lon_range = opts.domain_lon,
        lat_range = opts.domain_lat,
        z_range = opts.domain_z,
        z_faces = z_faces,
        grid_size = opts.grid_size,
        fallback_to_cpu = opts.fallback_to_cpu
    )

    println("Constructing immersed boundary with 2D bilinear regridding...")
    immersed_grid = build_immersed_grid_from_real_data(
        base_grid,
        target_bathy,
        min_water_depth = abs(opts.inshore_depth),
        mask_bay_of_fundy = opts.mask_bay_of_fundy
    )

    println("Grid summary:")
    println("  Longitude: $(opts.domain_lon[1])°E to $(opts.domain_lon[2])°E (Nx=$(opts.grid_size[1]))")
    println("  Latitude:  $(opts.domain_lat[1])°N to $(opts.domain_lat[2])°N (Ny=$(opts.grid_size[2]))")
    println("  Depth:     $(opts.domain_z[1]) m to $(opts.domain_z[2]) m (Nz=$(opts.grid_size[3]))")
    if !isnothing(z_faces)
        src = isempty(opts.vertical_grid_file) ? "stretched" : basename(opts.vertical_grid_file)
        println("  Vertical:  $src (surface dz=$(round(z_faces[end]-z_faces[end-1], digits=1))m, bed dz=$(round(z_faces[2]-z_faces[1], digits=1))m)")
    else
        println("  Vertical:  Oceananigans default (nz=$(opts.grid_size[3]) evenly spaced faces)")
    end

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
    tau_y::Real = 0.0,
    wind_file::AbstractString = ""
)
    target_grid = if isnothing(immersed_grid)
        grid_res = run_segment_grid(opts = opts)
        grid_res.immersed_grid
    else
        immersed_grid
    end

    println("\n=================================================================")
    println(" [Segment 3/8] Hydrodynamic Model & Tidal Forcing Setup")
    println("=================================================================")

    # Astronomical tides, read from the global TPXO9 solution rather than imposed as a uniform
    # oscillation. The previous form applied `tidal_u_amp` (0.25 m/s) uniformly in x, y AND z;
    # on this shelf that configuration diverged to 8-13 m/s within six simulated hours, and
    # disabling tides returned the same setup to 0.4 m/s. The amplitudes in [tides] are
    # therefore no longer used to scale a body force: they describe a uniform tide that is not
    # what the shelf actually experiences, and the real field varies from 0.006 m/s on the
    # southern boundary to 0.713 m/s in the Bay of Fundy.
    #
    # The atlas is read from this run's containerised inputs directory; see `fetch_input` and
    # the :tides entry in DATA_SOURCES for how to obtain it and what the alternatives are.
    tidal_forcing = if opts.enable_tides
        src_str = lowercase(strip(String(opts.tides_source)))
        if !(src_str in ("analytic", "synthetic", "none", "off", "false"))
            tides_file = replace(data_source(:tides).cache, "<output_dir>" => opts.output_dir)
            if !isfile(tides_file)
                println("Tidal solution missing at $(tides_file); retrieving via fetch_input(:tides)...")
                fetch_input(:tides, opts.output_dir)
            end
            if !isfile(tides_file)
                error(
                    "[hydrodynamics] enable_tides is true and tides_source is '$(opts.tides_source)', " *
                    "but the tidal solution could not be retrieved:\n" *
                    "  $(tides_file)\n" *
                    "Obtain it with `fetch_input(:tides, \"$(opts.output_dir)\")` or `--fetch-tides`, " *
                    "or set enable_tides = false (or pass `--no-tides`). There is deliberately no fallback to a " *
                    "uniform body-force approximation: that was what made this configuration diverge, " *
                    "and a run that silently used a different tidal solution than its provenance " *
                    "file records is the failure this project is trying to prevent.")
            end
            # Honour `[tides] constituents`. An empty list in a config that disables tides is not an
            # error; it simply means no constituents were requested, and the branch that reads
            # harmonics is not reached.
            cons = opts.tidal_constituents
            println("Reading tidal harmonics from $(basename(tides_file)) " *
                    "(constituents: $(join(cons, ", ")))...")
            h = read_tidal_velocity_harmonics(
                tides_file;
                constituents = cons,
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
            )
            # The forcing is reduced to a small `isbits` coefficient grid BEFORE it reaches the
            # model. A closure that searches the atlas arrays at call time is not GPU-safe and
            # fails to compile in gpu_compute_hydrostatic_free_surface_Gu!.
            # Grid size is bounded by GPU kernel *parameter* memory: the coefficient tuple is
            # passed by value into every invocation, so 16x12x2x4 doubles (12 KB) is rejected
            # with "Kernel invocation uses too much parameter memory". 6x4 gives 192 doubles
            # (~1.5 KB), which fits comfortably. The tidal field varies on scales of degrees and
            # the domain is ~18 deg across, so 3 deg spacing still resolves the 0.007 m/s southern
            # boundary from the 0.45 m/s Fundy -- the contrast that caused the original blowup.
            coeff = TidalCoefficientGrid(h, (opts.domain_lon[1], opts.domain_lon[2]),
                                         (opts.domain_lat[1], opts.domain_lat[2]);
                                         nx = 6, ny = 4)
            f = tidal_velocity_coefficients(coeff)
            edge = boundary_tidal_forcing(h; edges = (:east, :south, :west), sample_stride = 4)
            amps = vcat(edge[:east].uamp[1, :], edge[:south].uamp[1, :], edge[:west].uamp[1, :])
            println("  M2 boundary speed amplitude: ", round(minimum(amps), digits = 4), " - ",
                    round(maximum(amps), digits = 4), " m/s")
            f
        else
            nothing
        end
    else
        nothing
    end

    coriolis_lat = 0.5 * (opts.domain_lat[1] + opts.domain_lat[2])

"""
    fetch_boundary_tracers(source::Symbol, opts, cfg) -> NamedTuple

Thin wrapper binding a configured source to the bounding box and month from the TOML.

Everything else -- acquisition, caching, which file holds which variable, what units the
temperature is in, and whether the file is actually usable data -- lives in
`fetch_boundary_hydrography`. The driver deliberately does not know how any source is
packaged, because that knowledge is what made two WOA23 readers, a truncated boundary file
and a units error possible in the first place.

The bounding box is the `[boundaries] embedding_*` region rather than the study domain. It
has to be larger: the model needs the water just *outside* its own edges, and a subset cut
to the study domain would put the interpolation's edge exactly on the sponge.
"""
function fetch_boundary_tracers(source::Symbol, opts)
    return fetch_boundary_hydrography(source;
        lon_range = opts.embedding_lon, lat_range = opts.embedding_lat,
        input_dir = joinpath(pwd(), "inputs"),
        month = opts.hydrography_month)
end

# Configure open boundary conditions (lateral sponge) if requested.
    #
    # Two independent things are decided here and they were previously conflated into one
    # boolean:
    #
    #   * whether the sponge band is active at all, and
    #   * whether observed boundary temperature and salinity are available to relax toward.
    #
    # The second is what the source name actually selects. Previously *any* non-GLORYS value
    # fell through to `nothing`, which meant a keyless config asking for observed boundary
    # hydrography silently got an unconstrained edge -- the request was accepted and ignored.
    # That is the worst possible failure for this input, because an unconstrained edge still
    # produces plausible-looking output.
    obc_src = opts.ocean_boundary_source
    boundary_tracer_data = nothing
    if obc_src in (:woa23, :glorys12v1, :glorys_climatology, :hycom)
        boundary_tracer_data = fetch_boundary_tracers(obc_src, opts)
        println("Boundary hydrography: observed temperature and salinity from $(obc_src).")
        println("  This constrains what the incoming water IS, not how fast it is GOING -- " *
                "the inflow speed still comes from [boundaries] u_inflow/v_inflow.")
    elseif obc_src == :synthetic
        println("Boundary hydrography: none. The open edges are unconstrained.")
    else
        error(
            "Unknown [boundaries] ocean_boundary_source \"$(obc_src)\". Choose \"synthetic\" " *
            "(unconstrained edges), \"woa23\" (keyless climatology, no currents), " *
            "\"glorys12v1\" (reanalysis with currents; needs a Copernicus Marine account), " *
            "or \"hycom\" (keyless with currents, but its server currently refuses data " *
            "requests). Nothing was substituted."
        )
    end

    obc_craft = if obc_src in (:glorys12v1, :glorys_climatology, :woa23, :hycom)
        println("Enabling lateral boundary relaxation...")
        true
    else
        nothing
    end

    # Atmospheric surface forcing, computed from the ingested ERA5 surface state rather than
    # carried as a constant. `surface_heat_flux` in the config is what decides whether the
    # flux is applied; the stress is always applied, since without it there is no momentum
    # input at all.
    #
    # `wind_file` is a parameter of this function now, not a free variable left over from the
    # segment-1 scope. It was referenced here but never in scope, so `--all` died with
    # "UndefVarError: wind_file not defined in Main" before ever reaching the flux.
    wind_path = isempty(wind_file) ? joinpath(opts.input_dir, "wind_active.nc") : wind_file
    atmo_craft = if opts.atmospheric_source in (:era5, :era5_climatology, :mhw)
        if !isfile(wind_path)
            error(
                "Atmospheric forcing is requested but the surface-wind file $(wind_path) is " *
                "missing. Run the data-ingestion segment, or disable atmospheric forcing.")
        end
        sx, sy, hx = build_bulk_surface_flux(wind_path)
        if isnothing(sx)
            error(
                "Atmospheric forcing is requested but $(basename(wind_path)) has no bulk " *
                "flux inputs. Re-fetch the wind file so the surface state is included.")
        end
        # Which heat flux to apply is decided by `[hydrodynamics] bulk_heat_flux`:
        #   true  -> the bulk flux `hx` computed from the ERA5 surface state
        #   false -> the constant `surface_heat_flux` (W m-2)
        #
        # This was previously decided by `surface_heat_flux != 0.0`, which was wrong twice
        # over. `surface_heat_flux` is a VALUE in W m-2, so treating it as a switch conflates
        # "no flux" with "a small real flux"; and the resulting `want_flux` was then used only
        # to pick a printed message, while `heat_flux = hx` was returned unconditionally. So
        # setting `surface_heat_flux = 0.0` did not disable the bulk flux -- it still drove the
        # run -- and `bulk_heat_flux`, which is parsed and stored in `HydrodynamicOptions` and
        # documented in the TOML, was never read at all. Both keys now do what they say.
        want_bulk = opts.bulk_heat_flux
        if want_bulk && isnothing(hx)
            error(
                "[hydrodynamics] bulk_heat_flux = true but the ingested surface state in " *
                "$(basename(wind_path)) is unavailable, so the bulk flux cannot be computed. " *
                "Re-fetch the wind file, or set bulk_heat_flux = false to apply the constant " *
                "`surface_heat_flux` instead."
            )
        end
        applied_flux = want_bulk ? hx : opts.surface_heat_flux
        if want_bulk
            println("Atmospheric forcing: time-varying stress from ERA5 surface state, " *
                    "with a bulk net heat flux.")
        elseif opts.surface_heat_flux != 0.0
            println("Atmospheric forcing: time-varying stress from ERA5; constant heat flux " *
                    "of $(opts.surface_heat_flux) W m-2 (bulk_heat_flux = false).")
        else
            println("Atmospheric forcing: time-varying stress from ERA5; heat flux disabled.")
        end
        (stress_x = sx, stress_y = sy, heat_flux = applied_flux)
    else
        nothing
    end

    wind_x = isnothing(atmo_craft) ? tau_x : atmo_craft.stress_x
    wind_y = isnothing(atmo_craft) ? tau_y : atmo_craft.stress_y
    surface_heat_flux_val = if isnothing(atmo_craft) || isnothing(atmo_craft.heat_flux)
        0.0
    else
        atmo_craft.heat_flux
    end

    println("Building HydrostaticFreeSurfaceModel (Coriolis at $(coriolis_lat)°N, summer surface heat flux)...")
    free_surf = ImplicitFreeSurface(maxiter = 2000, reltol = 1e-6)
    closure_choice = hasproperty(opts, :turbulence_closure) ?
        Symbol(lowercase(string(opts.turbulence_closure))) : :smagorinsky
    # Model assembly is delegated entirely to NumericalEarth.ocean_simulation, which builds
    # a `HydrostaticFreeSurfaceModel` with a spherical Coriolis (retaining the beta term),
    # a full TEOS-10 equation of state, and the split-explicit free-surface solver.
    println("Assembling model via NumericalEarth.ocean_simulation " *
            "(spherical Coriolis, TEOS-10 EOS, split-explicit free surface)...")
    model = build_hydrodynamic_model(
        target_grid,
        Δt = opts.sim_dt,
        coriolis_latitude = coriolis_lat,
        surface_wind_stress_x = wind_x,
        surface_wind_stress_y = wind_y,
        surface_heat_flux = surface_heat_flux_val,
        tidal_forcing = tidal_forcing,
        open_boundary_conditions = obc_craft,
        lateral_boundary_relaxation = opts.boundary_method != :none,
        sponge_width = opts.sponge_layer_width,
        sponge_tau = opts.sponge_timescale,
    # The prescribed tide is a VELOCITY and is relaxed toward, not forced as an acceleration.
    # See the note in `build_hydrodynamic_model`.
    tidal_tau = opts.tidal_relaxation_timescale,
        u_inflow = opts.u_inflow,
        v_inflow = opts.v_inflow,
        boundary_tracers = boundary_tracer_data,
        closure = closure_choice in (:nemotke, :catke) ? :catke : nothing,
        ν = 1e-2,
        κ = 1e-2,
        tracers = (:T, :S),
        lon_range = opts.domain_lon,
        lat_range = opts.domain_lat,
        mask_bay_of_fundy = opts.mask_bay_of_fundy,
        active_boundaries = (:east, :south, :west)
    )

    # Initial thermal and haline stratification.
    #
    # The source of the initial water column is `hydrography_source`, read from `[data]`.
    # It is deliberately independent of `ocean_boundary_source`: the lateral boundary
    # provider and the initial state are separate choices, and conflating them is what left
    # every "real data" run silently initialising from the constants `T = 10.0, S = 35.0`.
    hydro_src = opts.hydrography_source

    if hydro_src === :woa23
        month_str = lpad(string(opts.hydrography_month), 2, '0')
        t_nc = joinpath(opts.input_dir, "woa23_temperature_$(month_str)_0.25deg.nc")
        s_nc = joinpath(opts.input_dir, "woa23_salinity_$(month_str)_0.25deg.nc")
        o_nc = joinpath(opts.input_dir, "woa23_oxygen_$(month_str)_0.25deg.nc")

        if !isfile(t_nc) || !isfile(s_nc)
            println("Fetching WOA23 climatology (keyless NOAA NCEI THREDDS)...")
            # Request the study box: that is the grid extent, so no regridding extrapolation
            # is needed at the edges. `woa23_regridded_tracers` clamps in depth regardless,
            # since the model domain routinely extends below the deepest WOA level used.
            fetch_open_woa_climatology(
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                month = opts.hydrography_month,
                output_dir = opts.input_dir
            )
        end

        if isfile(t_nc) && isfile(s_nc)
            println("Applying WOA23 climatological hydrography " *
                    "(month=$(opts.hydrography_month == 0 ? "annual" : opts.hydrography_month))...")
            woa = woa23_regridded_tracers(
                model.grid, t_nc, s_nc;
                oxygen_file = isfile(o_nc) ? o_nc : nothing
            )
            # Separate `set!` calls: Oceananigans has no NamedTuple form of `set!`, and the
            # closures are only accepted as keyword values for the tracer fields.
            set!(model, T = woa.T, S = woa.S)
            set!(model, u = 0.0, v = 0.0)
            if haskey(woa, :O_2) && hasproperty(model.tracers, :O_2)
                set!(model, O_2 = woa.O_2)
            end
        else
            error(
                "hydrography_source = \"woa23\" but WOA23 files are missing from " *
                "$(opts.input_dir) (expected $(basename(t_nc)) and $(basename(s_nc))). " *
                "Download them with `fetch_open_woa_climatology` or check network access to www.ncei.noaa.gov. " *
                "To use the analytic profile instead, set [data] hydrography_source = \"synthetic\" in your config."
            )
        end
    elseif hydro_src in (:glorys12v1, :glorys_climatology, :copernicus, :cmesms)
        # Try to fetch from Copernicus Marine with dataset fallback chain
        copernicus_file = joinpath(opts.input_dir, "copernicus_ts.nc")
        try
            println("Fetching Copernicus Marine hydrography (trying GLOPHY-2 → GLO12 → GLORYS12 → NRT)...")
            fetch_copernicus_hydrography_with_fallback(
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                start_date = "$(opts.start_year)-01-01",
                end_date = "$(opts.start_year)-12-31",
                output_path = copernicus_file,
                verbose = true
            )
            println("Applying Copernicus Marine hydrography...")
            # The copernicus file has thetao/so variables; need to open and set
            # For now, error since we don't have a reader wired up
            error(
                "Copernicus file downloaded but no reader for thetao/so variables is implemented. " *
                "Implement NetCDF reading of thetao/so on the model grid, or " *
                "set [data] hydrography_source = \"synthetic\" to use the analytic profile."
            )
        catch err
            error(
                "hydrography_source = \"$(hydro_src)\" requires Copernicus Marine download " *
                "which failed: $err. " *
                "Ensure 'copernicusmarine' is installed and credentials configured. " *
                "To use the analytic profile instead, set [data] hydrography_source = \"synthetic\"."
            )
        end
    else
        println("Applying baseline thermal stratification (T_surf=15°C, dT/dz=0.01°C/m)...")
        set_initial_stratification!(
            model,
            surface_temperature = 15.0,
            temperature_gradient = 0.01,
            salinity = 35.0
        )
    end

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
    println("  CIL Temperature Anomaly:     +$(round(deltas.ΔT_cil, digits=2)) °C")
    println("  Deep Temperature Anomaly:    +$(round(deltas.ΔT_deep, digits=2)) °C")
    println("  Surface Salinity Anomaly:    $(round(deltas.ΔS_surface, digits=2)) PSU")
    println("  Atmospheric Wind Factor:     x$(round(deltas.Δwind_factor, digits=2))")

    if !isnothing(model) && opts.scenario ∉ (:baseline, :historical, :climatology) &&
       opts.ocean_boundary_source ∉ (:glorys12v1, :glorys_climatology)
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

    # (Snow Crab Larval Thermal Ecology printout removed as it is not relevant to hydrodynamics phase)

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
    default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
    jld2_path, jld2_filename = resolve_hydro_model_path(opts, default_jld2)

    # Inspect existing target hydrodynamic file.
    #
    # A run extended into several `_partN` archives is inspected across ALL of them. Reading
    # only the base path would report the span of the first increment, so a later extension
    # would look like it is resuming from far earlier than it is -- and a run whose parts are
    # complete would be judged incomplete, restarting work that was already done.
    #
    # The resolver requires the base archive to exist, which is right when it is actually
    # reassembling parts but wrong here: a run that failed before writing any output leaves
    # the grid file and possibly an empty `parts/` directory behind, and there is then no
    # output to inspect. Asking the resolver in that state aborted the run with "Simulation
    # output file not found" instead of starting fresh, so a failed run could not be re-run
    # without manually deleting the directory. The base is therefore checked first, and a
    # missing one is reported as the single non-existent path, which `inspect_hydrodynamic_file`
    # records as `exists = false` and the fresh-start path below handles.
    parts = isfile(jld2_path) ? resolve_hydro_model_paths(jld2_path) : [jld2_path]
    infos = [inspect_hydrodynamic_file(p, expected_stop_time = opts.sim_duration) for p in parts]
    file_info = if length(infos) == 1
        infos[1]
    else
        (exists = any(i -> i.exists, infos),
         is_complete = all(i -> i.is_complete, infos),
         last_time = maximum(i -> i.last_time, infos),
         n_timesteps = sum(i -> i.n_timesteps, infos))
    end
    if length(infos) > 1
        println("Existing simulation spans $(length(infos)) part(s), " *
                "reaching $(round(file_info.last_time / 3600.0, digits=2)) h " *
                "over $(file_info.n_timesteps) snapshots.")
    end

    # If track-only: require existing valid file with snapshots
    if opts.track_only
        println("\n=================================================================")
        println(" [Segment 5/8] Hydrodynamic Simulation Time Stepping")
        println("=================================================================")
        if file_info.exists && file_info.n_timesteps > 0
            println("Using existing hydrodynamic flow solution from: $(jld2_path)")
            return (simulation = nothing, jld2_output_path = jld2_path)
        else
            error("Cannot run in --track-only mode: hydrodynamic model file does not exist or has no time snapshots: $(jld2_path)\n" *
                  "Please run with --hydro-only first or specify an existing file with --hydro-model=<path>.")
        end
    end

    # If reuse-hydro requested or animate-hydro on existing hydro-model file: reuse it
    can_reuse = (opts.reuse_hydro || (opts.animate_hydro && !isempty(opts.hydro_model_file))) &&
                file_info.exists && file_info.n_timesteps > 0
    if can_reuse
        println("\n=================================================================")
        println(" [Segment 5/8] Hydrodynamic Simulation Time Stepping")
        println("=================================================================")
        if file_info.is_complete
            println("Reusing completed hydrodynamic flow solution from: $(jld2_path)")
        else
            println("Reusing hydrodynamic flow solution from: $(jld2_path)")
            println("  Notice: $(file_info.n_timesteps) snapshot(s) present (up to $(round(file_info.last_time / 3600.0, digits=2)) h).")
        end
        return (simulation = nothing, jld2_output_path = jld2_path)
    end

    target_model = if isnothing(model)
        model_res = run_segment_model(opts = opts)
        run_segment_climate(opts = opts, model = model_res.model)
        model_res.model
    else
        model
    end

    println("\n=================================================================")
    println(" [Segment 5/8] Hydrodynamic Simulation Time Stepping")
    println("=================================================================")
    mkpath(opts.output_dir)

    # Existing output: resume, extend, or discard.
    #
    # `sim_duration` is the END TIME measured from t = 0, not the length of the work still to
    # do. That single decision makes extending a run require nothing but editing the number:
    # the inputs are already cached under this scenario's input directory and are not
    # re-downloaded, the model state is picked up from wherever it stopped, and raising
    # `sim_duration` runs the extra time. Resuming an interrupted run and extending a finished
    # one are therefore the same operation, and neither needs a separate flag.
    #
    # Two conditions must not loop, and they are different:
    #
    #   * the archive holds NO progress (stored time still zero) -- resuming cannot help, so
    #     it is discarded immediately;
    #   * the archive holds progress at T, but a previous attempt to continue from T already
    #     failed. The file is unchanged by that attempt, so resuming from it reproduces the
    #     same failure indefinitely. A small sidecar records the stored time at the moment each
    #     resume was attempted; if it still matches, that continuation has been tried and the
    #     archive is discarded.
    #
    # The sidecar needs no clearing. A successful continuation advances the archive's stored
    # time past the recorded value, so the match simply stops happening; the next extension
    # overwrites it. That is why a plain "last attempted time" is sufficient state -- the
    # archive itself is the record of whether progress was made.
    out_dir_target = dirname(jld2_path)
    attempt_path = string(jld2_path, ".attempt")

    discard_run!(reason) = begin
        println(reason)
        println("  Discarding the archive and its checkpoints; starting from t = 0.")
        for f in (jld2_path, replace(jld2_path, r"\.jld2$" => "_grid.jld2"))
            isfile(f) && rm(f; force = true)
        end
        cpdir = joinpath(out_dir_target, "checkpoints")
        isdir(cpdir) && for f in readdir(cpdir; join = true)
            endswith(f, ".jld2") && rm(f; force = true)
        end
        isfile(attempt_path) && rm(attempt_path; force = true)
        nothing
    end

    # `attempted_at` is the stored time of the last resume attempt, or NaN if none was recorded
    # or the sidecar is unreadable. An unreadable sidecar is treated as "never attempted",
    # which resumes rather than discards: losing an archive on the strength of a corrupt
    # sidecar would be the worse failure.
    attempted_at = NaN
    if isfile(attempt_path)
        attempted_at = try
            parse(Float64, strip(read(attempt_path, String)))
        catch
            NaN
        end
    end
    failed_before = !isnan(attempted_at) && attempted_at == file_info.last_time

    resume = false
    if file_info.exists && !file_info.is_complete && file_info.last_time <= 0
        discard_run!("Existing output $(basename(jld2_path)) is incomplete with no recorded progress.")
    elseif file_info.exists && failed_before && opts.auto_restart
        discard_run!(
            "A previous attempt to continue this run from t = " *
            "$(round(file_info.last_time / 3600.0, digits=2)) h did not get any further " *
            "(the archive's stored time is unchanged).")
    elseif file_info.exists && !file_info.is_complete && opts.auto_restart
        println("Found interrupted hydrodynamic simulation at: $(jld2_path)")
        println("  (progress: $(round(file_info.last_time / 3600.0, digits=2)) of " *
                "$(round(opts.sim_duration / 3600.0, digits=2)) hours)")
        file_info.last_time >= opts.sim_duration &&
            println("  Resuming to the configured end time $(round(opts.sim_duration / 3600.0, digits=2)) h " *
                    "-- raise sim_duration to extend this run further.")
        println("Resuming simulation gracefully from latest checkpoint...")
        resume = true
    elseif file_info.exists && !file_info.is_complete
        println("Found incomplete simulation at $(jld2_path), but auto-restart is disabled.")
        println("Starting fresh simulation from t = 0...")
    elseif file_info.exists && file_info.is_complete && opts.auto_restart
        if file_info.last_time >= opts.sim_duration
            println("Existing run already reached $(round(file_info.last_time / 3600.0, digits=2)) h, " *
                    "which is at or beyond the configured end time " *
                    "$(round(opts.sim_duration / 3600.0, digits=2)) h.")
            println("Raise sim_duration in the TOML to extend this run; nothing to do at the " *
                    "current length.")
        else
            println("Extending the existing run from $(round(file_info.last_time / 3600.0, digits=2)) h " *
                    "to $(round(opts.sim_duration / 3600.0, digits=2)) h.")
            resume = true
        end
    end

    # Record the attempt only when we are about to continue an existing archive, and only
    # after the decision to resume has been made: a run that starts from t = 0 leaves the
    # archive behind, and a stale sidecar must not be able to discard it on the next call.
    if resume
        mkpath(out_dir_target)
        write(attempt_path, string(file_info.last_time))
    end

    out_dir_target = dirname(jld2_path)
    mkpath(out_dir_target)

    # Checkpoint storage directory
    cp_dir_target = isempty(opts.checkpoint_dir) ?
        joinpath(out_dir_target, "checkpoints") : opts.checkpoint_dir

    # Where this invocation writes its snapshots.
    #
    # Continuing an existing run gets a NEW `_partN` file rather than reopening the old one.
    # `JLD2Writer` re-serialises `serialized/coriolis`, `buoyancy`, `closure` and the grid
    # entries every time it opens a file, so reopening an archive it previously wrote raises
    # `a group or dataset named rotation_rate is already present`; `overwrite_existing` only
    # decides whether the file is deleted first, so there is no append mode to fall back on.
    # Parts are therefore written once and never reopened, and `resolve_hydro_model_paths`
    # stitches them on read.
    #
    # The boundary coincides with the resume point by construction: the run picks up from the
    # checkpoint and writes the interval that follows it.
    ext = splitext(jld2_filename)[2]
    parts_dir = joinpath(out_dir_target, "parts")
    write_dir = out_dir_target
    write_filename = jld2_filename
    if resume
        # Parts go in a `parts/` subdirectory so a run that accumulates a dozen multi-GB
        # files does not bury the inputs, checkpoints and figures beside them. The checkpoint
        # directory deliberately stays at the top level: it is working state that gets cleaned
        # up, not part of the record.
        mkpath(parts_dir)
        write_dir = parts_dir
        stem = splitext(jld2_filename)[1]
        n = 1
        while isfile(joinpath(parts_dir, "$(stem)_part$(n)$(ext)"))
            n += 1
        end
        println("Continuing into a new part: parts/$(stem)_part$(n)$(ext)")
        write_filename = "$(stem)_part$(n)$(ext)"
    end

    # Initial time step: start at 10% of the configured sim_dt to allow the
    # CFL wizard to safely ramp up from rest without violating barotropic CFL.
    init_Δt  = Float64(opts.sim_dt) * 0.10
    max_Δt   = Float64(opts.max_dt)

    println("Setting up simulation (stop_time=$(opts.sim_duration)s, Δt₀=$(init_Δt)s, max_Δt=$(max_Δt)s, min_Δt=$(opts.min_dt_seconds)s)...")
    sim = setup_hydrodynamic_simulation(
        target_model,
        Δt = init_Δt,
        stop_time = opts.sim_duration,
        adaptive_time_step = opts.adaptive_cfl,
        target_cfl = opts.target_cfl,
        max_Δt = max_Δt,
        min_Δt = opts.min_dt_seconds,
        output_path = joinpath(write_dir, write_filename),
        output_schedule = (opts.output_schedule_seconds > 0.0) ?
            opts.output_schedule_seconds : 21600.0,
        progress_schedule = 100,
        enable_checkpoint = opts.enable_checkpoint,
        checkpoint_dir = cp_dir_target,
        checkpoint_prefix = opts.checkpoint_prefix,
        checkpoint_schedule = opts.checkpoint_schedule > 0 ? opts.checkpoint_schedule : nothing,
        cleanup_checkpoints = opts.checkpoint_cleanup,
        pickup = resume ? :auto : false
    )

    println("Integrating hydrostatic primitive equations...")
    run_hydrodynamic_simulation!(
        sim,
        pickup = resume ? :auto : false,
        verbose = true
    )

    completed_sim = sim.model.clock.time >= (sim.stop_time - 1e-3)
    if completed_sim
        println("Simulation completed successfully; fields saved to: $(jld2_path)")
    else
        println("Simulation paused/interrupted at $(prettytime(sim.model.clock.time)); state saved to checkpoints.")
    end

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

            db = open_storage(opts.zarr_path)
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
                println("Hydrodynamic fields for '$(target_run_id)' archived in Zarr.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive hydrodynamic fields to Zarr: $(err)"
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
    # Span, in seconds, of the hydrodynamic record the cohort will be advected through. The cohort
    # is tracked for this whole span and stops early only once every larva is dead or settled, so
    # the drift is never shorter than the model it is driven by.
    hydro_span = [0.0, 0.0]

    flow_interpolator = if !isnothing(sim_jld2_path)
        actual_path = sim_jld2_candidates[sim_jld2_path]
        println("Loading 4D simulated hydrodynamic fields from $(actual_path)...")
        try
            # A production run stores thousands of snapshots; materialising all of them as Float64
            # for u, v, w and T is hundreds of gigabytes and exhausts memory, which previously made
            # this silently fall through to the analytical jet below -- producing recruitment numbers
            # from a synthetic flow rather than the simulated one. So: keep the full record's time
            # span, but bound how many snapshots are held resident via `max_flow_snapshots`
            # (evenly subsampled, endpoints preserved), and report the retained span back.
            interp = create_flow_interpolator_from_jld2(actual_path;
                                                        t_window = nothing,
                                                        max_snapshots = opts.max_flow_snapshots,
                                                        time_range_out = hydro_span)
            println("  retained span: $(round(hydro_span[1], digits=1)) .. " *
                    "$(round(hydro_span[2], digits=1)) s " *
                    "($(round((hydro_span[2] - hydro_span[1]) / 86400, digits=2)) days)")
            interp
        catch err
            # Falling back to an analytical current here would silently replace the simulated
            # hydrodynamics with a synthetic jet, turning a memory limit into a wrong scientific
            # result. Make the failure explicit and actionable instead.
            if opts.allow_analytical_fallback
                @warn "Failed to build the flow interpolator from $(actual_path): " *
                      "$(sprint(showerror, err)). Falling back to the analytical jet because " *
                      "allow_analytical_fallback is set -- results will NOT reflect the simulation."
                nothing
            else
                error(
                    "Failed to build the flow interpolator from $(actual_path):\n  " *
                    sprint(showerror, err) * "\n\n" *
                    "Tracking needs the simulated fields; substituting an analytical current would " *
                    "silently change the science. To proceed deliberately with a synthetic flow, " *
                    "re-run with --allow-analytical-fallback. To reduce memory, lower " *
                    "[hydrodynamics].max_flow_snapshots (currently $(opts.max_flow_snapshots))."
                )
            end
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

    # Track over the full span of the hydrodynamic record. The cohort stops early by itself once
    # every larva is dead or settled, so this is an upper bound the run only reaches if larvae
    # persist. When there is no hydro file (analytical fallback or synthetic quick run) the
    # configured pelagic larval duration is the only sensible horizon.
    track_horizon = if !isnothing(flow_interpolator) && hydro_span[2] > hydro_span[1]
        hydro_span[2] - hydro_span[1]
    else
        opts.track_duration
    end
    println("Tracking cohort over up to $(track_horizon / 86400.0) days " *
            "(dt=$(opts.track_dt)s, stops early when all larvae are dead or settled)...")
    trajectories = track_larval_cohort(
        larvae,
        velocity_fn = flow_field_fn,
        temperature_fn = temp_field_fn,
        bathymetry_fn = bathy_field_fn,
        total_duration = track_horizon,
        dt = opts.track_dt,
        max_current_speed = opts.max_current_speed,
        κ_h = opts.diffusivity_h,
        κ_v = opts.diffusivity_v,
        is_lat_lon = true,
        enable_tides = opts.enable_tides,
        tidal_u_amp = opts.tidal_u_amp,
        tidal_v_amp = opts.tidal_v_amp,
        # The constituents and their amplitude ratios come from the TOML via `opts`. Hardcoding
        # them here would run the tides on a set the operator did not choose, which is the
        # same class of silent substitution the data path was changed to avoid.
        tidal_constituents = opts.tidal_constituents,
        tidal_u_amplitudes = Dict(:M2 => opts.tidal_u_amp,
                                  :S2 => 0.44 * opts.tidal_u_amp,
                                  :N2 => 0.19 * opts.tidal_u_amp,
                                  :K1 => 0.18 * opts.tidal_u_amp,
                                  :O1 => 0.15 * opts.tidal_u_amp),
        tidal_v_amplitudes = Dict(:M2 => opts.tidal_v_amp,
                                  :S2 => 0.42 * opts.tidal_v_amp,
                                  :N2 => 0.19 * opts.tidal_v_amp,
                                  :K1 => 0.18 * opts.tidal_v_amp,
                                  :O1 => 0.15 * opts.tidal_v_amp),
        enable_molting = opts.enable_molting,
        # Stochasticity controls. These must be passed explicitly: `track_larval_cohort` defaults
        # all three coefficients of variation to 0.0, i.e. a fully deterministic cohort. Loading
        # them into `opts` and writing them to `resolved_config.toml` is not enough on its own --
        # if they are not forwarded here the provenance file records a stochastic model that the
        # run never actually used.
        cv_molt = opts.cv_molt,
        cv_mortality = opts.cv_mortality,
        cv_settlement = opts.cv_settlement,
        settlement_stochastic = opts.settlement_stochastic,
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
        # The cohort-level stage diagnostics are saved too, so `--viz` can redraw the
        # progression figures from a checkpoint without re-running Segment 6.
        cohort_molt_fraction = hasproperty(trajectories, :cohort_molt_fraction) ? trajectories.cohort_molt_fraction : nothing,
        molt_schedule = hasproperty(trajectories, :molt_schedule) ? trajectories.molt_schedule : nothing,
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

    # Save trajectories directly into Zarr
    if opts.enable_zarr
        try
            println("Archiving trajectories to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
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
                println("Trajectories for '$(target_run_id)' successfully archived in Zarr.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive trajectories to Zarr: $(err)"
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
    target_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories

    println("\n=================================================================")
    println(" [Segment 7/8] Empirical Movement, Recruitment & Connectivity")
    println("=================================================================")

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

    # 5. Depth-stratified Voronoi Tessellation & Multi-Resolution Connectivity
    voronoi_res = if opts.enable_voronoi
        println("\nComputing depth-stratified Voronoi tessellation ($(opts.voronoi_n_units) units)...")
        v_tess = generate_depth_stratified_voronoi_units(
            :numerical_earth;
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            n_units = opts.voronoi_n_units,
            prob_core = opts.voronoi_prob_core,
            prob_shallow = opts.voronoi_prob_shallow,
            prob_deep = opts.voronoi_prob_deep,
            min_res_core_km = opts.voronoi_min_res_core_km,
            min_res_shallow_km = opts.voronoi_min_res_shallow_km,
            min_res_deep_km = opts.voronoi_min_res_deep_km,
            # The TOML states band edges as positive depths below the surface; the tessellator
            # bins on signed z, so the deeper edge gets the sign.
            core_depth_range = (-opts.voronoi_core_depth_max, -opts.voronoi_core_depth_min),
            shallow_depth_range = (-opts.voronoi_shallow_depth_max, -opts.voronoi_shallow_depth_min),
            deep_depth_range = (-opts.voronoi_deep_depth_max, -opts.voronoi_deep_depth_min),
            slope_weighting = opts.voronoi_slope_weighting,
            slope_factor = opts.voronoi_slope_factor,
            seed = opts.seed
        )
        println("Generated $(length(v_tess.units)) Voronoi units across depth strata.")

        # Export Voronoi polygons to standard GeoJSON format for GIS integration
        geojson_export_path = joinpath(opts.output_dir, "voronoi_units.geojson")
        try
            export_voronoi_geojson(v_tess, geojson_export_path)
            println("Exported Voronoi polygons to GeoJSON -> $(geojson_export_path)")
        catch geo_err
            @warn "Failed exporting Voronoi GeoJSON: $(geo_err)"
        end

        v_metrics = compute_tesselated_connectivity_matrix(target_trajs, v_tess)
        println("Voronoi Macro-Strata Transition Matrix (3x3):")
        for (idx, st) in enumerate(v_metrics.strata_names)
            row_str = join([string(round(v_metrics.macro_matrix[idx, j], digits = 3)) for j in 1:length(v_metrics.strata_names)], ", ")
            println("  $(rpad(string(st), 12)) -> [$(row_str)]")
        end
        (tessellation = v_tess, metrics = v_metrics)
    else
        nothing
    end

    # 6. Export comprehensive multi-variable NetCDF and JLD2 archives
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

# 7. Archive simulation run, trajectories, metrics & connectivity in Zarr
    if opts.enable_zarr
        try
            println("Archiving simulation run and metrics to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
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
                if opts.enable_voronoi && !isnothing(voronoi_res)
                    try
                        # Zarr doesn't support raw SQL - use the storage API instead
                        if !haskey(db["metrics"], target_run_id)
                            Zarr.create_group(db["metrics"], target_run_id)
                        end
                        mgroup = db["metrics"][target_run_id]
                        for u in voronoi_res.units
                            row_id = u.id
                            Zarr.write(mgroup["unit_$(row_id)"], Dict(
                                "run_id" => target_run_id,
                                "unit_id" => u.id,
                                "lon" => u.lon,
                                "lat" => u.lat,
                                "depth" => u.depth,
                                "stratum" => string(u.stratum),
                                "area_km2" => u.area_km2,
                                "settlement_count" => u.settlement_count,
                            ))
                        end
                        println("Archived $(length(voronoi_res.units)) Voronoi units to Zarr.")
                    catch v_err
                        @warn "Failed to archive Voronoi units to Zarr: $(v_err)"
                    end
                end
                println("Zarr run '$(target_run_id)' successfully archived.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive simulation run and metrics to Zarr: $(err)"
        end
    end
                        DuckDB.flush(appender)
                        DuckDB.close(appender)
                        println("Archived $(length(v_t.units)) Voronoi units to DuckDB table 'voronoi_units'.")
                    catch v_err
                        @warn "Failed to archive Voronoi units to DuckDB: $(v_err)"
                    end
                end
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
        voronoi = voronoi_res,
        netcdf_path = nc_export_path,
        jld2_path = jld_export_path,
        duckdb_path = opts.enable_duckdb ? opts.duckdb_path : nothing
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Segment 8: Scientific Visualizations
# ─────────────────────────────────────────────────────────────────────────────

"""
    render_all_figures(; opts) -> NamedTuple

Regenerate every figure and visual for an existing run, from stored results only.

This is the implementation behind `--figures`: particle figures, the larval-biology diagnostics,
the Eulerian hydrodynamic set, the interactive map, and the animation. It reads the trajectory
checkpoint and the hydrodynamic archive rather than re-simulating, so it is the fast path to
iterate on a plotting function.

# What it draws
- **Larval** (7): trajectories, DVM profiles, settlement density, empirical movement field,
  connectivity matrix, thermal exposure, recruitment summary, plus the three biology diagnostics
  (molt progression, degree-day growth, survival curves).
- **Eulerian** (9): advection, tracers, stratification, diffusion, section, station time series,
  time–depth diagram, near-bed temperature map, T–S diagram.
- **Interactive**: the decimated HTML dashboard.
- **Animation**: the field/dashboard animation, when `--animate-hydro` is set.

Missing inputs are reported rather than fatal, so a run with, say, no hydrodynamic archive still
produces the larval figures.
"""
function render_all_figures(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Figures] Regenerating all figures and visuals")
    println("=================================================================")
    viz = run_segment_visualize(opts = opts)
    if opts.animate_hydro
        anim_path = if !isempty(opts.hydro_model_file) && isfile(opts.hydro_model_file)
            render_hydrodynamic_animation_if_requested(opts, opts.hydro_model_file)
        else
            println("  --animate-hydro set but no hydro model file available; skipping animation.")
            nothing
        end
        return merge(viz, (animation = anim_path,))
    end
    return viz
end

"""
    _ck_or(saved::Dict, key::AbstractString, default)

Read `key` from a loaded JLD2 checkpoint dictionary, falling back to `default` when the key is
absent **or** stored as `nothing`.

The tracking checkpoint writes optional fields as `nothing` when the trajectory lacks them, and
older checkpoints predate several fields entirely, so a plain `get(saved, key, default)` is not
enough — a stored `nothing` would flow into the figures as a missing value.
"""
_ck_or(saved, key::AbstractString, default) =
    (haskey(saved, key) && !isnothing(saved[key])) ? saved[key] : default

"""
    _trajectories_from_checkpoint(saved::Dict) -> NamedTuple

Build the trajectory NamedTuple from a loaded `larval_trajectories*.jld2` dictionary.

Fields that older checkpoints lack are filled with the same placeholders used before
(`temperatures` 4 °C, `degree_days_timeseries` 40 °C·d, `survival_probability` 0.95). The
stage diagnostics `cohort_molt_fraction` and `molt_schedule` are set to `nothing` when absent
and are then rebuilt by [`_complete_stage_diagnostics`](@ref); a zero placeholder would plot
as a cohort that never moults.
"""
function _trajectories_from_checkpoint(saved)
    n_p, n_t = size(saved["lons"])
    t_end = saved["times"][end]
    return (
        lons = saved["lons"],
        lats = saved["lats"],
        depths = saved["depths"],
        temperatures = _ck_or(saved, "temperatures", fill(4.0, n_p, n_t)),
        degree_days = saved["degree_days"],
        degree_days_timeseries = _ck_or(saved, "degree_days_timeseries",
                                        fill(40.0, n_p, n_t)),
        survival_probability = _ck_or(saved, "survival_probability", fill(0.95, n_p, n_t)),
        stages = saved["stages"],
        alive = saved["alive"],
        settlement_status = saved["settlement_status"],
        settlement_age = _ck_or(saved, "settlement_age", fill(t_end, n_p)),
        stage_survival = _ck_or(saved, "stage_survival", fill(1.0, n_p)),
        ascent_duration = _ck_or(saved, "ascent_duration", fill(0.0, n_p)),
        cohort_molt_fraction = _ck_or(saved, "cohort_molt_fraction", nothing),
        molt_schedule = _ck_or(saved, "molt_schedule", nothing),
        times = saved["times"],
        ids = saved["ids"]
    )
end

"""
    _complete_stage_diagnostics(trajs::NamedTuple, opts) -> NamedTuple

Add `molt_schedule` and `cohort_molt_fraction` to a trajectory set that lacks them, such as one
read back from DuckDB, which stores per-particle rows only.

Both are deterministic functions of the run configuration and the stored degree-days, so they
are recomputed rather than approximated:

- `molt_schedule = molt_schedule(cv_molt = opts.cv_molt)`. Segment 6 calls
  `track_larval_cohort` without degree-day thresholds, so the run used the shared defaults
  (65, 130, 200 °C·d); if Segment 6 starts passing thresholds, pass them here too.
- `cohort_molt_fraction[:, j] = cohort_molt_fraction(schedule, D̄_j)` with
  ``\\bar D_j = N_p^{-1} \\sum_p D_{p,j}``, the cohort-mean degree-day at step `j`. This is the
  formula `track_larval_cohort` applies, and `degree_days_timeseries[:, j]` holds the same
  per-larva values it averages (dead larvae carry their last value).

Fields already present are left unchanged. When molting is disabled the run recorded no
schedule, so nothing is added and the figure reports it as not recorded.
"""
function _complete_stage_diagnostics(trajs::NamedTuple, opts)
    has(k) = hasproperty(trajs, k) && !isnothing(getproperty(trajs, k))
    opts.enable_molting || return trajs

    sched = has(:molt_schedule) ? trajs.molt_schedule : molt_schedule(cv_molt = opts.cv_molt)
    extra = has(:molt_schedule) ? NamedTuple() : (molt_schedule = sched,)

    if !has(:cohort_molt_fraction)
        dds = trajs.degree_days_timeseries
        n_p, n_t = size(dds)
        cm = Matrix{Float64}(undef, 3, n_t)
        for j in 1:n_t
            cm[:, j] .= cohort_molt_fraction(sched, sum(view(dds, :, j)) / n_p)
        end
        extra = merge(extra, (cohort_molt_fraction = cm,))
    end
    return isempty(extra) ? trajs : merge(trajs, extra)
end

"""
    load_run_trajectories(opts::HydrodynamicOptions) -> NamedTuple

Load the trajectory set of the configured run for Segments 7 and 8.

Sources, in order:
1. DuckDB (`opts.duckdb_path`, run id `opts.run_id` or `run_<scenario>_<year>`). If the JLD2
   checkpoint describes the same cohort (equal `ids` and `times`), the checkpoint is used
   instead, because it also holds `stage_survival`, `ascent_duration` and the true
   `settlement_age`, none of which DuckDB stores.
2. The JLD2 checkpoint `larval_trajectories[_<run_id>].jld2`.
3. Neither present: Segment 6 is run.

The result is passed through [`_complete_stage_diagnostics`](@ref), so every source yields the
same set of stage fields.
"""
function load_run_trajectories(opts::HydrodynamicOptions)
    target_run_id = !isempty(opts.run_id) ? opts.run_id :
        "run_$(opts.scenario)_$(opts.projection_year)"
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    tagged_cp = joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")
    track_cp = isfile(tagged_cp) ? tagged_cp :
               joinpath(opts.output_dir, "larval_trajectories.jld2")

    from_db = nothing
    if opts.enable_zarr && isdir(opts.zarr_path)
        db = open_storage(opts.zarr_path; read_only = true)
        try
            println("Loading particle trajectories from Zarr for run '$(target_run_id)'...")
            from_db = load_trajectories_namedtuple(db, target_run_id)
        catch err
            println("  Zarr load failed ($(sprint(showerror, err))); trying the checkpoint.")
        finally
            close_storage(db)
        end
    end

    from_cp = isfile(track_cp) ? _trajectories_from_checkpoint(load(track_cp)) : nothing

    trajs = if !isnothing(from_db) && !isnothing(from_cp) &&
               Int.(from_cp.ids) == Int.(from_db.ids) &&
               Float64.(from_cp.times) == Float64.(from_db.times)
        println("  using the matching checkpoint $(track_cp) (holds the full stage record).")
        from_cp
    elseif !isnothing(from_db)
        from_db
    elseif !isnothing(from_cp)
        println("Loading particle trajectories from checkpoint: $(track_cp)...")
        from_cp
    else
        println("Trajectories checkpoint not found. Running Segment 6...")
        run_segment_tracking(opts = opts).trajectories
    end
    return _complete_stage_diagnostics(trajs, opts)
end

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
    target_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories

    println("\n=================================================================")
    println(" [Segment 8/8] Scientific Visualizations & Spatial Figures")
    println("=================================================================")
    mkpath(opts.output_dir)

 

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

    # 7b. Larval biology diagnostics: stage progression, thermal-time growth, and survival.
    # These were previously only embedded in the interactive HTML payload and never drawn, which
    # left the quantile/CDF stage model unverifiable from the run output.
    fig_molt_path = joinpath(opts.output_dir, "molt_progression.png")
    fig_dd_path = joinpath(opts.output_dir, "degree_day_growth.png")
    fig_surv_path = joinpath(opts.output_dir, "survival_curves.png")
    println("Rendering molt progression -> $(fig_molt_path)...")
    plot_molt_progression(target_trajs, output_path = fig_molt_path)
    println("Rendering degree-day growth -> $(fig_dd_path)...")
    plot_degree_day_growth(target_trajs, output_path = fig_dd_path)
    println("Rendering survival curves -> $(fig_surv_path)...")
    plot_survival_curves(target_trajs, output_path = fig_surv_path)

    # 8. Hydrodynamic Model Eulerian Advection & Tracers Figures
    fig8_path = joinpath(opts.output_dir, "hydrodynamic_advection.png")
    fig9_path = joinpath(opts.output_dir, "hydrodynamic_tracers.png")
    println("Rendering hydrodynamic advection velocity field -> $(fig8_path)...")
    active_hydro = if opts.enable_zarr && isdir(opts.zarr_path)
        try
            db_h = open_storage(opts.zarr_path; read_only = true)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                load_hydrodynamic_field(db_h, target_run_id)
            finally
                close_storage(db_h)
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

    # 8b. Eulerian time series at a representative station. This plot function existed but was
    # never called by the driver, so no run ever produced an Eulerian time series — the standard
    # check that a simulation is evolving plausibly at a fixed point.
    fig_ts_path = joinpath(opts.output_dir, "hydrodynamic_timeseries.png")
    if !isnothing(active_hydro)
        println("Rendering hydrodynamic station time series -> $(fig_ts_path)...")
        try
            plot_hydrodynamic_timeseries(
                active_hydro,
                station = (opts.domain_lon[1] + 0.55 * (opts.domain_lon[2] - opts.domain_lon[1]),
                           opts.domain_lat[1] + 0.55 * (opts.domain_lat[2] - opts.domain_lat[1])),
                variable = :temperature,
                output_path = fig_ts_path
            )
        catch err
            @warn "Station time-series plot failed (non-fatal): $(sprint(showerror, err))"
        end
    else
        fig_ts_path = nothing
    end

    # 8c. Additional Eulerian diagnostics. The time series above is now built from the real
    # snapshots; these three were absent entirely.
    fig_hov_path = joinpath(opts.output_dir, "hydrodynamic_hovmoller.png")
    fig_tbot_path = joinpath(opts.output_dir, "bottom_temperature_map.png")
    fig_ts_diagram_path = joinpath(opts.output_dir, "ts_diagram.png")
    if !isnothing(active_hydro)
        st_lon = opts.domain_lon[1] + 0.55 * (opts.domain_lon[2] - opts.domain_lon[1])
        st_lat = opts.domain_lat[1] + 0.55 * (opts.domain_lat[2] - opts.domain_lat[1])
        for (label, path, f) in (
            ("time-depth diagram", fig_hov_path,
             () -> plot_hydrodynamic_hovmoller(active_hydro; station = (st_lon, st_lat),
                                               variable = :temperature, output_path = fig_hov_path)),
            ("near-bed temperature map", fig_tbot_path,
             () -> plot_bottom_temperature_map(active_hydro; output_path = fig_tbot_path)),
            ("temperature-salinity diagram", fig_ts_diagram_path,
             () -> plot_temperature_salinity_diagram(active_hydro; output_path = fig_ts_diagram_path))
        )
            println("Rendering $label -> $(path)...")
            try
                f()
            catch err
                # A diagnostic that cannot be built from this archive must not abort the run.
                @warn "$label failed (non-fatal): $(first(split(sprint(showerror, err), "\n")))"
            end
        end
    else
        fig_hov_path = nothing
        fig_tbot_path = nothing
        fig_ts_diagram_path = nothing
    end

    # 9. Hydrodynamic Field & Dashboard Animation (Oceananigans / CairoMakie)
    anim_file_path = nothing
    if opts.animate_hydro
        anim_fmt = lowercase(opts.anim_format) == "gif" ? "gif" : "mp4"
        hydro_source = if !isnothing(active_hydro)
            active_hydro
        else
            default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
            resolved_jld2, _ = resolve_hydro_model_path(opts, default_jld2)
            isfile(resolved_jld2) ? resolved_jld2 : nothing
        end

        anim_file_path = if !isempty(strip(opts.anim_output_path))
            opts.anim_output_path
        elseif opts.anim_variable == :dashboard
            joinpath(opts.output_dir, "hydrodynamic_dashboard_animation.$(anim_fmt)")
        else
            joinpath(opts.output_dir, "hydrodynamic_$(opts.anim_variable)_animation.$(anim_fmt)")
        end

        if opts.anim_variable == :dashboard
            println("Rendering animated hydrodynamic dashboard -> $(anim_file_path)...")
            animate_hydrodynamic_dashboard(
                hydro_source;
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        else
            println("Rendering animated hydrodynamic $(opts.anim_variable) field -> $(anim_file_path)...")
            animate_hydrodynamic_field(
                hydro_source;
                variable = opts.anim_variable,
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        end
    end

    # 9. Query Archived Scenarios from Zarr for Cross-Scenario Comparison & Multi-Layer Interactive Map
    scenarios_bundle = Dict{String, Any}()
    if opts.enable_zarr && isdir(opts.zarr_path)
        try
            db = open_storage(opts.zarr_path; read_only = true)
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
                close_storage(db)
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
            # Same biology as the primary cohort, so the cross-scenario comparison is not
            # confounded by one branch running deterministic biology and the other stochastic.
            enable_molting = opts.enable_molting,
            cv_molt = opts.cv_molt,
            cv_mortality = opts.cv_mortality,
            cv_settlement = opts.cv_settlement,
            settlement_stochastic = opts.settlement_stochastic,
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

    # 11. Depth-Stratified Voronoi Units Spatial Distribution
    fig_voronoi_path = nothing
    if opts.enable_voronoi
        fig_voronoi_path = joinpath(opts.output_dir, "voronoi_units_distribution.png")
        println("Rendering Voronoi units spatial distribution -> $(fig_voronoi_path)...")
        try
            target_bathy_path = joinpath(opts.input_dir, "bathymetry_active.nc")
            if !isfile(target_bathy_path)
                real_b = joinpath(opts.input_dir, "real_bathymetry.nc")
                target_bathy_path = isfile(real_b) ? real_b : target_bathy_path
            end
            v_tess = generate_depth_stratified_voronoi_units(
                target_bathy_path;
                n_units = opts.voronoi_n_units,
                prob_core = opts.voronoi_prob_core,
                prob_shallow = opts.voronoi_prob_shallow,
                prob_deep = opts.voronoi_prob_deep,
                min_res_core_km = opts.voronoi_min_res_core_km,
                min_res_shallow_km = opts.voronoi_min_res_shallow_km,
                min_res_deep_km = opts.voronoi_min_res_deep_km,
                # TOML band edges are positive depths; the tessellator bins on signed z.
                core_depth_range = (-opts.voronoi_core_depth_max, -opts.voronoi_core_depth_min),
                shallow_depth_range = (-opts.voronoi_shallow_depth_max, -opts.voronoi_shallow_depth_min),
                deep_depth_range = (-opts.voronoi_deep_depth_max, -opts.voronoi_deep_depth_min),
                slope_weighting = opts.voronoi_slope_weighting,
                slope_factor = opts.voronoi_slope_factor,
                seed = opts.seed
            )
            fig_v = CairoMakie.Figure(size = (1000, 750), fontsize = 13)
            ax_v = CairoMakie.Axis(
                fig_v[1, 1],
                title = "Depth-Stratified Voronoi Units (N=$(length(v_tess.units)))",
                xlabel = "Longitude (°E)",
                ylabel = "Latitude (°N)"
            )
            lons = [u.lon for u in v_tess.units]
            lats = [u.lat for u in v_tess.units]
            strata_codes = [u.stratum == :core ? 1 : (u.stratum == :shallow ? 2 : 3) for u in v_tess.units]
            CairoMakie.scatter!(ax_v, lons, lats, color = strata_codes, colormap = :viridis, markersize = 6)
            CairoMakie.save(fig_voronoi_path, fig_v)
        catch v_err
            @warn "Failed to render Voronoi unit distribution figure: $(v_err)"
        end
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
        # Larval biology diagnostics
        fig_molt_progression = fig_molt_path,
        fig_degree_day_growth = fig_dd_path,
        fig_survival = fig_surv_path,
        # Eulerian
        fig_hydro_advection = fig8_path,
        fig_hydro_tracers = fig9_path,
        fig_hydro_timeseries = fig_ts_path,
        fig_hovmoller = fig_hov_path,
        fig_bottom_temperature = fig_tbot_path,
        fig_ts_diagram = fig_ts_diagram_path,
        fig_comparison = fig10_path,
        fig_voronoi = fig_voronoi_path,
        interactive_map = opts.interactive_map ? html_path : nothing,
        animation = anim_file_path
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Animation Helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    render_hydrodynamic_animation_if_requested(
        opts::HydrodynamicOptions,
        hydro_source::AbstractString
    ) -> Union{Nothing, String}

Render 2D slice or multi-panel dashboard hydrodynamic animation if requested
by command-line options (`opts.animate_hydro = true`).

# Inputs
- `opts::HydrodynamicOptions`: Workflow options containing animation settings.
- `hydro_source::AbstractString`: Path to hydrodynamic JLD2 solution file.

# Outputs
- `Union{Nothing, String}`: Path to generated animation file, or `nothing` if not requested.
"""
function render_hydrodynamic_animation_if_requested(
    opts::HydrodynamicOptions,
    hydro_source::AbstractString
)
    if !opts.animate_hydro
        return nothing
    end
    target_bathy_path = joinpath(opts.input_dir, "bathymetry_active.nc")
    bathy_data = isfile(target_bathy_path) ?
        load_bathymetry_from_netcdf(target_bathy_path) : nothing
    anim_fmt = lowercase(opts.anim_format) == "gif" ? "gif" : "mp4"
    anim_file_path = if !isempty(strip(opts.anim_output_path))
        opts.anim_output_path
    elseif opts.anim_variable == :dashboard
        joinpath(opts.output_dir, "hydrodynamic_dashboard_animation.$(anim_fmt)")
    else
        joinpath(opts.output_dir, "hydrodynamic_$(opts.anim_variable)_animation.$(anim_fmt)")
    end

    mkpath(dirname(anim_file_path))
    if opts.anim_variable == :dashboard
        println("Rendering hydrodynamic dashboard -> $(anim_file_path)...")
        animate_hydrodynamic_dashboard(
            hydro_source;
            depth = opts.anim_depth,
            trajectories = nothing,
            bathymetry_data = bathy_data,
            output_path = anim_file_path,
            framerate = opts.anim_fps,
            show_trajectories = false,
            domain_lon = opts.domain_lon,
            domain_lat = opts.domain_lat
        )
    else
        println("Rendering hydrodynamic $(opts.anim_variable) field -> $(anim_file_path)...")
        animate_hydrodynamic_field(
            hydro_source;
            variable = opts.anim_variable,
            depth = opts.anim_depth,
            trajectories = nothing,
            bathymetry_data = bathy_data,
            output_path = anim_file_path,
            framerate = opts.anim_fps,
            show_trajectories = false,
            domain_lon = opts.domain_lon,
            domain_lat = opts.domain_lat
        )
    end
    println("Hydrodynamic animation saved to: $(anim_file_path)")
    return anim_file_path
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
        println(" Skipping larval particle tracking and metrics.")
        
        viz_res = nothing
        if opts.animate_hydro
            println(" Generating requested hydrodynamic animations...")
            hydro_source = isnothing(sim_res) ? opts.hydro_model_file : sim_res.jld2_output_path
            anim_file = render_hydrodynamic_animation_if_requested(opts, hydro_source)
            viz_res = (animation = anim_file,)
        end

        if opts.enable_zarr
            close_all_storage!()
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
            visualizations = viz_res,
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
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Archived Simulation Runs in Zarr: $(opts.zarr_path)")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = list_simulation_runs(db)
    close_storage(db)

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

Query and print comparative metrics across archived climate scenarios in Zarr.
"""
function run_cli_compare_scenarios(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Comparative Analysis from Zarr")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = compare_scenarios(db)
    close_storage(db)

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
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Ensemble Model Averaging")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    runs_df = list_simulation_runs(db)

    if nrow(runs_df) == 0
        close_storage(db)
        println("No simulation data available for model averaging.")
        return
    end

    scens = unique(runs_df.scenario)
    println("Averaging across $(length(scens)) scenario models: $(join(scens, ", "))")
    ens = compute_ensemble_model_average(db, scens)
    close_storage(db)

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

    # `--data=manifest` is a query, not a run: it prints the provenance registry that names
    # every physical input, its keyless source, and the alternatives. Handled before anything
    # is loaded so it works with no network and no config.
    if any(a -> a == "--data=manifest" || a == "--data-manifest", args)
        describe_data_sources()
        return
    end

    # 1. Resolve configuration file path and load centralized configuration
    is_snowcrab_tesselated = "--tesselated" in args || "--snowcrab-tesselated" in args ||
                             "--tessellated" in args || "--snowcrab-tessellated" in args
    is_snowcrab = is_snowcrab_tesselated || "--snowcrab-settings" in args ||
                  "--snowcrab" in args || "--snowcrab-mode" in args
    # Configuration is selected solely by --config=<path>; scenario and species settings
    # live in the TOML, not in dedicated flags.
    config_file = joinpath("configs", "default.toml")
    for a in args
        if startswith(a, "--config=")
            config_file = String(split(a, "=", limit = 2)[2])
        end
    end
    cfg = load_configuration(config_file)

    # 2. Build baseline options from the parsed configuration.
    #
    # The TOML is interpreted in exactly one place: `configuration_to_options` in
    # `src/configuration.jl`. This driver used to re-read all ~80 keys itself, which meant
    # any key added to the canonical loader was silently ignored here unless it was also
    # added to the duplicate parser. That is not hypothetical: `boundaries.method`,
    # the three sponge parameters, both inflow velocities, `hydrodynamics.max_dt_seconds`
    # and `tides.s2_u_amp`/`s2_v_amp` were all read by the loader and dropped by the driver,
    # so `configs/default.toml`'s documented `method = "none"` baseline actually ran with the
    # sponge *enabled*. Do not reintroduce per-key reads here; if a key is missing, add it to
    # `configuration_to_options` so both the library and the driver see it.
    base_opts = configuration_to_options(cfg)

    if any(a -> a == "--fetch-tides" || a == "--fetch-tide", args)
        target_dir = base_opts.output_dir
        for a in args
            if startswith(a, "--output-dir=")
                target_dir = String(split(a, "=", limit = 2)[2])
            elseif startswith(a, "-o=")
                target_dir = String(split(a, "=", limit = 2)[2])
            end
        end
        println("Fetching TPXO9 tidal solution into $(target_dir)...")
        p = fetch_input(:tides, target_dir)
        println("Tidal solution successfully verified at: $(p)")
        return
    end

    lon_range = base_opts.domain_lon
    lat_range = base_opts.domain_lat
    depth_range = base_opts.domain_z
    buffer_km = base_opts.buffer_km
    grid_dim = base_opts.grid_size
    bathy_source = base_opts.bathy_source
    coastline_source = base_opts.coastline_source
    wind_source = base_opts.wind_source
    wind_time_iso = base_opts.wind_time_iso
    tides_source = base_opts.tides_source
    tidal_constituents = base_opts.tidal_constituents
    vertical_depths = base_opts.vertical_depths
    inshore_depth = base_opts.inshore_depth
    shelf_slope = base_opts.shelf_slope
    mask_bay_of_fundy = base_opts.mask_bay_of_fundy
    enable_tides = base_opts.enable_tides
    tidal_u = Float64(base_opts.tidal_u_amp)
    tidal_v = Float64(base_opts.tidal_v_amp)
    s2_u = Float64(base_opts.s2_u_amp)
    s2_v = Float64(base_opts.s2_v_amp)
    max_dt_val = Float64(base_opts.max_dt)
    min_dt_val = Float64(base_opts.min_dt_seconds)
    scenario = base_opts.scenario
    proj_year = base_opts.projection_year
    sim_dur = base_opts.sim_duration
    sim_dt = base_opts.sim_dt
    adaptive_cfl = base_opts.adaptive_cfl
    target_cfl = base_opts.target_cfl
    surface_heat_flux = base_opts.surface_heat_flux
    hydro_model_file = base_opts.hydro_model_file
    hydro_only = base_opts.hydro_only
    track_only = base_opts.track_only
    reuse_hydro = base_opts.reuse_hydro
    allow_analytical_fallback = base_opts.allow_analytical_fallback
    run_id_val = base_opts.run_id
    n_parts = base_opts.n_particles
    track_dur = base_opts.track_duration
    track_dt = base_opts.track_dt
    min_depth = base_opts.min_seabed_depth
    diff_h = base_opts.diffusivity_h
    diff_v = base_opts.diffusivity_v
    rel_mode = base_opts.release_depth_mode
    bot_off = base_opts.bottom_release_offset
    init_ascent = base_opts.enable_initial_ascent
    asc_spd = base_opts.ascent_speed
    asc_target = base_opts.ascent_target_depth
    enable_dvm = base_opts.enable_dvm
    enable_molting = base_opts.enable_molting
    enable_duckdb = base_opts.enable_duckdb
    db_path = base_opts.duckdb_path
    use_gpu = base_opts.use_gpu
    fallback_cpu = base_opts.fallback_to_cpu
    interactive = base_opts.interactive_map
    output_dir = base_opts.output_dir
    input_dir = base_opts.input_dir
    seed = base_opts.seed

    enable_cp = base_opts.enable_checkpoint
    cp_prefix_default = "checkpoint_$(resolve_config_name(config_file))"
    cp_prefix = String(base_opts.checkpoint_prefix)
    if isempty(strip(cp_prefix)) || cp_prefix == "checkpoint"
        cp_prefix = cp_prefix_default
    end
    cp_sched = base_opts.checkpoint_schedule
    out_sched = base_opts.output_schedule_seconds
    cp_dir = String(base_opts.checkpoint_dir)
    cp_clean = base_opts.checkpoint_cleanup
    auto_res = base_opts.auto_restart

    res_scale = base_opts.resolution_scale
    v_mode = base_opts.vertical_stretching_mode
    v_file = base_opts.vertical_grid_file
    atmo_src = base_opts.atmospheric_source
    atmo_bulk_heat_flux = base_opts.bulk_heat_flux
    obc_src = base_opts.ocean_boundary_source
    obc_tp = base_opts.obc_type
    hydro_src = base_opts.hydrography_source
    hydro_month = base_opts.hydrography_month

    # Lateral boundary / sponge settings. These were previously read by neither the driver
    # nor (for the sponge parameters) honoured downstream, so `method = "none"` in
    # configs/default.toml did not actually disable the sponge.
    boundary_method = base_opts.boundary_method
    sponge_width = Float64(base_opts.sponge_layer_width)
    sponge_tau = Float64(base_opts.sponge_timescale)
    tidal_tau = Float64(base_opts.tidal_relaxation_timescale)
    u_inflow_val = Float64(base_opts.u_inflow)
    v_inflow_val = Float64(base_opts.v_inflow)

    # Voronoi tessellation defaults from config
    enable_voronoi = base_opts.enable_voronoi || is_snowcrab_tesselated
    voronoi_n_units = base_opts.voronoi_n_units
    voronoi_p_core = base_opts.voronoi_prob_core
    voronoi_p_shallow = base_opts.voronoi_prob_shallow
    voronoi_p_deep = base_opts.voronoi_prob_deep
    voronoi_min_core = base_opts.voronoi_min_res_core_km
    voronoi_min_shallow = base_opts.voronoi_min_res_shallow_km
    voronoi_min_deep = base_opts.voronoi_min_res_deep_km
    voronoi_core_min = base_opts.voronoi_core_depth_min
    voronoi_core_max = base_opts.voronoi_core_depth_max
    voronoi_shallow_min = base_opts.voronoi_shallow_depth_min
    voronoi_shallow_max = base_opts.voronoi_shallow_depth_max
    voronoi_deep_min = base_opts.voronoi_deep_depth_min
    voronoi_deep_max = base_opts.voronoi_deep_depth_max
    voronoi_slope_wt = base_opts.voronoi_slope_weighting
    voronoi_slope_fac = base_opts.voronoi_slope_factor

    # Visualization and animation defaults
    anim_hydro = base_opts.animate_hydro
    anim_var = base_opts.anim_variable
    anim_fps = base_opts.anim_fps
    anim_fmt = base_opts.anim_format
    anim_depth = base_opts.anim_depth
    anim_overlay_parts = base_opts.anim_overlay_particles
    anim_out_path = base_opts.anim_output_path

    # 3. Parse modifier flags that override config defaults
    is_quick = "--quick" in args || "-q" in args
    if "--gpu" in args || "--cuda" in args
        use_gpu = true
    elseif "--cpu" in args
        use_gpu = false
    end
    if "--fallback-cpu" in args
        fallback_cpu = true
    end
    if "--stretched-z" in args || "--vertical-tanh" in args
        v_mode = :tanh
    elseif "--uniform-z" in args
        v_mode = :uniform
    end
    if "--obc" in args
        obc_src = :glorys12v1
    elseif "--no-obc" in args
        obc_src = :none
    end
    # Atmospheric forcing source, named after `[atmosphere] source` so the flag and the TOML key
    # say the same thing. An empty value, or `false`/`none`/`off`, means no atmospheric forcing.
    # The driver decides whether to attach fluxes by testing membership of the list of sources
    # that produce them (see the `atmo_craft` branch), so "off" is a source that is not in it.
    for a in args
        if startswith(a, "--wind-source=")
            raw = lowercase(strip(String(split(a, "=", limit = 2)[2])))
            atmo_src = raw in ("", "false", "none", "off", "no") ? :none : Symbol(raw)
        end
    end
    if "--interactive" in args
        interactive = true
    elseif "--no-interactive" in args
        interactive = false
    end
    if "--animate-hydro" in args || "--anim-hydro" in args
        anim_hydro = true
    elseif "--no-animate-hydro" in args || "--no-anim-hydro" in args
        anim_hydro = false
    end
    if "--anim-overlay-particles" in args || "--anim-particles" in args
        anim_overlay_parts = true
    elseif "--no-anim-overlay-particles" in args
        anim_overlay_parts = false
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
    for a in args
        if startswith(a, "--tides-source=")
            tides_source = strip(String(split(a, "=", limit = 2)[2]))
        elseif startswith(a, "--tide-source=")
            tides_source = strip(String(split(a, "=", limit = 2)[2]))
        end
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
    if "--tesselated" in args || "--tessellated" in args || "--voronoi" in args
        enable_voronoi = true
    elseif "--no-tesselated" in args || "--no-tessellated" in args || "--no-voronoi" in args
        enable_voronoi = false
    end
    if "--hydro-only" in args
        hydro_only = true
    end
    if "--track-only" in args
        track_only = true
    end
    if "--allow-analytical-fallback" in args
        allow_analytical_fallback = true
    end
    if "--reuse-hydro" in args
        reuse_hydro = true
    end
    if "--restart" in args
        auto_res = true
    elseif "--no-restart" in args || "--force-new" in args
        auto_res = false
    end
    if "--checkpoint" in args
        enable_cp = true
    elseif "--no-checkpoint" in args
        enable_cp = false
    end
    if "--checkpoint-cleanup" in args
        cp_clean = true
    elseif "--no-checkpoint-cleanup" in args
        cp_clean = false
    end
    if "--real-5yr" in args
        scenario = :historical
        proj_year = 2020
        sim_dur = is_quick ? 432000.0 : 157788000.0
        if isempty(run_id_val)
            run_id_val = "snowcrab_real_5yr"
        end
        if isempty(hydro_model_file)
            hydro_model_file = "hydrodynamics_real_5yr.jld2"
        end
    elseif "--climatology-2yr" in args
        scenario = :climatology
        proj_year = 2022
        sim_dur = is_quick ? 172800.0 : 63115200.0
        if isempty(run_id_val)
            run_id_val = "snowcrab_climatology_2yr"
        end
        if isempty(hydro_model_file)
            hydro_model_file = "hydrodynamics_climatology_2yr.jld2"
        end
    elseif "--climatology-1.5yr" in args || "--climatology-18mo" in args
        scenario = :climatology
        proj_year = 2022
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
        elseif startswith(a, "--checkpoint-interval=") || startswith(a, "--checkpoint-schedule=")
            raw_cp = split(a, "=")[2]
            cp_sched = if endswith(raw_cp, "h")
                parse(Float64, raw_cp[1:end-1]) * 3600.0
            elseif endswith(raw_cp, "m")
                parse(Float64, raw_cp[1:end-1]) * 60.0
            elseif endswith(raw_cp, "d")
                parse(Float64, raw_cp[1:end-1]) * 86400.0
            else
                parse(Float64, raw_cp)
            end
        elseif startswith(a, "--checkpoints-dir=") || startswith(a, "--checkpoint-dir=")
            cp_dir = String(split(a, "=")[2])
        elseif startswith(a, "--checkpoint-prefix=") || startswith(a, "--cp-prefix=")
            cp_prefix = String(split(a, "=")[2])
        elseif startswith(a, "--res-scale=") || startswith(a, "--resolution-scale=")
            res_scale = parse(Float64, split(a, "=")[2])
            grid_dim = (max(10, round(Int, grid_dim[1] / res_scale)),
                        max(10, round(Int, grid_dim[2] / res_scale)),
                        grid_dim[3])
        elseif startswith(a, "--z-file=") || startswith(a, "--vertical-grid-file=")
            v_file = String(split(a, "=")[2])
        elseif startswith(a, "--obc-source=")
            obc_src = Symbol(split(a, "=")[2])
        elseif startswith(a, "--obc-type=")
            obc_tp = Symbol(split(a, "=")[2])
        elseif a == "--mask-bay-of-fundy" || a == "--mask-fundy"
            mask_bay_of_fundy = true
        elseif a == "--no-mask-bay-of-fundy" || a == "--no-mask-fundy"
            mask_bay_of_fundy = false
        elseif startswith(a, "--seed=")
            seed = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--voronoi-units=") || startswith(a, "--n-units=")
            voronoi_n_units = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--voronoi-prob-core=")
            voronoi_p_core = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--voronoi-prob-shallow=")
            voronoi_p_shallow = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--voronoi-prob-deep=")
            voronoi_p_deep = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--voronoi-min-core=")
            voronoi_min_core = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--voronoi-min-shallow=")
            voronoi_min_shallow = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--voronoi-min-deep=")
            voronoi_min_deep = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--anim-variable=") || startswith(a, "--anim-var=")
            anim_var = Symbol(lowercase(split(a, "=")[2]))
        elseif startswith(a, "--anim-fps=") || startswith(a, "--anim-framerate=")
            anim_fps = parse(Int, split(a, "=")[2])
        elseif startswith(a, "--anim-format=") || startswith(a, "--anim-fmt=")
            anim_fmt = String(lowercase(split(a, "=")[2]))
        elseif startswith(a, "--anim-depth=")
            anim_depth = parse(Float64, split(a, "=")[2])
        elseif startswith(a, "--anim-output=") || startswith(a, "--anim-file=")
            anim_out_path = String(split(a, "=")[2])
        end
    end

    # Check for mutually exclusive flags
    if hydro_only && track_only
        error("Flags --hydro-only and --track-only are mutually exclusive. Choose one.")
    end

    # Fast override for quick prototyping.
    #
    # `nz` is deliberately NOT reduced. The horizontal can be coarsened for speed without
    # changing what the column resolves, but the number of vertical layers determines the
    # top cell thickness, which sets the time step and -- far more importantly -- decides
    # whether the shallowest cell resolves the nursery band at all. Coarsening it silently
    # was how `--quick` turned a configured nz = 10 into Nz = 5 and produced
    # "Custom z_faces length 11 must equal Nz + 1 = 6" one segment later, with the real
    # disagreement buried in a grid-construction error that named neither `nz` nor `--quick`.
    if is_quick
        n_parts = min(n_parts, is_snowcrab ? 50 : 25)
        grid_dim = is_snowcrab ? (40, 40, grid_dim[3]) : (15, 15, grid_dim[3])
        sim_dur = min(sim_dur, is_snowcrab ? 432000.0 : 3600.0)
        track_dur = min(track_dur, 86400.0 * 5)
        voronoi_n_units = min(voronoi_n_units, 200)
    end

    opts = HydrodynamicOptions(
        domain_lon = lon_range,
        domain_lat = lat_range,
        domain_z = depth_range,
        grid_size = grid_dim,
        bathy_source = bathy_source,
        coastline_source = coastline_source,
        wind_source = wind_source,
        wind_time_iso = wind_time_iso,
        tides_source = tides_source,
        tidal_constituents = tidal_constituents,
        vertical_depths = vertical_depths,
        inshore_depth = inshore_depth,
        shelf_slope = shelf_slope,
        mask_bay_of_fundy = mask_bay_of_fundy,
        enable_tides = enable_tides,
        tidal_u_amp = tidal_u,
        tidal_v_amp = tidal_v,
        s2_u_amp = s2_u,
        s2_v_amp = s2_v,
        max_dt = max_dt_val,
        min_dt_seconds = min_dt_val,
        # The options struct is REBUILT here rather than mutated, so any field not listed silently
        # reverts to the constructor default and the TOML value is lost. These were omitted, which
        # quietly reset the biology to cv = 0.0 (fully deterministic) while
        # `resolved_config.toml` recorded 0.25 -- the exact "config value never reaches its
        # consumer" failure this project has hit three times now.
        cv_molt = base_opts.cv_molt,
        cv_mortality = base_opts.cv_mortality,
        cv_settlement = base_opts.cv_settlement,
        settlement_stochastic = base_opts.settlement_stochastic,
        max_current_speed = base_opts.max_current_speed,
        max_flow_snapshots = base_opts.max_flow_snapshots,
        nz_above = base_opts.nz_above,
        vertical_break_depth = base_opts.vertical_break_depth,
        # `allow_analytical_fallback`, `reuse_hydro` and `track_only` are already set further
        # down this call from their CLI locals; do not repeat them here (a repeated keyword is a
        # syntax error, and the pair of "fixes" for them was duplicated by exactly this mistake).
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
    allow_analytical_fallback = allow_analytical_fallback,
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
        enable_checkpoint = enable_cp,
        checkpoint_prefix = cp_prefix,
        checkpoint_schedule = cp_sched,
        output_schedule_seconds = out_sched,
        checkpoint_dir = cp_dir,
        checkpoint_cleanup = cp_clean,
        auto_restart = auto_res,
        seed = seed,
        vertical_stretching_mode = v_mode,
        vertical_grid_file = v_file,
        resolution_scale = res_scale,
        atmospheric_source = atmo_src,
        # `bulk_heat_flux` decides whether the heat flux is computed from the ingested
        # atmospheric surface state or taken as the constant `surface_heat_flux`. It is read
        # from the `[atmosphere]` section of the TOML, alongside `atmospheric_source`. Without
        # this the field reverted to its default `true` and the bulk flux was applied even
        # when the config set it to false.
        bulk_heat_flux = atmo_bulk_heat_flux,
        ocean_boundary_source = obc_src,
        hydrography_source = hydro_src,
        hydrography_month = hydro_month,
        obc_type = obc_tp,
        boundary_method = boundary_method,
        sponge_layer_width = sponge_width,
        sponge_timescale = sponge_tau,
    tidal_relaxation_timescale = tidal_tau,
        u_inflow = u_inflow_val,
        v_inflow = v_inflow_val,
        enable_voronoi = enable_voronoi,
        voronoi_n_units = voronoi_n_units,
        voronoi_prob_core = voronoi_p_core,
        voronoi_prob_shallow = voronoi_p_shallow,
        voronoi_prob_deep = voronoi_p_deep,
        voronoi_min_res_core_km = voronoi_min_core,
        voronoi_min_res_shallow_km = voronoi_min_shallow,
        voronoi_min_res_deep_km = voronoi_min_deep,
        voronoi_core_depth_min = voronoi_core_min,
        voronoi_core_depth_max = voronoi_core_max,
        voronoi_shallow_depth_min = voronoi_shallow_min,
        voronoi_shallow_depth_max = voronoi_shallow_max,
        voronoi_deep_depth_min = voronoi_deep_min,
        voronoi_deep_depth_max = voronoi_deep_max,
        voronoi_slope_weighting = voronoi_slope_wt,
        voronoi_slope_factor = voronoi_slope_fac,
        animate_hydro = anim_hydro,
        anim_variable = anim_var,
        anim_fps = anim_fps,
        anim_format = anim_fmt,
        anim_depth = anim_depth,
        anim_overlay_particles = anim_overlay_parts,
        anim_output_path = anim_out_path
    )

    # Executive Hardware & Workflow Summary Banner
    println("\n=================================================================")
    println(" ParticleTracking Regional Hydrodynamic & Transport Engine")
    println("=================================================================")
    println(" Active Configuration: $(config_file)")
    if opts.use_gpu
        cuda_ok = false
        cuda_name = "NVIDIA CUDA Device"
        cuda_vram = ""
        try
            if CUDA.functional()
                cuda_ok = true
                cuda_name = CUDA.name(CUDA.device())
                cuda_vram = " ($(round(CUDA.totalmem(CUDA.device()) / 1024^3, digits=1)) GB VRAM)"
            end
        catch
        end
        if cuda_ok
            println(" Compute Architecture: GPU ($(cuda_name)$(cuda_vram))")
            println(" Host CPU Worker Pool: $(Threads.nthreads()) threads (active for NetCDF I/O & preprocessing)")
        else
            println(" Compute Architecture: GPU requested, but CUDA unavailable (fallback = $(opts.fallback_to_cpu))")
        end
    else
        println(" Compute Architecture: CPU ($(Threads.nthreads()) worker threads)")
    end
    println(" Domain Discretization: $(opts.grid_size[1]) × $(opts.grid_size[2]) × $(opts.grid_size[3]) cells")
    println(" Geographic Extent:     Lon [$(opts.domain_lon[1]), $(opts.domain_lon[2])]°E, Lat [$(opts.domain_lat[1]), $(opts.domain_lat[2])]°N")
    if opts.enable_checkpoint
        println(" State Checkpointing:   Enabled (prefix: $(opts.checkpoint_prefix))")
    end
    println("=================================================================\n")

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

    # Record the fully-resolved configuration alongside the outputs.
    #
    # The TOML on disk is only the starting point: CLI overrides (--particles, --duration,
    # --output-dir, ...) and the defaults applied above all feed into `opts`. Writing the
    # realized state means every output can be traced back to the exact parameters that
    # produced it, even when the run was launched with ad-hoc flags.
    mkpath(opts.output_dir)
    resolved_path = joinpath(opts.output_dir, "resolved_config.toml")
    try
        open(resolved_path, "w") do io
            TOML.print(io, options_to_configuration(opts))
        end
        println("Resolved configuration written to: $(resolved_path)")
    catch err
        @warn "Could not write resolved configuration to $(resolved_path): " *
              "$(typeof(err).name.name)"
    end

    # Record which physical inputs actually backed this run.
    #
    # `resolved_config.toml` above records what was *asked for*; this records what was
    # *actually read*, with the cached file's size and SHA-256. Keeping the two separate is
    # deliberate: a run that degraded to analytic forcing, or that reused a stale cache, will
    # show up as a divergence between them, and that divergence is the thing worth seeing.
    if any(a -> startswith(a, "--data="), args) || opts.enable_checkpoint
        try
            prov_path = data_provenance(
                opts.output_dir;
                config_path = config_file,
                extra = Dict{String, Any}(
                    "scenario" => string(opts.scenario),
                    "projection_year" => opts.projection_year,
                    "seed" => opts.seed,
                ),
            )
            println("Data provenance written to: $(prov_path)")
        catch err
            @warn "Could not write data provenance: $(err)"
        end
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

        default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
        jld2_path, _ = resolve_hydro_model_path(opts, default_jld2)
        f_info = inspect_hydrodynamic_file(jld2_path, expected_stop_time = opts.sim_duration)

        s_res = if (opts.reuse_hydro || opts.animate_hydro) && f_info.exists && f_info.n_timesteps > 0
            if f_info.is_complete
                println("Reusing completed hydrodynamic flow solution from: $(jld2_path)")
            else
                println("Reusing hydrodynamic flow solution from: $(jld2_path)")
                println("  Notice: $(f_info.n_timesteps) snapshot(s) present (up to $(round(f_info.last_time / 3600.0, digits=2)) h).")
            end
            (simulation = nothing, jld2_output_path = jld2_path)
        else
            d_res = run_segment_data(opts = opts)
            g_res = run_segment_grid(opts = opts, bathy_file = d_res.bathy_file)
            m_res = run_segment_model(opts = opts, immersed_grid = g_res.immersed_grid,
                                      tau_x = d_res.tau_x, tau_y = d_res.tau_y)
            run_segment_climate(opts = opts, model = m_res.model)
            run_segment_simulation(opts = opts, model = m_res.model)
        end

        println("\nHydrodynamic solution ready: $(s_res.jld2_output_path)")

        if opts.animate_hydro
            render_hydrodynamic_animation_if_requested(opts, s_res.jld2_output_path)
        end

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
        s_res = run_segment_simulation(opts = opts)
        if opts.animate_hydro
            render_hydrodynamic_animation_if_requested(opts, s_res.jld2_output_path)
        end
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

    # Explicit "redraw every figure and visual from stored results" entry point. Distinct from
    # --viz in intent only: it names the intent (regenerate figures) and is what a user reaches
    # for after changing a plotting function, with no need to remember which segment owns it.
    if "--figures" in args || "--regenerate-figures" in args
        println("Regenerating all figures and visuals from stored run results...")
        render_all_figures(opts = opts)
        has_run = true
    end

    if !has_run
        if opts.animate_hydro && !isempty(opts.hydro_model_file) && isfile(opts.hydro_model_file)
            render_hydrodynamic_animation_if_requested(opts, opts.hydro_model_file)
        else
            println("No recognized execution flag provided.")
            display_help()
        end
    end
end

# Execute if run directly as script from command line
if is_cli_invocation
    main(ARGS)
end