"""
    cli_parser.jl

ArgParse-based command-line interface parser for ParticleTracking.
Replaces manual argument parsing with proper, validated argument handling.
"""

module CLIParser

using ArgParse

"""
    create_argparse_settings() -> ArgParseSettings

Create the complete ArgParse settings for ParticleTracking CLI.
"""
function create_argparse_settings()::ArgParseSettings
    s = ArgParseSettings(
        description = "Hydrodynamic Model & Particle Tracking CLI Interface",
        prog = "ParticleTrackingRun.jl",
        allow_ambiguous_args = false,
        error_on_conflict = true,
        add_version = true,
        version = "0.1.0"
    )

    # =============================================================================
    # Execution Mode Arguments
    # =============================================================================
    @add_arg_table! s begin
        "--all"
            help = "Execute entire end-to-end workflow pipeline"
            action = :store_true
        "--quick", "-q"
            help = "Fast debug mode with coarse resolution and short duration"
            action = :store_true
        "--segment"
            help = "Run a single workflow segment"
            arg_type = String
            default = ""
            choices = ["data", "grid", "model", "climate", "sim", "track", "metrics", "viz", "all"]
    end

    # =============================================================================
    # Decoupled Hydrodynamics & Multi-Cohort Tracking
    # =============================================================================
    @add_arg_table! s begin
        "--hydro-model"
            help = "Target hydrodynamic JLD2 model file (output or input)"
            arg_type = String
            default = ""
        "--hydro-only"
            help = "Run hydrodynamics only (Segments 1-5) and save to --hydro-model"
            action = :store_true
        "--track-only"
            help = "Run larval tracking only (Segments 6-8) using --hydro-model"
            action = :store_true
        "--reuse-hydro"
            help = "Reuse existing --hydro-model flow file if completed; else simulate"
            action = :store_true
        "--restart"
            help = "Gracefully resume hydrodynamic simulation from latest checkpoint (default: true)"
            action = :store_true
        "--no-restart", "--force-new"
            help = "Do not resume from checkpoint; overwrite and start fresh from t=0"
            action = :store_true
        "--checkpoint"
            help = "Enable prognostic state checkpoints (default: true)"
            action = :store_true
        "--no-checkpoint"
            help = "Disable prognostic state checkpoints"
            action = :store_true
        "--checkpoint-interval"
            help = "State checkpoint interval in seconds, hours (e.g. 6h), or days (e.g. 1d)"
            arg_type = String
            default = ""
        "--checkpoints-dir"
            help = "Custom directory for storing restart checkpoints"
            arg_type = String
            default = ""
        "--run-id"
            help = "Unique cohort run identifier for DuckDB persistence and figures"
            arg_type = String
            default = ""
    end

    # =============================================================================
    # Individual Segment Flags
    # =============================================================================
    @add_arg_table! s begin
        "--data"
            help = "Run environmental data ingestion (Segment 1)"
            action = :store_true
        "--grid"
            help = "Run grid and immersed boundary construction (Segment 2)"
            action = :store_true
        "--model"
            help = "Run hydrodynamic model setup with tidal forcing (Segment 3)"
            action = :store_true
        "--climate"
            help = "Run CMIP6 climate scenario & larval thermal biology (Segment 4)"
            action = :store_true
        "--sim", "--simulation"
            help = "Run Oceananigans hydrodynamic time integration (Segment 5)"
            action = :store_true
        "--track", "--tracking"
            help = "Run Lagrangian larval particle tracking (Segment 6)"
            action = :store_true
        "--metrics"
            help = "Compute retention, empirical diffusion & connectivity (Segment 7)"
            action = :store_true
        "--viz", "--visualize"
            help = "Generate CairoMakie spatial figures and interactive HTML map (Segment 8)"
            action = :store_true
    end

    # =============================================================================
    # DuckDB Analytical Storage & Model Averaging
    # =============================================================================
    @add_arg_table! s begin
        "--duckdb"
            help = "Enable DuckDB storage archiving (default: true)"
            action = :store_true
        "--no-duckdb"
            help = "Disable DuckDB storage archiving"
            action = :store_true
        "--db-path"
            help = "Custom DuckDB database path"
            arg_type = String
            default = "outputs/particle_tracking.duckdb"
        "--list-runs"
            help = "Query and display all simulation runs archived in DuckDB"
            action = :store_true
        "--compare-scenarios"
            help = "Query and display multi-scenario comparison metrics"
            action = :store_true
        "--model-average"
            help = "Compute ensemble model-averaged connectivity and recruitment"
            action = :store_true
    end

    # =============================================================================
    # Centralized Configuration
    # =============================================================================
    @add_arg_table! s begin
        "--config"
            help = "Path to centralized .toml file"
            arg_type = String
            default = "inputs/ParticleTracking.toml"
        "--save-config"
            help = "Export active parameter options to .toml file and exit"
            nargs = '?'
            arg_type = String
            default = ""
    end

    # =============================================================================
    # Climate Scenarios & Thermal Biology
    # =============================================================================
    @add_arg_table! s begin
        "--scenario"
            help = "Climate scenario: historical, ssp126, ssp245, ssp585, mhw, climatology"
            arg_type = String
            default = "ssp245"
        "--year"
            help = "Climate projection horizon year"
            arg_type = Int
            default = 2050
        "--heat-flux"
            help = "Summer atmospheric heat flux in W/m²"
            arg_type = Float64
            default = 50.0
    end

    # =============================================================================
    # Computational Architecture
    # =============================================================================
    @add_arg_table! s begin
        "--gpu", "--cuda"
            help = "Enable NVIDIA CUDA GPU acceleration for hydrodynamics"
            action = :store_true
        "--cpu"
            help = "Execute on multi-threaded CPU"
            action = :store_true
        "--fallback-cpu"
            help = "Automatically fall back to CPU if CUDA GPU is not functional"
            action = :store_true
    end

    # =============================================================================
    # Visualization & Hydrodynamic Animation Options
    # =============================================================================
    @add_arg_table! s begin
        "--interactive"
            help = "Export standalone interactive HTML5 Leaflet map (default)"
            action = :store_true
        "--no-interactive"
            help = "Disable interactive HTML map generation"
            action = :store_true
        "--animate-hydro", "--anim-hydro"
            help = "Render animated MP4/GIF video of hydrodynamic fields"
            action = :store_true
        "--no-animate-hydro"
            help = "Disable hydrodynamic animation generation (default)"
            action = :store_true
        "--anim-variable"
            help = "Field variable: dashboard, temperature (T), advection (speed), diffusion (kappa), salinity (S), viscosity (nu), stratification (N2), density (rho), vorticity (zeta), elevation (eta), w, richardson (Ri)"
            arg_type = String
            default = "dashboard"
        "--anim-fps"
            help = "Animation playback framerate in frames/sec"
            arg_type = Int
            default = 10
        "--anim-format"
            help = "Video container format: mp4, gif"
            arg_type = String
            default = "mp4"
        "--anim-depth"
            help = "Depth slice in meters for 2D fields"
            arg_type = Float64
            default = -2.5
        "--anim-output", "--anim-file"
            help = "Designate custom output video path"
            arg_type = String
            default = ""
        "--anim-overlay-particles", "--anim-particles"
            help = "Synchronously overlay Lagrangian larvae drifting with currents"
            action = :store_true
    end

    # =============================================================================
    # Spatial Domain & Grid Discretization
    # =============================================================================
    @add_arg_table! s begin
        "--lon"
            help = "Longitude bounding range in degrees East (format: min,max)"
            arg_type = String
            default = ""
        "--lat"
            help = "Latitude bounding range in degrees North (format: min,max)"
            arg_type = String
            default = ""
        "--depth-range", "--z"
            help = "Vertical depth range in meters (format: min,max)"
            arg_type = String
            default = ""
        "--grid"
            help = "Grid cell dimensions (format: nx,ny,nz)"
            arg_type = String
            default = ""
        "--nx"
            help = "Zonal grid cells"
            arg_type = Int
            default = 0
        "--ny"
            help = "Meridional grid cells"
            arg_type = Int
            default = 0
        "--nz"
            help = "Vertical grid layers"
            arg_type = Int
            default = 0
        "--res-scale"
            help = "Resolution scaling factor (e.g. 2.5 for ~2.5km calibrated grid)"
            arg_type = Float64
            default = 1.0
        "--stretched-z"
            help = "Enable hyperbolic tangent vertical layer stretching"
            action = :store_true
        "--uniform-z"
            help = "Use uniform vertical layer thicknesses"
            action = :store_true
        "--z-file", "--vertical-grid-file"
            help = "Path to vertical layer coordinate CSV file"
            arg_type = String
            default = "inputs/scotian_shelf_vertical_grid.csv"
    end

    # =============================================================================
    # Environmental Data & Forcing
    # =============================================================================
    @add_arg_table! s begin
        "--real"
            help = "Fetch real-world NOAA ERDDAP bathymetry and winds"
            action = :store_true
        "--synthetic"
            help = "Generate idealized synthetic shelf data"
            action = :store_true
        "--tides"
            help = "Enable astronomical tidal body forcing"
            action = :store_true
        "--no-tides"
            help = "Disable tidal body forcing"
            action = :store_true
        "--tidal-u"
            help = "Semi-major tidal current amplitude in m/s"
            arg_type = Float64
            default = 0.25
        "--tidal-v"
            help = "Semi-minor tidal current amplitude in m/s"
            arg_type = Float64
            default = 0.12
        "--era5-forcing", "--era5"
            help = "Enable NumericalEarth ERA5 atmospheric surface forcing"
            action = :store_true
        "--no-era5"
            help = "Use analytical wind stress forcing"
            action = :store_true
        "--obc"
            help = "Enable NumericalEarth GLORYS12V1 open boundary conditions"
            action = :store_true
        "--no-obc"
            help = "Disable open boundary conditions"
            action = :store_true
    end

    # =============================================================================
    # Climate Scenarios & Thermal Biology
    # =============================================================================
    @add_arg_table! s begin
        "--scenario"
            help = "Climate scenario: historical, ssp126, ssp245, ssp585, mhw, climatology"
            arg_type = String
            default = "ssp245"
        "--year"
            help = "Climate projection horizon year"
            arg_type = Int
            default = 2050
    end

    # =============================================================================
    # Hydrodynamic Simulation
    # =============================================================================
    @add_arg_table! s begin
        "--duration"
            help = "Hydrodynamic simulation duration in hours"
            arg_type = Float64
            default = 12.0
        "--sim-duration"
            help = "Hydrodynamic simulation duration in seconds"
            arg_type = Float64
            default = 43200.0
        "--sim-dt"
            help = "Initial hydrodynamic time step in seconds"
            arg_type = Float64
            default = 120.0
        "--adaptive-cfl"
            help = "Enable adaptive CFL time stepping"
            action = :store_true
        "--no-adaptive-cfl"
            help = "Disable adaptive CFL time stepping"
            action = :store_true
        "--target-cfl"
            help = "Target advective Courant-Friedrichs-Lewy limit"
            arg_type = Float64
            default = 0.2
    end

    # =============================================================================
    # Lagrangian Particle Tracking & Larval Ecology
    # =============================================================================
    @add_arg_table! s begin
        "--particles"
            help = "Number of larvae to initialize and track"
            arg_type = Int
            default = 100
        "--track-duration"
            help = "Cohort tracking duration in days"
            arg_type = Float64
            default = 5.0
        "--track-dt"
            help = "Lagrangian integration time step in seconds"
            arg_type = Float64
            default = 300.0
        "--min-depth"
            help = "Minimum water depth for larval placement"
            arg_type = Float64
            default = 100.0
        "--buffer-km", "--buffer", "--buf"
            help = "Spatial buffer beyond stratum boundaries in km"
            arg_type = Float64
            default = 100.0
        "--dvm"
            help = "Enable stage-dependent Diel Vertical Migration"
            action = :store_true
        "--no-dvm"
            help = "Disable Diel Vertical Migration"
            action = :store_true
        "--molting"
            help = "Enable degree-day thermal molting & mortality"
            action = :store_true
        "--no-molting"
            help = "Disable thermal molting"
            action = :store_true
        "--diff-h", "--diffusivity-h"
            help = "Horizontal turbulent diffusivity in m^2/s"
            arg_type = Float64
            default = 10.0
        "--diff-v", "--diffusivity-v"
            help = "Vertical turbulent diffusivity in m^2/s"
            arg_type = Float64
            default = 1e-4
        "--release-mode"
            help = "Release depth mode: bottom, range, surface"
            arg_type = String
            default = "bottom"
        "--ascent"
            help = "Enable post-hatch vertical ascent toward surface"
            action = :store_true
        "--no-ascent"
            help = "Disable initial vertical ascent"
            action = :store_true
        "--ascent-speed"
            help = "Vertical ascent swimming speed in m/s"
            arg_type = Float64
            default = 0.010
        "--ascent-target"
            help = "Target depth in meters for ascent completion"
            arg_type = Float64
            default = -10.0
    end

    # =============================================================================
    # Voronoi Tessellation Options
    # =============================================================================
    @add_arg_table! s begin
        "--voronoi-units", "--n-units"
            help = "Number of Voronoi units/centroids to generate"
            arg_type = Int
            default = 5000
        "--voronoi-prob-core"
            help = "Sampling probability weight for core depth stratum"
            arg_type = Float64
            default = 0.8
        "--voronoi-prob-shallow"
            help = "Sampling probability weight for shallow stratum"
            arg_type = Float64
            default = 0.1
        "--voronoi-prob-deep"
            help = "Sampling probability weight for deep stratum"
            arg_type = Float64
            default = 0.1
        "--voronoi-min-core"
            help = "Minimum Poisson-disc separation distance in core zone (km)"
            arg_type = Float64
            default = 1.5
        "--voronoi-min-shallow"
            help = "Minimum separation distance in shallow zone (km)"
            arg_type = Float64
            default = 5.0
        "--voronoi-min-deep"
            help = "Minimum separation distance in deep zone (km)"
            arg_type = Float64
            default = 10.0
    end

    # =============================================================================
    # I/O & Environment
    # =============================================================================
    @add_arg_table! s begin
        "--output-dir"
            help = "Directory for output figures and datasets"
            arg_type = String
            default = "outputs"
        "--input-dir"
            help = "Directory for input NetCDF caches"
            arg_type = String
            default = "inputs"
        "--seed"
            help = "Random number generator seed"
            arg_type = Int
            default = 42
        "--help", "-h"
            help = "Display this help documentation"
            action = :store_true
    end

    # =============================================================================
    # Checkpoint Options
    # =============================================================================
    @add_arg_table! s begin
        "--checkpoint-prefix", "--cp-prefix"
            help = "Checkpoint filename prefix"
            arg_type = String
            default = ""
        "--checkpoint-cleanup"
            help = "Retain only latest checkpoint on disk"
            action = :store_true
        "--no-checkpoint-cleanup"
            help = "Keep all intermediate checkpoints"
            action = :store_true
        "--obc-source"
            help = "Open boundary condition source"
            arg_type = String
            default = "glorys12v1"
        "--obc-type"
            help = "Open boundary condition type"
            arg_type = String
            default = "flather_chapman"
    end

    return s
end

"""
    parse_cli_args(args::Vector{String} = ARGS) -> Dict{String, Any}

Parse command-line arguments and return a dictionary of parsed values.
"""
function parse_cli_args(args::Vector{String} = ARGS)::Dict{String, Any}
    s = create_argparse_settings()
    parsed = parse_args(args, s)
    return parsed
end

"""
    build_options_from_parsed(parsed::Dict{String, Any}) -> HydrodynamicOptions

Build HydrodynamicOptions from parsed ArgParse arguments.
This function consolidates all the logic from the manual parsing.
"""
function build_options_from_parsed(parsed::Dict{String, Any})::Dict{String, Any}
    # This will be implemented to convert parsed args to HydrodynamicOptions
    # For now, return the parsed dict for inspection
    return parsed
end

# Export
export create_argparse_settings, parse_cli_args, build_options_from_parsed

end # module CLIParser