"""
    configuration.jl

Centralized configuration manager for `ParticleTracking.jl`.
Provides robust parsing, serialization, validation, and conversion between
structured TOML/INI configuration files (such as `inputs/ParticleTracking.config`),
nested dictionaries, metadata records, and runtime `HydrodynamicOptions` structs.
"""

using TOML

"""
    find_default_config_path() -> String

Locate the default configuration file in the project workspace, checking
`inputs/ParticalTracking.config`, `inputs/ParticleTracking.config`,
`ParticalTracking.config`, or `ParticleTracking.config`.
"""
function find_default_config_path()::String
    candidates = [
        joinpath("inputs", "ParticalTracking.config"),
        joinpath("inputs", "ParticleTracking.config"),
        "ParticalTracking.config",
        "ParticleTracking.config"
    ]
    for p in candidates
        if isfile(p)
            return p
        end
    end
    return joinpath("inputs", "ParticleTracking.config")
end

"""
    HydrodynamicOptions

Runtime parameter specification for regional hydrodynamic circulation,
tidal harmonics, CMIP6 climate anomalies, and Lagrangian snow crab
(Chionoecetes opilio) larval transport modeling.

# Fields
- `domain_lon::Tuple{Float64, Float64}`: Longitude bounding box (min_lon, max_lon) in °E.
- `domain_lat::Tuple{Float64, Float64}`: Latitude bounding box (min_lat, max_lat) in °N.
- `domain_z::Tuple{Float64, Float64}`: Vertical depth range (z_min, z_max) in meters.
- `grid_size::Tuple{Int, Int, Int}`: Grid cell counts (Nx, Ny, Nz).
- `data_mode::Symbol`: `:real` (NOAA ERDDAP) or `:synthetic` (idealized shelf).
- `enable_tides::Bool`: Whether to include astronomical tidal body forcing (M2).
- `tidal_u_amp::Float64`: Zonal tidal velocity amplitude in m/s.
- `tidal_v_amp::Float64`: Meridional tidal velocity amplitude in m/s.
- `scenario::Symbol`: Climate scenario (:historical, :ssp126, :ssp245, :ssp585, :mhw).
- `projection_year::Int`: Target climate scenario year (e.g. 2050).
- `sim_dt::Float64`: Initial hydrodynamic integration time step in seconds.
- `sim_duration::Float64`: Total hydrodynamic simulation duration in seconds.
- `adaptive_cfl::Bool`: Whether to adaptively modulate hydrodynamic time stepping.
- `target_cfl::Float64`: Target CFL limit for the numerical wizard.
- `n_particles::Int`: Number of Lagrangian larval particles to track.
- `track_duration::Float64`: Total particle tracking drift duration in seconds.
- `track_dt::Float64`: Particle integration step in seconds.
- `diffusivity_h::Float64`: Horizontal turbulent eddy diffusivity in m^2/s.
- `diffusivity_v::Float64`: Vertical turbulent eddy diffusivity in m^2/s.
- `enable_dvm::Bool`: Whether larvae undergo active Diel Vertical Migration.
- `enable_molting::Bool`: Whether degree-day accumulation triggers stage molting.
- `min_seabed_depth::Float64`: Minimum bathymetric water depth for larval release (meters).
- `buffer_km::Float64`: Spatial buffer distance in kilometers (default 100.0 km) extending beyond administrative CFAs.
- `use_gpu::Bool`: Whether to run hydrodynamic equations on NVIDIA CUDA GPU.
- `fallback_to_cpu::Bool`: Whether to fall back to CPU if CUDA GPU is not functional.
- `interactive_map::Bool`: Whether to export interactive HTML5 Leaflet map.
- `enable_duckdb::Bool`: Whether to persist simulation data to DuckDB.
- `duckdb_path::String`: DuckDB file path.
- `config_file::String`: Source configuration file path.
- `output_dir::String`: Directory for simulation artifacts and figures.
- `input_dir::String`: Directory for raw and processed bathymetry/wind files.
- `seed::Int`: Random number generator seed.
"""
struct HydrodynamicOptions
    domain_lon        :: Tuple{Float64, Float64}
    domain_lat        :: Tuple{Float64, Float64}
    domain_z          :: Tuple{Float64, Float64}
    grid_size         :: Tuple{Int, Int, Int}
    data_mode         :: Symbol
    enable_tides      :: Bool
    tidal_u_amp       :: Float64
    tidal_v_amp       :: Float64
    scenario          :: Symbol
    projection_year   :: Int
    sim_dt            :: Float64
    sim_duration      :: Float64
    adaptive_cfl      :: Bool
    target_cfl        :: Float64
    n_particles       :: Int
    track_duration    :: Float64
    track_dt          :: Float64
    diffusivity_h     :: Float64
    diffusivity_v     :: Float64
    enable_dvm        :: Bool
    enable_molting    :: Bool
    min_seabed_depth  :: Float64
    buffer_km         :: Float64
    use_gpu           :: Bool
    fallback_to_cpu   :: Bool
    interactive_map   :: Bool
    enable_duckdb     :: Bool
    duckdb_path       :: String
    config_file       :: String
    output_dir        :: String
    input_dir         :: String
    seed              :: Int
end

function HydrodynamicOptions(;
    domain_lon       :: Tuple{Real, Real} = (-71.0, -53.0),
    domain_lat       :: Tuple{Real, Real} = (40.0, 48.5),
    domain_z         :: Tuple{Real, Real} = (-3500.0, 0.0),
    grid_size        :: Tuple{Int, Int, Int} = (50, 50, 10),
    data_mode        :: Symbol = :synthetic,
    enable_tides     :: Bool = true,
    tidal_u_amp      :: Real = 0.25,
    tidal_v_amp      :: Real = 0.12,
    scenario         :: Symbol = :ssp245,
    projection_year  :: Int = 2050,
    sim_dt           :: Real = 120.0,
    sim_duration     :: Real = 43200.0,
    adaptive_cfl     :: Bool = true,
    target_cfl       :: Real = 0.2,
    n_particles      :: Int = 100,
    track_duration   :: Real = 86400.0 * 5,
    track_dt         :: Real = 300.0,
    diffusivity_h    :: Real = 10.0,
    diffusivity_v    :: Real = 1e-4,
    enable_dvm       :: Bool = true,
    enable_molting   :: Bool = true,
    min_seabed_depth :: Real = 100.0,
    buffer_km        :: Real = 100.0,
    use_gpu          :: Bool = false,
    fallback_to_cpu  :: Bool = false,
    interactive_map  :: Bool = true,
    enable_duckdb    :: Bool = true,
    duckdb_path      :: AbstractString = joinpath("outputs", "particle_tracking.duckdb"),
    config_file      :: AbstractString = find_default_config_path(),
    output_dir       :: AbstractString = "outputs",
    input_dir        :: AbstractString = "inputs",
    seed             :: Int = 42
)
    return HydrodynamicOptions(
        (Float64(domain_lon[1]), Float64(domain_lon[2])),
        (Float64(domain_lat[1]), Float64(domain_lat[2])),
        (Float64(domain_z[1]), Float64(domain_z[2])),
        grid_size,
        data_mode,
        enable_tides,
        Float64(tidal_u_amp),
        Float64(tidal_v_amp),
        scenario,
        projection_year,
        Float64(sim_dt),
        Float64(sim_duration),
        adaptive_cfl,
        Float64(target_cfl),
        n_particles,
        Float64(track_duration),
        Float64(track_dt),
        Float64(diffusivity_h),
        Float64(diffusivity_v),
        enable_dvm,
        enable_molting,
        Float64(min_seabed_depth),
        Float64(buffer_km),
        use_gpu,
        fallback_to_cpu,
        interactive_map,
        enable_duckdb,
        String(duckdb_path),
        String(config_file),
        String(output_dir),
        String(input_dir),
        seed
    )
end

"""
    load_configuration(config_path::AbstractString = find_default_config_path()) -> Dict{String, Any}

Read and parse a centralized `ParticleTracking.config` file into a nested Julia Dictionary.
If the requested file does not exist, returns the default parameter configuration dictionary.

# Inputs
- `config_path::AbstractString`: Path to the `.config` (TOML format) file.

# Outputs
- `Dict{String, Any}`: Nested dictionary containing all sectioned parameter settings.
"""
function load_configuration(config_path::AbstractString = find_default_config_path())::Dict{String, Any}
    if isfile(config_path)
        try
            return TOML.parsefile(config_path)
        catch err
            @warn "Failed to parse configuration file at $(config_path): $(err). Using defaults."
            return get_default_configuration()
        end
    else
        return get_default_configuration()
    end
end

"""
    save_configuration(
        config_dict::AbstractDict,
        config_path::AbstractString = joinpath("inputs", "ParticalTracking.config")
    ) -> String

Serialize a nested configuration dictionary to a centralized `.config` file in TOML format.

# Inputs
- `config_dict::AbstractDict`: Dictionary of configuration sections and key-values.
- `config_path::AbstractString`: Target destination file path.

# Outputs
- `String`: Path to the written configuration file.
"""
function save_configuration(
    config_dict::AbstractDict,
    config_path::AbstractString = joinpath("inputs", "ParticalTracking.config")
)::String
    out_dir = dirname(abspath(config_path))
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    open(config_path, "w") do io
        TOML.print(io, config_dict; sorted = true)
    end
    return config_path
end

"""
    get_default_configuration() -> Dict{String, Any}

Generate a comprehensive dictionary of all default parameters across all modeling domains.
"""
function get_default_configuration()::Dict{String, Any}
    return Dict{String, Any}(
        "domain" => Dict{String, Any}(
            "lon_min" => -71.0,
            "lon_max" => -53.0,
            "lat_min" => 40.0,
            "lat_max" => 48.5,
            "z_min" => -3500.0,
            "z_max" => 0.0,
            "buffer_km" => 100.0
        ),
        "grid" => Dict{String, Any}(
            "nx" => 50,
            "ny" => 50,
            "nz" => 10
        ),
        "data" => Dict{String, Any}(
            "data_mode" => "synthetic",
            "bathy_dataset_id" => "etopo180",
            "wind_dataset_id" => "erdBSwinds1day",
            "wind_time_iso" => "2023-06-01T00:00:00Z",
            "inshore_depth" => -100.0,
            "shelf_slope" => 500.0
        ),
        "tides" => Dict{String, Any}(
            "enable_tides" => true,
            "constituents" => ["M2"],
            "tidal_u_amp" => 0.25,
            "tidal_v_amp" => 0.12,
            "tidal_period" => 44712.0,
            "tidal_phase" => 0.0
        ),
        "climate" => Dict{String, Any}(
            "scenario" => "ssp245",
            "projection_year" => 2050,
            "baseline_year" => 2015,
            "horizon_year" => 2050
        ),
        "hydrodynamics" => Dict{String, Any}(
            "sim_duration_hours" => 12.0,
            "sim_dt_seconds" => 120.0,
            "adaptive_cfl" => true,
            "target_cfl" => 0.2,
            "divergence_velocity_limit" => 20.0,
            "coriolis_latitude" => 44.5
        ),
        "biology" => Dict{String, Any}(
            "n_particles" => 100,
            "track_duration_days" => 5.0,
            "track_dt_seconds" => 300.0,
            "min_seabed_depth" => 100.0,
            "buffer_km" => 100.0,
            "diffusivity_h" => 10.0,
            "diffusivity_v" => 1e-4
        ),
        "dvm" => Dict{String, Any}(
            "enable_dvm" => true,
            "megalopa_day_depth" => -120.0,
            "megalopa_night_depth" => -60.0,
            "zoea2_depth_factor" => 1.2,
            "megalopa_swim_factor" => 1.5
        ),
        "molting_and_settlement" => Dict{String, Any}(
            "enable_molting" => true,
            "t_base" => 0.0,
            "dd_zoea1_to_zoea2" => 150.0,
            "dd_zoea2_to_megalopa" => 310.0,
            "dd_megalopa_to_settle" => 510.0,
            "mortality_base" => 0.02,
            "mortality_thermal_threshold" => 10.0,
            "mortality_thermal_sensitivity" => 0.015,
            "settlement_min_depth" => -250.0,
            "settlement_max_depth" => -50.0,
            "settlement_max_temp" => 6.0
        ),
        "storage" => Dict{String, Any}(
            "enable_duckdb" => true,
            "duckdb_path" => "outputs/particle_tracking.duckdb",
            "export_parquet" => false
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => false,
            "fallback_to_cpu" => true
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => true,
            "title" => "Scotian Shelf Snow Crab Larval Dispersal & Demographic Connectivity"
        ),
        "paths" => Dict{String, Any}(
            "output_dir" => "outputs",
            "input_dir" => "inputs",
            "seed" => 42
        )
    )
end

"""
    configuration_to_options(config_dict::AbstractDict; overrides...) -> HydrodynamicOptions

Convert a configuration dictionary into a validated `HydrodynamicOptions` instance,
allowing optional keyword parameter overrides.

# Inputs
- `config_dict::AbstractDict`: Parsed configuration dictionary.
- `overrides...`: Additional keyword arguments to override dictionary entries.

# Outputs
- `HydrodynamicOptions`: Constructed runtime options instance.
"""
function configuration_to_options(config_dict::AbstractDict; overrides...)
    # Helper to safely extract nested values with defaults
    function get_val(section::String, key::String, default_val)
        if haskey(config_dict, section) && haskey(config_dict[section], key)
            return config_dict[section][key]
        end
        return default_val
    end

    lon_min = Float64(get_val("domain", "lon_min", -71.0))
    lon_max = Float64(get_val("domain", "lon_max", -53.0))
    lat_min = Float64(get_val("domain", "lat_min", 40.0))
    lat_max = Float64(get_val("domain", "lat_max", 48.5))
    z_min   = Float64(get_val("domain", "z_min", -3500.0))
    z_max   = Float64(get_val("domain", "z_max", 0.0))

    buffer_km = Float64(get_val("domain", "buffer_km", get_val("biology", "buffer_km", 100.0)))

    nx = Int(get_val("grid", "nx", 50))
    ny = Int(get_val("grid", "ny", 50))
    nz = Int(get_val("grid", "nz", 10))

    data_mode = Symbol(get_val("data", "data_mode", "synthetic"))

    enable_tides = Bool(get_val("tides", "enable_tides", true))
    tidal_u = Float64(get_val("tides", "tidal_u_amp", 0.25))
    tidal_v = Float64(get_val("tides", "tidal_v_amp", 0.12))

    scenario = Symbol(get_val("climate", "scenario", "ssp245"))
    proj_year = Int(get_val("climate", "projection_year", 2050))

    sim_hours = Float64(get_val("hydrodynamics", "sim_duration_hours", 12.0))
    sim_dt    = Float64(get_val("hydrodynamics", "sim_dt_seconds", 120.0))
    adapt_cfl = Bool(get_val("hydrodynamics", "adaptive_cfl", true))
    tgt_cfl   = Float64(get_val("hydrodynamics", "target_cfl", 0.2))

    n_parts   = Int(get_val("biology", "n_particles", 100))
    track_days = Float64(get_val("biology", "track_duration_days", 5.0))
    track_dt  = Float64(get_val("biology", "track_dt_seconds", 300.0))
    min_depth = Float64(get_val("biology", "min_seabed_depth", 100.0))
    diff_h    = Float64(get_val("biology", "diffusivity_h", 10.0))
    diff_v    = Float64(get_val("biology", "diffusivity_v", 1e-4))

    enable_dvm = Bool(get_val("dvm", "enable_dvm", true))
    enable_molting = Bool(get_val("molting_and_settlement", "enable_molting", true))

    use_gpu      = Bool(get_val("hardware", "use_gpu", false))
    fallback_cpu = Bool(get_val("hardware", "fallback_to_cpu", true))
    interactive  = Bool(get_val("visualization", "interactive_map", true))

    enable_duckdb = Bool(get_val("storage", "enable_duckdb", true))
    duckdb_path   = String(get_val("storage", "duckdb_path", "outputs/particle_tracking.duckdb"))

    output_dir = String(get_val("paths", "output_dir", "outputs"))
    input_dir  = String(get_val("paths", "input_dir", "inputs"))
    seed       = Int(get_val("paths", "seed", 42))

    # Construct HydrodynamicOptions with overrides applied
    return HydrodynamicOptions(;
        domain_lon = (lon_min, lon_max),
        domain_lat = (lat_min, lat_max),
        domain_z = (z_min, z_max),
        grid_size = (nx, ny, nz),
        data_mode = data_mode,
        enable_tides = enable_tides,
        tidal_u_amp = tidal_u,
        tidal_v_amp = tidal_v,
        scenario = scenario,
        projection_year = proj_year,
        sim_dt = sim_dt,
        sim_duration = sim_hours * 3600.0,
        adaptive_cfl = adapt_cfl,
        target_cfl = tgt_cfl,
        n_particles = n_parts,
        track_duration = track_days * 86400.0,
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
        duckdb_path = duckdb_path,
        output_dir = output_dir,
        input_dir = input_dir,
        seed = seed,
        overrides...
    )
end

"""
    options_to_configuration(opts::HydrodynamicOptions) -> Dict{String, Any}

Extract a complete configuration dictionary from a `HydrodynamicOptions` instance.
"""
function options_to_configuration(opts::HydrodynamicOptions)::Dict{String, Any}
    return Dict{String, Any}(
        "domain" => Dict{String, Any}(
            "lon_min" => opts.domain_lon[1],
            "lon_max" => opts.domain_lon[2],
            "lat_min" => opts.domain_lat[1],
            "lat_max" => opts.domain_lat[2],
            "z_min" => opts.domain_z[1],
            "z_max" => opts.domain_z[2],
            "buffer_km" => opts.buffer_km
        ),
        "grid" => Dict{String, Any}(
            "nx" => opts.grid_size[1],
            "ny" => opts.grid_size[2],
            "nz" => opts.grid_size[3]
        ),
        "data" => Dict{String, Any}(
            "data_mode" => string(opts.data_mode),
            "inshore_depth" => -100.0,
            "shelf_slope" => 500.0
        ),
        "tides" => Dict{String, Any}(
            "enable_tides" => opts.enable_tides,
            "tidal_u_amp" => opts.tidal_u_amp,
            "tidal_v_amp" => opts.tidal_v_amp
        ),
        "climate" => Dict{String, Any}(
            "scenario" => string(opts.scenario),
            "projection_year" => opts.projection_year
        ),
        "hydrodynamics" => Dict{String, Any}(
            "sim_duration_hours" => opts.sim_duration / 3600.0,
            "sim_dt_seconds" => opts.sim_dt,
            "adaptive_cfl" => opts.adaptive_cfl,
            "target_cfl" => opts.target_cfl
        ),
        "biology" => Dict{String, Any}(
            "n_particles" => opts.n_particles,
            "track_duration_days" => opts.track_duration / 86400.0,
            "track_dt_seconds" => opts.track_dt,
            "min_seabed_depth" => opts.min_seabed_depth,
            "buffer_km" => opts.buffer_km,
            "diffusivity_h" => opts.diffusivity_h,
            "diffusivity_v" => opts.diffusivity_v
        ),
        "dvm" => Dict{String, Any}(
            "enable_dvm" => opts.enable_dvm
        ),
        "molting_and_settlement" => Dict{String, Any}(
            "enable_molting" => opts.enable_molting
        ),
        "storage" => Dict{String, Any}(
            "enable_duckdb" => opts.enable_duckdb,
            "duckdb_path" => opts.duckdb_path
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => opts.use_gpu,
            "fallback_to_cpu" => opts.fallback_to_cpu
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => opts.interactive_map
        ),
        "paths" => Dict{String, Any}(
            "output_dir" => opts.output_dir,
            "input_dir" => opts.input_dir,
            "seed" => opts.seed
        )
    )
end
