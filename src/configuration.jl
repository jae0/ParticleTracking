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
        joinpath("inputs", "ParticleTracking.config"),
        "ParticleTracking.config",
        joinpath("inputs", "ParticalTracking.config"),
        "ParticalTracking.config"
    ]
    for p in candidates
        if isfile(p)
            return p
        end
    end
    return joinpath("inputs", "ParticleTracking.config")
end

"""
    resolve_config_name(config_file::AbstractString = "") -> String

Extract the base configuration name from a config file path, or return
the default configuration name (`"ParticleTracking"`) if unspecified or empty.

# Inputs
- `config_file::AbstractString`: Path to a `.config` file or empty string.

# Outputs
- `String`: Clean configuration identifier without directory or extension.
"""
function resolve_config_name(config_file::AbstractString = "")::String
    cfg_path = isempty(strip(config_file)) ? find_default_config_path() : String(config_file)
    base = basename(cfg_path)
    stem = splitext(base)[1]
    return isempty(stem) ? "ParticleTracking" : stem
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
- `buffer_km::Float64`: Spatial buffer distance in km (default 100.0 km) extending beyond CFAs.
- `release_depth_mode::Symbol`: Vertical release mode (:bottom, :range, :surface).
- `bottom_release_offset::Tuple{Float64, Float64}`: Elevation range above seabed (meters).
- `enable_initial_ascent::Bool`: Whether larvae perform active post-hatch vertical ascent.
- `ascent_speed::Float64`: Directed upward swimming speed during ascent (m/s).
- `ascent_target_depth::Float64`: Target epipelagic depth for ascent completion (meters).
- `use_gpu::Bool`: Whether to run hydrodynamic equations on NVIDIA CUDA GPU.
- `fallback_to_cpu::Bool`: Whether to fall back to CPU if CUDA GPU is not functional.
- `interactive_map::Bool`: Whether to export interactive HTML5 Leaflet map.
- `enable_duckdb::Bool`: Whether to persist simulation data to DuckDB.
- `duckdb_path::String`: DuckDB file path.
- `config_file::String`: Source configuration file path.
- `output_dir::String`: Directory for simulation artifacts and figures.
- `input_dir::String`: Directory for raw and processed bathymetry/wind files.
- `seed::Int`: Random number generator seed.
- `hydro_model_file::String`: Target hydrodynamic JLD2 checkpoint path (input/output).
- `hydro_only::Bool`: Run hydrodynamics only and persist to JLD2 checkpoint.
- `track_only::Bool`: Track larvae only using flow fields from hydro_model_file.
- `reuse_hydro::Bool`: Reuse existing hydro_model_file if present on disk; else run.
- `enable_checkpoint::Bool`: Whether to attach Oceananigans Checkpointer for state serialization.
- `checkpoint_schedule::Float64`: State checkpoint frequency in seconds (0.0 = match output schedule).
- `checkpoint_dir::String`: Directory where state checkpoints are archived.
- `checkpoint_cleanup::Bool`: Whether older intermediate checkpoints are automatically deleted.
- `auto_restart::Bool`: Automatically detect and pick up from existing checkpoints if available.
- `run_id::String`: Unique cohort identifier for DuckDB persistence and figures.
- `vertical_stretching_mode::Symbol`: Vertical grid stretching mode (:tanh, :uniform, or :csv).
- `vertical_grid_file::String`: Path to CSV file specifying vertical layer boundary faces.
- `resolution_scale::Float64`: Horizontal grid resolution scale multiplier / divisor.
- `atmospheric_source::Symbol`: Atmospheric forcing provider (:era5, :synthetic, :climatology).
- `ocean_boundary_source::Symbol`: Open boundary hydrographic data provider (:glorys12v1, :synthetic).
- `obc_type::Symbol`: Lateral open boundary formulation (:flather_chapman, :radiation, :clamped).
- `enable_voronoi::Bool`: Whether to compute multi-resolution Voronoi tessellation analysis.
- `voronoi_n_units::Int`: Number of Voronoi units/centroids to generate across strata.
- `voronoi_prob_core::Float64`: Sampling probability weight for core depth stratum (50 to 350 m).
- `voronoi_prob_shallow::Float64`: Sampling probability weight for shallow stratum (0 to 50 m).
- `voronoi_prob_deep::Float64`: Sampling probability weight for deep stratum (> 350 m).
- `voronoi_min_res_core_km::Float64`: Minimum Poisson-disc separation distance in core zone (km).
- `voronoi_min_res_shallow_km::Float64`: Minimum separation distance in shallow zone (km).
- `voronoi_min_res_deep_km::Float64`: Minimum separation distance in deep zone (km).
"""
struct HydrodynamicOptions
    domain_lon               :: Tuple{Float64, Float64}
    domain_lat               :: Tuple{Float64, Float64}
    domain_z                 :: Tuple{Float64, Float64}
    grid_size                :: Tuple{Int, Int, Int}
    data_mode                :: Symbol
    enable_tides             :: Bool
    tidal_u_amp              :: Float64
    tidal_v_amp              :: Float64
    scenario                 :: Symbol
    projection_year          :: Int
    sim_dt                   :: Float64
    sim_duration             :: Float64
    adaptive_cfl             :: Bool
    target_cfl               :: Float64
    surface_heat_flux        :: Float64
    n_particles              :: Int
    track_duration           :: Float64
    track_dt                 :: Float64
    diffusivity_h            :: Float64
    diffusivity_v            :: Float64
    enable_dvm               :: Bool
    enable_molting           :: Bool
    min_seabed_depth         :: Float64
    buffer_km                :: Float64
    release_depth_mode       :: Symbol
    bottom_release_offset    :: Tuple{Float64, Float64}
    enable_initial_ascent    :: Bool
    ascent_speed             :: Float64
    ascent_target_depth      :: Float64
    use_gpu                  :: Bool
    fallback_to_cpu          :: Bool
    interactive_map          :: Bool
    enable_duckdb            :: Bool
    duckdb_path              :: String
    config_file              :: String
    output_dir               :: String
    input_dir                :: String
    seed                     :: Int
    hydro_model_file         :: String
    hydro_only               :: Bool
    track_only               :: Bool
    reuse_hydro              :: Bool
    enable_checkpoint        :: Bool
    checkpoint_prefix        :: String
    checkpoint_schedule      :: Float64
    checkpoint_dir           :: String
    checkpoint_cleanup       :: Bool
    auto_restart             :: Bool
    run_id                   :: String
    vertical_stretching_mode :: Symbol
    vertical_grid_file       :: String
    resolution_scale         :: Float64
    atmospheric_source       :: Symbol
    ocean_boundary_source    :: Symbol
    obc_type                 :: Symbol
    enable_voronoi           :: Bool
    voronoi_n_units          :: Int
    voronoi_prob_core        :: Float64
    voronoi_prob_shallow     :: Float64
    voronoi_prob_deep        :: Float64
    voronoi_min_res_core_km  :: Float64
    voronoi_min_res_shallow_km :: Float64
    voronoi_min_res_deep_km  :: Float64
    animate_hydro            :: Bool
    anim_variable            :: Symbol
    anim_fps                 :: Int
    anim_format              :: String
    anim_depth               :: Float64
    anim_overlay_particles   :: Bool
end

function HydrodynamicOptions(;
    domain_lon            :: Tuple{Real, Real} = (-71.0, -53.0),
    domain_lat            :: Tuple{Real, Real} = (40.0, 48.5),
    domain_z              :: Tuple{Real, Real} = (-3500.0, 0.0),
    grid_size             :: Tuple{Int, Int, Int} = (50, 50, 10),
    data_mode             :: Symbol = :synthetic,
    enable_tides          :: Bool = true,
    tidal_u_amp           :: Real = 0.25,
    tidal_v_amp           :: Real = 0.12,
    scenario              :: Symbol = :ssp245,
    projection_year       :: Int = 2050,
    sim_dt                :: Real = 120.0,
    sim_duration          :: Real = 43200.0,
    adaptive_cfl          :: Bool = true,
    target_cfl            :: Real = 0.2,
    surface_heat_flux     :: Real = 50.0,
    n_particles           :: Int = 100,
    track_duration        :: Real = 86400.0 * 5,
    track_dt              :: Real = 300.0,
    diffusivity_h         :: Real = 10.0,
    diffusivity_v         :: Real = 1e-4,
    enable_dvm            :: Bool = true,
    enable_molting        :: Bool = true,
    min_seabed_depth      :: Real = 100.0,
    buffer_km             :: Real = 100.0,
    release_depth_mode    :: Symbol = :bottom,
    bottom_release_offset :: Tuple{Real, Real} = (0.5, 3.0),
    enable_initial_ascent :: Bool = true,
    ascent_speed          :: Real = 0.010,
    ascent_target_depth   :: Real = -10.0,
    use_gpu               :: Bool = false,
    fallback_to_cpu       :: Bool = false,
    interactive_map       :: Bool = true,
    enable_duckdb         :: Bool = true,
    duckdb_path           :: AbstractString = joinpath("outputs", "particle_tracking.duckdb"),
    config_file           :: AbstractString = find_default_config_path(),
    output_dir            :: AbstractString = "outputs",
    input_dir             :: AbstractString = "inputs",
    seed                  :: Int = 42,
    hydro_model_file      :: AbstractString = "",
    hydro_only            :: Bool = false,
    track_only            :: Bool = false,
    reuse_hydro           :: Bool = false,
    enable_checkpoint     :: Bool = true,
    checkpoint_prefix     :: AbstractString = "",
    checkpoint_schedule   :: Real = 0.0,
    checkpoint_dir        :: AbstractString = "",
    checkpoint_cleanup    :: Bool = true,
    auto_restart             :: Bool = true,
    run_id                   :: AbstractString = "",
    vertical_stretching_mode :: Symbol = :tanh,
    vertical_grid_file       :: AbstractString = joinpath("inputs", "scotian_shelf_vertical_grid.csv"),
    resolution_scale         :: Real = 1.0,
    atmospheric_source       :: Symbol = :era5,
    ocean_boundary_source    :: Symbol = :glorys12v1,
    obc_type                 :: Symbol = :flather_chapman,
    enable_voronoi           :: Bool = false,
    voronoi_n_units          :: Int = 5000,
    voronoi_prob_core        :: Real = 0.8,
    voronoi_prob_shallow     :: Real = 0.1,
    voronoi_prob_deep        :: Real = 0.1,
    voronoi_min_res_core_km  :: Real = 1.5,
    voronoi_min_res_shallow_km :: Real = 5.0,
    voronoi_min_res_deep_km  :: Real = 10.0,
    animate_hydro            :: Bool = false,
    anim_variable            :: Symbol = :dashboard,
    anim_fps                 :: Int = 10,
    anim_format              :: AbstractString = "mp4",
    anim_depth               :: Real = -2.5,
    anim_overlay_particles   :: Bool = false
)
    resolved_cp_prefix = if !isempty(strip(checkpoint_prefix)) && checkpoint_prefix != "checkpoint"
        String(checkpoint_prefix)
    else
        cfg_name = resolve_config_name(config_file)
        "checkpoint_$(cfg_name)"
    end

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
        Float64(surface_heat_flux),
        n_particles,
        Float64(track_duration),
        Float64(track_dt),
        Float64(diffusivity_h),
        Float64(diffusivity_v),
        enable_dvm,
        enable_molting,
        Float64(min_seabed_depth),
        Float64(buffer_km),
        release_depth_mode,
        (Float64(bottom_release_offset[1]), Float64(bottom_release_offset[2])),
        enable_initial_ascent,
        Float64(ascent_speed),
        Float64(ascent_target_depth),
        use_gpu,
        fallback_to_cpu,
        interactive_map,
        enable_duckdb,
        String(duckdb_path),
        String(config_file),
        String(output_dir),
        String(input_dir),
        seed,
        String(hydro_model_file),
        hydro_only,
        track_only,
        reuse_hydro,
        enable_checkpoint,
        resolved_cp_prefix,
        Float64(checkpoint_schedule),
        String(checkpoint_dir),
        checkpoint_cleanup,
        auto_restart,
        String(run_id),
        vertical_stretching_mode,
        String(vertical_grid_file),
        Float64(resolution_scale),
        atmospheric_source,
        ocean_boundary_source,
        obc_type,
        enable_voronoi,
        voronoi_n_units,
        Float64(voronoi_prob_core),
        Float64(voronoi_prob_shallow),
        Float64(voronoi_prob_deep),
        Float64(voronoi_min_res_core_km),
        Float64(voronoi_min_res_shallow_km),
        Float64(voronoi_min_res_deep_km),
        animate_hydro,
        anim_variable,
        anim_fps,
        String(anim_format),
        Float64(anim_depth),
        anim_overlay_particles
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
function load_configuration(
    config_path::AbstractString = find_default_config_path()
)::Dict{String, Any}
    if isfile(config_path)
        try
            return TOML.parsefile(config_path)
        catch err
            @warn "Failed to parse configuration file at $(config_path): $(err). Using defaults."
            if occursin("snowcrab_tesselated", lowercase(config_path))
                return get_snowcrab_tesselated_configuration()
            elseif occursin("snowcrab", lowercase(config_path))
                return get_snowcrab_configuration()
            else
                return get_default_configuration()
            end
        end
    else
        if occursin("snowcrab_tesselated", lowercase(config_path))
            return get_snowcrab_tesselated_configuration()
        elseif occursin("snowcrab", lowercase(config_path))
            return get_snowcrab_configuration()
        else
            return get_default_configuration()
        end
    end
end

"""
    save_configuration(
        config_dict::AbstractDict,
        config_path::AbstractString = joinpath("inputs", "ParticleTracking.config")
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
    config_path::AbstractString = joinpath("inputs", "ParticleTracking.config")
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
            "nz" => 10,
            "resolution_scale" => 1.0,
            "vertical_stretching_mode" => "tanh",
            "vertical_grid_file" => "inputs/scotian_shelf_vertical_grid.csv"
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
            "diffusivity_v" => 1e-4,
            "release_depth_mode" => "bottom",
            "bottom_release_offset" => [0.5, 3.0],
            "enable_initial_ascent" => true,
            "ascent_speed" => 0.010,
            "ascent_target_depth" => -10.0
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
            "export_parquet" => false,
            "enable_checkpoint" => true,
            "checkpoint_prefix" => "checkpoint_ParticleTracking",
            "checkpoint_schedule_seconds" => 21600.0,
            "checkpoint_dir" => "outputs/checkpoints",
            "checkpoint_cleanup" => true
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => false,
            "fallback_to_cpu" => true
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => true,
            "animate_hydro" => false,
            "anim_variable" => "dashboard",
            "anim_fps" => 10,
            "anim_format" => "mp4",
            "anim_depth" => -2.5,
            "anim_overlay_particles" => false,
            "title" => "Regional Marine Lagrangian Particle Tracking & Dispersion"
        ),
        "paths" => Dict{String, Any}(
            "output_dir" => "outputs",
            "input_dir" => "inputs",
            "seed" => 42
        )
    )
end

"""
    get_snowcrab_configuration() -> Dict{String, Any}

Generate a configuration dictionary calibrated specifically for snow crab (*Chionoecetes opilio*)
larval transport modeling across the Scotian Shelf and Northwest Atlantic.

# Biophysical Calibration & Parameters
- **Domain**: Scotian Shelf & continental slope (`lon`: [-68.0, -57.0] °E, `lat`: [42.0, 47.5] °N, `z`: [-3500.0, 0.0] m).
- **Grid**: High-resolution 100 × 100 × 20 cells.
- **Data**: Real bathymetry and atmospheric forcing (`data_mode`: "real").
- **Tides**: M2 (0.25 m/s) + S2 (0.11 m/s) tidal current forcing with quadratic bottom drag.
- **Biology**: 500 larvae released from commercial nursery depths (≥ 100 m) in benthic boundary layer ([0.5, 3.0] m off bed).
- **Behaviors**: Active vertical ascent (10 mm/s to -10 m), stage-dependent DVM, and degree-day thermal molting
  (base temperature \$T_{base} = -1.5\$ °C; stage thresholds: 65, 130, 200 DD).
- **Duration**: 60.0 days pelagic larval duration (PLD) advected at 300 s time steps.
- **Storage**: DuckDB target configured as `outputs/snowcrab_tracking.duckdb`.

# Outputs
- `Dict{String, Any}`: Nested dictionary containing snow crab calibrated parameters.
"""
function get_snowcrab_configuration()::Dict{String, Any}
    return Dict{String, Any}(
        "domain" => Dict{String, Any}(
            "lon_min" => -68.0,
            "lon_max" => -57.0,
            "lat_min" => 42.0,
            "lat_max" => 47.5,
            "z_min" => -3500.0,
            "z_max" => 0.0,
            "buffer_km" => 100.0
        ),
        "grid" => Dict{String, Any}(
            "nx" => 345,
            "ny" => 245,
            "nz" => 20,
            "resolution_scale" => 2.5,
            "vertical_stretching_mode" => "tanh",
            "vertical_grid_file" => "inputs/scotian_shelf_vertical_grid.csv"
        ),
        "data" => Dict{String, Any}(
            "data_mode" => "real",
            "bathy_dataset_id" => "etopo180",
            "wind_dataset_id" => "erdBSwinds1day",
            "wind_time_iso" => "2020-06-01T00:00:00Z",
            "inshore_depth" => -100.0,
            "shelf_slope" => 500.0
        ),
        "tides" => Dict{String, Any}(
            "enable_tides" => true,
            "constituents" => ["M2", "S2"],
            "tidal_u_amp" => 0.25,
            "tidal_v_amp" => 0.12,
            "s2_u_amp" => 0.11,
            "s2_v_amp" => 0.05,
            "tidal_period" => 44712.0,
            "tidal_phase" => 0.0
        ),
        "climate" => Dict{String, Any}(
            "scenario" => "historical",
            "projection_year" => 2020,
            "baseline_year" => 2020,
            "horizon_year" => 2050
        ),
        "hydrodynamics" => Dict{String, Any}(
            "sim_duration_hours" => 120.0,
            "sim_dt_seconds" => 120.0,
            "adaptive_cfl" => true,
            "target_cfl" => 0.2,
            "divergence_velocity_limit" => 20.0,
            "coriolis_latitude" => 44.5,
            "surface_heat_flux" => 50.0
        ),
        "biology" => Dict{String, Any}(
            "n_particles" => 500,
            "track_duration_days" => 60.0,
            "track_dt_seconds" => 300.0,
            "min_seabed_depth" => 100.0,
            "buffer_km" => 100.0,
            "diffusivity_h" => 10.0,
            "diffusivity_v" => 1e-4,
            "release_depth_mode" => "bottom",
            "bottom_release_offset" => [0.5, 3.0],
            "enable_initial_ascent" => true,
            "ascent_speed" => 0.010,
            "ascent_target_depth" => -10.0
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
            "t_base" => -1.5,
            "dd_zoea1_to_zoea2" => 65.0,
            "dd_zoea2_to_megalopa" => 130.0,
            "dd_megalopa_to_settle" => 200.0,
            "mortality_base" => 0.02,
            "mortality_thermal_threshold" => 7.0,
            "mortality_thermal_sensitivity" => 0.35,
            "mortality_cold_threshold" => -1.5,
            "mortality_cold_sensitivity" => 0.02,
            "settlement_min_depth" => -250.0,
            "settlement_max_depth" => -50.0,
            "settlement_max_temp" => 6.0
        ),
        "storage" => Dict{String, Any}(
            "enable_duckdb" => true,
            "duckdb_path" => "outputs/snowcrab_tracking.duckdb",
            "export_parquet" => false,
            "enable_checkpoint" => true,
            "checkpoint_prefix" => "checkpoint_snowcrab",
            "checkpoint_schedule_seconds" => 21600.0,
            "checkpoint_dir" => "outputs/checkpoints",
            "checkpoint_cleanup" => true
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => true,
            "fallback_to_cpu" => true
        ),
        "atmosphere" => Dict{String, Any}(
            "source" => "era5",
            "drag_formulation" => "garratt_1977",
            "bulk_heat_flux" => true,
            "climatology" => false
        ),
        "boundaries" => Dict{String, Any}(
            "ocean_boundary_source" => "glorys12v1",
            "obc_type" => "flather_chapman",
            "sponge_layer_width" => 0.25,
            "sponge_timescale" => 3600.0
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
    get_snowcrab_tesselated_configuration() -> Dict{String, Any}

Generate a configuration dictionary calibrated specifically for snow crab (*Chionoecetes opilio*)
larval transport with multi-resolution, depth-stratified Voronoi tessellation analysis.

# Biophysical Calibration & Parameters
- Inherits all baseline Scotian Shelf biophysical calibrations from `get_snowcrab_configuration()`.
- Incorporates a dedicated `[tessellation]` block configuring \$N = 5000\$ Voronoi units:
  - **Core Nursery Strata** (50 to 350 m depth): Probability weight \$p = 0.8\$, minimum
    separation distance 1.5 km (resolving complex banks, gullies, and nursery depressions).
  - **Shallow Coastal Strata** (0 to 50 m depth): Probability weight \$p = 0.1\$, minimum
    separation distance 5.0 km.
  - **Deep Slope & Basin Strata** (> 350 m depth): Probability weight \$p = 0.1\$, minimum
    separation distance 10.0 km.
- Persistence configured to `outputs/snowcrab_tesselated.duckdb`.

# Outputs
- `Dict{String, Any}`: Nested dictionary containing snow crab tessellated parameter settings.
"""
function get_snowcrab_tesselated_configuration()::Dict{String, Any}
    cfg = get_snowcrab_configuration()
    cfg["storage"]["duckdb_path"] = "outputs/snowcrab_tesselated.duckdb"
    cfg["storage"]["checkpoint_prefix"] = "checkpoint_snowcrab_tesselated"
    cfg["visualization"]["title"] = "Scotian Shelf Snow Crab Multi-Resolution Voronoi Dispersal & Connectivity"
    cfg["tessellation"] = Dict{String, Any}(
        "enable_voronoi" => true,
        "n_units" => 5000,
        "prob_core" => 0.8,
        "prob_shallow" => 0.1,
        "prob_deep" => 0.1,
        "min_res_core_km" => 1.5,
        "min_res_shallow_km" => 5.0,
        "min_res_deep_km" => 10.0,
        "core_depth_min" => -350.0,
        "core_depth_max" => -50.0
    )
    return cfg
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
    heat_flux = Float64(get_val("hydrodynamics", "surface_heat_flux", 50.0))

    n_parts   = Int(get_val("biology", "n_particles", 100))
    track_days = Float64(get_val("biology", "track_duration_days", 5.0))
    track_dt  = Float64(get_val("biology", "track_dt_seconds", 300.0))
    min_depth = Float64(get_val("biology", "min_seabed_depth", 100.0))
    diff_h    = Float64(get_val("biology", "diffusivity_h", 10.0))
    diff_v    = Float64(get_val("biology", "diffusivity_v", 1e-4))

    rel_mode_str = String(get_val("biology", "release_depth_mode", "bottom"))
    rel_mode = Symbol(lowercase(rel_mode_str))
    raw_offset = get_val("biology", "bottom_release_offset", [0.5, 3.0])
    b_offset = if raw_offset isa Vector && length(raw_offset) >= 2
        (Float64(raw_offset[1]), Float64(raw_offset[2]))
    elseif raw_offset isa Tuple && length(raw_offset) >= 2
        (Float64(raw_offset[1]), Float64(raw_offset[2]))
    else
        (0.5, 3.0)
    end
    init_ascent = Bool(get_val("biology", "enable_initial_ascent", true))
    asc_spd     = Float64(get_val("biology", "ascent_speed", 0.010))
    asc_target  = Float64(get_val("biology", "ascent_target_depth", -10.0))

    enable_dvm = Bool(get_val("dvm", "enable_dvm", true))
    enable_molting = Bool(get_val("molting_and_settlement", "enable_molting", true))

    use_gpu      = Bool(get_val("hardware", "use_gpu", false))
    fallback_cpu = Bool(get_val("hardware", "fallback_to_cpu", true))
    interactive  = Bool(get_val("visualization", "interactive_map", true))
    anim_hydro   = Bool(get_val("visualization", "animate_hydro", false))
    anim_var     = Symbol(lowercase(String(get_val("visualization", "anim_variable", "dashboard"))))
    anim_fps     = Int(get_val("visualization", "anim_fps", 10))
    anim_fmt     = String(get_val("visualization", "anim_format", "mp4"))
    anim_depth   = Float64(get_val("visualization", "anim_depth", -2.5))
    anim_overlay = Bool(get_val("visualization", "anim_overlay_particles", false))

    enable_duckdb = Bool(get_val("storage", "enable_duckdb", true))
    duckdb_path   = String(get_val("storage", "duckdb_path", "outputs/particle_tracking.duckdb"))

    output_dir = String(get_val("paths", "output_dir", "outputs"))
    input_dir  = String(get_val("paths", "input_dir", "inputs"))
    seed       = Int(get_val("paths", "seed", 42))

    hydro_file = String(get_val("hydrodynamics", "hydro_model_file",
                                get_val("storage", "output_filename", "")))
    run_id_val = String(get_val("storage", "run_id", ""))

    enable_cp = Bool(get_val("storage", "enable_checkpoint",
                             get_val("hydrodynamics", "enable_checkpoint", true)))
    cp_pfx_raw = String(get_val("storage", "checkpoint_prefix",
                                get_val("hydrodynamics", "checkpoint_prefix", "")))
    cp_sched  = Float64(get_val("storage", "checkpoint_schedule_seconds",
                                get_val("hydrodynamics", "checkpoint_schedule_seconds", 0.0)))
    cp_dir    = String(get_val("storage", "checkpoint_dir",
                               get_val("paths", "checkpoint_dir", "")))
    cp_clean  = Bool(get_val("storage", "checkpoint_cleanup", true))
    auto_res  = Bool(get_val("hydrodynamics", "auto_restart", true))

    res_scale = Float64(get_val("grid", "resolution_scale", 1.0))
    v_mode_str = String(get_val("grid", "vertical_stretching_mode", "tanh"))
    v_mode = Symbol(lowercase(v_mode_str))
    v_file = String(get_val("grid", "vertical_grid_file",
                            joinpath("inputs", "scotian_shelf_vertical_grid.csv")))
    atmo_src = Symbol(lowercase(String(get_val("atmosphere", "source", "era5"))))
    obc_src = Symbol(lowercase(String(get_val("boundaries", "ocean_boundary_source", "glorys12v1"))))
    obc_tp = Symbol(lowercase(String(get_val("boundaries", "obc_type", "flather_chapman"))))

    # Voronoi tessellation parameters
    enable_voronoi = Bool(get_val("tessellation", "enable_voronoi", false))
    voronoi_n_units = Int(get_val("tessellation", "n_units", 5000))
    voronoi_p_core = Float64(get_val("tessellation", "prob_core", 0.8))
    voronoi_p_shallow = Float64(get_val("tessellation", "prob_shallow", 0.1))
    voronoi_p_deep = Float64(get_val("tessellation", "prob_deep", 0.1))
    voronoi_min_core = Float64(get_val("tessellation", "min_res_core_km", 1.5))
    voronoi_min_shallow = Float64(get_val("tessellation", "min_res_shallow_km", 5.0))
    voronoi_min_deep = Float64(get_val("tessellation", "min_res_deep_km", 10.0))

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
        surface_heat_flux = heat_flux,
        n_particles = n_parts,
        track_duration = track_days * 86400.0,
        track_dt = track_dt,
        diffusivity_h = diff_h,
        diffusivity_v = diff_v,
        enable_dvm = enable_dvm,
        enable_molting = enable_molting,
        min_seabed_depth = min_depth,
        buffer_km = buffer_km,
        release_depth_mode = rel_mode,
        bottom_release_offset = b_offset,
        enable_initial_ascent = init_ascent,
        ascent_speed = asc_spd,
        ascent_target_depth = asc_target,
        use_gpu = use_gpu,
        fallback_to_cpu = fallback_cpu,
        interactive_map = interactive,
        enable_duckdb = enable_duckdb,
        duckdb_path = duckdb_path,
        output_dir = output_dir,
        input_dir = input_dir,
        seed = seed,
        hydro_model_file = hydro_file,
        enable_checkpoint = enable_cp,
        checkpoint_prefix = cp_pfx_raw,
        checkpoint_schedule = cp_sched,
        checkpoint_dir = cp_dir,
        checkpoint_cleanup = cp_clean,
        auto_restart = auto_res,
        run_id = run_id_val,
        vertical_stretching_mode = v_mode,
        vertical_grid_file = v_file,
        resolution_scale = res_scale,
        atmospheric_source = atmo_src,
        ocean_boundary_source = obc_src,
        obc_type = obc_tp,
        enable_voronoi = enable_voronoi,
        voronoi_n_units = voronoi_n_units,
        voronoi_prob_core = voronoi_p_core,
        voronoi_prob_shallow = voronoi_p_shallow,
        voronoi_prob_deep = voronoi_p_deep,
        voronoi_min_res_core_km = voronoi_min_core,
        voronoi_min_res_shallow_km = voronoi_min_shallow,
        voronoi_min_res_deep_km = voronoi_min_deep,
        animate_hydro = anim_hydro,
        anim_variable = anim_var,
        anim_fps = anim_fps,
        anim_format = anim_fmt,
        anim_depth = anim_depth,
        anim_overlay_particles = anim_overlay,
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
            "nz" => opts.grid_size[3],
            "resolution_scale" => opts.resolution_scale,
            "vertical_stretching_mode" => string(opts.vertical_stretching_mode),
            "vertical_grid_file" => opts.vertical_grid_file
        ),
        "atmosphere" => Dict{String, Any}(
            "source" => string(opts.atmospheric_source)
        ),
        "boundaries" => Dict{String, Any}(
            "ocean_boundary_source" => string(opts.ocean_boundary_source),
            "obc_type" => string(opts.obc_type)
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
            "target_cfl" => opts.target_cfl,
            "surface_heat_flux" => opts.surface_heat_flux,
            "hydro_model_file" => opts.hydro_model_file
        ),
        "biology" => Dict{String, Any}(
            "n_particles" => opts.n_particles,
            "track_duration_days" => opts.track_duration / 86400.0,
            "track_dt_seconds" => opts.track_dt,
            "min_seabed_depth" => opts.min_seabed_depth,
            "buffer_km" => opts.buffer_km,
            "diffusivity_h" => opts.diffusivity_h,
            "diffusivity_v" => opts.diffusivity_v,
            "release_depth_mode" => string(opts.release_depth_mode),
            "bottom_release_offset" => [opts.bottom_release_offset[1], opts.bottom_release_offset[2]],
            "enable_initial_ascent" => opts.enable_initial_ascent,
            "ascent_speed" => opts.ascent_speed,
            "ascent_target_depth" => opts.ascent_target_depth
        ),
        "dvm" => Dict{String, Any}(
            "enable_dvm" => opts.enable_dvm
        ),
        "molting_and_settlement" => Dict{String, Any}(
            "enable_molting" => opts.enable_molting
        ),
        "storage" => Dict{String, Any}(
            "enable_duckdb" => opts.enable_duckdb,
            "duckdb_path" => opts.duckdb_path,
            "run_id" => opts.run_id,
            "enable_checkpoint" => opts.enable_checkpoint,
            "checkpoint_prefix" => opts.checkpoint_prefix,
            "checkpoint_schedule_seconds" => opts.checkpoint_schedule,
            "checkpoint_dir" => opts.checkpoint_dir,
            "checkpoint_cleanup" => opts.checkpoint_cleanup
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => opts.use_gpu,
            "fallback_to_cpu" => opts.fallback_to_cpu
        ),
        "tessellation" => Dict{String, Any}(
            "enable_voronoi" => opts.enable_voronoi,
            "n_units" => opts.voronoi_n_units,
            "prob_core" => opts.voronoi_prob_core,
            "prob_shallow" => opts.voronoi_prob_shallow,
            "prob_deep" => opts.voronoi_prob_deep,
            "min_res_core_km" => opts.voronoi_min_res_core_km,
            "min_res_shallow_km" => opts.voronoi_min_res_shallow_km,
            "min_res_deep_km" => opts.voronoi_min_res_deep_km
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => opts.interactive_map,
            "animate_hydro" => opts.animate_hydro,
            "anim_variable" => string(opts.anim_variable),
            "anim_fps" => opts.anim_fps,
            "anim_format" => opts.anim_format,
            "anim_depth" => opts.anim_depth,
            "anim_overlay_particles" => opts.anim_overlay_particles
        ),
        "paths" => Dict{String, Any}(
            "output_dir" => opts.output_dir,
            "input_dir" => opts.input_dir,
            "seed" => opts.seed
        )
    )
end

"""
    SnowCrabRunOptions(; kwargs...) -> HydrodynamicOptions

Construct a `HydrodynamicOptions` instance pre-configured with calibrated physical
and biophysical parameters for snow crab (*Chionoecetes opilio*) larval transport
across the Scotian Shelf and Northwest Atlantic continental slope.

# Mathematical & Biophysical Formulations
Snow crab larvae hatch from benthic nurseries (\$z_{\\text{bed}} \\le -100\\text{ m}\$)
and actively ascend to the epipelagic zone:
- **Larval Cohort**: 500 larvae tracked for 60.0 days pelagic larval duration (PLD)
  at \$\\Delta t = 300\\text{ s}\$.
- **Spatial Domain**: Scotian Shelf & slope (\$\\lambda \\in [-68.0, -57.0]^\\circ\\text{E}\$,
  \$\\phi \\in [42.0, 47.5]^\\circ\\text{N}\$, \$z \\in [-3500.0, 0.0]\\text{ m}\$).
- **Grid Resolution**: \$100 \\times 100 \\times 20\$ grid cells with \$100\\text{ km}\$
  CFA buffer.
- **Benthic Release & Active Ascent**: Initial placement in near-bottom boundary layer
  ([0.5, 3.0] m off bed) with directed upward swimming (\$w = 0.010\\text{ m/s}\$)
  relaxing toward \$-10.0\\text{ m}\$, smoothly transitioning to stage-dependent DVM.
- **Degree-Day Thermal Molting**: Base temperature \$T_{\\text{base}} = -1.5^\\circ\\text{C}\$
  with cumulative degree-day thresholds (Zoea I \$\\to\$ Zoea II: 65 DD,
  Zoea II \$\\to\$ Megalopa: 130 DD, Megalopa \$\\to\$ Benthic Settlement: 200 DD).
- **Environmental Forcing & Storage**: Real bathymetry and winds, \$M_2 + S_2\$ tidal
  harmonics, summer heat flux (\$50\\text{ W/m}^2\$), and DuckDB storage
  (`outputs/snowcrab_tracking.duckdb`).

# Inputs
- `kwargs...`: Optional keyword overrides for any field of `HydrodynamicOptions`
  (e.g., `n_particles = 200`, `ascent_speed = 0.015`, `hydro_model_file = "..."`,
  `track_only = true`, `run_id = "cohort_2020"`).

# Outputs
- `HydrodynamicOptions`: Validated runtime options instance with snow crab defaults.

# References
- Sainte-Marie, G., & Sainte-Marie, B. (1999). Hatching and larval release in the
  snow crab, *Chionoecetes opilio*. *Journal of Crustacean Biology*, 19(4), 743-754.
- Lovrich, G. A., & Sainte-Marie, B. (1997). Cannibalism in the snow crab:
  implications for larval recruitment. *Marine Ecology Progress Series*, 148, 85-99.
- Kuhn, P. S., & Choi, J. S. (2011). Influence of temperature on snow crab larval
  development. *ICES Journal of Marine Science*, 68(8), 1673-1681.
"""
function SnowCrabRunOptions(; kwargs...)::HydrodynamicOptions
    cfg = get_snowcrab_configuration()
    return configuration_to_options(cfg; kwargs...)
end

"""
    SnowCrabTesselatedRunOptions(; kwargs...) -> HydrodynamicOptions

Construct a `HydrodynamicOptions` instance pre-configured with calibrated physical
and biophysical parameters for snow crab (*Chionoecetes opilio*) and multi-resolution
depth-stratified Voronoi tessellation analysis.

# Inputs
- `kwargs...`: Optional keyword overrides for any field of `HydrodynamicOptions`.

# Outputs
- `HydrodynamicOptions`: Validated runtime options instance with snow crab tessellated defaults.
"""
function SnowCrabTesselatedRunOptions(; kwargs...)::HydrodynamicOptions
    cfg = get_snowcrab_tesselated_configuration()
    return configuration_to_options(cfg; kwargs...)
end

# ─────────────────────────────────────────────────────────────────────────────
# Formally Decoupled Configuration Architecture
# ─────────────────────────────────────────────────────────────────────────────

"""
    HydrodynamicConfig

Physical parameters specifying Eulerian ocean circulation, lateral boundary
conditions, atmospheric fluxes, and numerical time integration.
"""
struct HydrodynamicConfig
    domain_lon               :: Tuple{Float64, Float64}
    domain_lat               :: Tuple{Float64, Float64}
    domain_z                 :: Tuple{Float64, Float64}
    grid_size                :: Tuple{Int, Int, Int}
    vertical_stretching_mode :: Symbol
    vertical_grid_file       :: String
    resolution_scale         :: Float64
    bathymetry_source        :: Symbol
    bathy_dataset_id         :: String
    inshore_depth            :: Float64
    shelf_slope              :: Float64
    atmospheric_source       :: Symbol
    drag_formulation         :: Symbol
    bulk_heat_flux           :: Bool
    climatology              :: Bool
    mhw_temp_anomaly         :: Float64
    ocean_boundary_source    :: Symbol
    obc_type                 :: Symbol
    sponge_layer_width       :: Float64
    sponge_timescale         :: Float64
    enable_tides             :: Bool
    tides_source             :: Symbol
    tidal_constituents       :: Vector{Symbol}
    tidal_u_amp              :: Float64
    tidal_v_amp              :: Float64
    cd_drag                  :: Float64
    bottom_drag              :: Float64
    bbl_mixing_closure       :: Symbol
    sim_duration_seconds     :: Float64
    sim_dt_seconds           :: Float64
    adaptive_cfl             :: Bool
    target_cfl               :: Float64
    target_wave_cfl          :: Float64
    max_dt_seconds           :: Float64
    min_dt_seconds           :: Float64
    coriolis_latitude        :: Float64
    divergence_limit         :: Float64
    output_dir               :: String
    output_filename          :: String
    output_schedule_seconds  :: Float64
    enable_checkpoint        :: Bool
    checkpoint_prefix        :: String
    checkpoint_schedule      :: Float64
    checkpoint_dir           :: String
    checkpoint_cleanup       :: Bool
    auto_restart             :: Bool
    use_gpu                  :: Bool
    fallback_to_cpu          :: Bool
end

"""
    LarvalDispersalConfig

Biological parameters specifying Lagrangian particle tracking, active vertical
locomotion, degree-day molting, thermal mortality, and benthic settlement.
"""
struct LarvalDispersalConfig
    n_particles              :: Int
    track_duration_seconds   :: Float64
    track_dt_seconds         :: Float64
    diffusivity_h            :: Float64
    diffusivity_v            :: Float64
    min_seabed_depth         :: Float64
    buffer_km                :: Float64
    release_depth_mode       :: Symbol
    bottom_release_offset    :: Tuple{Float64, Float64}
    enable_initial_ascent    :: Bool
    ascent_speed             :: Float64
    ascent_target_depth      :: Float64
    enable_dvm               :: Bool
    megalopa_day_depth       :: Float64
    megalopa_night_depth     :: Float64
    zoea2_depth_factor       :: Float64
    megalopa_swim_factor     :: Float64
    enable_molting           :: Bool
    t_base                   :: Float64
    dd_zoea1_to_zoea2        :: Float64
    dd_zoea2_to_megalopa     :: Float64
    dd_megalopa_to_settle    :: Float64
    mortality_base           :: Float64
    mortality_thermal_thresh :: Float64
    mortality_thermal_sens   :: Float64
    mortality_cold_thresh    :: Float64
    mortality_cold_sens      :: Float64
    settlement_min_depth     :: Float64
    settlement_max_depth     :: Float64
    settlement_max_temp      :: Float64
    enable_duckdb            :: Bool
    duckdb_path              :: String
    run_id                   :: String
    seed                     :: Int
end

"""
    CoupledSimulationConfig

Unified configuration aggregating physical `HydrodynamicConfig` and biological
`LarvalDispersalConfig` for integrated end-to-end simulations.
"""
struct CoupledSimulationConfig
    hydro  :: HydrodynamicConfig
    larval :: LarvalDispersalConfig
end

"""
    to_hydrodynamic_config(opts::HydrodynamicOptions; kwargs...) -> HydrodynamicConfig

Extract and convert a `HydrodynamicOptions` instance into a pure `HydrodynamicConfig`.
"""
function to_hydrodynamic_config(opts::HydrodynamicOptions; kwargs...)::HydrodynamicConfig
    d = Dict{Symbol, Any}(
        :domain_lon               => opts.domain_lon,
        :domain_lat               => opts.domain_lat,
        :domain_z                 => opts.domain_z,
        :grid_size                => opts.grid_size,
        :vertical_stretching_mode => opts.vertical_stretching_mode,
        :vertical_grid_file       => opts.vertical_grid_file,
        :resolution_scale         => opts.resolution_scale,
        :bathymetry_source        => opts.data_mode == :real ? :gebco : :synthetic,
        :bathy_dataset_id         => "etopo180",
        :inshore_depth            => -20.0,
        :shelf_slope              => 500.0,
        :atmospheric_source       => opts.atmospheric_source,
        :drag_formulation         => :garratt_1977,
        :bulk_heat_flux           => true,
        :climatology              => false,
        :mhw_temp_anomaly         => opts.scenario == :mhw ? 3.5 : 0.0,
        :ocean_boundary_source    => opts.ocean_boundary_source,
        :obc_type                 => opts.obc_type,
        :sponge_layer_width       => 0.25,
        :sponge_timescale         => 3600.0,
        :enable_tides             => opts.enable_tides,
        :tides_source             => :tpxo9_atlas,
        :tidal_constituents       => [:M2, :S2],
        :tidal_u_amp              => opts.tidal_u_amp,
        :tidal_v_amp              => opts.tidal_v_amp,
        :cd_drag                  => 0.0025,
        :bottom_drag              => 0.0001,
        :bbl_mixing_closure       => :tke,
        :sim_duration_seconds     => opts.sim_duration,
        :sim_dt_seconds           => opts.sim_dt,
        :adaptive_cfl             => opts.adaptive_cfl,
        :target_cfl               => opts.target_cfl,
        :target_wave_cfl          => 0.20,
        :max_dt_seconds           => 90.0,
        :min_dt_seconds           => 2.0,
        :coriolis_latitude        => 44.5,
        :divergence_limit         => 20.0,
        :output_dir               => opts.output_dir,
        :output_filename          => isempty(opts.hydro_model_file) ?
                                     "hydrodynamics_output.jld2" :
                                     basename(opts.hydro_model_file),
        :output_schedule_seconds  => 21600.0,
        :enable_checkpoint        => opts.enable_checkpoint,
        :checkpoint_prefix        => opts.checkpoint_prefix,
        :checkpoint_schedule      => opts.checkpoint_schedule,
        :checkpoint_dir           => opts.checkpoint_dir,
        :checkpoint_cleanup       => opts.checkpoint_cleanup,
        :auto_restart             => opts.auto_restart,
        :use_gpu                  => opts.use_gpu,
        :fallback_to_cpu          => opts.fallback_to_cpu
    )
    for (k, v) in kwargs
        d[k] = v
    end

    return HydrodynamicConfig(
        d[:domain_lon], d[:domain_lat], d[:domain_z], d[:grid_size],
        Symbol(d[:vertical_stretching_mode]), String(d[:vertical_grid_file]),
        Float64(d[:resolution_scale]), Symbol(d[:bathymetry_source]),
        String(d[:bathy_dataset_id]), Float64(d[:inshore_depth]),
        Float64(d[:shelf_slope]), Symbol(d[:atmospheric_source]),
        Symbol(d[:drag_formulation]), Bool(d[:bulk_heat_flux]),
        Bool(d[:climatology]), Float64(d[:mhw_temp_anomaly]),
        Symbol(d[:ocean_boundary_source]), Symbol(d[:obc_type]),
        Float64(d[:sponge_layer_width]), Float64(d[:sponge_timescale]),
        Bool(d[:enable_tides]), Symbol(d[:tides_source]),
        Vector{Symbol}(d[:tidal_constituents]), Float64(d[:tidal_u_amp]),
        Float64(d[:tidal_v_amp]), Float64(d[:cd_drag]),
        Float64(d[:bottom_drag]), Symbol(d[:bbl_mixing_closure]),
        Float64(d[:sim_duration_seconds]), Float64(d[:sim_dt_seconds]),
        Bool(d[:adaptive_cfl]), Float64(d[:target_cfl]),
        Float64(d[:target_wave_cfl]), Float64(d[:max_dt_seconds]),
        Float64(d[:min_dt_seconds]), Float64(d[:coriolis_latitude]),
        Float64(d[:divergence_limit]), String(d[:output_dir]),
        String(d[:output_filename]), Float64(d[:output_schedule_seconds]),
        Bool(d[:enable_checkpoint]), String(d[:checkpoint_prefix]),
        Float64(d[:checkpoint_schedule]), String(d[:checkpoint_dir]),
        Bool(d[:checkpoint_cleanup]), Bool(d[:auto_restart]),
        Bool(d[:use_gpu]), Bool(d[:fallback_to_cpu])
    )
end

"""
    to_larval_config(opts::HydrodynamicOptions; kwargs...) -> LarvalDispersalConfig

Extract and convert a `HydrodynamicOptions` instance into a pure `LarvalDispersalConfig`.
"""
function to_larval_config(opts::HydrodynamicOptions; kwargs...)::LarvalDispersalConfig
    d = Dict{Symbol, Any}(
        :n_particles              => opts.n_particles,
        :track_duration_seconds   => opts.track_duration,
        :track_dt_seconds         => opts.track_dt,
        :diffusivity_h            => opts.diffusivity_h,
        :diffusivity_v            => opts.diffusivity_v,
        :min_seabed_depth         => opts.min_seabed_depth,
        :buffer_km                => opts.buffer_km,
        :release_depth_mode       => opts.release_depth_mode,
        :bottom_release_offset    => opts.bottom_release_offset,
        :enable_initial_ascent    => opts.enable_initial_ascent,
        :ascent_speed             => opts.ascent_speed,
        :ascent_target_depth      => opts.ascent_target_depth,
        :enable_dvm               => opts.enable_dvm,
        :megalopa_day_depth       => -120.0,
        :megalopa_night_depth     => -60.0,
        :zoea2_depth_factor       => 1.2,
        :megalopa_swim_factor     => 1.5,
        :enable_molting           => opts.enable_molting,
        :t_base                   => -1.5,
        :dd_zoea1_to_zoea2        => 65.0,
        :dd_zoea2_to_megalopa     => 130.0,
        :dd_megalopa_to_settle    => 200.0,
        :mortality_base           => 0.02,
        :mortality_thermal_thresh => 7.0,
        :mortality_thermal_sens   => 0.35,
        :mortality_cold_thresh    => -1.5,
        :mortality_cold_sens      => 0.02,
        :settlement_min_depth     => -250.0,
        :settlement_max_depth     => -50.0,
        :settlement_max_temp      => 6.0,
        :enable_duckdb            => opts.enable_duckdb,
        :duckdb_path              => opts.duckdb_path,
        :run_id                   => opts.run_id,
        :seed                     => opts.seed
    )
    for (k, v) in kwargs
        d[k] = v
    end

    return LarvalDispersalConfig(
        Int(d[:n_particles]), Float64(d[:track_duration_seconds]),
        Float64(d[:track_dt_seconds]), Float64(d[:diffusivity_h]),
        Float64(d[:diffusivity_v]), Float64(d[:min_seabed_depth]),
        Float64(d[:buffer_km]), Symbol(d[:release_depth_mode]),
        (Float64(d[:bottom_release_offset][1]), Float64(d[:bottom_release_offset][2])),
        Bool(d[:enable_initial_ascent]), Float64(d[:ascent_speed]),
        Float64(d[:ascent_target_depth]), Bool(d[:enable_dvm]),
        Float64(d[:megalopa_day_depth]), Float64(d[:megalopa_night_depth]),
        Float64(d[:zoea2_depth_factor]), Float64(d[:megalopa_swim_factor]),
        Bool(d[:enable_molting]), Float64(d[:t_base]),
        Float64(d[:dd_zoea1_to_zoea2]), Float64(d[:dd_zoea2_to_megalopa]),
        Float64(d[:dd_megalopa_to_settle]), Float64(d[:mortality_base]),
        Float64(d[:mortality_thermal_thresh]), Float64(d[:mortality_thermal_sens]),
        Float64(d[:mortality_cold_thresh]), Float64(d[:mortality_cold_sens]),
        Float64(d[:settlement_min_depth]), Float64(d[:settlement_max_depth]),
        Float64(d[:settlement_max_temp]), Bool(d[:enable_duckdb]),
        String(d[:duckdb_path]), String(d[:run_id]), Int(d[:seed])
    )
end


