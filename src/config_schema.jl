"""
    config_schema.jl

Type-safe configuration schema using Configurations.jl for validation and parsing.
Replaces manual TOML parsing with structured, validated configuration objects.
"""

module ConfigSchema

using Configurations
using Configurations: @option

# Domain configuration
@option struct DomainConfig
    lon_min::Float64 = -68.0
    lon_max::Float64 = -57.0
    lat_min::Float64 = 42.0
    lat_max::Float64 = 47.5
    z_min::Float64 = -5000.0
    z_max::Float64 = 0.0
    buffer_km::Float64 = 50.0
end

# Grid configuration
@option struct GridConfig
    nx::Int = 345
    ny::Int = 245
    nz::Int = 20
    resolution_scale::Float64 = 2.5
    vertical_stretching_mode::String = "tanh"
    vertical_grid_file::String = "inputs/scotian_shelf_vertical_grid.csv"
end

# Data configuration
@option struct DataConfig
    data_mode::String = "real"
    bathy_dataset_id::String = "nceiEtopo2022"
    wind_dataset_id::String = "erdBSwinds1day"
    hydrography_source::String = "woa23"
    wind_time_iso::String = "2020-06-01T00:00:00Z"
    inshore_depth::Float64 = -20.0
    shelf_slope::Float64 = 500.0
end

# Bathymetry configuration
@option struct BathymetryConfig
    source::String = "numerical_earth"
    provider::String = "etopo2022"
    resolution_arcsec::Int = 15
    fallback_dataset_id::String = "etopo180"
    inshore_depth::Float64 = -20.0
    shelf_slope::Float64 = 500.0
end

# Boundaries configuration
@option struct BoundariesConfig
    method::String = "relaxation"
    parent_ocean::String = "glorys_climatology"
    ocean_boundary_source::String = "glorys12v1"
    climatology::Bool = true
    obc_type::String = "flather_chapman"
    active_boundaries::Vector{String} = String["east", "south", "west"]
    sponge_layer_width::Float64 = 0.25
    sponge_timescale_seconds::Float64 = 3600.0
    u_inflow::Float64 = -0.15
    v_inflow::Float64 = 0.05
end

# Atmosphere configuration
@option struct AtmosphereConfig
    source::String = "era5_climatology"
    climatology::Bool = true
    variables::Vector{String} = String["u10", "v10", "t2m", "ssrd", "strd"]
    drag_formulation::String = "garratt_1977"
    bulk_heat_flux::Bool = true
end

# Tides configuration
@option struct TidesConfig
    enable_tides::Bool = true
    tides_source::String = "tpxo9_atlas"
    constituents::Vector{String} = String["M2", "S2"]
    tidal_u_amp::Float64 = 0.25
    tidal_v_amp::Float64 = 0.12
    s2_u_amp::Float64 = 0.11
    s2_v_amp::Float64 = 0.05
    tidal_period::Float64 = 44712.0
    tidal_phase::Float64 = 0.0
end

# Bottom boundary layer configuration
@option struct BottomBoundaryLayerConfig
    cd_drag::Float64 = 0.0025
    bottom_drag::Float64 = 0.0001
    bbl_mixing_closure::String = "tke"
end

# Hydrodynamics configuration
@option struct HydrodynamicsConfig
    sim_duration_hours::Float64 = 17520.0
    sim_dt_seconds::Float64 = 180.0
    adaptive_cfl::Bool = true
    target_cfl::Float64 = 0.20
    target_wave_cfl::Float64 = 0.20
    max_dt_seconds::Float64 = 300.0
    min_dt_seconds::Float64 = 2.0
    coriolis_latitude::Float64 = 44.5
    divergence_velocity_limit::Float64 = 20.0
    surface_heat_flux::Float64 = 50.0
    turbulence_closure::String = "nemotke"
    enable_sea_ice::Bool = true
    enable_oxygen::Bool = true
    output_schedule_seconds::Float64 = 21600.0
    hydro_model_file::String = "hydrodynamics_climatology_2yr.jld2"
    enable_checkpoint::Bool = true
    checkpoint_schedule::Float64 = 3600.0
    checkpoint_dir::String = ""
    checkpoint_cleanup::Bool = true
    auto_restart::Bool = true
end

# Tessellation configuration
@option struct TessellationConfig
    enable_voronoi::Bool = true
    n_units::Int = 5000
    bathymetry_source::String = "numerical_earth"
    core_depth_min::Float64 = 50.0
    core_depth_max::Float64 = 350.0
    core_prob::Float64 = 0.8
    prob_core::Float64 = 0.8
    core_min_res_km::Float64 = 1.5
    min_res_core_km::Float64 = 1.5
    shallow_depth_min::Float64 = 0.0
    shallow_depth_max::Float64 = 50.0
    shallow_prob::Float64 = 0.1
    prob_shallow::Float64 = 0.1
    shallow_min_res_km::Float64 = 5.0
    min_res_shallow_km::Float64 = 5.0
    deep_depth_min::Float64 = 350.0
    deep_depth_max::Float64 = 6000.0
    deep_prob::Float64 = 0.1
    prob_deep::Float64 = 0.1
    deep_min_res_km::Float64 = 10.0
    min_res_deep_km::Float64 = 10.0
    slope_weighting::Bool = true
    slope_factor::Float64 = 15.0
end

# Biology configuration
@option struct BiologyConfig
    n_particles::Int = 500
    track_duration_days::Float64 = 60.0
    track_dt_seconds::Float64 = 300.0
    min_seabed_depth::Float64 = 50.0
    diffusivity_h::Float64 = 10.0
    diffusivity_v::Float64 = 1e-4
    release_depth_mode::String = "bottom"
    bottom_release_offset::Vector{Float64} = Float64[0.5, 3.0]
    enable_initial_ascent::Bool = true
    ascent_speed::Float64 = 0.010
    ascent_target_depth::Float64 = -10.0
    buffer_km::Float64 = 50.0
end

# DVM configuration
@option struct DVMConfig
    enable_dvm::Bool = true
    zoea1_day_depth::Float64 = -50.0
    zoea1_night_depth::Float64 = -10.0
    zoea2_day_depth::Float64 = -55.0
    zoea2_night_depth::Float64 = -8.0
    megalopa_day_depth::Float64 = -120.0
    megalopa_night_depth::Float64 = -60.0
    zoea2_depth_factor::Float64 = 1.0
    megalopa_swim_factor::Float64 = 1.0
end

# Molting and settlement configuration
@option struct MoltingAndSettlementConfig
    enable_molting::Bool = true
    t_base::Float64 = -1.5
    dd_zoea1_to_zoea2::Float64 = 65.0
    dd_zoea2_to_megalopa::Float64 = 130.0
    dd_megalopa_to_settle::Float64 = 200.0
    settlement_min_depth::Float64 = -350.0
    settlement_max_depth::Float64 = -50.0
    settlement_max_temp::Float64 = 6.0
    mortality_base::Float64 = 0.01
    mortality_thermal_threshold::Float64 = 15.0
    mortality_thermal_sensitivity::Float64 = 0.02
    mortality_cold_threshold::Float64 = 0.0
    mortality_cold_sensitivity::Float64 = 0.01
end

# Climate configuration
@option struct ClimateConfig
    scenario::String = "climatology"
    projection_year::Int = 2022
    baseline_year::Int = 2020
    horizon_year::Int = 2050
end

# Storage configuration
@option struct StorageConfig
    output_dir::String = "outputs"
    output_filename::String = "hydrodynamics_climatology_2yr.jld2"
    output_schedule_seconds::Float64 = 21600.0
    checkpoint_fields::Vector{String} = String["u", "v", "w", "T", "S", "η"]
    enable_duckdb::Bool = true
    duckdb_path::String = "outputs/snowcrab_tracking.duckdb"
    enable_checkpoint::Bool = true
    checkpoint_prefix::String = "checkpoint_hydrodynamics_climatology_2yr"
    checkpoint_schedule::Float64 = 21600.0
    checkpoint_schedule_seconds::Float64 = 21600.0
    checkpoint_dir::String = "outputs/checkpoints"
    cleanup_checkpoints::Bool = true
    auto_restart::Bool = true
end

# Hardware configuration
@option struct HardwareConfig
    use_gpu::Bool = true
    fallback_to_cpu::Bool = true
end

# Visualization configuration
@option struct VisualizationConfig
    interactive_map::Bool = true
    title::String = "Regional Marine Lagrangian Particle Tracking & Dispersion"
end

# Paths configuration
@option struct PathsConfig
    output_dir::String = "outputs"
    input_dir::String = "inputs"
    seed::Int = 42
end

# Metadata configuration
@option struct MetadataConfig
    name::String = "Scotian Shelf 2-Year Climatological Ocean Circulation"
    description::String = "Equilibrium multi-year simulation driven by repeating annual climatologies and tides"
    author::String = "Jae Choi"
    version::String = "2.0.0"
end

# Main configuration struct combining all sections
@option struct ParticleTrackingConfig
    metadata::MetadataConfig = MetadataConfig()
    data::DataConfig = DataConfig()
    climate::ClimateConfig = ClimateConfig()
    domain::DomainConfig = DomainConfig()
    grid::GridConfig = GridConfig()
    bathymetry::BathymetryConfig = BathymetryConfig()
    atmosphere::AtmosphereConfig = AtmosphereConfig()
    boundaries::BoundariesConfig = BoundariesConfig()
    tides::TidesConfig = TidesConfig()
    bottom_boundary_layer::BottomBoundaryLayerConfig = BottomBoundaryLayerConfig()
    hydrodynamics::HydrodynamicsConfig = HydrodynamicsConfig()
    tessellation::TessellationConfig = TessellationConfig()
    biology::BiologyConfig = BiologyConfig()
    dvm::DVMConfig = DVMConfig()
    molting_and_settlement::MoltingAndSettlementConfig = MoltingAndSettlementConfig()
    storage::StorageConfig = StorageConfig()
    hardware::HardwareConfig = HardwareConfig()
    visualization::VisualizationConfig = VisualizationConfig()
    paths::PathsConfig = PathsConfig()
end

# Convenience functions
"""
    load_config(config_path::AbstractString) -> ParticleTrackingConfig

Load and validate configuration from a TOML file using Configurations.jl.
"""
function load_config(config_path::AbstractString)::ParticleTrackingConfig
    if isfile(config_path)
        dict = Configurations.from_toml(config_path)
        return Configurations.parse_from_dict(ParticleTrackingConfig, dict)
    else
        @warn "Config file not found at $(config_path), using defaults"
        return ParticleTrackingConfig()
    end
end

"""
    save_config(config::ParticleTrackingConfig, config_path::AbstractString) -> String

Save configuration to a TOML file.
"""
function save_config(config::ParticleTrackingConfig, config_path::AbstractString)::String
    out_dir = dirname(abspath(config_path))
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    Configurations.to_toml(config, config_path)
    return config_path
end

# Conversion to HydrodynamicOptions (for backward compatibility)
# This will be implemented once we update the HydrodynamicOptions structure

export ParticleTrackingConfig, load_config, save_config

end # module ConfigSchema