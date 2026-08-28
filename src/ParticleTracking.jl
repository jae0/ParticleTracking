"""
    ParticleTracking

A Julia module coupling ocean hydrodynamic circulation modeling
(via Oceananigans.jl) with Lagrangian particle tracking for planktonic
and semiplanktonic snow crab (*Chionoecetes opilio*) larvae on coastal shelves.

# Mathematical & Architectural Foundation
The module provides:
1. Data generation, ingestion, and inspection of bathymetry and surface forcing (`synthetic_data.jl`).
2. Discretization of coastal shelf geometry with immersed boundary topography (`grid_bathymetry.jl`).
3. Hydrostatic Boussinesq primitive equation modeling with Coriolis, buoyancy, and wind stress (`hydrodynamic_model.jl`).
4. Simulation orchestration, time-stepping diagnostics, and field output writers (`simulation.jl`).
5. Larval stage parameterizations, Diel Vertical Migration (DVM), and stochastic transport (`larval_behavior.jl`).

# References
- Ramadhan, A., Marshall, J., Hill, C., Campin, J. M., Bischoff, T., & Wagner, G. L. (2020).
  Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs.
  *Journal of Open Source Software*, 5(53), 2018. DOI: 10.21105/joss.02018
- Incze, L. S., Armstrong, D. A., & Smith, S. L. (1987). Abundance of filter-feeding
  and pelagic stages of crab larvae in the southeastern Bering Sea.
  *Marine Biology*, 95(2), 195-200. DOI: 10.1007/BF00409006
- Lovrich, G. A., Sainte-Marie, B., & Smith, B. D. (1995). Depth distribution and
  seasonal movements of *Chionoecetes opilio* in Baie Sainte-Marguerite.
  *Canadian Journal of Fisheries and Aquatic Sciences*, 52(4), 903-913.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). Growth, developmental stages, and
  vertical distribution of snow crab larvae (*Chionoecetes opilio*).
  *Canadian Journal of Fisheries and Aquatic Sciences*, 56(11), 2181-2193.
"""
module ParticleTracking

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using CairoMakie
using NCDatasets
using Downloads
using Random
using DuckDB
using DataFrames
using DBInterface
using Dates
using Statistics
using LinearAlgebra
using TOML

# Sub-components
include("configuration.jl")
include("open_data.jl")
include("synthetic_data.jl")
include("architecture.jl")
include("grid_bathymetry.jl")
include("hydrodynamic_model.jl")
include("tides.jl")
include("climate_scenarios.jl")
include("simulation.jl")
include("larval_behavior.jl")
include("empirical_analysis.jl")
include("storage_duckdb.jl")
include("visualization.jl")

# Exported APIs
export
    # Centralized configuration management
    HydrodynamicOptions,
    load_configuration,
    save_configuration,
    configuration_to_options,
    options_to_configuration,
    get_default_configuration,
    find_default_config_path,

    # Architecture and device resolution
    resolve_architecture,

    # Open real-world data and regridding
    fetch_open_bathymetry,
    fetch_open_surface_winds,
    wind_speed_to_kinematic_stress,
    regrid_2d_field,

    # Synthetic data generation and inspection
    download_sample_data,
    inspect_netcdf,
    generate_synthetic_bathymetry,
    generate_synthetic_forcing,

    # Grid, bathymetry & coastline geometry
    REGIONAL_COASTLINE,
    build_shelf_grid,
    load_bathymetry_from_netcdf,
    get_bathymetry_interpolator,
    load_coastline_polygons,
    save_coastline_polygons,
    is_point_on_land,
    is_marine_water,
    extract_marine_cells,
    sample_marine_coordinates,
    buffer_distance_to_degrees,
    expand_domain_with_buffer,
    get_strata_buffered_envelope,
    build_immersed_grid,
    build_immersed_grid_from_real_data,

    # Hydrodynamic model
    build_hydrodynamic_model,
    set_initial_stratification!,
    set_initial_conditions!,

    # Tidal forcing & harmonic synthesis
    get_tidal_frequency,
    build_tidal_body_forcing,
    tidal_velocity_vector,
    simpson_hunter_parameter,

    # Climate scenarios and thermal biology
    get_climate_scenario_deltas,
    apply_climate_scenario!,
    temperature_dependent_pld,
    larval_thermal_mortality_rate,

    # Simulation and execution
    compute_advective_cfl,
    setup_hydrodynamic_simulation,
    run_hydrodynamic_simulation!,
    create_flow_interpolator_from_jld2,

    # Larval behavior and particle tracking
    diel_vertical_migration_velocity,
    superpose_tidal_velocity,
    update_larval_stage,
    evaluate_settlement_suitability,
    larval_transport_step,
    initialize_larval_particles,
    track_larval_cohort,

    # Empirical movement, recruitment & connectivity analysis
    point_in_polygon,
    intersect_polygon_with_coastline,
    load_cfa_polygons,
    estimate_empirical_movement,
    compute_gridded_recruitment_metrics,
    compute_gridded_thermal_metrics,
    compute_empirical_connectivity,
    export_larval_dispersal_netcdf,
    export_larval_dispersal_jld2,

    # DuckDB analytical storage, scenario management & model averaging
    open_duckdb_storage,
    close_duckdb_storage,
    initialize_duckdb_schema!,
    save_simulation_run!,
    load_run_configuration,
    list_simulation_runs,
    load_trajectories_df,
    load_trajectories_namedtuple,
    load_all_scenario_trajectories,
    save_hydrodynamic_field!,
    load_hydrodynamic_field,
    load_gridded_dispersal,
    load_connectivity_matrix,
    compare_scenarios,
    compute_ensemble_model_average,
    export_duckdb_to_parquet,

    # Visualization, interactive maps & scenario comparison
    plot_particle_trajectories,
    plot_vertical_migration_profiles,
    compare_scenario_dispersal,
    plot_larval_dispersal_density,
    plot_empirical_movement_field,
    plot_connectivity_matrix,
    plot_thermal_exposure_map,
    plot_recruitment_summary,
    extract_hydrodynamic_dataset,
    plot_hydrodynamic_advection,
    plot_hydrodynamic_tracers,
    export_interactive_tracks_html,
    plot_interactive_trajectories_map

end # module ParticleTracking
