"""
    configuration.jl

Centralized configuration manager for `ParticleTracking.jl`.
Provides robust parsing, serialization, validation, and conversion between
structured TOML configuration files (such as `inputs/ParticleTracking.toml`),
nested dictionaries, metadata records, and runtime `HydrodynamicOptions` structs.
"""

using TOML

"""
    find_default_config_path() -> String

Locate the configuration file in the project workspace, or return `""` when there is none.

Only project-root locations are searched. A configuration is an input the operator supplies;
it is never discovered inside a previous run's containerised `inputs/` directory, and it is
never guessed at. In particular there is no fallback to `configs/default.toml`: a run that
silently adopts a different scenario than the one intended produces output that looks
complete and answers a question nobody asked. Callers must fail loudly instead, which is what
`load_configuration` does with an empty path.

To run, pass one explicitly:

    julia --project=. ParticleTrackingRun.jl --config=configs/snowcrab.toml

or place `ParticleTracking.toml` at the project root.
"""
function find_default_config_path()::String
    for p in ("ParticleTracking.toml", "ParticleTracking.config")
        isfile(p) && return p
    end
    ""
end

const NO_CONFIG_ADVICE = """
No configuration file was found.

A configuration is required, not optional. It names the datasets to ingest, the domain and
grid, the forcing sources, and the biology parameters, and those choices are what a run is.
Without one the package would have to invent all of it, and the resulting output would carry
no statement of what it actually simulated -- it would be indistinguishable from a run that
was configured deliberately.

Available configurations in this project:
  configs/default.toml     small, fast baseline
  configs/opendata.toml   keyless observed products, moderate grid
  configs/snowcrab.toml    Scotian Shelf, high resolution, with tidal forcing

Run with one of:

    julia --project=. ParticleTrackingRun.jl --config=configs/snowcrab.toml

or copy one to ParticleTracking.toml at the project root, which the loader picks up
automatically. For a study of your own, copy the closest config and edit the sections the
run depends on: [data] for datasets, [domain] and [grid] for the box and resolution,
[tides] for tidal forcing, and [biology] for larval parameters."""

"""
    study_domain_ranges(cfg) -> NamedTuple

Geographic bounds of the *study domain* — the computational/analysis grid — read from a
parsed configuration.

Authoritative source is the TOML `[domain]` section (`lon_min`, `lon_max`, `lat_min`,
`lat_max`). No domain is compiled into the package: the configuration file is the single
source of truth, so a run against a different shelf only requires editing TOML.

# Inputs
- `cfg`: Parsed configuration `AbstractDict` (e.g. from [`load_configuration`](@ref)).

# Outputs
- `NamedTuple`: `(lon = (lon_min, lon_max), lat = (lat_min, lat_max))`, in degrees.
"""
function study_domain_ranges(cfg)::NamedTuple
    d = get(cfg, "domain", Dict{String, Any}())
    return (
        lon = (Float64(get(d, "lon_min", -68.0)), Float64(get(d, "lon_max", -57.0))),
        lat = (Float64(get(d, "lat_min", 42.0)), Float64(get(d, "lat_max", 47.5)))
    )
end

"""
    embedding_domain_ranges(cfg) -> NamedTuple

Geographic bounds of the *embedding domain* — the larger region from which open boundary
and atmospheric forcing data are extracted — read from a parsed configuration.

Authoritative source is the TOML `[boundaries]` section (`embedding_lon_min`,
`embedding_lon_max`, `embedding_lat_min`, `embedding_lat_max`). The study domain nests
inside this region.

# Inputs
- `cfg`: Parsed configuration `AbstractDict`.

# Outputs
- `NamedTuple`: `(lon = (lon_min, lon_max), lat = (lat_min, lat_max))`, in degrees.
"""
function embedding_domain_ranges(cfg)::NamedTuple
    b = get(cfg, "boundaries", Dict{String, Any}())
    return (
        lon = (Float64(get(b, "embedding_lon_min", -71.0)), Float64(get(b, "embedding_lon_max", -53.0))),
        lat = (Float64(get(b, "embedding_lat_min", 40.0)), Float64(get(b, "embedding_lat_max", 48.5)))
    )
end

"""
    resolve_config_name(config_file::AbstractString = "") -> String

Extract the base configuration name from a config file path, or return
the default configuration name (`"ParticleTracking"`) if unspecified or empty.

# Inputs
- `config_file::AbstractString`: Path to a `.toml` file or empty string.

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
- `bathy_source::String`: Which bathymetry product to ingest, read verbatim from
  `[data] bathy_source` (for example `"etopo2022"`, `"synthetic"`). There is no separate
  real/synthetic switch: acquisition always reads the named product, and falls back to analytic
  fields with a warning if a download fails.
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
- `fallback_to_cpu::Bool`: Whether a GPU request that cannot be honoured may degrade to CPU.
  Defaults to `false`, so the run fails loudly instead of quietly running somewhere other than
  where `use_gpu` said it would.
- `interactive_map::Bool`: Whether to export interactive HTML5 Leaflet map.
- `enable_duckdb::Bool`: Whether to persist simulation data to DuckDB.
- `duckdb_path::String`: DuckDB file path.
- `config_file::String`: Source configuration file path.
- `bathy_source::String`: Which bathymetry product to ingest; the run always reads the product
  named here rather than choosing between a real and a synthetic mode.
- `inshore_depth::Float64`: Water depth at the landward edge of the domain (m, negative).
- `shelf_slope::Float64`: Total seabed rise across the shelf (m).
- `output_schedule_seconds::Float64`: Field-output cadence. Independent of
  `checkpoint_schedule`, which is the state-checkpoint cadence; the two may legitimately differ.
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
  See the field comments for the full accepted set.
- `obc_type::Symbol`: Lateral open boundary formulation (:flather_chapman, :radiation, :clamped).
- `enable_voronoi::Bool`: Whether to compute multi-resolution Voronoi tessellation analysis.
- `voronoi_n_units::Int`: Number of Voronoi units/centroids to generate across strata.
- `voronoi_prob_core::Float64`: Sampling probability weight for core depth stratum (50 to 350 m).
- `voronoi_prob_shallow::Float64`: Sampling probability weight for shallow stratum (0 to 50 m).
- `voronoi_prob_deep::Float64`: Sampling probability weight for deep stratum (> 350 m).
- `voronoi_min_res_core_km::Float64`: Minimum Poisson-disc separation distance in core zone (km).
- `voronoi_min_res_shallow_km::Float64`: Minimum separation distance in shallow zone (km).
- `voronoi_min_res_deep_km::Float64`: Minimum separation distance in deep zone (km).
- `voronoi_core_depth_min`/`_max`, `voronoi_shallow_depth_min`/`_max`,
  `voronoi_deep_depth_min`/`_max`::Float64: The edges of the three depth strata, as positive
  depths below the surface. These decide which seabed area each stratum samples, so they change
  the connectivity result, not just its resolution. The sign is applied when the bands reach the
  tessellator, which bins on signed `z`.
- `voronoi_slope_weighting::Bool`: Whether minimum spacing is tightened on steep bathymetry.
- `voronoi_slope_factor::Float64`: How aggressively spacing tightens with slope; only used when
  `voronoi_slope_weighting` is true.
"""
struct HydrodynamicOptions
    domain_lon               :: Tuple{Float64, Float64}
    domain_lat               :: Tuple{Float64, Float64}
    domain_z                 :: Tuple{Float64, Float64}
    grid_size                :: Tuple{Int, Int, Int}
    bathy_source             :: String
    coastline_source         :: String
    wind_source              :: String
    wind_time_iso            :: String
    inshore_depth            :: Float64
    shelf_slope              :: Float64
    enable_tides             :: Bool
    tides_source             :: String
    tidal_constituents       :: Vector{Symbol}
    tidal_u_amp              :: Float64
    tidal_v_amp              :: Float64
    s2_u_amp                 :: Float64
    s2_v_amp                 :: Float64
    scenario                 :: Symbol
    projection_year          :: Int
    sim_dt                   :: Float64
    sim_duration             :: Float64
    min_dt_seconds           :: Float64
    adaptive_cfl             :: Bool
    target_cfl               :: Float64
      surface_heat_flux        :: Float64
      bulk_heat_flux           :: Bool
    max_flow_snapshots       :: Int
    allow_analytical_fallback :: Bool
    max_current_speed        :: Union{Nothing, Float64}
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
    output_schedule_seconds  :: Float64
    checkpoint_dir           :: String
    checkpoint_cleanup       :: Bool
    auto_restart             :: Bool
    run_id                   :: String
    vertical_stretching_mode :: Symbol
    nz_above :: Int
    vertical_break_depth :: Float64
    vertical_grid_file       :: String
    vertical_depths          :: Vector{Float64}
    resolution_scale         :: Float64
    atmospheric_source :: Symbol
    ocean_boundary_source :: Symbol
    # The larger region the open-boundary data is cut from. These are separate from
    # `domain_lon`/`domain_lat` because the boundary needs the water just *outside* the study
    # area; a subset cut to the study domain would put the interpolation's edge exactly on the
    # sponge. They are fields rather than re-read from the raw config dict because a value
    # fetched by key from somewhere downstream is a value nothing checks.
    embedding_lon :: Tuple{Float64, Float64}
    embedding_lat :: Tuple{Float64, Float64}
    # Where each piece of the model's outside world comes from.
    #
    # The model is a window cut out of the ocean, so three of its edges are not real coastline
    # but imaginary lines. Along those lines the model has to be told what the water is doing
    # just outside: how warm, how salty, and which way it is moving. That is the
    # "lateral boundary" below. If it is left unset, the model still runs, but water can drift
    # in or out of the open edges unchecked, and results that depend on the flow across an
    # edge are then unreliable.
    #
    # The values are stored as Symbols, so a misspelling stays visible in a comparison rather
    # than quietly becoming an empty field. They are read straight from the TOML with no
    # normalisation, so the accepted spellings are exactly as listed.
    #
    #   ocean_boundary_source  -- what sets the imaginary edges
    #     "synthetic"    Nothing is downloaded. The open edges are simply not constrained, and
    #                    the model behaves as if they are out of sight. This is a real choice
    #                    for a run that only studies the interior; it is not a stand-in for
    #                    data, and it writes no file.
    #     "glorys12v1"   A global real-world reconstruction, GLORYS12V1, at 1/12 degree (about
    #                    9 km). Comes from Copernicus Marine, which is free but does require
    #                    registering an account and putting a username and key in a file
    #                    called ~/.copernicusmarine/credentials. Without that file the download
    #                    is rejected.
    #     "hycom"        HYCOM GOFS at the same 1/12 degree (about 9 km), and unlike the
    #                    options below it also carries how fast the water is moving, which is
    #                    the part that matters most for an open edge. No account, no key. Read
    #                    the note on `fetch_hycom_boundary` before relying on it: its data is
    #                    currently not downloadable even though the catalogue works.
    #     "woa23"        A monthly average of the whole ocean, free and no account. But it
    #                    holds only temperature and salt, not movement, so it fixes what the
    #                    water is like but not how much of it flows through the edge.
    #
    #   atmospheric_source  -- what drives the sea surface from the sky
    #     "open_meteo"        Real weather, served for free by an open website. It also
    #                         reports air temperature, pressure, humidity, cloud cover and
    #                         sunlight, which is what the model needs to work out how much
    #                         heat the sea gains or loses each hour. Currently the only choice
    #                         that makes the heat flux change over the run. It gives one
    #                         representative spot rather than a picture of the sky, so the same
    #                         weather is used across the whole domain.
    #     "era5_climatology"  A fixed typical month, repeating. No download. Use it only to
    #                         debug; the sea will warm monotonically and never vary.
    #     "synthetic"         A single fixed number. No download, no variation. Diagnostics only.
    #
    #   tides_source  -- what sets the rise and fall of the tides
    #     "tpxo9_atlas"   Measured tide heights and currents, free and no account. Supplies the
    #                     real size of the M2 and S2 tides, which is what makes the spring-neap
    #                     cycle (big tides every two weeks) come out right.
    #     "tpxo9"         The same product, named more briefly.
    #     "analytic"      Tide shapes made up from a single number you supply. No download, and
    #                     no detail across space, so tides are the same everywhere in the
    #                     domain. Fine for checking the machinery, not for a real tide.
    #
    #   bathy_source  -- where the sea floor shape comes from
    #     "etopo2022"  A global map of the sea floor, free and no account, detailed to about
    #                  450 m. "etopo" and "open" are other spellings of the same thing.
    #
    #   hydrography_source  -- what the water starts out like at t = 0
    #     "woa23"      Average temperature and salt for a chosen month, free and no account.
    #     "synthetic"  A made-up smooth profile, no download. For testing.
    #     "glorys12v1" A real reconstruction; needs the Copernicus account described above.
    #
    #   coastline_source  -- where the outline of the land comes from
    #     "natural_earth_10m"  Very detailed, public domain, no account. "50m" and "110m" are
    #                         coarser outlines for when the detail is not needed.
    #
    # None of these falls back quietly. If a value is not recognised the run stops and says
    # what it could have used instead. That is deliberate: a run that asked for one product
    # and quietly received a different one would have different temperature or flow at its
    # edges, and nothing in the output would reveal it.
    # Source for the 3-D hydrographic (T, S) initial state. This is deliberately
    # separate from `bathy_source`/`ocean_boundary_source`, which govern bathymetry,
    # wind and lateral boundaries: the initial water column is its own choice.
    #   :synthetic — analytic profile from `set_initial_stratification!`
    #   :woa23     — World Ocean Atlas 2023 climatology (keyless, NOAA NCEI)
    #   :glorys12v1— GLORYS12V1 reanalysis (requires copernicusmarine credentials)
    hydrography_source       :: Symbol
    hydrography_month        :: Int      # WOA23 month selector (0 = annual climatology)
    obc_type                 :: Symbol
    enable_voronoi           :: Bool
    voronoi_n_units          :: Int
    voronoi_prob_core        :: Float64
    voronoi_prob_shallow     :: Float64
    voronoi_prob_deep        :: Float64
    voronoi_min_res_core_km  :: Float64
    voronoi_min_res_shallow_km :: Float64
    voronoi_min_res_deep_km  :: Float64
    voronoi_core_depth_min   :: Float64
    voronoi_core_depth_max   :: Float64
    voronoi_shallow_depth_min :: Float64
    voronoi_shallow_depth_max :: Float64
    voronoi_deep_depth_min   :: Float64
    voronoi_deep_depth_max   :: Float64
    voronoi_slope_weighting  :: Bool
    voronoi_slope_factor     :: Float64
    animate_hydro            :: Bool
    anim_variable            :: Symbol
    anim_fps                 :: Int
    anim_format              :: String
    anim_depth               :: Float64
    anim_overlay_particles   :: Bool
    anim_output_path         :: String
    # Lateral boundary relaxation parameters (Price & Aumont 2011)
    boundary_method          :: Symbol   # :relaxation | :none
    sponge_layer_width       :: Float64  # sponge buffer width (degrees)
    sponge_timescale         :: Float64  # relaxation timescale (s)
    tidal_relaxation_timescale :: Float64  # relaxation toward the prescribed tide (s)
    u_inflow                 :: Float64  # reference zonal inflow velocity (m/s)
    v_inflow                 :: Float64  # reference meridional inflow velocity (m/s)
    # Time integration ceiling
    max_dt                   :: Float64  # maximum allowable time step (s)
    # Stochastic biology parameters
    settlement_stochastic    :: Bool     # whether settlement uses probabilistic HSI draw
    cv_molt                  :: Float64  # coefficient of variation for molt thresholds (default 0.25)
    cv_mortality             :: Float64  # coefficient of variation for mortality (default 0.25)
    cv_settlement            :: Float64  # coefficient of variation for settlement HSI (default 0.25)
end

function HydrodynamicOptions(;
    domain_lon            :: Tuple{Real, Real} = (-71.0, -53.0),
    domain_lat            :: Tuple{Real, Real} = (40.0, 48.5),
    domain_z              :: Tuple{Real, Real} = (-3500.0, 0.0),
    grid_size             :: Tuple{Int, Int, Int} = (50, 50, 10),
    bathy_source           :: AbstractString = "etopo2022",
    coastline_source       :: AbstractString = "natural_earth_10m",
    wind_source            :: AbstractString = "open_meteo",
    wind_time_iso          :: AbstractString = "2020-06-01T00:00:00Z",
    inshore_depth          :: Real = -20.0,
    shelf_slope            :: Real = 500.0,
    enable_tides          :: Bool = true,
    tides_source          :: AbstractString = "tpxo9",
    tidal_constituents    :: AbstractVector{Symbol} = [:M2, :S2],
    tidal_u_amp           :: Real = 0.25,
    tidal_v_amp           :: Real = 0.12,
    s2_u_amp              :: Real = 0.11,
    s2_v_amp              :: Real = 0.05,
    scenario              :: Symbol = :ssp245,
    projection_year       :: Int = 2050,
    sim_dt                :: Real = 120.0,
    sim_duration          :: Real = 43200.0,
    min_dt_seconds        :: Real = 2.0,
    adaptive_cfl          :: Bool = true,
    target_cfl            :: Real = 0.2,
      surface_heat_flux     :: Real = 50.0,
      bulk_heat_flux        :: Bool = true,
    max_flow_snapshots    :: Int = 400,
    allow_analytical_fallback :: Bool = false,
    max_current_speed     :: Union{Nothing, Real} = 3.0,
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
    output_schedule_seconds :: Real = 21600.0,
    checkpoint_dir        :: AbstractString = "",
    checkpoint_cleanup    :: Bool = true,
    auto_restart             :: Bool = true,
    run_id                   :: AbstractString = "",
    vertical_stretching_mode :: Symbol = :oceananigans,
    nz_above :: Int = 0,
    vertical_break_depth :: Real = -400.0,
    vertical_grid_file       :: AbstractString = "",
    vertical_depths          :: AbstractVector{<:Real} = Float64[],
    resolution_scale         :: Real = 1.0,
    atmospheric_source       :: Symbol = :era5,
    ocean_boundary_source    :: Symbol = :glorys12v1,
    embedding_lon            :: Tuple{Float64, Float64} = (-71.0, -53.0),
    embedding_lat            :: Tuple{Float64, Float64} = (40.0, 48.5),
    hydrography_source       :: Symbol = :synthetic,
    hydrography_month        :: Int     = 0,
    obc_type                 :: Symbol = :flather_chapman,
    enable_voronoi           :: Bool = false,
    voronoi_n_units          :: Int = 5000,
    voronoi_prob_core        :: Real = 0.8,
    voronoi_prob_shallow     :: Real = 0.1,
    voronoi_prob_deep        :: Real = 0.1,
    voronoi_min_res_core_km  :: Real = 1.5,
    voronoi_min_res_shallow_km :: Real = 5.0,
    voronoi_min_res_deep_km  :: Real = 10.0,
    voronoi_core_depth_min   :: Real = 50.0,
    voronoi_core_depth_max   :: Real = 350.0,
    voronoi_shallow_depth_min :: Real = 0.0,
    voronoi_shallow_depth_max :: Real = 50.0,
    voronoi_deep_depth_min   :: Real = 350.0,
    voronoi_deep_depth_max   :: Real = 6000.0,
    voronoi_slope_weighting  :: Bool = true,
    voronoi_slope_factor     :: Real = 15.0,
    animate_hydro            :: Bool = false,
    anim_variable            :: Symbol = :dashboard,
    anim_fps                 :: Int = 10,
    anim_format              :: AbstractString = "mp4",
    anim_depth               :: Real = -2.5,
    anim_overlay_particles   :: Bool = false,
    anim_output_path         :: AbstractString = "",
    boundary_method          :: Symbol = :relaxation,
    sponge_layer_width       :: Real = 0.35,
    sponge_timescale         :: Real = 3600.0,
    tidal_relaxation_timescale :: Real = 3600.0,
    u_inflow                 :: Real = -0.15,
    v_inflow                 :: Real = 0.05,
    max_dt                   :: Real = 600.0,
    settlement_stochastic    :: Bool = true,
    cv_molt                  :: Real = 0.25,
    cv_mortality             :: Real = 0.25,
    cv_settlement            :: Real = 0.25
)
    resolved_cp_prefix = if !isempty(strip(checkpoint_prefix)) && checkpoint_prefix != "checkpoint"
        String(checkpoint_prefix)
    else
        cfg_name = resolve_config_name(config_file)
        "checkpoint_$(cfg_name)"
    end

    # The open-boundary embedding region arrives as the `embedding_lon`/`embedding_lat`
    # keywords rather than being re-derived from a config dict here: this constructor takes
    # individual keywords and has no parsed `cfg` to read. `configuration_to_options`, which
    # does have the parsed TOML, is what supplies them -- so the value has exactly one path
    # from the file into the struct.

    return HydrodynamicOptions(
        (Float64(domain_lon[1]), Float64(domain_lon[2])),
        (Float64(domain_lat[1]), Float64(domain_lat[2])),
        (Float64(domain_z[1]), Float64(domain_z[2])),
        grid_size,
        bathy_source,
        coastline_source,
        wind_source,
        wind_time_iso,
        Float64(inshore_depth),
        Float64(shelf_slope),
        enable_tides,
        tides_source,
        tidal_constituents,
        Float64(tidal_u_amp),
        Float64(tidal_v_amp),
        Float64(s2_u_amp),
        Float64(s2_v_amp),
        scenario,
        projection_year,
        Float64(sim_dt),
        Float64(sim_duration),
        Float64(min_dt_seconds),
        adaptive_cfl,
        Float64(target_cfl),
  Float64(surface_heat_flux),
  bulk_heat_flux,
  Int(max_flow_snapshots),
        allow_analytical_fallback,
        (max_current_speed isa Real && max_current_speed > 0) ? Float64(max_current_speed) : nothing,
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
        Float64(output_schedule_seconds),
        String(checkpoint_dir),
        checkpoint_cleanup,
        auto_restart,
        String(run_id),
        vertical_stretching_mode,
        Int(nz_above),
        Float64(vertical_break_depth),
        String(vertical_grid_file),
        Float64[Float64(d) for d in vertical_depths],
        Float64(resolution_scale),
        atmospheric_source,
        ocean_boundary_source,
        embedding_lon,
        embedding_lat,
        hydrography_source,
        hydrography_month,
        obc_type,
        enable_voronoi,
        voronoi_n_units,
        Float64(voronoi_prob_core),
        Float64(voronoi_prob_shallow),
        Float64(voronoi_prob_deep),
        Float64(voronoi_min_res_core_km),
        Float64(voronoi_min_res_shallow_km),
        Float64(voronoi_min_res_deep_km),
        Float64(voronoi_core_depth_min),
        Float64(voronoi_core_depth_max),
        Float64(voronoi_shallow_depth_min),
        Float64(voronoi_shallow_depth_max),
        Float64(voronoi_deep_depth_min),
        Float64(voronoi_deep_depth_max),
        voronoi_slope_weighting,
        Float64(voronoi_slope_factor),
        animate_hydro,
        anim_variable,
        anim_fps,
        String(anim_format),
        Float64(anim_depth),
        anim_overlay_particles,
        String(anim_output_path),
        boundary_method,
        Float64(sponge_layer_width),
        Float64(sponge_timescale),
        Float64(tidal_relaxation_timescale),
        Float64(u_inflow),
        Float64(v_inflow),
        Float64(max_dt),
        settlement_stochastic,
        Float64(cv_molt),
        Float64(cv_mortality),
        Float64(cv_settlement)
    )
end

"""
    load_configuration(config_path::AbstractString = find_default_config_path()) -> Dict{String, Any}

Read and parse a centralized `ParticleTracking.toml` file into a nested Julia Dictionary.
A missing or unparseable file is an error, not a silent fall back to defaults.

# Inputs
- `config_path::AbstractString`: Path to the `.toml` (TOML format) file.

# Outputs
- `Dict{String, Any}`: Nested dictionary containing all sectioned parameter settings.
"""
function load_configuration(
    config_path::AbstractString = find_default_config_path()
)::Dict{String, Any}
    # A run must never proceed on defaults it did not ask for: a silently substituted
    # configuration produces output that looks real and is not. `get_default_configuration()`
    # remains available for callers that genuinely want the defaults.
    isempty(config_path) && error(strip(NO_CONFIG_ADVICE, '\n'))
    isfile(config_path) || error(
        "Configuration file not found: $(config_path). Pass --config=<path>."
    )
    try
        return TOML.parsefile(config_path)
    catch err
        error("Failed to parse configuration file at $(config_path): $(err)")
    end
end

"""
    save_configuration(
        config_dict::AbstractDict,
        config_path::AbstractString = "ParticleTracking.toml"
    ) -> String

Serialize a nested configuration dictionary to a centralized `.toml` file in TOML format.

# Inputs
- `config_dict::AbstractDict`: Dictionary of configuration sections and key-values.
- `config_path::AbstractString`: Target destination file path.

# Outputs
- `String`: Path to the written configuration file.
"""
function save_configuration(
    config_dict::AbstractDict,
    config_path::AbstractString = "ParticleTracking.toml"
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
            "vertical_stretching_mode" => "oceananigans",
            "vertical_grid_file" => ""
        ),
        "data" => Dict{String, Any}(
            "bathy_source" => "etopo2022",
            "coastline_source" => "natural_earth_10m",
            "wind_source" => "open_meteo",
            "hydrography_source" => "woa23",
            "wind_time_iso" => "2023-06-01T00:00:00Z"
        ),
        "bathymetry" => Dict{String, Any}(
            "source" => "numerical_earth",
            "resolution_arcsec" => 15,
            "inshore_depth" => -20.0,
            "shelf_slope" => 500.0
        ),
        "boundaries" => Dict{String, Any}(
            "method" => "relaxation",
            "ocean_boundary_source" => "glorys12v1",
            "obc_type" => "flather_chapman",
            "sponge_layer_width" => 0.35,
            "sponge_timescale" => 3600.0,
            "u_inflow" => -0.15,
            "v_inflow" => 0.05
        ),
        "atmosphere" => Dict{String, Any}(
            "source" => "era5",
            "drag_formulation" => "garratt_1977",
            "bulk_heat_flux" => true
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
            "coriolis_latitude" => 44.5,
            "surface_heat_flux" => 50.0,
            "max_flow_snapshots" => 0,
            "allow_analytical_fallback" => false,
            "turbulence_closure" => "nemotke",
            "enable_sea_ice" => true,
            "enable_oxygen" => true
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
            "enable_checkpoint" => true,
            "checkpoint_prefix" => "checkpoint_ParticleTracking",
            "checkpoint_schedule_seconds" => 21600.0,
            "checkpoint_dir" => "outputs/checkpoints",
            "checkpoint_cleanup" => true
        ),
        "tessellation" => Dict{String, Any}(
            "enable_voronoi" => false,
            "n_units" => 5000,
            "prob_core" => 0.8,
            "prob_shallow" => 0.1,
            "prob_deep" => 0.1,
            "min_res_core_km" => 1.5,
            "min_res_shallow_km" => 5.0,
            "min_res_deep_km" => 10.0,
            # Depth band edges as positive depths below the surface. Core is the mid-shelf
            # nursery, shallow the nearshore strip, deep everything below the shelf break.
            "core_depth_min" => 50.0,
            "core_depth_max" => 350.0,
            "shallow_depth_min" => 0.0,
            "shallow_depth_max" => 50.0,
            "deep_depth_min" => 350.0,
            "deep_depth_max" => 6000.0,
            "slope_weighting" => true,
            "slope_factor" => 15.0
        ),
        "hardware" => Dict{String, Any}(
            "use_gpu" => false,
            # Off by default, and deliberately. A GPU request that cannot be honoured must
            # fail loudly: silently running the CPU path instead would leave
            # `resolved_config.toml` recording `use_gpu = true` for a run that never touched
            # a GPU, which is the same silent lie the provenance file exists to prevent. Set
            # it to true only on a machine where degrading is genuinely acceptable.
            "fallback_to_cpu" => false
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => true,
            "animate_hydro" => false,
            "anim_variable" => "dashboard",
            "anim_fps" => 10,
            "anim_format" => "mp4",
            "anim_depth" => -2.5,
            "anim_overlay_particles" => false,
            "anim_output_path" => "",
            "title" => "Regional Marine Lagrangian Particle Tracking & Dispersion"
        ),
        "paths" => Dict{String, Any}(
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

    bathy_source = String(get_val("data", "bathy_source", "etopo2022"))
    coastline_source = String(get_val("data", "coastline_source", "natural_earth_10m"))
    wind_source = String(get_val("data", "wind_source", "open_meteo"))
    wind_time_iso = String(get_val("data", "wind_time_iso", "2020-06-01T00:00:00Z"))

    # Shelf geometry, from [bathymetry]. `inshore_depth` is the water depth at the landward
    # edge of the domain; `shelf_slope` is the total rise across the shelf. Both are read
    # from here rather than hardcoded, so a config can actually shape its own seabed.
    inshore_depth = Float64(get_val("bathymetry", "inshore_depth", -20.0))
    shelf_slope   = Float64(get_val("bathymetry", "shelf_slope", 500.0))

    # Tidal body forcing. `constituents` is an empty list in a config that disables tides, so
    # an empty read is not an error; it simply means no constituents were requested.
    tides_src = String(get_val("tides", "tides_source", "tpxo9"))
    raw_cons = get_val("tides", "constituents", Any[:M2, :S2])
    # Normalised to upper case because `get_tidal_frequency` dispatches on `:M2`/`:S2`
    # case-sensitively, while the atlas file lookup is case-insensitive. Normalising here
    # means a TOML writing "m2" and one writing "M2" behave identically.
    tidal_cons = Symbol[Symbol(uppercase(String(c))) for c in raw_cons]

    enable_tides = Bool(get_val("tides", "enable_tides", true))
    tides_source = tides_src
    tidal_constituents = tidal_cons
    tidal_u = Float64(get_val("tides", "tidal_u_amp", 0.25))
    tidal_v = Float64(get_val("tides", "tidal_v_amp", 0.12))
    s2_u = Float64(get_val("tides", "s2_u_amp", 0.11))
    s2_v = Float64(get_val("tides", "s2_v_amp", 0.05))

    scenario = Symbol(get_val("climate", "scenario", "ssp245"))
    proj_year = Int(get_val("climate", "projection_year", 2050))

    sim_hours = Float64(get_val("hydrodynamics", "sim_duration_hours", 12.0))
    sim_dt    = Float64(get_val("hydrodynamics", "sim_dt_seconds", 120.0))
    # Floor on the adaptive time step, forwarded to the simulation as `min_Δt`. The surface cell
    # sets the admissible step (Δt ≈ target_cfl × Δz_min / w), so a floor above that value stops
    # the wizard stabilising the surface layer.
    min_dt_seconds = Float64(get_val("hydrodynamics", "min_dt_seconds", 2.0))
    adapt_cfl = Bool(get_val("hydrodynamics", "adaptive_cfl", true))
    tgt_cfl   = Float64(get_val("hydrodynamics", "target_cfl", 0.2))
    heat_flux = Float64(get_val("hydrodynamics", "surface_heat_flux", 50.0))
    # `bulk_heat_flux` selects WHICH heat flux is applied: true computes it from the ingested
    # atmospheric surface state, false applies the constant `surface_heat_flux` above. It was
    # parsed into `HydrodynamicConfig` and documented in the TOML but never carried onto
    # `HydrodynamicOptions`, so nothing downstream could read it and the bulk flux was applied
    # unconditionally.
    #
    # It is read from `[atmosphere]`, which is where the TOML actually carries it (next to
    # `atmospheric_source`, `wind_source` and the rest of the atmospheric configuration). The
    # `[hydrodynamics]` section is accepted as a fallback so the key is not silently ignored
    # if it is written there instead.
    bulk_heat_flux = if haskey(config_dict, "atmosphere") &&
                        haskey(config_dict["atmosphere"], "bulk_heat_flux")
        Bool(get_val("atmosphere", "bulk_heat_flux", true))
    else
        Bool(get_val("hydrodynamics", "bulk_heat_flux", true))
    end
    # Cap on how many hydrodynamic snapshots are held in memory when building the 4D flow
    # interpolator for particle tracking. A production run stores thousands of snapshots and
    # materialising all of them for u, v, w and T is hundreds of gigabytes; the snapshots are
    # evenly subsampled down to this count. `0` means "no cap" (load every snapshot).
    max_flow_snapshots = Int(get_val("hydrodynamics", "max_flow_snapshots", 400))
    allow_analytical_fallback = Bool(get_val("hydrodynamics", "allow_analytical_fallback", false))
    # Cap on the background current a larva is advected by. A diverged hydrodynamics cell can
    # otherwise teleport a particle across the domain in a single step; 3 m/s is well above any
    # real Scotian Shelf current. Use a negative value to disable the clamp.
    _mcs = Float64(get_val("hydrodynamics", "max_current_speed", 3.0))
    max_current_speed = _mcs > 0 ? _mcs : nothing

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

    # Stochastic biology parameters
    settlement_stochastic = Bool(get_val("biology", "settlement_stochastic", true))
    cv_molt = Float64(get_val("biology", "cv_molt", 0.25))
    cv_mortality = Float64(get_val("biology", "cv_mortality", 0.25))
    cv_settlement = Float64(get_val("biology", "cv_settlement", 0.25))

    use_gpu      = Bool(get_val("hardware", "use_gpu", false))
    fallback_cpu = Bool(get_val("hardware", "fallback_to_cpu", false))
    interactive  = Bool(get_val("visualization", "interactive_map", true))
    anim_hydro   = Bool(get_val("visualization", "animate_hydro", false))
    anim_var     = Symbol(lowercase(String(get_val("visualization", "anim_variable", "dashboard"))))
    anim_fps     = Int(get_val("visualization", "anim_fps", 10))
    anim_fmt     = String(get_val("visualization", "anim_format", "mp4"))
    anim_depth   = Float64(get_val("visualization", "anim_depth", -2.5))
    anim_overlay = Bool(get_val("visualization", "anim_overlay_particles", false))
    anim_out_path = String(get_val("visualization", "anim_output_path", ""))

    enable_duckdb = Bool(get_val("storage", "enable_duckdb", true))
    duckdb_path   = String(get_val("storage", "duckdb_path", "outputs/particle_tracking.duckdb"))

    output_dir = String(get_val("storage", "output_dir", "outputs"))
    # Ingested data lives in `<output_dir>/inputs`, beside this run's outputs, provenance and
    # resolved config. It is deliberately NOT configurable and NOT shared between runs: a
    # single top-level `inputs/` let a run with a changed domain silently reuse another
    # domain's bathymetry, which is exactly the quiet mismatch `resolved_config.toml` and
    # `data_provenance.json` exist to expose. The cost is disk -- a dataset is fetched once
    # per scenario rather than once per machine -- and that is the intended trade, because it
    # makes each scenario a self-contained, shareable directory.
    input_dir  = joinpath(output_dir, "inputs")
    seed       = Int(get_val("paths", "seed", 42))

    hydro_file = String(get_val("hydrodynamics", "hydro_model_file",
                                get_val("storage", "output_filename", "")))
    hydro_only = Bool(get_val("hydrodynamics", "hydro_only", false))
    track_only = Bool(get_val("hydrodynamics", "track_only", false))
    reuse_hydro = Bool(get_val("hydrodynamics", "reuse_hydro", false))
    run_id_val = String(get_val("storage", "run_id", ""))

    enable_cp = Bool(get_val("storage", "enable_checkpoint",
                             get_val("hydrodynamics", "enable_checkpoint", true)))
    cp_pfx_raw = String(get_val("storage", "checkpoint_prefix",
                                get_val("hydrodynamics", "checkpoint_prefix", "")))
    cp_sched  = Float64(get_val("storage", "checkpoint_schedule_seconds",
                                get_val("hydrodynamics", "checkpoint_schedule_seconds", 0.0)))
    # Field-output cadence, distinct from the checkpoint cadence. Read from [storage] rather
    # than hardcoded downstream, so the two can legitimately differ.
    out_sched = Float64(get_val("storage", "output_schedule_seconds", 21600.0))
    cp_dir    = String(get_val("storage", "checkpoint_dir",
                               get_val("paths", "checkpoint_dir", "")))
    cp_clean  = Bool(get_val("storage", "checkpoint_cleanup", true))
    auto_res  = Bool(get_val("hydrodynamics", "auto_restart", true))

    res_scale = Float64(get_val("grid", "resolution_scale", 1.0))
    v_mode_str = String(get_val("grid", "vertical_stretching_mode", "oceananigans"))
    # Two-segment shelf grid parameters. `vertical_stretching_mode = "two_segment"` splits `nz`
    # into a surface-refined upper segment above `vertical_break_depth` and a coarse lower segment,
    # which is the only way to resolve the 10-400 m active layer while still reaching the deep basin.
    # `nz_above = 0` means "derive it as 60% of nz".
    nz_above_cfg = Int(get_val("grid", "nz_above", 0))
    vertical_break_depth = Float64(get_val("grid", "vertical_break_depth", -400.0))
    v_mode = Symbol(lowercase(v_mode_str))
    v_file = strip(String(get_val("grid", "vertical_grid_file", "")))

    # `vertical_grid_file` is the operator's own file: a relative path is resolved against the
    # working directory they ran from, not against a previous run's inputs. If they name a
    # file, it has to exist -- silently falling back to a generated grid would run the whole
    # simulation on a vertical grid nobody chose. The file is then copied into this run's
    # containerised input directory so the run is self-contained and reproducible.
    #
    # `[grid] depths` is the discretised vertical axis written inline in the config, as
    # positive-downward face depths starting at the surface: depths = [0, 10, 20, ...].
    # It is the alternative to naming a CSV, and the two are mutually exclusive so a config
    # can never silently prefer one over the other.
    raw_depths = get_val("grid", "depths", Any[])
    v_depths = Float64[Float64(d) for d in raw_depths]
    if !isempty(v_depths)
        isempty(v_file) || error(
            "[grid] sets both vertical_grid_file and depths. Use one: an explicit depth list " *
            "is easier to review in the config than a side-car CSV.")
        abs(first(v_depths)) <= 1e-6 || error(
            "[grid] depths must start at 0.0 (the surface); got $(first(v_depths)).")
        all(diff(v_depths) .> 0) || error(
            "[grid] depths must increase monotonically (positive downward): $(v_depths).")
        length(v_depths) == nz + 1 || error(
            "[grid] depths has $(length(v_depths)) entries, but nz = $nz requires " *
            "$(nz + 1) faces.")
        # The deepest face has to reach the bottom of the domain, or the grid silently stops
        # short of the seabed and the deepest water column is left unresolvable.
        abs(last(v_depths) - abs(z_min)) <= 1.0 || error(
            "[grid] depths ends at $(last(v_depths)) m but [domain] z_min is $(z_min) m. " *
            "The deepest face must reach the domain bottom.")
    end
    if !isempty(v_file)
        src = abspath(v_file)
        isfile(src) || error(
            "[grid] vertical_grid_file names '$v_file', which does not exist (resolved to " *
            "$src). Supply the file, or remove the key to derive the grid from the depth " *
            "range and nz instead."
        )
        dst = joinpath(input_dir, basename(v_file))
        isfile(dst) || cp(src, dst, force = true)
        v_file = dst
    end
    atmo_src = Symbol(lowercase(String(get_val("atmosphere", "forcing", "era5_climatology"))))
    obc_src = Symbol(lowercase(String(get_val("boundaries", "ocean_boundary_source", "glorys12v1"))))
    obc_tp = Symbol(lowercase(String(get_val("boundaries", "obc_type", "flather_chapman"))))

    # Initial 3-D hydrography (T, S). Read from `[data]`, independent of the
    # bathymetry/wind `bathy_source` and the lateral `ocean_boundary_source`, because the
    # initial water column is a separate choice from the forcing and the boundaries.
    hydro_src = Symbol(lowercase(String(get_val("data", "hydrography_source", "synthetic"))))
    hydro_month = Int(get_val("data", "hydrography_month", 0))

    # Lateral boundary relaxation (Price & Aumont 2011) parameters
    boundary_method = Symbol(lowercase(String(get_val("boundaries", "method", "relaxation"))))
    sponge_width_val = Float64(get_val("boundaries", "sponge_layer_width", 0.35))
    sponge_tau_val   = Float64(get_val("boundaries", "sponge_timescale_seconds", 3600.0))
    u_inflow_val     = Float64(get_val("boundaries", "u_inflow", -0.15))
    v_inflow_val     = Float64(get_val("boundaries", "v_inflow", 0.05))

    # Maximum allowable time step ceiling
    max_dt_val = Float64(get_val("hydrodynamics", "max_dt_seconds", sim_dt))

    # Voronoi tessellation parameters
    enable_voronoi = Bool(get_val("tessellation", "enable_voronoi", false))
    voronoi_n_units = Int(get_val("tessellation", "n_units", 5000))
    voronoi_p_core = Float64(get_val("tessellation", "prob_core", 0.8))
    voronoi_p_shallow = Float64(get_val("tessellation", "prob_shallow", 0.1))
    voronoi_p_deep = Float64(get_val("tessellation", "prob_deep", 0.1))
    voronoi_min_core = Float64(get_val("tessellation", "min_res_core_km", 1.5))
    voronoi_min_shallow = Float64(get_val("tessellation", "min_res_shallow_km", 5.0))
    voronoi_min_deep = Float64(get_val("tessellation", "min_res_deep_km", 10.0))
    # Depth band edges, as positive depths below the surface (the TOML convention). The
    # tessellator works in signed z (negative downwards), so the sign is applied when the
    # bands are passed to it; see the `generate_depth_stratified_voronoi_units` call.
    voronoi_core_min     = Float64(get_val("tessellation", "core_depth_min", 50.0))
    voronoi_core_max     = Float64(get_val("tessellation", "core_depth_max", 350.0))
    voronoi_shallow_min  = Float64(get_val("tessellation", "shallow_depth_min", 0.0))
    voronoi_shallow_max  = Float64(get_val("tessellation", "shallow_depth_max", 50.0))
    voronoi_deep_min     = Float64(get_val("tessellation", "deep_depth_min", 350.0))
    voronoi_deep_max     = Float64(get_val("tessellation", "deep_depth_max", 6000.0))
    voronoi_slope_wt     = Bool(get_val("tessellation", "slope_weighting", true))
    voronoi_slope_fac    = Float64(get_val("tessellation", "slope_factor", 15.0))

    # Construct HydrodynamicOptions with overrides applied
    return HydrodynamicOptions(;
        domain_lon = (lon_min, lon_max),
        domain_lat = (lat_min, lat_max),
        domain_z = (z_min, z_max),
        grid_size = (nx, ny, nz),
        bathy_source = bathy_source,
    coastline_source = coastline_source,
    wind_source = wind_source,
    wind_time_iso = wind_time_iso,
        inshore_depth = inshore_depth,
        shelf_slope = shelf_slope,
        enable_tides = enable_tides,
    tides_source = tides_source,
    tidal_constituents = tidal_constituents,
        tidal_u_amp = tidal_u,
        tidal_v_amp = tidal_v,
        s2_u_amp = s2_u,
        s2_v_amp = s2_v,
        scenario = scenario,
        projection_year = proj_year,
        sim_dt = sim_dt,
        sim_duration = sim_hours * 3600.0,
        min_dt_seconds = min_dt_seconds,
        adaptive_cfl = adapt_cfl,
        target_cfl = tgt_cfl,
          surface_heat_flux = heat_flux,
          bulk_heat_flux = Bool(bulk_heat_flux),
    max_flow_snapshots = max_flow_snapshots,
    allow_analytical_fallback = allow_analytical_fallback,
    max_current_speed = (max_current_speed isa Real && max_current_speed > 0) ?
                        Float64(max_current_speed) : nothing,
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
        hydro_only = hydro_only,
        track_only = track_only,
        reuse_hydro = reuse_hydro,
        enable_checkpoint = enable_cp,
        checkpoint_prefix = cp_pfx_raw,
        checkpoint_schedule = cp_sched,
        output_schedule_seconds = out_sched,
        checkpoint_dir = cp_dir,
        checkpoint_cleanup = cp_clean,
        auto_restart = auto_res,
        run_id = run_id_val,
        vertical_stretching_mode = v_mode,
        nz_above = nz_above_cfg,
        vertical_break_depth = vertical_break_depth,
        vertical_grid_file = v_file,
    vertical_depths = v_depths,
        resolution_scale = res_scale,
    atmospheric_source = atmo_src,
    ocean_boundary_source = obc_src,
    embedding_lon = embedding_domain_ranges(config_dict).lon,
    embedding_lat = embedding_domain_ranges(config_dict).lat,
    hydrography_source = hydro_src,
    hydrography_month = hydro_month,
        obc_type = obc_tp,
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
        anim_overlay_particles = anim_overlay,
        anim_output_path = anim_out_path,
        boundary_method = boundary_method,
        sponge_layer_width = sponge_width_val,
        sponge_timescale = sponge_tau_val,
        tidal_relaxation_timescale = Float64(get_val("tides", "tidal_relaxation_timescale_seconds", 3600.0)),
        u_inflow = u_inflow_val,
        v_inflow = v_inflow_val,
        max_dt = max_dt_val,
        settlement_stochastic = settlement_stochastic,
        cv_molt = cv_molt,
        cv_mortality = cv_mortality,
        cv_settlement = cv_settlement,
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
            "nz_above" => opts.nz_above,
            "vertical_break_depth" => opts.vertical_break_depth,
            "vertical_grid_file" => opts.vertical_grid_file,
            "depths" => opts.vertical_depths
        ),
        "atmosphere" => Dict{String, Any}(
            "source" => string(opts.atmospheric_source)
        ),
        "boundaries" => Dict{String, Any}(
            "ocean_boundary_source" => string(opts.ocean_boundary_source),
            "obc_type" => string(opts.obc_type),
            "method" => string(opts.boundary_method),
            "sponge_layer_width" => opts.sponge_layer_width,
            "sponge_timescale_seconds" => opts.sponge_timescale,
            "u_inflow" => opts.u_inflow,
            "v_inflow" => opts.v_inflow
        ),
        "data" => Dict{String, Any}(
            "bathy_source" => opts.bathy_source,
            "coastline_source" => opts.coastline_source,
            "wind_source" => opts.wind_source,
            "wind_time_iso" => opts.wind_time_iso,
            "hydrography_source" => string(opts.hydrography_source),
            "hydrography_month" => opts.hydrography_month
        ),
        "bathymetry" => Dict{String, Any}(
            "inshore_depth" => opts.inshore_depth,
            "shelf_slope" => opts.shelf_slope
        ),
        "tides" => Dict{String, Any}(
            "enable_tides" => opts.enable_tides,
            "tides_source" => opts.tides_source,
            "tidal_constituents" => String.(opts.tidal_constituents),
            "tidal_u_amp" => opts.tidal_u_amp,
            "tidal_v_amp" => opts.tidal_v_amp,
            "s2_u_amp" => opts.s2_u_amp,
            "s2_v_amp" => opts.s2_v_amp
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
            "max_dt_seconds" => opts.max_dt,
            "min_dt_seconds" => opts.min_dt_seconds,
            "surface_heat_flux" => opts.surface_heat_flux,
            "max_flow_snapshots" => something(opts.max_flow_snapshots, 0),
            "allow_analytical_fallback" => opts.allow_analytical_fallback,
            "max_current_speed" => something(opts.max_current_speed, 0.0),
            "hydro_model_file" => opts.hydro_model_file,
            "hydro_only" => opts.hydro_only,
            "track_only" => opts.track_only,
            "reuse_hydro" => opts.reuse_hydro
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
            "ascent_target_depth" => opts.ascent_target_depth,
            "settlement_stochastic" => opts.settlement_stochastic,
            "cv_molt" => opts.cv_molt,
            "cv_mortality" => opts.cv_mortality,
            "cv_settlement" => opts.cv_settlement
        ),
        "dvm" => Dict{String, Any}(
            "enable_dvm" => opts.enable_dvm
        ),
        "molting_and_settlement" => Dict{String, Any}(
            "enable_molting" => opts.enable_molting,
            "cv_molt" => opts.cv_molt,
            "cv_mortality" => opts.cv_mortality
        ),
        "storage" => Dict{String, Any}(
            "output_dir" => opts.output_dir,
            "enable_duckdb" => opts.enable_duckdb,
            "duckdb_path" => opts.duckdb_path,
            "run_id" => opts.run_id,
            "enable_checkpoint" => opts.enable_checkpoint,
            "checkpoint_prefix" => opts.checkpoint_prefix,
            "checkpoint_schedule_seconds" => opts.checkpoint_schedule,
            "output_schedule_seconds" => opts.output_schedule_seconds,
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
            "min_res_deep_km" => opts.voronoi_min_res_deep_km,
            "core_depth_min" => opts.voronoi_core_depth_min,
            "core_depth_max" => opts.voronoi_core_depth_max,
            "shallow_depth_min" => opts.voronoi_shallow_depth_min,
            "shallow_depth_max" => opts.voronoi_shallow_depth_max,
            "deep_depth_min" => opts.voronoi_deep_depth_min,
            "deep_depth_max" => opts.voronoi_deep_depth_max,
            "slope_weighting" => opts.voronoi_slope_weighting,
            "slope_factor" => opts.voronoi_slope_factor
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => opts.interactive_map,
            "animate_hydro" => opts.animate_hydro,
            "anim_variable" => string(opts.anim_variable),
            "anim_fps" => opts.anim_fps,
            "anim_format" => opts.anim_format,
            "anim_depth" => opts.anim_depth,
            "anim_overlay_particles" => opts.anim_overlay_particles,
            "anim_output_path" => opts.anim_output_path
        ),
        "paths" => Dict{String, Any}(
            "seed" => opts.seed
        )
    )
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
    nz_above :: Int
    vertical_break_depth :: Float64
    vertical_grid_file       :: String
    vertical_depths          :: Vector{Float64}
    resolution_scale         :: Float64
    bathymetry_source        :: Symbol
    bathy_source             :: String
    coastline_source         :: String
    wind_source              :: String
    wind_time_iso            :: String
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
    tidal_relaxation_timescale :: Float64
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
    settlement_stochastic    :: Bool
    cv_molt                  :: Float64
    cv_mortality             :: Float64
    cv_settlement            :: Float64
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
        :vertical_depths          => opts.vertical_depths,
        :resolution_scale         => opts.resolution_scale,
        :bathymetry_source        => Symbol(opts.bathy_source),
        :bathy_source              => opts.bathy_source,
        :coastline_source          => opts.coastline_source,
        :wind_source               => opts.wind_source,
        :wind_time_iso             => opts.wind_time_iso,
        :inshore_depth            => opts.inshore_depth,
        :shelf_slope              => opts.shelf_slope,
        :atmospheric_source       => opts.atmospheric_source,
        :drag_formulation         => :garratt_1977,
        :bulk_heat_flux           => true,
        :climatology              => false,
        :mhw_temp_anomaly         => opts.scenario == :mhw ? 3.5 : 0.0,
        :ocean_boundary_source    => opts.ocean_boundary_source,
        :hydrography_source       => opts.hydrography_source,
        :hydrography_month        => opts.hydrography_month,
        :obc_type                 => opts.obc_type,
        :sponge_layer_width       => opts.sponge_layer_width,
        :sponge_timescale         => opts.sponge_timescale,
        :enable_tides             => opts.enable_tides,
        :tides_source             => Symbol(opts.tides_source),
        :tidal_constituents       => opts.tidal_constituents,
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
        :max_dt_seconds           => opts.max_dt,
        # Was hardcoded to 2.0, so the second options struct (built from this dict at the call
        # site below) always carried 2.0 regardless of what the TOML asked for. Use the value
        # that was actually parsed.
        :min_dt_seconds           => opts.min_dt_seconds,
        :coriolis_latitude        => 44.5,
        :divergence_limit         => 20.0,
        :output_dir               => opts.output_dir,
        :output_filename          => isempty(opts.hydro_model_file) ?
                                     "hydrodynamics_output.jld2" :
                                     basename(opts.hydro_model_file),
        :output_schedule_seconds  => opts.output_schedule_seconds,
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
        Symbol(d[:vertical_stretching_mode]),         String(d[:vertical_grid_file]),
        Float64[Float64(x) for x in d[:vertical_depths]],
        Float64(d[:resolution_scale]), Symbol(d[:bathymetry_source]),
        String(d[:bathy_source]), String(d[:coastline_source]), String(d[:wind_source]), String(d[:wind_time_iso]), Float64(d[:inshore_depth]),
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
        :settlement_stochastic    => opts.settlement_stochastic,
        :cv_molt                  => opts.cv_molt,
        :cv_mortality             => opts.cv_mortality,
        :cv_settlement            => opts.cv_settlement,
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
        Float64(d[:settlement_max_temp]), Bool(d[:settlement_stochastic]),
        Float64(d[:cv_molt]), Float64(d[:cv_mortality]), Float64(d[:cv_settlement]),
        Bool(d[:enable_duckdb]),
        String(d[:duckdb_path]), String(d[:run_id]), Int(d[:seed])
    )
end

