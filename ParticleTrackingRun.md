# Hydrodynamic Model for Snow Crab Larval Particle Tracking

## Overview

This workflow establishes a regional hydrodynamic model of the Scotian Shelf and
adjacent slope waters using [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).
The hydrodynamic solution provides time-evolving advection ($\boldsymbol{u} = (u, v, w)$)
and turbulent diffusion ($\kappa_h, \kappa_v$) fields. These physical fields drive
Lagrangian particle tracking of planktonic (Zoea I, Zoea II) and semiplanktonic (Megalopa)
snow crab (*Chionoecetes opilio*) larvae undergoing stage-specific vertical migration,
temperature-dependent development, and spatial dispersal across coastal nursery habitats.

The underlying routines are modularized in `src/` under the `ParticleTracking` package:
- [`src/open_data.jl`](file:///c:/home/jae/projects/ParticleTracking/src/open_data.jl): Real-world data ingestion (NOAA ERDDAP ETOPO/GEBCO, NOAA Blended Sea Winds), coordinate normalization, 2D bilinear regridding, and Large & Pond (1981) drag formulation.
- [`src/tides.jl`](file:///c:/home/jae/projects/ParticleTracking/src/tides.jl): Astronomical tidal constituents ($M_2, S_2, N_2, K_1, O_1$), barotropic tidal momentum body forcing, and Simpson-Hunter tidal mixing front diagnostics.
- [`src/climate_scenarios.jl`](file:///c:/home/jae/projects/ParticleTracking/src/climate_scenarios.jl): CMIP6 climate change scenarios (SSP1-2.6, SSP2-4.5, SSP5-8.5, Marine Heatwaves), thermal stratification adjustments, Pelagic Larval Duration (PLD), and larval thermal mortality modeling.
- [`src/synthetic_data.jl`](file:///c:/home/jae/projects/ParticleTracking/src/synthetic_data.jl): Synthetic bathymetry and surface wind stress generation and inspection.
- [`src/grid_bathymetry.jl`](file:///c:/home/jae/projects/ParticleTracking/src/grid_bathymetry.jl): Spherical curvilinear grid and immersed boundary configuration (including automated regridding from real bathymetry).
- [`src/hydrodynamic_model.jl`](file:///c:/home/jae/projects/ParticleTracking/src/hydrodynamic_model.jl): Hydrostatic Boussinesq primitive equation model construction with optional tidal momentum body forcing.
- [`src/simulation.jl`](file:///c:/home/jae/projects/ParticleTracking/src/simulation.jl): Time integration orchestration, adaptive CFL time stepping, stability watchdogs, and 4D JLD2 spatiotemporal output interpolation.
- [`src/larval_behavior.jl`](file:///c:/home/jae/projects/ParticleTracking/src/larval_behavior.jl): Diel Vertical Migration (DVM), $M_2$ tidal velocity superposition, in situ thermal degree-day molting, and benthic settlement suitability filtering.
- [`src/empirical_analysis.jl`](file:///c:/home/jae/projects/ParticleTracking/src/empirical_analysis.jl): Empirical advection & diffusivity estimation, gridded recruitment/retention metrics, thermal exposure accounting, transition probability connectivity matrices, and multi-layer NetCDF/JLD2 exports.
- [`src/visualization.jl`](file:///c:/home/jae/projects/ParticleTracking/src/visualization.jl): CairoMakie spatial trajectory maps, DVM depth profiles, 2D settlement density kernels, empirical velocity quivers, annotated connectivity heatmaps, and cross-scenario comparisons.

---

## 1. Setup & Package Loading

```julia
cd("c:/home/jae/projects/ParticleTracking")

# Load all dependencies and functional modules via unified initializer
include("src/init.jl")
```

---

## 2. Environmental Data Ingestion (Real vs. Synthetic)

### Option 2A: Fetch Real-World Data from Open Repositories

Download real-world seafloor bathymetry from NOAA CoastWatch / ERDDAP (`etopo180` or GEBCO)
and observed 10m surface vector winds (`erdBSwinds1day`), converting winds to kinematic
surface wind stress via the Large & Pond (1981) parameterization:

```julia
# 1. Fetch real bathymetry from NOAA ERDDAP (Scotian Shelf bounding box)
real_bathy_file = fetch_open_bathymetry(lon_range=(-68.0, -57.0),
                                        lat_range=(42.0, 47.0),
                                        output_path="inputs/real_bathymetry.nc",
                                        dataset_id="etopo180")

# 2. Fetch real 10m surface winds for a target date
real_wind_file = fetch_open_surface_winds(lon_range=(-68.0, -57.0),
                                          lat_range=(42.0, 47.0),
                                          time_iso="2023-06-01T00:00:00Z",
                                          output_path="inputs/real_surface_winds.nc")

# 3. Inspect the retrieved NetCDF datasets
inspect_netcdf(real_bathy_file)
inspect_netcdf(real_wind_file)

# 4. Example: Compute kinematic wind stress from 10m wind velocity
u10, v10 = 7.5, 3.2 # Zonal and meridional 10m wind speed (m/s)
tau_x, tau_y = wind_speed_to_kinematic_stress(u10, v10)
println("Kinematic wind stress: tau_x = $(round(tau_x, digits=6)), tau_y = $(round(tau_y, digits=6)) m^2/s^2")
```

### Option 2B: Generate Synthetic Benchmark Data

Alternatively, generate synthetic bathymetry and idealized cyclic surface wind forcing:

```julia
synth_bathy_file   = generate_synthetic_bathymetry("inputs/nova_scotia_bathymetry.nc",
                                                   lon_range=(-68.0, -57.0),
                                                   lat_range=(42.0, 47.0),
                                                   n_lon=50, n_lat=50,
                                                   inshore_depth=-100.0,
                                                   shelf_slope=500.0)

synth_forcing_file = generate_synthetic_forcing("inputs/surface_forcing.nc",
                                                lon_range=(-68.0, -57.0),
                                                lat_range=(42.0, 47.0),
                                                time_range=(0.0, 86400.0),
                                                n_lon=50, n_lat=50, n_time=24,
                                                tau_x_amplitude=0.1)
```

---

## 3. Grid & Immersed Boundary Topography

Discretize the regional domain using a `LatitudeLongitudeGrid`. For real-world bathymetry,
`build_immersed_grid_from_real_data` automatically normalizes coordinates and regrids
the elevation matrix onto the target model resolution via 2D bilinear interpolation:

```julia
# 1. Define base spherical grid
underlying_grid = build_shelf_grid(lon_range=(-68.0, -57.0),
                                   lat_range=(42.0, 47.0),
                                   z_range=(-1000.0, 0.0),
                                   grid_size=(50, 50, 10))

# 2. Build immersed boundary (using real bathymetry with automatic 2D regridding)
grid = build_immersed_grid_from_real_data(underlying_grid, real_bathy_file)

println("Immersed boundary grid initialized:")
println(grid)
```

---

## 4. Hydrodynamic Model with Optional Astronomical Tidal Forcing

Construct a `HydrostaticFreeSurfaceModel` with Coriolis acceleration, SeawaterBuoyancy,
active temperature and salinity tracers, surface wind stress boundary conditions, and
optional barotropic tidal body forcing ($M_2, S_2$):

```julia
# 1. (Optional) Configure barotropic M2 semi-diurnal tidal momentum body forcing
tides = build_tidal_body_forcing(constituents=[:M2],
                                 u_amplitudes=Dict(:M2 => 0.25),
                                 v_amplitudes=Dict(:M2 => 0.12))

# 2. Instantiate hydrostatic model with wind stress and tidal forcing
model = build_hydrodynamic_model(grid,
                                 coriolis_latitude=45.0,
                                 surface_wind_stress_x=tau_x,
                                 surface_wind_stress_y=tau_y,
                                 tidal_forcing=tides,
                                 ν=1e-2, κ=1e-2,
                                 tracers=(:T, :S))

# 3. Set initial thermal and haline stratification (baseline)
set_initial_stratification!(model,
                            surface_temperature=15.0,
                            temperature_gradient=0.01,
                            salinity=35.0)

# 4. Diagnose Simpson-Hunter tidal mixing index (χ = log10(h / U³)) across shelf banks
chi_bank  = simpson_hunter_parameter(40.0, 1.1)  # Shallow bank (well-mixed if χ < 1.5)
chi_shelf = simpson_hunter_parameter(150.0, 0.2) # Deep shelf (stratified if χ > 2.0)
println("Simpson-Hunter Mixing Parameter: Bank = $(round(chi_bank, digits=2)), Shelf = $(round(chi_shelf, digits=2))")
```

---

## 5. Optional Climate Forcing Scenarios & Larval Thermal Ecology

Simulate future ocean conditions under IPCC CMIP6 Shared Socioeconomic Pathways
(`:ssp126`, `:ssp245`, `:ssp370`, `:ssp585`) or transient Marine Heatwaves (`:marine_heatwave`).
Climate anomalies alter vertical stratification, pelagic larval drift duration, and survival rates:

```julia
# 1. Inspect climate scenario deltas for Scotian Shelf (e.g. SSP2-4.5 in 2050)
deltas = get_climate_scenario_deltas(:ssp245, year=2050)
println("Climate Scenario: ", deltas.description)
println("  Surface warming: +$(deltas.ΔT_surface) °C")
println("  Deep warming:    +$(deltas.ΔT_deep) °C")
println("  Freshening:      $(deltas.ΔS_surface) PSU")

# 2. Apply climate anomalies to the model's thermal and haline fields
apply_climate_scenario!(model, scenario=:ssp245, year=2050)

# 3. Calculate temperature-dependent Pelagic Larval Duration (PLD)
t_ambient = 5.5 # Mean larval exposure temperature (°C)
pld_days  = temperature_dependent_pld(t_ambient)
mortality = larval_thermal_mortality_rate(t_ambient)

println("Larval Development at $(t_ambient) °C:")
println("  Estimated PLD:      $(round(pld_days, digits=1)) days")
println("  Daily Mortality:    $(round(mortality * 100, digits=2)) %/day")
```

---

## 6. Simulation Setup & Execution

Orchestrate time stepping with adaptive CFL monitoring, stability watchdogs, and JLD2 field output writers:

```julia
# 1. Configure simulation integration parameters and JLD2 outputs
simulation = setup_hydrodynamic_simulation(model,
                                           Δt=2minutes,
                                           stop_time=12hours,
                                           adaptive_time_step=true,
                                           target_cfl=0.2,
                                           output_dir="outputs",
                                           output_filename="nova_scotia_hydrodynamics.jld2",
                                           output_schedule=100)

# 2. Execute hydrodynamic model integration
run_hydrodynamic_simulation!(simulation)
```

---

## 7. Coupling with Snow Crab Larval Particle Tracking & Life History

Snow crab larvae transition through planktonic zoeal stages (Zoea I and II) before
entering the semiplanktonic megalopal stage. In-situ thermal degree-days ($DD = \int \max(0, T - T_0) dt$)
trigger ontogenetic molting, while settled Instar I juveniles are evaluated against
cold-water benthic nursery habitat criteria:

```julia
# 1. Initialize larval particle cohort across the nursery spawning grounds
n_larvae = 200
larvae = initialize_larval_particles(n_larvae,
                                    lon_range=(-64.0, -62.0),
                                    lat_range=(43.5, 45.0),
                                    depth_range=(-40.0, -10.0),
                                    stage=:zoea1)

# 2. Multi-day Lagrangian trajectory simulation with tidal oscillations and degree-day molting
flow_field_fn(lon, lat, z, t) = (0.05 + 0.02 * sin(t / 43200.0), 0.02, 0.0001)
temp_field_fn(lon, lat, z, t) = 4.0 + 0.01 * z # In-situ temperature
bathy_field_fn(lon, lat)      = -180.0         # Local seabed elevation

trajectories = track_larval_cohort(larvae,
                                   velocity_fn=flow_field_fn,
                                   temperature_fn=temp_field_fn,
                                   bathymetry_fn=bathy_field_fn,
                                   total_duration=86400.0 * 10, # 10 days
                                   dt=300.0,                    # 5 min step
                                   κ_h=10.0, κ_v=1e-4,
                                   is_lat_lon=true,
                                   enable_tides=true,
                                   enable_molting=true)

println("Tracked $(n_larvae) larvae over $(length(trajectories.times)) time steps.")
println("Settlement summary: ", count(==( :settled_successful), trajectories.settlement_status), " settled successfully.")
```

---

## 8. Empirical Movement, Recruitment, & Demographic Connectivity Analysis

Subject the simulated Lagrangian particle tracks to empirical movement analysis to estimate
spatially-resolved velocity and turbulent diffusivity fields, recruitment metrics, thermal exposure,
and demographic connectivity matrices:

```julia
# 1. Estimate empirical advection velocity (u_emp, v_emp) and turbulent diffusivity (D_emp)
emp_mov = estimate_empirical_movement(trajectories,
                                      lon_bins=range(-68.0, -57.0, length=30),
                                      lat_bins=range(42.0, 47.0, length=30))
println("Empirical Movement Analysis:")
println("  Mean empirical speed: ", round(nanmean(emp_mov.speed_mean), digits=4), " m/s")
println("  Mean empirical diffusivity: ", round(nanmean(emp_mov.diffusivity), digits=2), " m^2/s")

# 2. Compute gridded larval recruitment, settlement density, and self-retention
rec_metrics = compute_gridded_recruitment_metrics(trajectories,
                                                  lon_bins=range(-68.0, -57.0, length=30),
                                                  lat_bins=range(42.0, 47.0, length=30))
println("Total settled recruits: ", sum(rec_metrics.settlement_density))
println("Successful recruits in nursery: ", sum(rec_metrics.successful_settlement_density))

# 3. Compute gridded thermal exposure and degree-days
therm_metrics = compute_gridded_thermal_metrics(trajectories,
                                                lon_bins=range(-68.0, -57.0, length=30),
                                                lat_bins=range(42.0, 47.0, length=30))

# 4. Derive macro-regional connectivity matrix between Crab Fishing Areas (CFAs)
cfa_definitions = [
    (name="CFA 20-22 (Eastern NS)", lon=(-62.0, -57.0), lat=(44.5, 47.5)),
    (name="CFA 23-24 (Middle Shelf)", lon=(-64.5, -60.0), lat=(43.0, 45.5)),
    (name="CFA 4X (Southwest NS)", lon=(-68.0, -64.0), lat=(42.0, 44.5)),
    (name="Offshore / Slope", lon=(-68.0, -57.0), lat=(40.0, 43.0))
]
conn = compute_empirical_connectivity(trajectories, strata_definitions=cfa_definitions)
println("\nTransition Probability Matrix (P_ij):")
display(round.(conn.matrix, digits=3))

# 5. Export comprehensive multi-variable NetCDF and JLD2 archives
export_larval_dispersal_netcdf("outputs/larval_dispersal_analysis.nc",
                               trajectories=trajectories,
                               strata_definitions=cfa_definitions)
export_larval_dispersal_jld2("outputs/larval_dispersal_analysis.jld2",
                             trajectories=trajectories,
                             strata_definitions=cfa_definitions)
```

---

## 9. Visualization & Cross-Scenario Comparisons

Render particle drift trajectories, empirical velocity quivers, annotated connectivity matrices,
and multi-scenario climate comparisons:

```julia
# 1. Plot 2D trajectories over Scotian Shelf bathymetry
bathy_grid_data = load_bathymetry_from_netcdf(bathy_path)
fig_tracks = plot_particle_trajectories(trajectories,
                                        bathymetry_data=bathy_grid_data,
                                        output_path="outputs/larval_trajectories.png")

# 2. Plot Diel Vertical Migration (DVM) depth profiles over time
fig_dvm = plot_vertical_migration_profiles(trajectories,
                                           sample_indices=1:5,
                                           output_path="outputs/dvm_depth_profiles.png")

# 3. Plot 2D settlement nursery density distribution
fig_density = plot_larval_dispersal_density(trajectories,
                                            output_path="outputs/settlement_density.png")

# 4. Plot empirical velocity vector arrows over turbulent diffusivity
fig_emp = plot_empirical_movement_field(emp_mov,
                                        output_path="outputs/empirical_movement_field.png")

# 5. Plot annotated macro-regional connectivity matrix
fig_conn = plot_connectivity_matrix(conn,
                                     output_path="outputs/regional_connectivity_matrix.png")

# 6. Plot larval thermal exposure and cumulative degree-days
fig_therm = plot_thermal_exposure_map(therm_metrics,
                                      output_path="outputs/thermal_exposure_map.png")

# 7. Plot recruitment and settlement summary
fig_rec = plot_recruitment_summary(rec_metrics,
                                   output_path="outputs/recruitment_summary.png")

# 8. Compare dispersal patterns across multiple climate scenarios
flow_ssp585(lon, lat, z, t) = (0.07 + 0.03 * sin(t / 43200.0), 0.03, 0.0001)
traj_ssp585 = track_larval_cohort(larvae, velocity_fn=flow_ssp585,
                                  total_duration=86400.0 * 10, dt=300.0)

scenario_comp = Dict(
    :historical => trajectories,
    :ssp585_2050 => traj_ssp585
)
fig_comp = compare_scenario_dispersal(scenario_comp,
                                      output_path="outputs/climate_scenario_comparison.png")
```

---

## 10. Streamlined Climate- & Tide-Coupled Production Pipeline

The entire real-data hydrodynamic modeling, tidal forcing, climate scenario integration,
particle tracking, empirical movement analysis, and visualization workflow can be executed in a unified script:

```julia
cd("c:/home/jae/projects/ParticleTracking")

# Load environment via unified initializer
include("src/init.jl")

# 1. Real bathymetry retrieval & wind stress calculation
bathy_path   = fetch_open_bathymetry(lon_range=(-68.0, -57.0), lat_range=(42.0, 47.0))
tau_x, tau_y = wind_speed_to_kinematic_stress(8.0, 2.5)

# 2. Grid & immersed real seafloor configuration with automated 2D regridding
base_grid = build_shelf_grid(lon_range=(-68.0, -57.0), lat_range=(42.0, 47.0),
                             z_range=(-1000.0, 0.0), grid_size=(50, 50, 10))
grid      = build_immersed_grid_from_real_data(base_grid, bathy_path)

# 3. Model instantiation with M2 tidal body forcing & SSP2-4.5 climate stratification
tides = build_tidal_body_forcing(constituents=[:M2], u_amplitudes=Dict(:M2 => 0.3))
model = build_hydrodynamic_model(grid, coriolis_latitude=45.0,
                                 surface_wind_stress_x=tau_x,
                                 surface_wind_stress_y=tau_y,
                                 tidal_forcing=tides)
apply_climate_scenario!(model, scenario=:ssp245, year=2050)

# 4. Simulation & execution
sim = setup_hydrodynamic_simulation(model, Δt=2minutes, stop_time=12hours,
                                    adaptive_time_step=true, target_cfl=0.2,
                                    output_dir="outputs",
                                    output_filename="nova_scotia_ssp245_tides.jld2")
run_hydrodynamic_simulation!(sim)

# 5. Initialize snow crab larvae cohort & track dispersal with tides and molting
larval_cohort = initialize_larval_particles(500, stage=:zoea1)
flow_fn(lon, lat, z, t) = (0.05, 0.02, 0.0001)
trajs = track_larval_cohort(larval_cohort, velocity_fn=flow_fn,
                            total_duration=86400.0 * 5, enable_tides=true, enable_molting=true)

# 6. Empirical movement, recruitment & connectivity analysis + NetCDF export
nc_out = export_larval_dispersal_netcdf("outputs/production_larval_dispersal.nc",
                                        trajectories=trajs)

# 7. Generate summary visualization maps
bathy_data = load_bathymetry_from_netcdf(bathy_path)
plot_particle_trajectories(trajs, bathymetry_data=bathy_data,
                           output_path="outputs/production_larval_tracks.png")
plot_larval_dispersal_density(trajs, output_path="outputs/production_settlement_density.png")

emp_mov = estimate_empirical_movement(trajs)
plot_empirical_movement_field(emp_mov, output_path="outputs/production_empirical_movement.png")

conn = compute_empirical_connectivity(trajs)
plot_connectivity_matrix(conn, output_path="outputs/production_connectivity_matrix.png")

# 8. Export interactive HTML5 + Leaflet oceanographic dashboard
export_interactive_tracks_html(
    "outputs/production_interactive_tracks.html",
    trajectories = trajs,
    bathymetry = bathy_data,
    strata_definitions = conn.strata_names,
    title = "Scotian Shelf Snow Crab Larval Dispersal & Connectivity Dashboard"
)

println("Climate- and tide-coupled production pipeline complete with diagnostics and interactive map exported.")
```

---

## 11. Command-Line Interface (CLI) Complete Reference

The workflow script [`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl) provides a comprehensive, flag-driven command-line interface. It allows granular execution of specific workflow segments, runtime parameter customization, hardware selection, and DuckDB analytical queries.

### Basic Syntax
```bash
julia --project=. ParticleTrackingRun.jl [FLAGS...]
```

### Complete Options Reference

| Category | Flag / Option | Argument Type | Default Value | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Execution Modes** | `--all` | Flag | — | Execute entire 8-segment end-to-end pipeline. |
| | `--quick`, `-q` | Flag | — | Fast debug mode (coarse grid $15\times 15\times 5$, $1\text{ h}$ sim, $2\text{ d}$ track). |
| | `--segment=<name>` | String | `all` | Target segment: `data`, `grid`, `model`, `climate`, `sim`, `track`, `metrics`, `viz`, `all`. |
| **Segment Flags** | `--data` | Flag | — | Run environmental data ingestion (Option 2A/2B). |
| | `--grid` | Flag | — | Run grid and immersed boundary construction. |
| | `--model` | Flag | — | Run hydrodynamic model setup with tidal forcing. |
| | `--climate` | Flag | — | Run CMIP6 climate scenario & larval thermal biology. |
| | `--sim`, `--simulation` | Flag | — | Run Oceananigans hydrodynamic time integration. |
| | `--track`, `--tracking` | Flag | — | Run Lagrangian larval particle tracking. |
| | `--metrics` | Flag | — | Compute retention, empirical diffusion & connectivity. |
| | `--viz`, `--visualize` | Flag | — | Generate CairoMakie spatial figures and Leaflet map. |
| **DuckDB Analytics** | `--duckdb` | Flag | `true` | Enable DuckDB analytical storage archiving. |
| | `--no-duckdb` | Flag | — | Disable DuckDB storage archiving. |
| | `--db-path=<path>` | String | `outputs/particle_tracking.duckdb` | Custom file path for the DuckDB analytical database. |
| | `--list-runs` | Flag | — | Query and print all archived simulation runs table. |
| | `--compare-scenarios` | Flag | — | Query and print multi-scenario comparative metrics. |
| | `--model-average` | Flag | — | Compute ensemble model-averaged connectivity ($P_{ij} \pm \sigma$). |
| | `--export-parquet` | Flag | — | Export all DuckDB tables to Apache Parquet files. |
| **Architecture** | `--gpu`, `--cuda` | Flag | — | Enable NVIDIA CUDA GPU acceleration for hydrodynamics. |
| | `--cpu` | Flag | `true` | Execute on multi-threaded CPU. |
| | `--fallback-cpu` | Flag | — | Automatically fall back to CPU if CUDA GPU is absent. |
| **Visualization** | `--interactive` | Flag | `true` | Export standalone interactive HTML5 Leaflet map. |
| | `--no-interactive` | Flag | — | Disable interactive HTML map generation. |
| **Spatial Domain** | `--lon=<min,max>` | Real,Real | `-68.0,-57.0` | Longitude bounding range in degrees East. |
| | `--lat=<min,max>` | Real,Real | `42.0,47.0` | Latitude bounding range in degrees North. |
| | `--depth-range=<min,max>`| Real,Real | `-1000.0,0.0` | Vertical depth range in meters ($z_{\text{bottom}}, z_{\text{surface}}$). |
| | `--grid=<nx,ny,nz>` | Int,Int,Int | `50,50,10` | Grid cell resolution ($N_\lambda, N_\phi, N_z$). |
| | `--nx=<int>`, `--ny=<int>`, `--nz=<int>` | Int | `50`, `50`, `10` | Individual grid axis dimensions. |
| **Data & Tides** | `--real` | Flag | — | Fetch real-world NOAA ERDDAP bathymetry and winds. |
| | `--synthetic` | Flag | `true` | Generate idealized synthetic shelf topography & winds. |
| | `--tides` | Flag | `true` | Enable astronomical tidal body forcing ($M_2$). |
| | `--no-tides` | Flag | — | Disable tidal body forcing. |
| | `--tidal-u=<val>` | Float (m/s) | `0.25` | Semi-major barotropic tidal current amplitude. |
| | `--tidal-v=<val>` | Float (m/s) | `0.12` | Semi-minor barotropic tidal current amplitude. |
| **Climate Scenarios** | `--scenario=<name>` | Symbol | `ssp245` | Climate scenario: `historical`, `ssp126`, `ssp245`, `ssp585`, `mhw`. |
| | `--year=<int>` | Int | `2050` | Climate projection benchmark horizon year. |
| **Hydro Simulation** | `--duration=<hours>` | Float (hrs) | `12.0` | Hydrodynamic simulation duration in hours. |
| | `--sim-duration=<sec>` | Float (sec) | `43200.0` | Hydrodynamic simulation duration in seconds. |
| | `--sim-dt=<sec>` | Float (sec) | `120.0` | Initial hydrodynamic time integration step. |
| | `--adaptive-cfl` | Flag | `true` | Enable advective CFL-limited adaptive time stepping. |
| | `--no-adaptive-cfl` | Flag | — | Disable adaptive CFL time stepping. |
| | `--target-cfl=<val>` | Float | `0.2` | Target advective Courant-Friedrichs-Lewy limit. |
| **Lagrangian Tracking** | `--particles=<int>` | Int | `100` | Number of larval particles to initialize and track. |
| | `--track-duration=<days>`| Float (days)| `5.0` | Larval cohort tracking duration in days. |
| | `--track-dt=<sec>` | Float (sec) | `300.0` | Lagrangian advection time step in seconds. |
| | `--min-depth=<meters>` | Float (m) | `100.0` | Minimum seabed depth for marine water placement. |
| | `--dvm` | Flag | `true` | Enable stage-dependent Diel Vertical Migration. |
| | `--no-dvm` | Flag | — | Disable Diel Vertical Migration. |
| | `--molting` | Flag | `true` | Enable degree-day thermal molting & mortality. |
| | `--no-molting` | Flag | — | Disable thermal molting calculations. |
| | `--diff-h=<val>` | Float ($m^2/s$)| `10.0` | Horizontal turbulent diffusivity ($\kappa_h$). |
| | `--diff-v=<val>` | Float ($m^2/s$)| `1e-4` | Vertical turbulent diffusivity ($\kappa_v$). |
| **I/O & Environment** | `--output-dir=<path>` | String | `outputs` | Target directory for outputs and figures. |
| | `--input-dir=<path>` | String | `inputs` | Cache directory for input NetCDF datasets. |
| | `--seed=<int>` | Int | `42` | Pseudorandom number generator seed. |
| | `--help`, `-h` | Flag | — | Print command-line help and flag reference. |

---

## 12. Hardware Acceleration & Interactive Visualization

### GPU / CUDA Hardware Acceleration
Oceananigans natively supports NVIDIA CUDA GPU execution using unified memory and streaming multiprocessors. The workflow integrates [`resolve_architecture`](file:///c:/home/jae/projects/ParticleTracking/src/architecture.jl):
```bash
# Execute on NVIDIA CUDA GPU with automatic fallback to CPU
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu

# Explicit CPU execution
julia --project=. ParticleTrackingRun.jl --all --cpu
```

### Standalone Interactive HTML5 Map Dashboard
The visualizer exports an interactive Leaflet.js dashboard (`outputs/interactive_larval_tracks.html`) that can be opened in any browser:
- **Base Maps**: ESRI Ocean Basemap, CartoDB Dark Matter, OpenStreetMap.
- **Dynamic Playback**: Simulation time scrubber ($t = 0 \to 5\text{ days}$) with Play/Pause animation.
- **Stage Coloration**: Zoea I–IV, Megalopa, Instar I (nursery recruit), and thermal mortality.
- **Live Telemetry HUD**: Interactive popups with depth ($z$), temperature ($T$), degree-days ($DD$), and survival probability ($S(t)$).
- **Management Strata**: Bounding overlays for CFAs 20–22, 23–24, and 4X.

---

## 13. DuckDB Analytical Storage, Scenario Management & Ensemble Model Averaging

Similar to the BSTM modeling framework, `ParticleTracking.jl` includes a high-performance **DuckDB** analytical storage backend (`src/storage_duckdb.jl`) for persisting multi-scenario simulation runs, millions of trajectory steps, demographic transition matrices, and gridded dispersal fields into a single relational database (`outputs/particle_tracking.duckdb`).

### Database Relational Schema
1. **`simulation_runs`**: Metadata for each run (scenario name, projection year, $N_{\text{particles}}$, duration, time step, physical & biological options, seed, timestamps).
2. **`particle_trajectories`**: Columnar time series of particle coordinates $(\lambda, \phi, z)$, ambient temperature $T$, degree-days $DD$, survival probability $S(t)$, developmental stage, and settlement state.
3. **`recruitment_metrics`**: Summary outcomes per cohort (settlement success rate, mean PLD, mean degree-days, thermal mortality, mean dispersal displacement).
4. **`connectivity_transitions`**: Demographic transition probabilities $P_{ij}$ and raw transit counts between spatial management strata (e.g. CFAs 20–22, 23–24, 4X).
5. **`gridded_dispersal_summary`**: Spatial matrices of empirical velocities $(\bar{u}, \bar{v})$, turbulent diffusivity $D_{\text{emp}}$, settlement density, and thermal exposure.

### CLI Querying & Model Averaging Commands
```bash
# Execute workflow and archive results to DuckDB (enabled by default)
julia --project=. ParticleTrackingRun.jl --all --scenario=ssp245 --year=2050
julia --project=. ParticleTrackingRun.jl --all --scenario=ssp585 --year=2050

# Query and display all archived simulation runs in DuckDB
julia --project=. ParticleTrackingRun.jl --list-runs

# Multi-scenario comparative analytics
julia --project=. ParticleTrackingRun.jl --compare-scenarios

# Multi-scenario Bayesian / ensemble model averaging (connectivity & recruitment)
julia --project=. ParticleTrackingRun.jl --model-average

# Export DuckDB tables to Apache Parquet format (Python/R/BSTM interoperability)
julia --project=. ParticleTrackingRun.jl --export-parquet
```

### Julia API Usage
```julia
using ParticleTracking

# Open DuckDB analytical database connection
db = open_duckdb_storage("outputs/particle_tracking.duckdb")

# 1. Query table of all simulation runs
runs_df = list_simulation_runs(db; scenario = "ssp245")

# 2. Extract trajectory DataFrame with particle/stage filtering
traj_df = load_trajectories_df(db, "run_ssp245_2050", stage = :megalopa)

# 3. Retrieve connectivity transition matrix
conn = load_connectivity_matrix(db, "run_ssp245_2050")

# 4. Multi-scenario comparative summary
comparison_df = compare_scenarios(db)

# 5. Compute weighted ensemble model-averaged connectivity and recruitment
ens = compute_ensemble_model_average(
    db,
    ["historical", "ssp245", "ssp585"],
    weights = [0.2, 0.5, 0.3]
)
println("Ensemble mean connectivity matrix: ", ens.mean_connectivity)
println("Ensemble connectivity uncertainty (std): ", ens.std_connectivity)

# 6. Export to Apache Parquet
export_duckdb_to_parquet(db, "outputs/parquet")

close_duckdb_storage(db)
```

---

## 14. Centralized Configuration Management (`inputs/ParticalTracking.config`)

To streamline workflow reproducibility and parameter configuration across multiple scenarios, all user-configurable defaults are declared in a centralized TOML configuration file at `inputs/ParticalTracking.config`.

### Configuration Sections
- `[domain]`: Bounding coordinates ($\lambda_{\min}, \lambda_{\max}, \phi_{\min}, \phi_{\max}, z_{\min}, z_{\max}$).
- `[grid]`: Numerical grid dimensions ($N_x, N_y, N_z$).
- `[data]`: Environmental data mode (`synthetic` vs `real`), dataset IDs, and synthetic shelf parameters.
- `[tides]`: Barotropic tidal forcing options ($M_2$ constituent amplitude, period, phase).
- `[climate]`: Climate scenario selection (`ssp126`, `ssp245`, `ssp585`, `mhw`), projection and baseline years.
- `[hydrodynamics]`: Oceananigans integration duration, initial time step, adaptive CFL target.
- `[biology]`: Larval cohort size ($N$), pelagic drift duration, tracking time step, minimum seabed depth ($100\text{ m}$), turbulent diffusivities ($\kappa_h, \kappa_v$).
- `[dvm]`: Stage-specific Diel Vertical Migration daytime/nighttime target depths and swimming speeds.
- `[molting_and_settlement]`: Degree-day thresholds ($150, 310, 510\text{ DD}$), thermal mortality parameters, and benthic nursery suitability windows ($-250\text{ m} \le z \le -50\text{ m}$, $T \le 6^\circ\text{C}$).
- `[storage]`: DuckDB analytical database persistence and Parquet export.
- `[hardware]`: NVIDIA CUDA GPU hardware acceleration and automatic CPU fallback.
- `[visualization]`: Interactive HTML5 Leaflet map export and dashboard options.
- `[paths]`: File system directories (`inputs`, `outputs`) and pseudorandom seed.

### Julia Configuration APIs
```julia
using ParticleTracking

# Load centralized configuration
cfg = load_configuration("inputs/ParticalTracking.config")

# Convert configuration dict to validated HydrodynamicOptions
opts = configuration_to_options(cfg, scenario = :ssp585, n_particles = 250)

# Save active options back to a configuration file
save_configuration(options_to_configuration(opts), "inputs/custom_run.config")
```

---

## 15. Administrative Crab Fishing Area (CFA) Boundary Polygons

The regional modeling platform ingests official administrative management boundaries defined in CSV/DAT format (`lon,lat`):
- `inputs/cfa4x.dat`: Southwest Nova Scotia (CFA 4X, 30 vertices).
- `inputs/cfanorth.dat`: Eastern Cape Breton & Glace Bay (CFA North / 20–22, 17 vertices).
- `inputs/cfasouth.dat`: Middle Scotian Shelf & Halifax (CFA South / 23–24, 28 vertices).

### Point-in-Polygon Classification
Spatial connectivity calculations use the Jordan Curve (ray-casting) theorem via `point_in_polygon(lon, lat, poly_lons, poly_lats)` to assign particle coordinates $(\lambda_p, \phi_p)$ directly to irregular administrative polygons rather than rectangular bounding boxes.

### Leaflet Visualization
Interactive Leaflet dashboards automatically load all available `inputs/cfa*.dat` polygons, rendering vector overlays via `L.polygon(stratum.polygon, ...)` with stage-dependent color highlights and management zone tooltips.

---

## 16. References

### Oceanic Dispersion & Empirical Movement
- **Okubo, A.** (1971). Oceanic diffusion diagrams. *Deep Sea Research and Oceanographic Abstracts*, 18(8), 789-802. DOI: [10.1016/0011-7471(71)90046-5](https://doi.org/10.1016/0011-7471(71)90046-5)

### Tidal Dynamics & Shelf Mixing
- **Egbert, G. D., & Erofeeva, S. Y.** (2002). Efficient inverse modeling of barotropic ocean tides. *Journal of Atmospheric and Oceanic Technology*, 19(2), 183-204. DOI: [10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2](https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2)
- **Pugh, D., & Woodworth, P.** (2014). *Sea-Level Science: Understanding Tides, Surges, Tsunamis and Mean Sea-Level Changes*. Cambridge University Press. DOI: [10.1017/CBO9781139235778](https://doi.org/10.1017/CBO9781139235778)
- **Simpson, J. H., & Hunter, J. R.** (1974). Fronts in the Irish Sea. *Nature*, 250(5465), 404-406. DOI: [10.1038/250404a0](https://doi.org/10.1038/250404a0)

### Climate Forcing & Regional Projections
- **Brickman, D., Wang, Z., & DeTracey, B.** (2018). Variability and trends in the Scotian Shelf and Gulf of Maine region from a high-resolution regional ocean climate model. *Progress in Oceanography*, 164, 49-64. DOI: [10.1016/j.pocean.2018.04.004](https://doi.org/10.1016/j.pocean.2018.04.004)
- **Hobday, A. J., et al.** (2016). A hierarchical approach to defining marine heatwaves. *Progress in Oceanography*, 141, 227-238. DOI: [10.1016/j.pocean.2015.12.014](https://doi.org/10.1016/j.pocean.2015.12.014)
- **Loder, J. W., van der Baaren, A., & Yashayaev, I.** (2015). Climate change trends and projections for the Canadian Northwest Atlantic. *Canadian Technical Report of Hydrography and Ocean Sciences*, 305, 142 pp.
- **O'Neill, B. C., et al.** (2016). The Scenario Model Intercomparison Project (ScenarioMIP) for CMIP6. *Geoscientific Model Development*, 9(9), 3461-3482. DOI: [10.5194/gmd-9-3461-2016](https://doi.org/10.5194/gmd-9-3461-2016)
- **Saba, V. S., et al.** (2016). Enhanced warming of the Northwest Atlantic Ocean under climate change. *Journal of Geophysical Research: Oceans*, 121(1), 118-132. DOI: [10.1002/2015JC011346](https://doi.org/10.1002/2015JC011346)

### Ocean Hydrodynamics & Numerical Methods
- **Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A.** (2007). *Spectral Methods: Evolution to Complex Geometries and Applications to Fluid Dynamics*. Springer-Verlag, Berlin.
- **Courant, R., Friedrichs, K., & Lewy, H.** (1928). Über die partiellen Differenzengleichungen der mathematischen Physik. *Mathematische Annalen*, 100(1), 32-74. DOI: [10.1007/BF01448839](https://doi.org/10.1007/BF01448839)
- **Marshall, J., Adcroft, A., Hill, C., Perelman, L., & Heisey, C.** (1997). A finite-volume, incompressible Navier Stokes model for studies of the ocean on parallel computers. *Journal of Geophysical Research: Oceans*, 102(C3), 5753-5766. DOI: [10.1029/96JC02775](https://doi.org/10.1029/96JC02775)
- **Ramadhan, A., Marshall, J., Hill, C., Campin, J. M., Bischoff, T., & Wagner, G. L.** (2020). Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs. *Journal of Open Source Software*, 5(53), 2018. DOI: [10.21105/joss.02018](https://doi.org/10.21105/joss.02018)
- **Vallis, G. K.** (2017). *Atmospheric and Oceanic Fluid Dynamics: Fundamentals and Large-Scale Circulation*. 2nd Edition. Cambridge University Press, Cambridge. DOI: [10.1017/9781107588417](https://doi.org/10.1017/9781107588417)
- **Verzicco, R.** (2023). Immersed boundary methods for ocean modeling. *Annual Review of Fluid Mechanics*, 55, 305-333. DOI: [10.1146/annurev-fluid-030322-040713](https://doi.org/10.1146/annurev-fluid-030322-040713)

### Atmospheric Forcing & Open Datasets
- **GEBCO Compilation Group** (2023). *GEBCO 2023 Grid*. DOI: [10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e](https://doi.org/10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e)
- **Large, W. G., & Pond, S.** (1981). Open ocean momentum flux measurements in moderate to strong winds. *Journal of Physical Oceanography*, 11(3), 324-336. DOI: [10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2](https://doi.org/10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2)
- **NOAA National Centers for Environmental Information** (2022). *NOAA ETOPO 2022 15 Arc-Second Global Relief Model*. NOAA NCEI. DOI: [10.25921/fd1h-fy81](https://doi.org/10.25921/fd1h-fy81)
- **Simons, R. A.** (2019). *ERDDAP: The Environmental Research Division's Data Access Program*. NOAA CoastWatch / SWFSC.
- **Wu, J.** (1982). Wind-stress coefficients over sea surface from breeze to hurricane. *Journal of Geophysical Research: Oceans*, 87(C12), 9704-9706. DOI: [10.1029/JC087iC12p09704](https://doi.org/10.1029/JC087iC12p09704)
- **Zhang, H.-M., Bates, J. J., & Reynolds, R. W.** (2006). Assessment of composite global sampling: Sea surface wind speed. *Geophysical Research Letters*, 33(17), L17714. DOI: [10.1029/2006GL027086](https://doi.org/10.1029/2006GL027086)

### Snow Crab Larval Ecology & Lagrangian Transport
- **Epifanio, C. E., & Cohen, J. H.** (2016). Behavioral adaptations in larvae of brachyuran crabs: a review. *Journal of Experimental Marine Biology and Ecology*, 482, 85-105. DOI: [10.1016/j.jembe.2016.05.006](https://doi.org/10.1016/j.jembe.2016.05.006)
- **Incze, L. S., Armstrong, D. A., & Smith, S. L.** (1987). Abundance of filter-feeding and pelagic stages of crab larvae in the southeastern Bering Sea. *Marine Biology*, 95(2), 195-200. DOI: [10.1007/BF00409006](https://doi.org/10.1007/BF00409006)
- **Kloeden, P. E., & Platen, E.** (1992). *Numerical Solution of Stochastic Differential Equations*. Springer-Verlag, Berlin. DOI: [10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)
- **Kuhn, P. S., & Choi, J. S.** (2011). Influence of temperature on embryo incubation and larval development in snow crab (*Chionoecetes opilio*). *Fisheries Research*, 107(1-3), 81-87. DOI: [10.1016/j.fishres.2010.10.011](https://doi.org/10.1016/j.fishres.2010.10.011)
- **Lovrich, G. A., Sainte-Marie, B., & Smith, B. D.** (1995). Depth distribution and seasonal movements of *Chionoecetes opilio* (Brachyura: Majidae) in Baie Sainte-Marguerite, Gulf of Saint Lawrence. *Canadian Journal of Fisheries and Aquatic Sciences*, 52(4), 903-913. DOI: [10.1139/f95-090](https://doi.org/10.1139/f95-090)
- **North, E. W., Gallego, A., & Petitgas, P. (Eds.)** (2009). Manual of recommended practices for modelling physical - biological interactions during fish early life. *ICES Cooperative Research Report*, No. 295, 111 pp.
- **Sainte-Marie, G., & Sainte-Marie, B.** (1999). Growth, developmental stages, and vertical distribution of snow crab larvae (*Chionoecetes opilio*) in the northwestern Gulf of St. Lawrence. *Canadian Journal of Fisheries and Aquatic Sciences*, 56(11), 2181-2193. DOI: [10.1139/f99-151](https://doi.org/10.1139/f99-151)