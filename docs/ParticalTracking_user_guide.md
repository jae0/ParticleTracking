# ParticleTracking.jl: Comprehensive User Guide & Reference 

**A High-Performance Biophysical Ocean Modeling, Lagrangian Particle Tracking, and Demographic Population Connectivity Framework**

---

## 1. Foundations

`ParticleTracking.jl` is a high-performance Julia biophysical modeling system designed to simulate 3D ocean hydrodynamics, Lagrangian particle transport, individual-based larval physiology, and macro-regional population connectivity across continental shelf environments.

Originally parameterized for the **Scotian Shelf snow crab (*Chionoecetes opilio*)** ecosystem in the Northwest Atlantic (Fisheries and Oceans Canada Crab Fishing Areas CFAs 20–22, 23–24, and 4X), the framework is fully generalizable to any marine species or regional ocean basin worldwide.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     ParticleTracking.jl End-to-End System                   │
└─────────────────────────────────────────────────────────────────────────────┘
  1. Open Data Ingestion (NOAA ERDDAP ETOPO/GEBCO & Blended Sea Winds)
                                │
                                ▼
  2. 3D Hydrodynamic Model (Oceananigans.jl on CUDA GPU / Multi-Threaded CPU)
     - Spherical Curvilinear Grid + Immersed Boundary Topography
     - M2/S2 Tidal Harmonic Body Forcing + Coriolis + Seawater Buoyancy
     - CMIP6 Climate Scenarios (Historical, SSP1-2.6, SSP2-4.5, SSP5-8.5, MHW)
                                │
                                ▼
  3. Individual-Based Lagrangian Particle Tracking (Euler-Maruyama SDE)
     - Marine Bathymetric Rejection Sampling (Strict Water Placement >= 100m)
     - Stage-Dependent Diel Vertical Migration (DVM) Active Swimming
     - In Situ Thermal Degree-Day Ontogenetic Molting (Zoea I -> II -> Megalopa)
     - Cumulative Temperature-Dependent PLD & Thermal Stress Mortality
     - Cold Intermediate Layer (CIL) Nursery Settlement Criteria
                                │
                                ▼
  4. Spatial Demographics & Administrative CFA Polygon Ray-Casting
     - CFA North (20-22), CFA South (23-24), CFA 4X (SW Nova Scotia)
     - Transition Probability Matrices (P_ij), Self-Retention & Export
     - Taylor Dispersion Empirical Velocity & Turbulent Diffusivity Fields
                                │
                                ▼
  5. Multi-Layer NetCDF / JLD2 / DuckDB Analytical Engine
     - Embedded Columnar DuckDB Database (Outputs/particle_tracking.duckdb)
     - Multi-Scenario Intercomparison & Bayesian / Ensemble Model Averaging
     - Zero-Copy Apache Parquet Export (BSTM / R / Python Interoperability)
                                │
                                ▼
  6. Visualizations & Standalone Interactive Leaflet.js Dashboard
     - 8 Publication-Ready CairoMakie Diagnostic Charts
     - Interactive HTML5 Map with Playback Scrubber, Live Telemetry & CFA Polygons
```

---

## 2. Software Architecture & Module Structure

The project follows modular software engineering standards. Core algorithms are encapsulated in `src/`, with shared exports declared in `src/ParticleTracking.jl` which is the single entry point for all module loading.


```
ParticleTracking/
├── Project.toml                      # Project dependencies & UUID declarations
├── inputs/
│   ├── ParticleTracking.config       # Master centralized configuration file
│   ├── coastline.dat                 # Nova Scotia coastline boundary polygon vertices
│   ├── cfa4x.dat                     # CFA 4X boundary polygon vertices
│   ├── cfanorth.dat                  # CFA North (20-22) boundary polygon vertices
│   └── cfasouth.dat                  # CFA South (23-24) boundary polygon vertices
├── outputs/                          # Generated NetCDF, JLD2, DuckDB & figures
├── docs/
│   ├── ParticalTracking_user_guide.md # Comprehensive technical manual
│   └── snow_crab_larval_connectivity_paper.md # Peer-reviewed paper manuscript
├── src/
│   ├── ParticleTracking.jl           # Main module: entry point, exports & all sub-includes
│   ├── configuration.jl              # TOML config manager & HydrodynamicOptions struct
│   ├── architecture.jl               # NVIDIA CUDA GPU / CPU hardware resolution
│   ├── open_data.jl                  # NOAA ERDDAP data fetcher & Large & Pond drag law
│   ├── synthetic_data.jl             # Idealized bathymetry & wind benchmark generator
│   ├── grid_bathymetry.jl            # Spherical grid & cut-cell immersed boundary setup
│   ├── hydrodynamic_model.jl         # Hydrostatic primitive equation model setup & heat flux
│   ├── tides.jl                      # Astronomical tidal forcing & Simpson-Hunter metric
│   ├── climate_scenarios.jl          # CMIP6 anomalies, PLD formulas, thermal mortality
│   ├── simulation.jl                 # Oceananigans time integration & adaptive CFL
│   ├── larval_behavior.jl            # DVM swimming, BBL shear, sinking, drift, tracking
│   ├── empirical_analysis.jl         # Taylor dispersion, CFA polygons & connectivity
│   ├── storage_duckdb.jl             # DuckDB analytical backend & ensemble averaging
│   └── visualization.jl              # CairoMakie plots & Leaflet interactive map
├── test/
│   └── runtests.jl                   # Comprehensive 17-testset unit test suite (618 unit tests)
├── ParticleTrackingRun.jl            # Production-grade flag-driven CLI application
├── ParticleTrackingRun.md            # Execution guide & CLI quickstart
└── ParticleTracking issues.md        # Quality assurance audit log & bug resolutions
```

### Module Responsibilities & Key Exports

| Source File                                                                                            | Primary Purpose                                               | Exported Functions / Types                                                                                                                                                                                                                                                                                                |
| :-------------------------------------------------------------------------------------------------------| :--------------------------------------------------------------| :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`src/configuration.jl`](file:///c:/home/jae/projects/ParticleTracking/src/configuration.jl)           | Centralized TOML configuration parsing and options conversion | `HydrodynamicOptions`, `load_configuration`, `save_configuration`, `configuration_to_options`, `options_to_configuration`, `get_default_configuration`, `find_default_config_path`                                                                                                                                        |
| [`src/architecture.jl`](file:///c:/home/jae/projects/ParticleTracking/src/architecture.jl)             | Hardware platform detection and device allocation             | `resolve_architecture`                                                                                                                                                                                                                                                                                                    |
| [`src/open_data.jl`](file:///c:/home/jae/projects/ParticleTracking/src/open_data.jl)                   | NOAA ERDDAP bathymetry & wind ingestion, drag laws            | `fetch_open_bathymetry`, `fetch_open_surface_winds`, `wind_speed_to_kinematic_stress`, `regrid_2d_field`                                                                                                                                                                                                                  |
| [`src/synthetic_data.jl`](file:///c:/home/jae/projects/ParticleTracking/src/synthetic_data.jl)         | Synthetic idealized benchmarks & NetCDF inspection            | `download_sample_data`, `inspect_netcdf`, `generate_synthetic_bathymetry`, `generate_synthetic_forcing`                                                                                                                                                                                                                   |
| [`src/grid_bathymetry.jl`](file:///c:/home/jae/projects/ParticleTracking/src/grid_bathymetry.jl)       | Spherical grid & cut-cell immersed boundary setup             | `build_shelf_grid`, `load_bathymetry_from_netcdf`, `apply_immersed_bathymetry`, `get_bathymetry_interpolator`, `load_coastline_polygons`, `is_point_on_land`                                                                                                                                                            |
| [`src/hydrodynamic_model.jl`](file:///c:/home/jae/projects/ParticleTracking/src/hydrodynamic_model.jl) | Hydrostatic Boussinesq primitive equation & heat flux         | `build_hydrodynamic_model`, `set_initial_stratification!`                                                                                                                                                                                                                                                                 |
| [`src/tides.jl`](file:///c:/home/jae/projects/ParticleTracking/src/tides.jl)                           | Astronomical tidal forcing, spring-neap & Simpson-Hunter      | `build_tidal_body_forcing`, `tidal_velocity_vector`, `get_tidal_frequency`, `simpson_hunter_parameter`                                                                                                                                                                                                                     |
| [`src/climate_scenarios.jl`](file:///c:/home/jae/projects/ParticleTracking/src/climate_scenarios.jl)   | CMIP6 downscaling, PLD formulas, thermal mortality            | `ClimateScenarioDelta`, `get_climate_scenario`, `apply_climate_anomaly_to_stratification`, `temperature_dependent_pld`, `larval_thermal_mortality_rate`                                                                                                                                                                     |
| [`src/simulation.jl`](file:///c:/home/jae/projects/ParticleTracking/src/simulation.jl)                 | Time integration, adaptive CFL, 4D output interpolation       | `setup_hydrodynamic_simulation`, `run_hydrodynamic_simulation!`, `create_flow_interpolator_from_jld2`, `compute_advective_cfl`                                                                                                                                                                                            |
| [`src/larval_behavior.jl`](file:///c:/home/jae/projects/ParticleTracking/src/larval_behavior.jl)       | DVM swimming, BBL shear, sinking, drift, tracking             | `initialize_larval_particles`, `larval_ascent_velocity`, `diel_vertical_migration_velocity`, `superpose_tidal_velocity`, `bbl_velocity_factor`, `larval_passive_sinking_velocity`, `update_larval_stage`, `evaluate_settlement_suitability`, `larval_transport_step`, `track_larval_cohort`, `canonicalize_trajectories`                  |
| [`src/empirical_analysis.jl`](file:///c:/home/jae/projects/ParticleTracking/src/empirical_analysis.jl) | Taylor dispersion, CFA polygons, recruitment connectivity     | `estimate_empirical_movement`, `compute_gridded_recruitment_metrics`, `compute_gridded_thermal_metrics`, `point_in_polygon`, `load_cfa_polygons`, `compute_empirical_connectivity`, `connectivity_transitions`, `export_larval_dispersal_netcdf`, `export_larval_dispersal_jld2`                                     |
| [`src/storage_duckdb.jl`](file:///c:/home/jae/projects/ParticleTracking/src/storage_duckdb.jl)         | DuckDB analytical backend & ensemble averaging                | `open_duckdb_storage`, `close_duckdb_storage`, `save_simulation_run!`, `load_run_configuration`, `list_simulation_runs`, `load_trajectories_df`, `load_connectivity_matrix`, `compare_scenarios`, `compute_ensemble_model_average`, `export_duckdb_to_parquet`                                                            |
| [`src/visualization.jl`](file:///c:/home/jae/projects/ParticleTracking/src/visualization.jl)           | CairoMakie figures & interactive Leaflet HTML dashboard       | `plot_particle_trajectories`, `plot_dvm_depth_profiles`, `plot_larval_dispersal_density`, `plot_empirical_movement_field`, `plot_connectivity_matrix`, `plot_thermal_exposure_map`, `plot_recruitment_summary`, `plot_climate_scenario_comparison`, `export_interactive_tracks_html`, `plot_interactive_trajectories_map` |

---

## 3. Biophysical Formulations

### 3.1 3D Hydrodynamic Circulation & Boussinesq Equations

The regional hydrodynamic model ([Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl)) integrates the 3D hydrostatic Boussinesq primitive equations on a spherical latitude-longitude coordinate grid $(\lambda, \phi, z)$:

$$\frac{\partial \boldsymbol{u}_h}{\partial t} + (\boldsymbol{u} \cdot \nabla) \boldsymbol{u}_h + f \hat{\boldsymbol{k}} \times \boldsymbol{u}_h = -\frac{1}{\rho_0} \nabla_h p + \nabla_h \cdot (\nu_h \nabla_h \boldsymbol{u}_h) + \frac{\partial}{\partial z}\left(\nu_v \frac{\partial \boldsymbol{u}_h}{\partial z}\right) + \boldsymbol{F}_{\text{tide}}$$

$$\frac{\partial p}{\partial z} = -\rho g = b \rho_0$$

$$\nabla \cdot \boldsymbol{u} = \frac{1}{R_E \cos\phi}\frac{\partial u}{\partial \lambda} + \frac{1}{R_E}\frac{\partial v}{\partial \phi} + \frac{\partial w}{\partial z} = 0$$

$$\frac{\partial T}{\partial t} + \boldsymbol{u} \cdot \nabla T = \nabla \cdot (\boldsymbol{\kappa}_T \nabla T)$$

$$\frac{\partial S}{\partial t} + \boldsymbol{u} \cdot \nabla S = \nabla \cdot (\boldsymbol{\kappa}_S \nabla S)$$

where:
- $\boldsymbol{u}_h = (u, v)$ is the horizontal velocity vector (zonal, meridional);
- $w$ is the diagnostic vertical velocity;
- $f = 2\Omega \sin\phi$ is the Coriolis parameter evaluated at latitude $\phi$;
- $b = -g(\rho - \rho_0)/\rho_0$ is buoyancy computed via `SeawaterBuoyancy`;
- $\boldsymbol{\kappa}_T, \boldsymbol{\kappa}_S$ are turbulent thermal and salt diffusivities;
- $\boldsymbol{\kappa}_T, \boldsymbol{\kappa}_S$ are turbulent thermal and salt diffusivities;
- $\boldsymbol{F}_{\text{tide}}$ is the astronomical tidal momentum body forcing;
- Bottom momentum drag incorporates linear Rayleigh and quadratic drag: $\boldsymbol{F}_{\text{drag}} = -(r_{\text{drag}} + C_d |\boldsymbol{u}_h|) \boldsymbol{u}_h$.

---

### 3.2 Atmospheric Wind Drag & Air–Sea Surface Heat Flux

Surface momentum boundary conditions are parameterized from 10-meter wind vectors $(u_{10}, v_{10})$ via the empirical non-linear drag law of Large & Pond (1981) modified by Wu (1982):

$$\boldsymbol{\tau}_{\text{kinematic}} = \frac{\boldsymbol{\tau}}{\rho_{\text{water}}} = \frac{\rho_{\text{air}}}{\rho_{\text{water}}} C_d |\boldsymbol{u}_{10}| \boldsymbol{u}_{10}$$

$$C_d(U_{10}) = \begin{cases} 1.2 \times 10^{-3}, & U_{10} \le 11.0\text{ m s}^{-1} \\ (0.49 + 0.065 U_{10}) \times 10^{-3}, & U_{10} > 11.0\text{ m s}^{-1} \end{cases}$$

Net atmospheric heat exchange $Q_{\text{net}}$ ($\text{W m}^{-2}$, positive downward warming) enters as a kinematic temperature flux top boundary condition:

$$J_T = -\frac{Q_{\text{net}}}{\rho_0 c_p} \quad [^\circ\text{C}\text{ m s}^{-1}]$$

where volumetric heat capacity $\rho_0 c_p \approx 1025.0 \times 3990.0 \approx 4.09 \times 10^6\text{ J m}^{-3}\;^\circ\text{C}^{-1}$.

---

### 3.3 Astronomical Tidal Body Forcing, Spring–Neap Cycles & Mixing Fronts

Barotropic tidal harmonics (e.g. $M_2$ semi-diurnal, $\tau = 12.42\text{ h}$; $S_2$ solar semi-diurnal, $\tau = 12.00\text{ h}$) are introduced via bottom-drag-compensated momentum body forcing accelerations $\boldsymbol{F}_{\text{tide}} = (F_u, F_v)$ (Egbert & Erofeeva, 2002):

$$F_u(\lambda, \phi, z, t) = \sum_{k} U_k \sqrt{\omega_k^2 + r_{\text{drag}}^2} \cos(\omega_k t + \phi_{k,u})$$

$$F_v(\lambda, \phi, z, t) = \sum_{k} V_k \sqrt{\omega_k^2 + r_{\text{drag}}^2} \sin(\omega_k t + \phi_{k,v})$$

Superposition of $M_2$ and $S_2$ generates the classic 14.77-day spring-neap beat envelope:

$$T_{\text{spring-neap}} = \frac{T_{M2} T_{S2}}{|T_{M2} - T_{S2}|} \approx 14.77\text{ days}$$

To diagnose tidal mixing fronts and boundary layer separation over shallow banks under tidal and atmospheric forcing, we compute the generalized Simpson-Hunter parameter $\chi$ (Garrett, Keeley & Greenberg 1978; Loder & Greenberg 1986):

$$\chi = \log_{10}\left( \frac{h}{U_{\text{tide}}^3 + \gamma U_{\text{wind}}^3} \right)$$

- $\chi < 1.5$: Vertically well-mixed water column (energetic bank crests).
- $1.5 \le \chi \le 2.0$: Tidal mixing front (biophysical retention zone).
- $\chi > 2.0$: Stratified water column (intermediate shelf basins).

---

### 3.4 CMIP6 Regional Climate Scenarios & Thermal Structure

Climate change perturbations are implemented through vertical anomaly profiles downscaled for the Scotian Shelf (Brickman et al., 2018; Saba et al., 2016):

$$\Delta T(z) = \Delta T_{\text{deep}} + (\Delta \text{SST} - \Delta T_{\text{deep}}) \exp\left( \frac{z}{H_{\text{mix}}} \right)$$

$$\Delta S(z) = \Delta S_{\text{surface}} \exp\left( \frac{z}{H_{\text{mix}}} \right)$$

where $H_{\text{mix}} = 30.0\text{ m}$ is the epipelagic mixed layer depth.

| Scenario     | Description                             | $\Delta \text{SST}$ (°C) | $\Delta T_{\text{deep}}$ (°C) | $\Delta S_{\text{surface}}$ (PSU) | Wind Multiplier |
| :-------------| :----------------------------------------| :-------------------------| :------------------------------| :----------------------------------| :----------------|
| `historical` | Historical baseline climatology         | $0.0$                    | $0.0$                         | $0.0$                             | $1.00 \times$   |
| `ssp126`     | CMIP6 SSP1-2.6 (2050 Low Emissions)     | $+1.1$                   | $+0.5$                        | $-0.25$                           | $1.05 \times$   |
| `ssp245`     | CMIP6 SSP2-4.5 (2050 Intermediate)      | $+1.8$                   | $+0.9$                        | $-0.45$                           | $1.10 \times$   |
| `ssp585`     | CMIP6 SSP5-8.5 (2050 High Emissions)    | $+3.5$                   | $+1.9$                        | $-0.85$                           | $1.20 \times$   |
| `mhw`        | Marine Heatwave (Transient Category IV) | $+3.5$                   | $+0.2$                        | $-0.10$                           | $0.80 \times$   |

---

### 3.5 Individual-Based Lagrangian Particle Tracking SDE

Larval trajectories are governed by the stochastic Langevin equation integrated via the Euler-Maruyama scheme (North et al., 2009), enhanced with bottom boundary layer shear, passive gravitational sinking, and diffusive pseudo-drift (Hunter 1993; Visser 1997):

$$\boldsymbol{x}_h^{n+1} = \boldsymbol{x}_h^n + \left[ f_{\text{bbl}}(z^n) \boldsymbol{u}_h(\boldsymbol{x}^n, t^n) + \boldsymbol{u}_{\text{tide}}(\boldsymbol{x}^n, t^n) \right] \Delta t + \sqrt{2 \kappa_h \Delta t} \, \boldsymbol{\xi}_h^n$$

$$z^{n+1} = z^n + \left[ w(\boldsymbol{x}^n, t^n) + w_{\text{swim}}(z^n, t^n) + w_{\text{sink}}(\text{stage}) + \frac{\partial \kappa_v}{\partial z} \right] \Delta t + \sqrt{2 \kappa_v\left(z^n + \frac{1}{2}\frac{\partial \kappa_v}{\partial z}\Delta t\right) \Delta t} \, \xi_z^n$$

where:
- $f_{\text{bbl}}(z) = \text{clamp}\left( \frac{\ln(\max(z_0, z - z_{\text{bed}}) / z_0)}{\ln(h_{\text{bbl}} / z_0)}, 0.0, 1.0 \right)$ attenuates horizontal velocity logarithmically near the seabed ($h_{\text{bbl}} = 10\text{ m}, z_0 = 10^{-3}\text{ m}$);
- $w_{\text{sink}}$ is passive gravitational settling velocity (Zoea I: $-0.5\text{ mm s}^{-1}$, Zoea II: $-1.0\text{ mm s}^{-1}$, Megalopa: $-2.5\text{ mm s}^{-1}$);
- $\frac{\partial \kappa_v}{\partial z}$ is the Visser (1997) deterministic drift correction preventing unphysical particle trapping in low-diffusivity pycnoclines;
- Vertical boundaries are absorbing: $z^{n+1} \in [z_{\text{bed}}, 0.0]$;
- Shorelines apply polygon ray-casting with alongshore tangential slip.

---

### 3.6 Larval Biological Behaviors: Benthic Release, Ascent, DVM, Molting & Mortality

#### A. Benthic Larval Release & Directed Post-Hatch Vertical Ascent
Adult female snow crabs dwell strictly on the benthic continental shelf floor. Ovigerous
females release egg clutches directly into the near-bottom boundary layer
($z_{\text{init}} \in [z_{\text{bed}} + 0.5, z_{\text{bed}} + 3.0]\text{ m}$) during
spring (Lovrich et al., 1995; Sainte-Marie & Sainte-Marie, 1999). Freshly hatched
Stage I zoeae exhibit high swimming motility characterized by strong negative geotaxis
and positive phototaxis (Sulkin, 1984; Forward, 1988), initiating an active vertical
ascent through the water column toward the epipelagic surface mixed layer
($z_{\text{target}} = -10.0\text{ m}$):

$$w_{\text{ascent}}(z) = w_{\text{ascent,max}} \tanh\left( \frac{z_{\text{target}} - z}{L_{\text{relax}}} \right)$$

where $w_{\text{ascent,max}} = 10.0\text{ mm s}^{-1}$ ($0.010\text{ m s}^{-1}$) and
$L_{\text{relax}} = 10.0\text{ m}$. This directed upward swimming velocity readily
overcomes passive gravitational settling ($w_{\text{sink}} = -0.5\text{ mm s}^{-1}$),
allowing larvae to traverse the deep Cold Intermediate Layer and reach the surface
mixed layer within $\approx 2\text{--}4\text{ hours}$. Once larvae attain the surface
mixed layer ($z \ge z_{\text{target}}$) or exceed the ascent time limit, they
transition smoothly into stage-specific Diel Vertical Migration (DVM).

| Developmental Stage    | Daytime Depth ($z_{\text{day}}$) | Nighttime Depth ($z_{\text{night}}$) | Max Speed ($w_{\max}$) | Ecological Niche                                |
| :-----------------------| :---------------------------------| :-------------------------------------| :-----------------------| :-----------------------------------------------|
| **Zoea I (Ascent)**    | Surface Target ($-10.0\text{ m}$) | Surface Target ($-10.0\text{ m}$)     | $10.0\text{ mm s}^{-1}$| Post-hatch benthic ascent                       |
| **Zoea I (DVM)**       | $-50.0\text{ m}$                 | $-10.0\text{ m}$                     | $5.0\text{ mm s}^{-1}$ | Epipelagic warm mixed layer ($T > 4^\circ\text{C}$)|
| **Zoea II**            | $-55.0\text{ m}$                 | $-8.0\text{ m}$                      | $6.0\text{ mm s}^{-1}$ | Pycnocline feeding / shallow nighttime refuge   |
| **Megalopa**           | $-120.0\text{ m}$                | $-60.0\text{ m}$                     | $7.5\text{ mm s}^{-1}$ | CIL benthic staging & settlement search         |
| **Instar I (Recruit)** | Benthic Bed ($z_{\text{bed}}$)   | Benthic Bed ($z_{\text{bed}}$)       | $0.0\text{ mm s}^{-1}$ | Substrate-anchored benthic recruit              |

#### B. Context-Dependent Diel Vertical Migration (DVM)
Larvae modulate their vertical swimming toward stage-specific diurnal target depths, attenuated by turbidity $\alpha_{\text{turb}}$ and constrained by the Cold Intermediate Layer (CIL) thermocline boundary $z_{\text{cil}}$ (Incze et al., 1987; Sainte-Marie & Sainte-Marie, 1999):

$$z_{\text{amp}} = \frac{z_{\text{night}} - z_{\text{day}}}{2} \cdot \text{clamp}(\alpha_{\text{turb}}, 0.0, 1.0)$$

$$z_{\text{target}}(t) = \frac{z_{\text{day}} + z_{\text{night}}}{2} - z_{\text{amp}} \cos\left( \frac{2\pi t}{86400} \right)$$

$$w_{\text{swim}}(z, t) = w_{\max} \tanh\left( \frac{z_{\text{target}}(t) - z}{L_{\text{relax}}} \right)$$

#### C. In Situ Thermal Degree-Days & Calibrated Ontogenetic Molting
Thermal degree-days accumulate along 3D trajectories relative to the physiological baseline $T_0 = -1.5^\circ\text{C}$ (Kuhn & Choi, 2011):

$$DD(t) = \int_0^t \max\left( 0.0, T(x(\tau), y(\tau), z(\tau), \tau) - T_0 \right) d\tau \quad (T_0 = -1.5^\circ\text{C})$$

- **Zoea I $\to$ Zoea II**: $DD \ge 65.0^\circ\text{C} \cdot \text{days}$
- **Zoea II $\to$ Megalopa**: $DD \ge 130.0^\circ\text{C} \cdot \text{days}$
- **Megalopa $\to$ Instar I (Competent Settlement)**: $DD \ge 200.0^\circ\text{C} \cdot \text{days}$

#### D. Temperature-Dependent PLD & Thermal Stress Mortality
Pelagic Larval Duration scales with temperature: $\text{PLD}(T) = 135.0 \cdot (T - T_0)^{-0.75}\text{ days}$.

Instantaneous daily mortality rate $\mu(T)$ includes baseline mortality, stage-specific scaling (Zoea I: $1.1\times$, Zoea II: $1.0\times$, Megalopa: $0.8\times$), exponential warm stress, and linear cold stress:

$$\mu(T) = \mu_0 s_{\text{stage}} \exp\left( \beta \max(0, T - T_{\text{crit}}) \right) + \mu_{\text{cold}} \max(0, T_{\text{cold}} - T)$$

where $T_{\text{crit}} = 7.0^\circ\text{C}$, $\beta = 0.35^\circ\text{C}^{-1}$, and $T_{\text{cold}} = -1.5^\circ\text{C}$. Stage survival is recorded independently across each ontogenetic molt.

#### E. Benthic Nursery Settlement Suitability with Subgrid Tidal Filtering
In energetic tidal regions (CFAs 20–22), semi-diurnal $M_2$ tides displace the thermocline by $\pm 2^\circ\text{C}$ over 12.42 hours. A low-pass exponential filter extracts the true benthic thermal regime:

$$\bar{T}_{\text{bed}}(t + \Delta t) = (1 - \alpha) \bar{T}_{\text{bed}}(t) + \alpha T_{\text{bed}}(t), \quad \alpha = \text{clamp}\left(\frac{\Delta t}{44712}, 0.005, 1.0\right)$$

Competent Megalopae settle when:
1. **Depth Suitability**: $-250.0\text{ m} \le z_{\text{bed}} \le -50.0\text{ m}$;
2. **Thermal Suitability**: Tidally averaged bottom temperature $\bar{T}_{\text{bed}} \le 6.0^\circ\text{C}$.

---

### 3.7 Administrative CFA Boundary Polygons & Ray-Casting Algorithm

Spatial classification utilizes true administrative boundary vertices from `inputs/cfa*.dat`:
- **CFA North (CFAs 20–22)**: Eastern Cape Breton & Glace Bay (`inputs/cfanorth.dat`, 17 vertices).
- **CFA South (CFAs 23–24)**: Middle Scotian Shelf & Halifax (`inputs/cfasouth.dat`, 28 vertices).
- **CFA 4X**: Southwest Nova Scotia & Browns Bank (`inputs/cfa4x.dat`, 30 vertices).

#### Jordan Curve (Ray-Casting) Theorem
For any particle coordinate $(\lambda_p, \phi_p)$, a horizontal ray is cast eastward: $\{\lambda \ge \lambda_p, \phi = \phi_p\}$. For each polygon edge from $(\lambda_i, \phi_i)$ to $(\lambda_{i+1}, \phi_{i+1})$, an intersection occurs if $\phi_p \in [\min(\phi_i, \phi_{i+1}), \max(\phi_i, \phi_{i+1}))$ and:

$$\lambda_{\text{int}} = \lambda_i + (\phi_p - \phi_i) \frac{\lambda_{i+1} - \lambda_i}{\phi_{i+1} - \phi_i} > \lambda_p$$

The particle is inside the polygon if the total number of intersections is odd.

---

### 3.8 Demographic Recruitment Connectivity Matrix

The demographic recruitment transition matrix $P_{ij}$ weights individual settlement by cumulative survival probability $S_p(t_{\text{end}})$:

$$P_{ij} = \frac{\sum_{p \in (i \to j)} S_p(t_{\text{end}})}{N_{\text{released}, i}}$$

- **Self-Retention Probability**: $S_i = P_{ii}$ (demographic recruitment retained in natal stratum).
- **Export Probability**: $E_i = \sum_{j \ne i} P_{ij}$ (recruitment subsidizing downstream strata).
- **Pelagic Mortality Loss**: $L_i = 1 - \sum_{j} P_{ij}$ (unsettled or dead larvae).

---

## 4. Centralized Configuration System (`inputs/ParticleTracking.config`)

All user-tunable parameters and physical/biological defaults are centralized in `inputs/ParticleTracking.config` (TOML format).

```toml
[domain]
lon_min = -68.0                     # Western longitude limit (°E)
lon_max = -57.0                     # Eastern longitude limit (°E)
lat_min = 42.0                      # Southern latitude limit (°N)
lat_max = 47.0                      # Northern latitude limit (°N)
z_min = -1000.0                     # Seabed bottom limit (meters)
z_max = 0.0                         # Sea surface limit (meters)

[grid]
nx = 50                             # Zonal grid cell resolution
ny = 50                             # Meridional grid cell resolution
nz = 10                             # Vertical grid layer resolution

[data]
mode = "synthetic"                  # "synthetic" (idealized shelf) or "real" (NOAA ERDDAP)
elevation_offset = -50.0            # Minimum offshore bank depth (meters)
elevation_slope = 1.8               # Continental shelf-to-slope gradient (m/km)
ref_wind_u = 8.5                    # Baseline 10m zonal wind speed (m/s)
ref_wind_v = 3.2                    # Baseline 10m meridional wind speed (m/s)

[tides]
enable_tides = true                 # Astronomical tidal forcing (M2 + S2)
constituents = ["M2", "S2"]         # Tidal constituents (:M2, :S2 for spring-neap envelope)
tidal_u_amp = 0.25                  # Semi-major M2 tidal amplitude (m/s)
tidal_v_amp = 0.12                  # Semi-minor M2 tidal amplitude (m/s)
s2_u_amp = 0.11                     # Semi-major S2 tidal amplitude (m/s)
s2_v_amp = 0.05                     # Semi-minor S2 tidal amplitude (m/s)
tidal_period = 44712.0              # M2 tidal period (seconds, ~12.42 hours)
tidal_phase = 0.0                   # Initial tidal phase (radians)

[climate]
scenario = "historical"             # "historical", "ssp126", "ssp245", "ssp585", "mhw"
projection_year = 2020              # Target climate projection horizon year
baseline_year = 2020                # Baseline historical benchmark year
horizon_year = 2050                 # Standard CMIP6 horizon year

[hydrodynamics]
sim_duration_hours = 120.0          # Hydrodynamic model integration time (hours)
sim_dt_seconds = 120.0              # Initial time step (seconds)
adaptive_cfl = true                 # Dynamic CFL adaptive time-stepping
target_cfl = 0.2                    # Target advective CFL limit
divergence_velocity_limit = 20.0    # Maximum velocity threshold for stability watchdog (m/s)
coriolis_latitude = 44.5            # Reference latitude for Coriolis parameter (degrees N)
surface_heat_flux = 50.0            # Net surface heat flux (W/m², positive downward warming)

[biology]
n_particles = 500                   # Larval cohort size
track_duration_days = 60.0          # Lagrangian drift duration (days)
track_dt_seconds = 300.0            # Particle tracking time step (seconds)
min_seabed_depth = 100.0            # Minimum release water depth (meters)
diffusivity_h = 10.0                # Horizontal eddy diffusivity (m^2/s)
diffusivity_v = 1e-4                # Vertical eddy diffusivity (m^2/s)
enable_bbl = true                   # Logarithmic bottom boundary layer shear
h_bbl = 10.0                        # BBL thickness (meters)
enable_sinking = true               # Stage-dependent passive gravitational sinking
release_depth_mode = "bottom"       # Larval release depth mode ("bottom", "range", "surface")
bottom_release_offset = [0.5, 3.0]  # Release height interval [min, max] above seafloor (m)
enable_initial_ascent = true        # Active upward swimming post-hatch toward surface
ascent_speed = 0.010                # Initial post-hatch ascent speed (m/s, ~10 mm/s)
ascent_target_depth = -10.0         # Target epipelagic depth for ascent completion (meters)

[dvm]
enable_dvm = true                   # Stage-specific Diel Vertical Migration
zoea1_day_depth = -50.0             # Zoea I daytime target depth (meters)
zoea1_night_depth = -10.0           # Zoea I nighttime target depth (meters)
zoea2_day_depth = -55.0             # Zoea II daytime target depth (meters, empirical niche)
zoea2_night_depth = -8.0            # Zoea II nighttime target depth (meters, shallow refuge)
megalopa_day_depth = -120.0         # Megalopa daytime target depth (meters)
megalopa_night_depth = -60.0        # Megalopa nighttime target depth (meters)
turbidity_attenuation = 1.0         # Turbidity vertical migration attenuation factor

[molting_and_settlement]
enable_molting = true               # Degree-day molting and thermal mortality
t_base = -1.5                       # Base physiological development temperature (°C)
dd_zoea1_to_zoea2 = 65.0            # Zoea I -> Zoea II molt threshold (DD)
dd_zoea2_to_megalopa = 130.0        # Zoea II -> Megalopa molt threshold (DD)
dd_megalopa_to_settle = 200.0       # Megalopa settlement competency threshold (DD)
mortality_base = 0.02               # Natural baseline mortality (day^-1)
mortality_thermal_threshold = 7.0   # Thermal stress threshold (°C)
mortality_thermal_sensitivity = 0.015 # Excess thermal mortality rate
settlement_min_depth = -250.0       # Nursery maximum depth limit (meters)
settlement_max_depth = -50.0        # Nursery minimum depth limit (meters)
settlement_max_temp = 6.0           # Maximum benthic temperature for settlement (°C)

[storage]
enable_duckdb = true                # Persist all simulation runs into DuckDB
duckdb_path = "outputs/particle_tracking.duckdb" # DuckDB file path
export_parquet = false              # Zero-copy export to Apache Parquet

[hardware]
use_gpu = false                     # NVIDIA CUDA GPU hardware acceleration
fallback_to_cpu = true              # Fallback to multi-threaded CPU if GPU is absent

[visualization]
interactive_map = true              # Export interactive HTML5 Leaflet map
title = "Scotian Shelf Snow Crab Larval Dispersal & Demographic Connectivity"

[paths]
output_dir = "outputs"              # Target directory for artifacts and figures
input_dir = "inputs"                # Directory for bathymetry, winds, polygons
seed = 42                           # RNG seed for reproducible stochastic runs
```

---

---

## 1. Setup & Package Loading

```julia
cd("c:/home/jae/work/ParticleTracking")  # your work directory

# Load module and all dependencies (this may require multiple runs and Julia restarts)  
include("c:/home/jae/projects/ParticleTracking/src/ParticleTracking.jl")

using .ParticleTracking # the leading "." is not a typo

# define scope
lon_range=(-68.0, -57.0)
lat_range=(42.0, 47.0)
z_range=(-600.0, -10.0)
grid_size=(100, 100, 20)
time_range=(0.0, 86400.0)
year = "2023-06-01T00:00:00Z"

```

---

## 2. Environmental Data Ingestion (Real vs. Synthetic)

### Fetch Real-World Data from Open Repositories

Download real-world seafloor bathymetry from NOAA CoastWatch / ERDDAP (`etopo180` or GEBCO)
and observed 10m surface vector winds (`erdBSwinds1day`), converting winds to kinematic
surface wind stress via the Large & Pond (1981) parameterization:

```julia


# 1. Fetch real bathymetry from NOAA ERDDAP (Scotian Shelf bounding box)
real_bathy_file = fetch_open_bathymetry(
  lon_range=lon_range,
  lat_range=lat_range,
  output_path="inputs/real_bathymetry.nc",
  dataset_id="etopo180")

# 2. Fetch real 10m surface winds for a target date
real_wind_file = fetch_open_surface_winds(
  lon_range=lon_range,
  lat_range=lat_range,
  time_iso="2023-06-01T00:00:00Z",
  output_path="inputs/real_surface_winds.nc")  # noaa sources not working

real_wind_file = fetch_copernicus_surface_winds(  # requires: pip install copernicusmarine .. api issues
  lon_range=lon_range,
  lat_range=lat_range,
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

### Alternative: Generate Synthetic Benchmark Data

Alternatively, generate synthetic bathymetry and idealized cyclic surface wind forcing:

```julia
synth_bathy_file   = generate_synthetic_bathymetry(
  "inputs/nova_scotia_bathymetry.nc",
  lon_range=lon_range,
  lat_range=lat_range,
  n_lon=grid_size[1], 
  n_lat=grid_size[2],
  inshore_depth=z_range[2],
  shelf_slope=z_range[1])
  
synth_forcing_file = generate_synthetic_forcing(
  "inputs/surface_forcing.nc",
  lon_range=lon_range,
  lat_range=lat_range,
  time_range=time_range,
  n_lon=grid_size[1], 
  n_lat=grid_size[2],
  n_time=24,
  tau_x_amplitude=0.1)
```

---

## 3. Grid & Immersed Boundary Topography

Discretize the regional domain using a `LatitudeLongitudeGrid`. For real-world bathymetry,
`build_immersed_grid_from_real_data` automatically normalizes coordinates and regrids
the elevation matrix onto the target model resolution via 2D bilinear interpolation:

```julia
# 1. Define base spherical grid
underlying_grid = build_shelf_grid(
  lon_range=lon_range,
  lat_range=lat_range,
  z_range=z_range,
  grid_size=grid_size)
  
# 2. Build immersed boundary (using real bathymetry with automatic 2D regridding)
grid = build_immersed_grid_from_real_data(
  underlying_grid, 
  real_bathy_file)

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
tides = build_tidal_body_forcing(
  constituents=[:M2],
  u_amplitudes=Dict(:M2 => 0.25),
  v_amplitudes=Dict(:M2 => 0.12))
  
# 2. Instantiate hydrostatic model with wind stress and tidal forcing
model = build_hydrodynamic_model(
  grid,
  coriolis_latitude=45.0,
  surface_wind_stress_x=tau_x,
  surface_wind_stress_y=tau_y,
  tidal_forcing=tides,
  ν=1e-2, κ=1e-2,
  tracers=(:T, :S))
  
# 3. Set initial thermal and haline stratification (baseline)
set_initial_stratification!(
  model,
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
simulation = setup_hydrodynamic_simulation(
  model,
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
# 1. Initialize larval particle cohort from benthic egg release near seafloor
bathy_field_fn(lon, lat)      = -180.0         # Local seabed elevation (meters)
n_larvae = 200
larvae = initialize_larval_particles(
  n_larvae,
  lon_range=(-64.0, -62.0),
  lat_range=(43.5, 45.0),
  bathymetry_fn=bathy_field_fn,
  release_depth_mode=:bottom,
  bottom_offset=(0.5, 3.0),
  stage=:zoea1)

# 2. Multi-day Lagrangian trajectory simulation with ascent, DVM, and molting
flow_field_fn(lon, lat, z, t) = (0.05 + 0.02 * sin(t / 43200.0), 0.02, 0.0001)
temp_field_fn(lon, lat, z, t) = 4.0 + 0.01 * z # In-situ temperature

trajectories = track_larval_cohort(
  larvae,
  velocity_fn=flow_field_fn,
  temperature_fn=temp_field_fn,
  bathymetry_fn=bathy_field_fn,
  total_duration=86400.0 * 10, # 10 days
  dt=300.0,                    # 5 min step
  κ_h=10.0, κ_v=1e-4,
  is_lat_lon=true,
  enable_tides=true,
  enable_molting=true,
  enable_initial_ascent=true,
  ascent_speed=0.010,
  ascent_target_depth=-10.0)
  
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
emp_mov = estimate_empirical_movement(
  trajectories,
  lon_bins=range(-68.0, -57.0, length=30),
  lat_bins=range(42.0, 47.0, length=30))

nanmean(x) = mean(filter(!isnan, x))
println("Empirical Movement Analysis:")
println("  Mean empirical speed: ", round(nanmean(emp_mov.speed_mean), digits=4), " m/s")
println("  Mean empirical diffusivity: ", round(nanmean(emp_mov.diffusivity), digits=2), " m^2/s")

# 2. Compute gridded larval recruitment, settlement density, and self-retention
rec_metrics = compute_gridded_recruitment_metrics(
  trajectories,
  lon_bins=range(-68.0, -57.0, length=30),
  lat_bins=range(42.0, 47.0, length=30))

println("Total settled recruits: ", sum(rec_metrics.settlement_density))
println("Successful recruits in nursery: ", sum(rec_metrics.successful_settlement_density))

# 3. Compute gridded thermal exposure and degree-days
therm_metrics = compute_gridded_thermal_metrics(
  trajectories,
  lon_bins=range(-68.0, -57.0, length=30),
  lat_bins=range(42.0, 47.0, length=30))

# 4. Derive macro-regional connectivity matrix between Crab Fishing Areas (CFAs)
cfa_definitions = [
    (name="CFA 20-22 (Eastern NS)", lon=(-62.0, -57.0), lat=(44.5, 47.5)),
    (name="CFA 23-24 (Middle Shelf)", lon=(-64.5, -60.0), lat=(43.0, 45.5)),
    (name="CFA 4X (Southwest NS)", lon=(-68.0, -64.0), lat=(42.0, 44.5)),
    (name="Offshore / Slope", lon=(-68.0, -57.0), lat=(40.0, 43.0))
]

conn = compute_empirical_connectivity(
  trajectories, 
  strata_definitions=cfa_definitions)

println("\nTransition Probability Matrix (P_ij):")
display(round.(conn.matrix, digits=3))

# 5. Export comprehensive multi-variable NetCDF and JLD2 archives
export_larval_dispersal_netcdf(
  "outputs/larval_dispersal_analysis.nc",
  trajectories=trajectories,
  strata_definitions=cfa_definitions)
  
export_larval_dispersal_jld2(
  "outputs/larval_dispersal_analysis.jld2",
  trajectories=trajectories,
  strata_definitions=cfa_definitions)
```

---

## 9. Visualization & Cross-Scenario Comparisons

Render particle drift trajectories, empirical velocity quivers, annotated connectivity matrices,
and multi-scenario climate comparisons:

```julia
# 1. Plot 2D trajectories over Scotian Shelf bathymetry
bathy_grid_data = load_bathymetry_from_netcdf(real_bathy_file)

fig_tracks = plot_particle_trajectories(
  trajectories,
  bathymetry_data=bathy_grid_data,
  output_path="outputs/larval_trajectories.png")

# 2. Plot Diel Vertical Migration (DVM) depth profiles over time
fig_dvm = plot_vertical_migration_profiles(
  trajectories,
  sample_indices=1:5,
  output_path="outputs/dvm_depth_profiles.png")
  
# 3. Plot 2D settlement nursery density distribution
fig_density = plot_larval_dispersal_density(
  trajectories,
  output_path="outputs/settlement_density.png")

# 4. Plot empirical velocity vector arrows over turbulent diffusivity
fig_emp = plot_empirical_movement_field(
  emp_mov,
  output_path="outputs/empirical_movement_field.png")

# 5. Plot annotated macro-regional connectivity matrix
fig_conn = plot_connectivity_matrix(
  conn,
  output_path="outputs/regional_connectivity_matrix.png")

# 6. Plot larval thermal exposure and cumulative degree-days
fig_therm = plot_thermal_exposure_map(
  therm_metrics,
  output_path="outputs/thermal_exposure_map.png")

# 7. Plot recruitment and settlement summary
fig_rec = plot_recruitment_summary(
  rec_metrics,
  output_path="outputs/recruitment_summary.png")

# 8. Compare dispersal patterns across multiple climate scenarios
flow_ssp585(lon, lat, z, t) = (0.07 + 0.03 * sin(t / 43200.0), 0.03, 0.0001)
traj_ssp585 = track_larval_cohort(
  larvae, velocity_fn=flow_ssp585,
  total_duration=86400.0 * 10, 
  dt=300.0)

scenario_comp = Dict(
  :historical => trajectories,
  :ssp585_2050 => traj_ssp585
)
fig_comp = compare_scenario_dispersal(
  scenario_comp,
  output_path="outputs/climate_scenario_comparison.png")
```

---

## 10. Dashboard
 
```julia 
 
# Export interactive HTML5 + Leaflet oceanographic dashboard
export_interactive_tracks_html(
    "outputs/production_interactive_tracks.html",
    trajectories = trajectories,
    bathymetry = bathy_grid_data,
    strata_definitions = conn.strata_names,
    title = "Scotian Shelf Snow Crab Larval Dispersal & Connectivity Dashboard"
)
 
```

---

## 11. Hardware Acceleration & Interactive Visualization

### GPU / CUDA Hardware Acceleration
Oceananigans natively supports NVIDIA CUDA GPU execution using unified memory and streaming multiprocessors. The workflow integrates [`resolve_architecture`](file:///c:/home/jae/projects/ParticleTracking/src/architecture.jl).

---

### Standalone Interactive HTML5 Map Dashboard
The visualizer exports an interactive Leaflet.js dashboard (`outputs/interactive_larval_tracks.html`) that can be opened in any browser:
- **Base Maps**: ESRI Ocean Basemap, CartoDB Dark Matter, OpenStreetMap.
- **Dynamic Playback**: Simulation time scrubber ($t = 0 \to 5\text{ days}$) with Play/Pause animation.
- **Stage Coloration**: Zoea I–IV, Megalopa, Instar I (nursery recruit), and thermal mortality.
- **Live Telemetry HUD**: Interactive popups with depth ($z$), temperature ($T$), degree-days ($DD$), and survival probability ($S(t)$).
- **Management Strata**: Bounding overlays for CFAs 20–22, 23–24, and 4X.

---

## 12. DuckDB Analytical Storage, Scenario Management & Ensemble Model Averaging

Similar to the BSTM modeling framework, `ParticleTracking.jl` includes a high-performance **DuckDB** analytical storage backend (`src/storage_duckdb.jl`) for persisting multi-scenario simulation runs, millions of trajectory steps, demographic transition matrices, and gridded dispersal fields into a single relational database (`outputs/particle_tracking.duckdb`).

### Database Relational Schema
1. **`simulation_runs`**: Metadata for each run (scenario name, projection year, $N_{\text{particles}}$, duration, time step, physical & biological options, seed, timestamps).
2. **`particle_trajectories`**: Columnar time series of particle coordinates $(\lambda, \phi, z)$, ambient temperature $T$, degree-days $DD$, survival probability $S(t)$, developmental stage, and settlement state.
3. **`recruitment_metrics`**: Summary outcomes per cohort (settlement success rate, mean PLD, mean degree-days, thermal mortality, mean dispersal displacement).
4. **`connectivity_transitions`**: Demographic transition probabilities $P_{ij}$ and raw transit counts between spatial management strata (e.g. CFAs 20–22, 23–24, 4X).
5. **`gridded_dispersal_summary`**: Spatial matrices of empirical velocities $(\bar{u}, \bar{v})$, turbulent diffusivity $D_{\text{emp}}$, settlement density, and thermal exposure.


### API Usage

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

## 13. Centralized Configuration Management (`inputs/ParticleTracking.config`)

To streamline workflow reproducibility and parameter configuration across multiple
scenarios, all user-configurable defaults are declared in a centralized TOML
configuration file at [`inputs/ParticleTracking.config`](file:///c:/home/jae/projects/ParticleTracking/inputs/ParticleTracking.config).

### Configuration Sections
- `[domain]`: Bounding coordinates ($\lambda_{\min}, \lambda_{\max}, \phi_{\min}, \phi_{\max}, z_{\min}, z_{\max}$).
- `[grid]`: Numerical grid dimensions ($N_x, N_y, N_z$).
- `[data]`: Environmental data mode (`synthetic` vs `real`), dataset IDs, and synthetic shelf parameters.
- `[tides]`: Barotropic tidal forcing options ($M_2$ constituent amplitude, period, phase).
- `[climate]`: Climate scenario selection (`ssp126`, `ssp245`, `ssp585`, `mhw`), projection and baseline years.
- `[hydrodynamics]`: Oceananigans integration duration, initial time step, adaptive CFL target.
- `[biology]`: Larval cohort size ($N$), drift duration, tracking step, minimum depth ($100\text{ m}$),
  diffusivities ($\kappa_h, \kappa_v$), release depth mode (`:bottom`, `:range`, `:surface`), bottom offset
  ($[0.5, 3.0]\text{ m}$), and active vertical ascent swimming.
- `[dvm]`: Stage-specific Diel Vertical Migration daytime/nighttime target depths and swimming speeds.
- `[molting_and_settlement]`: Degree-day thresholds ($150, 310, 510\text{ DD}$), thermal mortality
  parameters, and benthic nursery suitability windows ($-250\text{ m} \le z \le -50\text{ m}$, $T \le 6^\circ\text{C}$).
- `[storage]`: DuckDB analytical database persistence and Parquet export.
- `[hardware]`: NVIDIA CUDA GPU hardware acceleration and automatic CPU fallback.
- `[visualization]`: Interactive HTML5 Leaflet map export and dashboard options.
- `[paths]`: File system directories (`inputs`, `outputs`) and pseudorandom seed.

```julia
# Load centralized configuration
cfg = load_configuration("inputs/ParticleTracking.config")

# Convert configuration dict to validated HydrodynamicOptions
opts = configuration_to_options(cfg, scenario = :ssp585, n_particles = 250)

# Save active options back to a configuration file
save_configuration(options_to_configuration(opts), "inputs/custom_run.config")
```

---

## 14. Administrative Area (CFA) Polygons

The regional modeling platform ingests official administrative management boundaries
defined in CSV/DAT format (`lon,lat`):
- [`inputs/cfa4x.dat`](file:///c:/home/jae/projects/ParticleTracking/inputs/cfa4x.dat): Southwest Nova Scotia (CFA 4X, 30 vertices).
- [`inputs/cfanorth.dat`](file:///c:/home/jae/projects/ParticleTracking/inputs/cfanorth.dat): Eastern Cape Breton & Glace Bay (CFA North / 20–22, 17 vertices).
- [`inputs/cfasouth.dat`](file:///c:/home/jae/projects/ParticleTracking/inputs/cfasouth.dat): Middle Scotian Shelf & Halifax (CFA South / 23–24, 28 vertices).

### Point-in-Polygon Classification
Spatial connectivity calculations use the Jordan Curve (ray-casting) theorem via
`point_in_polygon(lon, lat, poly_lons, poly_lats)` to assign particle coordinates
$(\lambda_p, \phi_p)$ directly to irregular administrative polygons rather than rectangular bounding boxes.

### Leaflet Visualization
Interactive Leaflet dashboards automatically load all available `inputs/cfa*.dat` polygons,
rendering vector overlays via `L.polygon(stratum.polygon, ...)` with stage-dependent color
highlights and management zone tooltips.

---

## 15. Operational Workflows: Decoupled Multi-Cohort Tracking & Segmented Execution

### 15.1 Architectural Separation of Hydrodynamics and Lagrangian Tracking

A critical design feature of `ParticleTracking.jl` is the complete mathematical and computational
separation between **Eulerian hydrodynamics** (solving 3D Navier-Stokes momentum, continuity, and
tracer equations on a grid) and **Lagrangian particle tracking** (solving Euler-Maruyama stochastic
differential equations for discrete individuals):

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 1. 3D Hydrodynamic Simulation (Oceananigans.jl on GPU / CPU)               │
│    Solves Navier-Stokes, Coriolis, Stratification & Tidal Equations          │
│    Execution time: Minutes to Hours (Heavy PDE Solvers)                     │
│    Output: 4D Flow & Temperature Fields -> outputs/simulation_flow.jld2     │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │ Saved checkpoint archive
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ 2. Sub-Microsecond 4D Flow Interpolation Wrapper                           │
│    create_flow_interpolator_from_jld2("outputs/simulation_flow.jld2")       │
│    Zero-allocation in-memory cubic/linear spatio-temporal lookup             │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │ Rapid re-use across cohorts
         ┌─────────────────────────────┼─────────────────────────────┐
         ▼                             ▼                             ▼
┌─────────────────┐           ┌─────────────────┐           ┌─────────────────┐
│ Cohort A:       │           │ Cohort B:       │           │ Cohort C:       │
│ Spring Bloom    │           │ Peak Hatch      │           │ Late Summer     │
│ Benthic Ascent  │           │ Higher Speed    │           │ Surface Release │
│ (t0 = 0 d)      │           │ (t0 = 15 d)     │           │ (t0 = 30 d)     │
└────────┬────────┘           └────────┬────────┘           └────────┬────────┘
         │                             │                             │
         ▼                             ▼                             ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ 3. Unified Relational DuckDB Storage (outputs/particle_tracking.duckdb)      │
│    save_simulation_run!(db, "cohort_spring_benthic", opts; ...)             │
│    save_simulation_run!(db, "cohort_peak_ascent15", opts; ...)              │
│    save_simulation_run!(db, "cohort_late_surface", opts; ...)               │
│    Tagged metadata, SQL cohort filtering, scenario comparison & Parquet     │
└─────────────────────────────────────────────────────────────────────────────┘
```

Because Lagrangian tracking is thousands of times faster than solving fluid equations,
an entire season of hydrodynamic circulation can be computed once and archived. Then,
dozens or hundreds of distinct larval cohorts—varying in hatch dates ($t_0$), release depth
modes, vertical ascent velocities, DVM behaviors, or stochastic realization seeds—can be
simulated and tagged without re-running the hydrodynamic model.

---

### 15.2 Workflow 1: Scripted Multi-Cohort Batching via Julia API

This workflow demonstrates how to load an existing hydrodynamic JLD2 checkpoint, instantiate
a fast 4D flow interpolator, loop across multiple larval cohorts with different hatch dates
and ascent parameters, persist each into DuckDB with unique tags, and perform cross-cohort
comparative analytics.

```julia
using ParticleTracking
using Random
using Dates
using Statistics

# Step 1: Open DuckDB analytical database
db_path = "outputs/particle_tracking.duckdb"
db = open_duckdb_storage(db_path)

# Step 2: Ingest 4D flow field from completed hydrodynamic simulation
flow_checkpoint = "outputs/simulation_flow.jld2"
bathy_checkpoint = "inputs/bathymetry_active.nc"

if !isfile(flow_checkpoint)
    error("Hydrodynamic flow checkpoint not found at $(flow_checkpoint). " *
          "Run hydrodynamic simulation once before cohort tracking.")
end

# Build high-performance in-memory 4D interpolator: (lon, lat, z, t) -> (u, v, w, T)
flow_interpolator = create_flow_interpolator_from_jld2(flow_checkpoint)

velocity_fn = (lon, lat, z, t) -> begin
    f = flow_interpolator(lon, lat, z, t)
    (f.u, f.v, f.w)
end

temperature_fn = (lon, lat, z, t) -> begin
    f = flow_interpolator(lon, lat, z, t)
    f.T
end

# Bathymetry elevation interpolator: (lon, lat) -> z_bed
bathy_fn = if isfile(bathy_checkpoint)
    get_bathymetry_interpolator(bathy_checkpoint)
else
    (lon, lat) -> -120.0 - 150.0 * (lat - 42.0) / 5.0
end

# Administrative management strata for demographic connectivity
cfa_definitions = [
    (name = "CFA 20-22 (North)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
    (name = "CFA 23-24 (South)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
    (name = "CFA 4X (SW Nova)",  lon = (-68.0, -64.0), lat = (42.0, 44.5))
]

# Step 3: Define experimental cohort parameter matrix
# Test different hatch dates (temporal offsets t0), release modes, and ascent speeds
cohort_configs = [
    (
        tag = "cohort_2020_spring_benthic",
        notes = "Spring bloom hatch (day 0), benthic release with active ascent (10 mm/s)",
        t0 = 0.0,
        n_particles = 200,
        release_mode = :bottom,
        enable_ascent = true,
        ascent_speed = 0.010,
        seed = 101
    ),
    (
        tag = "cohort_2020_spring_ascent_fast",
        notes = "Spring bloom hatch (day 0), fast benthic ascent (15 mm/s)",
        t0 = 0.0,
        n_particles = 200,
        release_mode = :bottom,
        enable_ascent = true,
        ascent_speed = 0.015,
        seed = 102
    ),
    (
        tag = "cohort_2020_summer_benthic",
        notes = "Mid-summer hatch (day 15 offset), benthic release with active ascent",
        t0 = 15.0 * 86400.0,
        n_particles = 200,
        release_mode = :bottom,
        enable_ascent = true,
        ascent_speed = 0.010,
        seed = 201
    ),
    (
        tag = "cohort_2020_summer_surface_ctrl",
        notes = "Control: Mid-summer hatch (day 15 offset), direct surface release without ascent",
        t0 = 15.0 * 86400.0,
        n_particles = 200,
        release_mode = :surface,
        enable_ascent = false,
        ascent_speed = 0.0,
        seed = 202
    )
]

# Base hydrodynamic options container for metadata logging
base_opts = HydrodynamicOptions(
    scenario = :ssp245,
    projection_year = 2020,
    track_duration = 10.0 * 86400.0,
    track_dt = 300.0
)

# Step 4: Batch execute Lagrangian tracking across all cohorts
for cfg in cohort_configs
    println("\n>>> Processing Cohort: $(cfg.tag)")
    println("    $(cfg.notes)")

    rng = MersenneTwister(cfg.seed)

    # Initialize larvae according to release mode
    larvae = initialize_larvae(
        cfg.n_particles,
        domain_lon = (-68.0, -57.0),
        domain_lat = (42.0, 47.0),
        release_depth = cfg.release_mode,
        bottom_offset = (0.5, 3.0),
        bathymetry_fn = bathy_fn,
        rng = rng
    )

    # Track cohort trajectory through time using shared flow field
    trajectories = track_larval_cohort(
        larvae,
        velocity_fn = velocity_fn,
        temperature_fn = temperature_fn,
        bathymetry_fn = bathy_fn,
        total_duration = base_opts.track_duration,
        dt = base_opts.track_dt,
        t0 = cfg.t0,
        κ_h = 10.0,
        κ_v = 1e-4,
        is_lat_lon = true,
        enable_tides = true,
        enable_molting = true,
        enable_bbl = true,
        enable_sinking = true,
        enable_initial_ascent = cfg.enable_ascent,
        ascent_speed = cfg.ascent_speed,
        ascent_target_depth = -10.0,
        rng = rng
    )

    # Compute demographic connectivity matrix
    conn = compute_empirical_connectivity(
        trajectories,
        strata_definitions = cfa_definitions
    )

    # Compute empirical movement and retention fields
    emp_mov = estimate_empirical_movement(
        trajectories,
        lon_bins = range(-68.0, -57.0, length = 25),
        lat_bins = range(42.0, 47.0, length = 25)
    )

    rec_metrics = compute_gridded_recruitment_metrics(
        trajectories,
        lon_bins = range(-68.0, -57.0, length = 25),
        lat_bins = range(42.0, 47.0, length = 25)
    )

    therm_metrics = compute_gridded_thermal_metrics(
        trajectories,
        lon_bins = range(-68.0, -57.0, length = 25),
        lat_bins = range(42.0, 47.0, length = 25)
    )

    # Tag and persist run into DuckDB with custom run_id and metadata
    cohort_opts = HydrodynamicOptions(
        scenario = base_opts.scenario,
        projection_year = base_opts.projection_year,
        n_particles = cfg.n_particles,
        track_duration = base_opts.track_duration,
        track_dt = base_opts.track_dt,
        release_depth_mode = cfg.release_mode,
        enable_initial_ascent = cfg.enable_ascent,
        ascent_speed = cfg.ascent_speed,
        seed = cfg.seed
    )

    save_simulation_run!(
        db,
        cfg.tag,
        cohort_opts;
        trajectories = trajectories,
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
        notes = cfg.notes
    )

    println("    Archived $(cfg.tag) in DuckDB with $(cfg.n_particles) particles.")
end

# Step 5: Post-hoc comparative analytics across all tagged cohorts
println("\n=================================================================")
println(" Multi-Cohort Comparative Analytics")
println("=================================================================")

df_runs = list_simulation_runs(db)
println("Archived runs in database:")
display(df_runs[:, [:run_id, :scenario, :n_particles, :notes]])

# Compare recruitment success and dispersal distances across cohorts
df_comparison = compare_scenarios(db)
display(df_comparison)

# Compute ensemble model average across benthic release cohorts
benthic_run_ids = [
    "cohort_2020_spring_benthic",
    "cohort_2020_spring_ascent_fast",
    "cohort_2020_summer_benthic"
]

ens = compute_ensemble_model_average(
    db,
    benthic_run_ids,
    weights = [0.4, 0.2, 0.4]
)

println("\nEnsemble-Averaged Connectivity Matrix (P_ij):")
display(round.(ens.mean_connectivity, digits = 3))
println("\nEnsemble Uncertainty (Standard Deviation):")
display(round.(ens.std_connectivity, digits = 3))

# Export database tables to Apache Parquet for downstream analysis
export_duckdb_to_parquet(db, "outputs/parquet_cohorts")

close_duckdb_storage(db)
println("\nMulti-cohort batching and analysis complete.")
```

---

### 15.3 Workflow 2: Production Segmented CLI Pipeline (`ParticleTrackingRun.jl`)

In production and High-Performance Computing (HPC) environments, the CLI runner
[`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl) allows decoupling hydrodynamics from tracking directly from the command
line using segment flags:

#### Step 1: Execute 3D Hydrodynamics Once
Run environmental data ingestion, grid generation, model configuration, and hydrodynamic
time integration to generate the flow checkpoint `outputs/baseline/simulation_flow.jld2`:

```bash
# Run Hydrodynamics (Segments 1 through 5)
julia --project=. ParticleTrackingRun.jl \
    --data --grid --model --sim \
    --duration=120.0 \
    --grid=50,50,10 \
    --scenario=ssp245 \
    --year=2050 \
    --output-dir=outputs/baseline
```

#### Step 2: Run Lagrangian Tracking for Cohort A (Spring Benthic Release)
Track larvae using the existing flow checkpoint without re-running hydrodynamics,
archiving results to DuckDB under a dedicated tag:

```bash
# Run Lagrangian tracking & analytics for Cohort A
julia --project=. ParticleTrackingRun.jl \
    --track --metrics --viz \
    --output-dir=outputs/baseline \
    --particles=500 \
    --track-duration=10.0 \
    --release-mode=bottom \
    --ascent \
    --ascent-speed=0.010 \
    --ascent-target=-10.0 \
    --seed=101 \
    --duckdb \
    --db-path=outputs/particle_tracking.duckdb
```

#### Step 3: Run Lagrangian Tracking for Cohort B (Summer Surface Control)
Track an alternate cohort using the same hydrodynamic solution, varying the biological
behavior:

```bash
# Run Lagrangian tracking & analytics for Cohort B
julia --project=. ParticleTrackingRun.jl \
    --track --metrics --viz \
    --output-dir=outputs/baseline \
    --particles=500 \
    --track-duration=10.0 \
    --release-mode=surface \
    --no-ascent \
    --seed=202 \
    --duckdb \
    --db-path=outputs/particle_tracking.duckdb
```

#### Step 4: Query Database and Compare Cohorts
Inspect the archived cohorts and compare metrics without launching any simulations:

```bash
# List all completed runs in DuckDB
julia --project=. ParticleTrackingRun.jl --list-runs --db-path=outputs/particle_tracking.duckdb

# Output side-by-side recruitment and connectivity metrics
julia --project=. ParticleTrackingRun.jl --compare-scenarios --db-path=outputs/particle_tracking.duckdb

# Compute weighted ensemble model average
julia --project=. ParticleTrackingRun.jl --model-average --db-path=outputs/particle_tracking.duckdb

# Export all tables to Apache Parquet for external Python / R / BSTM analysis
julia --project=. ParticleTrackingRun.jl --export-parquet --db-path=outputs/particle_tracking.duckdb
```

---

### 15.4 Workflow 3: Dedicated Snow Crab Dispersal Platform (`--snowcrab-settings`)

For dedicated snow crab (*Chionoecetes opilio*) assessment workflows, the platform is integrated
directly into [`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl)
via the `--snowcrab-settings` flag (shorthand: `--snowcrab`). Passing `--snowcrab-settings` automatically
pre-configures all snow crab calibrated parameters (500 larvae, 60-day PLD, bottom boundary release
with active vertical ascent, -3500 m to 0 m depth domain, -1.5°C base molting temperature, and DuckDB
target database `outputs/snowcrab_tracking.duckdb`), while any additional CLI arguments override
those defaults.

For programmatic Julia scripting, the exported function `SnowCrabRunOptions(; kwargs...)` returns
a calibrated `HydrodynamicOptions` instance with identical parameters.

#### Step 1: Hydrodynamic Simulation Only (`--hydro-only`)
Integrate Oceananigans hydrostatic equations and persist flow fields into a chosen
`--hydro-model` file:

```bash
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --real-5yr --hydro-only --hydro-model=hydrodynamics1.jld2
```

#### Step 2: Multi-Cohort Lagrangian Tracking Reusing Hydrodynamics (`--track-only`)
Instantly load the pre-computed flow fields from `--hydro-model` and track cohorts with
different biological parameters (which override the snow crab defaults), tagging each in DuckDB:

```bash
# Track Cohort 1: Spring benthic release with active ascent (~10 mm/s)
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_baseline --particles=500 --ascent

# Track Cohort 2: Fast ascent swimming (~15 mm/s, overriding default 10 mm/s)
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_fast --particles=500 --ascent-speed=0.015

# Track Cohort 3: Surface release control (no ascent)
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_surface --particles=500 --release-mode=surface --no-ascent
```

#### Step 3: Comparative Analytics Across Cohorts
Query DuckDB and generate side-by-side demographic and recruitment summaries:

```bash
julia --project=. ParticleTrackingRun.jl --compare
```

---

### 15.5 Standalone Interactive HTML5 + Leaflet.js Dashboard

Exported automatically via `export_interactive_tracks_html` to
[`outputs/interactive_larval_tracks.html`](file:///c:/home/jae/projects/ParticleTracking/outputs/interactive_larval_tracks.html):
- **Multiple Basemaps**: Seamless switching between ESRI Ocean Basemap, CartoDB Dark Matter, and OpenStreetMap.
- **Dynamic Playback Control**: Time slider scrubber ($t = 0 \to 5\text{ days}$) with Play/Pause animation and speed toggles ($1\times, 5\times, 20\times$).
- **Ontogenetic Color-Coding**: Dynamic color transitions for Zoea I, Zoea II, Megalopa, Instar I (settled recruit), and thermal mortality.
- **Live Telemetry HUD**: Hover and click popups displaying real-time depth ($z$), in-situ temperature ($T$), cumulative degree-days ($DD$), and survival probability ($S(t)$).
- **Administrative CFA Overlays**: True boundary vector polygons (`L.polygon`) for CFA North, CFA South, and CFA 4X with interactive tooltips.

---

## 16. Developing & Extending the Framework

### 16.1 Running the Full Test Suite
The continuous integration test suite validates all 17 component test sets with **618 unit tests**:

```bash
julia --project=. test/runtests.jl
```

The test sets cover:
1. Drag law and open data ingestion utilities.
2. Synthetic data generation and NetCDF dataset inspection.
3. Immersed boundary grid geometry and Coriolis parameters.
4. Hydrodynamic model configuration, boundary conditions, and CMIP6 scenarios.
5. Tidal dynamics, harmonic body forcing, and Simpson-Hunter stratification.
6. Larval DVM, tidal velocity superposition, and settlement environmental filters.
7. Lagrangian particle tracking, numerical integration, and marine land exclusion.
8. Empirical movement velocity and turbulent diffusivity estimation.
9. Gridded recruitment metrics and thermal exposure analysis.
10. Empirical demographic connectivity matrix computation and NetCDF / JLD2 export.
11. CairoMakie spatial figures, vertical profiles, and diagnostics rendering.
12. GPU architecture detection and automatic CPU fallback.
13. Standalone interactive HTML5 Leaflet dashboard generation.
14. DuckDB analytical database storage, multi-scenario querying, and ensemble averaging.
15. Centralized TOML configuration parsing, validation, and scenario metadata.
16. Coastline geometry, 0% land seeding, and CFA polygon intersection.
17. Enhanced physical and biophysical processes (bottom release, active ascent, surface heat flux).

### 16.2 Adding a New Climate Scenario
To add a new CMIP6 scenario (e.g. `ssp370`):
1. In [`src/climate_scenarios.jl`](file:///c:/home/jae/projects/ParticleTracking/src/climate_scenarios.jl), register a new `ClimateScenarioDelta` entry in `get_climate_scenario`:
   ```julia
   elseif sc_sym == :ssp370
       return ClimateScenarioDelta(:ssp370, "CMIP6 SSP3-7.0", 2.6, 1.4, -0.65, 1.15)
   ```
2. Update the CLI choices in [`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl) and configuration parsing in [`src/configuration.jl`](file:///c:/home/jae/projects/ParticleTracking/src/configuration.jl).

### 16.3 Adding a New Biological Species
To model a different marine species (e.g. American Lobster *Homarus americanus*):
1. Adjust degree-day molting thresholds and swimming target depths in [`inputs/ParticleTracking.config`](file:///c:/home/jae/projects/ParticleTracking/inputs/ParticleTracking.config).
2. Specify nursery depth and bottom temperature limits under `[molting_and_settlement]`.
3. Configure stage-specific vertical migration targets and phototaxis responses under `[dvm]`.

---

## 17. References

### Snow Crab Larval Ecology & Lagrangian Transport
- **Boudreau, S. A., Anderson, S. C., & Worm, B.** (2011). Top-down interactions and temperature constraints on large-scale patterns of biomass and distribution in snow crab (*Chionoecetes opilio*). *Marine Ecology Progress Series*, 429, 169–183. DOI: [10.3354/meps09082](https://doi.org/10.3354/meps09082)
- **Choi, J. S., & Zisserson, B.** (2012). Assessment of Scotian Shelf snow crab (*Chionoecetes opilio*) in 2011. *DFO Canadian Science Advisory Secretariat Research Document*, 2012/025, 88 pp.
- **Epifanio, C. E., & Cohen, J. H.** (2016). Behavioral adaptations in larvae of brachyuran crabs: a review. *Journal of Experimental Marine Biology and Ecology*, 482, 85–105. DOI: [10.1016/j.jembe.2016.05.006](https://doi.org/10.1016/j.jembe.2016.05.006)
- **Incze, L. S., Armstrong, D. A., & Smith, S. L.** (1987). Abundance of filter-feeding and pelagic stages of crab larvae in the southeastern Bering Sea. *Marine Biology*, 95(2), 195–200. DOI: [10.1007/BF00409006](https://doi.org/10.1007/BF00409006)
- **Kloeden, P. E., & Platen, E.** (1992). *Numerical Solution of Stochastic Differential Equations*. Springer-Verlag, Berlin. DOI: [10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)
- **Kuhn, P. S., & Choi, J. S.** (2011). Influence of temperature on embryo incubation and larval development in snow crab (*Chionoecetes opilio*). *Fisheries Research*, 107(1-3), 81–87. DOI: [10.1016/j.fishres.2010.10.011](https://doi.org/10.1016/j.fishres.2010.10.011)
- **Lovrich, G. A., Sainte-Marie, B., & Smith, B. D.** (1995). Depth distribution and seasonal movements of *Chionoecetes opilio* (Brachyura: Majidae) in Baie Sainte-Marguerite, Gulf of Saint Lawrence. *Canadian Journal of Fisheries and Aquatic Sciences*, 52(4), 903–913. DOI: [10.1139/f95-090](https://doi.org/10.1139/f95-090)
- **North, E. W., Gallego, A., & Petitgas, P. (Eds.)** (2009). Manual of recommended practices for modelling physical - biological interactions during fish early life. *ICES Cooperative Research Report*, No. 295, 111 pp.
- **Okubo, A.** (1971). Oceanic diffusion diagrams. *Deep Sea Research and Oceanographic Abstracts*, 18(8), 789–802. DOI: [10.1016/0011-7471(71)90046-5](https://doi.org/10.1016/0011-7471(71)90046-5)
- **Sainte-Marie, B., Raymond, S., & Brêthes, J. C.** (1996). Determinants of size at morphometric maturity and fecundity in female snow crab, *Chionoecetes opilio*, in the Gulf of St. Lawrence. *Canadian Journal of Fisheries and Aquatic Sciences*, 53(11), 2419–2426. DOI: [10.1139/f96-200](https://doi.org/10.1139/f96-200)
- **Sainte-Marie, G., & Sainte-Marie, B.** (1999). Growth, developmental stages, and vertical distribution of snow crab larvae (*Chionoecetes opilio*) in the northwestern Gulf of St. Lawrence. *Canadian Journal of Fisheries and Aquatic Sciences*, 56(11), 2181–2193. DOI: [10.1139/f99-151](https://doi.org/10.1139/f99-151)
- **Tremblay, M. J.** (1997). Snow crab (*Chionoecetes opilio*) distribution and abundance in the Eastern Nova Scotia area. *DFO Canadian Science Advisory Secretariat Research Document*, 97/80, 24 pp.

### Ocean Hydrodynamics & Numerical Methods
- **Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A.** (2007). *Spectral Methods: Evolution to Complex Geometries and Applications to Fluid Dynamics*. Springer-Verlag, Berlin.
- **Courant, R., Friedrichs, K., & Lewy, H.** (1928). Über die partiellen Differenzengleichungen der mathematischen Physik. *Mathematische Annalen*, 100(1), 32–74. DOI: [10.1007/BF01448839](https://doi.org/10.1007/BF01448839)
- **Hannah, C. G., Shore, J. A., Loder, J. W., & Xu, Z.** (2001). Seasonal circulation on the Western and Central Scotian Shelf. *Journal of Physical Oceanography*, 31(2), 591–615. DOI: [10.1175/1520-0485(2001)031<0591:SCOTWA>2.0.CO;2](https://doi.org/10.1175/1520-0485(2001)031<0591:SCOTWA>2.0.CO;2)
- **Loder, J. W., Han, G., Hannah, C. G., Greenberg, D. A., & Smith, P. C.** (1997). Hydrography and circulation in the Scotian Shelf–Gulf of Maine region. *Canadian Journal of Fisheries and Aquatic Sciences*, 54(S1), 95–113. DOI: [10.1139/f97-158](https://doi.org/10.1139/f97-158)
- **Marshall, J., Adcroft, A., Hill, C., Perelman, L., & Heisey, C.** (1997). A finite-volume, incompressible Navier Stokes model for studies of the ocean on parallel computers. *Journal of Geophysical Research: Oceans*, 102(C3), 5753–5766. DOI: [10.1029/96JC02775](https://doi.org/10.1029/96JC02775)
- **Ramadhan, A., Marshall, J., Hill, C., Campin, J. M., Bischoff, T., & Wagner, G. L.** (2020). Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs. *Journal of Open Source Software*, 5(53), 2018. DOI: [10.21105/joss.02018](https://doi.org/10.21105/joss.02018)
- **Vallis, G. K.** (2017). *Atmospheric and Oceanic Fluid Dynamics: Fundamentals and Large-Scale Circulation*. 2nd Edition. Cambridge University Press, Cambridge. DOI: [10.1017/9781107588417](https://doi.org/10.1017/9781107588417)
- **Verzicco, R.** (2023). Immersed boundary methods for ocean modeling. *Annual Review of Fluid Mechanics*, 55, 305–333. DOI: [10.1146/annurev-fluid-030322-040713](https://doi.org/10.1146/annurev-fluid-030322-040713)

### Climate Forcing & Regional Projections
- **Brickman, D., Wang, Z., & DeTracey, B.** (2018). Variability and trends in the Scotian Shelf and Gulf of Maine region from a high-resolution regional ocean climate model. *Progress in Oceanography*, 164, 49–64. DOI: [10.1016/j.pocean.2018.04.004](https://doi.org/10.1016/j.pocean.2018.04.004)
- **Hobday, A. J., et al.** (2016). A hierarchical approach to defining marine heatwaves. *Progress in Oceanography*, 141, 227–238. DOI: [10.1016/j.pocean.2015.12.014](https://doi.org/10.1016/j.pocean.2015.12.014)
- **Loder, J. W., van der Baaren, A., & Yashayaev, I.** (2015). Climate change trends and projections for the Canadian Northwest Atlantic. *Canadian Technical Report of Hydrography and Ocean Sciences*, 305, 142 pp.
- **O'Neill, B. C., et al.** (2016). The Scenario Model Intercomparison Project (ScenarioMIP) for CMIP6. *Geoscientific Model Development*, 9(9), 3461–3482. DOI: [10.5194/gmd-9-3461-2016](https://doi.org/10.5194/gmd-9-3461-2016)
- **Saba, V. S., et al.** (2016). Enhanced warming of the Northwest Atlantic Ocean under climate change. *Journal of Geophysical Research: Oceans*, 121(1), 118–132. DOI: [10.1002/2015JC011346](https://doi.org/10.1002/2015JC011346)

### Atmospheric Forcing & Open Datasets
- **GEBCO Compilation Group** (2023). *GEBCO 2023 Grid*. DOI: [10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e](https://doi.org/10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e)
- **Large, W. G., & Pond, S.** (1981). Open ocean momentum flux measurements in moderate to strong winds. *Journal of Physical Oceanography*, 11(3), 324–336. DOI: [10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2](https://doi.org/10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2)
- **NOAA National Centers for Environmental Information** (2022). *NOAA ETOPO 2022 15 Arc-Second Global Relief Model*. NOAA NCEI. DOI: [10.25921/fd1h-fy81](https://doi.org/10.25921/fd1h-fy81)
- **Simons, R. A.** (2019). *ERDDAP: The Environmental Research Division's Data Access Program*. NOAA CoastWatch / SWFSC.
- **Wu, J.** (1982). Wind-stress coefficients over sea surface from breeze to hurricane. *Journal of Geophysical Research: Oceans*, 87(C12), 9704–9706. DOI: [10.1029/JC087iC12p09704](https://doi.org/10.1029/JC087iC12p09704)
- **Zhang, H.-M., Bates, J. J., & Reynolds, R. W.** (2006). Assessment of composite global sampling: Sea surface wind speed. *Geophysical Research Letters*, 33(17), L17714. DOI: [10.1029/2006GL027086](https://doi.org/10.1029/2006GL027086)

### Tidal Dynamics & Shelf Mixing
- **Egbert, G. D., & Erofeeva, S. Y.** (2002). Efficient inverse modeling of barotropic ocean tides. *Journal of Atmospheric and Oceanic Technology*, 19(2), 183–204. DOI: [10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2](https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2)
- **Pugh, D., & Woodworth, P.** (2014). *Sea-Level Science: Understanding Tides, Surges, Tsunamis and Mean Sea-Level Changes*. Cambridge University Press. DOI: [10.1017/CBO9781139235778](https://doi.org/10.1017/CBO9781139235778)
- **Simpson, J. H., & Hunter, J. R.** (1974). Fronts in the Irish Sea. *Nature*, 250(5465), 404–406. DOI: [10.1038/250404a0](https://doi.org/10.1038/250404a0)

