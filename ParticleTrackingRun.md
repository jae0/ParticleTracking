# Hydrodynamic Model for Snow Crab Larval Particle Tracking

## Overview

This workflow establishes a regional hydrodynamic model of the Scotian Shelf and
adjacent slope waters using [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).
The hydrodynamic solution provides time-evolving advection ($\boldsymbol{u} = (u, v, w)$)
and turbulent diffusion ($\kappa_h, \kappa_v$) fields. These physical fields drive
Lagrangian particle tracking of planktonic (Zoea I, Zoea II) and semiplanktonic (Megalopa)
snow crab (*Chionoecetes opilio*) larvae undergoing stage-specific vertical migration,
temperature-dependent development, and spatial dispersal across coastal nursery habitats.



## Command-Line Interface (CLI) Reference (`ParticleTrackingRun.jl`)

The command-line runner [`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl) provides full operational control over segmented execution, parameter overrides, DuckDB queries, and scenario comparisons.

### Syntax
```bash
julia --project=. ParticleTrackingRun.jl [OPTIONS...]
```

### Complete CLI Flags Reference

| Category              | Flag / Option                            | Argument Type   | Default Value                      | Description                                                                               |
| :----------------------| :-----------------------------------------| :----------------| :-----------------------------------| :------------------------------------------------------------------------------------------|
| **Execution Modes**   | `--all`                                  | Flag            | —                                  | Execute complete 8-segment production pipeline.                                           |
|                       | `--quick`, `-q`                          | Flag            | —                                  | Fast debug mode ($15\times 15\times 5$, $1\text{ h}$ hydro, $2\text{ d}$ track).          |
|                       | `--segment=<name>`                       | String          | `all`                              | Run segment: `data`, `grid`, `model`, `climate`, `sim`, `track`, `metrics`, `viz`, `all`. |
| **Decoupled Hydro**   | `--hydro-model=<path>`                   | String          | `hydrodynamics_<scen>_<yr>.jld2`   | Target hydrodynamic JLD2 model file path (input or output).                               |
|                       | `--hydro-only`                           | Flag            | —                                  | Execute hydrodynamic simulation only (Segments 1–5), saving solution to `--hydro-model`.  |
|                       | `--track-only`                           | Flag            | —                                  | Execute larval tracking only (Segments 6–8), reading flow fields from `--hydro-model`.    |
|                       | `--reuse-hydro`                          | Flag            | —                                  | Reuse existing `--hydro-model` if present; otherwise run hydrodynamic integration.        |
|                       | `--run-id=<string>`                      | String          | `run_<scenario>_<year>`            | Unique cohort run identifier for DuckDB persistence and figure naming.                    |
| **Segment Execution** | `--data`                                 | Flag            | —                                  | Run environmental data ingestion (Segment 1).                                             |
|                       | `--grid`                                 | Flag            | —                                  | Run grid & immersed boundary construction (Segment 2).                                    |
|                       | `--model`                                | Flag            | —                                  | Run hydrodynamic model setup & tidal forcing (Segment 3).                                 |
|                       | `--climate`                              | Flag            | —                                  | Run CMIP6 climate scenario & thermal biology (Segment 4).                                 |
|                       | `--sim`, `--simulation`                  | Flag            | —                                  | Run Oceananigans hydrodynamic time stepping (Segment 5).                                  |
|                       | `--track`, `--tracking`                  | Flag            | —                                  | Run Lagrangian particle tracking & DVM (Segment 6).                                       |
|                       | `--metrics`                              | Flag            | —                                  | Compute retention, diffusion & connectivity (Segment 7).                                  |
|                       | `--viz`, `--visualize`                   | Flag            | —                                  | Generate Makie figures & Leaflet map (Segment 8).                                         |
| **DuckDB Analytics**  | `--duckdb`                               | Flag            | `true`                             | Enable DuckDB scenario storage archiving.                                                 |
|                       | `--no-duckdb`                            | Flag            | —                                  | Disable DuckDB storage archiving.                                                         |
|                       | `--db-path=<path>`                       | String          | `outputs/particle_tracking.duckdb` | Custom file path for the DuckDB analytical database.                                      |
|                       | `--list-runs`                            | Flag            | —                                  | Query and print all archived simulation runs.                                             |
|                       | `--compare-scenarios`                    | Flag            | —                                  | Query and print multi-scenario comparative analytics.                                     |
|                       | `--model-average`                        | Flag            | —                                  | Compute ensemble model-averaged connectivity ($P_{ij} \pm \sigma$).                       |
|                       | `--export-parquet`                       | Flag            | —                                  | Export all DuckDB tables to Apache Parquet format.                                        |
| **Configuration**     | `--config=<path>`                        | String          | `inputs/ParticleTracking.config`   | Load parameter settings from custom `.config` file.                                       |
|                       | `--save-config[=<path>]`                 | String          | `inputs/ParticleTracking.config`   | Export active CLI options to `.config` file and exit.                                     |
| **Species Calibration**| `--snowcrab-settings`                   | Flag            | —                                  | Load calibrated Snow Crab defaults (500 larvae, 60d PLD, ascent, 100x100x20, -3500m to 0m)|
|                       | `--snowcrab`, `--snowcrab-mode`          | Flag            | —                                  | Aliases for `--snowcrab-settings`.                                                        |
|                       | `--real-5yr`                             | Flag            | —                                  | Execute 5-Year physical hydrodynamic cycle scenario (`historical`, 2020).                 |
|                       | `--climatology-2yr`                      | Flag            | —                                  | Execute 2-Year climatological average cycle scenario (`ssp245`, 2022).                    |
|                       | `--compare`                              | Flag            | —                                  | Query DuckDB database and print comparative analytics across runs.                        |
| **Hardware**          | `--gpu`, `--cuda`                        | Flag            | —                                  | Enable NVIDIA CUDA GPU acceleration for hydrodynamics.                                    |
|                       | `--cpu`                                  | Flag            | `true`                             | Execute on multi-threaded CPU.                                                            |
|                       | `--fallback-cpu`                         | Flag            | —                                  | Automatically fall back to CPU if CUDA GPU is absent.                                     |
| **Visualization**     | `--interactive`                          | Flag            | `true`                             | Export standalone interactive HTML5 Leaflet dashboard.                                    |
|                       | `--no-interactive`                       | Flag            | —                                  | Disable interactive HTML map generation.                                                  |
| **Spatial Domain**    | `--lon=<min,max>`                        | Real,Real       | `-68.0,-57.0`                      | Longitude bounding range in degrees East.                                                 |
|                       | `--lat=<min,max>`                        | Real,Real       | `42.0,47.0`                        | Latitude bounding range in degrees North.                                                 |
|                       | `--depth-range=<min,max>`                | Real,Real       | `-1000.0,0.0`                      | Vertical depth range in meters ($z_{\text{bottom}}, z_{\text{surface}}$).                 |
|                       | `--grid=<nx,ny,nz>`                      | Int,Int,Int     | `50,50,10`                         | Grid cell resolution ($N_\lambda, N_\phi, N_z$).                                          |
|                       | `--nx=<int>`, `--ny=<int>`, `--nz=<int>` | Int             | `50`, `50`, `10`                   | Individual grid axis dimensions.                                                          |
| **Data & Tides**      | `--real`                                 | Flag            | —                                  | Fetch real-world NOAA ERDDAP bathymetry and winds.                                        |
|                       | `--synthetic`                            | Flag            | `true`                             | Generate idealized synthetic shelf data.                                                  |
|                       | `--tides`                                | Flag            | `true`                             | Enable astronomical tidal body forcing ($M_2 + S_2$).                                      |
|                       | `--no-tides`                             | Flag            | —                                  | Disable tidal body forcing.                                                               |
|                       | `--tidal-u=<val>`                        | Float (m/s)     | `0.25`                             | Semi-major barotropic tidal current amplitude.                                            |
|                       | `--tidal-v=<val>`                        | Float (m/s)     | `0.12`                             | Semi-minor barotropic tidal current amplitude.                                            |
|                       | `--heat-flux=<val>`                      | Float ($W/m^2$) | `50.0`                             | Net atmospheric surface heat flux (positive warming).                                     |
| **Climate Scenarios** | `--scenario=<name>`                      | Symbol          | `historical`                       | Climate scenario: `historical`, `ssp126`, `ssp245`, `ssp585`, `mhw`.                      |
|                       | `--year=<int>`                           | Int             | `2020`                             | Climate projection horizon year.                                                          |
| **Simulation**        | `--duration=<hours>`                     | Float (hrs)     | `120.0`                            | Hydrodynamic simulation duration in hours.                                                |
|                       | `--sim-duration=<sec>`                   | Float (sec)     | `432000.0`                         | Hydrodynamic simulation duration in seconds.                                              |
|                       | `--sim-dt=<sec>`                         | Float (sec)     | `120.0`                            | Initial hydrodynamic time integration step.                                               |
|                       | `--adaptive-cfl`                         | Flag            | `true`                             | Enable advective CFL-limited adaptive time stepping.                                      |
|                       | `--no-adaptive-cfl`                      | Flag            | —                                  | Disable adaptive CFL time stepping.                                                       |
|                       | `--target-cfl=<val>`                     | Float           | `0.2`                              | Target advective Courant-Friedrichs-Lewy limit.                                           |
| **Lagrangian Drift**  | `--particles=<int>`                      | Int             | `100`                              | Number of larvae to initialize and track.                                                 |
|                       | `--track-duration=<days>`                | Float (days)    | `5.0`                              | Larval cohort tracking duration in days.                                                  |
|                       | `--track-dt=<sec>`                       | Float (sec)     | `300.0`                            | Lagrangian advection time step in seconds.                                                |
|                       | `--min-depth=<meters>`                   | Float (m)       | `100.0`                            | Minimum water depth for larval placement.                                                 |
|                       | `--release-mode=<mode>`                  | Symbol          | `bottom`                           | Larval release depth mode (`bottom`, `range`, `surface`).                                 |
|                       | `--ascent`                               | Flag            | `true`                             | Enable post-hatch vertical ascent toward surface.                                         |
|                       | `--no-ascent`                            | Flag            | —                                  | Disable initial vertical ascent.                                                          |
|                       | `--ascent-speed=<val>`                   | Float (m/s)     | `0.010`                            | Vertical ascent swimming speed (~10 mm/s).                                                |
|                       | `--ascent-target=<val>`                  | Float (m)       | `-10.0`                            | Target depth for ascent completion.                                                       |
|                       | `--dvm`                                  | Flag            | `true`                             | Enable stage-dependent Diel Vertical Migration.                                           |
|                       | `--no-dvm`                               | Flag            | —                                  | Disable Diel Vertical Migration.                                                          |
|                       | `--molting`                              | Flag            | `true`                             | Enable degree-day thermal molting & mortality.                                            |
|                       | `--no-molting`                           | Flag            | —                                  | Disable thermal molting calculations.                                                     |
|                       | `--diff-h=<val>`                         | Float ($m^2/s$) | `10.0`                             | Horizontal turbulent diffusivity ($\kappa_h$).                                            |
|                       | `--diff-v=<val>`                         | Float ($m^2/s$) | `1e-4`                             | Vertical turbulent diffusivity ($\kappa_v$).                                              |
| **I/O & Environment** | `--output-dir=<path>`                    | String          | `outputs`                          | Target directory for outputs and figures.                                                 |
|                       | `--input-dir=<path>`                     | String          | `inputs`                           | Cache directory for input NetCDF datasets.                                                |
|                       | `--seed=<int>`                           | Int             | `42`                               | Pseudorandom number generator seed.                                                       |
|                       | `--help`, `-h`                           | Flag            | —                                  | Print command-line help and flag reference.                                               |

---

### Common Operational CLI Workflows

```bash
# 1. Fast end-to-end debug verification
julia --project=. ParticleTrackingRun.jl --all --quick

# 2. Production run under CMIP6 SSP5-8.5 warming (2050) with 500 larvae
julia --project=. ParticleTrackingRun.jl --all --scenario=ssp585 --year=2050 --particles=500

# 3. Production run with custom configuration file
julia --project=. ParticleTrackingRun.jl --all --config=inputs/ParticleTracking.config

# 4. GPU-accelerated run with automatic CPU fallback
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu

# 5. Query and list all archived simulation runs in DuckDB
julia --project=. ParticleTrackingRun.jl --list-runs

# 6. Multi-scenario comparative analytics across climate projections
julia --project=. ParticleTrackingRun.jl --compare-scenarios

# 7. Compute ensemble model-averaged connectivity matrix (P_ij ± σ)
julia --project=. ParticleTrackingRun.jl --model-average

# 8. Export all DuckDB tables to Apache Parquet format
julia --project=. ParticleTrackingRun.jl --export-parquet

# 9. Modify parameters and export new configuration file
julia --project=. ParticleTrackingRun.jl --particles=1000 --min-depth=120.0 --save-config=inputs/deep_shelf.config
```


---

## Decoupled Multi-Cohort Operational Workflow via CLI

Because Lagrangian tracking is thousands of times faster than solving 3D Navier-Stokes
equations, `ParticleTrackingRun.jl` allows decoupling the Eulerian hydrodynamic model
from larval transport. You can compute regional circulation once, archive the flow field
to `outputs/baseline/simulation_flow.jld2`, and then execute multiple Lagrangian cohorts
with different biological parameters, hatch dates, release modes, or random seeds.

### 1. Execute Hydrodynamics Once (--hydro-only)
Compute the 3D physical ocean state (velocity $\boldsymbol{u}$, temperature $T$, and
bathymetry) and write the checkpoint directly to `--hydro-model`:

```bash
julia --project=. ParticleTrackingRun.jl \
    --hydro-only \
    --hydro-model=outputs/baseline/hydrodynamics_ssp245_2050.jld2 \
    --duration=120.0 \
    --grid=50,50,10 \
    --scenario=ssp245 \
    --year=2050 \
    --output-dir=outputs/baseline
```

### 2. Track Cohort A: Spring Benthic Release with Active Ascent (--track-only)
Using the pre-computed flow field, track 500 larvae released from
the seabed that actively swim upward toward the surface:

```bash
julia --project=. ParticleTrackingRun.jl \
    --track-only \
    --hydro-model=outputs/baseline/hydrodynamics_ssp245_2050.jld2 \
    --run-id=cohort_spring_ascent \
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

### 3. Track Cohort B: Alternate Ascent Speed (15 mm/s)
Simulate larvae with a faster vertical ascent velocity ($w_{\text{ascent}} = 15\text{ mm/s}$):

```bash
julia --project=. ParticleTrackingRun.jl \
    --track-only \
    --hydro-model=outputs/baseline/hydrodynamics_ssp245_2050.jld2 \
    --run-id=cohort_fast_ascent \
    --output-dir=outputs/baseline \
    --particles=500 \
    --track-duration=10.0 \
    --release-mode=bottom \
    --ascent \
    --ascent-speed=0.015 \
    --ascent-target=-10.0 \
    --seed=102 \
    --duckdb \
    --db-path=outputs/particle_tracking.duckdb
```

### 4. Track Cohort C: Surface Release Control
Simulate a control cohort released directly into the surface mixed layer without initial
vertical ascent:

```bash
julia --project=. ParticleTrackingRun.jl \
    --track-only \
    --hydro-model=outputs/baseline/hydrodynamics_ssp245_2050.jld2 \
    --run-id=cohort_surface_control \
    --output-dir=outputs/baseline \
    --particles=500 \
    --track-duration=10.0 \
    --release-mode=surface \
    --no-ascent \
    --seed=201 \
    --duckdb \
    --db-path=outputs/particle_tracking.duckdb
```

### 5. Cross-Cohort Analysis & DuckDB Querying
All runs automatically log their trajectories, demographic connectivity matrices ($P_{ij}$),
recruitment metrics, and empirical movement fields into `outputs/particle_tracking.duckdb`.
Inspect and compare cohorts directly from the command line:

```bash
# List all completed runs in DuckDB
julia --project=. ParticleTrackingRun.jl --list-runs --db-path=outputs/particle_tracking.duckdb

# Multi-cohort comparative metrics (retention, PLD, degree-days, displacement)
julia --project=. ParticleTrackingRun.jl --compare-scenarios --db-path=outputs/particle_tracking.duckdb

# Bayesian / ensemble model-averaged demographic connectivity (P_ij ± σ)
julia --project=. ParticleTrackingRun.jl --model-average --db-path=outputs/particle_tracking.duckdb

# Export all DuckDB tables to Apache Parquet for external Python / R / BSTM analysis
julia --project=. ParticleTrackingRun.jl --export-parquet --db-path=outputs/particle_tracking.duckdb
```

---

## Centralized Configuration Management

All physical, numerical, and biological parameters can be saved or loaded from centralized
TOML configuration files:

```bash
# Execute run using settings from custom configuration file
julia --project=. ParticleTrackingRun.jl --all --config=inputs/ParticleTracking.config

# Modify parameters on CLI and export new validated configuration file
julia --project=. ParticleTrackingRun.jl \
    --particles=1000 \
    --min-depth=120.0 \
    --ascent-speed=0.012 \
    --save-config=inputs/deep_shelf.config
```

---

## Hardware Acceleration

Oceananigans automatically targets NVIDIA CUDA GPUs when available, falling back
gracefully to multi-threaded CPU architectures:

```bash
# GPU execution with automatic CPU fallback
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu

# Explicit CPU execution
julia --project=. ParticleTrackingRun.jl --all --cpu
```

---

## Unified Multi-Year Snow Crab Dispersal Platform

The snow crab modeling platform is integrated directly into [`ParticleTrackingRun.jl`](file:///c:/home/jae/projects/ParticleTracking/ParticleTrackingRun.jl) via the `--snowcrab-settings` flag (or shorthand `--snowcrab`). Passing this flag loads all calibrated physical and biological parameters for snow crab (*Chionoecetes opilio*):
- **Cohort Scale**: 500 larvae (quick: 50).
- **Pelagic Duration**: 60.0 days (quick: 5 days).
- **Seabed Placement & Depth**: Commercial nursery grounds ($\ge 100\text{ m}$) within $100\text{ km}$ CFA buffer.
- **Vertical Behavior**: Active post-hatch ascent ($10\text{ mm/s}$ to $-10\text{ m}$), stage-dependent DVM, and BBL shear attenuation.
- **Thermal Biology**: Molting base $T_{\text{base}} = -1.5^\circ\text{C}$ with degree-day thresholds (65, 130, 200 DD).
- **Domain & Grid**: High-resolution $100\times 100\times 20$ grid across Scotian Shelf and continental slope ($z \in [-3500, 0]\text{ m}$).
- **Analytical Database**: Persists cohorts to `outputs/snowcrab_tracking.duckdb`.

Any additional CLI flags passed alongside `--snowcrab-settings` directly **override** those defaults (e.g. `--particles=200`, `--ascent-speed=0.015`, `--hydro-model=custom_hydro.jld2`).

For programmatic Julia scripting, the exported function `SnowCrabRunOptions(; kwargs...)` generates the identical runtime configuration struct.

### Invocation Examples

```bash
# Production Snow Crab run with 5-year hydrodynamic simulation:
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --real-5yr --hydro-only --hydro-model=hydrodynamics1.jld2

# Quick test run:
julia --project=. ParticleTrackingRun.jl --snowcrab --all --quick
```

| Flag / Option | Argument Type | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `--all` | Flag | — | Run both Real 5-Year and 2-Year Climatological models + comparison. |
| `--real-5yr` | Flag | — | Run only the high-resolution 5-year real physical cycle. |
| `--climatology-2yr` | Flag | — | Run only the 2-year recent climatological average model. |
| `--compare` | Flag | — | Query DuckDB database and print comparative analytics across runs. |
| `--quick`, `-q` | Flag | — | Rapid testing mode with reduced domain resolution and duration. |
| `--hydro-model=<path>` | String | `outputs/hydrodynamics_...` | Hydrodynamic JLD2 snapshot file to save into or load from. |
| `--hydro-only` | Flag | — | Run hydrodynamics only and save snapshot to `--hydro-model`. |
| `--track-only` | Flag | — | Skip hydrodynamics and track larvae using flow fields from `--hydro-model`. |
| `--reuse-hydro` | Flag | — | Reuse `--hydro-model` if file exists on disk; otherwise run simulation. |
| `--run-id=<string>` | String | Auto-generated | Custom run identifier for DuckDB persistence and figure filenames. |
| `--particles=<int>` | Int | `500` (`50` quick) | Number of larvae per released cohort. |
| `--grid=<nx,ny,nz>` | Int,Int,Int | `100,100,20` | Spatial grid cell dimensions. |
| `--sim-dt=<sec>` | Float | `120.0` | Hydrodynamic integration time step. |
| `--heat-flux=<val>` | Float | `50.0` | Summer atmospheric downward heat flux ($W/m^2$). |
| `--release-mode=<mode>`| Symbol | `bottom` | Larval release depth mode (`bottom`, `range`, `surface`). |
| `--ascent` | Flag | `true` | Enable post-hatch active vertical ascent toward the surface. |
| `--no-ascent` | Flag | — | Disable post-hatch vertical ascent. |
| `--ascent-speed=<val>` | Float (m/s) | `0.010` | Vertical swimming speed during post-hatch ascent. |
| `--ascent-target=<val>`| Float (m) | `-10.0` | Target depth in meters for ascent completion. |
| `--db-path=<path>` | String | `outputs/snowcrab_tracking.duckdb` | Custom DuckDB analytical database path. |
| `--output-dir=<dir>` | String | `outputs` | Target directory for outputs and figures. |
| `--gpu`, `--cuda` | Flag | — | Enable GPU acceleration. |
| `--fallback-cpu` | Flag | `true` | Fallback to CPU if CUDA GPU is not detected. |

### Decoupled Hydrodynamics & Multi-Cohort Larval Tracking Workflow

Lagrangian particle tracking runs significantly faster than 3D Navier-Stokes integration.
Using `--hydro-only` and `--track-only` allows generating a hydrodynamic solution once and
reusing it across multiple larval dispersal cohorts:

```bash
# Step 1: Run 5-year hydrodynamic simulation only and save to a chosen file:
julia --project=. ParticleTrackingRun.jl --snowcrab --real-5yr --hydro-only --hydro-model=hydrodynamics1.jld2

# Step 2: Track Cohort 1 (Spring benthic release with active ascent):
julia --project=. ParticleTrackingRun.jl --snowcrab --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_baseline --particles=500 --release-mode=bottom --ascent

# Step 3: Track Cohort 2 (Fast vertical ascent at 15 mm/s):
julia --project=. ParticleTrackingRun.jl --snowcrab --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_fast_ascent --particles=500 --ascent-speed=0.015

# Step 4: Track Cohort 3 (Surface release control, no initial ascent):
julia --project=. ParticleTrackingRun.jl --snowcrab --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_surface_ctrl --particles=500 --release-mode=surface --no-ascent

# Step 5: Compare recruitment and connectivity across cohorts in DuckDB:
julia --project=. ParticleTrackingRun.jl --compare
```

