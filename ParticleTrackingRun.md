# ParticleTrackingRun.jl — CLI & Workflow Execution Guide

## Overview

This workflow establishes a regional 3D hydrodynamic model using
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) coupled to individual-based
Lagrangian particle tracking. The hydrodynamic solution provides time-evolving advection
($\boldsymbol{u} = (u, v, w)$) and turbulent diffusion ($\kappa_h, \kappa_v$) fields, which
drive Lagrangian tracking of pelagic and semiplanktonic larvae undergoing stage-specific
vertical migration, temperature-dependent development, and spatial dispersal across
coastal nursery habitats.

The engine is **species-agnostic and region-agnostic**: all biological, ecological, and
domain-specific parameters are supplied through TOML configuration files rather than
compiled-in defaults or species-specific CLI flags.

---

## Quick Start

```bash
# Instantiate dependencies
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Fast end-to-end debug verification (coarse grid, short durations)
julia --project=. ParticleTrackingRun.jl --all --quick

# Production run driven by a TOML configuration
julia --project=. ParticleTrackingRun.jl --all --config=work/snowcrab/snowcrab.toml

# Print full CLI help
julia --project=. ParticleTrackingRun.jl --help
```

---

## Architecture: Eight-Segment Pipeline

The runner decomposes the modelling workflow into eight independently addressable
segments. Each may be run alone, or as a complete pipeline via `--all`.

| # | Segment | Purpose |
|---|---------|---------|
| 1 | `data`     | Environmental data ingestion (bathymetry, winds, hydrography) and drag processing |
| 2 | `grid`     | Spherical grid construction and immersed boundary (bathymetry) regridding |
| 3 | `model`    | Hydrodynamic model assembly, tidal forcing, OBC and atmospheric fluxes |
| 4 | `climate`  | Climate scenario integration and thermal ecology |
| 5 | `sim`      | Oceananigans hydrodynamic time integration with checkpointing |
| 6 | `track`    | Lagrangian particle tracking with DVM, molting, and settlement |
| 7 | `metrics`  | Empirical movement, recruitment, connectivity, Voronoi tessellation |
| 8 | `viz`      | Scientific visualizations and spatial figures |

---

## CLI Flags Reference

### Execution Modes

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--all` | Flag | — | Execute the complete 8-segment production pipeline |
| `--quick`, `-q` | Flag | — | Fast debug mode (coarse grid, ~1 h hydrodynamics, ~2 d tracking) |
| `--segment=<name>` | String | `all` | Run a single segment: `data`, `grid`, `model`, `climate`, `sim`, `track`, `metrics`, `viz`, `all` |
| `--help`, `-h` | Flag | — | Print CLI help |

### Decoupled Hydrodynamics & Multi-Cohort Tracking

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--hydro-model=<path>` | String | `<output_dir>/hydrodynamics.jld2` | Hydrodynamic JLD2 file (output for `--hydro-only`, input for `--track-only`) |
| `--hydro-only` | Flag | — | Run Segments 1–5 only; save flow fields to `--hydro-model` |
| `--track-only` | Flag | — | Run Segments 6–8 only; read flow fields from `--hydro-model` |
| `--reuse-hydro` | Flag | — | Reuse existing `--hydro-model` if complete; otherwise integrate |
| `--run-id=<string>` | String | `run_<scenario>_<year>` | Cohort identifier for DuckDB persistence and figure naming |
| `--restart` | Flag | `true` | Resume from the latest checkpoint if an incomplete run exists |
| `--no-restart`, `--force-new` | Flag | — | Ignore existing checkpoints and restart from $t=0$ |
| `--checkpoint` | Flag | `true` | Enable prognostic state checkpointing |
| `--no-checkpoint` | Flag | — | Disable prognostic state checkpointing |
| `--checkpoint-interval=<spec>` | String | output cadence | Checkpoint interval: seconds (`3600`), hours (`6h`), or days (`1d`) |
| `--checkpoints-dir=<path>` | String | `<output_dir>/checkpoints` | Checkpoint archive directory |
| `--checkpoint-prefix=<name>` | String | `checkpoint` | Checkpoint filename prefix |
| `--checkpoint-cleanup` | Flag | — | Retain only the latest checkpoint on disk |
| `--no-checkpoint-cleanup` | Flag | `true` | Keep all intermediate checkpoints |

### Configuration

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--config=<path>` | String | `inputs/ParticleTracking.toml` | Load parameters from a TOML configuration file |
| `--save-config[=<path>]` | String | — | Export the active resolved configuration to TOML and exit |

**All species-specific and domain-specific parameters belong in the TOML file**, not the CLI.
The runner has no `--snowcrab*`, `--real-5yr`, or `--climatology-2yr` flags; select a
scenario by pointing `--config` at the corresponding TOML file.

### DuckDB Analytics

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--duckdb` | Flag | `true` | Enable DuckDB archiving |
| `--no-duckdb` | Flag | — | Disable DuckDB archiving |
| `--db-path=<path>` | String | `outputs/particle_tracking.duckdb` | DuckDB database path |
| `--list-runs` | Flag | — | Print all archived simulation runs |
| `--compare-scenarios` | Flag | — | Print multi-scenario comparative analytics |
| `--model-average` | Flag | — | Compute ensemble model-averaged connectivity and recruitment |

### Hardware

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--gpu`, `--cuda` | Flag | — | Enable NVIDIA CUDA GPU acceleration |
| `--cpu` | Flag | `true` | Execute on multi-threaded CPU |
| `--fallback-cpu` | Flag | — | Fall back to CPU automatically if CUDA is unavailable |

### Spatial Domain & Grid

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--lon=<min,max>` | String | from config | Longitude bounds (°E) |
| `--lat=<min,max>` | String | from config | Latitude bounds (°N) |
| `--depth-range=<min,max>` | String | from config | Vertical depth range (m) |
| `--grid=<nx,ny,nz>` | String | from config | Grid cell dimensions |
| `--nx=<int>`, `--ny=<int>`, `--nz=<int>` | Int | from config | Individual grid axis dimensions |
| `--res-scale=<float>` | Float | `1.0` | Resolution scaling factor |
| `--stretched-z` | Flag | — | Enable hyperbolic-tangent vertical stretching |
| `--uniform-z` | Flag | — | Use uniform vertical layer thicknesses |
| `--z-file`, `--vertical-grid-file` | String | `inputs/scotian_shelf_vertical_grid.csv` | Vertical coordinate CSV |

### Environmental Data & Forcing

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--real` | Flag | — | Fetch real bathymetry and atmospheric forcing |
| `--synthetic` | Flag | `true` | Generate idealized synthetic shelf data |
| `--tides` | Flag | `true` | Enable astronomical tidal body forcing |
| `--no-tides` | Flag | — | Disable tidal forcing |
| `--tidal-u=<val>` | Float (m/s) | `0.25` | M2 semi-major tidal velocity amplitude |
| `--tidal-v=<val>` | Float (m/s) | `0.12` | M2 semi-minor tidal velocity amplitude |
| `--era5-forcing`, `--era5` | Flag | — | Enable ERA5 atmospheric surface forcing |
| `--no-era5` | Flag | — | Use analytical wind stress forcing |
| `--obc` | Flag | — | Enable GLORYS open boundary conditions |
| `--no-obc` | Flag | — | Disable open boundary conditions |
| `--obc-source` | String | `glorys12v1` | Open boundary data source |
| `--obc-type` | String | `flather_chapman` | Open boundary condition scheme |

### Climate Scenarios

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--scenario=<name>` | String | `ssp245` | `historical`, `ssp126`, `ssp245`, `ssp585`, `mhw`, `climatology` |
| `--year=<int>` | Int | `2050` | Climate projection horizon year |
| `--heat-flux=<val>` | Float (W/m²) | `50.0` | Net atmospheric surface heat flux |

### Hydrodynamic Simulation

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--duration=<hours>` | Float | from config | Simulation duration in hours |
| `--sim-duration=<sec>` | Float | from config | Simulation duration in seconds |
| `--sim-dt=<sec>` | Float | `120.0` | Initial hydrodynamic time step |
| `--adaptive-cfl` | Flag | `true` | Enable CFL-limited adaptive time stepping |
| `--no-adaptive-cfl` | Flag | — | Disable adaptive CFL stepping |
| `--target-cfl=<val>` | Float | `0.2` | Target advective Courant–Friedrichs–Lewy number |

### Lagrangian Particle Tracking & Larval Ecology

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--particles=<int>` | Int | `100` | Cohort size |
| `--track-duration=<days>` | Float | `5.0` | Tracking duration in days |
| `--track-dt=<sec>` | Float | `300.0` | Lagrangian integration time step |
| `--min-depth=<m>` | Float | `100.0` | Minimum water depth for larval placement |
| `--buffer-km`, `--buffer`, `--buf` | Float | `100.0` | Spatial buffer beyond stratum boundaries (km) |
| `--dvm` | Flag | `true` | Enable Diel Vertical Migration |
| `--no-dvm` | Flag | — | Disable DVM |
| `--molting` | Flag | `true` | Enable degree-day molting and mortality |
| `--no-molting` | Flag | — | Disable thermal molting |
| `--diff-h`, `--diffusivity-h` | Float (m²/s) | `10.0` | Horizontal turbulent diffusivity |
| `--diff-v`, `--diffusivity-v` | Float (m²/s) | `1e-4` | Vertical turbulent diffusivity |
| `--release-mode` | String | `bottom` | `bottom`, `range`, `surface` |
| `--ascent` | Flag | `true` | Enable post-hatch vertical ascent |
| `--no-ascent` | Flag | — | Disable initial vertical ascent |
| `--ascent-speed` | Float (m/s) | `0.010` | Ascent swimming speed |
| `--ascent-target` | Float (m) | `-10.0` | Target depth for ascent completion |

### Voronoi Tessellation

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--voronoi-units`, `--n-units` | Int | `5000` | Number of depth-stratified Voronoi units |
| `--voronoi-prob-core` | Float | `0.8` | Core stratum sampling probability |
| `--voronoi-prob-shallow` | Float | `0.1` | Shallow stratum sampling probability |
| `--voronoi-prob-deep` | Float | `0.1` | Deep stratum sampling probability |
| `--voronoi-min-core` | Float (km) | `1.5` | Minimum separation, core stratum |
| `--voronoi-min-shallow` | Float (km) | `5.0` | Minimum separation, shallow stratum |
| `--voronoi-min-deep` | Float (km) | `10.0` | Minimum separation, deep stratum |

### Visualization & Animation

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--interactive` | Flag | `true` | Export interactive HTML5 Leaflet map |
| `--no-interactive` | Flag | — | Disable interactive HTML map |
| `--animate-hydro`, `--anim-hydro` | Flag | — | Render MP4/GIF animation of hydrodynamic fields |
| `--no-animate-hydro` | Flag | `true` | Disable hydrodynamic animation |
| `--anim-variable=<name>` | String | `dashboard` | Diagnostic field (see table below) |
| `--anim-fps=<int>` | Int | `10` | Playback framerate |
| `--anim-format` | String | `mp4` | `mp4` or `gif` |
| `--anim-depth=<m>` | Float | `-2.5` | Target depth for 2D field slices |
| `--anim-output`, `--anim-file` | String | `""` | Custom output video path |
| `--anim-overlay-particles`, `--anim-particles` | Flag | — | Overlay Lagrangian particles on animation |

### I/O & Environment

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--output-dir=<path>` | String | from config | Output directory for figures and datasets |
| `--input-dir=<path>` | String | from config | Input cache directory |
| `--seed=<int>` | Int | `42` | Random number generator seed |

---

## Configuration Files (TOML)

Configuration files use [TOML](https://toml.io) syntax and are validated against a typed
schema (`src/config_schema.jl`) via **Configurations.jl**. Files are located by
`--config=<path>`; the default is resolved by `find_default_config_path()`.

### Available Configurations

| File | Description |
| :--- | :--- |
| `inputs/ParticleTracking.toml` | Default regional modelling configuration |
| `inputs/hydrodynamics_climatology_2yr.toml` | 2-year climatological circulation (repeating annual cycle) |
| `inputs/hydrodynamics_2020_2yr.toml` | 2020–2022 real-data hindcast |
| `inputs/hydrodynamics_mhw_2yr.toml` | Marine heat wave scenario with thermal anomaly |
| `work/snowcrab/snowcrab.toml` | Snow crab larval dispersal; all outputs under `work/snowcrab/` |

### Section Reference

| Section | Purpose |
| :--- | :--- |
| `[metadata]` | Run name, description, author, version |
| `[data]` | Data ingestion mode and source dataset identifiers |
| `[climate]` | Scenario, projection/baseline/horizon years |
| `[domain]` | Study domain bounds (`lon_*`, `lat_*`, `z_*`) and buffer |
| `[grid]` | Grid dimensions, resolution scale, vertical stretching |
| `[bathymetry]` | Bathymetry provider, resolution, fallbacks |
| `[atmosphere]` | Atmospheric forcing source, variables, drag formulation |
| `[boundaries]` | OBC method, source, sponge parameters, embedding domain |
| `[tides]` | Tidal constituents, amplitudes, period, phase |
| `[bottom_boundary_layer]` | Quadratic drag and benthic mixing closure |
| `[hydrodynamics]` | Duration, time step, CFL control, tracers, model file |
| `[tessellation]` | Depth-stratified Voronoi unit parameters |
| `[biology]` | Cohort size, tracking duration, diffusivity, release, ascent |
| `[dvm]` | Stage-specific diel vertical migration depths |
| `[molting_and_settlement]` | Degree-day thresholds, settlement criteria, mortality |
| `[storage]` | Output filename, DuckDB path, checkpoint cadence and directory |
| `[hardware]` | GPU preference and CPU fallback |
| `[visualization]` | Interactive map toggle and figure title |
| `[paths]` | `output_dir`, `input_dir`, RNG `seed` |

### Two Domain Scales

Two distinct spatial scales are used and are declared separately in the TOML:

- **Study domain** (`[domain]`): the computational grid and analysis region — here the
  Scotian Shelf, `lon ∈ [-68.0, -57.0]°E`, `lat ∈ [42.0, 47.5]°N`.
- **Embedding domain** (`[boundaries]`): the larger Northwest Atlantic region
  (`lon ∈ [-71.0, -53.0]°E`, `lat ∈ [40.0, 48.5]°N`) from which open boundary and
  atmospheric forcing data are extracted. The study domain nests inside it.

### Path Resolution

`[paths] output_dir` is the single source of truth for all generated artefacts. Within it:

- `output_filename` — hydrodynamic field archive (`.jld2`)
- `duckdb_path` — analytical database
- `checkpoint_dir` — restart checkpoint archive

Set `checkpoint_dir = ""` to inherit `<output_dir>/checkpoints` automatically.

### Exporting Configuration

```bash
# Resolve CLI overrides and export to TOML
julia --project=. ParticleTrackingRun.jl \
    --particles=1000 --min-depth=120.0 --ascent-speed=0.012 \
    --save-config=inputs/deep_shelf.toml
```

---

## Decoupled Multi-Cohort Workflow

Lagrangian tracking is orders of magnitude cheaper than solving the 3D primitive
equations. Compute the circulation once, archive it, then run many cohorts against it.

```bash
# 1. Solve hydrodynamics once
julia --project=. ParticleTrackingRun.jl \
    --config=work/snowcrab/snowcrab.toml \
    --hydro-only \
    --hydro-model=work/snowcrab/hydrodynamics_base.jld2

# 2. Cohort A — benthic release with active ascent
julia --project=. ParticleTrackingRun.jl \
    --config=work/snowcrab/snowcrab.toml \
    --track-only --hydro-model=work/snowcrab/hydrodynamics_base.jld2 \
    --run-id=cohort_spring --ascent --ascent-speed=0.010 --seed=101

# 3. Cohort B — faster ascent
julia --project=. ParticleTrackingRun.jl \
    --config=work/snowcrab/snowcrab.toml \
    --track-only --hydro-model=work/snowcrab/hydrodynamics_base.jld2 \
    --run-id=cohort_fast --ascent --ascent-speed=0.015 --seed=202

# 4. Cohort C — surface release control
julia --project=. ParticleTrackingRun.jl \
    --config=work/snowcrab/snowcrab.toml \
    --track-only --hydro-model=work/snowcrab/hydrodynamics_base.jld2 \
    --run-id=cohort_control --release-mode=surface --no-ascent --seed=303

# 5. Compare all cohorts
julia --project=. ParticleTrackingRun.jl --compare-scenarios --db-path=work/snowcrab/snowcrab.duckdb
```

---

## Checkpointing and Graceful Restart

1. **Automatic periodic checkpointing** — Oceananigans `Checkpointer` serializes the full
   prognostic state ($u, v, w, T, S, \eta$, clock) into
   `<checkpoint_dir>/<prefix>_iteration<N>.jld2`. With `cleanup_checkpoints = true`, only
   the newest file is retained, preventing storage exhaustion on long integrations.
   On Windows platforms, file removal traps transient OS sharing locks (`EBUSY` / `IOError`)
   with explicit garbage collection (`GC.gc()`), deferring removal across iterations rather
   than crashing multi-day integrations.
2. **Emergency checkpointing on interruption** — `run_hydrodynamic_simulation!` traps
   `InterruptException` (`Ctrl+C`, SLURM walltime) and flushes state before exiting cleanly.
3. **Seamless resumption** — with `auto_restart = true` the runner inspects the diagnostic
   timeseries and checkpoint directory; if $t_{\text{last}} < t_{\text{stop}}$ it picks up
   from the newest checkpoint and appends new output without overwriting existing records.
   `--no-restart` / `--force-new` forces a fresh integration from $t=0$.

---

## Numerical Stability, Adaptive CFL & Boundary Forcing Assumptions

### 1. Layer-Resolved Vertical Courant–Friedrichs–Lewy (CFL) Condition
In vertically stretched coordinates (e.g. surface layer $\Delta z_1 \approx 10\text{ m}$ vs deep basin
$\Delta z_k \approx 300\text{ m}$), evaluating the vertical advective Courant number against a single
surface layer thickness $\min(\Delta z)$ introduces an artificial $30\times$ penalty for deep vertical
velocities. The adaptive CFL controller evaluates layer-resolved vertical grid spacing:

$$\text{CFL}_z = \max_{i, j, k} \left( \frac{|w_{i,j,k}| \Delta t}{\Delta z_k} \right)$$

where $\Delta z_k = z_{F, k+1} - z_{F, k}$ is the true vertical grid spacing at vertical cell $k$.
The overall 3D advective CFL number is computed as:

$$\text{CFL} = \left( \max \frac{|u|}{\Delta x} + \max \frac{|v|}{\Delta y} + \max_{k} \frac{|w_k|}{\Delta z_k} \right) \Delta t \le \text{CFL}_{\text{target}}$$

### 2. Barotropic Tidal Forcing vs Boundary Relaxation
Astronomical tides are implemented as momentum body accelerations in the interior,
and target velocity vectors in open boundary relaxation sponge layers:
- **Interior Momentum Acceleration**: $\boldsymbol{F}_{\text{tide}} = (F_u, F_v)$ in $\text{m s}^{-2}$,
  compensating bottom friction $r_{\text{drag}}$.
- **Selective Open Boundary Sponge Masking**: Sponge damping layers are applied strictly to open
  ocean boundaries (e.g., `:east`, `:south`, `:west`), masking land-adjacent boundaries (e.g., `:north`)
  to prevent artificial inflow generation and numerical boundary reflection.

### 3. Hydrodynamic Field Serialization & Part Sizing
Hydrodynamic fields are written to JLD2 archives according to `output_schedule_seconds` in the `[hydrodynamics]` section (defaulting to 21,600 s / 6 h). Intermediate segment runs serialize into segmented files under `parts/` (e.g. `hydrodynamics_..._part1.jld2`). Non-advective diagnostic fields ($\nu, \kappa, N^2, \zeta$) can be selectively included via `include_diagnostics` in `setup_hydrodynamic_simulation` to optimize storage for particle tracking.

```bash
# Resume after interruption (default behaviour)
julia --project=. ParticleTrackingRun.jl --config=inputs/hydrodynamics_2020_2yr.toml \
    --hydro-only --restart

# Force a clean re-run
julia --project=. ParticleTrackingRun.jl --config=inputs/hydrodynamics_2020_2yr.toml \
    --hydro-only --force-new
```

---

## Depth-Stratified Voronoi Tessellation

Refines demographic resolution in core nursery habitats without incurring prohibitive
Eulerian CFL penalties, by post-processing Lagrangian endpoints into adaptive areal units.

### Stratum Sampling Strategy

| Stratum | Depth range | Probability | Min. separation |
| :--- | :--- | :--- | :--- |
| Core nursery | 50–350 m | $P = 0.8$ | 1.5 km |
| Shallow / inshore | 0–50 m | $P = 0.1$ | 5.0 km |
| Deep slope / basin | > 350 m | $P = 0.1$ | 10.0 km |

$N$ units (default 5000; 200 in `--quick`) are sampled by Poisson-disc process across the
regional bathymetry. `slope_weighting` refines spacing along steep shelf breaks and canyons
by `slope_factor`.

### Usage

```bash
# Tessellation parameters come from the [tessellation] section of the TOML
julia --project=. ParticleTrackingRun.jl --config=work/snowcrab/snowcrab.toml --all

# Override on the CLI
julia --project=. ParticleTrackingRun.jl --config=work/snowcrab/snowcrab.toml --all \
    --voronoi-units=5000 --voronoi-prob-core=0.85 --voronoi-min-core=1.5
```

### Outputs

- **$N \times N$ connectivity matrix** $P_{ij}$ across Voronoi units
- **$3 \times 3$ macro-strata matrix** (shallow ↔ core ↔ deep)
- **GeoJSON polygons** — `<output_dir>/voronoi_units.geojson`
- **DuckDB table** `voronoi_units` — coordinates, stratum, depth, area, settlement density
- **Figure** `<output_dir>/voronoi_tessellation.png` — stratum-coloured cells with settlement density

---

## Visualization System

Segment 8 renders a full set of publication figures. All functions respect the resolved
`output_dir` from the TOML and use a shared masking helper so land and below-seafloor cells
render transparently rather than producing degenerate colour ranges.

### 2D Figures (CairoMakie)

| Artefact | Content |
| :--- | :--- |
| `larval_trajectories.png` | Trajectory map coloured by developmental stage |
| `dvm_depth_profiles.png` | Stage-specific diel depth distributions |
| `settlement_density.png` | Gridded benthic settlement density |
| `empirical_movement_field.png` | Quiver field of empirical displacement and diffusivity |
| `regional_connectivity_matrix.png` | Macro-stratum transition matrix $P_{ij}$ |
| `thermal_exposure_map.png` | Degree-day thermal exposure |
| `recruitment_summary.png` | Survival and recruitment bar summary |
| `hydrodynamic_advection.png` | Surface current speed with vector overlay |
| `hydrodynamic_tracers.png` | Temperature and salinity fields |
| `hydrodynamic_stratification.png` | Buoyancy frequency $N^2$ and salinity gradient |
| `hydrodynamic_diffusion.png` | Eddy diffusivity $\kappa_v$ and viscosity $\nu_v$ |
| `hydrodynamic_section.png` | Vertical cross-section with bathymetry masking |
| `multi_panel_dashboard.png` | 6-panel combined hydrodynamic + Lagrangian dashboard |
| `voronoi_tessellation.png` | Voronoi units by stratum with settlement density |
| `particle_fate_summary.png` | Stage progression, fate distribution, thermal exposure |

### 3D Visualizations (GLMakie, GPU-accelerated)

Exported through `ParticleTracking.GLVisualization`:

| Function | Content |
| :--- | :--- |
| `plot_3d_hydrodynamic_field` | Volume rendering with isosurfaces, bathymetry terrain mesh, optional particle overlay |
| `plot_3d_particle_trajectories` | 3D trajectory tubes with stage colour gradients and fading tails |
| `plot_3d_connectivity` | 3D connectivity network between management strata with arc weighting |

### Animation Diagnostic Variables

| Key | Aliases | Physical metric | Colormap |
| :--- | :--- | :--- | :--- |
| `temperature` | `T`, `temp`, `theta` | Potential temperature (°C) | `:thermal` |
| `advection` | `speed`, `velocity`, `u_h` | Horizontal current speed $\sqrt{u^2+v^2}$ (cm s⁻¹) | `:viridis` |
| `diffusion` | `diffusivity`, `kappa` | Vertical eddy diffusivity $\kappa_v$ ($10^{-4}$ m² s⁻¹) | `:turbid` |
| `salinity` | `S`, `sal` | Practical salinity (PSU) | `:haline` |
| `viscosity` | `nu`, `eddy_viscosity` | Vertical eddy viscosity $\nu_v$ ($10^{-4}$ m² s⁻¹) | `:deep` |
| `stratification` | `N2` | Buoyancy frequency squared $N^2$ ($10^{-4}$ s⁻²) | `:ice` |
| `density` | `rho` | Potential density (kg m⁻³) | `:dense` |
| `vorticity` | `zeta` | Relative vorticity $\zeta$ ($10^{-5}$ s⁻¹) | `:balance` |
| `elevation` | `eta`, `ssh` | Free surface height $\eta$ (cm) | `:delta` |
| `w` | `vertical_velocity` | Vertical velocity $w$ (mm s⁻¹) | `:curl` |
| `richardson` | `Ri` | Gradient Richardson number (dimensionless) | `:spectral` |
| `dashboard` | — | Synchronized multi-diagnostic dashboard | Multiple |

Depth targeting is continuous via nearest-vertical-coordinate matching
($k = \arg\min_k |z_k - z_{\text{target}}|$), so `--anim-depth=-50.0` targets the cold
intermediate layer and `--anim-depth=0.0` the surface mixed layer.

```bash
# 4-panel dashboard with larval drift overlay
julia --project=. ParticleTrackingRun.jl --config=work/snowcrab/snowcrab.toml \
    --animate-hydro --anim-variable=dashboard --anim-depth=-2.5 --anim-fps=12 \
    --anim-output=work/snowcrab/dashboard.mp4

# CIL temperature evolution
julia --project=. ParticleTrackingRun.jl --config=work/snowcrab/snowcrab.toml \
    --animate-hydro --anim-variable=temperature --anim-depth=-50.0 \
    --anim-output=work/snowcrab/cil_temperature.mp4

# Surface advection with particle overlay
julia --project=. ParticleTrackingRun.jl --config=work/snowcrab/snowcrab.toml \
    --animate-hydro --anim-variable=advection --anim-depth=0.0 --anim-particles
```

---

## Programmatic API

```julia
using ParticleTracking

# Load and validate a TOML configuration
cfg = load_configuration("work/snowcrab/snowcrab.toml")

# Resolve the runtime options struct
opts = configuration_to_options(cfg)

# Inspect or override individual options
opts.domain_lon, opts.grid_size, opts.n_particles

# Run individual segments
run_segment_data(; opts = opts)
run_segment_sim(; opts = opts)
run_segment_track(; opts = opts, trajectories = nothing)
```

The configuration layer is split into two coordinated representations:

- **`configuration.jl`** — nested `Dict` parsing, `HydrodynamicOptions` construction, and
  `HydrodynamicConfig` / `LarvalDispersalConfig` / `CoupledSimulationConfig` decoupled views.
- **`config_schema.jl`** — `Configurations.jl` `@option` structs providing typed,
  validated schema defaults; `load_config` / `save_config` / `schema_to_options` bridge
  the schema to runtime options.

---

## Common Workflows

```bash
# 1. Fast end-to-end debug verification
julia --project=. ParticleTrackingRun.jl --all --quick

# 2. Production climatological run (all output under work/snowcrab/)
julia --project=. ParticleTrackingRun.jl --all --config=work/snowcrab/snowcrab.toml

# 3. Marine heat wave scenario
julia --project=. ParticleTrackingRun.jl --all --config=inputs/hydrodynamics_mhw_2yr.toml

# 4. Real 2020–2022 hindcast
julia --project=. ParticleTrackingRun.jl --all --config=inputs/hydrodynamics_2020_2yr.toml

# 5. GPU-accelerated with automatic CPU fallback
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu \
    --config=work/snowcrab/snowcrab.toml

# 6. Query archived runs and comparative analytics
julia --project=. ParticleTrackingRun.jl --list-runs
julia --project=. ParticleTrackingRun.jl --compare-scenarios
julia --project=. ParticleTrackingRun.jl --model-average

# 7. Export a modified configuration
julia --project=. ParticleTrackingRun.jl --particles=1000 --min-depth=120.0 \
    --save-config=inputs/deep_shelf.toml
```

---

## License

Released under the [MIT License](LICENSE).
