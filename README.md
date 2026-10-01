# ParticleTracking.jl

**Individual-based larval transport, biophysical ocean circulation, and demographic connectivity on the Scotian Shelf**

[![Julia](https://img.shields.io/badge/Julia-1.13-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Backend: Oceananigans.jl](https://img.shields.io/badge/Physics-Oceananigans.jl-informational.svg)](https://github.com/CliMA/Oceananigans.jl)
[![Storage: DuckDB](https://img.shields.io/badge/Storage-DuckDB-yellow.svg)](https://duckdb.org)
[![Visualization: CairoMakie](https://img.shields.io/badge/Visualization-Makie-purple.svg)](https://makie.juliaplots.org)

> Developed and verified on Julia 1.13. Oceananigans currently targets Julia 1.12; it emits a
> warning on newer versions. Everything here is tested on 1.13.1.

---

## Overview

`ParticleTracking.jl` couples regional 3D hydrostatic Boussinesq ocean circulation (Oceananigans.jl)
with individual-based stochastic Lagrangian particle tracking, snow crab larval thermal ecology, and
demographic connectivity analytics.

Parameterised for the **Scotian Shelf snow crab (*Chionoecetes opilio*)** fishery across Crab
Fishing Areas (CFAs 4X, 20–22, 23–24), the platform is modular and generalisable to other species
or shelf seas.

```
WOA23 / ETOPO 2022 / ERA5–CDS wind          (keyless, cached under inputs/)
                   │
                   ▼
3D Hydrostatic Circulation — Oceananigans.jl, CUDA GPU or CPU
  ├── Immersed boundary shelf bathymetry, Large & Pond (1981) wind drag
  ├── Air–sea surface heat flux & bottom drag (linear Rayleigh + quadratic)
  ├── Astronomical M2 + S2 spring–neap tides, generalized Simpson–Hunter fronts
  ├── Two-segment vertically stretched grid (surface-refined 0–400 m, coarse below)
  └── Climate scenarios: historical, SSP1-2.6, SSP2-4.5, SSP5-8.5, marine heatwave
                   │
                   ▼
Individual-Based Lagrangian Tracking (Euler–Maruyama SDE)
  ├── Strict marine placement, benthic release with bottom offset
  ├── Post-hatch vertical ascent toward the surface mixed layer
  ├── Logarithmic bottom boundary layer shear & passive gravitational sinking
  ├── Visser (1997) diffusive pseudo-drift correction at the pycnocline
  ├── Stage-specific DVM with CIL boundaries & turbidity attenuation
  ├── Degree-day ontogenetic molting via a per-larva developmental quantile
  ├── Tidally filtered benthic settlement suitability & thermal mortality
  └── Three persistent per-larva traits: developmental rate, vigour, settlement readiness
                   │
                   ▼
Demographic Connectivity & Analytics — DuckDB
  ├── CFA polygon boundary classification (ray casting)
  ├── Survival-weighted connectivity matrices (P_ij = Σ S_p / N_released)
  ├── Embedded DuckDB storage, scenario SQL queries & ensemble averaging
  └── CairoMakie figures, time–depth diagnostics & a self-contained Leaflet HTML map
```

---

## Quickstart

### Install

```bash
git clone <this-repo> && cd ParticleTracking
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### Run the test suite

```bash
julia --project=. test/runtests.jl          # 445 tests, ~3 min
```

### Run the pipeline

```bash
# Fast debug run: coarse grid, short hydrodynamics, small cohort (~3 min end to end)
julia --project=. ParticleTrackingRun.jl --all --quick --config=configs/default.toml

# Production run against real data
julia --project=. ParticleTrackingRun.jl --all --config=configs/snowcrab.toml

# Keyless open-data run (no credentials required)
julia --project=. ParticleTrackingRun.jl --all --config=configs/opendata.toml

# GPU, falling back to CPU if unavailable
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu --config=configs/snowcrab.toml
```

### Redraw every figure without re-simulating

```bash
julia --project=. ParticleTrackingRun.jl --figures --config=configs/snowcrab.toml
```

Reads the trajectory checkpoint and hydrodynamics archive and redraws the whole figure set —
particle figures, the larval-biology diagnostics, the Eulerian set, the interactive HTML map and the
animation. This is the fast path when iterating on a plotting function.

### Run individual segments

```bash
# Hydrodynamics only, then reuse for several cohorts
julia --project=. ParticleTrackingRun.jl --sim  --config=configs/snowcrab.toml
julia --project=. ParticleTrackingRun.jl --track --config=configs/snowcrab.toml --reuse-hydro
julia --project=. ParticleTrackingRun.jl --track --config=configs/snowcrab.toml --reuse-hydro

# Tracking + analytics against an existing hydrodynamics archive
julia --project=. ParticleTrackingRun.jl --track-only --config=configs/snowcrab.toml
```

Segment flags: `--data`, `--grid`, `--model`, `--climate`, `--sim`, `--track`, `--metrics`, `--viz`.

Useful modifiers: `--gpu` / `--cpu`, `--quick`, `--no-tides`, `--no-obc`, `--no-molting`,
`--animate-hydro`, `--interactive`, `--no-duckdb`, `--force-new`, `--restart`, `--allow-analytical-fallback`.

Run `julia --project=. ParticleTrackingRun.jl --help` for the full list, or see
**[CLI & Workflow Guide](ParticleTrackingRun.md)**.

### Analytics

```bash
julia --project=. ParticleTrackingRun.jl --list-runs
julia --project=. ParticleTrackingRun.jl --compare-scenarios
julia --project=. ParticleTrackingRun.jl --model-average
```

---

## Configuration

All physical, biological and numerical parameters live in **TOML** files. Nothing
species-specific is hard-coded in the core or exposed as a CLI flag.

| File | Purpose |
|---|---|
| `configs/default.toml` | Small synthetic case for development and tests |
| `configs/snowcrab.toml` | Production Scotian Shelf configuration, real data |
| `configs/opendata.toml` | Keyless open-data configuration (WOA23, ETOPO, wind) |

A configuration has 19 sections: `metadata`, `data`, `climate`, `domain`, `grid`, `bathymetry`,
`atmosphere`, `boundaries`, `tides`, `bottom_boundary_layer`, `hydrodynamics`, `tessellation`,
`biology`, `dvm`, `molting_and_settlement`, `storage`, `hardware`, `visualization`, `paths`.

Before each run the fully resolved configuration is written to
`<output_dir>/resolved_config.toml`, so every result carries its own provenance.

**Stochasticity.** Three coefficients of variation control the per-larva dispersion, each as a
*dispersion-only* knob that never shifts an expected outcome:

| Key | Meaning |
|---|---|
| `cv_molt` | spread of the developmental quantile (lognormal on degree-day thresholds) |
| `cv_mortality` | spread of the mortality-rate multiplier (lognormal frailty) |
| `cv_settlement` | spread of the settlement propensity (log-odds perturbation of the HSI) |

The `0.25` defaults are placeholders, not fitted values.

**Vertical grid & adaptive CFL.** `vertical_stretching_mode = "two_segment"` splits the column at
`vertical_break_depth`, giving a surface-refined upper segment (e.g. $\Delta z \approx 10\text{ m}$)
and a coarse lower segment (e.g. $\Delta z \approx 300\text{ m}$). Adaptive CFL time stepping evaluates
vertical stability locally across each layer ($\text{CFL}_z = \max_k [|w_k| \Delta t / \Delta z_k]$)
to prevent artificial time step collapse from deep vertical motion.

**Boundary forcing & sponge relaxation.** Interior astronomical tides are forced as momentum
accelerations $\boldsymbol{F}_{\text{tide}}$ compensating bottom drag, while open boundary relaxation
targets physical tidal velocity vectors $\boldsymbol{u}_{\text{target}}$ ($\text{m s}^{-1}$) across
active ocean boundaries, with landward boundaries masked to avoid spurious wave reflection.

---

## Data

| Data | Source | Credentials |
|---|---|---|
| Hydrography (T/S/O₂) | WOA23 0.25° (T/S) and 1.00° (O₂) | none |
| Bathymetry | ETOPO 2022 via NOAA CoastWatch ERDDAP (`ETOPO_2022_v1_15s`) | none |
| Surface wind | Open-Meteo → NOAA PSL → CDS | none for the first two |
| Reanalysis fields | Copernicus CDS (ERA5) | CDS API key in `~/.cdsapirc` |
| GLORYS reanalysis | Copernicus Marine ARCO Zarr | none (anonymous) |

Fetched files are cached under `inputs/`, so repeat runs are offline.

---

## Visualisation

**Larval** (`ParticleTrackingRun.jl --all`)
trajectory maps · DVM depth profiles · settlement density · empirical movement field · connectivity
matrix · thermal exposure · recruitment summary · **molt progression** (cohort transition CDF,
realised stage composition, molting degree-day vs developmental quantile) · **degree-day growth**
(per-larva thermal time and temperature history against the 65/130/200 thresholds) · **survival
curves** (survival probability, final cohort disposition, stage-transition mortality)

**Eulerian** (hydrodynamics)
advection · tracers · stratification · diffusion · vertical section · station time series ·
**time–depth diagram** (Hovmöller) · **near-bed temperature map** · **T–S diagram**

**Interactive** — a self-contained Leaflet HTML map with trajectory playback, layer toggles and a
time slider. The payload is decimated (default 400 tracks × 900 points) so it opens in a browser at
production scale; summary statistics are computed from the full-resolution arrays.

**Animation** — `plot_hydrodynamic_field` and `plot_hydrodynamic_dashboard`, GIF or MP4.

---

## Project layout

```
ParticleTrackingRun.jl   CLI driver; 8 segments, one function each
src/
  config/                TOML loading, validation, resolved-config provenance
  data/                  WOA23 / ETOPO / wind / GLORYS acquisition, grid construction
  model/                 hydrodynamics, climate scenarios, simulation loop
  biology/               DVM, molting, mortality, settlement
  output/                figures, interactive map, animation, exports
  utils/                 drag laws, coordinate helpers
test/runtests.jl         445 tests in 21 testsets
configs/                 default.toml, snowcrab.toml, opendata.toml
inputs/                  cached downloaded data
work/                    per-run output (DuckDB, hydrodynamics archive, figures)
```

---

## References

- 📖 **[CLI & Workflow Guide](ParticleTrackingRun.md)** — full command reference.
- 📝 **[TODO / engineering notes](todo.md)** — current state, known limitations and open work.

Method references: Oceananigans.jl · ClimaOcean.jl · NumericalEarth.jl ·
Loder & Petrie (1991) shelf currents · Large & Pond (1981) wind drag · Visser (1997) pycnocline
drift · Simpson & Hunter (2004) fronts · Sainte-Marie & Sainte-Marie (1999) stage durations ·
Kuhn & Choi (2011) degree-day thresholds.

---

## License

Released under the [MIT License](LICENSE).
