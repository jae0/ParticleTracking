# ParticleTracking.jl

**A High-Performance Biophysical Ocean Modeling, Individual-Based Larval Transport, and Demographic Population Connectivity Framework**

[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Backend: Oceananigans.jl](https://img.shields.io/badge/Physics-Oceananigans.jl-informational.svg)](https://github.com/CliMA/Oceananigans.jl)
[![Storage: DuckDB](https://img.shields.io/badge/Storage-DuckDB-yellow.svg)](https://duckdb.org)

---

## Overview

`ParticleTracking.jl` couples regional 3D hydrostatic Boussinesq ocean circulation with individual-based stochastic Lagrangian particle tracking, larval thermal bioenergetics, and demographic connectivity analytics. 

Parameterized for the **Scotian Shelf snow crab (*Chionoecetes opilio*)** fishery ecosystem across Crab Fishing Areas (CFAs 20–22, 23–24, 4X), the platform is fully modular and generalizable to any marine species or regional shelf sea worldwide.

```
NOAA ERDDAP Data / Synthetic Benchmarks
                   │
                   ▼
3D Hydrodynamic Circulation (Oceananigans.jl on CUDA GPU / CPU)
  ├── Immersed boundary shelf bathymetry & Large & Pond (1981) wind drag
  ├── Air-sea surface heat flux (50 W/m²) & bottom drag (linear Rayleigh + quadratic)
  ├── Astronomical M2 + S2 spring-neap tidal forcing & generalized Simpson-Hunter fronts
  └── CMIP6 climate warming scenarios (Historical, SSP1-2.6, SSP2-4.5, SSP5-8.5, MHW)
                   │
                   ▼
Individual-Based Lagrangian Particle Tracking (Euler-Maruyama SDE)
  ├── Strict marine bathymetric placement & benthic release with bottom offset (0.5–3.0 m)
  ├── Active post-hatch vertical ascent swimming (10 mm/s) toward surface mixed layer
  ├── Logarithmic bottom boundary layer (BBL) shear & larval passive gravitational sinking
  ├── Visser (1997) diffusive pseudo-drift correction for stratified pycnoclines
  ├── Stage-specific Diel Vertical Migration (DVM) with CIL boundaries & turbidity attenuation
  ├── Calibrated thermal degree-day molting (T_base = -1.5°C; Zoea I -> II -> Megalopa)
  └── Tidally filtered benthic settlement suitability & exponential thermal mortality
                   │
                   ▼
Demographic Connectivity & Analytical Engine (DuckDB & Leaflet)
  ├── Administrative CFA polygon boundary classification (Jordan Curve ray-casting)
  ├── Stochastic survival-weighted recruitment connectivity matrices (P_ij = Σ S_p / N_released)
  ├── Embedded DuckDB analytical storage, scenario SQL querying & ensemble averaging
  └── CairoMakie publication charts & interactive HTML5 Leaflet.js dashboard
```

---

## Quickstart

### 1. Installation & Environment Setup
Clone the repository and instantiate Julia dependencies:
```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### 2. Run the Full Test Suite
Verify all 17 test sets (**618 unit tests**):
```bash
julia --project=. test/runtests.jl
```

### 3. Run the Production Pipeline
Execute an end-to-end simulation using the command-line interface:
```bash
# Fast debug mode (coarse grid, 1 hr hydro, 2 day track)
julia --project=. ParticleTrackingRun.jl --all --quick

# Full production run under CMIP6 SSP5-8.5 warming (Year 2050) with 500 larvae
julia --project=. ParticleTrackingRun.jl --all --scenario=ssp585 --year=2050 --particles=500

# GPU-accelerated run with automatic CPU fallback
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu
```

### 4. Decoupled Multi-Cohort Batching
Run the heavy 3D hydrodynamics once, then track multiple distinct cohorts (different
hatch dates, release depths, swimming speeds, or seeds) without re-solving fluid equations:
```bash
# Step 1: Solve hydrodynamics once and archive flow field
julia --project=. ParticleTrackingRun.jl --data --grid --model --sim --output-dir=outputs/baseline

# Step 2: Track Cohort A (Spring hatch, benthic release with ascent)
julia --project=. ParticleTrackingRun.jl --track --metrics --viz --output-dir=outputs/baseline \
    --particles=500 --release-mode=bottom --ascent --ascent-speed=0.010 --seed=101

# Step 3: Track Cohort B (Summer hatch, alternate ascent speed)
julia --project=. ParticleTrackingRun.jl --track --metrics --viz --output-dir=outputs/baseline \
    --particles=500 --release-mode=bottom --ascent --ascent-speed=0.015 --seed=201
```

---

## DuckDB Analytics & Scenario Management

All runs, trajectory time series, recruitment metrics, and demographic transition matrices
are archived in `outputs/particle_tracking.duckdb`:

```bash
# List all archived simulation runs
julia --project=. ParticleTrackingRun.jl --list-runs

# Multi-scenario comparative analytics
julia --project=. ParticleTrackingRun.jl --compare-scenarios

# Bayesian / ensemble model-averaged demographic connectivity (P_ij ± σ)
julia --project=. ParticleTrackingRun.jl --model-average

# Export DuckDB tables to Apache Parquet format
julia --project=. ParticleTrackingRun.jl --export-parquet
```

---

## Multi-Year Snow Crab Hydrodynamics & Dispersal (`--snowcrab-settings`)

The platform provides a dedicated snow crab modeling workflow integrated directly into
[`ParticleTrackingRun.jl --snowcrab-settings`](ParticleTrackingRun.jl) (shorthand: `--snowcrab`):
- **Real 5-Year Physical Cycle**: High-resolution regional circulation driven by real
  NOAA bathymetry/winds, $M_2+S_2$ spring-neap tides, and surface heat flux.
- **2-Year Climatological Average**: Continuous multi-year seasonal stratification
  cycle capturing settlement patterns and cohort carryover.
- **Coupled Larval Biology**: Context-dependent DVM, $T_0 = -1.5^\circ\text{C}$
  degree-day molting, benthic release with vertical ascent, and BBL shear.
- **Decoupled Workflow**: Decouple 3D hydrodynamic integration from particle tracking
  via `--hydro-model=<file>`, `--hydro-only`, `--track-only`, and `--reuse-hydro`.
- **DuckDB Comparison**: Multi-cohort storage, tagging (`--run-id`), and side-by-side
  SQL comparative analytics (`--compare`).
- **Overrides**: Any subsequent CLI flag directly overrides the baseline snow crab settings.
- **Scripting API**: Call `SnowCrabRunOptions(; kwargs...)` for programmatic Julia execution.

```bash
# 1. Run hydrodynamics only and save flow field to a chosen model file:
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --real-5yr --hydro-only --hydro-model=hydrodynamics1.jld2

# 2. Track larval cohort reusing the saved hydrodynamics:
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_spring_2020 --particles=500 --ascent

# 3. Track a second cohort under alternate biological behavior (e.g. faster ascent):
julia --project=. ParticleTrackingRun.jl --snowcrab-settings --track-only --hydro-model=hydrodynamics1.jld2 \
    --run-id=cohort_summer_fast --particles=500 --ascent-speed=0.015

# 4. Compare all cohorts side-by-side in DuckDB:
julia --project=. ParticleTrackingRun.jl --compare

# 5. Quick end-to-end test run:
julia --project=. ParticleTrackingRun.jl --snowcrab --all --quick
```

---

## Centralized Configuration

All physical, biological, and numerical parameters are declared in standardized configuration files (TOML format) spanning 13 sections:
- [`inputs/ParticleTracking.config`](inputs/ParticleTracking.config): Default regional modeling configuration.
- [`inputs/snowcrab.config`](inputs/snowcrab.config): Calibrated snow crab (*Chionoecetes opilio*) baseline configuration.

---

## Documentation & References

- 📖 **[Comprehensive User & Reference Manual](docs/ParticalTracking_user_guide.md)**: Full mathematical derivations, biophysical equations, schema definitions, and developer guide.
- 🚀 **[CLI & Workflow Execution Guide](ParticleTrackingRun.md)**: Complete command-line options reference and workflow recipes.
- 📄 **[Scientific Research Paper](docs/snow_crab_larval_connectivity_paper.md)**: Peer-reviewed paper manuscript describing larval transport mechanisms across the Scotian Shelf.

---

## License

Released under the [MIT License](LICENSE).
