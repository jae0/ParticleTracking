# ParticleTracking.jl

**A High-Performance Biophysical Ocean Modeling, Individual-Based Larval Transport, and Demographic Population Connectivity Framework**

[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Backend: Oceananigans.jl](https://img.shields.io/badge/Physics-Oceananigans.jl-informational.svg)](https://github.com/CliMA/Oceananigans.jl)
[![Storage: DuckDB](https://img.shields.io/badge/Storage-DuckDB-yellow.svg)](https://duckdb.org)
[![Visualization: CairoMakie + GLMakie](https://img.shields.io/badge/Visualization-Makie-purple.svg)](https://makie.juliaplots.org)

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
  └── CairoMakie 2D publication charts, GLMakie 3D GPU visualizations & interactive HTML5 Leaflet.js dashboard
```

---

## Quickstart

### 1. Installation & Environment Setup
Clone the repository and instantiate Julia dependencies:
```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### 2. Run the Full Test Suite
Verify all test sets:
```bash
julia --project=. test/runtests.jl
```

### 3. Run the Production Pipeline
Execute an end-to-end simulation using the command-line interface:
```bash
# Fast debug mode (coarse grid, 1 hr hydro, 2 day track)
julia --project=. ParticleTrackingRun.jl --all --quick

# Full production run with custom TOML configuration
julia --project=. ParticleTrackingRun.jl --all --config=inputs/snowcrab.toml

# GPU-accelerated run with automatic CPU fallback
julia --project=. ParticleTrackingRun.jl --all --gpu --fallback-cpu
```

### 4. Decoupled Multi-Cohort Batching
Run the heavy 3D hydrodynamics once, then track multiple distinct cohorts without re-solving fluid equations:
```bash
# Step 1: Solve hydrodynamics once and archive flow field
julia --project=. ParticleTrackingRun.jl --data --grid --model --sim --output-dir=outputs/baseline --config=inputs/snowcrab.toml

# Step 2: Track Cohort A (Spring hatch, benthic release with ascent)
julia --project=. ParticleTrackingRun.jl --track --metrics --viz --output-dir=outputs/baseline --config=inputs/snowcrab.toml \
    --particles=500 --release-mode=bottom --ascent --ascent-speed=0.010 --seed=101

# Step 3: Track Cohort B (Summer hatch, alternate ascent speed)
julia --project=. ParticleTrackingRun.jl --track --metrics --viz --output-dir=outputs/baseline --config=inputs/snowcrab.toml \
    --particles=500 --release-mode=bottom --ascent --ascent-speed=0.015 --seed=201
```

---

## DuckDB Analytics & Scenario Management

All runs, trajectory time series, recruitment metrics, and demographic transition matrices are archived in `outputs/particle_tracking.duckdb`:

```bash
# List all archived simulation runs
julia --project=. ParticleTrackingRun.jl --list-runs

# Multi-scenario comparative analytics
julia --project=. ParticleTrackingRun.jl --compare-scenarios

# Bayesian / ensemble model-averaged demographic connectivity (P_ij ± σ)
julia --project=. ParticleTrackingRun.jl --model-average
```

---

## Configuration System

All physical, biological, and numerical parameters are declared in standardized **TOML configuration files** spanning 13 sections:

| Config File | Purpose |
|-------------|---------|
| `inputs/ParticleTracking.toml` | Default regional modeling configuration |
| `inputs/snowcrab.toml` | Calibrated snow crab (*Chionoecetes opilio*) baseline (Scotian Shelf) |
| `inputs/hydrodynamics_climatology_2yr.toml` | 2-year climatological circulation |
| `inputs/hydrodynamics_2020_2yr.toml` | 2020-2022 real hindcast |
| `inputs/hydrodynamics_mhw_2yr.toml` | Marine heat wave scenario |
| `work/snowcrab/snowcrab.toml` | Snow crab run with output to `work/snowcrab/` |

All species-specific parameters (domain bounds, larval biology, thermal ecology, tessellation) are configured **exclusively via TOML files** — no snowcrab-specific CLI flags remain. The core code is fully species-agnostic.

---

## Documentation & References

- 📖 **[CLI & Workflow Execution Guide](ParticleTrackingRun.md)**: Complete command-line options reference and workflow recipes.
- 📄 **Scientific Research Paper** (`docs/snow_crab_larval_connectivity_paper.md`): Peer-reviewed paper manuscript describing larval transport mechanisms across the Scotian Shelf.

---

## Visualization Capabilities

### 2D Publication Figures (CairoMakie)
- Particle trajectory maps & DVM depth profiles
- Settlement density & thermal exposure maps
- Hydrodynamic fields (velocity, T/S, stratification, diffusion, vorticity)
- Vertical cross-sections with bathymetry masking
- Connectivity matrices & recruitment summaries
- Multi-panel dashboards (`plot_multi_panel_dashboard`)
- Voronoi tessellation visualization (`plot_voronoi_tessellation`)
- Particle fate summaries (`plot_particle_fate_summary`)

### 3D GPU-Accelerated Visualizations (GLMakie)
- **Volume rendering** of hydrodynamic fields with isosurface extraction
- **3D trajectory tubes** with stage-based color gradients
- **3D connectivity networks** between management strata
- Bathymetry terrain mesh & sea surface elevation
- Publication-ready 1920×1080 output

### Interactive Dashboards (Leaflet.js)
- Standalone HTML5 maps with trajectory playback
- Time-slider for temporal exploration
- Layer toggles for bathymetry, currents, particles

---

## License

Released under the [MIT License](LICENSE).