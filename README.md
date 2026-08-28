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
  ├── Astronomical M2 tidal body forcing & Simpson-Hunter mixing fronts
  └── CMIP6 climate warming scenarios (Historical, SSP1-2.6, SSP2-4.5, SSP5-8.5, MHW)
                   │
                   ▼
Individual-Based Lagrangian Particle Tracking (Euler-Maruyama SDE)
  ├── Strict marine water placement (depth >= 100 m)
  ├── Stage-specific Diel Vertical Migration (DVM) active swimming
  ├── In situ thermal degree-day ontogenetic molting (Zoea I -> II -> Megalopa)
  └── Cumulative temperature-dependent PLD, thermal mortality & CIL settlement
                   │
                   ▼
Demographic Connectivity & Analytical Engine (DuckDB & Leaflet)
  ├── Administrative CFA polygon boundary classification (Jordan Curve ray-casting)
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
Verify all 15 test sets (471 unit tests):
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

# Export DuckDB tables to Apache Parquet format
julia --project=. ParticleTrackingRun.jl --export-parquet
```

---

## Centralized Configuration

All physical, biological, and numerical parameters are declared in [`inputs/ParticalTracking.config`](inputs/ParticalTracking.config) (TOML format) spanning 13 sections (`[domain]`, `[grid]`, `[data]`, `[tides]`, `[climate]`, `[hydrodynamics]`, `[biology]`, `[dvm]`, `[molting_and_settlement]`, `[storage]`, `[hardware]`, `[visualization]`, `[paths]`).

---

## Documentation & References

- 📖 **[Comprehensive User & Reference Manual](docs/ParticalTracking_manual.md)**: Full mathematical derivations, biophysical equations, schema definitions, and developer guide.
- 🚀 **[CLI & Workflow Execution Guide](ParticleTrackingRun.md)**: Complete command-line options reference and workflow recipes.
- 📄 **[Scientific Research Paper](snow_crab_larval_connectivity_paper.md)**: Peer-reviewed paper manuscript describing larval transport mechanisms across the Scotian Shelf.

---

## License

Released under the [MIT License](LICENSE).
