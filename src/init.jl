"""
    init.jl

Unified initialization script for the particle tracking and hydrodynamic modeling
environment. Loads all required Julia dependencies and includes all functional
sub-components directly into the caller's workspace.
"""

# 1. Load external package dependencies
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using CairoMakie
using NCDatasets
using Downloads
using JLD2
using Random
using DuckDB
using DataFrames
using DBInterface
using Dates
using Statistics
using LinearAlgebra
using TOML

# 2. Include all functional modules directly into the active workspace
# Determine base directory of this script to allow robust relative inclusions
const _SRC_DIR = @__DIR__

include(joinpath(_SRC_DIR, "configuration.jl"))
include(joinpath(_SRC_DIR, "open_data.jl"))
include(joinpath(_SRC_DIR, "synthetic_data.jl"))
include(joinpath(_SRC_DIR, "architecture.jl"))
include(joinpath(_SRC_DIR, "grid_bathymetry.jl"))
include(joinpath(_SRC_DIR, "hydrodynamic_model.jl"))
include(joinpath(_SRC_DIR, "tides.jl"))
include(joinpath(_SRC_DIR, "climate_scenarios.jl"))
include(joinpath(_SRC_DIR, "simulation.jl"))
include(joinpath(_SRC_DIR, "larval_behavior.jl"))
include(joinpath(_SRC_DIR, "empirical_analysis.jl"))
include(joinpath(_SRC_DIR, "storage_duckdb.jl"))
include(joinpath(_SRC_DIR, "visualization.jl"))

println("Particle tracking and hydrodynamic modeling environment successfully initialized.")
