"""
    simulation.jl

Simulation setup, adaptive time integration, numerical stability diagnostics,
and output writing for regional hydrodynamic modeling.
"""

using Oceananigans
using Oceananigans.AbstractOperations: ∂z
using Oceananigans.Fields: indices
using SeawaterPolynomials
using SeawaterPolynomials: TEOS10EquationOfState
import SeawaterPolynomials.ρ
using SpecialFunctions: erf
using Oceanostics: KineticEnergyDissipationRate
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using Oceananigans.OutputWriters: JLD2Writer, Checkpointer, checkpoint
import Oceananigans.OutputWriters: cleanup_checkpoints
using JLD2

# Import domain constants for fallback coordinate generation

"""
    compute_advective_cfl(
        model,
        Δt::Real;
        bottom_drag::Real = 1e-4,
        cd_drag::Real = 2.5e-3
    ) -> Float64

Estimate the composite Courant-Friedrichs-Lewy (CFL) stability parameter
across the computational domain, evaluating advective velocity transport
in interior cells and the explicit bottom friction / drag Courant limit.

# Mathematical Formulation
```math
\\text{CFL}_{\\text{adv}} = \\max \\left(
    \\frac{|u| \\Delta t}{\\Delta x},
    \\frac{|v| \\Delta t}{\\Delta y},
    \\frac{|w| \\Delta t}{\\Delta z}
\\right)
```
```math
\\text{CFL}_{\\text{drag}} = (r_{\\text{drag}} + C_d |\\boldsymbol{u}_h|) \\Delta t
```
The composite stability index returned is:
```math
\\text{CFL} = \\max(\\text{CFL}_{\\text{adv}}, 0.5 \\, \\text{CFL}_{\\text{drag}})
```
ensuring that adaptive time-step controllers (e.g. `TimeStepWizard`) automatically
throttle \$\\Delta t\$ during intense spring-neap tidal acceleration episodes.

# References
- Canuto, C., et al. (2007). *Spectral Methods*. Springer-Verlag.
"""
function compute_advective_cfl(
    model,
    Δt::Real;
    bottom_drag::Real = 1e-4,
    cd_drag::Real = 2.5e-3
)
    u_max = maximum(abs, interior(model.velocities.u))
    v_max = maximum(abs, interior(model.velocities.v))
    w_max = maximum(abs, interior(model.velocities.w))

    # Estimate horizontal grid spacing (in meters)
    base_g = model.grid isa ImmersedBoundaryGrid ? model.grid.underlying_grid : model.grid
    dx_approx = if hasproperty(base_g, :radius)
        # Spherical grid: compute minimum zonal spacing at highest absolute latitude
        r_earth = Float64(base_g.radius)
        dlon_deg = minimum(diff(collect(base_g.λᶠᵃᵃ[1:base_g.Nx + 1])))
        lat_max_abs = max(abs(base_g.φᵃᶠᵃ[1]), abs(base_g.φᵃᶠᵃ[base_g.Ny + 1]))
        r_earth * cosd(lat_max_abs) * deg2rad(dlon_deg)
    elseif hasproperty(base_g, :Lx)
        Float64(base_g.Lx) / base_g.Nx
    else
        1000.0
    end
    dy_approx = if hasproperty(base_g, :radius)
        r_earth = Float64(base_g.radius)
        dlat_deg = minimum(diff(collect(base_g.φᵃᶠᵃ[1:base_g.Ny + 1])))
        r_earth * deg2rad(dlat_deg)
    elseif hasproperty(base_g, :Ly)
        Float64(base_g.Ly) / base_g.Ny
    else
        1000.0
    end
    dz_approx = Float64(base_g.Lz) / base_g.Nz

    cfl_x = (u_max * Δt) / max(1.0, dx_approx)
    cfl_y = (v_max * Δt) / max(1.0, dy_approx)
    cfl_z = (w_max * Δt) / max(1.0, dz_approx)
    cfl_adv = max(cfl_x, cfl_y, cfl_z)

    # Frictional / quadratic drag Courant number
    speed_h = sqrt(u_max^2 + v_max^2)
    cfl_drag = (bottom_drag + cd_drag * speed_h) * Δt

    return Float64(max(cfl_adv, 0.5 * cfl_drag))
end

"""
    find_latest_checkpoint(
        checkpoint_dir::AbstractString;
        prefix::AbstractString = "checkpoint_\$(resolve_config_name())"
    ) -> Union{Nothing, String}

Search `checkpoint_dir` for serialized Oceananigans model checkpoints matching
`\$(prefix)_iteration*.jld2` and return the path with the highest iteration index.
Falls back to the most recently modified matching checkpoint file if iterations
cannot be parsed from filenames. Returns `nothing` if no checkpoints exist.

# Inputs
- `checkpoint_dir::AbstractString`: Directory containing checkpoint archives.
- `prefix::AbstractString`: Checkpoint filename prefix (default `"checkpoint_\$(resolve_config_name())"`).

# Outputs
- `Union{Nothing, String}`: Path to the latest checkpoint, or `nothing` if none found.
"""
function find_latest_checkpoint(
    checkpoint_dir::AbstractString;
    prefix::AbstractString = "checkpoint_$(resolve_config_name())"
)::Union{Nothing, String}
    if !isdir(checkpoint_dir)
        return nothing
    end
    resolved_prefix = if isempty(strip(prefix)) || prefix == "checkpoint"
        "checkpoint_$(resolve_config_name())"
    else
        String(prefix)
    end
    files = readdir(checkpoint_dir)
    esc_pfx = replace(resolved_prefix, "." => "\\.")
    pattern = Regex("^" * esc_pfx * raw"_iteration(\d+)\.jld2$")
    matches = Tuple{Int, String}[]
    for f in files
        m = match(pattern, f)
        if !isnothing(m)
            iter = tryparse(Int, m.captures[1])
            if !isnothing(iter)
                push!(matches, (iter, joinpath(checkpoint_dir, f)))
            end
        end
    end
    if isempty(matches)
        fallback = filter(f -> startswith(f, resolved_prefix) && endswith(f, ".jld2"), files)
        if isempty(fallback)
            return nothing
        end
        full_paths = joinpath.(checkpoint_dir, fallback)
        return full_paths[argmax(mtime.(full_paths))]
    end
    sort!(matches, by = first)
    return last(matches)[2]
end

"""
    inspect_hydrodynamic_checkpoint(filepath::AbstractString) -> NamedTuple

Inspect an Oceananigans model checkpoint archive on disk, extracting its discrete
iteration number, simulation clock time, serialized prognostic fields, and validity.

# Inputs
- `filepath::AbstractString`: Path to the `.jld2` checkpoint file.

# Outputs
- `NamedTuple`:
  - `exists::Bool`: Whether the target file exists on disk.
  - `valid::Bool`: Whether the file is a readable JLD2 archive with valid metadata.
  - `iteration::Int`: Discrete numerical iteration number recorded at checkpoint.
  - `time::Float64`: Elapsed simulation time in seconds at checkpoint.
  - `fields::Vector{Symbol}`: Prognostic variable symbols found in the archive.
  - `filepath::String`: Normalized file path.
"""
function inspect_hydrodynamic_checkpoint(filepath::AbstractString)
    if !isfile(filepath)
        return (
            exists = false,
            valid = false,
            iteration = 0,
            time = 0.0,
            fields = Symbol[],
            filepath = String(filepath)
        )
    end
    try
        jldopen(filepath, "r") do f
            # Locate the model root group (nested under simulation/model or model or root)
            model_root = if haskey(f, "simulation") && haskey(f["simulation"], "model")
                f["simulation"]["model"]
            elseif haskey(f, "model")
                f["model"]
            else
                f
            end

            iter = 0
            t_val = 0.0

            if haskey(model_root, "clock")
                clk = model_root["clock"]
                if hasproperty(clk, :iteration)
                    iter = Int(clk.iteration)
                end
                if hasproperty(clk, :time)
                    t_val = Float64(clk.time)
                end
            elseif haskey(f, "iteration")
                iter = Int(f["iteration"])
            elseif haskey(f, "clock/iteration")
                iter = Int(f["clock/iteration"])
            end

            if haskey(f, "time") && t_val == 0.0
                t_val = Float64(f["time"])
            elseif haskey(f, "clock/time") && t_val == 0.0
                t_val = Float64(f["clock/time"])
            end

            # Filename fallback for iteration if unparsed
            if iter == 0
                m = match(r"_iteration(\d+)\.jld2$", basename(filepath))
                if !isnothing(m)
                    parsed_iter = tryparse(Int, m.captures[1])
                    if !isnothing(parsed_iter)
                        iter = parsed_iter
                    end
                end
            end

            f_keys = Symbol[]
            for group_name in ("velocities", "tracers")
                if haskey(model_root, group_name)
                    for k in keys(model_root[group_name])
                        push!(f_keys, Symbol(k))
                    end
                elseif haskey(f, group_name)
                    for k in keys(f[group_name])
                        push!(f_keys, Symbol(k))
                    end
                end
            end
            unique!(f_keys)

            return (
                exists = true,
                valid = true,
                iteration = iter,
                time = t_val,
                fields = f_keys,
                filepath = String(filepath)
            )
        end
    catch err
        @warn "Corrupted or unreadable checkpoint archive at $(filepath): $(err)"
        return (
            exists = true,
            valid = false,
            iteration = 0,
            time = 0.0,
            fields = Symbol[],
            filepath = String(filepath)
        )
    end
end

"""
    verify_checkpoint_compatibility(
        filepath::AbstractString,
        expected_u_dim::Tuple;
        divergence_velocity_limit::Real = 20.0
    ) -> Bool

Inspect a JLD2 checkpoint archive for grid compatibility and numerical stability
prior to resuming simulation execution. Validates spatial array dimensions and
bounds checks the physical interior domain velocities (excluding boundary halos).

# Inputs
- `filepath::AbstractString`: Path to the JLD2 checkpoint archive on disk.
- `expected_u_dim::Tuple`: Expected tuple dimensions of the `u` velocity field.
- `divergence_velocity_limit::Real`: Velocity threshold (m/s) indicating numerical instability.

# Outputs
- `Bool`: `true` if checkpoint exists, dimensions match, and interior velocity is stable.
"""
function verify_checkpoint_compatibility(
    filepath::AbstractString,
    expected_u_dim::Tuple;
    divergence_velocity_limit::Real = 20.0
)
    !isfile(filepath) && return false
    is_compat = true
    try
        jldopen(filepath, "r") do f
            if haskey(f, "simulation") && haskey(f["simulation"], "model") &&
               haskey(f["simulation"]["model"], "velocities") &&
               haskey(f["simulation"]["model"]["velocities"], "u") &&
               haskey(f["simulation"]["model"]["velocities"]["u"], "data")

                u_data = f["simulation"]["model"]["velocities"]["u"]["data"]
                cp_dim = size(u_data)
                if cp_dim != expected_u_dim
                    @warn "Checkpoint u dimensions $(cp_dim) differ from model $(expected_u_dim). " *
                          "Ignoring incompatible checkpoint: $(filepath)"
                    is_compat = false
                    return
                end

                # Exclude boundary halos (standard 3 cells on each boundary)
                u_int = (ndims(u_data) == 3 && all(size(u_data) .> 6)) ?
                    @view(u_data[4:end-3, 4:end-3, 4:end-3]) : u_data
                max_u = maximum(abs, u_int)

                max_v = 0.0
                if haskey(f["simulation"]["model"]["velocities"], "v") &&
                   haskey(f["simulation"]["model"]["velocities"]["v"], "data")
                    v_data = f["simulation"]["model"]["velocities"]["v"]["data"]
                    v_int = (ndims(v_data) == 3 && all(size(v_data) .> 6)) ?
                        @view(v_data[4:end-3, 4:end-3, 4:end-3]) : v_data
                    max_v = maximum(abs, v_int)
                end

                spd_max = max(max_u, max_v)
                if isnan(spd_max) || isinf(spd_max) || spd_max > divergence_velocity_limit
                    @warn "Checkpoint $(filepath) contains divergent interior velocities " *
                          "(max|u,v| = $(round(spd_max, digits=2)) m/s > limit $(divergence_velocity_limit) m/s). " *
                          "Ignoring corrupt checkpoint."
                    is_compat = false
                end
            else
                is_compat = false
            end
        end
    catch err
        @warn "Failed to inspect checkpoint $(filepath): $(err)"
        is_compat = false
    end
    return is_compat
end

# Windows-safe checkpoint cleanup:
# On Windows, recently accessed or written JLD2 archives may encounter sharing violations (EBUSY).
# Rather than allowing an unhandled IOError to abort a multi-day simulation, we safely retry
# after garbage collection and gracefully defer removal to subsequent checkpoint intervals.
function Oceananigans.OutputWriters.cleanup_checkpoints(checkpointer::Checkpointer)
    prefix = Oceananigans.OutputWriters.checkpoint_superprefix(checkpointer.prefix)
    filepaths = Oceananigans.OutputWriters.glob(prefix * "*.jld2", checkpointer.dir)
    latest_checkpoint_filepath = Oceananigans.OutputWriters.latest_checkpoint(checkpointer, filepaths)
    for filepath in filepaths
        if filepath != latest_checkpoint_filepath
            try
                rm(filepath; force = true)
            catch err
                if err isa IOError
                    GC.gc()
                    try
                        rm(filepath; force = true)
                    catch retry_err
                        @warn "Unable to prune intermediate checkpoint $(filepath) (resource locked): $(retry_err)"
                    end
                else
                    rethrow(err)
                end
            end
        end
    end
    return nothing
end

"""
    inspect_hydrodynamic_file(
        filepath::AbstractString;
        expected_stop_time::Union{Nothing, Real} = nothing
    ) -> NamedTuple

Inspect an Oceananigans JLD2 diagnostic timeseries archive. Verifies structural
integrity, extracts recorded timestamps, and evaluates whether the simulation
reached the planned integration stop time.

# Mathematical Formulation
A simulation archive is marked complete when the latest recorded timestamp
`t_last` satisfies:
```math
t_{\\text{last}} \\ge t_{\\text{stop}} - \\epsilon
```
where `\\epsilon = 10^{-3}` seconds accounts for discrete floating point
time-stepping rounding.

# Inputs
- `filepath::AbstractString`: Path to the `.jld2` timeseries output file.
- `expected_stop_time::Union{Nothing, Real}`: Optional target stop duration in seconds.

# Outputs
- `NamedTuple`:
  - `exists::Bool`: Whether the output file exists on disk.
  - `valid::Bool`: Whether the archive contains readable timeseries groups.
  - `is_complete::Bool`: Whether the archive reached `expected_stop_time`.
  - `last_time::Float64`: Most advanced simulated timestamp in seconds.
  - `n_timesteps::Int`: Total count of recorded diagnostic time levels.
  - `times::Vector{Float64}`: Chronologically ordered timestamps.
  - `filepath::String`: Normalized file path.
"""
function inspect_hydrodynamic_file(
    filepath::AbstractString;
    expected_stop_time::Union{Nothing, Real} = nothing,
    stop_time::Union{Nothing, Real} = nothing
)
    tgt_stop = isnothing(expected_stop_time) ? stop_time : expected_stop_time

    if !isfile(filepath)
        return (
            exists = false,
            valid = false,
            is_complete = false,
            last_time = 0.0,
            n_timesteps = 0,
            n_records = 0,
            times = Float64[],
            filepath = String(filepath)
        )
    end

    try
        jldopen(filepath, "r") do f
            if !haskey(f, "timeseries/t") || !haskey(f, "timeseries/u")
                return (
                    exists = true,
                    valid = false,
                    is_complete = false,
                    last_time = 0.0,
                    n_timesteps = 0,
                    n_records = 0,
                    times = Float64[],
                    filepath = String(filepath)
                )
            end

            t_obj = f["timeseries/t"]
            t_vals = if t_obj isa JLD2.Group
                num_keys = filter(k -> !isnothing(tryparse(Float64, k)), collect(keys(t_obj)))
                sorted_k = sort(num_keys, by = k -> parse(Float64, k))
                [Float64(t_obj[k]) for k in sorted_k]
            else
                collect(Float64, t_obj)
            end

            last_t = isempty(t_vals) ? 0.0 : last(t_vals)
            complete = if isnothing(tgt_stop)
                !isempty(t_vals)
            else
                last_t >= (Float64(tgt_stop) - 1e-3)
            end

            return (
                exists = true,
                valid = true,
                is_complete = complete,
                last_time = last_t,
                n_timesteps = length(t_vals),
                n_records = length(t_vals),
                times = t_vals,
                filepath = String(filepath)
            )
        end
    catch err
        @warn "Failed to inspect hydrodynamic timeseries file at $(filepath): $(err)"
        return (
            exists = true,
            valid = false,
            is_complete = false,
            last_time = 0.0,
            n_timesteps = 0,
            n_records = 0,
            times = Float64[],
            filepath = String(filepath)
        )
    end
end

"""
    fill_field_const(grid, value) -> Field

A centre-point `Field` on `grid` filled with the constant `value` (m² s⁻¹), halo filled.

Used to publish constant closure diffusivity/viscosity as an ordinary gridded output so
consumers can treat it exactly like any other diagnostic field.
"""
function fill_field_const(grid, value::Real)
    f = Field{Center, Center, Center}(grid)
    f .= value
    return f
end

"""
    closure_diffusivity_constants(model) -> Tuple{Union{Nothing, Real}, Union{Nothing, Real}}

Extract the scalar vertical diffusivity `κ` and viscosity `ν` actually used by the model's
closure, if they are constants.

Handles `ScalarDiffusivity` (scalar `κ`/`ν`); returns `(nothing, nothing)` for closures whose
mixing is a prognostic `Field`, which is handled by [`native_mixing_diagnostics`](@ref).
"""
function closure_diffusivity_constants(model)
    hasproperty(model, :closure) || return (nothing, nothing)

    κ, ν = nothing, nothing
    # `model.closure` is a Tuple when several closures are supplied, but a single closure
    # object when the model normalises one (which is what `ocean_simulation` produces).
    closures = model.closure isa Tuple ? model.closure : (model.closure,)
    for c in closures
        (c isa Tuple || c isa NamedTuple) && continue
        hasproperty(c, :κ) || continue
        hasproperty(c, :ν) || continue
        k_raw = getproperty(c, :κ)
        n_raw = getproperty(c, :ν)
        κ_cand = scalar_or_constant(k_raw)
        ν_cand = scalar_or_constant(n_raw)
        κ_cand === nothing && continue
        ν_cand === nothing && continue
        κ = κ_cand
        ν = ν_cand
    end
    return (κ, ν)
end

"""
    scalar_or_constant(x) -> Union{Nothing, Real}

Return `x` as a `Real` when it is a constant diffusivity/viscosity, or when it is a
`NamedTuple` (as `ClimateDataWrangler`/`ScalarDiffusivity` closures present mixing
coefficients) whose entries are all constants. Returns `nothing` when the coefficient is a
grid `Field`, i.e. genuinely spatially varying, since that case is published as a field
rather than a constant.
"""
function scalar_or_constant(x)
    x isa Real && return Float64(x)
    if x isa NamedTuple
        vals = collect(values(x))
        isempty(vals) && return nothing
        all(v -> v isa Real, vals) || return nothing
        return Float64(first(vals))
    end
    return nothing
end

"""
    closure_field_smagarinsky_ν(model) -> Union{Nothing, Any}

Return the Smagorinsky grid eddy-viscosity field `νₑ` from `model.closure_fields` if the
closure carries one, otherwise `nothing`.
"""
function closure_field_smagarinsky_ν(model)
    hasproperty(model, :closure_fields) || return nothing
    cfs = model.closure_fields isa Tuple ? model.closure_fields : (model.closure_fields,)
    for cf in cfs
        isnothing(cf) && continue
        cf isa NamedTuple || continue
        haskey(cf, :νₑ) || continue
        return cf.νₑ
    end
    return nothing
end

"""
    native_mixing_diagnostics(model) -> Dict{Symbol, Any}

Model-native vertical mixing fields (`κ`, `ν`) promoted to `Field`s on the model's grid.

The hydrodynamics uses `ScalarDiffusivity` closures, so the vertical eddy diffusivity `κ`
and eddy viscosity `ν` that actually enter the tracer/velocity equations are *constants*.
Reporting those constants is the only way for downstream figures to reflect the mixing the
model really applied; re-deriving a spatially varying Richardson-number proxy in the
visualization layer produces values that contradict the integrated physics.

If the closure exposes a prognostic grid `Field` (e.g. `NEMOTKE`, or the Smagorinsky `νₑ`
in `closure_fields`), that `Field` is preferred, since it is then the genuinely spatially
varying quantity the model used.

# Inputs
- `model`: Oceananigans `AbstractModel` (e.g. `HydrostaticFreeSurfaceModel`).

# Outputs
- `Dict{Symbol, Any}`: `:κ` and/or `:ν` mapped to grid `Field`s, restricted to whatever the
  closure can supply. Empty when the closure exposes nothing usable.
"""
function native_mixing_diagnostics(model)
    out = Dict{Symbol, Any}()
    grid = model.grid

    κ_val, ν_val = closure_diffusivity_constants(model)
    isnothing(κ_val) || (out[:κ] = fill_field_const(grid, κ_val))
    isnothing(ν_val) || (out[:ν] = fill_field_const(grid, ν_val))

    # A Smagorinsky closure carries a prognostic eddy viscosity `νₑ` in `closure_fields`.
    # It is published under its own key rather than overwriting the constant background `ν`,
    # so downstream consumers can see both the prescribed and the turbulent contribution
    # instead of having one silently replace the other.
    νₑ = closure_field_smagarinsky_ν(model)
    if !isnothing(νₑ) && occursin("Field", string(typeof(νₑ).name.wrapper))
        out[:νₑ] = νₑ
    end

    return out
end

# Reference seawater density and gravitational acceleration used to convert the
# potential-density gradient into buoyancy frequency: N² = ∂z b.
# (Buoyancy is defined as b = g * (ρ₀ - ρ) / ρ₀, so N² = ∂z b = -(g/ρ₀) ∂ρ/∂z.)
const g_ref = 9.80665  # m s⁻²

"""
    model_equation_of_state(model) -> Union{Nothing, BoussinesqEquationOfState}

Recover the Boussinesq equation of state the model is actually integrating with.

`HydrostaticFreeSurfaceModel` stores it at `model.buoyancy.formulation.equation_of_state`, where
the formulation is a `SeawaterBuoyancy`. Every model field that needs a density — `:ρ` and the
`:N2` derived from it — is therefore evaluated with exactly the EOS the physics used, rather
than a hard-coded polynomial that could drift away from it.

Returns `nothing` when the model carries no such formulation (e.g. a `ConstantBuoyancy`
diagnostic model), so callers can skip density diagnostics instead of guessing.
"""
function model_equation_of_state(model)
    hasproperty(model, :buoyancy) || return nothing
    buoyancy = model.buoyancy
    hasproperty(buoyancy, :formulation) || return nothing
    formulation = buoyancy.formulation
    hasproperty(formulation, :equation_of_state) || return nothing
    return formulation.equation_of_state
end

"""
    model_gravitational_acceleration(model) -> Float64

Gravitational acceleration the model is integrating with, taken from its buoyancy formulation
so the published N² is consistent with the dynamics. Falls back to [`g_ref`](@ref) only when
the model does not advertise one.
"""
function model_gravitational_acceleration(model)
    candidates = if hasproperty(model, :buoyancy)
        buoyancy = model.buoyancy
        (hasproperty(buoyancy, :gravitational_acceleration) ? buoyancy.gravitational_acceleration : nothing,
         hasproperty(buoyancy, :formulation) && hasproperty(buoyancy.formulation, :gravitational_acceleration) ?
             buoyancy.formulation.gravitational_acceleration : nothing,
         hasproperty(buoyancy, :formulation) && hasproperty(buoyancy.formulation, :g) ?
             buoyancy.formulation.g : nothing)
    else
        (nothing, nothing, nothing)
    end

    for candidate in candidates
        candidate isa Real && return Float64(candidate)
    end

    return g_ref
end

"""
    reference_density(model) -> Float64

Reference density ρ₀ used to scale the potential-density gradient into buoyancy frequency.
Prefers the equation of state's own `reference_density`; falls back to `1025.0` kg m⁻³.
"""
function reference_density(model)
    eos = model_equation_of_state(model)
    if !isnothing(eos) && hasproperty(eos, :reference_density)
        ρ₀ = eos.reference_density
        ρ₀ isa Real && return Float64(ρ₀)
    end
    return 1025.0
end

"""
    fill_seawater_density!(ρ_field, model)

Fill `ρ_field` with the in-situ seawater density (kg m⁻³) implied by the model's own `T`/`S`
tracers and its equation of state, evaluated on `ρ_field`'s own grid.

The polynomial is `SeawaterPolynomials.ρ(Θ, Sᴬ, Z, eos)` — the real implementation behind
NumericalEarth's `total_density` alias, which is declared without methods in this build and
therefore cannot be called. Note the argument order: conservative temperature `Θ` first, absolute
salinity `Sᴬ` second. Swapping them yields a plausible-looking but physically wrong density
(≈ 1004 kg m⁻³ for shelf conditions that should read ≈ 1027 kg m⁻³), so the order is fixed here
rather than inferred at each call site. Depth `Z` is taken from the grid's own vertical coordinate
nodes, consistent with the convention used by [`write_grid_coordinates`](@ref).

Halos are filled alongside the interior because `:N2` differentiates vertically and so reads one
column beyond the interior at the top and bottom. Two details matter here:

* The loop runs over the *data* axes, not `axes(ρ_field, k)`. `axes` on a `Field` reports the
  interior (`1:Nz`), so iterating it leaves the halo columns holding whatever `new_data`
  initialised them to; the vertical difference at the first and last interior cells then reads
  uninitialised memory and reports N² off by ~1000×.
* The depth index is *clamped* into `1:Nz` when sampling the vertical coordinate. `znode` is
  undefined outside the grid, and extrapolating there produced density values orders of magnitude
  off. Clamping makes the halo a constant extension of the nearest interior cell, so the boundary
  gradient collapses to zero instead of being garbage.

The fill is always performed on the **host** and copied to the destination array, including on a
GPU run. Two reasons: the TEOS-10 polynomial carries a coefficient vector, so `eos` is not `isbits`
and cannot be captured inside a GPU kernel; and calling `ρ` elementwise over `T`/`S` scalar-indexes
device arrays, which CUDA disallows outright ("Scalar indexing is disallowed"). Since ρ is only
refreshed on the output schedule, the host round-trip is negligible next to a time step.
"""
function fill_seawater_density!(ρ_field, model)
    eos = model_equation_of_state(model)
    isnothing(eos) && return nothing

    T = model.tracers.T
    S = model.tracers.S

    grid = ρ_field.grid
    base_grid = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    cpu_grid = architecture(base_grid) isa GPU ? on_architecture(CPU(), base_grid) : base_grid
    Nz = base_grid.Nz

    data = ρ_field.data
    size(data) == size(T.data) == size(S.data) || error(
        "Seawater density field $(size(data)) does not align with the model tracers " *
        "(T=$(size(T.data)), S=$(size(S.data))).")

    # Transfer unconditionally: `on_architecture(CPU(), …)` is a no-op for host arrays,
    # and testing the destination with `isa CuArray` is not reliable because a `Field`'s
    # `.data` can be an OffsetArray *wrapping* the device buffer rather than a CuArray.
    T_in = on_architecture(CPU(), T.data)
    S_in = on_architecture(CPU(), S.data)

    # Halo cells of a freshly assembled model hold *uninitialised* memory, which is not
    # merely garbage: reading it and feeding it to the TEOS-10 polynomial can overflow to
    # Inf/NaN and then throw a `DomainError` from a `sqrt` of a huge negative discriminant,
    # or corrupt the interpreter outright. Neither is recoverable, and the halo is not
    # something the diagnostics need — they are read from the interior and the halo is
    # masked by the immersed boundary. So the buffer starts life filled with the equation
    # of state's reference density and only the interior is ever evaluated.
    host = fill(Float64(eos.reference_density), size(data))

    # Interior index range per dimension. `indices(field, k)` reports a Colon, so the
    # halo width is derived from the data axes against the underlying grid's interior
    # size; the halo is symmetric, so half the excess is the halo on each side.
    ax1, ax2, ax3 = axes(data, 1), axes(data, 2), axes(data, 3)
    Nx, Ny = base_grid.Nx, base_grid.Ny
    ix = (first(ax1) + (length(ax1) - Nx) ÷ 2) : (first(ax1) + (length(ax1) - Nx) ÷ 2 + Nx - 1)
    iy = (first(ax2) + (length(ax2) - Ny) ÷ 2) : (first(ax2) + (length(ax2) - Ny) ÷ 2 + Ny - 1)
    iz = (first(ax3) + (length(ax3) - Nz) ÷ 2) : (first(ax3) + (length(ax3) - Nz) ÷ 2 + Nz - 1)

    # A halo column is still at a real depth, so the interior index is clamped into 1:Nz
    # when the vertical coordinate is sampled.
    k_phys(k) = clamp(k - first(iz) + 1, 1, Nz)

    @inbounds for k in iz, j in iy, i in ix
        z = Float64(Oceananigans.Grids.znode(k_phys(k), cpu_grid, Center()))
        Ti = Float64(T_in[i, j, k])
        Si = Float64(S_in[i, j, k])

        # Clamp into the polynomial's validity box; land cells legitimately hold T = S = 0,
        # and any residual uninitialised memory is caught by the isfinite guard.
        value = (isfinite(Ti) && isfinite(Si)) ?
                ρ(clamp(Ti, -2.0, 40.0), clamp(Si, 0.0, 42.0), z, eos) :
                Float64(eos.reference_density)

        host[i, j, k] = isfinite(value) ? value : Float64(eos.reference_density)
    end

    # Copy back in place so the caller's field object is updated, on CPU or device alike.
    copyto!(data, host)

    return ρ_field
end

"""
    native_stratification_diagnostics(model) -> Dict{Symbol, Any}

Model-native stratification and vorticity fields, ready to be written by a JLD2 writer.

* `:ρ` — in-situ seawater density (kg m⁻³) on the centre grid, from the model's own EOS.
* `:N2` — buoyancy frequency squared N² = −(g/ρ₀) ∂ρ/∂z (s⁻²), the scaled `∂z` operation
  on `:ρ`, using the model's own gravitational acceleration and reference density.
* `:ζ` — relative vertical vorticity ζ = ∂v/∂x − ∂u/∂y (s⁻¹), likewise on the centre grid.

These are `Field`/`AbstractOperation` objects that Oceananigans evaluates with the model's own
grid, stencils and stencilling, replacing the post-hoc finite-difference re-derivation
previously performed in the visualization layer.

`:ρ` is state derived from the prognostic tracers rather than a static field, so it must be
refreshed as the simulation advances. Call [`refresh_stratification_diagnostics!`](@ref) on each
tick of the output schedule; `setup_hydrodynamic_simulation` registers that as a callback, which
Oceananigans runs before the output writers on the same tick.

# Inputs
- `model`: Oceananigans `AbstractModel` carrying `T` and `S` tracers and `u`/`v` velocities.

# Outputs
- `Dict{Symbol, Any}`: keys `:ρ`, `:N2`, `:ζ` mapped to fields/operations. `:ρ` and `:N2` are
  omitted when the model carries no `T`/`S` tracers or no Boussinesq equation of state; `:ζ` is
  omitted when the model has no horizontal velocities.
"""
function native_stratification_diagnostics(model)
    out = Dict{Symbol, Any}()

    # Relative vertical vorticity ζ = ∂v/∂x − ∂u/∂y on the centre grid.
    # ∂x/∂y are the abstract operations exported by Oceananigans; the writer calls
    # `compute!` on them at each output time.
    if hasproperty(model, :velocities)
        u = model.velocities.u
        v = model.velocities.v
        out[:ζ] = ∂x(v) - ∂y(u)
    end

    # Seawater density and the stratification it implies. Both need the `T` and `S`
    # tracers plus the model's Boussinesq equation of state; without them there is no
    # physically meaningful density to publish, so the keys are omitted rather than
    # filled with a placeholder.
    #
    # NOTE: Native ρ/N² diagnostics are currently CPU-only because the TEOS-10 equation
    # of state is not isbits and cannot be called from GPU kernels. On GPU we skip these
    # diagnostics entirely (the GPU path has deeper stability issues in any case).
    if hasproperty(model, :tracers) && haskey(model.tracers, :T) && haskey(model.tracers, :S) &&
       !isnothing(model_equation_of_state(model)) &&
       !(architecture(model) isa GPU)
        cpu_model = on_architecture(CPU(), model)
        ρ_field_cpu = CenterField(cpu_model.grid)
        fill_seawater_density!(ρ_field_cpu, cpu_model)
        
        ρ_field = on_architecture(architecture(model), ρ_field_cpu)
        out[:ρ] = ρ_field

        # N² = ∂z b with buoyancy b = g (ρ₀ − ρ) / ρ₀, so N² = −(g/ρ₀) ∂ρ/∂z.
        # `∂z(ρ_field)` alone carries units of kg m⁻⁴, so the constant is applied here to
        # publish genuine buoyancy frequency in s⁻². Keeping it as a scalar-multiply
        # operation (rather than a pre-scaled Field) means the writer re-evaluates the whole
        # expression against the current ρ at every write.
        scale = -(model_gravitational_acceleration(model) / reference_density(model))
        out[:N2] = scale * ∂z(ρ_field)
    end

    return out
end

"""
    refresh_stratification_diagnostics!(out::Dict{Symbol, Any}, model)

Recompute the state-derived entries of [`native_stratification_diagnostics`](@ref) from the
model's current tracers, in place.

Only `:ρ` needs refreshing: `:N2` and `:ζ` are `∂z`/`∂x`/`∂y` operations over it or over the
velocities, and the output writer re-evaluates those with `compute!` at every write. A missing
`:ρ` (model without `T`/`S` or without an equation of state) is a no-op.
"""
function refresh_stratification_diagnostics!(out::Dict{Symbol, Any}, model)
    ρ_field = get(out, :ρ, nothing)
    isnothing(ρ_field) && return out

    # Native ρ/N² diagnostics are CPU-only (TEOS-10 not isbits). On GPU this callback
    # is a no-op since the keys were never added in native_stratification_diagnostics.
    architecture(model) isa GPU && return out

    cpu_model = on_architecture(CPU(), model)
    cpu_ρ = on_architecture(CPU(), ρ_field)
    fill_seawater_density!(cpu_ρ, cpu_model)
    
    copyto!(ρ_field.data, cpu_ρ.data)
    return out
end

"""
    write_grid_coordinates

Write the model's cell-centre geographic coordinates to a sidecar JLD2 file next to the
simulation output.

`write_grid_coordinates(model, output_path) -> String`

The main simulation file carries a `serialized/grid` entry, but deserialising it yields an
opaque `JLD2.ReconstructedStatic` rather than a usable `LatitudeLongitudeGrid`, because the
grid's type parameters cannot be resolved in the reading session. Anything that needs the
domain geometry — most importantly the Lagrangian flow interpolator — therefore reads these
plain vectors instead of attempting type reconstruction.

# Inputs
- `model`: Model whose grid coordinates are persisted.
- `output_path`: Path of the main simulation JLD2 file; the sidecar is derived from it.

# Outputs
- `String`: Path of the sidecar file written.

# Notes
Coordinates are geographic degrees for `lons`/`lats` and metres (negative down) for
`depths`, matching the conventions used throughout the package.
"""
function write_grid_coordinates(model, output_path::AbstractString)
    base, _ = splitext(output_path)
    sidecar = string(base, "_grid.jld2")

    base_g = model.grid isa ImmersedBoundaryGrid ? model.grid.underlying_grid : model.grid
    cpu_g = architecture(base_g) isa GPU ? on_architecture(CPU(), base_g) : base_g

    lons = collect(Float64, cpu_g.λᶜᵃᵃ)[cpu_g.Hx+1 : cpu_g.Hx+cpu_g.Nx]
    lats = collect(Float64, cpu_g.φᵃᶜᵃ)[cpu_g.Hy+1 : cpu_g.Hy+cpu_g.Ny]
    depths = [Float64(Oceananigans.Grids.znode(k, cpu_g, Center())) for k in 1:cpu_g.Nz]
    halo = (Int(cpu_g.Hx), Int(cpu_g.Hy), Int(cpu_g.Hz))
    grid_size = (Int(cpu_g.Nx), Int(cpu_g.Ny), Int(cpu_g.Nz))

    jldopen(sidecar, "w") do file
        file["lons"] = lons
        file["lats"] = lats
        file["depths"] = depths
        file["halo"] = halo
        file["size"] = grid_size
    end

    return sidecar
end



"""
    setup_hydrodynamic_simulation(
        model;
        Δt::Real = 2minutes,
        stop_time::Real = 12hours,
        adaptive_time_step::Bool = false,
        target_cfl::Real = 0.2,
        max_Δt::Real = 12.0,
        min_Δt::Real = 10.0,
        progress_schedule::Union{Real, Int} = 20,
        enable_output::Bool = true,
        output_dir::AbstractString = "outputs",
        output_filename::AbstractString = "nova_scotia_hydrodynamics.jld2",
        output_schedule::Union{Real, Int} = 100,
        overwrite_existing::Union{Nothing, Bool} = nothing,
        enable_checkpoint::Bool = true,
        checkpoint_dir::Union{Nothing, AbstractString} = nothing,
        checkpoint_prefix::AbstractString = "checkpoint_$(resolve_config_name())",
        checkpoint_schedule::Union{Nothing, Real, Int} = nothing,
        cleanup_checkpoints::Bool = true,
        pickup::Union{Bool, Symbol, String, Int} = :auto,
        watchdog::Bool = true,
        divergence_velocity_limit::Real = 20.0
    )

Configure an Oceananigans `Simulation` with progress callbacks, CFL monitoring,
numerical stability watchdogs, periodic prognostic state checkpointing, and JLD2
output writers for background hydrodynamic fields.

# Mathematical Formulation & Time Integration
Advances discrete time \$t^{n+1} = t^n + \\Delta t\$ using fractional step or
multistage schemes subject to the CFL condition:
```math
\\text{CFL} = \\max \\left( \\frac{|u| \\Delta t}{\\Delta x},
\\frac{|v| \\Delta t}{\\Delta y} \\right) \\le C_{\\text{target}}
```

# Inputs
- `model`: Configured `HydrostaticFreeSurfaceModel`.
- `Δt::Real`: Initial numerical time step in seconds (default 2 minutes).
- `stop_time::Real`: Simulation duration in seconds (default 12 hours).
- `adaptive_time_step::Bool`: Whether to dynamically scale `Δt` based on CFL.
- `target_cfl::Real`: Target CFL number when `adaptive_time_step` is enabled.
- `max_Δt::Real`: Maximum allowed time step in seconds (default `12.0` to ensure
  high-order WENO stability across immersed step bathymetry).
- `min_Δt::Real`: Minimum allowed time step in seconds (default `0.1`). If the wizard attempts to scale below this floor to preserve CFL, it implies spatial divergence.
- `progress_schedule::Union{Real, Int}`: Logging frequency (iterations or seconds).
- `enable_output::Bool`: Whether to attach a JLD2 timeseries output writer.
- `output_dir::AbstractString`: Directory where outputs will be saved.
- `output_filename::AbstractString`: Output file name.
- `output_schedule::Union{Real, Int}`: Output saving frequency (iterations or seconds).
- `overwrite_existing::Union{Nothing, Bool}`: Whether to overwrite pre-existing output.
  When `nothing` (default) and `pickup` is active, automatically preserves and appends.
- `enable_checkpoint::Bool`: Whether to attach Oceananigans `Checkpointer`.
- `checkpoint_dir::Union{Nothing, AbstractString}`: Directory for checkpoint archives.
  Defaults to `joinpath(output_dir, "checkpoints")`.
- `checkpoint_prefix::AbstractString`: Checkpoint filename prefix (default `"checkpoint_\$(resolve_config_name())"`).
- `checkpoint_schedule::Union{Nothing, Real, Int}`: Checkpoint frequency (defaults to
  `output_schedule` if unspecified).
- `cleanup_checkpoints::Bool`: Automatically prune older intermediate checkpoints.
- `pickup::Union{Bool, Symbol, String, Int}`: Checkpoint pickup specification.
  `:auto` (default) detects the latest existing checkpoint in `checkpoint_dir`.
- `watchdog::Bool`: Whether to attach a NaN/divergence early detection callback.
- `divergence_velocity_limit::Real`: Maximum permissible fluid velocity threshold before
  raising divergence error (default 20.0 m/s).

# Outputs
- `Simulation`: Ready-to-run Oceananigans simulation instance.

# References
- Ramadhan, A., et al. (2020). Oceananigans.jl: Fast and friendly geophysical
  fluid dynamics on GPUs. *Journal of Open Source Software*, 5(53), 2018.
"""
function setup_hydrodynamic_simulation(
    model;
    Δt::Real = 2minutes,
    stop_time::Real = 12hours,
    adaptive_time_step::Bool = false,
    target_cfl::Real = 0.2,
    max_Δt::Real = 12.0,
    min_Δt::Real = 0.1,
    progress_schedule::Union{Real, Int} = 20,
    enable_output::Bool = true,
    output_dir::AbstractString = "outputs",
    output_filename::AbstractString = "nova_scotia_hydrodynamics.jld2",
    output_schedule::Union{Real, Int} = 100,
    overwrite_existing::Union{Nothing, Bool} = nothing,
    enable_checkpoint::Bool = true,
    checkpoint_dir::Union{Nothing, AbstractString} = nothing,
    checkpoint_prefix::AbstractString = "checkpoint_$(resolve_config_name())",
    checkpoint_schedule::Union{Nothing, Real, Int} = nothing,
    cleanup_checkpoints::Bool = true,
    pickup::Union{Bool, Symbol, String, Int} = :auto,
    watchdog::Bool = true,
    divergence_velocity_limit::Real = 20.0
)
    mkpath(output_dir)
    full_output_path = joinpath(output_dir, output_filename)

    # Resolve checkpoint directory
    cp_dir = if isnothing(checkpoint_dir) || isempty(checkpoint_dir)
        joinpath(output_dir, "checkpoints")
    else
        String(checkpoint_dir)
    end
    if enable_checkpoint
        mkpath(cp_dir)
    end

    resolved_cp_prefix = if isempty(strip(checkpoint_prefix)) || checkpoint_prefix == "checkpoint"
        "checkpoint_$(resolve_config_name())"
    else
        String(checkpoint_prefix)
    end

    latest_cp = find_latest_checkpoint(cp_dir, prefix = resolved_cp_prefix)

    # Determine resolved pickup target
    resolved_pickup = if pickup === :auto
        if !isnothing(latest_cp)
            m_dim = size(model.velocities.u.data)
            if verify_checkpoint_compatibility(
                latest_cp,
                m_dim;
                divergence_velocity_limit = divergence_velocity_limit
            )
                @info "Detected existing hydrodynamic checkpoint for pickup: $(latest_cp)"
                latest_cp
            else
                false
            end
        else
            false
        end
    else
        pickup
    end

    # Determine whether to overwrite or append to existing diagnostic output file
    resolved_overwrite = if !isnothing(overwrite_existing)
        overwrite_existing
    elseif resolved_pickup != false && isfile(full_output_path)
        @info "Resuming simulation: appending to existing timeseries at $(full_output_path)."
        false
    else
        true
    end

    sim = Simulation(model, Δt = Δt, stop_time = stop_time)

    # 1. Progress reporting callback
    prog_sched = if progress_schedule isa Int
        IterationInterval(progress_schedule)
    else
        TimeInterval(progress_schedule)
    end

    progress_fn(s) = begin
        u_max = maximum(abs, interior(s.model.velocities.u))
        v_max = maximum(abs, interior(s.model.velocities.v))
        cfl_val = compute_advective_cfl(s.model, s.Δt)
        @info(
            "Iter: $(iteration(s)) | Time: $(prettytime(s)) | " *
            "Δt: $(round(s.Δt, digits=1))s | max(|u|): $(round(u_max, digits=4)) m/s | " *
            "max(|v|): $(round(v_max, digits=4)) m/s | CFL: $(round(cfl_val, digits=3))"
        )
    end
    sim.callbacks[:progress] = Callback(progress_fn, prog_sched)

    # 2. Adaptive time step management
    if adaptive_time_step
        wizard = TimeStepWizard(
            cfl = target_cfl,
            max_change = 1.1,
            min_change = 0.5,
            max_Δt = max_Δt,
            min_Δt = min_Δt
        )
        sim.callbacks[:wizard] = Callback(wizard, IterationInterval(5))
    end

    # 3. Numerical stability watchdog callback
    if watchdog
        stability_check(s) = begin
            u_max = maximum(abs, interior(s.model.velocities.u))
            v_max = maximum(abs, interior(s.model.velocities.v))
            spd_max = max(u_max, v_max)
            if isnan(spd_max) || isinf(spd_max) || spd_max > divergence_velocity_limit
                error(
                    "Numerical divergence detected at iteration $(iteration(s)), " *
                    "time $(prettytime(s)). Velocity magnitude u_max = $(spd_max) m/s " *
                    "(limit: $(divergence_velocity_limit) m/s)."
                )
            elseif spd_max > 5.0 && (iteration(s) % 50 == 0)
                @warn(
                    "Elevated interior velocity magnitude: max(|u|,|v|) = " *
                    "$(round(spd_max, digits=2)) m/s at iteration $(iteration(s)) ($(prettytime(s)))."
                )
            end
        end
        sim.callbacks[:watchdog] = Callback(stability_check, IterationInterval(10))
    end

    # 4. Field output writing
    if enable_output
        outputs_dict = Dict{Symbol, Any}(
            :u => model.velocities.u,
            :v => model.velocities.v,
            :w => model.velocities.w
        )
        for tracer_name in keys(model.tracers)
            outputs_dict[tracer_name] = model.tracers[tracer_name]
        end
        if hasproperty(model, :free_surface) && hasproperty(model.free_surface, :η)
            outputs_dict[:η] = model.free_surface.η
        end

        # Model-native vertical mixing, written on the tracer (centre) grid so the
        # consumer never has to de-stagger or guess at halos for these fields.
        merge!(outputs_dict, native_mixing_diagnostics(model))

        # Model-native stratification (ρ, N²) and vorticity (ζ), evaluated by
        # Oceananigans on the model's own grid rather than re-derived downstream.
        strat_dict = native_stratification_diagnostics(model)
        merge!(outputs_dict, strat_dict)

        out_sched = if output_schedule isa Int
            IterationInterval(output_schedule)
        else
            TimeInterval(output_schedule)
        end

        # ρ is derived from the prognostic `T`/`S` tracers, so it has to be recomputed as
        # the simulation advances; `∂z(ρ)` is an operation and is recomputed by the writer
        # itself. Oceananigans runs `TimeStepCallsite` callbacks before the output writers on
        # each tick (see `Simulations/run.jl`), so sharing the writer's schedule guarantees ρ
        # is current at the moment it is written rather than one interval behind.
        if haskey(strat_dict, :ρ)
            sim.callbacks[:stratification] = Callback(
                s -> refresh_stratification_diagnostics!(strat_dict, s.model),
                out_sched
            )
        end

        writer = JLD2Writer(
            model,
            outputs_dict,
            filename = full_output_path,
            schedule = out_sched,
            overwrite_existing = resolved_overwrite
        )
        sim.output_writers[:fields] = writer

        # Persist the grid coordinate vectors as plain numeric arrays.
        #
        # The writer also stores a `serialized/grid` entry, but deserialising it yields an
        # opaque `JLD2.ReconstructedStatic` rather than a usable grid (the grid's type
        # parameters cannot be resolved in the reading session), so consumers such as
        # `create_flow_interpolator_from_jld2` cannot read coordinates off it. Writing the
        # vectors explicitly is version-independent and removes the dependency on JLD2 type
        # reconstruction entirely.
        write_grid_coordinates(model, full_output_path)
    end

    # 5. Prognostic state checkpointing
    if enable_checkpoint
        cp_sched_val = (isnothing(checkpoint_schedule) || (checkpoint_schedule isa Real && checkpoint_schedule <= 0)) ?
            output_schedule : checkpoint_schedule
        cp_sched = if cp_sched_val isa Int
            IterationInterval(cp_sched_val)
        else
            TimeInterval(cp_sched_val)
        end

        checkpointer = Checkpointer(
            model,
            schedule = cp_sched,
            dir = cp_dir,
            prefix = resolved_cp_prefix,
            overwrite_existing = true,
            cleanup = cleanup_checkpoints
        )
        sim.output_writers[:checkpointer] = checkpointer
    end

    return sim
end

"""
    run_hydrodynamic_simulation!(
        simulation;
        pickup::Union{Bool, Symbol, String, Int} = :auto,
        checkpoint_at_end::Bool = true,
        verbose::Bool = true
    ) -> Simulation

Execute time integration of the ocean hydrodynamic simulation with error diagnostics,
automated checkpoint pickup, and graceful emergency checkpointing on interruption.

# Mathematical Formulation & Time Integration
Advances prognostic primitive equation variables:
```math
\\boldsymbol{q}^{n+1} = \\boldsymbol{q}^n + \\int_{t^n}^{t^{n+1}} \\mathcal{F}(\\boldsymbol{q}, t) \\, dt
```
When `pickup` is specified or auto-detected, restores prognostic state `q^k`
from discrete checkpoint iteration `k` without recomputing preceding time steps.

# Inputs
- `simulation::Simulation`: Oceananigans simulation to execute.
- `pickup::Union{Bool, Symbol, String, Int}`: Checkpoint pickup specification.
  - `:auto` (default): Picks up from latest checkpoint in `:checkpointer` if present.
  - `true`: Resumes from most recently modified checkpoint in `:checkpointer`.
  - `filepath::String`: Resumes from specified checkpoint archive path.
  - `false`: Integrates from model's current initialized state without checkpoint pickup.
- `checkpoint_at_end::Bool`: Automatically serialize final state upon stopping (default `true`).
- `verbose::Bool`: Whether to log progress, elapsed wall-clock time, and completion.

# Outputs
- `Simulation`: Completed or gracefully paused simulation instance.

# References
- Ramadhan, A., et al. (2020). *JOSS*, 5(53), 2018.
"""
function run_hydrodynamic_simulation!(
    simulation;
    pickup::Union{Bool, Symbol, String, Int} = :auto,
    checkpoint_at_end::Bool = true,
    verbose::Bool = true
)
    # Resolve pickup target
    actual_pickup = if pickup === :auto
        if haskey(simulation.output_writers, :checkpointer)
            cp_writer = simulation.output_writers[:checkpointer]
            latest = find_latest_checkpoint(cp_writer.dir, prefix = cp_writer.prefix)
            if !isnothing(latest)
                m_dim = size(simulation.model.velocities.u.data)
                if verify_checkpoint_compatibility(latest, m_dim)
                    @info "Picking up hydrodynamic simulation from checkpoint: $(latest)"
                    latest
                else
                    false
                end
            else
                false
            end
        else
            false
        end
    else
        pickup
    end

    if verbose
        @info "Starting hydrodynamic simulation..."
        @info "Initial Δt: $(simulation.Δt) s, Stop time: $(prettytime(simulation.stop_time))"
        if actual_pickup != false
            @info "Pickup enabled: $(actual_pickup)"
        end
    end

    start_wall_time = time()
    try
        run!(simulation; pickup = actual_pickup, checkpoint_at_end = checkpoint_at_end)
    catch err
        if err isa InterruptException
            @warn(
                "Hydrodynamic simulation interrupted at iteration $(iteration(simulation)), " *
                "time: $(prettytime(simulation))."
            )
            # Perform emergency state checkpoint
            if haskey(simulation.output_writers, :checkpointer)
                try
                    checkpoint(simulation)
                    cp_w = simulation.output_writers[:checkpointer]
                    @info "Emergency restart checkpoint safely saved to directory: $(cp_w.dir)"
                catch cp_err
                    @error "Failed to save emergency checkpoint during interruption: $(cp_err)"
                end
            end
            @info "Simulation state safely preserved on disk. Resume anytime with --restart."
            return simulation
        else
            @error(
                "Simulation failed at iteration $(iteration(simulation)), " *
                "time: $(prettytime(simulation))",
                exception = (err, catch_backtrace())
            )
            rethrow(err)
        end
    end

    elapsed_time = time() - start_wall_time
    if verbose
        @info "Simulation completed successfully in $(round(elapsed_time, digits=2)) seconds."
    end
    return simulation
end

"""
    create_flow_interpolator_from_jld2(
        jld2_filepath::AbstractString;
        variables::Tuple = (:u, :v, :w, :T)
    )

Construct a high-performance 4D spatiotemporal interpolator function `(lon, lat, z, t)`
from saved Oceananigans JLD2 simulation outputs.

# Mathematical Formulation
Given discrete grid nodes \$(x_i, y_j, z_k)\$ and time levels \$t_m\$, queries at
arbitrary Lagrangian particle positions \$(\\lambda, \\phi, z, t)\$ are evaluated via
trilinear spatial interpolation followed by linear temporal interpolation:
```math
f(\\boldsymbol{x}, t) = (1 - \\theta) f(\\boldsymbol{x}, t_m) + \\theta f(\\boldsymbol{x}, t_{m+1})
```
where \$\\theta = (t - t_m) / (t_{m+1} - t_m)\$.

# Inputs
- `jld2_filepath::AbstractString`: Path to the simulation `.jld2` output file.
- `variables::Tuple`: Tuple of variable symbols to load (`:u, :v, :w, :T`).
- `t_window::Union{Nothing, Tuple{<:Real, <:Real}}`: Optional `(t_start, t_end)` in **seconds**
  restricting which snapshots are loaded. Queries outside the window are clamped to its edges.
  Particle tracking usually covers a small fraction of a long simulation (e.g. a 60-day larval
  window inside a 2-year run), and loading only that window is what keeps the interpolator within
  memory. `nothing` (the default) loads every snapshot.
- `max_snapshots::Union{Nothing, Integer}`: Optional hard cap on the number of snapshots retained,
  applied by evenly subsampling the selected window if the window still exceeds this many. A final
  backstop against running out of memory on a very long or very high-cadence record.
- `time_range_out::Union{Nothing, Vector{Float64}}`: Optional two-element output vector, resized to
  `[t_start, t_end]` (seconds) to report the span actually retained. Use this to size a run from the
  hydrodynamics rather than assuming a duration: queries outside the retained span are clamped to
  the edge snapshot, so tracking past `t_end` would integrate a frozen field.

# Outputs
- `Function`: Callable `(lon, lat, z, t) -> NamedTuple` returning interpolated field values.

# Memory
A production `snowcrab` simulation stores thousands of snapshots on a 345x245x20 grid. Materialising
all of them as `Float64` for `u, v, w, T` needs hundreds of gigabytes, so an unrestricted load
fails with `OutOfMemoryError` and the caller silently substitutes an analytical flow — which is a
scientifically wrong result, not a graceful degradation. Pass `t_window` (and/or `max_snapshots`) so
the resident set matches the part of the record actually being interrogated.

# References
- Marshall, J., et al. (1997). *J. Geophys. Res. Oceans*, 102(C3), 5753-5766.
"""
function create_flow_interpolator_from_jld2(
    jld2_filepath::AbstractString;
    variables::Tuple = (:u, :v, :w, :T),
    domain_lon::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    domain_lat::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    t_window::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    max_snapshots::Union{Nothing, Integer} = nothing,
    time_range_out::Union{Nothing, Vector{Float64}} = nothing
)
    if !isfile(jld2_filepath)
        error("Simulation output file not found at: $(jld2_filepath)")
    end

    # Load all grid coordinates and field arrays into local arrays.
    # The file is closed after this block; the closure captures in-memory arrays.
    local lons_vec, lats_vec, deps_vec, t_vec
    local u_arr, v_arr, w_arr, T_arr
    local has_eta, η_arr

    jldopen(jld2_filepath, "r") do file
        if !haskey(file, "timeseries/u")
            error("JLD2 file at $(jld2_filepath) is missing required 'timeseries/u' group.")
        end

        u_group = file["timeseries/u"]
        u_raw_keys = collect(keys(u_group))
        if isempty(u_raw_keys)
            error("JLD2 timeseries group 'timeseries/u' contains no time snapshots in $(jld2_filepath).")
        end

        # Robust chronological ordering: filter out non-numeric metadata keys (e.g. 'serialized')
        numeric_keys = filter(k -> !isnothing(tryparse(Float64, k)), u_raw_keys)
        if isempty(numeric_keys)
            error("No numeric time snapshot keys found in 'timeseries/u' in $(jld2_filepath).")
        end
        sorted_keys = sort(numeric_keys, by = k -> parse(Float64, k))

        nt_total = length(sorted_keys)

        # Resolve the true simulation times BEFORE selecting snapshots.
        #
        # The snapshot *keys* are not times. They are elapsed-second stamps written by the
        # checkpointer ("0", "50", "100", ... "227100"), while the physical time of each snapshot
        # lives in `timeseries/t` and spans the full simulated period (here 0 -> 6.31e7 s, i.e.
        # 730 days, at ~3.86 h spacing). Selecting on the keys therefore selects on checkpointer
        # indices, not on the model's time axis, and a `t_window` in seconds would select the wrong
        # snapshots -- or none. The interpolator's own temporal axis is `t_vec`, so the window must
        # be applied to `t_vec` and the two must be kept in lockstep.
        t_all = if haskey(file, "timeseries/t")
            t_obj = file["timeseries/t"]
            if t_obj isa JLD2.Group
                [Float64(t_obj[k]) for k in sorted_keys]
            else
                collect(Float64, t_obj)
            end
        else
            # No time group: fall back to the key stamps, which are at least monotonically ordered.
            [parse(Float64, k) for k in sorted_keys]
        end
        length(t_all) == nt_total || (t_all = t_all[1:min(length(t_all), nt_total)])

        # Restrict the retained snapshots BEFORE allocating any 4D field array. This is the only
        # place the selection can happen cheaply: the arrays below are sized from `nt`, and their
        # product with the grid is what exhausts memory on a long or high-cadence run.
        t_lo, t_hi = t_all[1], t_all[end]

        keep = trues(nt_total)
        if !isnothing(t_window)
            lo, hi = Float64(t_window[1]), Float64(t_window[2])
            keep .= (t_all .>= lo) .& (t_all .<= hi)
            # If the window missed every snapshot (e.g. it lies outside the simulated period),
            # fall back to the nearest snapshot rather than producing an empty interpolator.
            if !any(keep)
                keep[argmin(abs.(t_all .- clamp(lo, t_lo, t_hi)))] = true
            end
        end

        selected = findall(keep)
        # `max_snapshots <= 0` (and `nothing`) mean "no cap": load every snapshot in the window.
        # Without this guard a cap of 0 would subsample to zero snapshots and produce an empty
        # interpolator, because `length(selected) > 0` holds for any non-empty selection.
        if !isnothing(max_snapshots) && max_snapshots > 0 && length(selected) > max_snapshots
            # Evenly subsample the selection, always keeping the first and last so the temporal
            # extent of the requested window is preserved.
            idx = unique!(round.(Int, range(1, length(selected), length = max_snapshots)))
            selected = selected[idx]
        end
        isempty(selected) && error(
            "No snapshots retained from $(jld2_filepath): the selection is empty " *
            "(t_window=$(t_window), max_snapshots=$(max_snapshots)).")

        sorted_keys = sorted_keys[selected]
        t_vec = t_all[selected]
        nt = length(sorted_keys)
        if nt < nt_total
            span_days = (t_vec[end] - t_vec[1]) / 86400
            @info "create_flow_interpolator_from_jld2: retained $nt of $nt_total snapshots " *
                  "from $(basename(jld2_filepath)) (t_window=$(t_window), " *
                  "max_snapshots=$(max_snapshots)); covering " *
                  "$(round(t_vec[1], digits=1)) .. $(round(t_vec[end], digits=1)) s " *
                  "($(round(span_days, digits=2)) days)"
        end
        # Report the retained span to the caller, so a run can be sized from the hydrodynamics
        # rather than from an assumed duration.
        if !isnothing(time_range_out)
            resize!(time_range_out, 2)
            time_range_out[1] = t_vec[1]
            time_range_out[2] = t_vec[end]
        end

        # Grid coordinates come from the sidecar written alongside the simulation
        # (`<output>_grid.jld2`), which holds plain numeric vectors. We deliberately do not
        # read `serialized/grid`: deserialising it yields an opaque
        # `JLD2.ReconstructedStatic` rather than a usable grid, because the grid's type
        # parameters cannot be resolved in the reading session.
        sidecar = string(first(splitext(jld2_filepath)), "_grid.jld2")
        coords = if isfile(sidecar)
            JLD2.jldopen(sidecar, "r") do cf
                (lons = collect(Float64, cf["lons"]),
                 lats = collect(Float64, cf["lats"]),
                 depths = collect(Float64, cf["depths"]),
                 halo = Int.(collect(cf["halo"])),
                 size = Int.(collect(cf["size"])))
            end
        else
            nothing
        end

        ug = nothing   # grid type reconstruction is not used; see `coords` above

        sample_u = file["timeseries/u/$(first(sorted_keys))"]
        Nx = (!isnothing(coords)) ? coords.size[1] : size(sample_u, 1)
        Ny = (!isnothing(coords)) ? coords.size[2] : size(sample_u, 2)
        Nz = (!isnothing(coords)) ? coords.size[3] : size(sample_u, 3)

        Hx = (!isnothing(coords)) ? coords.halo[1] : 0
        Hy = (!isnothing(coords)) ? coords.halo[2] : 0
        Hz = (!isnothing(coords)) ? coords.halo[3] : 0

        i_c = (1 + Hx):(Nx + Hx)
        j_c = (1 + Hy):(Ny + Hy)
        k_c = (1 + Hz):(Nz + Hz)

        extract_coords(obj) = try
            if hasproperty(obj, :parent)
                collect(Float64, obj.parent)
            else
                collect(Float64, obj)
            end
        catch
            nothing
        end

        # Horizontal coordinates come from the sidecar when present. Otherwise fall back to a
        # caller-supplied domain, which originates from the TOML `[domain]` section.
        #
        # The fallback is resolved lazily: `domain_lon`/`domain_lat` are only consulted when
        # `coords === nothing`, so a missing argument is an error only when there is genuinely
        # no sidecar to read. Resolving it eagerly made a perfectly valid sidecar appear to
        # fail whenever the caller passed no explicit domain.
        missing_sidecar_hint(key) =
            "Simulation file $(jld2_filepath) has no sidecar grid coordinates " *
            "($(first(splitext(jld2_filepath)))_grid.jld2, written by " *
            "write_grid_coordinates) and no `$(key)` was supplied. Pass the configured " *
            "study-domain range (e.g. opts.$(key) from the TOML `[domain]` section)."

        lon_raw = if isnothing(coords)
            isnothing(domain_lon) && error(missing_sidecar_hint("domain_lon"))
            collect(range(Float64(domain_lon[1]), Float64(domain_lon[2]), length = Nx))
        else
            coords.lons
        end

        lats_vec_raw = if isnothing(coords)
            isnothing(domain_lat) && error(missing_sidecar_hint("domain_lat"))
            collect(range(Float64(domain_lat[1]), Float64(domain_lat[2]), length = Ny))
        else
            coords.lats
        end

        lons_vec = length(lon_raw) >= (Nx + 2 * Hx) ? lon_raw[i_c] :
                   (length(lon_raw) >= Nx ? lon_raw[1:Nx] :
                    collect(range(first(lon_raw), last(lon_raw), length = Nx)))

        lat_raw = lats_vec_raw
        lats_vec = length(lat_raw) >= (Ny + 2 * Hy) ? lat_raw[j_c] :
                   (length(lat_raw) >= Ny ? lat_raw[1:Ny] :
                    collect(range(first(lat_raw), last(lat_raw), length = Ny)))

        dep_raw = if !isnothing(coords)
            coords.depths
        elseif !isnothing(ug) && hasproperty(ug, :zᵃᵃᶜ)
            something(extract_coords(ug.zᵃᵃᶜ), collect(range(-1000.0, 0.0, length = Nz)))
        elseif !isnothing(ug) && hasproperty(ug, :z) && hasproperty(ug.z, :cᵃᵃᶜ)
            something(extract_coords(ug.z.cᵃᵃᶜ), collect(range(-1000.0, 0.0, length = Nz)))
        else
            collect(range(-1000.0, 0.0, length = Nz))
        end
        deps_vec = length(dep_raw) >= (Nz + 2 * Hz) ? dep_raw[k_c] :
                   (length(dep_raw) >= Nz ? dep_raw[1:Nz] :
                    collect(range(first(dep_raw), last(dep_raw), length = Nz)))

        # Load all time snapshots for each variable into 4D arrays (Nx, Ny, Nz, nt)
        u_arr = Array{Float64}(undef, Nx, Ny, Nz, nt)
        v_arr = Array{Float64}(undef, Nx, Ny, Nz, nt)
        w_arr = Array{Float64}(undef, Nx, Ny, Nz, nt)
        T_arr = Array{Float64}(undef, Nx, Ny, Nz, nt)

        has_w   = haskey(file, "timeseries/w")
        has_T   = haskey(file, "timeseries/T")
        has_eta = haskey(file, "timeseries/η")

        η_arr = has_eta ? Array{Float64}(undef, Nx, Ny, nt) : Array{Float64}(undef, 0, 0, 0)

        for (m, key) in enumerate(sorted_keys)
            raw_u = Float64.(file["timeseries/u/$(key)"])
            u_arr[:, :, :, m] = if size(raw_u) == (Nx, Ny, Nz)
                raw_u
            elseif size(raw_u) >= (Nx + 1 + 2 * Hx, Ny + 2 * Hy, Nz + 2 * Hz)
                u_face = raw_u[(1 + Hx):(Nx + 1 + Hx), j_c, k_c]
                0.5 .* (u_face[1:Nx, :, :] .+ u_face[2:Nx + 1, :, :])
            elseif size(raw_u, 1) == Nx + 1
                0.5 .* (raw_u[1:Nx, 1:Ny, 1:Nz] .+ raw_u[2:Nx + 1, 1:Ny, 1:Nz])
            else
                raw_u[1:Nx, 1:Ny, 1:Nz]
            end

            raw_v = Float64.(file["timeseries/v/$(key)"])
            v_arr[:, :, :, m] = if size(raw_v) == (Nx, Ny, Nz)
                raw_v
            elseif size(raw_v) >= (Nx + 2 * Hx, Ny + 1 + 2 * Hy, Nz + 2 * Hz)
                v_face = raw_v[i_c, (1 + Hy):(Ny + 1 + Hy), k_c]
                0.5 .* (v_face[:, 1:Ny, :] .+ v_face[:, 2:Ny + 1, :])
            elseif size(raw_v, 2) == Ny + 1
                0.5 .* (raw_v[1:Nx, 1:Ny, 1:Nz] .+ raw_v[1:Nx, 2:Ny + 1, 1:Nz])
            else
                raw_v[1:Nx, 1:Ny, 1:Nz]
            end

            if has_w
                raw_w = Float64.(file["timeseries/w/$(key)"])
                w_arr[:, :, :, m] = if size(raw_w) == (Nx, Ny, Nz)
                    raw_w
                elseif size(raw_w) >= (Nx + 2 * Hx, Ny + 2 * Hy, Nz + 1 + 2 * Hz)
                    w_face = raw_w[i_c, j_c, (1 + Hz):(Nz + 1 + Hz)]
                    0.5 .* (w_face[:, :, 1:Nz] .+ w_face[:, :, 2:Nz + 1])
                elseif size(raw_w, 3) == Nz + 1
                    0.5 .* (raw_w[1:Nx, 1:Ny, 1:Nz] .+ raw_w[1:Nx, 1:Ny, 2:Nz + 1])
                else
                    raw_w[1:Nx, 1:Ny, 1:Nz]
                end
            else
                w_arr[:, :, :, m] = zeros(Nx, Ny, Nz)
            end

            if has_T
                raw_T = Float64.(file["timeseries/T/$(key)"])
                T_arr[:, :, :, m] = if size(raw_T) == (Nx, Ny, Nz)
                    raw_T
                elseif size(raw_T) >= (Nx + 2 * Hx, Ny + 2 * Hy, Nz + 2 * Hz)
                    raw_T[i_c, j_c, k_c]
                else
                    raw_T[1:Nx, 1:Ny, 1:Nz]
                end
            else
                T_arr[:, :, :, m] = fill(4.5, Nx, Ny, Nz)
            end

            if has_eta
                raw_eta = Float64.(file["timeseries/η/$(key)"])
                η_arr[:, :, m] = if size(raw_eta) == (Nx, Ny)
                    raw_eta
                elseif size(raw_eta) >= (Nx + 2 * Hx, Ny + 2 * Hy)
                    raw_eta[i_c, j_c]
                else
                    raw_eta[1:Nx, 1:Ny]
                end
            end
        end
    end

    nx, ny, nz, nt = size(u_arr)

    function flow_interpolator(lon::Real, lat::Real, z::Real, t::Real)
        # Find spatial indices via binary search with boundary-safe clamping
        i_f = searchsortedlast(lons_vec, Float64(lon))
        j_f = searchsortedlast(lats_vec, Float64(lat))
        k_f = searchsortedlast(deps_vec, Float64(z))
        i = nx > 1 ? clamp(i_f, 1, nx - 1) : 1
        j = ny > 1 ? clamp(j_f, 1, ny - 1) : 1
        k = nz > 1 ? clamp(k_f, 1, nz - 1) : 1

        # Spatial fractional coordinates
        sx = (nx > 1 && (lons_vec[i+1] - lons_vec[i]) != 0.0) ?
             clamp((Float64(lon) - lons_vec[i]) / (lons_vec[i+1] - lons_vec[i]), 0.0, 1.0) : 0.0
        sy = (ny > 1 && (lats_vec[j+1] - lats_vec[j]) != 0.0) ?
             clamp((Float64(lat) - lats_vec[j]) / (lats_vec[j+1] - lats_vec[j]), 0.0, 1.0) : 0.0
        sz = (nz > 1 && (deps_vec[k+1] - deps_vec[k]) != 0.0) ?
             clamp((Float64(z) - deps_vec[k])  / (deps_vec[k+1]  - deps_vec[k]),  0.0, 1.0) : 0.0

        # Find temporal bracket
        m_f = searchsortedlast(t_vec, Float64(t))
        m   = nt > 1 ? clamp(m_f, 1, nt - 1) : 1
        θ   = (nt > 1 && (t_vec[m+1] - t_vec[m]) != 0.0) ?
              clamp((Float64(t) - t_vec[m]) / (t_vec[m+1] - t_vec[m]), 0.0, 1.0) : 0.0

        function sample_at_time(arr, time_idx)
            i_next = nx > 1 ? i + 1 : i
            j_next = ny > 1 ? j + 1 : j
            k_next = nz > 1 ? k + 1 : k

            v000 = arr[i,      j,      k,      time_idx]
            v100 = arr[i_next, j,      k,      time_idx]
            v010 = arr[i,      j_next, k,      time_idx]
            v110 = arr[i_next, j_next, k,      time_idx]
            v001 = arr[i,      j,      k_next, time_idx]
            v101 = arr[i_next, j,      k_next, time_idx]
            v011 = arr[i,      j_next, k_next, time_idx]
            v111 = arr[i_next, j_next, k_next, time_idx]

            return (1.0 - sx) * (1.0 - sy) * (1.0 - sz) * v000 +
                   sx         * (1.0 - sy) * (1.0 - sz) * v100 +
                   (1.0 - sx) * sy         * (1.0 - sz) * v010 +
                   sx         * sy         * (1.0 - sz) * v110 +
                   (1.0 - sx) * (1.0 - sy) * sz         * v001 +
                   sx         * (1.0 - sy) * sz         * v101 +
                   (1.0 - sx) * sy         * sz         * v011 +
                   sx         * sy         * sz         * v111
        end

        function sample_eta_at_time(time_idx)
            i_next = nx > 1 ? i + 1 : i
            j_next = ny > 1 ? j + 1 : j
            e00 = η_arr[i,      j,      time_idx]
            e10 = η_arr[i_next, j,      time_idx]
            e01 = η_arr[i,      j_next, time_idx]
            e11 = η_arr[i_next, j_next, time_idx]
            return (1.0 - sx) * (1.0 - sy) * e00 +
                   sx         * (1.0 - sy) * e10 +
                   (1.0 - sx) * sy         * e01 +
                   sx         * sy         * e11
        end

        function interp(arr)
            val0 = sample_at_time(arr, m)
            if nt > 1 && θ > 0.0
                val1 = sample_at_time(arr, m + 1)
                return (1.0 - θ) * val0 + θ * val1
            else
                return val0
            end
        end

        function interp_eta()
            if !has_eta
                return 0.0
            end
            val0 = sample_eta_at_time(m)
            if nt > 1 && θ > 0.0
                val1 = sample_eta_at_time(m + 1)
                return (1.0 - θ) * val0 + θ * val1
            else
                return val0
            end
        end

        return (
            u = interp(u_arr),
            v = interp(v_arr),
            w = interp(w_arr),
            T = interp(T_arr),
            η = interp_eta()
        )
    end

    return flow_interpolator
end
