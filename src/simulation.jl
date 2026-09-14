"""
    simulation.jl

Simulation setup, adaptive time integration, numerical stability diagnostics,
and output writing for regional hydrodynamic modeling.
"""

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using Oceananigans.OutputWriters: JLD2Writer, Checkpointer, checkpoint
import Oceananigans.OutputWriters: cleanup_checkpoints
using JLD2

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
    setup_hydrodynamic_simulation(
        model;
        Δt::Real = 2minutes,
        stop_time::Real = 12hours,
        adaptive_time_step::Bool = false,
        target_cfl::Real = 0.2,
        max_Δt::Real = 5minutes,
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
- `max_Δt::Real`: Maximum allowed time step in seconds.
- `min_Δt::Real`: Minimum allowed time step before raising divergence error.
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
    max_Δt::Real = 5minutes,
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

        out_sched = if output_schedule isa Int
            IterationInterval(output_schedule)
        else
            TimeInterval(output_schedule)
        end

        writer = JLD2Writer(
            model,
            outputs_dict,
            filename = full_output_path,
            schedule = out_sched,
            overwrite_existing = resolved_overwrite
        )
        sim.output_writers[:fields] = writer
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

# Outputs
- `Function`: Callable `(lon, lat, z, t) -> NamedTuple` returning interpolated field values.

# References
- Marshall, J., et al. (1997). *J. Geophys. Res. Oceans*, 102(C3), 5753-5766.
"""
function create_flow_interpolator_from_jld2(
    jld2_filepath::AbstractString;
    variables::Tuple = (:u, :v, :w, :T)
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

        nt = length(sorted_keys)

        # Extract time vector: handle both JLD2.Group and flat Vector representations
        t_vec = if haskey(file, "timeseries/t")
            t_obj = file["timeseries/t"]
            if t_obj isa JLD2.Group
                [Float64(t_obj[k]) for k in sorted_keys]
            else
                collect(Float64, t_obj)
            end
        else
            [something(tryparse(Float64, k), Float64(idx)) for (idx, k) in enumerate(sorted_keys)]
        end

        # Extract grid coordinate vectors from stored metadata or infer from sample
        sample_u = file["timeseries/u/$(first(sorted_keys))"]
        nx, ny, nz = size(sample_u)

        grid_obj = if haskey(file, "grid")
            file["grid"]
        elseif haskey(file, "serialized/grid")
            file["serialized/grid"]
        else
            nothing
        end

        extract_coords(obj, n) = try
            raw = if hasproperty(obj, :parent)
                collect(Float64, obj.parent)
            else
                collect(Float64, obj)
            end
            if length(raw) >= n
                raw[1:n]
            else
                collect(range(first(raw), last(raw), length = n))
            end
        catch
            nothing
        end

        lons_vec = if !isnothing(grid_obj) && hasproperty(grid_obj, :λᶜᵃᵃ)
            something(extract_coords(grid_obj.λᶜᵃᵃ, nx), collect(range(-68.0, -57.0, length = nx)))
        elseif !isnothing(grid_obj) && hasproperty(grid_obj, :xᶜᵃᵃ)
            something(extract_coords(grid_obj.xᶜᵃᵃ, nx), collect(range(-68.0, -57.0, length = nx)))
        else
            collect(range(-68.0, -57.0, length = nx))
        end

        lats_vec = if !isnothing(grid_obj) && hasproperty(grid_obj, :φᵃᶜᵃ)
            something(extract_coords(grid_obj.φᵃᶜᵃ, ny), collect(range(42.0, 47.0, length = ny)))
        elseif !isnothing(grid_obj) && hasproperty(grid_obj, :yᵃᶜᵃ)
            something(extract_coords(grid_obj.yᵃᶜᵃ, ny), collect(range(42.0, 47.0, length = ny)))
        else
            collect(range(42.0, 47.0, length = ny))
        end

        deps_vec = if !isnothing(grid_obj) && hasproperty(grid_obj, :zᵃᵃᶜ)
            something(extract_coords(grid_obj.zᵃᵃᶜ, nz), collect(range(-1000.0, 0.0, length = nz)))
        elseif !isnothing(grid_obj) && hasproperty(grid_obj, :z) &&
               hasproperty(grid_obj.z, :cᵃᵃᶜ)
            something(extract_coords(grid_obj.z.cᵃᵃᶜ, nz), collect(range(-1000.0, 0.0, length = nz)))
        else
            collect(range(-1000.0, 0.0, length = nz))
        end

        # Load all time snapshots for each variable into 4D arrays (nx, ny, nz, nt)
        u_arr = Array{Float64}(undef, nx, ny, nz, nt)
        v_arr = Array{Float64}(undef, nx, ny, nz, nt)
        w_arr = Array{Float64}(undef, nx, ny, nz, nt)
        T_arr = Array{Float64}(undef, nx, ny, nz, nt)

        has_w   = haskey(file, "timeseries/w")
        has_T   = haskey(file, "timeseries/T")
        has_eta = haskey(file, "timeseries/η")

        η_arr = has_eta ? Array{Float64}(undef, nx, ny, nt) : Array{Float64}(undef, 0, 0, 0)

        for (m, key) in enumerate(sorted_keys)
            raw_u = Float64.(file["timeseries/u/$(key)"])
            raw_v = Float64.(file["timeseries/v/$(key)"])
            u_arr[:, :, :, m] = (size(raw_u) == (nx, ny, nz)) ? raw_u : raw_u[1:nx, 1:ny, 1:nz]
            v_arr[:, :, :, m] = (size(raw_v) == (nx, ny, nz)) ? raw_v : raw_v[1:nx, 1:ny, 1:nz]
            if has_w
                raw_w = Float64.(file["timeseries/w/$(key)"])
                if size(raw_w) == (nx, ny, nz)
                    w_arr[:, :, :, m] = raw_w
                elseif size(raw_w, 3) == nz + 1
                    # Average vertical face values to cell centers
                    w_arr[:, :, :, m] = 0.5 .* (raw_w[1:nx, 1:ny, 1:nz] .+ raw_w[1:nx, 1:ny, 2:nz+1])
                else
                    w_arr[:, :, :, m] = raw_w[1:nx, 1:ny, 1:nz]
                end
            else
                w_arr[:, :, :, m] = zeros(nx, ny, nz)
            end
            if has_T
                raw_T = Float64.(file["timeseries/T/$(key)"])
                T_arr[:, :, :, m] = (size(raw_T) == (nx, ny, nz)) ? raw_T : raw_T[1:nx, 1:ny, 1:nz]
            else
                T_arr[:, :, :, m] = fill(4.5, nx, ny, nz)
            end
            if has_eta
                raw_eta = Float64.(file["timeseries/η/$(key)"])
                η_arr[:, :, m] = (size(raw_eta) == (nx, ny)) ? raw_eta : raw_eta[1:nx, 1:ny]
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
