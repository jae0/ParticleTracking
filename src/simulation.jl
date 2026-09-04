"""
    simulation.jl

Simulation setup, adaptive time integration, numerical stability diagnostics,
and output writing for regional hydrodynamic modeling.
"""

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using Oceananigans.OutputWriters: JLD2Writer
using JLD2

"""
    compute_advective_cfl(model, Δt::Real)

Compute the maximum instantaneous horizontal and vertical advective Courant-
Friedrichs-Lewy (CFL) number across active fluid cells:
```math
\\text{CFL} = \\max \\left( \\frac{|u| \\Delta t}{\\Delta x},
\\frac{|v| \\Delta t}{\\Delta y}, \\frac{|w| \\Delta t}{\\Delta z} \\right)
```

# Inputs
- `model`: Configured `HydrostaticFreeSurfaceModel`.
- `Δt::Real`: Current time step in seconds.

# Outputs
- `Float64`: Maximum advective CFL number.

# References
- Courant, R., Friedrichs, K., & Lewy, H. (1928). Über die partiellen
  Differenzengleichungen der mathematischen Physik. *Mathematische Annalen*,
  100(1), 32-74. DOI: 10.1007/BF01448839
- Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A. (2007).
  *Spectral Methods: Evolution to Complex Geometries and Applications to Fluid
  Dynamics*. Springer-Verlag, Berlin.
"""
function compute_advective_cfl(model, Δt::Real)
    u_max = maximum(abs, model.velocities.u)
    v_max = maximum(abs, model.velocities.v)
    w_max = maximum(abs, model.velocities.w)

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

    return Float64(max(cfl_x, cfl_y, cfl_z))
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
        overwrite_existing::Bool = true,
        watchdog::Bool = true
    )

Configure an Oceananigans `Simulation` with progress callbacks, CFL monitoring,
numerical stability watchdogs, and JLD2 output writers for background hydrodynamic fields.

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
- `adaptive_time_step::Bool`: Whether to dynamically scale \$\\Delta t\$ based on CFL.
- `target_cfl::Real`: Target CFL number when `adaptive_time_step` is enabled.
- `max_Δt::Real`: Maximum allowed time step in seconds.
- `min_Δt::Real`: Minimum allowed time step before raising numerical divergence error.
- `progress_schedule::Union{Real, Int}`: Logging frequency (iterations if Int, seconds if Real).
- `enable_output::Bool`: Whether to attach a JLD2 output writer.
- `output_dir::AbstractString`: Directory where outputs will be saved.
- `output_filename::AbstractString`: Output file name.
- `output_schedule::Union{Real, Int}`: Output saving frequency (iterations or seconds).
- `overwrite_existing::Bool`: Whether to overwrite pre-existing output files.
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
    overwrite_existing::Bool = true,
    watchdog::Bool = true,
    divergence_velocity_limit::Real = 20.0
)
    mkpath(output_dir)
    full_output_path = joinpath(output_dir, output_filename)

    sim = Simulation(model, Δt = Δt, stop_time = stop_time)

    # 1. Progress reporting callback
    prog_sched = if progress_schedule isa Int
        IterationInterval(progress_schedule)
    else
        TimeInterval(progress_schedule)
    end

    progress_fn(s) = begin
        u_max = maximum(abs, s.model.velocities.u)
        v_max = maximum(abs, s.model.velocities.v)
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
            u_max = maximum(abs, s.model.velocities.u)
            if isnan(u_max) || isinf(u_max) || u_max > divergence_velocity_limit
                error(
                    "Numerical divergence detected at iteration $(iteration(s)), " *
                    "time $(prettytime(s)). Velocity magnitude u_max = $(u_max) m/s (limit: $(divergence_velocity_limit) m/s)."
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
            overwrite_existing = overwrite_existing
        )
        sim.output_writers[:fields] = writer
    end

    return sim
end

"""
    run_hydrodynamic_simulation!(simulation; verbose::Bool=true)

Execute time integration of the ocean hydrodynamic simulation with error diagnostics.

# Inputs
- `simulation::Simulation`: Oceananigans simulation to execute.
- `verbose::Bool`: Whether to log start, elapsed wall-clock time, and completion.

# Outputs
- `Simulation`: Completed simulation instance.
"""
function run_hydrodynamic_simulation!(simulation; verbose::Bool = true)
    if verbose
        @info "Starting hydrodynamic simulation..."
        @info "Initial Δt: $(simulation.Δt) s, Stop time: $(prettytime(simulation.stop_time))"
    end

    start_wall_time = time()
    try
        run!(simulation)
    catch err
        @error(
            "Simulation failed at iteration $(iteration(simulation)), " *
            "time: $(prettytime(simulation))",
            exception = (err, catch_backtrace())
        )
        rethrow(err)
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

        # Robust chronological ordering: parse numeric timestamps or iterations
        sorted_keys = sort(u_raw_keys, by = k -> begin
            parsed = tryparse(Float64, k)
            !isnothing(parsed) ? parsed : 0.0
        end)

        nt = length(sorted_keys)

        # Extract time vector
        t_vec = if haskey(file, "timeseries/t")
            collect(Float64, file["timeseries/t"])
        else
            [something(tryparse(Float64, k), Float64(idx)) for (idx, k) in enumerate(sorted_keys)]
        end

        # Extract grid coordinate vectors from stored metadata or infer from sample
        sample_u = file["timeseries/u/$(first(sorted_keys))"]
        nx, ny, nz = size(sample_u)

        if haskey(file, "grid")
            g = file["grid"]
            lons_vec = haskey(g, "λᶜᵃᵃ") ? collect(Float64, g["λᶜᵃᵃ"][1:nx]) :
                       collect(range(-68.0, -57.0, length = nx))
            lats_vec = haskey(g, "φᵃᶜᵃ") ? collect(Float64, g["φᵃᶜᵃ"][1:ny]) :
                       collect(range(42.0, 47.0, length = ny))
            deps_vec = haskey(g, "zᵃᵃᶜ") ? collect(Float64, g["zᵃᵃᶜ"][1:nz]) :
                       collect(range(-1000.0, 0.0, length = nz))
        else
            @warn "No 'grid' key in JLD2 file; using uniform coordinate fallback."
            lons_vec = collect(range(-68.0, -57.0, length = nx))
            lats_vec = collect(range(42.0, 47.0, length = ny))
            deps_vec = collect(range(-1000.0, 0.0, length = nz))
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
            u_arr[:, :, :, m] = Float64.(file["timeseries/u/$(key)"])
            v_arr[:, :, :, m] = Float64.(file["timeseries/v/$(key)"])
            w_arr[:, :, :, m] = has_w ? Float64.(file["timeseries/w/$(key)"]) : zeros(nx, ny, nz)
            T_arr[:, :, :, m] = has_T ? Float64.(file["timeseries/T/$(key)"]) : fill(4.5, nx, ny, nz)
            if has_eta
                η_arr[:, :, m] = Float64.(file["timeseries/η/$(key)"])
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
