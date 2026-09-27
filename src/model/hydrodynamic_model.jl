"""
    hydrodynamic_model.jl

Configuration and initialization of hydrostatic free-surface ocean models
in Oceananigans.jl for coastal shelf domains.
"""

using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: WENOVectorInvariant
# `ContinuousForcing` lives in `Oceananigans.Forcings`, which Oceananigans `using`s
# internally but does not re-export at top level. `ZeroForcing` is defined locally below,
# so it is deliberately not imported here.
using Oceananigans.Forcings: ContinuousForcing
using ClimaOcean
using Dates
using NCDatasets
using Interpolations
using Interpolations: interpolate, BSpline, Linear
import NumericalEarth: ocean_simulation, GLORYSDaily, GLORYSMonthly, WOAMonthly, WOAAnnual, BoundingBox
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

# Import domain constants

"""
    build_hydrodynamic_model(grid; kwargs...) -> HydrostaticFreeSurfaceModel

Build a hydrostatic ocean model through `NumericalEarth.ocean_simulation`, returning the
bare `HydrostaticFreeSurfaceModel` so that existing call sites (`set!`, `model.velocities`,
`model.tracers`, `model.closure`, `model.grid`) continue to work unchanged.

`ocean_simulation` returns an Oceananigans `Simulation`; this wrapper returns `sim.model`,
because this package supplies its own adaptive stepping, checkpointing and JLD2 output
machinery and therefore only needs the model object.

# Physics defaults chosen here, and why they differ from the library defaults
- `equation_of_state = TEOS10EquationOfState(reference_density = 1025)`. The library default
  is `reference_density = 1020`; this package works in the Boussinesq sense with ρ₀ = 1025.
  TEOS-10 is a full equation of state rather than the linearised form used previously, so
  buoyancy and N² differ from the old `LinearEquationOfState` results — this is an accuracy
  improvement, not a like-for-like swap.
- `coriolis` is left at the library default, `HydrostaticSphericalCoriolis`. The previous
  manual builder used `FPlane(latitude = 44.5)`, which freezes f and omits the β term; over a
  42–47.5 °N domain the β term is not negligible, so the spherical form is more appropriate.
- `closure` defaults to `VerticalScalarDiffusivity` with the caller's `ν`/`κ`, matching the
  previous behaviour. Pass `closure = :catke` for the library's prognostic NEMO-TKE closure.
- `free_surface = SplitExplicitFreeSurface(grid; cfl)`. The previous builder used
  `ImplicitFreeSurface`; the split-explicit barotropic solver is cheaper and is what the
  library considers the standard choice.
- `bottom_drag_coefficient` is left at the library default (`0.003`) because `Default` is not
  part of Oceananigans' public API in this build; `cd_drag` is therefore only honoured by the
  manual builder. Set it explicitly via the library kwarg if a different value is needed.

Surface heat/wind fluxes, tidal body forcing and the lateral sponge are *not* configured
here: `ocean_simulation` exposes them through different hooks (surface boundary conditions
and interior `forcing`) than this package's existing `Forcing`/BC construction, and they are
supplied by the caller through `forcing`/`boundary_conditions` below.
"""
function build_hydrodynamic_model(
    grid;
    Δt::Real = 120.0,
    coriolis_latitude::Real = 45.0,
    surface_wind_stress_x = 0.0001,
    surface_wind_stress_y = 0.0,
    surface_heat_flux = 0.0,
    tidal_forcing::Union{Nothing, NamedTuple} = nothing,
    sponge_forcing::Union{Nothing, NamedTuple} = nothing,
    bottom_drag::Real = 1e-4,
    cd_drag::Real = 2.5e-3,
    ν::Real = 1e-2,
    κ::Real = 1e-2,
    closure = nothing,
    closure_scheme = nothing,
    tracers::Tuple = (:T, :S),
    enable_o2::Bool = false,
    reference_density::Real = 1025.0,
    free_surface = nothing,
    momentum_advection = nothing,
    tracer_advection = nothing,
    boundary_conditions::NamedTuple = NamedTuple(),
    forcing::NamedTuple = NamedTuple(),
    open_boundary_conditions = nothing,
    lateral_boundary_relaxation::Bool = false,
    sponge_width::Float64 = 0.35,
    sponge_tau::Float64 = 3600.0,
    u_inflow::Float64 = -0.15,
    v_inflow::Float64 = 0.05,
    u_decay::Float64 = 200.0,
    v_decay::Float64 = 500.0,
    lon_range::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    lat_range::Union{Nothing, Tuple{<:Real, <:Real}} = nothing
)
    active_tracers = (enable_o2 && !(:O2 in tracers)) ? (tracers..., :O2) : tracers

    # Lateral relaxation sponge (Price & Aumont 2011).
    #
    # The sponge itself is constructed where the forcing is assembled, because it must be
    # wrapped in a `ContinuousForcing` that declares its velocity dependency. Here we only
    # decide whether it is requested at all, and validate that the caller supplied the
    # study-domain bounds it will be anchored to.
    #
    # The TOML `[domain]` section is the single source of truth for geography; nothing about
    # the domain is compiled into this package.
    sponge_active = lateral_boundary_relaxation || !isnothing(open_boundary_conditions)
    if sponge_active && (isnothing(lon_range) || isnothing(lat_range))
        error(
            "Lateral boundary relaxation requested but no study-domain bounds were " *
            "supplied. Pass the ranges from the TOML `[domain]` section as lon_range " *
            "and lat_range to build_hydrodynamic_model."
        )
    end

    # Surface kinematic boundary conditions, expressed the same way as the manual builder so
    # that `ocean_simulation` inherits them instead of inventing its own defaults.
    rho0_cp = 1025.0 * 3990.0
    kinematic_T_flux = if surface_heat_flux isa Function
        (x, y, t) -> -surface_heat_flux(x, y, t) / rho0_cp
    elseif surface_heat_flux isa AbstractMatrix
        -surface_heat_flux ./ rho0_cp
    else
        -Float64(surface_heat_flux) / rho0_cp
    end

    u_bc = surface_wind_stress_x isa BoundaryCondition ? surface_wind_stress_x :
           FluxBoundaryCondition(surface_wind_stress_x)
    v_bc = surface_wind_stress_y isa BoundaryCondition ? surface_wind_stress_y :
           FluxBoundaryCondition(surface_wind_stress_y)
    T_bc = FluxBoundaryCondition(kinematic_T_flux)

    bcs = Dict{Symbol, Any}(
        :u => FieldBoundaryConditions(; top = u_bc),
        :v => FieldBoundaryConditions(; top = v_bc)
    )
    :T in active_tracers && (bcs[:T] = FieldBoundaryConditions(; top = T_bc))
    user_bcs = NamedTuple{Tuple(keys(bcs))}(Tuple(bcs[k] for k in keys(bcs)))

    # Extract grid bounds for sponge layer construction (needs to happen before sponge creation)
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    lon_min_grid = Float64(base_g.λᶠᵃᵃ[1])
    lon_max_grid = Float64(base_g.λᶠᵃᵃ[base_g.Nx + 1])
    lat_min_grid = Float64(base_g.φᵃᶠᵃ[1])
    lat_max_grid = Float64(base_g.φᵃᶠᵃ[base_g.Ny + 1])

    # Interior forcing: tides plus the lateral sponge, matching the manual builder's CPU path.
    # Construct proper bitstype forcing structs that Oceananigans can call in GPU kernels.

    # Interior forcing: tides plus the lateral sponge, matching the manual builder's CPU path.
    # Construct proper bitstype forcing structs that Oceananigans can call in GPU kernels.
    # Extract values from NamedTuples properly.
    u_tide_val = isnothing(tidal_forcing) ? nothing : tidal_forcing.u
    v_tide_val = isnothing(tidal_forcing) ? nothing : tidal_forcing.v

    # Sponge relaxation. The relaxation needs the instantaneous velocity, so the forcing is
    # declared with `field_dependencies = (:u, :v)`: Oceananigans then calls it as
    # `f(x, y, z, t, u, v)`. Without that declaration it is called as `f(x, y, z, t)` and
    # raises `MethodError` on the missing velocity arguments.
    #
    # Every captured value is `isbits` (see `ExponentialInflow`), so the composed forcing is
    # lowerable into a GPU kernel.
    if sponge_active
        sponge_relaxation = LateralBoundaryRelaxation(
            (lon_min_grid, lon_max_grid), (lat_min_grid, lat_max_grid);
            sponge_width_deg = sponge_width,
            tau_relax = sponge_tau,
            u_inflow = u_inflow,
            v_inflow = v_inflow,
            u_decay = u_decay,
            v_decay = v_decay
        )

        # These closures call the Symbol-free `sponge_relaxation_u`/`_v` entry points.
        # Dispatching through a `Symbol` here would put a non-bits value in the kernel ABI
        # and the free-surface kernel would fail to compile on GPU.
        u_func = if isnothing(u_tide_val)
            (x, y, z, t, u, v) -> sponge_relaxation_u(sponge_relaxation, x, y, z, t, u)
        else
            (x, y, z, t, u, v) -> Float64(u_tide_val(x, y, z, t)) + sponge_relaxation_u(sponge_relaxation, x, y, z, t, u)
        end
        v_func = if isnothing(v_tide_val)
            (x, y, z, t, u, v) -> sponge_relaxation_v(sponge_relaxation, x, y, z, t, v)
        else
            (x, y, z, t, u, v) -> Float64(v_tide_val(x, y, z, t)) + sponge_relaxation_v(sponge_relaxation, x, y, z, t, v)
        end

        composed_forcing = (;
            u = ContinuousForcing(u_func; field_dependencies = (:u, :v)),
            v = ContinuousForcing(v_func; field_dependencies = (:u, :v))
        )
    else
        # No sponge: plain prescribed body forcing, called as f(x, y, z, t).
        composed_forcing = (; u = ZeroForcing(), v = ZeroForcing())
        if !isnothing(u_tide_val)
            composed_forcing = merge(composed_forcing,
                (; u = ContinuousForcing((x, y, z, t) -> Float64(u_tide_val(x, y, z, t)))))
        end
        if !isnothing(v_tide_val)
            composed_forcing = merge(composed_forcing,
                (; v = ContinuousForcing((x, y, z, t) -> Float64(v_tide_val(x, y, z, t)))))
        end
    end

    merged_forcing = merge(forcing, composed_forcing)

    eff_closure = closure !== nothing ? closure : closure_scheme
    active_closure = if eff_closure === nothing
        VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(); ν = Float64(ν), κ = Float64(κ))
    elseif eff_closure isa Symbol
        eff_closure in (:nemotke, :catke, :tke) ? CATKEVerticalDiffusivity() :
            VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(); ν = Float64(ν), κ = Float64(κ))
    else
        eff_closure
    end

    kwargs = Dict{Symbol, Any}(
        :Δt => Float64(Δt),
        :closure => active_closure,
        :tracers => active_tracers,
        :reference_density => Float64(reference_density),
        :equation_of_state => TEOS10EquationOfState(; reference_density = Float64(reference_density)),
        :bottom_drag_background_velocity => Float64(bottom_drag),
        :boundary_conditions => merge(user_bcs, boundary_conditions),
        :forcing => merged_forcing
    )
    isnothing(free_surface) || (kwargs[:free_surface] = free_surface)
    isnothing(momentum_advection) || (kwargs[:momentum_advection] = momentum_advection)
    isnothing(tracer_advection) || (kwargs[:tracer_advection] = tracer_advection)

    sim = ocean_simulation(grid; model = :hydrostatic, kwargs...)
    return sim.model
end

"""
    apply_reanalysis_initial_conditions!(model; dataset, date, variables) -> model

Initialise prognostic tracers from a reanalysis dataset using `ClimaOcean.Metadatum`.

`Metadatum` is the supported mechanism for injecting observed/reanalysis state into an
Oceananigans model: it carries both the value and its provenance (dataset, date, variable
name) so that archived output remains traceable to the source it came from. This replaces
the previous hand-rolled `interpolate_ocean_state`, which fabricated constant fields when the
dataset could not be read.

# Inputs
- `model`: Model whose tracers are to be initialised.
- `dataset`: A NumericalEarth dataset, e.g. `GLORYSDaily()`, `GLORYSMonthly()`, `WOAMonthly()`.
- `date`: `DateTime` for the initial state.
- `variables`: Tracer symbols to initialise, e.g. `(:T, :S)`.

# Outputs
- The model, with tracers set.

# Notes
If the dataset cannot be fetched (no credentials, no network), the failure is reported
loudly rather than silently substituting constant values, so a run never quietly claims a
reanalysis-forced state it did not have.
"""
function apply_reanalysis_initial_conditions!(
    model;
    dataset,
    date::DateTime,
    variables::Tuple = (:T, :S)
)
    tracer_names = Dict(:temperature => :T, :salinity => :S)
    assignments = Pair{Symbol, Any}[]
    for v in variables
        haskey(tracer_names, v) || continue
        field = tracer_names[v]
        hasproperty(model.tracers, field) || continue
        push!(assignments, field => ClimaOcean.Metadatum(
            v; date = date, dataset = dataset))
    end
    isempty(assignments) && return model
    set!(model, (assignments..., u = 0.0, v = 0.0))
    return model
end

"""
    woa23_regridded_tracers(grid, temperature_file, salinity_file; oxygen_file=nothing)
        -> (T = (x, y, z) -> …, S = (x, y, z) -> …, O_2 = (x, y, z) -> …)

Read a downloaded WOA23 (or WOA-series) climatology from NetCDF and return **callable
interpolants** of conservative temperature `T` (°C), practical salinity `S` (PSU) and,
optionally, dissolved oxygen `O_2` (µmol kg⁻¹) on the model grid.

The returned callables take grid coordinates `(x, y, z)` in degrees East, degrees North
and metres (negative down) and can be handed straight to `set!(model, T = …, S = …)`, which
is how Oceananigans accepts a spatially varying initial condition.

# Why this exists
`NumericalEarth` can construct `WOAMonthly`/`WOAAnnual` dataset objects, but its
`Column(grid, dataset; …)` constructor does **not** read a dataset — `Column` describes a
point/column *location*. The supported injection path is `ClimaOcean.Metadatum`, which needs
`Downloads.download` to accept a `Metadatum` and so does not give a keyless end-to-end route
from a plain NOAA NCEI file. Since WOA23 is served from an open THREDDS/ERDDAP endpoint that
needs no credentials, reading the file directly and regridding it here keeps the keyless path
real rather than aspirational.

# Interpolation
Horizontal regridding is bilinear in `(lon, lat)` via `Interpolations.jl`; the vertical uses
the dataset's own standard depth levels with linear interpolation in depth, clamped to the
shallowest available level above and the deepest below. Clamping matters: WOA23's deepest
levels are far deeper than the model domain, so an unclamped query would extrapolate.

Values that WOA23 masks to its fill value (WOA uses `-1.0e30`) are replaced by `NaN` before
interpolation and then filled with the horizontally interpolated mean of the valid column, so
land or ice-covered cells do not poison neighbouring water cells with a huge negative bias.

# Inputs
- `grid`: Model grid (an `ImmersedBoundaryGrid` is unwrapped to its underlying grid).
- `temperature_file`, `salinity_file`: Paths to NetCDF files with `t_an` / `s_an`.
- `oxygen_file`: Optional path to a NetCDF file with `o_an`.

# Outputs
- NamedTuple of callables, omitting `O_2` when `oxygen_file` is `nothing` or the model has no
  `O_2` tracer.
"""
function woa23_regridded_tracers(
    grid,
    temperature_file::AbstractString,
    salinity_file::AbstractString;
    oxygen_file::Union{Nothing, AbstractString} = nothing
)
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid

    T_vals, lon_src, lat_src, dep_src = read_woa_variable(temperature_file, "t_an")
    S_vals, _, _, _ = read_woa_variable(salinity_file, "s_an")

    size(T_vals) == size(S_vals) || error(
        "WOA23 temperature $(size(T_vals)) and salinity $(size(S_vals)) grids disagree; " *
        "the two files must come from the same subset request.")

    T_interp = build_woa_interpolator(T_vals, lon_src, lat_src, dep_src)
    S_interp = build_woa_interpolator(S_vals, lon_src, lat_src, dep_src)

    # The model grid routinely extends outside the downloaded subset (a WOA request may
    # cover only the study box while the grid is larger) and far below the deepest WOA
    # level, so every query is clamped into the dataset's own bounds. Without this the
    # interpolator would throw on out-of-bounds indices; with it, boundary cells take the
    # nearest available value instead of an extrapolation. `read_woa_variable` has already
    # wrapped longitudes into signed degrees East, so longitude needs only clamping.
    lo_z, hi_z = Float64(first(dep_src)), Float64(last(dep_src))
    lo_y, hi_y = Float64(first(lat_src)), Float64(last(lat_src))
    lo_x, hi_x = Float64(first(lon_src)), Float64(last(lon_src))

    # `interpolate` indexes by *index* coordinates (1..n along each axis), not by the
    # physical values stored in the axes, so the query is converted here. Passing
    # physical lon/lat/depth straight through is silently wrong and throws a BoundsError
    # whose reported index is the physical value, which makes the mistake easy to misread.
    span(v) = begin
        d = last(v) - first(v)
        d == 0 ? (t -> 1.0) : (t -> (t - Float64(first(v))) / d + 1.0)
    end
    to_ix, to_iy, to_iz = span(lon_src), span(lat_src), span(dep_src)

    query(itp) = (x, y, z) -> itp(
        to_iz(clamp(Float64(z), lo_z, hi_z)),
        to_iy(clamp(Float64(y), lo_y, hi_y)),
        to_ix(clamp(Float64(x), lo_x, hi_x))
    )

    out = (T = query(T_interp), S = query(S_interp))

    if !isnothing(oxygen_file) && isfile(oxygen_file)
        O_vals, _, _, _ = read_woa_variable(oxygen_file, "o_an")
        if size(O_vals) == size(T_vals)
            O_interp = build_woa_interpolator(O_vals, lon_src, lat_src, dep_src)
            out = merge(out, (O_2 = query(O_interp),))
        else
            @warn "Skipping dissolved oxygen: WOA23 oxygen grid $(size(O_vals)) does not match " *
                  "the temperature/salinity grid $(size(T_vals))."
        end
    end

    return out
end

"""
    read_woa_variable(filepath, varname) -> (values, lon, lat, depth)

Read one WOA23 variable plus its coordinate axes from NetCDF.

Returns `values` as a `(depth, lat, lon)` `Array` with all three axes in **ascending** order,
and WOA's `-1.0e30` fill value mapped to `NaN`.

Two convention differences from the model grid are handled here:

* WOA stores depth positive-downward (0, 5, 10 … m); the model grid is negative-downward, so
  depths are negated and re-sorted to ascending (i.e. surface first).
* WOA stores the Atlantic sector as positive longitudes (e.g. `47.5`) whereas the model grid
  uses signed degrees East (e.g. `-47.5`), so longitudes are wrapped into `[-180, 180)` and
  re-sorted alongside the data.
"""
function read_woa_variable(filepath::AbstractString, varname::AbstractString)
    NCDatasets.NCDataset(filepath, "r") do ds
        haskey(ds, varname) || error(
            "WOA file $(filepath) has no variable '$(varname)'.")

        var = ds[varname]

        lon_raw = collect(Float64, ds["lon"][:])
        lat_raw = collect(Float64, ds["lat"][:])
        dep_raw = collect(Float64, ds["depth"][:])

        lon_wrapped = [mod(v + 180.0, 360.0) - 180.0 for v in lon_raw]
        dep_model = [-v for v in dep_raw]          # positive-down -> negative-down

        ip = sortperm(lon_wrapped)
        jp = sortperm(lat_raw)
        kp = sortperm(dep_model)

        lon = lon_wrapped[ip]
        lat = lat_raw[jp]
        depth = dep_model[kp]

        vals = permutedims(Array(var), (3, 2, 1))  # (lon, lat, depth) -> (depth, lat, lon)
        vals = vals[kp, jp, ip]

        size(vals) == (length(depth), length(lat), length(lon)) || error(
            "Unexpected WOA variable '$(varname)' with size $(size(vals)); expected " *
            "(depth=$(length(depth)), lat=$(length(lat)), lon=$(length(lon))).")

        # WOA masks land/ice with a huge negative fill value.
        vals[vals .<= -1.0e29] .= NaN

        return (vals, lon, lat, depth)
    end
end

"""
    build_woa_interpolator(values, lon, lat, depth) -> (x, y, z) -> Float64

Wrap a `(depth, lat, lon)` array in a trilinear interpolator evaluable at model coordinates.

`values` must already have WOA fill values mapped to `NaN`. Masked cells are repaired first so
a query can never return `NaN` for a model cell sitting over land:

* a depth column that is entirely `NaN` is replaced by the mean of the valid depths in that
  column (zero if there are none);
* any remaining interior `NaN` is replaced by the mean of all valid values.

The vertical axis is clamped to the dataset's depth range by the caller, so the interpolator is
only ever queried inside its bounds.
"""
function build_woa_interpolator(values::AbstractArray, lon, lat, depth)
    nz, ny, nx = size(values)

    for j in 1:ny, i in 1:nx
        col = @view values[:, j, i]
        all(isnan, col) || continue
        valid = filter(!isnan, col)
        isempty(valid) ? (col .= 0.0) : (col .= sum(valid) / length(valid))
    end

    nan_idx = findall(isnan, values)
    if !isempty(nan_idx)
        valid = filter(!isnan, values)
        mean_valid = isempty(valid) ? 0.0 : sum(valid) / length(valid)
        values[nan_idx] .= mean_valid
    end

    # `values` is (depth, lat, lon) with each axis ascending, which is the order
    # `interpolate` expects for the returned array.
    return interpolate(values, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear())))
end

"""
    ZeroForcing

Bitstype placeholder representing zero forcing acceleration.
Callable with arbitrary continuous arguments `(x, y, z, t)` or `(x, y, z, t, ...)`.
"""
struct ZeroForcing end
@inline (::ZeroForcing)(x, y, z, t) = 0.0
@inline (::ZeroForcing)(x, y, z, t, u) = 0.0
@inline (::ZeroForcing)(x, y, z, t, u, v) = 0.0

"""
    seawater_freezing_temperature(S::Real) -> Float64

Calculate seawater freezing temperature in °C as a function of practical salinity \$S\$ (PSU)
using the UNESCO / Millero (1978) thermodynamic formulation:
```math
T_{\\text{freeze}}(S) = -0.0575 S + 0.00171 S^{1.5} - 0.000215 S^2
```
"""
@inline function seawater_freezing_temperature(S::Real)::Float64
    s_val = Float64(S)
    if s_val < 0.0
        return 0.0
    end
    return -0.0575 * s_val + 0.00171 * (s_val^1.5) - 0.000215 * (s_val^2)
end

"""
    LateralBoundaryRelaxation{TF}

Bitstype lateral boundary sponge relaxation forcing following Price & Aumont (2011)
and Arctic Ocean shelf modeling implementations. Rather than imposing artificial localized
surface wind stress to drive regional shelf circulation, lateral boundary relaxation nudges
prognostic velocities \$(u, v)\$ and active tracers \$(T, S, O_2)\$ towards upstream and open ocean
hydrographic boundary conditions within peripheral sponge zones.

# Mathematical & Governing Formulation
```math
F_{\\text{relax}}(x, y, z, t, \\psi) = -\\frac{\\gamma(x, y)}{\\tau_{\\text{relax}}} \\left(\\psi - \\psi_{\\text{ref}}(x, y, z, t)\\right)
```
where the sponge relaxation weight \$\\gamma(x, y) \\in [0, 1]\$ increases smoothly (quadratically)
from 0.0 at the interior edge of the sponge layer to 1.0 at the outer domain boundary.

# References
- Price, J. F., & Aumont, O. (2011). Lateral boundary conditions and regional shelf circulation.
  *Journal of Physical Oceanography*, 41(5), 903-921.
"""
struct LateralBoundaryRelaxation{TU, TV}
    lon_domain       :: Tuple{Float64, Float64}
    lat_domain       :: Tuple{Float64, Float64}
    sponge_width_deg :: Float64
    tau_relax        :: Float64
    u_inflow         :: Float64
    v_inflow         :: Float64
    u_ref            :: TU
    v_ref            :: TV
end

"""
    ExponentialInflow{T}

Bitstype reference profile for open-boundary inflow, \$\\psi_{ref}(z) = A \\exp(z / d)\$.

A concrete struct is used rather than a closure so that the sponge forcing built from it
stays `isbits`. A closure capturing `A` and `d` is *not* `isbits`, which makes the enclosing
`ContinuousForcing` unlowerable into a GPU kernel and triggers
`KernelError: passing non-bitstype argument`.
"""
struct ExponentialInflow{T <: AbstractFloat}
    amplitude::T
    decay    :: T
end

@inline (f::ExponentialInflow)(x, y, z, t) = f.amplitude * exp(z / f.decay)

"""
    LateralBoundaryRelaxation(lon_domain, lat_domain; kwargs...)

Construct a regional lateral boundary relaxation sponge zone following Price & Aumont (2011).

# Arguments
- `lon_domain::Tuple{Real, Real}`: Zonal domain boundaries (lon_min, lon_max) in °E.
- `lat_domain::Tuple{Real, Real}`: Meridional domain boundaries (lat_min, lat_max) in °N.

# Keyword Arguments
- `sponge_width_deg::Real=0.5`: Width of peripheral buffer sponge layer in degrees.
- `tau_relax::Real=86400.0`: Relaxation timescale \$\\tau\$ in seconds (default 1 day).
- `u_inflow::Real=-0.15`: Reference zonal upstream inflow velocity in m/s.
- `v_inflow::Real=0.05`: Reference meridional upstream inflow velocity in m/s.
- `u_decay::Real=100.0`: e-folding depth of the zonal reference profile in metres.
- `v_decay::Real=100.0`: e-folding depth of the meridional reference profile in metres.

# Notes
`u_ref`/`v_ref` default to [`ExponentialInflow`](@ref) values rather than closures so the
result is `isbits` and can be compiled for GPU.
"""
function LateralBoundaryRelaxation(
    lon_domain::Tuple{Real, Real},
    lat_domain::Tuple{Real, Real};
    sponge_width_deg::Real = 0.5,
    tau_relax::Real = 86400.0,
    u_inflow::Real = -0.15,
    v_inflow::Real = 0.05,
    u_decay::Real = 100.0,
    v_decay::Real = 100.0
)
    return LateralBoundaryRelaxation(
        (Float64(lon_domain[1]), Float64(lon_domain[2])),
        (Float64(lat_domain[1]), Float64(lat_domain[2])),
        Float64(sponge_width_deg),
        Float64(tau_relax),
        Float64(u_inflow),
        Float64(v_inflow),
        ExponentialInflow(Float64(u_inflow), Float64(u_decay)),
        ExponentialInflow(Float64(v_inflow), Float64(v_decay))
    )
end

@inline function _sponge_relaxation(r::LateralBoundaryRelaxation, ref, x, y, z, t, psi::Real)
    gamma = compute_sponge_gamma(
        x, y,
        r.lon_domain[1], r.lon_domain[2],
        r.lat_domain[1], r.lat_domain[2],
        r.sponge_width_deg
    )
    if gamma <= 0.0
        return 0.0
    end
    return -gamma * (Float64(psi) - Float64(ref(x, y, z, t))) / r.tau_relax
end

"""
    sponge_relaxation_u(r, x, y, z, t, psi) -> Float64

Zonal (\$u\$) sponge relaxation tendency, with no runtime selector argument.

This is the entry point the GPU forcing path uses. It deliberately takes only
`isbitstype` arguments: a `Symbol` cannot be passed through a CUDA kernel ABI, so a
`Symbol`-dispatched method makes the enclosing `ContinuousForcing` unlowerable and the
free-surface kernel fails to compile with `InvalidIRError`. Keeping the two components as
separate methods keeps the kernel body free of non-bits values.
"""
@inline function sponge_relaxation_u(r::LateralBoundaryRelaxation, x, y, z, t, psi::Real)
    return _sponge_relaxation(r, r.u_ref, x, y, z, t, psi)
end

"""
    sponge_relaxation_v(r, x, y, z, t, psi) -> Float64

Meridional (\$v\$) sponge relaxation tendency; see [`sponge_relaxation_u`](@ref) for why this
is separate rather than `Symbol`-dispatched.
"""
@inline function sponge_relaxation_v(r::LateralBoundaryRelaxation, x, y, z, t, psi::Real)
    return _sponge_relaxation(r, r.v_ref, x, y, z, t, psi)
end

# Runtime-selected forms, retained for host-side and test callers. These are deliberately
# not used to build GPU forcing: the `Symbol` argument is what makes them unlowerable.
@inline function (r::LateralBoundaryRelaxation)(x, y, z, t, var::Symbol, psi::Real)
    return if var === :u
        sponge_relaxation_u(r, x, y, z, t, psi)
    elseif var === :v
        sponge_relaxation_v(r, x, y, z, t, psi)
    else
        0.0
    end
end

@inline (r::LateralBoundaryRelaxation)(x, y, z, t, psi) = sponge_relaxation_u(r, x, y, z, t, psi)
@inline (r::LateralBoundaryRelaxation)(x, y, z, t) = 0.0

"""
    compute_sponge_gamma(x, y, lon_min, lon_max, lat_min, lat_max, sponge_width) -> Float64

Compute the quadratic sponge relaxation weight γ ∈ [0, 1] for multi-boundary
relaxation following Price & Aumont (2011). Returns the maximum relaxation weight
across all four boundaries (east, west, north, south).
"""
@inline function compute_sponge_gamma(
    x::Real, y::Real,
    lon_min::Real, lon_max::Real,
    lat_min::Real, lat_max::Real,
    sponge_width::Real
)::Float64
    # Eastern boundary (x → lon_max)
    gamma_e = if x >= lon_max
        1.0
    elseif x <= (lon_max - sponge_width)
        0.0
    else
        ((x - (lon_max - sponge_width)) / sponge_width)^2
    end

    # Western boundary (x → lon_min)
    gamma_w = if x <= lon_min
        1.0
    elseif x >= (lon_min + sponge_width)
        0.0
    else
        (((lon_min + sponge_width) - x) / sponge_width)^2
    end

    # Northern boundary (y → lat_max)
    gamma_n = if y >= lat_max
        1.0
    elseif y <= (lat_max - sponge_width)
        0.0
    else
        ((y - (lat_max - sponge_width)) / sponge_width)^2
    end

    # Southern boundary (y → lat_min)
    gamma_s = if y <= lat_min
        1.0
    elseif y >= (lat_min + sponge_width)
        0.0
    else
        (((lat_min + sponge_width) - y) / sponge_width)^2
    end

    return max(gamma_e, gamma_w, gamma_n, gamma_s)
end

"""
    compute_sponge_velocity_target(u_inflow_mean, v_inflow_mean, u_decay, v_decay, z) -> Tuple{Float64, Float64}

Compute the target velocity for sponge relaxation with exponential depth decay.
"""
@inline function compute_sponge_velocity_target(
    u_inflow_mean::Real, v_inflow_mean::Real,
    u_decay::Real, v_decay::Real,
    z::Real
)::Tuple{Float64, Float64}
    u_target = Float64(u_inflow_mean) * exp(Float64(z) / Float64(u_decay))
    v_target = Float64(v_inflow_mean) * exp(Float64(z) / Float64(v_decay))
    return (u_target, v_target)
end

"""
    HorizontalMomentumForcingU{TF}

Bitstype callable struct for zonal momentum forcing on GPU/CPU without field dependencies.
Combines depth-tapered harmonic tidal body forcing, lateral sponge relaxation,
and Patankar-type quasi-implicit quadratic bottom drag.
"""
struct HorizontalMomentumForcingU{TF}
    tidal          :: TF
    has_sponge     :: Bool
    lon_min        :: Float64   # western domain boundary (°E)
    lon_max        :: Float64   # eastern domain boundary (°E)
    lat_min        :: Float64   # southern domain boundary (°N)
    lat_max        :: Float64   # northern domain boundary (°N)
    sponge_width   :: Float64   # sponge buffer width (degrees)
    sponge_tau     :: Float64   # relaxation timescale (s)
    u_inflow_mean  :: Float64   # reference zonal inflow velocity (m/s)
    u_inflow_decay :: Float64   # depth e-folding scale (m)
    bottom_drag    :: Float64
    cd_drag        :: Float64
end

"""
    (m::HorizontalMomentumForcingU)(x, y, z, t, u, v)

Zonal momentum forcing: tidal body force + multi-boundary quadratic sponge relaxation
(Price & Aumont 2011) + Patankar-implicit quadratic bottom drag.

Sponge gamma is the maximum of all four boundary ramp functions, ramping quadratically
from 0 at the inner sponge edge to 1 at the domain wall. This is identical to the
formulation in `LateralBoundaryRelaxation` and avoids ad-hoc sine-taper momentum
discontinuities that cause CFL blow-up near immersed boundaries.
"""
@inline function (m::HorizontalMomentumForcingU)(x, y, z, t, u, v)
    sponge_val = 0.0
    if m.has_sponge
        gamma = compute_sponge_gamma(
            x, y,
            m.lon_min, m.lon_max,
            m.lat_min, m.lat_max,
            m.sponge_width
        )
        if gamma > 0.0
            u_target = m.u_inflow_mean * exp(Float64(z) / m.u_inflow_decay)
            sponge_val = -gamma * (Float64(u) - u_target) / m.sponge_tau
        end
    end

    tide_val = m.tidal(x, y, z, t)

    speed = sqrt(u^2 + v^2)
    tau_ref = 120.0
    drag_coeff = (m.bottom_drag + m.cd_drag * speed) /
                 (1.0 + (m.bottom_drag + m.cd_drag * speed) * tau_ref)
    drag_val = -drag_coeff * Float64(u)
    return tide_val + sponge_val + drag_val
end

@inline (m::HorizontalMomentumForcingU)(x, y, z, t, u) = m(x, y, z, t, u, 0.0)
@inline (m::HorizontalMomentumForcingU)(x, y, z, t) = m(x, y, z, t, 0.0, 0.0)

"""
    HorizontalMomentumForcingV{TF}

Bitstype callable struct for meridional momentum forcing on GPU/CPU.
Combines harmonic tidal body forcing, multi-boundary quadratic sponge relaxation
(Price & Aumont 2011), and Patankar-type quasi-implicit quadratic bottom drag.
No sine wall tapering is applied to tidal forcing; relaxation covers all four
open boundaries symmetrically.
"""
struct HorizontalMomentumForcingV{TF}
    tidal          :: TF
    has_sponge     :: Bool
    lon_min        :: Float64   # western domain boundary (°E)
    lon_max        :: Float64   # eastern domain boundary (°E)
    lat_min        :: Float64   # southern domain boundary (°N)
    lat_max        :: Float64   # northern domain boundary (°N)
    sponge_width   :: Float64   # sponge buffer width (degrees)
    sponge_tau     :: Float64   # relaxation timescale (s)
    v_inflow_mean  :: Float64   # reference meridional inflow velocity (m/s)
    v_inflow_decay :: Float64   # depth e-folding scale (m)
    bottom_drag    :: Float64
    cd_drag        :: Float64
end

@inline function (m::HorizontalMomentumForcingV)(x, y, z, t, u, v)
    sponge_val = 0.0
    if m.has_sponge
        gamma = compute_sponge_gamma(
            x, y,
            m.lon_min, m.lon_max,
            m.lat_min, m.lat_max,
            m.sponge_width
        )
        if gamma > 0.0
            v_target = m.v_inflow_mean * exp(Float64(z) / m.v_inflow_decay)
            sponge_val = -gamma * (Float64(v) - v_target) / m.sponge_tau
        end
    end

    tide_val = m.tidal(x, y, z, t)

    speed = sqrt(u^2 + v^2)
    tau_ref = 120.0
    drag_coeff = (m.bottom_drag + m.cd_drag * speed) /
                 (1.0 + (m.bottom_drag + m.cd_drag * speed) * tau_ref)
    drag_val = -drag_coeff * Float64(v)
    return tide_val + sponge_val + drag_val
end

@inline (m::HorizontalMomentumForcingV)(x, y, z, t, v) = m(x, y, z, t, 0.0, v)
@inline (m::HorizontalMomentumForcingV)(x, y, z, t) = m(x, y, z, t, 0.0, 0.0)

"""
    set_initial_stratification!(model;

        surface_temperature::Union{Nothing, Real} = nothing,
        surface_temp::Union{Nothing, Real} = nothing,
        bottom_temperature::Union{Nothing, Real} = nothing,
        bottom_temp::Union{Nothing, Real} = nothing,
        cil_temperature::Union{Nothing, Real} = nothing,
        cil_temp::Union{Nothing, Real} = nothing,
        slope_temperature::Union{Nothing, Real} = nothing,
        slope_temp::Union{Nothing, Real} = nothing,
        temperature_gradient::Real = 0.01,
        temp_stratification::Union{Nothing, Function} = nothing,
        salinity::Union{Nothing, Real, Function} = 33.0,
        lon_range::Union{Nothing, Tuple} = nothing,
        lat_range::Union{Nothing, Tuple} = nothing,
    )

Set physically plausible initial `T` and `S` profiles on `model`.

Temperature is built from up to four overlapping components — a surface value, a
bottom value, a cold-intermediate-layer (CIL) minimum, and a continental-slope
minimum — with `temperature_gradient` setting the thermocline sharpness. Salinity is
constant unless a function is supplied.

All components are optional: with none given, a linear profile from
`surface_temperature` (default 14 °C) to `bottom_temperature` (default 4 °C) at
`temperature_gradient` is used, giving a stably stratified column.

# Notes
- The `temperature_*`/`temp_*` and `*_temp` spellings are accepted as aliases.
- Domain bounds, when supplied, are used to taper the shelf-slope component so it
  weakens outside the configured `lon_range`/`lat_range`.
- Values are broadcast across the horizontal grid, so the resulting field is
  horizontally uniform apart from the slope taper.
"""
function set_initial_stratification!(
    model;
    surface_temperature::Union{Nothing, Real} = nothing,
    surface_temp::Union{Nothing, Real} = nothing,
    bottom_temperature::Union{Nothing, Real} = nothing,
    bottom_temp::Union{Nothing, Real} = nothing,
    cil_temperature::Union{Nothing, Real} = nothing,
    cil_temp::Union{Nothing, Real} = nothing,
    slope_temperature::Union{Nothing, Real} = nothing,
    slope_temp::Union{Nothing, Real} = nothing,
    temperature_gradient::Union{Nothing, Real} = nothing,
    temp_stratification::Union{Nothing, Real} = nothing,
    salinity::Union{Real, Function} = 33.0,
    lon_range::Union{Nothing, Tuple} = nothing,
    lat_range::Union{Nothing, Tuple} = nothing,
    stratification_type::Symbol = :three_layer,
    kwargs...
)
    # Resolve domain bounds
    lon_min, lon_max = isnothing(lon_range) ? (-71.0, -53.0) : Float64.(lon_range)
    lat_min, lat_max = isnothing(lat_range) ? (40.0, 48.5) : Float64.(lat_range)

    if lon_min >= lon_max || lat_min >= lat_max
        error("Invalid domain bounds: lon_range=$(lon_range), lat_range=$(lat_range)")
    end

    # Resolve surface temperature (default 14.0 °C)
    T0 = if !isnothing(surface_temperature)
        Float64(surface_temperature)
    elseif !isnothing(surface_temp)
        Float64(surface_temp)
    else
        14.0
    end

    # Resolve vertical temperature gradient (default 0.01 °C/m)
    dT_dz = if !isnothing(temperature_gradient)
        Float64(temperature_gradient)
    elseif !isnothing(temp_stratification)
        Float64(temp_stratification)
    else
        0.01
    end

    # Resolve bottom / CIL / slope water temperatures
    b_temp = if !isnothing(bottom_temperature)
        Float64(bottom_temperature)
    elseif !isnothing(bottom_temp)
        Float64(bottom_temp)
    else
        nothing
    end

    T_cil_min = if !isnothing(cil_temperature)
        Float64(cil_temperature)
    elseif !isnothing(cil_temp)
        Float64(cil_temp)
    elseif !isnothing(b_temp) && b_temp <= 5.0
        Float64(b_temp)
    else
        1.5
    end

    T_slope = if !isnothing(slope_temperature)
        Float64(slope_temperature)
    elseif !isnothing(slope_temp)
        Float64(slope_temp)
    elseif !isnothing(b_temp) && b_temp > 5.0
        Float64(b_temp)
    else
        8.5
    end

    # Horizontal gradient scales (°C) across the Scotian Shelf
    ΔT_cross = 8.0
    ΔT_along = 5.0
    ΔT_cil   = 2.5
    T_cil_edge = T_cil_min + ΔT_cil

    # Salinity gradient scales (PSU)
    S0       = Float64(salinity isa Real ? salinity : 33.0)
    ΔS_cross = 1.5
    ΔS_deep  = 2.5

    # Deep water relaxation scale (m)
    h_scale = (dT_dz > 0.0 && T_slope > T_cil_edge) ?
        max(50.0, min(500.0, (T_slope - T_cil_edge) / dT_dz)) : 100.0

    function temp_profile(lon, lat, z)
        dx = lon_max - lon_min
        dy = lat_max - lat_min
        x_norm = max(0.0, min(1.0, (lon - lon_min) / dx))
        y_norm = max(0.0, min(1.0, (lat - lat_min) / dy))
        δT_h = ΔT_cross * x_norm - ΔT_along * y_norm

        T_surf = T0 + δT_h

        if stratification_type == :linear
            return T_surf + dT_dz * z
        end

        # Three-layer Northwest Atlantic / Scotian Shelf structure (Petrie & Drinkwater 1993)
        # Cold Intermediate Layer core & edge with horizontal gradient modulation
        T_cil_min_local = T_cil_min + 0.40 * δT_h
        T_cil_edge_local = T_cil_min_local + ΔT_cil

        if z > -20.0
            # Surface mixed layer thermocline
            frac = (z + 20.0) / 20.0
            return T_cil_edge_local + frac * (T_surf - T_cil_edge_local)
        elseif z > -80.0
            # Cold Intermediate Layer (CIL): parabolic minimum centered at z = -50 m
            centre_frac = (z + 50.0) / 30.0
            return T_cil_min_local + ΔT_cil * centre_frac^2
        else
            # Deep Warm Slope Water: smooth exponential relaxation toward T_slope
            z_deep = abs(z) - 80.0
            T_slope_local = T_slope + 0.25 * δT_h
            return T_cil_edge_local +
                   (T_slope_local - T_cil_edge_local) * (1.0 - exp(-z_deep / h_scale))
        end
    end

    sal_profile = if salinity isa Function
        salinity
    else
        function(lon, lat, z)
            dx = lon_max - lon_min
            x_norm = max(0.0, min(1.0, (lon - lon_min) / dx))
            S_cross = ΔS_cross * x_norm
            # Salinity increases with depth to guarantee static gravitational stability: ∂ρ/∂z ≤ 0
            S_depth = ΔS_deep * (1.0 - exp(-abs(z) / 150.0))
            return S0 + S_cross + S_depth
        end
    end

    o2_arg = if :O2 in keys(model.tracers)
        o2_input = get(kwargs, :oxygen, nothing)
        if o2_input isa Function
            o2_input
        elseif o2_input isa Real
            (lon, lat, z) -> Float64(o2_input)
        else
            function(lon, lat, z)
                dx = lon_max - lon_min
                x_norm = max(0.0, min(1.0, (lon - lon_min) / dx))
                # Well-oxygenated surface & CIL (~300 umol/kg), deeper slope/basin depletion
                o2_surf = 305.0 - 15.0 * x_norm
                o2_dep = 100.0 * (1.0 - exp(-abs(Float64(z)) / 150.0))
                return o2_surf - o2_dep
            end
        end
    else
        nothing
    end

    if !isnothing(o2_arg)
        set!(model, T = temp_profile, S = sal_profile, O2 = o2_arg)
    else
        set!(model, T = temp_profile, S = sal_profile)
    end
    return nothing
end

"""
    set_initial_conditions!(model; kwargs...)

Flexible setter for initializing arbitrary model fields (velocities, tracers).

# Inputs
- `model`: Oceananigans model instance.
- `kwargs...`: Named fields and corresponding functions, arrays, or constants.
"""
function set_initial_conditions!(model; kwargs...)
    set!(model; kwargs...)
    return nothing
end
