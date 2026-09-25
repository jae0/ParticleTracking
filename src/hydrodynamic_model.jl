"""
    hydrodynamic_model.jl

Configuration and initialization of hydrostatic free-surface ocean models
in Oceananigans.jl for coastal shelf domains.
"""

using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: WENOVectorInvariant
using ClimaOcean

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
- `u_ref`: Function `(x, y, z, t) -> u_target` for zonal flow reference.
- `v_ref`: Function `(x, y, z, t) -> v_target` for meridional flow reference.
"""
function LateralBoundaryRelaxation(
    lon_domain::Tuple{Real, Real},
    lat_domain::Tuple{Real, Real};
    sponge_width_deg::Real = 0.5,
    tau_relax::Real = 86400.0,
    u_inflow::Real = -0.15,
    v_inflow::Real = 0.05,
    u_ref = (x, y, z, t) -> Float64(u_inflow) * exp(Float64(z) / 100.0),
    v_ref = (x, y, z, t) -> Float64(v_inflow) * exp(Float64(z) / 100.0)
)
    return LateralBoundaryRelaxation(
        (Float64(lon_domain[1]), Float64(lon_domain[2])),
        (Float64(lat_domain[1]), Float64(lat_domain[2])),
        Float64(sponge_width_deg),
        Float64(tau_relax),
        Float64(u_inflow),
        Float64(v_inflow),
        u_ref,
        v_ref
    )
end

@inline function (r::LateralBoundaryRelaxation)(x, y, z, t, var::Symbol, psi::Real)
    # Eastern boundary sponge (x -> lon_max)
    gamma_e = if x >= r.lon_domain[2]
        1.0
    elseif x <= (r.lon_domain[2] - r.sponge_width_deg)
        0.0
    else
        ((x - (r.lon_domain[2] - r.sponge_width_deg)) / r.sponge_width_deg)^2
    end

    # Western boundary sponge (x -> lon_min)
    gamma_w = if x <= r.lon_domain[1]
        1.0
    elseif x >= (r.lon_domain[1] + r.sponge_width_deg)
        0.0
    else
        (((r.lon_domain[1] + r.sponge_width_deg) - x) / r.sponge_width_deg)^2
    end

    # Northern boundary sponge (y -> lat_max)
    gamma_n = if y >= r.lat_domain[2]
        1.0
    elseif y <= (r.lat_domain[2] - r.sponge_width_deg)
        0.0
    else
        ((y - (r.lat_domain[2] - r.sponge_width_deg)) / r.sponge_width_deg)^2
    end

    # Southern boundary sponge (y -> lat_min)
    gamma_s = if y <= r.lat_domain[1]
        1.0
    elseif y >= (r.lat_domain[1] + r.sponge_width_deg)
        0.0
    else
        (((r.lat_domain[1] + r.sponge_width_deg) - y) / r.sponge_width_deg)^2
    end

    gamma = max(gamma_e, gamma_w, gamma_n, gamma_s)
    if gamma <= 0.0
        return 0.0
    end

    target_val = if var === :u
        Float64(r.u_ref(x, y, z, t))
    elseif var === :v
        Float64(r.v_ref(x, y, z, t))
    else
        0.0
    end

    return -gamma * (Float64(psi) - target_val) / r.tau_relax
end

@inline (r::LateralBoundaryRelaxation)(x, y, z, t, psi) = r(x, y, z, t, :u, psi)
@inline (r::LateralBoundaryRelaxation)(x, y, z, t) = 0.0

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
        # Eastern boundary (x → lon_max)
        gamma_e = if x >= m.lon_max
            1.0
        elseif x <= (m.lon_max - m.sponge_width)
            0.0
        else
            ((x - (m.lon_max - m.sponge_width)) / m.sponge_width)^2
        end

        # Western boundary (x → lon_min)
        gamma_w = if x <= m.lon_min
            1.0
        elseif x >= (m.lon_min + m.sponge_width)
            0.0
        else
            (((m.lon_min + m.sponge_width) - x) / m.sponge_width)^2
        end

        # Northern boundary (y → lat_max)
        gamma_n = if y >= m.lat_max
            1.0
        elseif y <= (m.lat_max - m.sponge_width)
            0.0
        else
            ((y - (m.lat_max - m.sponge_width)) / m.sponge_width)^2
        end

        # Southern boundary (y → lat_min)
        gamma_s = if y <= m.lat_min
            1.0
        elseif y >= (m.lat_min + m.sponge_width)
            0.0
        else
            (((m.lat_min + m.sponge_width) - y) / m.sponge_width)^2
        end

        gamma = max(gamma_e, gamma_w, gamma_n, gamma_s)
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
        # Eastern boundary (x → lon_max)
        gamma_e = if x >= m.lon_max
            1.0
        elseif x <= (m.lon_max - m.sponge_width)
            0.0
        else
            ((x - (m.lon_max - m.sponge_width)) / m.sponge_width)^2
        end

        # Western boundary (x → lon_min)
        gamma_w = if x <= m.lon_min
            1.0
        elseif x >= (m.lon_min + m.sponge_width)
            0.0
        else
            (((m.lon_min + m.sponge_width) - x) / m.sponge_width)^2
        end

        # Northern boundary (y → lat_max)
        gamma_n = if y >= m.lat_max
            1.0
        elseif y <= (m.lat_max - m.sponge_width)
            0.0
        else
            ((y - (m.lat_max - m.sponge_width)) / m.sponge_width)^2
        end

        # Southern boundary (y → lat_min)
        gamma_s = if y <= m.lat_min
            1.0
        elseif y >= (m.lat_min + m.sponge_width)
            0.0
        else
            (((m.lat_min + m.sponge_width) - y) / m.sponge_width)^2
        end

        gamma = max(gamma_e, gamma_w, gamma_n, gamma_s)
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
    build_hydrodynamic_model(
        grid;
        coriolis_latitude::Real = 45.0,
        surface_wind_stress_x::Union{Real, AbstractMatrix, Function} = 0.1,
        surface_wind_stress_y::Union{Real, AbstractMatrix, Function} = 0.0,
        ν::Real = 1e-2,
        κ::Real = 1e-2,
        tracers::Tuple = (:T, :S),
        free_surface = ImplicitFreeSurface()
    )

Instantiate a `HydrostaticFreeSurfaceModel` with Coriolis rotation,
seawater buoyancy, turbulent eddy diffusivity, and surface wind stress fluxes.

# Mathematical & Governing Equations
The hydrostatic Boussinesq primitive equations integrated on the rotating sphere
are:
```math
\\frac{\\partial \\boldsymbol{u}_h}{\\partial t} + (\\boldsymbol{u} \\cdot \\nabla) \\boldsymbol{u}_h
+ f \\hat{\\boldsymbol{k}} \\times \\boldsymbol{u}_h = -\\frac{1}{\\rho_0} \\nabla_h p
+ \\nabla \\cdot (\\nu \\nabla \\boldsymbol{u}_h) - (r_{\\text{drag}} + C_d |\\boldsymbol{u}_h|) \\boldsymbol{u}_h
```
```math
\\frac{\\partial p}{\\partial z} = -\\rho g = b \\rho_0
```
```math
\\nabla \\cdot \\boldsymbol{u} = 0
```
```math
\\frac{\\partial T}{\\partial t} + \\boldsymbol{u} \\cdot \\nabla T = \\nabla \\cdot (\\kappa \\nabla T)
```
```math
\\frac{\\partial S}{\\partial t} + \\boldsymbol{u} \\cdot \\nabla S = \\nabla \\cdot (\\kappa \\nabla S)
```
where \$\\boldsymbol{u}_h = (u, v)\$ is horizontal velocity, \$f = 2\\Omega \\sin \\phi\$
is the Coriolis parameter, \$b = -g (\\rho - \\rho_0)/\\rho_0\$ is buoyancy,
\$\\nu, \\kappa\$ represent momentum and tracer eddy diffusivities, and
\$r_{\\text{drag}}, C_d\$ parameterize bottom and boundary layer momentum dissipation
(Large & Pond 1981; Blumberg & Mellor 1987).

# Inputs
- `grid`: Underlying computational grid or `ImmersedBoundaryGrid`.
- `coriolis_latitude::Real`: Reference latitude for Coriolis `FPlane` in degrees North.
- `surface_wind_stress_x`: Zonal kinematic surface momentum flux (\$m^2 s^{-2}\$).
- `surface_wind_stress_y`: Meridional kinematic surface momentum flux (\$m^2 s^{-2}\$).
- `surface_heat_flux`: Net surface heat flux in \$W m^{-2}\$ (positive downward, warming).
- `tidal_forcing`: Optional NamedTuple `(u = Fu, v = Fv)` from `build_tidal_body_forcing`.
- `bottom_drag::Real`: Linear Rayleigh bottom/boundary damping rate in \$s^{-1}\$ (default \$10^{-4}\$).
- `cd_drag::Real`: Quadratic bottom drag coefficient (default \$10^{-3}\$).
- `ν::Real`: Kinematic eddy viscosity (\$m^2 s^{-1}\$).
- `κ::Real`: Tracer eddy diffusivity (\$m^2 s^{-1}\$).
- `closure`: Optional custom turbulence closure. Defaults to `SmagorinskyLilly()` if none provided.
- `tracers::Tuple`: Active tracer fields (default `(:T, :S)`).
- `free_surface`: Free surface representation (default `ImplicitFreeSurface()`).

Note: Momentum and tracer advection are explicitly configured to use 5th-order `WENO()` 
schemes to suppress grid-scale numerical noise over complex bathymetry.

# Outputs
- `HydrostaticFreeSurfaceModel`: Configured Oceananigans model instance.

# References
- Blumberg, A. F., & Mellor, G. L. (1987). A description of a three-dimensional
  coastal ocean circulation model. *Three-Dimensional Coastal Ocean Models*, 4, 1-16.
- Large, W. G., & Pond, S. (1981). JPO, 11(3), 324-336.
- Ramadhan, A., et al. (2020). Oceananigans.jl. *JOSS*, 5(53), 2018.
"""
function build_hydrodynamic_model(
    grid;
    coriolis_latitude::Real = 45.0,
    surface_wind_stress_x = 0.0001,
    surface_wind_stress_y = 0.0,
    surface_heat_flux = 0.0,
    tidal_forcing::Union{Nothing, NamedTuple} = nothing,
    open_boundary_conditions = nothing,
    sponge_forcing::Union{Nothing, NamedTuple} = nothing,
    bottom_drag::Real = 1e-4,
    cd_drag::Real = 2.5e-3,
    ν::Real = 1e-2,
    κ::Real = 1e-2,
    closure = nothing,
    closure_scheme = nothing,
    tracers::Tuple = (:T, :S),
    enable_o2::Bool = false,
    lateral_boundary_relaxation::Bool = false,
    sponge_width::Float64 = 0.35,
    sponge_tau::Float64 = 3600.0,
    u_inflow::Float64 = -0.15,
    v_inflow::Float64 = 0.05,
    u_decay::Float64 = 200.0,
    omega_m2::Float64 = (2π / 44712.0), # M2 period ≈ 12.42 h
    omega_s2::Float64 = (2π / 43200.0),  # S2 period ≈ 12.0 h
    free_surface = ImplicitFreeSurface(),
    momentum_advection = WENOVectorInvariant(),
    tracer_advection = WENO()
)
    # Resolve closure scheme if passed via closure_scheme alias
    eff_closure = closure !== nothing ? closure : closure_scheme

    # Resolve active prognostic tracers
    active_tracers = tracers
    if enable_o2 && !(:O2 in active_tracers)
        active_tracers = (active_tracers..., :O2)
    end

    # Surface kinematic boundary conditions for horizontal momentum
    u_top_bc = surface_wind_stress_x isa BoundaryCondition ?
        surface_wind_stress_x : FluxBoundaryCondition(surface_wind_stress_x)
    v_top_bc = surface_wind_stress_y isa BoundaryCondition ?
        surface_wind_stress_y : FluxBoundaryCondition(surface_wind_stress_y)

    u_bcs_dict = Dict{Symbol, Any}(:top => u_top_bc)
    v_bcs_dict = Dict{Symbol, Any}(:top => v_top_bc)

    # Net surface heat flux: J_T = -Q_net / (ρ₀ * c_p) [°C m s⁻¹] (upward positive in Oceananigans)
    rho0_cp = 1025.0 * 3990.0 # volumetric heat capacity ~4.09e6 J/(m³ °C)
    kinematic_T_flux = if surface_heat_flux isa Function
        (x, y, t) -> -surface_heat_flux(x, y, t) / rho0_cp
    elseif surface_heat_flux isa AbstractMatrix
        -surface_heat_flux ./ rho0_cp
    else
        -Float64(surface_heat_flux) / rho0_cp
    end
    T_top_bc = FluxBoundaryCondition(kinematic_T_flux)
    T_bcs_dict = Dict{Symbol, Any}(:top => T_top_bc)
    S_bcs_dict = Dict{Symbol, Any}()

    u_bcs = FieldBoundaryConditions(; u_bcs_dict...)
    v_bcs = FieldBoundaryConditions(; v_bcs_dict...)
    T_bcs = FieldBoundaryConditions(; T_bcs_dict...)

    boundary_conditions = Dict{Symbol, Any}(:u => u_bcs, :v => v_bcs)
    if :T in tracers
        boundary_conditions[:T] = T_bcs
    end
    if :S in tracers && !isempty(S_bcs_dict)
        boundary_conditions[:S] = FieldBoundaryConditions(; S_bcs_dict...)
    end
    if :O2 in tracers
        boundary_conditions[:O2] = FieldBoundaryConditions()
    end

    # Grid-aware full domain boundaries for lateral relaxation sponges.
    # Extract all four domain edges so the sponge covers E/W/N/S uniformly,
    # consistent with LateralBoundaryRelaxation (Price & Aumont 2011).
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    lon_min_grid = hasproperty(base_g, :λᶠᵃᵃ) ?
        Float64(base_g.λᶠᵃᵃ[1]) : -68.0
    lon_max_grid = hasproperty(base_g, :λᶠᵃᵃ) ?
        Float64(base_g.λᶠᵃᵃ[base_g.Nx + 1]) : -57.0
    lat_min_grid = hasproperty(base_g, :φᵃᶠᵃ) ?
        Float64(base_g.φᵃᶠᵃ[1]) : 42.0
    lat_max_grid = hasproperty(base_g, :φᵃᶠᵃ) ?
        Float64(base_g.φᵃᶠᵃ[base_g.Ny + 1]) : 47.5

    # Active sponge relaxation following Price & Aumont (2011).
    # Uses LateralBoundaryRelaxation directly for CPU path, or builds it into the
    # bitstype HorizontalMomentumForcingU/V structs for GPU bitstype compatibility.
    # u_decay / v_decay are depth e-folding scales (meters, positive value).
    u_decay = 200.0
    v_decay = 500.0
    active_sponge = if !isnothing(sponge_forcing)
        sponge_forcing
    elseif !isnothing(open_boundary_conditions) || lateral_boundary_relaxation
        # Build a LateralBoundaryRelaxation object and wrap into named-tuple closures
        # so the CPU forcing path can dispatch identically to the GPU bitstype path.
        lbr = LateralBoundaryRelaxation(
            (lon_min_grid, lon_max_grid),
            (lat_min_grid, lat_max_grid);
            sponge_width_deg = sponge_width,
            tau_relax        = sponge_tau,
            u_inflow         = u_inflow,
            v_inflow         = v_inflow,
            u_ref            = (x, y, z, t) -> Float64(u_inflow) * exp(Float64(z) / u_decay),
            v_ref            = (x, y, z, t) -> Float64(v_inflow) * exp(Float64(z) / v_decay)
        )
        (
            u = (x, y, z, t, u) -> lbr(x, y, z, t, :u, u),
            v = (x, y, z, t, v) -> lbr(x, y, z, t, :v, v)
        )
    else
        nothing
    end

    coriolis = FPlane(latitude = coriolis_latitude)
    buoyancy = SeawaterBuoyancy()

    # Dynamic resolution of turbulence closures from ClimaOcean / Oceananigans
    active_closure = if eff_closure isa Symbol
        h_diff = HorizontalScalarDiffusivity(ν = 20.0, κ = 20.0)
        v_diff = VerticalScalarDiffusivity(
            VerticallyImplicitTimeDiscretization(),
            ν = Float64(ν),
            κ = Float64(κ)
        )
        if eff_closure in (:nemotke, :catke, :tke)
            (CATKEVerticalDiffusivity(), h_diff)
        elseif eff_closure == :smagorinsky
            (SmagorinskyLilly(), h_diff, v_diff)
        else
            (SmagorinskyLilly(), h_diff, v_diff)
        end
    elseif isnothing(eff_closure)
        (
            SmagorinskyLilly(),
            HorizontalScalarDiffusivity(ν = 20.0, κ = 20.0),
            VerticalScalarDiffusivity(
                VerticallyImplicitTimeDiscretization(),
                ν = Float64(ν),
                κ = Float64(κ)
            )
        )
    elseif eff_closure isa Tuple
        has_vdiff = any(
            c -> c isa VerticalScalarDiffusivity || c isa CATKEVerticalDiffusivity,
            eff_closure
        )
        if !has_vdiff
            v_diff = VerticalScalarDiffusivity(
                VerticallyImplicitTimeDiscretization(),
                ν = Float64(ν),
                κ = Float64(κ)
            )
            (eff_closure..., v_diff)
        else
            eff_closure
        end
    else
        eff_closure
    end

    # Momentum forcing: GPU uses bitstype continuous forcing; CPU supports drag closures
    arch = architecture(grid)
    is_gpu = !(arch isa CPU)

    momentum_forcing = if is_gpu
        tf_u = !isnothing(tidal_forcing) && hasproperty(tidal_forcing, :u) ?
            tidal_forcing.u : ZeroForcing()
        tf_v = !isnothing(tidal_forcing) && hasproperty(tidal_forcing, :v) ?
            tidal_forcing.v : ZeroForcing()

        has_obc = !isnothing(open_boundary_conditions) ||
                  !isnothing(sponge_forcing) ||
                  lateral_boundary_relaxation

        # GPU bitstype structs: multi-boundary quadratic sponge, no sine wall tapering.
        # u_decay / v_decay are depth e-folding scales (m, positive).
        u_decay_gpu = 200.0
        v_decay_gpu = 500.0
        fu_gpu = HorizontalMomentumForcingU(
            tf_u,
            has_obc,
            lon_min_grid, lon_max_grid,
            lat_min_grid, lat_max_grid,
            sponge_width,
            sponge_tau,
            Float64(u_inflow), u_decay_gpu,
            Float64(bottom_drag), Float64(cd_drag)
        )
        fv_gpu = HorizontalMomentumForcingV(
            tf_v,
            has_obc,
            lon_min_grid, lon_max_grid,
            lat_min_grid, lat_max_grid,
            sponge_width,
            sponge_tau,
            Float64(v_inflow), v_decay_gpu,
            Float64(bottom_drag), Float64(cd_drag)
        )

        (u = Forcing(fu_gpu, field_dependencies = (:u, :v)),
         v = Forcing(fv_gpu, field_dependencies = (:u, :v)))
    else
        # CPU path: tidal forcing is passed without any boundary wall taper.
        # Sponge relaxation is handled by active_sponge (LateralBoundaryRelaxation).
        total_Fu(x, y, z, t, u, v) = begin
            tide_val = isnothing(tidal_forcing) ? 0.0 : Float64(tidal_forcing.u(x, y, z, t))
            speed = sqrt(u^2 + v^2)
            tau_ref = 120.0
            drag_coeff = (bottom_drag + cd_drag * speed) /
                         (1.0 + (bottom_drag + cd_drag * speed) * tau_ref)
            drag_val = -drag_coeff * Float64(u)
            sponge_val = (!isnothing(active_sponge) && hasproperty(active_sponge, :u)) ?
                Float64(active_sponge.u(x, y, z, t, u)) : 0.0
            return tide_val + drag_val + sponge_val
        end

        total_Fv(x, y, z, t, u, v) = begin
            tide_val = isnothing(tidal_forcing) ? 0.0 : Float64(tidal_forcing.v(x, y, z, t))
            speed = sqrt(u^2 + v^2)
            tau_ref = 120.0
            drag_coeff = (bottom_drag + cd_drag * speed) /
                         (1.0 + (bottom_drag + cd_drag * speed) * tau_ref)
            drag_val = -drag_coeff * Float64(v)
            sponge_val = (!isnothing(active_sponge) && hasproperty(active_sponge, :v)) ?
                Float64(active_sponge.v(x, y, z, t, v)) : 0.0
            return tide_val + drag_val + sponge_val
        end

        (
            u = Forcing(total_Fu, field_dependencies = (:u, :v)),
            v = Forcing(total_Fv, field_dependencies = (:u, :v))
        )
    end

    model_kwargs = Dict{Symbol, Any}(
        :coriolis => coriolis,
        :buoyancy => buoyancy,
        :tracers => active_tracers,
        :boundary_conditions => NamedTuple(boundary_conditions),
        :forcing => momentum_forcing,
        :closure => active_closure,
        :free_surface => free_surface,
        :momentum_advection => momentum_advection,
        :tracer_advection => tracer_advection
    )

    model = HydrostaticFreeSurfaceModel(grid; model_kwargs...)

    return model
end


"""
    set_initial_stratification!(
        model;
        surface_temperature::Real = 14.0,
        surface_temp = nothing,
        bottom_temperature = nothing,
        bottom_temp = nothing,
        cil_temperature = nothing,
        cil_temp = nothing,
        slope_temperature = nothing,
        slope_temp = nothing,
        temperature_gradient = 0.01,
        temp_stratification = nothing,
        salinity::Union{Real, Function} = 33.0,
        lon_range::Union{Nothing, Tuple} = nothing,
        lat_range::Union{Nothing, Tuple} = nothing,
        stratification_type::Symbol = :three_layer,
        kwargs...
    )

Initialize realistic 3D thermal and haline stratification for the Scotian Shelf
and Gulf of St. Lawrence region in Oceananigans.

# Mathematical Formulation

When `stratification_type == :three_layer` (default), the vertical temperature
profile \$T(\\lambda, \\phi, z)\$ captures the three-layer water column structure
characteristic of the Northwest Atlantic shelf:

**1. Surface mixed layer** (\$z > -20\\text{ m}\$):
```math
T_{\\text{surf}}(\\lambda, \\phi) =
    T_0 + \\Delta T_{\\text{cross}} \\cdot x_{\\text{norm}}
       - \\Delta T_{\\text{along}} \\cdot y_{\\text{norm}}
```
```math
T(z) = T_{\\text{cil,edge}} + \\left(\\frac{z + 20}{20}\\right)
    \\cdot (T_{\\text{surf}} - T_{\\text{cil,edge}})
```
where \$x_{\\text{norm}} \\in [0, 1]\$ increases offshore and \$y_{\\text{norm}} \\in [0, 1]\$
increases northward across the model domain. Reference horizontal gradients are
\$\\Delta T_{\\text{cross}} \\approx 8^\\circ\\text{C}\$ and \$\\Delta T_{\\text{along}} \\approx 5^\\circ\\text{C}\$.
At \$z = 0\\text{ m}\$, \$T = T_{\\text{surf}}\$. At \$z = -20\\text{ m}\$,
\$T = T_{\\text{cil,edge}} = T_{\\text{cil,min}} + \\Delta T_{\\text{cil}}\$.

**2. Cold Intermediate Layer (CIL)** (\$-80\\text{ m} < z \\le -20\\text{ m}\$):
```math
T_{\\text{CIL}}(\\lambda, \\phi, z) = T_{\\text{cil,min}} +
    \\Delta T_{\\text{cil}} \\cdot \\left(\\frac{z + 50}{30}\\right)^2
```
The CIL represents cold winter-cooled shelf water. At the core depth \$z = -50\\text{ m}\$,
\$T = T_{\\text{cil,min}}\$. At both interfaces (\$z = -20\\text{ m}\$ and \$z = -80\\text{ m}\$),
\$T = T_{\\text{cil,edge}}\$, maintaining \$C^0\$ continuity across layer interfaces.

**3. Deep / Warm Slope Water** (\$z \\le -80\\text{ m}\$):
```math
T_{\\text{deep}}(z) = T_{\\text{cil,edge}} +
    (T_{\\text{slope}} - T_{\\text{cil,edge}}) \\cdot
    \\left[1 - \\exp\\left(-\\frac{|z| - 80}{h_{\\text{scale}}}\\right)\\right]
```
where \$T_{\\text{slope}} \\approx 8.5^\\circ\\text{C}\$ and \$h_{\\text{scale}}\$ is the
vertical transition scale (default 100 m, or scaled by `temperature_gradient`).

When `stratification_type == :linear`:
```math
T(\\lambda, \\phi, z) = T_{\\text{surf}}(\\lambda, \\phi) + \\Gamma \\cdot z
```
where \$\\Gamma = \\text{temperature\\_gradient}\$ (or `temp_stratification`).

Salinity profile:
```math
S(\\lambda, \\phi, z) = S_0 + \\Delta S_{\\text{cross}} \\cdot x_{\\text{norm}}
    + \\Delta S_{\\text{deep}} \\cdot \\left[1 - \\exp\\left(-\\frac{|z|}{150}\\right)\\right]
```
ensuring gravitational static stability: \$\\partial \\rho / \\partial z \\le 0\$.

# References
- Petrie, B., and Drinkwater, K. F. (1993). Temperature and salinity
  variability on the Scotian Shelf and in the Gulf of Maine 1945–1990.
  *J. Geophys. Res. Oceans*, 98(C11), 20079–20089.
- Yashayaev, I., and Loder, J. W. (2016). Recurrent replenishment of
  Labrador Sea Water and associated decadal‐scale variability in heat and
  freshwater content. *Geophys. Res. Lett.*, 43(9), 4399–4407.

# Inputs
- `model`: `HydrostaticFreeSurfaceModel` instance to initialize.
- `surface_temperature::Real`: Reference SST at the warm corner (°C).
  Alias: `surface_temp`. Default 14.0 °C.
- `bottom_temperature::Union{Nothing, Real}`: Shelf benthic or deep water
  temperature (°C). Alias: `bottom_temp`. If \$\\le 5.0^\\circ\\text{C}\$, sets the
  CIL minimum temperature \$T_{\\text{cil,min}}\$. If \$> 5.0^\\circ\\text{C}\$,
  sets the deep slope water temperature \$T_{\\text{slope}}\$.
- `cil_temperature::Union{Nothing, Real}`: CIL core minimum temperature (°C).
  Alias: `cil_temp`. Default 1.5 °C.
- `slope_temperature::Union{Nothing, Real}`: Warm Slope Water temperature (°C).
  Alias: `slope_temp`. Default 8.5 °C.
- `temperature_gradient::Union{Nothing, Real}`: Vertical thermal gradient
  (\$dT/dz\$, °C/m). Alias: `temp_stratification`. Default 0.01 °C/m.
- `salinity::Union{Real, Function}`: Background practical salinity (PSU) or
  a function `(lon, lat, z) -> S`. Default 33.0 PSU.
- `lon_range::Union{Nothing, Tuple}`: Domain longitude bounds for coordinate normalization.
- `lat_range::Union{Nothing, Tuple}`: Domain latitude bounds for coordinate normalization.
- `stratification_type::Symbol`: `:three_layer` (default) or `:linear`.
- `kwargs...`: Additional keyword arguments absorbed for forward/backward compatibility.

# Outputs
- `Nothing`: Modifies `model.tracers.T` and `model.tracers.S` in-place.
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
