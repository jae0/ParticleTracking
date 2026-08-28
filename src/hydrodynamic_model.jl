"""
    hydrodynamic_model.jl

Configuration and initialization of hydrostatic free-surface ocean models
in Oceananigans.jl for coastal shelf domains.
"""

using Oceananigans
using Oceananigans.Units

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
+ \\nabla \\cdot (\\nu \\nabla \\boldsymbol{u}_h)
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
and \$\\nu, \\kappa\$ represent momentum and tracer eddy diffusivities.

# Inputs
- `grid`: Underlying computational grid or `ImmersedBoundaryGrid`.
- `coriolis_latitude::Real`: Reference latitude for Coriolis `FPlane` in degrees North.
- `surface_wind_stress_x`: Zonal kinematic surface momentum flux (scalar, 2D array, or function).
- `surface_wind_stress_y`: Meridional kinematic surface momentum flux.
- `ν::Real`: Kinematic eddy viscosity (\$m^2 s^{-1}\$).
- `κ::Real`: Tracer eddy diffusivity (\$m^2 s^{-1}\$).
- `tracers::Tuple`: Active tracer fields (default `(:T, :S)`).
- `free_surface`: Free surface representation (default `ImplicitFreeSurface()`).

# Outputs
- `HydrostaticFreeSurfaceModel`: Configured Oceananigans model instance.

# References
- Ramadhan, A., Marshall, J., Hill, C., Campin, J. M., Bischoff, T., & Wagner,
  G. L. (2020). Oceananigans.jl: Fast and friendly geophysical fluid dynamics on
  GPUs. *Journal of Open Source Software*, 5(53), 2018. DOI: 10.21105/joss.02018
- Marshall, J., Adcroft, A., Hill, C., Perelman, L., & Heisey, C. (1997).
  A finite-volume, incompressible Navier Stokes model for studies of the ocean on
  parallel computers. *Journal of Geophysical Research: Oceans*, 102(C3), 5753-5766.
- Vallis, G. K. (2017). *Atmospheric and Oceanic Fluid Dynamics: Fundamentals
  and Large-Scale Circulation*. 2nd Edition. Cambridge University Press.
"""
function build_hydrodynamic_model(
    grid;
    coriolis_latitude::Real = 45.0,
    surface_wind_stress_x::Union{Real, AbstractMatrix, Function} = 0.1,
    surface_wind_stress_y::Union{Real, AbstractMatrix, Function} = 0.0,
    tidal_forcing::Union{Nothing, NamedTuple} = nothing,
    ν::Real = 1e-2,
    κ::Real = 1e-2,
    tracers::Tuple = (:T, :S),
    free_surface = ImplicitFreeSurface()
)
    # Surface kinematic boundary conditions for horizontal momentum
    u_top_bc = FluxBoundaryCondition(surface_wind_stress_x)
    v_top_bc = FluxBoundaryCondition(surface_wind_stress_y)

    u_bcs = FieldBoundaryConditions(top = u_top_bc)
    v_bcs = FieldBoundaryConditions(top = v_top_bc)

    coriolis = FPlane(latitude = coriolis_latitude)
    buoyancy = SeawaterBuoyancy()
    closure  = ScalarDiffusivity(ν = ν, κ = κ)

    model_kwargs = Dict{Symbol, Any}(
        :coriolis => coriolis,
        :buoyancy => buoyancy,
        :tracers => tracers,
        :boundary_conditions => (u = u_bcs, v = v_bcs),
        :closure => closure,
        :free_surface => free_surface
    )

    if !isnothing(tidal_forcing)
        model_kwargs[:forcing] = tidal_forcing
    end

    model = HydrostaticFreeSurfaceModel(grid; model_kwargs...)

    return model
end

"""
    set_initial_stratification!(
        model;
        surface_temperature::Real = 15.0,
        temperature_gradient::Real = 0.01,
        salinity::Union{Real, Function} = 35.0
    )

Initialize thermal and haline stratification across the vertical water column.

# Mathematical Formulation
Linear temperature profile with depth \$z \\le 0\$:
```math
T(x, y, z) = T_{\\text{surface}} + \\left(\\frac{dT}{dz}\\right) z
```
Salinity is set uniformly or via a depth-dependent profile:
```math
S(x, y, z) = S_0
```

# Inputs
- `model::HydrostaticFreeSurfaceModel`: Model instance to initialize.
- `surface_temperature::Real`: Sea surface temperature in °C.
- `temperature_gradient::Real`: Vertical thermal gradient \$dT/dz\$ in °C/m.
- `salinity::Union{Real, Function}`: Background practical salinity in PSU.

# Outputs
- `Nothing`: Modifies `model.tracers` in-place.
"""
function set_initial_stratification!(
    model;
    surface_temperature::Real = 15.0,
    temperature_gradient::Real = 0.01,
    salinity::Union{Real, Function} = 35.0
)
    temp_profile(x, y, z) = surface_temperature + temperature_gradient * z

    if salinity isa Real
        sal_profile = salinity
    else
        sal_profile = salinity
    end

    set!(model, T = temp_profile, S = sal_profile)
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
