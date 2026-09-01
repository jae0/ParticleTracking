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
        surface_temperature::Real = 14.0,
        salinity::Union{Real, Function} = 33.0,
        lon_range::Tuple{Real, Real} = (-71.0, -53.0),
        lat_range::Tuple{Real, Real} = (40.0, 48.5)
    )

Initialize realistic 3D thermal and haline stratification for the Scotian Shelf
and Gulf of St. Lawrence region.

# Mathematical Formulation

Temperature is a composite of three regimes along the vertical coordinate z ≤ 0:

**Surface mixed layer** (z > -20 m):
```math
T_{\\text{surf}}(\\lambda, \\phi) =
    T_0 + \\Delta T_{\\text{cross}} \\cdot x_{\\text{norm}}
       - \\Delta T_{\\text{along}} \\cdot y_{\\text{norm}}
```
where x_norm ∈ [0,1] increases offshore (westward/deeper) and y_norm ∈ [0,1]
increases northward. ΔT_cross ≈ 8°C (cold inshore → warm offshore slope);
ΔT_along ≈ 5°C (warm SW → cold NE, Labrador Current influence).

**Cold Intermediate Layer (CIL)** (-80 m < z ≤ -20 m):
```math
T_{\\text{CIL}}(\\lambda, \\phi, z) = T_{\\text{CIL,min}} +
    \\Delta T_{\\text{CIL}} \\cdot \\left(\\frac{z + 50}{30}\\right)^2
```
The CIL is a diagnostic feature of the Scotian Shelf, formed by winter
convection of surface water that persists as a lens of cold water (1–3°C)
between the warm surface layer and the deep slope water.

**Slope water** (z ≤ -80 m):
```math
T_{\\text{deep}}(z) = T_{\\text{CIL,min}} +
    \\frac{|z| - 80}{100} \\cdot (T_{\\text{slope}} - T_{\\text{CIL,min}})
```
where T_slope ≈ 8°C at z = -180 m (Warm Slope Water).

Salinity:
```math
S(\\lambda, \\phi, z) = S_0 + \\Delta S_{\\text{cross}} \\cdot x_{\\text{norm}}
    + \\frac{|z|}{100} \\cdot \\Delta S_{\\text{deep}}
```
Fresher coastal/northern waters (Gulf of St. Lawrence outflow), saltier
offshore slope water (North Atlantic Deep Water influence).

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
- `salinity::Union{Real, Function}`: Background salinity (PSU) or a function
  `(x, y, z) -> S` that returns practical salinity.
- `lon_range::Tuple`: Domain longitude bounds used for normalizing coordinates.
- `lat_range::Tuple`: Domain latitude bounds used for normalizing coordinates.

# Outputs
- `Nothing`: Modifies `model.tracers` in-place.
"""
function set_initial_stratification!(
    model;
    surface_temperature::Real = 15.0,
    temperature_gradient::Real = 0.01,
    salinity::Union{Real, Function} = 35.0,
    lon_range::Union{Nothing, Tuple} = nothing,
    lat_range::Union{Nothing, Tuple} = nothing
)

    lon_min, lon_max = Float64.(lon_range)
    lat_min, lat_max = Float64.(lat_range)
    T0               = Float64(surface_temperature)

    # Horizontal gradient scales (°C)
    ΔT_cross = 8.0   # cross-shelf: cold inshore/NE → warm offshore/SW
    ΔT_along = 5.0   # along-shelf: warm SW → cold NE (Labrador influence)

    # CIL (Cold Intermediate Layer) parameters
    T_cil_min = 1.5  # Coldest CIL core temperature (°C)
    T_slope   = 8.5  # Warm Slope Water below -180 m (°C)

    # Salinity gradient scales (PSU)
    S0        = Float64(salinity isa Real ? salinity : 33.0)
    ΔS_cross  = 2.0  # saltier offshore slope water
    ΔS_deep   = 1.5  # slightly saltier at depth

    function temp_profile(lon, lat, z)
        x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
        y_norm = clamp((lat - lat_min) / (lat_max - lat_min), 0.0, 1.0)

        # Surface temperature: warmer offshore/southwest, cooler inshore/northeast
        T_surf = T0 + ΔT_cross * x_norm - ΔT_along * y_norm

        if z > -20.0
            # Surface mixed layer: light linear cline to CIL top
            frac = (z + 20.0) / 20.0       # 0 at z=-20, 1 at z=0
            return clamp(T_cil_min + frac * (T_surf - T_cil_min), -1.5, 22.0)
        elseif z > -80.0
            # Cold Intermediate Layer: parabolic minimum centred at -50 m
            centre_frac = (z + 50.0) / 30.0
            T_cil = T_cil_min + 2.5 * centre_frac^2
            return clamp(T_cil, -1.5, T_surf)
        else
            # Warm Slope Water: linear recovery toward T_slope below CIL
            depth_factor = (abs(z) - 80.0) / 100.0
            T_deep = T_cil_min + depth_factor * (T_slope - T_cil_min)
            return clamp(T_deep, T_cil_min, T_slope)
        end
    end

    sal_profile = if salinity isa Function
        salinity
    else
        function(lon, lat, z)
            x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
            S_cross = ΔS_cross * x_norm
            S_depth = ΔS_deep * min(abs(z), 200.0) / 200.0
            clamp(S0 + S_cross + S_depth, 28.0, 36.0)
        end
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
