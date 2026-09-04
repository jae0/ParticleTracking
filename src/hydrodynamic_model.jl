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
- `closure`: Optional custom turbulence closure (e.g. `RiBasedVerticalDiffusivity()`).
- `tracers::Tuple`: Active tracer fields (default `(:T, :S)`).
- `free_surface`: Free surface representation (default `ImplicitFreeSurface()`).

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
    surface_wind_stress_x::Union{Real, AbstractMatrix, Function} = 0.0001,
    surface_wind_stress_y::Union{Real, AbstractMatrix, Function} = 0.0,
    surface_heat_flux::Union{Real, AbstractMatrix, Function} = 0.0,
    tidal_forcing::Union{Nothing, NamedTuple} = nothing,
    bottom_drag::Real = 1e-4,
    cd_drag::Real = 1e-3,
    ν::Real = 1e-2,
    κ::Real = 1e-2,
    closure = nothing,
    tracers::Tuple = (:T, :S),
    free_surface = ImplicitFreeSurface()
)
    # Surface kinematic boundary conditions for horizontal momentum
    u_top_bc = FluxBoundaryCondition(surface_wind_stress_x)
    v_top_bc = FluxBoundaryCondition(surface_wind_stress_y)

    u_bcs = FieldBoundaryConditions(top = u_top_bc)
    v_bcs = FieldBoundaryConditions(top = v_top_bc)

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
    T_bcs = FieldBoundaryConditions(top = T_top_bc)

    boundary_conditions = Dict{Symbol, Any}(:u => u_bcs, :v => v_bcs)
    if :T in tracers
        boundary_conditions[:T] = T_bcs
    end

    coriolis = FPlane(latitude = coriolis_latitude)
    buoyancy = SeawaterBuoyancy()
    active_closure = isnothing(closure) ? ScalarDiffusivity(ν = ν, κ = κ) : closure

    # Momentum forcing combining tidal oscillations and bottom boundary layer drag
    total_Fu(x, y, z, t, u, v) = begin
        tide_val = isnothing(tidal_forcing) ? 0.0 : tidal_forcing.u(x, y, z, t)
        speed = sqrt(u^2 + v^2)
        drag_val = -(bottom_drag + cd_drag * speed) * u
        return tide_val + drag_val
    end

    total_Fv(x, y, z, t, u, v) = begin
        tide_val = isnothing(tidal_forcing) ? 0.0 : tidal_forcing.v(x, y, z, t)
        speed = sqrt(u^2 + v^2)
        drag_val = -(bottom_drag + cd_drag * speed) * v
        return tide_val + drag_val
    end

    momentum_forcing = (
        u = Forcing(total_Fu, field_dependencies = (:u, :v)),
        v = Forcing(total_Fv, field_dependencies = (:u, :v))
    )

    model_kwargs = Dict{Symbol, Any}(
        :coriolis => coriolis,
        :buoyancy => buoyancy,
        :tracers => tracers,
        :boundary_conditions => NamedTuple(boundary_conditions),
        :forcing => momentum_forcing,
        :closure => active_closure,
        :free_surface => free_surface
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
    h_scale = dT_dz > 0.0 ? clamp((T_slope - T_cil_edge) / dT_dz, 50.0, 500.0) : 100.0

    function temp_profile(lon, lat, z)
        x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
        y_norm = clamp((lat - lat_min) / (lat_max - lat_min), 0.0, 1.0)
        T_surf = T0 + ΔT_cross * x_norm - ΔT_along * y_norm

        if stratification_type == :linear
            return T_surf + dT_dz * z
        end

        # Three-layer Northwest Atlantic / Scotian Shelf structure
        if z > -20.0
            # Surface mixed layer thermocline
            frac = (z + 20.0) / 20.0
            return clamp(T_cil_edge + frac * (T_surf - T_cil_edge), -1.5, 25.0)
        elseif z > -80.0
            # Cold Intermediate Layer (CIL): parabolic minimum at z = -50 m
            centre_frac = (z + 50.0) / 30.0
            T_cil = T_cil_min + ΔT_cil * centre_frac^2
            return clamp(T_cil, -1.5, T_surf)
        else
            # Deep Warm Slope Water: smooth exponential relaxation toward T_slope
            z_deep = abs(z) - 80.0
            T_deep = T_cil_edge + (T_slope - T_cil_edge) * (1.0 - exp(-z_deep / h_scale))
            return clamp(T_deep, min(T_cil_edge, T_slope), max(T_cil_edge, T_slope))
        end
    end

    sal_profile = if salinity isa Function
        salinity
    else
        function(lon, lat, z)
            x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
            S_cross = ΔS_cross * x_norm
            # Salinity increases with depth to guarantee static gravitational stability: ∂ρ/∂z ≤ 0
            S_depth = ΔS_deep * (1.0 - exp(-abs(z) / 150.0))
            clamp(S0 + S_cross + S_depth, 28.0, 36.5)
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
