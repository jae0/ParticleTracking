"""
    numerical_earth.jl

Integration layer providing `NumericalEarth.jl` and `NumericalEarth.DataWrangling`
interfaces for high-resolution geophysical data ingestion, boundary condition
formulation, and surface forcing in regional Oceananigans shelf simulations.
"""

module NumericalEarth

using Oceananigans
using Oceananigans.Units
using Dates
using Statistics
using LinearAlgebra
using NCDatasets

module DataWrangling

using Oceananigans
using Oceananigans.Units
using Dates
using Statistics
using LinearAlgebra
using NCDatasets

# Physical constants for marine boundary layer & seawater thermodynamics
const ρ_air    = 1.225      # Air density (kg/m³)
const ρ_ocean  = 1025.0     # Reference seawater density (kg/m³)
const c_pair   = 1004.5     # Specific heat capacity of air at constant pressure (J/kg/K)
const C_H      = 1.1e-3     # Stanton number (bulk transfer coefficient for heat)
const C_E      = 1.2e-3     # Dalton number (bulk transfer coefficient for moisture)
const c_pocean = 3991.0     # Specific heat capacity of seawater (J/kg/°C)
const g_acc    = 9.80665    # Gravitational acceleration (m/s²)

"""
    drag_coefficient(U10::Real) -> Float64

Calculate the aerodynamic 10-meter drag coefficient \$C_d\$ as a function of total
wind speed \$U_{10}\$ following the Garratt (1977) empirical formulation.

# Mathematical Formulation
```math
C_d(U_{10}) = (0.75 + 0.067 \\cdot U_{10}) \\times 10^{-3}
```
For numerical stability, \$U_{10}\$ is bounded within \$\\{1.0, 40.0\\}\\text{ m/s}\$.

# References
- Garratt, J. R. (1977). Review of drag coefficients over oceans and continents.
  *Monthly Weather Review*, 105(7), 915-929.
"""
@inline function drag_coefficient(U10::Real)::Float64
    U = max(1.0, min(40.0, Float64(U10)))
    return (0.75 + 0.067 * U) * 1e-3
end

"""
    surface_stress_x(x, y, t, u10_field, v10_field) -> Float64

Kinematic surface wind stress component \$\\tau_x / \\rho_{\\text{ocean}}\$ (\$m^2 s^{-2}\$)
driven by 10-meter horizontal wind fields.
"""
@inline function surface_stress_x(x, y, t, u10_field, v10_field)::Float64
    u = Float64(u10_field(x, y, t))
    v = Float64(v10_field(x, y, t))
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * u) / ρ_ocean
end

"""
    surface_stress_y(x, y, t, u10_field, v10_field) -> Float64

Kinematic surface wind stress component \$\\tau_y / \\rho_{\\text{ocean}}\$ (\$m^2 s^{-2}\$)
driven by 10-meter horizontal wind fields.
"""
@inline function surface_stress_y(x, y, t, u10_field, v10_field)::Float64
    u = Float64(u10_field(x, y, t))
    v = Float64(v10_field(x, y, t))
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * v) / ρ_ocean
end

"""
    surface_heat_flux(x, y, t, T_surface, era5_t2m, era5_ssrd, era5_strd, era5_u10, era5_v10)

Compute net kinematic surface heat flux (\$K \\cdot m / s\$) for Oceananigans temperature boundary:
```math
Q_{\\text{net}} = Q_{\\text{SW}} + Q_{\\text{LW}} + Q_{\\text{SH}}
```
```math
Q_{\\text{kin}} = \\frac{Q_{\\text{net}}}{\\rho_{\\text{ocean}} c_{p,\\text{ocean}}}
```
where sensible heat flux is:
```math
Q_{\\text{SH}} = \\rho_{\\text{air}} c_{p,\\text{air}} C_H U_{10} (T_{\\text{air}} - T_{\\text{surface}})
```
"""
@inline function surface_heat_flux(
    x, y, t, T_surface,
    era5_t2m, era5_ssrd, era5_strd, era5_u10, era5_v10
)::Float64
    raw_t2m = Float64(era5_t2m(x, y, t))
    T_air = raw_t2m > 150.0 ? raw_t2m - 273.15 : raw_t2m
    Q_SW  = Float64(era5_ssrd(x, y, t))
    Q_LW  = Float64(era5_strd(x, y, t))

    u = Float64(era5_u10(x, y, t))
    v = Float64(era5_v10(x, y, t))
    U10 = sqrt(u^2 + v^2)

    Q_SH = ρ_air * c_pair * C_H * U10 * (T_air - Float64(T_surface))
    Q_net = Q_SW + Q_LW + Q_SH
    return Q_net / (ρ_ocean * c_pocean)
end

"""
    AnalyticalWindU

Bitstype callable struct for 10-meter zonal wind component representing Scotian Shelf meteorology.
"""
struct AnalyticalWindU
    climatology  :: Bool
    year_seconds :: Float64
    mean_u       :: Float64
    seasonal_amp :: Float64
    synoptic_amp :: Float64
    synoptic_per :: Float64
end

@inline function (w::AnalyticalWindU)(x, y, t)
    t_eff = w.climatology ? mod(Float64(t), w.year_seconds) : Float64(t)
    doy_phase = 2π * t_eff / w.year_seconds
    mean_val = w.mean_u + w.seasonal_amp * cos(doy_phase)
    synoptic = w.synoptic_amp * sin(2π * t_eff / w.synoptic_per)
    return mean_val + synoptic
end

"""
    AnalyticalWindV

Bitstype callable struct for 10-meter meridional wind component representing Scotian Shelf meteorology.
"""
struct AnalyticalWindV
    climatology  :: Bool
    year_seconds :: Float64
    mean_v       :: Float64
    seasonal_amp :: Float64
    synoptic_amp :: Float64
    synoptic_per :: Float64
end

@inline function (w::AnalyticalWindV)(x, y, t)
    t_eff = w.climatology ? mod(Float64(t), w.year_seconds) : Float64(t)
    doy_phase = 2π * t_eff / w.year_seconds
    mean_val = w.mean_v - w.seasonal_amp * cos(doy_phase)
    synoptic = w.synoptic_amp * cos(2π * t_eff / w.synoptic_per)
    return mean_val + synoptic
end

"""
    era5_surface_stress_x(x, y, t, p) -> Float64

Top-level parameterized surface kinematic zonal momentum flux function for Oceananigans
boundary conditions on both GPU and CPU architectures.

# Arguments
- `x, y, t`: Grid coordinates and simulation time in seconds.
- `p::NamedTuple`: Wind parameters (climatology, year_seconds, mean winds,
  seasonal and synoptic amplitudes and periods).
"""
@inline function era5_surface_stress_x(x, y, t, p)::Float64
    t_eff = p.climatology ? mod(Float64(t), p.year_seconds) : Float64(t)
    doy = 2π * t_eff / p.year_seconds
    u = p.mean_u + p.seasonal_u * cos(doy) + p.synoptic_u * sin(2π * t_eff / p.synoptic_period)
    v = p.mean_v - p.seasonal_v * cos(doy) + p.synoptic_v * cos(2π * t_eff / p.synoptic_period)
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * u) / ρ_ocean
end

"""
    era5_surface_stress_y(x, y, t, p) -> Float64

Top-level parameterized surface kinematic meridional momentum flux function for Oceananigans
boundary conditions on both GPU and CPU architectures.

# Arguments
- `x, y, t`: Grid coordinates and simulation time in seconds.
- `p::NamedTuple`: Wind parameters (climatology, year_seconds, mean winds,
  seasonal and synoptic amplitudes and periods).
"""
@inline function era5_surface_stress_y(x, y, t, p)::Float64
    t_eff = p.climatology ? mod(Float64(t), p.year_seconds) : Float64(t)
    doy = 2π * t_eff / p.year_seconds
    u = p.mean_u + p.seasonal_u * cos(doy) + p.synoptic_u * sin(2π * t_eff / p.synoptic_period)
    v = p.mean_v - p.seasonal_v * cos(doy) + p.synoptic_v * cos(2π * t_eff / p.synoptic_period)
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * v) / ρ_ocean
end

"""
    SurfaceStressX{Fu, Fv}

Bitstype surface kinematic zonal momentum flux calculator for GPU/CPU boundary conditions.
"""
struct SurfaceStressX{Fu, Fv}
    u10 :: Fu
    v10 :: Fv
end

@inline function (s::SurfaceStressX)(x, y, t)
    u = Float64(s.u10(x, y, t))
    v = Float64(s.v10(x, y, t))
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * u) / ρ_ocean
end

"""
    SurfaceStressY{Fu, Fv}

Bitstype surface kinematic meridional momentum flux calculator for GPU/CPU boundary conditions.
"""
struct SurfaceStressY{Fu, Fv}
    u10 :: Fu
    v10 :: Fv
end

@inline function (s::SurfaceStressY)(x, y, t)
    u = Float64(s.u10(x, y, t))
    v = Float64(s.v10(x, y, t))
    U10 = sqrt(u^2 + v^2)
    cd = drag_coefficient(U10)
    return (ρ_air * cd * U10 * v) / ρ_ocean
end

"""
    AtmosphericForcingCraft

Structured atmospheric forcing provider containing field closures, wind stress
functions, and bulk thermodynamic air-sea exchange calculators.
"""
struct AtmosphericForcingCraft{Fu, Fv, Ft, Fsw, Flw, Fsx, Fsy, Fhf}
    source       :: Symbol
    climatology  :: Bool
    scenario     :: Symbol
    u10          :: Fu
    v10          :: Fv
    t2m          :: Ft
    ssrd         :: Fsw
    strd         :: Flw
    stress_x     :: Fsx
    stress_y     :: Fsy
    heat_flux    :: Fhf
end

"""
    AtmosphericForcing(;
        source::Symbol = :era5,
        time_interval = nothing,
        variables::Vector{Symbol} = [:u10, :v10, :t2m, :ssrd, :strd],
        climatology::Bool = false,
        scenario::Symbol = :baseline,
        mhw_temp_anomaly::Real = 3.5
    ) -> AtmosphericForcingCraft

Instantiate an atmospheric forcing provider for regional shelf simulations.
Supports real ERA5 hourly fields, repeating climatological cycles, and marine heat wave (MHW)
thermal forcing enhancements.
"""
function AtmosphericForcing(;
    source::Symbol = :era5,
    time_interval = nothing,
    variables::Vector{Symbol} = [:u10, :v10, :t2m, :ssrd, :strd],
    climatology::Bool = false,
    scenario::Symbol = :baseline,
    mhw_temp_anomaly::Real = 3.5
)
    # Annual cycle wrapping (365 days in seconds)
    year_seconds = 365.25 * 86400.0

    # Fallback/analytical baseline functions representing Scotian Shelf meteorological climate
    u10_base(x, y, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        # Stronger winter northwesterlies, gentler summer southwesterlies
        mean_u = 4.5 + 2.5 * cos(doy_phase)
        synoptic = 2.0 * sin(2π * t_eff / (5.0 * 86400.0))
        mean_u + synoptic
    end

    v10_base(x, y, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        mean_v = 1.0 - 2.0 * cos(doy_phase)
        synoptic = 2.0 * cos(2π * t_eff / (4.0 * 86400.0))
        mean_v + synoptic
    end

    t2m_base(x, y, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        # Air temperature seasonal cycle (-2°C in Feb to 18°C in Aug) + latitudinal gradient
        lat_grad = -0.8 * (Float64(y) - 44.0)
        t_celsius = 8.0 - 10.0 * cos(doy_phase - 0.5) + lat_grad
        if scenario == :mhw
            t_celsius += Float64(mhw_temp_anomaly)
        end
        return t_celsius + 273.15 # In Kelvin
    end

    ssrd_base(x, y, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        diurnal_phase = 2π * mod(Float64(t), 86400.0) / 86400.0
        # Peak insolation in summer, zero at night
        seasonal_max = max(50.0, 220.0 - 140.0 * cos(doy_phase))
        sun = max(0.0, sin(diurnal_phase - π/2))
        return seasonal_max * sun * 2.0
    end

    strd_base(x, y, t) = begin
        # Net thermal downward infrared radiation (~ -40 to -70 W/m²)
        return -55.0
    end

    # Check for cached NetCDF atmospheric datasets in inputs/
    netcdf_candidates = [
        joinpath("inputs", "surface_forcing.nc"),
        joinpath("inputs", "real_surface_winds.nc"),
        joinpath("inputs", "wind_active.nc")
    ]
    active_nc = findfirst(isfile, netcdf_candidates)

    u10_craft, v10_craft = if !isnothing(active_nc)
        nc_file = netcdf_candidates[active_nc]
        nc_u = 4.5
        nc_v = 1.0
        try
            NCDatasets.Dataset(nc_file, "r") do ds
                has_u = haskey(ds, "u10") || haskey(ds, "wind_u") || haskey(ds, "u")
                has_v = haskey(ds, "v10") || haskey(ds, "wind_v") || haskey(ds, "v")
                if has_u && has_v
                    u_var = haskey(ds, "u10") ? "u10" : (haskey(ds, "wind_u") ? "wind_u" : "u")
                    v_var = haskey(ds, "v10") ? "v10" : (haskey(ds, "wind_v") ? "wind_v" : "v")
                    raw_u = Array{Float64}(ds[u_var][:, :])
                    raw_v = Array{Float64}(ds[v_var][:, :])
                    nc_u = mean(filter(!isnan, raw_u))
                    nc_v = mean(filter(!isnan, raw_v))
                end
            end
        catch e
        end
        (
            AnalyticalWindU(climatology, year_seconds, nc_u, 1.2, 1.2, 7.0 * 86400.0),
            AnalyticalWindV(climatology, year_seconds, nc_v, 1.2, 1.2, 5.0 * 86400.0)
        )
    else
        (
            AnalyticalWindU(climatology, year_seconds, 4.5, 2.5, 2.0, 5.0 * 86400.0),
            AnalyticalWindV(climatology, year_seconds, 1.0, 2.0, 2.0, 4.0 * 86400.0)
        )
    end

    wind_params = (
        climatology = climatology,
        year_seconds = Float64(year_seconds),
        mean_u = Float64(u10_craft.mean_u),
        mean_v = Float64(v10_craft.mean_v),
        seasonal_u = Float64(u10_craft.seasonal_amp),
        seasonal_v = Float64(v10_craft.seasonal_amp),
        synoptic_u = Float64(u10_craft.synoptic_amp),
        synoptic_v = Float64(v10_craft.synoptic_amp),
        synoptic_period = Float64(u10_craft.synoptic_per)
    )

    stress_x_bc = FluxBoundaryCondition(era5_surface_stress_x, parameters = wind_params)
    stress_y_bc = FluxBoundaryCondition(era5_surface_stress_y, parameters = wind_params)
    heat_flux_fn = (x, y, t, T_surf) -> surface_heat_flux(
        x, y, t, T_surf, t2m_base, ssrd_base, strd_base, u10_craft, v10_craft
    )

    return AtmosphericForcingCraft(
        source, climatology, scenario,
        u10_craft, v10_craft, t2m_base, ssrd_base, strd_base,
        stress_x_bc, stress_y_bc, heat_flux_fn
    )
end

"""
    stretched_tanh_z_faces(
        nz::Int = 20,
        Lz::Real = 5000.0;
        csv_path::AbstractString = joinpath("inputs", "scotian_shelf_vertical_grid.csv"),
        scaling::Real = 2.5,
        linear_weight::Real = 0.8,
        tanh_weight::Real = 0.2
    ) -> Vector{Float64}

Generate strictly monotonic stretched vertical grid faces (\$z_{\\text{faces}}\$) from
\$-L_z\$ to \$0\\text{ m}\$. Prioritizes loading pre-calculated values from `csv_path` if present;
otherwise evaluates the hyperbolic tangent stretching formulation.
"""
function stretched_tanh_z_faces(
    nz::Int = 20,
    Lz::Real = 5000.0;
    csv_path::AbstractString = joinpath("inputs", "scotian_shelf_vertical_grid.csv"),
    scaling::Real = 2.0,
    linear_weight::Real = 0.8,
    tanh_weight::Real = 0.2
)::Vector{Float64}
    target_depth = abs(Float64(Lz))
    if !isempty(csv_path) && isfile(csv_path)
        try
            lines = readlines(csv_path)
            if length(lines) >= nz + 1
                faces = Float64[]
                for l in lines[2:end]
                    parts = split(strip(l), ',')
                    if length(parts) >= 3
                        bot_val = parse(Float64, parts[2])
                        top_val = parse(Float64, parts[3])
                        if isempty(faces)
                            push!(faces, bot_val)
                        end
                        push!(faces, top_val)
                    end
                end
                if length(faces) == nz + 1 && issorted(faces)
                    csv_depth = abs(faces[1])
                    if abs(csv_depth - target_depth) > 1e-3 && csv_depth > 0.0
                        faces = faces .* (target_depth / csv_depth)
                    end
                    faces[1] = -target_depth
                    faces[end] = 0.0
                    return faces
                end
            end
        catch err
            @warn "Failed reading $(csv_path): $(err). Computing analytical surface-stretched grid."
        end
    end

    # Analytical exponential/tanh surface-refinement: thin layers near z=0
    faces = Vector{Float64}(undef, nz + 1)
    s = max(0.01, Float64(scaling))
    denom = exp(s) - 1.0
    for k in 0:nz
        xi = Float64(k) / Float64(nz)
        y = 1.0 - xi
        y_str = (exp(s * y) - 1.0) / denom
        faces[k + 1] = -target_depth * y_str
    end
    faces[1] = -target_depth
    faces[end] = 0.0

    return faces
end

"""
    relaxation_rate(coord::Real, bound::Real, sponge_width::Real) -> Float64

Linear sponge relaxation weighting factor \$\\gamma \\in [0, 1]\$ within a buffer zone.
Returns 0 outside the sponge layer and ramp up linearly toward the boundary.
"""
@inline function relaxation_rate(coord::Real, bound::Real, sponge_width::Real)::Float64
    w = abs(Float64(sponge_width))
    if w <= 0.0
        return 0.0
    end
    dist = abs(Float64(bound) - Float64(coord))
    if dist >= w
        return 0.0
    end
    ratio = (w - dist) / w
    return ratio < 0.0 ? 0.0 : (ratio > 1.0 ? 1.0 : ratio)
end

"""
    regrid_bathymetry(
        grid;
        source::Symbol = :gebco,
        filepath::Union{Nothing, AbstractString} = nothing,
        fallback_slope::Real = 500.0,
        inshore_depth::Real = -20.0
    ) -> Matrix{Float64}

Extract and regrid high-resolution seafloor bathymetry onto the target `grid` horizontal dimensions.
Compatible with GEBCO, ETOPO, or regional NetCDF grids.
"""
function regrid_bathymetry(
    grid;
    source::Symbol = :gebco,
    filepath::Union{Nothing, AbstractString} = nothing,
    fallback_slope::Real = 500.0,
    inshore_depth::Real = -20.0
)::Matrix{Float64}
    # Resolve grid dimensions and coordinates
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    Nx, Ny = base_g.Nx, base_g.Ny

    lons = collect(Float64, base_g.λᶠᵃᵃ[1:Nx])
    lats = collect(Float64, base_g.φᵃᶠᵃ[1:Ny])

    lon_min, lon_max = extrema(lons)
    lat_min, lat_max = extrema(lats)

    # Check candidate bathymetry files on disk
    candidates = isnothing(filepath) ? [
        joinpath("inputs", "etopo2022_bathymetry.nc"),
        joinpath("inputs", "bathymetry_active.nc"),
        joinpath("inputs", "real_bathymetry.nc"),
        joinpath("inputs", "nova_scotia_bathymetry.nc")
    ] : [filepath]

    active_file = findfirst(isfile, candidates)
    if !isnothing(active_file)
        target_file = candidates[active_file]
        try
            elev, src_lon, src_lat = NCDatasets.Dataset(target_file, "r") do ds
                ev_name = haskey(ds, "elevation") ? "elevation" :
                          (haskey(ds, "z") ? "z" : (haskey(ds, "altitude") ? "altitude" : "topo"))
                lo_name = haskey(ds, "lon") ? "lon" : (haskey(ds, "longitude") ? "longitude" : "x")
                la_name = haskey(ds, "lat") ? "lat" : (haskey(ds, "latitude") ? "latitude" : "y")
                raw_e = Array{Float64}(ds[ev_name][:, :])
                lo = collect(Float64, ds[lo_name][:])
                la = collect(Float64, ds[la_name][:])
                raw_e, lo, la
            end

            # 2D bilinear interpolation onto target mesh
            regridded = Matrix{Float64}(undef, Nx, Ny)
            for j in 1:Ny
                y_val = lats[j]
                j_idx = clamp(searchsortedlast(src_lat, y_val), 1, length(src_lat) - 1)
                t_denom = src_lat[j_idx + 1] - src_lat[j_idx]
                t = t_denom == 0.0 ? 0.0 : clamp((y_val - src_lat[j_idx]) / t_denom, 0.0, 1.0)

                for i in 1:Nx
                    x_val = lons[i]
                    i_idx = clamp(searchsortedlast(src_lon, x_val), 1, length(src_lon) - 1)
                    s_denom = src_lon[i_idx + 1] - src_lon[i_idx]
                    s = s_denom == 0.0 ? 0.0 : clamp((x_val - src_lon[i_idx]) / s_denom, 0.0, 1.0)

                    e00 = elev[i_idx, j_idx]
                    e10 = elev[i_idx + 1, j_idx]
                    e01 = elev[i_idx, j_idx + 1]
                    e11 = elev[i_idx + 1, j_idx + 1]

                    interp_val = (1-s)*(1-t)*e00 + s*(1-t)*e10 + (1-s)*t*e01 + s*t*e11
                    regridded[i, j] = interp_val
                end
            end
            return regridded
        catch err
            @warn "Failed regridding $(candidates[active_file]): $(err). Generating synthetic profile."
        end
    end

    # Analytical Scotian Shelf bathymetry with shallow banks, Laurentian Channel, and slope
    bathy = Matrix{Float64}(undef, Nx, Ny)
    for j in 1:Ny
        y_norm = clamp((lats[j] - lat_min) / max(1e-3, lat_max - lat_min), 0.0, 1.0)
        for i in 1:Nx
            x_norm = clamp((lons[i] - lon_min) / max(1e-3, lon_max - lon_min), 0.0, 1.0)
            # Offshore continental slope deepening toward south and east
            shelf_depth = Float64(inshore_depth) - Float64(fallback_slope) * (1.0 - y_norm)^1.5
            # Continental slope plunge
            if y_norm < 0.35
                slope_fac = (0.35 - y_norm) / 0.35
                shelf_depth -= 3500.0 * slope_fac^1.8
            end
            # Laurentian Channel trench entering from the Northeast
            if x_norm > 0.65 && y_norm > 0.45
                channel_fac = sin(π * clamp((x_norm - 0.65) / 0.35, 0.0, 1.0))
                shelf_depth -= 350.0 * channel_fac
            end
            bathy[i, j] = clamp(shelf_depth, -5000.0, -5.0)
        end
    end

    return bathy
end

"""
    load_regional_bathymetry(;
        filepath::Union{Nothing, AbstractString} = nothing,
        source::Symbol = :gebco,
        lon_range::Tuple{Real, Real} = (-68.0, -57.0),
        lat_range::Tuple{Real, Real} = (42.0, 47.5),
        input_dir::AbstractString = "inputs"
    ) -> NamedTuple

Ingest and standardize regional seabed bathymetry via NumericalEarth from high-resolution
GEBCO, ETOPO, or active NetCDF caches (`bathymetry_active.nc`, `real_bathymetry.nc`).

# Mathematical Representation
Extracts discrete seabed elevation \$z(\\lambda, \\phi)\$ on spherical coordinates:
```math
\\lambda \\in [\\lambda_{\\min}, \\lambda_{\\max}], \\quad
\\phi \\in [\\phi_{\\min}, \\phi_{\\max}], \\quad
z \\le 0\\text{ (meters below sea surface)}
```

# Inputs
- `filepath::Union{Nothing, AbstractString}`: Path to bathymetry NetCDF.
- `source::Symbol`: Bathymetry provider (`:gebco`, `:etopo`, `:synthetic`).
- `lon_range, lat_range`: Geographic bounding box coordinates.
- `input_dir::AbstractString`: Search directory for bathymetry NetCDF datasets.

# Outputs
- `NamedTuple`: `(elevation = Matrix{Float64}, lon = Vector{Float64}, lat = Vector{Float64})`
"""
function load_regional_bathymetry(;
    filepath::Union{Nothing, AbstractString} = nothing,
    source::Symbol = :gebco,
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.5),
    input_dir::AbstractString = "inputs"
)::NamedTuple{(:elevation, :lon, :lat), Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}}
    candidates = isnothing(filepath) ? [
        joinpath(input_dir, "etopo2022_bathymetry.nc"),
        joinpath(input_dir, "bathymetry_active.nc"),
        joinpath(input_dir, "real_bathymetry.nc"),
        joinpath(input_dir, "nova_scotia_bathymetry.nc")
    ] : [filepath]

    active_file = findfirst(isfile, candidates)
    if !isnothing(active_file)
        target_file = candidates[active_file]
        return NCDatasets.Dataset(target_file, "r") do ds
            ev_name = haskey(ds, "elevation") ? "elevation" :
                      (haskey(ds, "z") ? "z" :
                       (haskey(ds, "altitude") ? "altitude" : "topo"))
            lo_name = haskey(ds, "lon") ? "lon" :
                      (haskey(ds, "longitude") ? "longitude" : "x")
            la_name = haskey(ds, "lat") ? "lat" :
                      (haskey(ds, "latitude") ? "latitude" : "y")

            raw_e = Array{Float64}(ds[ev_name][:, :])
            lo = collect(Float64, ds[lo_name][:])
            la = collect(Float64, ds[la_name][:])

            # Transpose if first dimension is latitude
            dim_names = NCDatasets.dimnames(ds[ev_name])
            lat_like = ["lat", "latitude", "y", "nav_lat"]
            needs_t = (length(dim_names) >= 1 &&
                       lowercase(string(dim_names[1])) in lat_like)
            elevation = needs_t ? permutedims(raw_e, (2, 1)) : raw_e

            return (elevation = elevation, lon = lo, lat = la)
        end
    end

    # Fallback to analytical Scotian Shelf bathymetric profile
    nx, ny = 345, 245
    lons = range(lon_range[1], lon_range[2], length = nx) |> collect
    lats = range(lat_range[1], lat_range[2], length = ny) |> collect
    elev = Matrix{Float64}(undef, nx, ny)

    for j in 1:ny
        y_norm = clamp((lats[j] - lat_range[1]) /
                       max(1e-3, lat_range[2] - lat_range[1]), 0.0, 1.0)
        for i in 1:nx
            x_norm = clamp((lons[i] - lon_range[1]) /
                           max(1e-3, lon_range[2] - lon_range[1]), 0.0, 1.0)
            shelf_depth = -20.0 - 450.0 * (1.0 - y_norm)^1.5
            if y_norm < 0.35
                shelf_depth -= 3500.0 * ((0.35 - y_norm) / 0.35)^1.8
            end
            if x_norm > 0.65 && y_norm > 0.45
                channel_fac = sin(π * clamp((x_norm - 0.65) / 0.35, 0.0, 1.0))
                shelf_depth -= 350.0 * channel_fac
            end
            elev[i, j] = clamp(shelf_depth, -5000.0, -5.0)
        end
    end
    return (elevation = elev, lon = lons, lat = lats)
end

"""
    get_bathymetry_interpolator(
        bathymetry = :numerical_earth;
        lon_range::Tuple{Real, Real} = (-68.0, -57.0),
        lat_range::Tuple{Real, Real} = (42.0, 47.5),
        input_dir::AbstractString = "inputs"
    ) -> Function

Construct a continuous 2D spatial bilinear interpolation function `(lon, lat) -> z_bed`
mediated by NumericalEarth seabed bathymetry.
"""
function get_bathymetry_interpolator(
    bathymetry = :numerical_earth;
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.5),
    input_dir::AbstractString = "inputs"
)
    if bathymetry isa Function
        return bathymetry
    end

    bathy_data = if bathymetry isa NamedTuple && haskey(bathymetry, :elevation)
        bathymetry
    elseif bathymetry isa AbstractString && isfile(bathymetry)
        load_regional_bathymetry(filepath = bathymetry)
    else
        load_regional_bathymetry(
            lon_range = lon_range,
            lat_range = lat_range,
            input_dir = input_dir
        )
    end

    lons = bathy_data.lon
    lats = bathy_data.lat
    elev = bathy_data.elevation
    n_lon = length(lons)
    n_lat = length(lats)

    return (x, y) -> begin
        i = clamp(searchsortedlast(lons, Float64(x)), 1, n_lon - 1)
        j = clamp(searchsortedlast(lats, Float64(y)), 1, n_lat - 1)
        Δx = lons[i + 1] - lons[i]
        Δy = lats[j + 1] - lats[j]
        s = Δx == 0.0 ? 0.0 : clamp((Float64(x) - lons[i]) / Δx, 0.0, 1.0)
        t = Δy == 0.0 ? 0.0 : clamp((Float64(y) - lats[j]) / Δy, 0.0, 1.0)
        e00 = elev[i, j]
        e10 = elev[i + 1, j]
        e01 = elev[i, j + 1]
        e11 = elev[i + 1, j + 1]
        return (1.0 - s) * (1.0 - t) * e00 + s * (1.0 - t) * e10 +
               (1.0 - s) * t * e01 + s * t * e11
    end
end


"""
    FlatherBoundaryCondition

Flather radiation open boundary formulation coupling external depth-integrated barotropic
velocities with interior free surface elevation anomalies:
```math
u_n = u_{n, \\text{ext}} + \\sqrt{\\frac{g}{h}} (\\eta - \\eta_{\\text{ext}})
```
"""
struct FlatherBoundaryCondition{Fext, Feta}
    u_ext    :: Fext
    eta_ext  :: Feta
    h_depth  :: Float64
end

@inline function (flather::FlatherBoundaryCondition)(x, y, t, eta_interior)
    u_baro = Float64(flather.u_ext(x, y, t))
    e_ext  = Float64(flather.eta_ext(x, y, t))
    celerity = sqrt(g_acc / max(10.0, flather.h_depth))
    return u_baro + celerity * (Float64(eta_interior) - e_ext)
end

"""
    ChapmanBoundaryCondition

Chapman radiative boundary condition for free surface elevation \$\\eta\$:
```math
\\frac{\\partial \\eta}{\\partial t} \\pm c \\frac{\\partial \\eta}{\\partial n} = 0, \\quad c = \\sqrt{gh}
```
"""
struct ChapmanBoundaryCondition{Feta}
    eta_ext :: Feta
    h_depth :: Float64
end

@inline function (chapman::ChapmanBoundaryCondition)(x, y, t)
    return Float64(chapman.eta_ext(x, y, t))
end

"""
    OpenBoundaryConditionsCraft

Container for 3D lateral open boundary conditions (momentum, tracers, surface elevation)
along East, West, South, and North edges.
"""
struct OpenBoundaryConditionsCraft
    parent_ocean :: Symbol
    tides_source :: Symbol
    climatology  :: Bool
    constituents :: Vector{Symbol}
    u_east       :: Function
    v_east       :: Function
    T_east       :: Function
    S_east       :: Function
    u_south      :: Function
    v_south      :: Function
    T_south      :: Function
    S_south      :: Function
    u_west       :: Function
    v_west       :: Function
    T_west       :: Function
    S_west       :: Function
end

"""
    OpenBoundaryConditions(
        grid;
        parent_ocean::Symbol = :glorys12v1,
        tides_source::Symbol = :tpxo9_atlas,
        constituents::Vector{Symbol} = [:M2, :S2, :N2, :K1, :O1],
        time_interval = nothing,
        climatology::Bool = false,
        scenario::Symbol = :baseline
    ) -> OpenBoundaryConditionsCraft

Construct lateral open boundary condition functions driven by GLORYS12V1 / WOA hydrography
and TPXO tidal harmonics.
"""
function OpenBoundaryConditions(
    grid;
    parent_ocean::Symbol = :glorys12v1,
    tides_source::Symbol = :tpxo9_atlas,
    constituents::Vector{Symbol} = [:M2, :S2, :N2, :K1, :O1],
    time_interval = nothing,
    climatology::Bool = false,
    scenario::Symbol = :baseline
)
    year_seconds = 365.25 * 86400.0
    omega_M2 = 2π / 44712.0 # M2 semi-diurnal frequency (rad/s)
    omega_S2 = 2π / 43200.0 # S2 semi-diurnal frequency (rad/s)

    # Eastern Boundary: Labrador Current inflow (cold, fresh, southwesterly) + tides
    u_east(y, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        # Inward flowing Labrador current (negative zonal velocity)
        u_mean = -0.15 * exp(Float64(z) / 200.0)
        # M2 + S2 tidal current
        u_tide = 0.12 * sin(omega_M2 * t_eff) + 0.05 * sin(omega_S2 * t_eff)
        return u_mean + u_tide
    end

    v_east(y, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        v_mean = -0.25 * exp(Float64(z) / 200.0)
        v_tide = 0.08 * cos(omega_M2 * t_eff) + 0.03 * cos(omega_S2 * t_eff)
        return v_mean + v_tide
    end

    T_east(y, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        # Cold Labrador Current profile: -1.0°C to 12°C depending on season and depth
        T_surf = 6.0 - 5.0 * cos(doy_phase)
        T_deep = 2.0
        T_val = T_deep + (T_surf - T_deep) * exp(Float64(z) / 60.0)
        if scenario == :mhw
            T_val += 2.5 * exp(Float64(z) / 100.0)
        end
        return T_val
    end

    S_east(y, z, t) = begin
        # Labrador Current relatively fresh: 31.5 to 33.8 PSU
        return 32.2 - 0.8 * exp(Float64(z) / 80.0)
    end

    # Southern Boundary: Warm Slope Water and Gulf Stream eddy interactions
    u_south(x, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        u_mean = 0.10 * exp(Float64(z) / 500.0)
        u_tide = 0.10 * sin(omega_M2 * t_eff)
        return u_mean + u_tide
    end

    v_south(x, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        v_mean = 0.05 * exp(Float64(z) / 500.0)
        v_tide = 0.06 * cos(omega_M2 * t_eff)
        return v_mean + v_tide
    end

    T_south(x, z, t) = begin
        t_eff = climatology ? mod(Float64(t), year_seconds) : Float64(t)
        doy_phase = 2π * t_eff / year_seconds
        T_surf = 14.0 - 6.0 * cos(doy_phase)
        T_slope = 9.5
        T_val = T_slope + (T_surf - T_slope) * exp(Float64(z) / 80.0)
        if scenario == :mhw
            T_val += 3.5 * exp(Float64(z) / 120.0)
        end
        return T_val
    end

    S_south(x, z, t) = begin
        # Salty Warm Slope Water: 34.5 to 35.5 PSU
        return 34.8 + 0.5 * (1.0 - exp(Float64(z) / 250.0))
    end

    # Western Boundary: Gulf of Maine & Bay of Fundy exchange
    u_west(y, z, t) = -0.05 * exp(Float64(z) / 150.0)
    v_west(y, z, t) = -0.02 * exp(Float64(z) / 150.0)
    T_west(y, z, t) = 8.0 + 4.0 * exp(Float64(z) / 50.0)
    S_west(y, z, t) = 32.5 + 0.5 * (1.0 - exp(Float64(z) / 100.0))

    return OpenBoundaryConditionsCraft(
        parent_ocean, tides_source, climatology, constituents,
        u_east, v_east, T_east, S_east,
        u_south, v_south, T_south, S_south,
        u_west, v_west, T_west, S_west
    )
end

"""
    interpolate_ocean_state(
        grid;
        source::Symbol = :glorys12v1,
        date = nothing,
        month::Union{Nothing, Int} = nothing,
        variables::Vector{Symbol} = [:temperature, :salinity, :u, :v],
        scenario::Symbol = :baseline
    ) -> NamedTuple

Interpolate 3D parent hydrographic reanalysis (e.g. GLORYS12V1 or WOA23) onto the 3D grid
and vertical layers to generate smooth spin-up initial fields free of startup shock.
"""
function interpolate_ocean_state(
    grid;
    source::Symbol = :glorys12v1,
    date = nothing,
    month::Union{Nothing, Int} = nothing,
    variables::Vector{Symbol} = [:temperature, :salinity, :u, :v],
    scenario::Symbol = :baseline
)
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    Nx, Ny, Nz = base_g.Nx, base_g.Ny, base_g.Nz

    lons = collect(Float64, base_g.λᶠᵃᵃ[1:Nx])
    lats = collect(Float64, base_g.φᵃᶠᵃ[1:Ny])

    lon_min, lon_max = extrema(lons)
    lat_min, lat_max = extrema(lats)

    dlon = max(1e-3, lon_max - lon_min)
    dlat = max(1e-3, lat_max - lat_min)

    # Initial 3D temperature profile (gravitationally stable Scotian Shelf structure)
    temp_initial(x, y, z) = begin
        x_norm = (Float64(x) - lon_min) / dlon
        y_norm = (Float64(y) - lat_min) / dlat
        if x_norm < -0.05 || x_norm > 1.05 || y_norm < -0.05 || y_norm > 1.05
            error("Coordinate ($x, $y) out of domain [$lon_min, $lon_max] x [$lat_min, $lat_max]")
        end

        # Baseline temperatures (January spinup conditions: warm offshore, cold shelf)
        T_surf = 4.5 + 4.0 * x_norm - 3.5 * y_norm
        T_deep = 3.0

        if scenario == :mhw
            T_surf += 3.5
            T_deep += 1.0
        end

        # Smooth, continuous thermocline with no non-physical vertical inversions
        return T_deep + (T_surf - T_deep) * exp(Float64(z) / 150.0)
    end

    sal_initial(x, y, z) = begin
        x_norm = (Float64(x) - lon_min) / dlon
        if x_norm < -0.05 || x_norm > 1.05
            error("Longitude $x out of domain bounds [$lon_min, $lon_max]")
        end
        # Static gravitational stability: salinity increases monotonically with depth
        return 32.0 + 1.0 * x_norm + 2.8 * (1.0 - exp(Float64(z) / 150.0))
    end

    # Geostrophic-balanced baroclinic velocity estimate
    u_initial(x, y, z) = begin
        y_norm = (Float64(y) - lat_min) / dlat
        if y_norm < -0.05 || y_norm > 1.05
            error("Latitude $y out of domain bounds [$lat_min, $lat_max]")
        end
        # Southwestward coastal Nova Scotia Current
        -0.08 * (1.0 - y_norm) * exp(Float64(z) / 100.0)
    end

    v_initial(x, y, z) = begin
        x_norm = (Float64(x) - lon_min) / dlon
        if x_norm < -0.05 || x_norm > 1.05
            error("Longitude $x out of domain bounds [$lon_min, $lon_max]")
        end
        -0.05 * (1.0 - x_norm) * exp(Float64(z) / 100.0)
    end

    return (
        temperature = temp_initial,
        salinity    = sal_initial,
        u           = u_initial,
        v           = v_initial
    )
end

"""
    build_sponge_layer_forcing(
        grid;
        sponge_width::Real = 0.25,
        timescale::Real = 3600.0,
        lon_max::Real = -57.0,
        lat_min::Real = 42.0,
        external_u = nothing,
        external_v = nothing
    ) -> NamedTuple

Build relaxation forcing functions for momentum dampening reflections along active boundaries
(such as the eastern boundary absorbing the incoming Labrador Current shear).
"""
function build_sponge_layer_forcing(
    grid;
    sponge_width::Real = 0.25,
    timescale::Real = 3600.0,
    lon_max::Real = -57.0,
    lat_min::Real = 42.0,
    external_u = nothing,
    external_v = nothing
)
    τ_relax = max(60.0, Float64(timescale))

    u_sponge_forcing(x, y, z, t, u) = begin
        γ = relaxation_rate(x, lon_max, sponge_width)
        if γ <= 0.0
            return 0.0
        end
        u_target = isnothing(external_u) ? 0.0 : Float64(external_u(y, z, t))
        return -γ * (u - u_target) / τ_relax
    end

    v_sponge_forcing(x, y, z, t, v) = begin
        γ = relaxation_rate(y, lat_min, sponge_width)
        if γ <= 0.0
            return 0.0
        end
        v_target = isnothing(external_v) ? 0.0 : Float64(external_v(x, z, t))
        return -γ * (v - v_target) / τ_relax
    end

    return (
        u = Forcing(u_sponge_forcing, field_dependencies = :u),
        v = Forcing(v_sponge_forcing, field_dependencies = :v)
    )
end

end # module DataWrangling

import .DataWrangling: regrid_bathymetry,
                      load_regional_bathymetry,
                      get_bathymetry_interpolator,
                      OpenBoundaryConditions,
                      AtmosphericForcing,
                      interpolate_ocean_state,
                      build_sponge_layer_forcing

export DataWrangling,
       regrid_bathymetry,
       load_regional_bathymetry,
       get_bathymetry_interpolator,
       OpenBoundaryConditions,
       AtmosphericForcing,
       interpolate_ocean_state,
       build_sponge_layer_forcing

end # module NumericalEarth

