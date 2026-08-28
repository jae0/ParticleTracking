"""
    open_data.jl

Open-access scientific data acquisition, coordinate standardization, 2D regridding,
and aerodynamic drag formulations for realistic ocean hydrodynamic modeling.
"""

using NCDatasets
using Downloads

"""
    wind_speed_to_kinematic_stress(
        u10::Real,
        v10::Real;
        ρ_air::Real = 1.225,
        ρ_water::Real = 1025.0
    )

Convert 10-meter atmospheric wind velocity components \$(u_{10}, v_{10})\$ to
kinematic surface wind stress components \$(\\tau_{x, \\text{kin}}, \\tau_{y, \\text{kin}})\$.

# Mathematical & Empirical Formulation
Following Large & Pond (1981) and Wu (1982), the aerodynamic surface wind stress
vector \$\\boldsymbol{\\tau} = (\\tau_x, \\tau_y)\$ is:
```math
\\boldsymbol{\\tau} = \\rho_{\\text{air}} C_d |\\boldsymbol{u}_{10}| \\boldsymbol{u}_{10}
```
The dimensionless drag coefficient \$C_d\$ depends on 10m wind speed
\$U_{10} = \\sqrt{u_{10}^2 + v_{10}^2}\$:
```math
C_d(U_{10}) = \\begin{cases}
1.2 \\times 10^{-3} & \\text{if } U_{10} \\le 11.0\\text{ m/s} \\\\
(0.49 + 0.065 U_{10}) \\times 10^{-3} & \\text{if } U_{10} > 11.0\\text{ m/s}
\\end{cases}
```
The kinematic wind stress used as upper boundary flux in ocean models is:
```math
\\boldsymbol{\\tau}_{\\text{kinematic}} = \\frac{\\boldsymbol{\\tau}}{\\rho_{\\text{water}}}
= \\left( \\frac{\\rho_{\\text{air}}}{\\rho_{\\text{water}}} \\right) C_d |\\boldsymbol{u}_{10}| \\boldsymbol{u}_{10}
```

# Inputs
- `u10::Real`: Zonal 10-meter wind speed in \$m s^{-1}\$ (positive eastward).
- `v10::Real`: Meridional 10-meter wind speed in \$m s^{-1}\$ (positive northward).
- `ρ_air::Real`: Air density in \$kg m^{-3}\$ (default 1.225 kg/m³).
- `ρ_water::Real`: Seawater reference density in \$kg m^{-3}\$ (default 1025.0 kg/m³).

# Outputs
- `Tuple{Float64, Float64}`: `(tau_x_kinematic, tau_y_kinematic)` in \$m^2 s^{-2}\$.

# References
- Large, W. G., & Pond, S. (1981). JPO, 11(3), 324-336.
- Wu, J. (1982). JGR: Oceans, 87(C12), 9704-9706.
"""
function wind_speed_to_kinematic_stress(
    u10::Real,
    v10::Real;
    ρ_air::Real = 1.225,
    ρ_water::Real = 1025.0
)
    speed = sqrt(u10^2 + v10^2)
    if speed == 0.0
        return (0.0, 0.0)
    end

    cd = if speed <= 11.0
        1.2e-3
    else
        (0.49 + 0.065 * speed) * 1e-3
    end

    factor = (ρ_air / ρ_water) * cd * speed
    tau_x = factor * u10
    tau_y = factor * v10

    return (Float64(tau_x), Float64(tau_y))
end

"""
    regrid_2d_field(
        src_lon::AbstractVector,
        src_lat::AbstractVector,
        src_field::AbstractMatrix,
        target_lon::AbstractVector,
        target_lat::AbstractVector
    )

Interpolate a 2D scalar field (e.g. bathymetry, temperature, wind stress) from
source coordinates onto target model coordinates using bilinear interpolation.

# Mathematical Formulation
For target location \$(\\lambda, \\phi)\$ bounded by source nodes
\$(\\lambda_i, \\lambda_{i+1})\$ and \$(\\phi_j, \\phi_{j+1})\$:
```math
s = \\frac{\\lambda - \\lambda_i}{\\lambda_{i+1} - \\lambda_i}, \\quad
t = \\frac{\\phi - \\phi_j}{\\phi_{j+1} - \\phi_j}
```
```math
f(\\lambda, \\phi) = (1-s)(1-t) f_{i, j} + s(1-t) f_{i+1, j}
                    + (1-s)t f_{i, j+1} + st f_{i+1, j+1}
```

# Inputs
- `src_lon::AbstractVector`: Source longitudes (monotonic).
- `src_lat::AbstractVector`: Source latitudes (monotonic).
- `src_field::AbstractMatrix`: Source data matrix of size `(length(src_lon), length(src_lat))`.
- `target_lon::AbstractVector`: Destination longitude grid coordinates.
- `target_lat::AbstractVector`: Destination latitude grid coordinates.

# Outputs
- `Matrix{Float64}`: Interpolated 2D matrix of size `(length(target_lon), length(target_lat))`.
"""
function regrid_2d_field(
    src_lon::AbstractVector,
    src_lat::AbstractVector,
    src_field::AbstractMatrix,
    target_lon::AbstractVector,
    target_lat::AbstractVector
)
    # Ensure source coordinate sorting (ascending)
    lon_perm = sortperm(collect(src_lon))
    lat_perm = sortperm(collect(src_lat))

    s_lon = collect(Float64, src_lon[lon_perm])
    s_lat = collect(Float64, src_lat[lat_perm])
    s_field = Float64.(src_field[lon_perm, lat_perm])

    n_src_x = length(s_lon)
    n_src_y = length(s_lat)
    n_tgt_x = length(target_lon)
    n_tgt_y = length(target_lat)

    if n_src_x < 2 || n_src_y < 2
        error("Source grid must have at least 2 points in each dimension.")
    end

    interpolated = Matrix{Float64}(undef, n_tgt_x, n_tgt_y)

    for (j_idx, y_val) in enumerate(target_lat)
        # Find bracket in latitude
        j = searchsortedlast(s_lat, y_val)
        j = max(1, min(j, n_src_y - 1))
        t_denom = s_lat[j + 1] - s_lat[j]
        t = t_denom == 0.0 ? 0.0 : (y_val - s_lat[j]) / t_denom
        t = clamp(t, 0.0, 1.0)

        for (i_idx, x_val) in enumerate(target_lon)
            # Find bracket in longitude
            i = searchsortedlast(s_lon, x_val)
            i = max(1, min(i, n_src_x - 1))
            s_denom = s_lon[i + 1] - s_lon[i]
            s = s_denom == 0.0 ? 0.0 : (x_val - s_lon[i]) / s_denom
            s = clamp(s, 0.0, 1.0)

            f00 = s_field[i, j]
            f10 = s_field[i + 1, j]
            f01 = s_field[i, j + 1]
            f11 = s_field[i + 1, j + 1]

            val = (1.0 - s) * (1.0 - t) * f00 +
                  s * (1.0 - t) * f10 +
                  (1.0 - s) * t * f01 +
                  s * t * f11

            interpolated[i_idx, j_idx] = val
        end
    end

    return interpolated
end

"""
    fetch_open_bathymetry(;
        lon_range = (-68.0, -57.0),
        lat_range = (42.0, 47.0),
        output_path = joinpath("inputs", "real_bathymetry.nc"),
        dataset_id = "etopo180",
        stride = 1,
        verbose = true
    )

Retrieve real ocean bathymetry from open scientific data repositories
(NOAA ERDDAP / CoastWatch ETOPO or GEBCO) for a designated regional bounding box.

# Inputs
- `lon_range::Tuple{Real, Real}`: `(min_lon, max_lon)` in degrees East [-180, 180].
- `lat_range::Tuple{Real, Real}`: `(min_lat, max_lat)` in degrees North [-90, 90].
- `output_path::AbstractString`: Destination NetCDF filepath.
- `dataset_id::AbstractString`: ERDDAP dataset ID (default "etopo180" or "nceiEtopo2022").
- `stride::Int`: Subsampling index stride (default 1).
- `verbose::Bool`: Whether to print status messages.

# Outputs
- `String`: Path to the downloaded and verified NetCDF bathymetry file.

# References
- NOAA National Centers for Environmental Information. (2022). NOAA ETOPO 2022
  15 Arc-Second Global Relief Model. NOAA NCEI. DOI: 10.25921/fd1h-fy81
- GEBCO Compilation Group. (2023). GEBCO 2023 Grid.
  DOI: 10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e
- Simons, R. A. (2019). ERDDAP: The Environmental Research Division's Data
  Access Program. NOAA CoastWatch / SWFSC.
"""
function fetch_open_bathymetry(;
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.0),
    output_path::AbstractString = joinpath("inputs", "real_bathymetry.nc"),
    dataset_id::AbstractString = "etopo180",
    stride::Int = 1,
    verbose::Bool = true
)
    mkpath(dirname(output_path))

    min_lat, max_lat = lat_range[1], lat_range[2]
    min_lon, max_lon = lon_range[1], lon_range[2]

    # Standard NOAA CoastWatch ERDDAP endpoint
    base_url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/$(dataset_id).nc"
    query_str = "?altitude[($(min_lat)):$(stride):($(max_lat))]" *
                "[($(min_lon)):$(stride):($(max_lon))]"
    target_url = base_url * query_str

    if verbose
        println("Fetching real bathymetry from NOAA ERDDAP ($(dataset_id))...")
        println("URL: $(target_url)")
    end

    try
        Downloads.download(target_url, output_path)
    catch err
        if verbose
            println("ERDDAP download failed ($(err)). Checking backup mirror...")
        end
        # Alternative NOAA NCEI ETOPO endpoint
        backup_url = "https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/" *
                     "DEM_global_mosaic/ImageServer/exportImage" *
                     "?bbox=$(min_lon),$(min_lat),$(max_lon),$(max_lat)" *
                     "&bboxSR=4326&imageSR=4326&format=tiff&f=json"
        @warn "If direct ERDDAP is offline, please check network or use backup: $(backup_url)"
        rethrow(err)
    end

    if verbose
        println("Bathymetry saved successfully to: $(output_path)")
    end
    return output_path
end

"""
    fetch_open_surface_winds(;
        lon_range = (-68.0, -57.0),
        lat_range = (42.0, 47.0),
        time_iso = "2023-06-01T00:00:00Z",
        output_path = joinpath("inputs", "real_surface_winds.nc"),
        verbose = true
    )

Retrieve real observed/reanalyzed surface winds from NOAA Blended Sea Winds or
NCEP via open ERDDAP services.

# Inputs
- `lon_range::Tuple{Real, Real}`: Longitude bounds in degrees East.
- `lat_range::Tuple{Real, Real}`: Latitude bounds in degrees North.
- `time_iso::AbstractString`: ISO-8601 timestamp (e.g. "2023-06-01T00:00:00Z").
- `output_path::AbstractString`: Destination NetCDF filepath.
- `verbose::Bool`: Whether to log connection progress.

# Outputs
- `String`: Path to the downloaded NetCDF wind file.

# References
- Zhang, H.-M., Bates, J. J., & Reynolds, R. W. (2006). Assessment of composite
  global sampling: Sea surface wind speed. *Geophysical Research Letters*,
  33(17), L17714. DOI: 10.1029/2006GL027086
- Large, W. G., & Pond, S. (1981). Open ocean momentum flux measurements in
  moderate to strong winds. *Journal of Physical Oceanography*, 11(3), 324-336.
"""
function fetch_open_surface_winds(;
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.0),
    time_iso::AbstractString = "2023-06-01T00:00:00Z",
    output_path::AbstractString = joinpath("inputs", "real_surface_winds.nc"),
    verbose::Bool = true
)
    mkpath(dirname(output_path))

    min_lat, max_lat = lat_range[1], lat_range[2]
    min_lon, max_lon = lon_range[1], lon_range[2]

    # NOAA ERDDAP Blended Sea Winds daily dataset (erdBSwinds1day)
    base_url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdBSwinds1day.nc"
    query_str = "?u[($(time_iso)):1:($(time_iso))]" *
                "[($(min_lat)):1:($(max_lat))]" *
                "[($(min_lon)):1:($(max_lon))]," *
                "v[($(time_iso)):1:($(time_iso))]" *
                "[($(min_lat)):1:($(max_lat))]" *
                "[($(min_lon)):1:($(max_lon))]"
    target_url = base_url * query_str

    if verbose
        println("Fetching real wind fields from NOAA ERDDAP (erdBSwinds1day)...")
        println("URL: $(target_url)")
    end

    try
        Downloads.download(target_url, output_path)
    catch err
        if verbose
            println("Wind download error ($(err)).")
        end
        rethrow(err)
    end

    if verbose
        println("Surface wind fields saved to: $(output_path)")
    end
    return output_path
end
