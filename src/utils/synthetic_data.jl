"""
    synthetic_data.jl

Data ingestion, synthetic generation, and NetCDF dataset inspection utilities
for hydrodynamic modeling and particle tracking workflows.
"""

using NCDatasets
using Downloads

"""
    download_sample_data(bathymetry_path="inputs/bathymetry.nc",
                         forcing_path="inputs/forcing.nc";
                         verbose=true)

Download sample Oceananigans bathymetry and surface forcing datasets with
automated fallbacks.

# Mathematical & Physical Context
Sample hydrodynamic datasets provide standard test fields for surface wind
stress and bottom topography to initialize idealized ocean simulations.

# Inputs
- `bathymetry_path::AbstractString`: Destination path for the bathymetry file.
- `forcing_path::AbstractString`: Destination path for the surface forcing file.
- `verbose::Bool`: Whether to print status messages.

# Outputs
- `Tuple{String, String}`: Paths to downloaded bathymetry and forcing files.
"""
function download_sample_data(
    bathymetry_path::AbstractString = joinpath("inputs", "bathymetry.nc"),
    forcing_path::AbstractString = joinpath("inputs", "forcing.nc");
    verbose::Bool = true
)
    mkpath(dirname(bathymetry_path))
    mkpath(dirname(forcing_path))

    bathy_primary = "https://github.com/CliMA/OceananigansArtifacts.jl/" *
                    "raw/main/bathymetries/bathymetry_six_degree.nc"
    bathy_backup  = "https://raw.githubusercontent.com/CliMA/" *
                    "OceananigansArtifacts.jl/main/bathymetries/" *
                    "bathymetry_six_degree.nc"

    forcing_primary = "https://github.com/CliMA/OceananigansArtifacts.jl/" *
                      "raw/main/forcing/surface_forcing.nc"
    forcing_backup  = "https://raw.githubusercontent.com/CliMA/" *
                      "OceananigansArtifacts.jl/main/forcing/surface_forcing.nc"

    if verbose
        println("Downloading sample bathymetry dataset...")
    end
    try
        Downloads.download(bathy_primary, bathymetry_path)
    catch err
        if verbose
            println("Primary bathymetry URL failed ($(err)). Trying backup...")
        end
        Downloads.download(bathy_backup, bathymetry_path)
    end

    if verbose
        println("Downloading sample forcing dataset...")
    end
    try
        Downloads.download(forcing_primary, forcing_path)
    catch err
        if verbose
            println("Primary forcing URL failed ($(err)). Trying backup...")
        end
        Downloads.download(forcing_backup, forcing_path)
    end

    if verbose
        println("Sample datasets downloaded successfully.")
    end
    return (bathymetry_path, forcing_path)
end

"""
    inspect_netcdf(filepath::AbstractString; verbose::Bool=true)

Inspect dimensions, coordinate ranges, and variables of a NetCDF file.

# Inputs
- `filepath::AbstractString`: Path to the NetCDF file.
- `verbose::Bool`: If true, prints summary to stdout.

# Outputs
- `Dict{Symbol, Any}`: Dictionary containing keys `:dimensions`, `:variables`,
  and `:attributes`.
"""
function inspect_netcdf(filepath::AbstractString; verbose::Bool = true)
    if !isfile(filepath)
        error("NetCDF file not found at: $(filepath)")
    end

    file_size = filesize(filepath)
    if file_size == 0
        error("NetCDF file at $(filepath) is empty (0 bytes).")
    end

    info = Dict{Symbol, Any}()
    NCDatasets.Dataset(filepath, "r") do ds
        dims_dict = Dict{String, Int}()
        for (dim_name, dim_len) in ds.dim
            dims_dict[dim_name] = dim_len
        end
        var_names = collect(keys(ds))
        attrs_dict = Dict{String, Any}()
        for (attr_name, attr_val) in ds.attrib
            attrs_dict[attr_name] = attr_val
        end

        info[:dimensions] = dims_dict
        info[:variables] = var_names
        info[:attributes] = attrs_dict
        info[:filesize_bytes] = file_size

        if verbose
            println("--- NetCDF Inspection: $(basename(filepath)) ---")
            println("Path:        $(filepath)")
            println("Dimensions:  $(dims_dict)")
            println("Variables:   $(var_names)")
            println("File size:   $(file_size) bytes")
            println("----------------------------------------------")
        end
    end

    return info
end

"""
    generate_synthetic_bathymetry(
        filepath="inputs/nova_scotia_bathymetry.nc";
        lon_range=(-68.0, -57.0),
        lat_range=(42.0, 47.0),
        n_lon=50,
        n_lat=50,
        inshore_depth=-80.0,
        shelf_slope=600.0,
        land_elevation=80.0,
        include_coastline::Bool=true
    )

Generate synthetic regional bathymetry embedding high-resolution coastline landmasses
(mainland Nova Scotia, Cape Breton, PEI, New Brunswick) with positive terrestrial
elevations (\$z > 0\$) and sloping shelf-to-basin ocean bathymetry (\$z < 0\$).

# Mathematical Formulation
For each grid node \$(\\lambda_i, \\phi_j)\$:
```math
z_b(\\lambda_i, \\phi_j) = \\begin{cases}
+h_{\\text{land}} & \\text{if } (\\lambda_i, \\phi_j) \\in \\mathcal{L} \\quad (\\text{terrestrial land}) \\\\
z_0 - \\Delta z \\left[ \\frac{\\lambda_i - \\lambda_{\\min}}{\\lambda_{\\max} - \\lambda_{\\min}} + \\frac{\\phi_j - \\phi_{\\min}}{\\phi_{\\max} - \\phi_{\\min}} \\right] & \\text{if } (\\lambda_i, \\phi_j) \\notin \\mathcal{L} \\quad (\\text{marine waters})
\\end{cases}
```
where \$\\mathcal{L}\$ is the regional coastline land polygon set, \$z_0 < 0\$ is the coastal
marine shelf datum, and \$\\Delta z\$ scales the offshore continental slope.

# Inputs
- `filepath::AbstractString`: Output NetCDF path.
- `lon_range::Tuple{Real, Real}`: `(min_lon, max_lon)` in degrees East.
- `lat_range::Tuple{Real, Real}`: `(min_lat, max_lat)` in degrees North.
- `n_lon::Int`: Number of longitude grid points.
- `n_lat::Int`: Number of latitude grid points.
- `inshore_depth::Float64`: Baseline coastal marine shelf elevation in meters (default -80.0 m).
- `shelf_slope::Float64`: Elevation difference across the domain in meters (default 600.0 m).
- `land_elevation::Float64`: Terrestrial topographic height in meters (default +80.0 m).
- `include_coastline::Bool`: Whether to embed true terrestrial land boundaries (default `true`).

# Outputs
- `String`: Path to the generated NetCDF file.

# References
- Eaton, B., et al. (2022). NetCDF Climate and Forecast (CF) Metadata Conventions. Version 1.10.
"""
function generate_synthetic_bathymetry(
    filepath::AbstractString = joinpath("inputs", "nova_scotia_bathymetry.nc");
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.0),
    n_lon::Int = 50,
    n_lat::Int = 50,
    inshore_depth::Float64 = -80.0,
    shelf_slope::Float64 = 600.0,
    land_elevation::Float64 = 80.0,
    include_coastline::Bool = true
)
    if lon_range[1] >= lon_range[2]
        error("Invalid longitude range: $(lon_range). lon_min must be < lon_max.")
    end
    if lat_range[1] >= lat_range[2]
        error("Invalid latitude range: $(lat_range). lat_min must be < lat_max.")
    end
    if n_lon <= 1 || n_lat <= 1
        error("Grid resolution (n_lon=$(n_lon), n_lat=$(n_lat)) must be > 1.")
    end

    mkpath(dirname(filepath))

    lon_coords = range(lon_range[1], lon_range[2], length = n_lon)
    lat_coords = range(lat_range[1], lat_range[2], length = n_lat)
    Δlon = lon_range[2] - lon_range[1]
    Δlat = lat_range[2] - lat_range[1]

    NCDatasets.Dataset(filepath, "c") do ds
        NCDatasets.defDim(ds, "lon", n_lon)
        NCDatasets.defDim(ds, "lat", n_lat)

        vlon = NCDatasets.defVar(ds, "lon", Float64, ("lon",),
            attrib = Dict("units" => "degrees_east", "standard_name" => "longitude"))
        vlat = NCDatasets.defVar(ds, "lat", Float64, ("lat",),
            attrib = Dict("units" => "degrees_north", "standard_name" => "latitude"))
        vtopo = NCDatasets.defVar(ds, "elevation", Float64, ("lon", "lat"),
            attrib = Dict("units" => "meters", "standard_name" => "bedrock_altitude"))

        vlon[:] = collect(lon_coords)
        vlat[:] = collect(lat_coords)

        elevation_matrix = Matrix{Float64}(undef, n_lon, n_lat)
        for (i, x) in enumerate(lon_coords), (j, y) in enumerate(lat_coords)
            on_land = include_coastline && is_point_on_land(x, y)
            if on_land
                # Positive elevation for terrestrial landmass
                elevation_matrix[i, j] = max(10.0, land_elevation)
            else
                # Sloping ocean bathymetry from coastal shelf to deep basin
                marine_z = inshore_depth - shelf_slope * (
                    (x - lon_range[1]) / Δlon + (y - lat_range[1]) / Δlat
                )
                elevation_matrix[i, j] = min(-10.0, marine_z)
            end
        end
        vtopo[:, :] = elevation_matrix

        ds.attrib["title"] = "Synthetic Scotian Shelf Bathymetry with Coastline"
        ds.attrib["description"] = "Bathymetric grid with embedded landmass topography for larval particle tracking"
    end

    return filepath
end

"""
    generate_synthetic_forcing(
        filepath="inputs/surface_forcing.nc";
        lon_range=(-68.0, -57.0),
        lat_range=(42.0, 47.0),
        time_range=(0.0, 86400.0),
        n_lon=50,
        n_lat=50,
        n_time=24,
        tau_x_amplitude=0.1,
        tau_y_amplitude=0.0
    )

Generate synthetic time-dependent surface wind stress NetCDF forcing file.

# Mathematical Formulation
Surface kinematic wind stress vector \$\\boldsymbol{\\tau} = (\\tau_x, \\tau_y)\$ is
applied at the upper ocean boundary:
```math
\\tau_x(t) = \\tau_{x, 0} \\sin\\left( \\frac{2\\pi t}{T_{\\text{cycle}}} \\right)
```
where \$\\tau_{x, 0}\$ is the peak wind stress in \$\\text{N m}^{-2}\$ (or kinematic
\$\\text{m}^2\\text{s}^{-2}\$) and \$T_{\\text{cycle}}\$ is the forcing period.

# Inputs
- `filepath::AbstractString`: Output NetCDF path.
- `lon_range::Tuple{Real, Real}`: `(min_lon, max_lon)`.
- `lat_range::Tuple{Real, Real}`: `(min_lat, max_lat)`.
- `time_range::Tuple{Real, Real}`: `(start_time, end_time)` in seconds.
- `n_lon::Int`: Number of longitude points.
- `n_lat::Int`: Number of latitude points.
- `n_time::Int`: Number of time steps.
- `tau_x_amplitude::Float64`: Zonal wind stress amplitude.
- `tau_y_amplitude::Float64`: Meridional wind stress amplitude.

# Outputs
- `String`: Path to the generated NetCDF file.

# References
- Large, W. G., & Pond, S. (1981). Open ocean momentum flux measurements in
  moderate to strong winds. *Journal of Physical Oceanography*, 11(3), 324-336.
"""
function generate_synthetic_forcing(
    filepath::AbstractString = joinpath("inputs", "surface_forcing.nc");
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.0),
    time_range::Tuple{Real, Real} = (0.0, 86400.0),
    n_lon::Int = 50,
    n_lat::Int = 50,
    n_time::Int = 24,
    tau_x_amplitude::Float64 = 0.1,
    tau_y_amplitude::Float64 = 0.0
)
    mkpath(dirname(filepath))

    lon_coords = range(lon_range[1], lon_range[2], length = n_lon)
    lat_coords = range(lat_range[1], lat_range[2], length = n_lat)
    time_coords = range(time_range[1], time_range[2], length = n_time)
    period = time_range[2] > time_range[1] ? (time_range[2] - time_range[1]) : 86400.0

    NCDatasets.Dataset(filepath, "c") do ds
        NCDatasets.defDim(ds, "lon", n_lon)
        NCDatasets.defDim(ds, "lat", n_lat)
        NCDatasets.defDim(ds, "time", n_time)

        vlon = NCDatasets.defVar(ds, "lon", Float64, ("lon",),
            attrib = Dict("units" => "degrees_east"))
        vlat = NCDatasets.defVar(ds, "lat", Float64, ("lat",),
            attrib = Dict("units" => "degrees_north"))
        vtime = NCDatasets.defVar(ds, "time", Float64, ("time",),
            attrib = Dict("units" => "seconds"))
        vtaux = NCDatasets.defVar(
            ds, "tau_x", Float64, ("lon", "lat", "time"),
            attrib = Dict(
                "units" => "N m-2",
                "standard_name" => "surface_downward_eastward_stress"
            )
        )
        vtauy = NCDatasets.defVar(
            ds, "tau_y", Float64, ("lon", "lat", "time"),
            attrib = Dict(
                "units" => "N m-2",
                "standard_name" => "surface_downward_northward_stress"
            )
        )

        vlon[:] = collect(lon_coords)
        vlat[:] = collect(lat_coords)
        vtime[:] = collect(time_coords)

        for (t_idx, t_val) in enumerate(time_coords)
            vtaux[:, :, t_idx] .= tau_x_amplitude * sin(2 * π * t_val / period)
            vtauy[:, :, t_idx] .= tau_y_amplitude * cos(2 * π * t_val / period)
        end

        ds.attrib["title"] = "Synthetic Surface Wind Stress Forcing"
    end

    return filepath
end
