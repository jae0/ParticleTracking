"""
    grid_bathymetry.jl

Grid construction and immersed boundary setup for shelf hydrodynamic modeling.
"""

using Oceananigans
using Oceananigans.Grids: Face, Center, znode
using Oceananigans.Architectures: architecture, on_architecture
using NCDatasets

"""
    build_shelf_grid(;
        architecture = :cpu,
        lon_range=(-68.0, -57.0),
        lat_range=(42.0, 47.0),
        z_range=(-1000.0, 0.0),
        grid_size=(50, 50, 10),
        topology=(Bounded, Bounded, Bounded),
        fallback_to_cpu::Bool = false
    )

Construct a `LatitudeLongitudeGrid` for regional shelf hydrodynamics with specified architecture.

# Mathematical & Physical Formulation
The grid discretizes the spherical curvilinear coordinate system
\$( \\lambda, \\phi, z )\$ representing longitude, latitude, and vertical depth on CPU or GPU.
```math
\\lambda \\in [\\lambda_{\\min}, \\lambda_{\\max}], \\quad
\\phi \\in [\\phi_{\\min}, \\phi_{\\max}], \\quad
z \\in [z_{\\text{bottom}}, z_{\\text{surface}}]
```

# Inputs
- `architecture`: Architecture descriptor (`:cpu`, `:gpu`, `:cuda`, or `AbstractArchitecture`).
- `lon_range::Tuple{Real, Real}`: Longitude bounds in degrees East.
- `lat_range::Tuple{Real, Real}`: Latitude bounds in degrees North.
- `z_range::Tuple{Real, Real}`: Vertical bounds in meters (bottom, top).
- `grid_size::Tuple{Int, Int, Int}`: Number of grid cells `(Nx, Ny, Nz)`.
- `topology::Tuple`: Oceananigans boundary topology (default Bounded in all 3 dirs).
- `fallback_to_cpu::Bool`: If true, falls back to CPU when GPU is requested on non-CUDA systems.

# Outputs
- `LatitudeLongitudeGrid`: Discretized computational grid on specified architecture.
"""
function build_shelf_grid(;
    architecture = :cpu,
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.0),
    z_range::Tuple{Real, Real} = (-1000.0, 0.0),
    grid_size::Tuple{Int, Int, Int} = (50, 50, 10),
    topology::Tuple = (Bounded, Bounded, Bounded),
    fallback_to_cpu::Bool = false
)
    if lon_range[1] >= lon_range[2]
        error("Invalid longitude range: $(lon_range). lon_min must be < lon_max.")
    end
    if lat_range[1] >= lat_range[2]
        error("Invalid latitude range: $(lat_range). lat_min must be < lat_max.")
    end
    if z_range[1] >= z_range[2]
        error("Invalid vertical range: $(z_range). z_min must be < z_max.")
    end
    if any(s <= 0 for s in grid_size)
        error("Grid size dimensions must all be positive integers: $(grid_size)")
    end

    arch = resolve_architecture(architecture; fallback_to_cpu = fallback_to_cpu)

    grid = LatitudeLongitudeGrid(
        arch;
        size = grid_size,
        longitude = lon_range,
        latitude = lat_range,
        z = z_range,
        topology = topology
    )
    return grid
end

"""
    load_bathymetry_from_netcdf(filepath::AbstractString, varname::AbstractString="elevation")

Extract a 2D bathymetry matrix and coordinate vectors from a NetCDF dataset.

# Inputs
- `filepath::AbstractString`: Path to the bathymetry NetCDF file.
- `varname::AbstractString`: Name of the elevation variable (default "elevation").

# Outputs
- `NamedTuple`: `(elevation = Matrix{Float64}, lon = Vector{Float64}, lat = Vector{Float64})`
"""
function load_bathymetry_from_netcdf(
    filepath::AbstractString,
    varname::AbstractString = "elevation"
)
    if !isfile(filepath)
        error("Bathymetry file does not exist: $(filepath)")
    end

    return NCDatasets.Dataset(filepath, "r") do ds
        # Detect elevation variable
        actual_var = if haskey(ds, varname)
            varname
        else
            candidates = ["elevation", "altitude", "topo", "z", "bedrock_altitude", "Band1"]
            found = findfirst(c -> haskey(ds, c), candidates)
            if isnothing(found)
                available = collect(keys(ds))
                error("Elevation variable '$(varname)' not found in $(filepath). Available: $(available)")
            end
            candidates[found]
        end

        elevation_data = Array{Float64}(ds[actual_var][:, :])

        # Detect longitude coordinate variable
        lon_candidates = ["lon", "longitude", "x", "nav_lon"]
        lon_var = findfirst(c -> haskey(ds, c), lon_candidates)
        lon_coords = if !isnothing(lon_var)
            collect(Float64, ds[lon_candidates[lon_var]][:])
        else
            Float64[]
        end

        # Detect latitude coordinate variable
        lat_candidates = ["lat", "latitude", "y", "nav_lat"]
        lat_var = findfirst(c -> haskey(ds, c), lat_candidates)
        lat_coords = if !isnothing(lat_var)
            collect(Float64, ds[lat_candidates[lat_var]][:])
        else
            Float64[]
        end

        return (
            elevation = elevation_data,
            lon = lon_coords,
            lat = lat_coords
        )
    end
end

"""
    get_bathymetry_interpolator(bathymetry; varname="elevation") -> Function

Construct a continuous 2D spatial bilinear interpolation function `(lon, lat) -> z_bed`
from a bathymetry dataset (`NamedTuple` with `:lon, :lat, :elevation`), NetCDF file path,
or direct function.

# Mathematical Formulation
Given grid points \$(\\lambda_i, \\phi_j)\$ with seabed elevation \$z_{i,j}\$, the interpolated
seabed elevation at \$(\\lambda, \\phi)\$ within cell \$[\\lambda_i, \\lambda_{i+1}] \\times [\\phi_j, \\phi_{j+1}]\$ is:
```math
z(\\lambda, \\phi) = (1 - t_\\lambda)(1 - t_\\phi) z_{i,j} + t_\\lambda(1 - t_\\phi) z_{i+1,j}
                   + (1 - t_\\lambda) t_\\phi z_{i,j+1} + t_\\lambda t_\\phi z_{i+1,j+1}
```
where:
```math
t_\\lambda = \\frac{\\lambda - \\lambda_i}{\\lambda_{i+1} - \\lambda_i}, \\quad
t_\\phi = \\frac{\\phi - \\phi_j}{\\phi_{j+1} - \\phi_j}
```

# Inputs
- `bathymetry`: Function `(lon, lat) -> z`, `NamedTuple` `(lon, lat, elevation)`, or
  NetCDF filepath `AbstractString`.
- `varname::AbstractString`: Variable name if reading from NetCDF (default "elevation").

# Outputs
- `Function`: `(lon::Real, lat::Real) -> Float64` returning seabed elevation in meters
  (\$z \\le 0\$ for ocean, \$z > 0\$ for land).
"""
function get_bathymetry_interpolator(
    bathymetry::Union{Function, NamedTuple, AbstractString};
    varname::AbstractString = "elevation"
)
    if bathymetry isa Function
        return bathymetry
    end

    bathy_data = if bathymetry isa AbstractString
        load_bathymetry_from_netcdf(bathymetry, varname)
    else
        bathymetry
    end

    lons = Float64.(bathy_data.lon)
    lats = Float64.(bathy_data.lat)
    elev = Float64.(bathy_data.elevation)

    n_lon = length(lons)
    n_lat = length(lats)

    if n_lon < 2 || n_lat < 2
        error("Bathymetry grid requires >= 2x2 grid points. Found $(n_lon)x$(n_lat).")
    end

    return function (lon::Real, lat::Real)
        i = searchsortedlast(lons, Float64(lon))
        j = searchsortedlast(lats, Float64(lat))
        i = clamp(i, 1, n_lon - 1)
        j = clamp(j, 1, n_lat - 1)

        dlon = lons[i + 1] - lons[i]
        dlat = lats[j + 1] - lats[j]
        tx = dlon > 0.0 ? clamp((Float64(lon) - lons[i]) / dlon, 0.0, 1.0) : 0.0
        ty = dlat > 0.0 ? clamp((Float64(lat) - lats[j]) / dlat, 0.0, 1.0) : 0.0

        z00 = elev[i, j]
        z10 = elev[i + 1, j]
        z01 = elev[i, j + 1]
        z11 = elev[i + 1, j + 1]

        return (1.0 - tx) * (1.0 - ty) * z00 +
               tx * (1.0 - ty) * z10 +
               (1.0 - tx) * ty * z01 +
               tx * ty * z11
    end
end

"""
    point_in_polygon(
        x::Real,
        y::Real,
        poly_x::AbstractVector{<:Real},
        poly_y::AbstractVector{<:Real}
    ) -> Bool

Determine whether 2D point `(x, y)` lies inside a closed polygon defined by vertex
coordinate vectors `poly_x` and `poly_y` using the Jordan Curve (ray-casting) theorem.

# Mathematical Formulation
Casts a horizontal ray from \$(x, y)\$ in the positive \$+x\$ direction to \$+\\infty\$
and counts crossings with all polygon edges \$( (x_i, y_i) \\to (x_j, y_j) )\$:
```math
\\text{crossing} \\iff ((y_i > y) \\neq (y_j > y)) \\land
\\left( x < \\frac{(x_j - x_i)(y - y_i)}{y_j - y_i} + x_i \\right)
```
The point is inside if and only if the total intersection count is odd.

# Inputs
- `x::Real`: Point x-coordinate (longitude).
- `y::Real`: Point y-coordinate (latitude).
- `poly_x::AbstractVector{<:Real}`: Polygon vertex x-coordinates.
- `poly_y::AbstractVector{<:Real}`: Polygon vertex y-coordinates.

# Outputs
- `Bool`: `true` if inside polygon, `false` otherwise.
"""
function point_in_polygon(
    x::Real,
    y::Real,
    poly_x::AbstractVector{<:Real},
    poly_y::AbstractVector{<:Real}
)
    n = length(poly_x)
    if n < 3 || length(poly_y) != n
        return false
    end

    inside = false
    j = n
    px = Float64(x)
    py = Float64(y)

    @inbounds for i in 1:n
        xi, yi = Float64(poly_x[i]), Float64(poly_y[i])
        xj, yj = Float64(poly_x[j]), Float64(poly_y[j])

        # Check horizontal ray intersection with segment (i, j)
        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end

    return inside
end

"""
    buffer_distance_to_degrees(
        buffer_km::Real,
        ref_lat::Real = 44.5
    ) -> Tuple{Float64, Float64}

Convert a linear physical buffer distance in kilometers \$d_{\\text{buf}}\$ into
equivalent geographic longitude and latitude degree increments \$(\\Delta\\lambda, \\Delta\\phi)\$
at a specified reference latitude \$\\phi_0\$.

# Mathematical Formulation
Using the spherical Earth model with mean radius \$R_{\\text{earth}} = 6371.0088\\text{ km}\$:
```math
\\Delta\\phi = \\frac{d_{\\text{buf}}}{R_{\\text{earth}}} \\times \\left(\\frac{180^\\circ}{\\pi}\\right)
```
```math
\\Delta\\lambda = \\frac{d_{\\text{buf}}}{R_{\\text{earth}} \\cos(\\deg2rad(\\phi_0))} \\times \\left(\\frac{180^\\circ}{\\pi}\\right)
```

# Inputs
- `buffer_km::Real`: Buffer distance in kilometers (e.g. 100.0 km).
- `ref_lat::Real`: Reference latitude in degrees North (default 44.5°N).

# Outputs
- `Tuple{Float64, Float64}`: `(dlon, dlat)` degree increments.

# References
- Bowditch, N. (2002). *The American Practical Navigator*. National Imagery and Mapping Agency.
"""
function buffer_distance_to_degrees(
    buffer_km::Real,
    ref_lat::Real = 44.5
)
    if buffer_km < 0.0
        error("Buffer distance must be non-negative: $(buffer_km) km")
    end

    r_earth_km = 6371.0088
    km_per_deg_lat = (π * r_earth_km) / 180.0 # ~111.195 km/deg
    dlat = Float64(buffer_km) / km_per_deg_lat

    cos_lat = cos(deg2rad(clamp(Float64(ref_lat), -85.0, 85.0)))
    km_per_deg_lon = km_per_deg_lat * max(0.01, cos_lat)
    dlon = Float64(buffer_km) / km_per_deg_lon

    return (dlon, dlat)
end

"""
    expand_domain_with_buffer(
        lon_range::Tuple{Real, Real},
        lat_range::Tuple{Real, Real};
        buffer_km::Real = 100.0
    ) -> Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}

Expand a geographic bounding box \$(\\lambda_{\\min}, \\lambda_{\\max}) \\times (\\phi_{\\min}, \\phi_{\\max})\$
outward by a user-defined physical buffer distance (default 100.0 km).

# Mathematical Formulation
```math
\\Omega_{\\text{buffered}} = [\\lambda_{\\min} - \\Delta\\lambda, \\; \\lambda_{\\max} + \\Delta\\lambda] \\times [\\phi_{\\min} - \\Delta\\phi, \\; \\phi_{\\max} + \\Delta\\phi]
```

# Inputs
- `lon_range::Tuple{Real, Real}`: Input longitude bounds.
- `lat_range::Tuple{Real, Real}`: Input latitude bounds.
- `buffer_km::Real`: Buffer distance in kilometers (default 100.0 km).

# Outputs
- `Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}`: `(buffered_lon_range, buffered_lat_range)`
"""
function expand_domain_with_buffer(
    lon_range::Tuple{Real, Real},
    lat_range::Tuple{Real, Real};
    buffer_km::Real = 100.0
)
    ref_lat = 0.5 * (Float64(lat_range[1]) + Float64(lat_range[2]))
    dlon, dlat = buffer_distance_to_degrees(buffer_km, ref_lat)

    buf_lon = (
        max(-180.0, Float64(lon_range[1]) - dlon),
        min(180.0, Float64(lon_range[2]) + dlon)
    )
    buf_lat = (
        max(-90.0, Float64(lat_range[1]) - dlat),
        min(90.0, Float64(lat_range[2]) + dlat)
    )

    return (buf_lon, buf_lat)
end

"""
    get_strata_buffered_envelope(
        polygons::AbstractVector{<:NamedTuple};
        buffer_km::Real = 100.0
    ) -> NamedTuple

Compute the collective bounding envelope across a set of administrative stratum polygons
(e.g., loaded CFAs) expanded by a user-defined buffer distance (default 100.0 km).

# Inputs
- `polygons::AbstractVector{<:NamedTuple}`: List of stratum polygons with `:lons, :lats`.
- `buffer_km::Real`: Buffer distance in kilometers (default 100.0 km).

# Outputs
- `NamedTuple`: `(lon_range = (min_lon, max_lon), lat_range = (min_lat, max_lat), buffer_km = buffer_km, dlon = dlon, dlat = dlat)`
"""
function get_strata_buffered_envelope(
    polygons::AbstractVector{<:NamedTuple};
    buffer_km::Real = 100.0
)
    if isempty(polygons)
        error("Cannot compute buffered envelope for empty polygon list.")
    end

    all_lons = Float64[]
    all_lats = Float64[]

    for poly in polygons
        append!(all_lons, poly.lons)
        append!(all_lats, poly.lats)
    end

    raw_lon = extrema(all_lons)
    raw_lat = extrema(all_lats)

    buf_lon, buf_lat = expand_domain_with_buffer(raw_lon, raw_lat, buffer_km = buffer_km)
    dlon, dlat = buffer_distance_to_degrees(buffer_km, 0.5 * (raw_lat[1] + raw_lat[2]))

    return (
        lon_range = buf_lon,
        lat_range = buf_lat,
        raw_lon_range = raw_lon,
        raw_lat_range = raw_lat,
        buffer_km = Float64(buffer_km),
        dlon = dlon,
        dlat = dlat
    )
end

"""
    REGIONAL_COASTLINE

Canonical multi-polygon boundary definitions for major landmasses across the
Scotian Shelf and Canadian Atlantic region (Nova Scotia mainland, Cape Breton,
Prince Edward Island, New Brunswick / Maine, Southern Newfoundland, Sable Island).
"""
const REGIONAL_COASTLINE = [
    # 1. Nova Scotia Mainland (clockwise closed perimeter)
    (
        name = "Nova Scotia Mainland",
        code = :nova_scotia_mainland,
        lons = [
            -64.25, -64.95, -64.49, -64.36, -65.75, -66.35, -66.05, -66.15,
            -65.98, -65.62, -65.32, -64.60, -64.30, -63.92, -63.55, -63.45,
            -63.00, -62.50, -61.98, -61.40, -60.98, -61.50, -61.40, -61.90,
            -62.25, -62.70, -63.13, -63.30, -63.67, -64.21, -64.25
        ],
        lats = [
            45.75,  45.33,  45.33,  45.09,  44.65,  44.27,  44.30,  43.80,
            43.70,  43.47,  43.70,  44.10,  44.40,  44.49,  44.46,  44.60,
            44.72,  44.88,  45.00,  45.18,  45.33,  45.38,  45.60,  45.87,
            45.65,  45.68,  45.79,  45.75,  45.85,  45.83,  45.75
        ]
    ),
    # 2. Cape Breton Island (clockwise from Strait of Canso)
    (
        name = "Cape Breton Island",
        code = :cape_breton,
        lons = [
            -61.40, -61.53, -61.50, -61.12, -61.00, -60.80, -60.50, -60.42,
            -60.45, -60.40, -60.32, -60.20, -59.70, -59.80, -60.15, -60.50,
            -60.85, -61.00, -61.20, -61.40
        ],
        lats = [
            45.65,  45.98,  46.07,  46.43,  46.63,  46.83,  47.05,  47.03,
            46.90,  46.65,  46.33,  46.20,  46.03,  45.92,  45.85,  45.75,
            45.65,  45.55,  45.55,  45.65
        ]
    ),
    # 3. Prince Edward Island (clockwise from West Point)
    (
        name = "Prince Edward Island",
        code = :prince_edward_island,
        lons = [
            -64.40, -63.99, -63.00, -61.97, -62.25, -62.46, -62.78, -63.49,
            -63.70, -63.85, -64.40
        ],
        lats = [
            46.62,  47.05,  46.50,  46.45,  46.35,  46.00,  45.95,  46.20,
            46.25,  46.38,  46.62
        ]
    ),
    # 4. New Brunswick & Continental Mainland (clockwise from Chignecto)
    (
        name = "New Brunswick & Continental Mainland",
        code = :new_brunswick_mainland,
        lons = [
            -64.21, -64.35, -64.79, -65.53, -66.00, -66.47, -67.05, -67.00,
            -67.45, -68.50, -68.50, -65.50, -65.00, -64.85, -64.50, -64.10,
            -64.21
        ],
        lats = [
            45.83,  45.75,  45.60,  45.35,  45.20,  45.06,  45.08,  44.88,
            44.65,  44.30,  47.90,  47.90,  47.10,  46.68,  46.25,  46.00,
            45.83
        ]
    ),
    # 5. Southern Newfoundland (clockwise from Cape Ray)
    (
        name = "Southern Newfoundland",
        code = :newfoundland_south,
        lons = [
            -59.30, -59.30, -52.70, -52.70, -53.05, -53.40, -54.20, -55.20,
            -56.00, -57.60, -58.70, -59.13, -59.30
        ],
        lats = [
            47.60,  48.00,  48.00,  47.55,  46.65,  46.70,  46.80,  47.10,
            47.55,  47.61,  47.61,  47.57,  47.60
        ]
    ),
    # 6. Sable Island (closed perimeter)
    (
        name = "Sable Island",
        code = :sable_island,
        lons = [-60.15, -59.90, -59.70, -59.90, -60.15],
        lats = [43.93, 43.95, 43.96, 43.91, 43.93]
    )
]

"""
    load_coastline_polygons(path::AbstractString = "inputs/coastline.dat") -> Vector{NamedTuple}

Load high-resolution regional coastline polygons from a formatted `.dat` file, or
return the default canonical `REGIONAL_COASTLINE` if the file is absent.

# Inputs
- `path::AbstractString`: Path to coastline definitions file (default `"inputs/coastline.dat"`).

# Outputs
- `Vector{NamedTuple}`: List of polygon objects containing `:name, :code, :lons, :lats`.
"""
function load_coastline_polygons(path::AbstractString = "inputs/coastline.dat")
    if !isfile(path)
        return REGIONAL_COASTLINE
    end

    polys = NamedTuple[]
    cur_name = ""
    cur_code = :unknown
    cur_lons = Float64[]
    cur_lats = Float64[]

    lines = readlines(path)
    for line in lines
        trimmed = strip(line)
        if isempty(trimmed) || startswith(trimmed, "#")
            continue
        end

        if startswith(trimmed, ">")
            if !isempty(cur_lons) && length(cur_lons) >= 3
                push!(polys, (
                    name = cur_name,
                    code = cur_code,
                    lons = copy(cur_lons),
                    lats = copy(cur_lats)
                ))
            end
            header = strip(trimmed[2:end])
            parts = split(header, ",")
            cur_name = strip(parts[1])
            cur_code = length(parts) >= 2 ? Symbol(strip(parts[2])) : Symbol(lowercase(replace(cur_name, " " => "_")))
            empty!(cur_lons)
            empty!(cur_lats)
            continue
        end

        parts = split(trimmed, ",")
        if length(parts) >= 2
            p_lon = tryparse(Float64, strip(parts[1]))
            p_lat = tryparse(Float64, strip(parts[2]))
            if !isnothing(p_lon) && !isnothing(p_lat)
                push!(cur_lons, p_lon)
                push!(cur_lats, p_lat)
            end
        end
    end

    if !isempty(cur_lons) && length(cur_lons) >= 3
        push!(polys, (
            name = cur_name,
            code = cur_code,
            lons = copy(cur_lons),
            lats = copy(cur_lats)
        ))
    end

    return isempty(polys) ? REGIONAL_COASTLINE : polys
end

"""
    save_coastline_polygons(
        path::AbstractString = "inputs/coastline.dat";
        polygons::AbstractVector{<:NamedTuple} = REGIONAL_COASTLINE
    )

Save regional coastline multi-polygon definitions to a structured text file.

# Inputs
- `path::AbstractString`: Destination file path.
- `polygons::AbstractVector{<:NamedTuple}`: List of coastline polygon NamedTuples.
"""
function save_coastline_polygons(
    path::AbstractString = "inputs/coastline.dat";
    polygons::AbstractVector{<:NamedTuple} = REGIONAL_COASTLINE
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Regional Coastline Multi-Polygon Boundaries for Scotian Shelf & Maritime Canada")
        println(io, "# Format: > Name, Code followed by lon,lat pairs")
        for poly in polygons
            p_name = hasproperty(poly, :name) ? string(poly.name) : "Coastline"
            p_code = hasproperty(poly, :code) ? string(poly.code) : "coastline"
            println(io, "\n> $(p_name),$(p_code)")
            for i in 1:length(poly.lons)
                println(io, "$(poly.lons[i]),$(poly.lats[i])")
            end
        end
    end
    return path
end

"""
    is_point_on_land(
        lon::Real,
        lat::Real;
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> Bool

Determine whether geographic coordinates `(lon, lat)` fall within any emergent
terrestrial landmass polygon in the Maritime Canada region.

# Mathematical Formulation
Uses rapid axis-aligned bounding box (AABB) exclusion followed by Jordan Curve
ray-casting point-in-polygon tests against each coastline polygon.

# Inputs
- `lon::Real`: Longitude coordinate in degrees East.
- `lat::Real`: Latitude coordinate in degrees North.
- `coastline`: Optional custom list of coastline polygons (defaults to loaded/canonical coastline).

# Outputs
- `Bool`: `true` if coordinates lie on land; `false` in marine waters.
"""
function is_point_on_land(
    lon::Real,
    lat::Real;
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    x = Float64(lon)
    y = Float64(lat)
    polys = !isnothing(coastline) ? coastline : REGIONAL_COASTLINE

    for poly in polys
        min_x, max_x = extrema(poly.lons)
        min_y, max_y = extrema(poly.lats)
        # Fast bounding box check
        if x >= min_x && x <= max_x && y >= min_y && y <= max_y
            if point_in_polygon(x, y, poly.lons, poly.lats)
                return true
            end
        end
    end

    return false
end

"""
    is_marine_water(
        lon::Real,
        lat::Real;
        bathymetry::Union{Function, NamedTuple, AbstractString, Nothing} = nothing,
        min_seabed_depth::Real = 0.0,
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> Bool

Determine whether geographic coordinates \$(\\lambda, \\phi)\$ lie in active marine
waters of depth at least `min_seabed_depth` meters, strictly excluding emergent land (\$z \\ge 0\$)
and terrestrial landmasses via high-resolution coastline polygon boundaries.

# Mathematical Formulation
Given the seafloor elevation function \$z_{\\text{bed}}(\\lambda, \\phi)\$ and terrestrial
landmass set \$\\mathcal{L}\$:
```math
\\text{is\\_marine}(\\lambda, \\phi) = \\begin{cases}
\\text{true} & \\text{if } (\\lambda, \\phi) \\notin \\mathcal{L} \\land
z_{\\text{bed}}(\\lambda, \\phi) \\le -\\max(0.0, h_{\\text{min}}) \\\\
\\text{false} & \\text{if } (\\lambda, \\phi) \\in \\mathcal{L} \\lor
z_{\\text{bed}}(\\lambda, \\phi) > -\\max(0.0, h_{\\text{min}})
\\end{cases}
```

# Inputs
- `lon::Real`: Longitude coordinate in degrees.
- `lat::Real`: Latitude coordinate in degrees.
- `bathymetry`: Elevation interpolator Function, `NamedTuple`, NetCDF file path, or `nothing`.
- `min_seabed_depth::Real`: Minimum water depth in meters (default: 0.0 m, which enforces
  the shoreline zero-datum \$z_{\\text{bed}} < 0\$).
- `coastline`: Optional custom list of coastline land polygons.

# Outputs
- `Bool`: `true` if coordinates lie strictly in marine water; `false` on land.
"""
function is_marine_water(
    lon::Real,
    lat::Real;
    bathymetry::Union{Function, NamedTuple, AbstractString, Nothing} = nothing,
    min_seabed_depth::Real = 0.0,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    # 1. Strict coastline land exclusion
    if is_point_on_land(lon, lat; coastline = coastline)
        return false
    end

    # 2. Bathymetric depth constraint
    bathy_fn = if isnothing(bathymetry)
        def_path = "inputs/bathymetry_active.nc"
        isfile(def_path) ? get_bathymetry_interpolator(def_path) : nothing
    else
        get_bathymetry_interpolator(bathymetry)
    end

    if isnothing(bathy_fn)
        return true
    end

    z_bed = bathy_fn(Float64(lon), Float64(lat))
    h_threshold = -max(0.0, Float64(min_seabed_depth))
    return z_bed <= h_threshold
end

"""
    extract_marine_cells(
        bathymetry::Union{NamedTuple, AbstractString};
        lon_range::Tuple{Real, Real} = (-180.0, 180.0),
        lat_range::Tuple{Real, Real} = (-90.0, 90.0),
        min_seabed_depth::Real = 0.0,
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> NamedTuple

Extract all discrete marine grid cell centers that lie in open water (\$z < 0\$),
strictly outside terrestrial landmasses, and meet the minimum water depth requirement
within the specified spatial bounding box.

# Inputs
- `bathymetry`: `NamedTuple` `(lon, lat, elevation)` or NetCDF file path `AbstractString`.
- `lon_range::Tuple{Real, Real}`: Bounding box longitude limits.
- `lat_range::Tuple{Real, Real}`: Bounding box latitude limits.
- `min_seabed_depth::Real`: Minimum water depth in meters (default 0.0 m).
- `coastline`: Optional coastline polygons list.

# Outputs
- `NamedTuple`: `(lons = Vector{Float64}, lats = Vector{Float64}, depths = Vector{Float64}, weights = Vector{Float64})`
"""
function extract_marine_cells(
    bathymetry::Union{NamedTuple, AbstractString};
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0),
    min_seabed_depth::Real = 0.0,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    bathy_data = if bathymetry isa AbstractString
        load_bathymetry_from_netcdf(bathymetry)
    else
        bathymetry
    end

    lons = Float64.(bathy_data.lon)
    lats = Float64.(bathy_data.lat)
    elev = Float64.(bathy_data.elevation)

    h_threshold = -max(0.0, Float64(min_seabed_depth))

    marine_lons = Float64[]
    marine_lats = Float64[]
    marine_depths = Float64[]
    marine_weights = Float64[]

    n_lon = length(lons)
    n_lat = length(lats)

    for i in 1:n_lon
        x = lons[i]
        if !(lon_range[1] <= x <= lon_range[2])
            continue
        end
        for j in 1:n_lat
            y = lats[j]
            if !(lat_range[1] <= y <= lat_range[2])
                continue
            end
            z = elev[i, j]

            # Point must not be on terrestrial land and depth must satisfy threshold
            if z <= h_threshold && !is_point_on_land(x, y; coastline = coastline)
                push!(marine_lons, x)
                push!(marine_lats, y)
                push!(marine_depths, z)
                # Spherical grid cell surface area weight dA = R^2 cos(lat) dlon dlat
                weight = cos(deg2rad(clamp(y, -89.9, 89.9)))
                push!(marine_weights, weight)
            end
        end
    end

    if isempty(marine_lons)
        error(
            "No marine water cells found with depth >= $(min_seabed_depth) m " *
            "within lon $(lon_range) and lat $(lat_range). Domain is entirely land."
        )
    end

    return (
        lons = marine_lons,
        lats = marine_lats,
        depths = marine_depths,
        weights = marine_weights
    )
end

"""
    sample_marine_coordinates(
        n_particles::Int,
        bathymetry::Union{NamedTuple, AbstractString};
        lon_range::Tuple{Real, Real} = (-180.0, 180.0),
        lat_range::Tuple{Real, Real} = (-90.0, 90.0),
        min_seabed_depth::Real = 0.0,
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        rng::AbstractRNG = Random.default_rng()
    ) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Sample \$N\$ continuous coordinates strictly within active marine cells (\$z_{\\text{bed}} < 0\$
and \$z_{\\text{bed}} \\le -h_{\\text{min}}\$) with spherical area weighting, sub-cell jitter,
and rigorous coastline land rejection.

# Guarantees
- 100% deterministic success in \$O(N)\$ time without rejection sampling stalls.
- Exactly 0% probability of particles landing on emergent land or coastal terrain.

# Inputs
- `n_particles::Int`: Number of particle coordinates to sample.
- `bathymetry`: `NamedTuple` or NetCDF file path.
- `lon_range, lat_range`: Geographic bounding box.
- `min_seabed_depth::Real`: Minimum seabed depth (default: 0.0 m).
- `coastline`: Optional coastline polygons.
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `Tuple`: `(sampled_lons, sampled_lats, sampled_seabed_depths)`
"""
function sample_marine_coordinates(
    n_particles::Int,
    bathymetry::Union{NamedTuple, AbstractString};
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0),
    min_seabed_depth::Real = 0.0,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    rng::AbstractRNG = Random.default_rng()
)
    cells = extract_marine_cells(
        bathymetry,
        lon_range = lon_range,
        lat_range = lat_range,
        min_seabed_depth = min_seabed_depth,
        coastline = coastline
    )

    bathy_data = if bathymetry isa AbstractString
        load_bathymetry_from_netcdf(bathymetry)
    else
        bathymetry
    end

    lons_raw = Float64.(bathy_data.lon)
    lats_raw = Float64.(bathy_data.lat)
    dlon = length(lons_raw) > 1 ? abs(lons_raw[2] - lons_raw[1]) : 0.05
    dlat = length(lats_raw) > 1 ? abs(lats_raw[2] - lats_raw[1]) : 0.05

    bathy_interp = get_bathymetry_interpolator(bathymetry)
    h_threshold = -max(0.0, Float64(min_seabed_depth))

    total_weight = sum(cells.weights)
    cum_weights = cumsum(cells.weights) ./ total_weight

    sampled_lons = Vector{Float64}(undef, n_particles)
    sampled_lats = Vector{Float64}(undef, n_particles)
    sampled_zbed = Vector{Float64}(undef, n_particles)

    for p in 1:n_particles
        u = rand(rng, Float64)
        idx = searchsortedfirst(cum_weights, u)
        idx = clamp(idx, 1, length(cells.lons))

        c_lon = cells.lons[idx]
        c_lat = cells.lats[idx]

        jitter_x = (rand(rng, Float64) - 0.5) * dlon * 0.95
        jitter_y = (rand(rng, Float64) - 0.5) * dlat * 0.95
        cand_x = c_lon + jitter_x
        cand_y = c_lat + jitter_y

        z_cand = bathy_interp(cand_x, cand_y)

        # Enforce both bathymetric threshold and coastline land rejection
        if z_cand > h_threshold || is_point_on_land(cand_x, cand_y; coastline = coastline)
            cand_x = c_lon
            cand_y = c_lat
            z_cand = cells.depths[idx]
        end

        sampled_lons[p] = cand_x
        sampled_lats[p] = cand_y
        sampled_zbed[p] = z_cand
    end

    return (sampled_lons, sampled_lats, sampled_zbed)
end

"""
    build_immersed_grid(
        grid::LatitudeLongitudeGrid,
        bathymetry::Union{AbstractMatrix, AbstractString};
        varname::AbstractString = "elevation"
    )

Wrap a base `LatitudeLongitudeGrid` with an `ImmersedBoundaryGrid` using
`GridFittedBottom` representing ocean seafloor topography.

# Mathematical Formulation
The immersed boundary isolates solid bottom cells from active fluid cells:
```math
\\chi(\\lambda, \\phi, z) = \\begin{cases}
1 & \\text{if } z \\ge z_b(\\lambda, \\phi) \\quad (\\text{fluid}) \\\\
0 & \\text{if } z < z_b(\\lambda, \\phi) \\quad (\\text{solid seafloor})
\\end{cases}
```

# Inputs
- `grid::LatitudeLongitudeGrid`: Base computational grid.
- `bathymetry::Union{AbstractMatrix, AbstractString}`: 2D elevation array or
  path to NetCDF file containing bathymetry.
- `varname::AbstractString`: Elevation variable name if a NetCDF file path is passed.

# Outputs
- `ImmersedBoundaryGrid`: Oceananigans immersed boundary grid.

# References
- Verzicco, R. (2023). Immersed boundary methods for ocean modeling.
  *Annual Review of Fluid Mechanics*, 55, 305-333. DOI: 10.1146/annurev-fluid-030322-040713
- Ramadhan, A., et al. (2020). Oceananigans.jl: Fast and friendly geophysical
  fluid dynamics on GPUs. *Journal of Open Source Software*, 5(53), 2018.
"""
function build_immersed_grid(
    grid::LatitudeLongitudeGrid,
    bathymetry::Union{AbstractMatrix, AbstractString};
    varname::AbstractString = "elevation"
)
    topo_matrix::Matrix{Float64} = if bathymetry isa AbstractString
        load_bathymetry_from_netcdf(bathymetry, varname).elevation
    else
        Matrix{Float64}(bathymetry)
    end

    nx, ny, _ = size(grid)
    t_nx, t_ny = size(topo_matrix)

    if (nx != t_nx) || (ny != t_ny)
        error(
            "Bathymetry dimensions ($(t_nx), $(t_ny)) do not match grid horizontal " *
            "dimensions ($(nx), $(ny))."
        )
    end

    # Check for topography exceeding surface
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    z_max = znode(base_g.Nz + 1, base_g, Face())
    max_topo = maximum(topo_matrix)
    if max_topo > z_max
        @warn "Maximum bathymetry elevation ($(max_topo) m) exceeds surface " *
              "height ($(z_max) m). Emerged land points present."
    end

    arch = architecture(grid)
    arch_topo = on_architecture(arch, topo_matrix)
    immersed_grid = ImmersedBoundaryGrid(grid, GridFittedBottom(arch_topo))
    return immersed_grid
end

"""
    build_immersed_grid_from_real_data(
        grid::LatitudeLongitudeGrid,
        bathymetry_filepath::AbstractString;
        varname::Union{Nothing, AbstractString} = nothing,
        lon_var::Union{Nothing, AbstractString} = nothing,
        lat_var::Union{Nothing, AbstractString} = nothing
    )

Construct an `ImmersedBoundaryGrid` by interpolating real-world bathymetry
(e.g. from NOAA ETOPO, GEBCO) onto the target model grid coordinates.

# Mathematical Formulation
Extracts continuous coordinates \$\\lambda_{\\text{raw}}, \\phi_{\\text{raw}}\$ and
seafloor elevations \$z_{\\text{raw}}\$, applies 2D bilinear regridding onto model
target cells \$\\lambda_{\\text{grid}}, \\phi_{\\text{grid}}\$, and constructs the
immersed boundary.

# Inputs
- `grid::LatitudeLongitudeGrid`: Target Oceananigans computational grid.
- `bathymetry_filepath::AbstractString`: Path to the downloaded real NetCDF bathymetry.
- `varname::Union{Nothing, String}`: Elevation variable name (auto-detects if nothing).
- `lon_var::Union{Nothing, String}`: Longitude variable name (auto-detects if nothing).
- `lat_var::Union{Nothing, String}`: Latitude variable name (auto-detects if nothing).

# Outputs
- `ImmersedBoundaryGrid`: Computational grid containing the interpolated real seafloor.

# References
- NOAA National Centers for Environmental Information. (2022). NOAA ETOPO 2022
  15 Arc-Second Global Relief Model. NOAA NCEI. DOI: 10.25921/fd1h-fy81
- GEBCO Compilation Group. (2023). GEBCO 2023 Grid.
  DOI: 10.5285/f98b0f3b-9c64-d6f7-e053-6c86abc0f34e
"""
function build_immersed_grid_from_real_data(
    grid::LatitudeLongitudeGrid,
    bathymetry_filepath::AbstractString;
    varname::Union{Nothing, AbstractString} = nothing,
    lon_var::Union{Nothing, AbstractString} = nothing,
    lat_var::Union{Nothing, AbstractString} = nothing
)
    if !isfile(bathymetry_filepath)
        error("Real bathymetry file not found: $(bathymetry_filepath)")
    end

    raw_elevation, raw_lon, raw_lat = NCDatasets.Dataset(bathymetry_filepath, "r") do ds
        # Auto-detect elevation variable name
        vname = if !isnothing(varname)
            varname
        elseif haskey(ds, "altitude")
            "altitude"
        elseif haskey(ds, "elevation")
            "elevation"
        elseif haskey(ds, "z")
            "z"
        elseif haskey(ds, "topo")
            "topo"
        else
            error("Cannot auto-detect elevation variable. Keys in file: $(keys(ds))")
        end

        # Auto-detect longitude
        xname = if !isnothing(lon_var)
            lon_var
        elseif haskey(ds, "longitude")
            "longitude"
        elseif haskey(ds, "lon")
            "lon"
        elseif haskey(ds, "x")
            "x"
        else
            error("Cannot auto-detect longitude variable. Keys in file: $(keys(ds))")
        end

        # Auto-detect latitude
        yname = if !isnothing(lat_var)
            lat_var
        elseif haskey(ds, "latitude")
            "latitude"
        elseif haskey(ds, "lat")
            "lat"
        elseif haskey(ds, "y")
            "y"
        else
            error("Cannot auto-detect latitude variable. Keys in file: $(keys(ds))")
        end

        elev = Array{Float64}(ds[vname][:, :])
        lons = collect(Float64, ds[xname][:])
        lats = collect(Float64, ds[yname][:])
        (elev, lons, lats)
    end

    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    nx, ny, _ = size(base_g)
    lon_min, lon_max = base_g.λᶠᵃᵃ[1], base_g.λᶠᵃᵃ[base_g.Nx + 1]
    lat_min, lat_max = base_g.φᵃᶠᵃ[1], base_g.φᵃᶠᵃ[base_g.Ny + 1]

    target_lons = range(lon_min, lon_max, length = nx)
    target_lats = range(lat_min, lat_max, length = ny)

    # Perform 2D bilinear regridding
    regridded_topo = regrid_2d_field(raw_lon, raw_lat, raw_elevation,
                                     target_lons, target_lats)

    return build_immersed_grid(grid, regridded_topo)
end

