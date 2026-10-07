"""
    grid_bathymetry.jl

Grid construction and immersed boundary setup for shelf hydrodynamic modeling.

ParticleTracking-specific functions requiring Oceananigans grid types.
All generic geospatial operations (bathymetry I/O, regridding, coastline, marine cell
extraction, smoothing, buffering) are in GeoData.Data.
"""

using Oceananigans
using Oceananigans.Grids: Face, Center, znode
using Oceananigans.Architectures: architecture, on_architecture, CPU, GPU
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom
using NCDatasets
using Interpolations
using NumericalEarth: ETOPO2022, regrid_bathymetry, smooth_topography!, BoundingBox
using GeoData
using GeoData.Data
using Random

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
"""
function build_shelf_grid(;
    architecture = :cpu,
    lon_range::Tuple{Real, Real} = (-71.0, -53.0),
    lat_range::Tuple{Real, Real} = (40.0, 48.5),
    z_range::Tuple{Real, Real} = (-3500.0, 0.0),
    z_faces::Union{Nothing, AbstractVector, Function} = nothing,
    grid_size::Tuple{Int, Int, Int} = (50, 50, 10),
    halo::Tuple{Int, Int, Int} = (7, 7, 5),
    topology::Tuple = (Bounded, Bounded, Bounded),
    fallback_to_cpu::Bool = false
)
    if lon_range[1] >= lon_range[2]
        error("Invalid longitude range: $(lon_range). lon_min must be < lon_max.")
    end
    if lat_range[1] >= lat_range[2]
        error("Invalid latitude range: $(lat_range). lat_min must be < lat_max.")
    end
    if any(s <= 0 for s in grid_size)
        error("Grid size dimensions must all be positive integers: $(grid_size)")
    end

    z_specification = if !isnothing(z_faces)
        if z_faces isa AbstractVector
            if length(z_faces) != grid_size[3] + 1
                error("Custom z_faces length $(length(z_faces)) must equal Nz + 1 = $(grid_size[3] + 1)")
            end
            if !issorted(z_faces)
                error("Custom z_faces must be strictly monotonically increasing.")
            end
        end
        z_faces
    else
        if z_range[1] >= z_range[2]
            error("Invalid vertical range: $(z_range). z_min must be < z_max.")
        end
        z_range
    end

    arch = resolve_architecture(architecture; fallback_to_cpu = fallback_to_cpu)
     
    grid = LatitudeLongitudeGrid(
           arch;
           size = grid_size,
           longitude = lon_range,
           latitude = lat_range,
           z = z_specification,
           topology = topology,
           halo = halo
       )
    return grid
end

"""
    load_vertical_grid_csv(filepath::AbstractString) -> Vector{Float64}

Load vertical layer face coordinates from a CSV file (e.g. `scotian_shelf_vertical_grid.csv`).
"""
function load_vertical_grid_csv(filepath::AbstractString)::Vector{Float64}
    if !isfile(filepath)
        error("Vertical grid file does not exist: $(filepath)")
    end
    lines = readlines(filepath)
    if length(lines) < 2
        error("Vertical grid file $(filepath) is empty or missing data.")
    end
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
    if !issorted(faces)
        if issorted(faces; rev = true)
            error(
                "Vertical grid file $(filepath) is written positive-downward (first face " *
                "$(first(faces)) m, last $(last(faces)) m). Oceananigans requires z faces " *
                "NEGATIVE and ascending, ending at 0 m at the surface: rewrite so the first " *
                "row is the deepest layer and the last row's top is 0.0.")
        end
        error("Parsed vertical grid faces from $(filepath) are not monotonically sorted.")
    end
    abs(last(faces)) <= 1e-6 || error(
        "Vertical grid faces from $(filepath) must end at 0.0 m (the surface); got " *
        "$(last(faces)) m.")
    first(faces) < 0 || error(
        "Vertical grid faces from $(filepath) must be negative below the surface; got " *
        "$(first(faces)) m as the first (deepest) face.")
    return faces
end

"""
    build_immersed_grid(
        grid::LatitudeLongitudeGrid,
        bathymetry::Union{AbstractMatrix, AbstractString};
        varname::AbstractString = "elevation"
    )

Wrap a base `LatitudeLongitudeGrid` with an `ImmersedBoundaryGrid` using
`GridFittedBottom` representing ocean seafloor topography.
"""
function build_immersed_grid(
    grid::LatitudeLongitudeGrid,
    bathymetry::Union{AbstractMatrix, AbstractString};
    varname::AbstractString = "elevation"
)
    topo_matrix::Matrix{Float64} = if bathymetry isa AbstractString
        GeoData.Data.load_bathymetry_geodata(bathymetry).elevation
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

    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    arch = architecture(grid)
    z_max = arch isa GPU ? 0.0 : znode(base_g.Nz + 1, base_g, Face())
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
    is_in_bay_of_fundy(lon::Real, lat::Real) -> Bool

Determine whether geographic coordinates `(lon, lat)` in degrees East and North fall
within the hypertidal Bay of Fundy, Minas Basin, or Chignecto Bay exclusion zone.
"""
@inline function is_in_bay_of_fundy(lon::Real, lat::Real)::Bool
    if lat < 44.4 || lat > 46.0
        return false
    end
    if lon < -67.3 || lon > -64.0
        return false
    end
    if lat <= 44.8 && lon > -65.6
        return false
    end
    if lat <= 45.1 && lon > -65.2
        return false
    end
    if lat >= 45.7 && lon > -64.2
        return false
    end
    return true
end

"""
    build_immersed_grid_from_real_data(
        grid::LatitudeLongitudeGrid,
        bathymetry_filepath::AbstractString;
        varname::Union{Nothing, AbstractString} = nothing,
        lon_var::Union{Nothing, AbstractString} = nothing,
        lat_var::Union{Nothing, AbstractString} = nothing,
        min_water_depth::Real = 10.0,
        mask_bay_of_fundy::Bool = false
    )

Construct an `ImmersedBoundaryGrid` by interpolating real-world bathymetry
(e.g. from NOAA ETOPO, GEBCO) onto the target model grid coordinates.
"""
function build_immersed_grid_from_real_data(
    grid::LatitudeLongitudeGrid,
    bathymetry_filepath::AbstractString;
    varname::Union{Nothing, AbstractString} = nothing,
    lon_var::Union{Nothing, AbstractString} = nothing,
    lat_var::Union{Nothing, AbstractString} = nothing,
    min_water_depth::Real = 10.0,
    mask_bay_of_fundy::Bool = false
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

        raw_elev = Array{Float64}(ds[vname][:, :])
        lons = collect(Float64, ds[xname][:])
        lats = collect(Float64, ds[yname][:])

        # Guarantee (n_lon, n_lat) layout
        dim_names   = NCDatasets.dimnames(ds[vname])
        lat_like    = ["lat", "latitude", "y", "nav_lat"]
        needs_tp    = (length(dim_names) >= 1 &&
                       lowercase(string(dim_names[1])) in lat_like)
        elev = needs_tp ? permutedims(raw_elev, (2, 1)) : raw_elev

        (elev, lons, lats)
    end

    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    nx, ny, _ = size(base_g)
    lon_min, lon_max = base_g.λᶠᵃᵃ[1], base_g.λᶠᵃᵃ[base_g.Nx + 1]
    lat_min, lat_max = base_g.φᵃᶠᵃ[1], base_g.φᵃᶠᵃ[base_g.Ny + 1]

    target_lons = range(lon_min, lon_max, length = nx)
    target_lats = range(lat_min, lat_max, length = ny)

    # Perform 2D bilinear regridding
    regridded_topo = GeoData.Data.regrid_2d_field(raw_lon, raw_lat, raw_elevation,
                                     target_lons, target_lats)

    if mask_bay_of_fundy
        n_masked = 0
        for j in 1:ny, i in 1:nx
            if is_in_bay_of_fundy(target_lons[i], target_lats[j])
                regridded_topo[i, j] = 0.0
                n_masked += 1
            end
        end
        @info "Bay of Fundy masking active: set $(n_masked) cells to solid land (0.0 m)."
    end

    # Condition bathymetry
    h_floor = Float64(min_water_depth)
    conditioned_topo = GeoData.Data.smooth_bathymetry(
        regridded_topo,
        passes = 4,
        alpha = 0.5,
        h_min = h_floor
    )

    if mask_bay_of_fundy
        for j in 1:ny, i in 1:nx
            if is_in_bay_of_fundy(target_lons[i], target_lats[j])
                conditioned_topo[i, j] = 0.0
            end
        end
    end

    return build_immersed_grid(grid, conditioned_topo)
end

"""
    extract_grid_coordinates(grid) -> NamedTuple

Extract 1D cell-center spatial coordinates `(lons, lats, depths)` from an
Oceananigans computational grid, returning geographic coordinates in **degrees**
and the vertical coordinate in **metres**.
"""
function extract_grid_coordinates(grid)
    base_g = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    cpu_g = architecture(base_g) isa GPU ? on_architecture(CPU(), base_g) : base_g
    nx, ny, nz = cpu_g.Nx, cpu_g.Ny, cpu_g.Nz

    if !(cpu_g isa LatitudeLongitudeGrid)
        error(
            "extract_grid_coordinates requires a LatitudeLongitudeGrid (or an " *
            "ImmersedBoundaryGrid wrapping one) to report geographic coordinates in " *
            "degrees; received $(nameof(typeof(cpu_g))). Metric-distance grids must be " *
            "converted to geographic coordinates before extraction."
        )
    end

    # Geographic horizontal coordinates (degrees), halo region stripped (interior indices 1:nx, 1:ny).
    lons = [Float64(cpu_g.λᶜᵃᵃ[i]) for i in 1:nx]
    lats = [Float64(cpu_g.φᵃᶜᵃ[j]) for j in 1:ny]

    # Vertical coordinate (metres) is genuinely metric, so `znode` is appropriate.
    depths = [Float64(Oceananigans.Grids.znode(k, cpu_g, Oceananigans.Grids.Center()))
              for k in 1:nz]

    return (lons = lons, lats = lats, depths = depths)
end

"""
    REGIONAL_COASTLINE

Canonical multi-polygon boundary definitions for major landmasses across the
Scotian Shelf, Gulf of St. Lawrence, and northwestern Atlantic region.
Covers the expanded simulation domain [-71, -53]°E / [40, 48.5]°N.

Polygons are closed (first == last vertex) and listed clockwise as viewed
from above (standard GIS convention for exterior rings). Used as fallback
when GeoData coastline is not available.
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
    # 4. New Brunswick, Quebec & Maine Continental Mainland
    # Extended west to -71°E and north to 48.5°N to cover the expanded domain.
    (
        name = "NB-Quebec-Maine Continental Mainland",
        code = :new_brunswick_mainland,
        lons = [
            -64.21, -64.35, -64.79, -65.53, -66.00, -66.47, -67.05, -67.00,
            -67.45, -67.80, -68.50, -69.10, -69.80, -70.20, -70.60, -71.00,
            -71.00, -70.50, -70.00, -69.50, -69.00, -68.50, -68.00, -67.60,
            -67.00, -66.50, -66.00, -65.50, -65.00, -64.85, -64.50, -64.10,
            -64.21
        ],
        lats = [
            45.83,  45.75,  45.60,  45.35,  45.20,  45.06,  45.08,  44.88,
            44.65,  44.45,  44.30,  44.10,  43.80,  43.58,  43.50,  43.00,
            48.50,  48.50,  48.50,  48.50,  48.50,  48.50,  48.50,  48.50,
            48.50,  48.50,  48.50,  47.90,  47.10,  46.68,  46.25,  46.00,
            45.83
        ]
    ),
    # 5. Southern & Western Newfoundland
    # Expanded to cover both Avalon Peninsula and western Newfoundland coasts
    # visible in the domain north of 46.5°N and east of -60°W.
    (
        name = "Southern-Western Newfoundland",
        code = :newfoundland_south,
        lons = [
            -59.13, -58.70, -57.60, -56.00, -55.20, -54.20, -53.40, -53.05,
            -52.70, -53.00, -53.50, -54.10, -54.80, -55.60, -56.40, -57.30,
            -58.00, -58.40, -58.80, -59.13
        ],
        lats = [
            47.57,  47.61,  47.61,  47.55,  47.10,  46.80,  46.70,  46.65,
            47.55,  48.00,  48.30,  48.50,  48.50,  48.50,  48.50,  48.50,
            48.50,  48.30,  48.00,  47.57
        ]
    ),
    # 6. Newfoundland Avalon Peninsula (eastern Newfoundland protruding into domain)
    (
        name = "Newfoundland Avalon Peninsula",
        code = :newfoundland_avalon,
        lons = [
            -52.70, -52.70, -53.30, -53.60, -53.20, -52.90, -52.80, -52.70
        ],
        lats = [
            47.55,  48.00,  48.00,  47.50,  46.65,  46.80,  47.10,  47.55
        ]
    ),
    # 7. Anticosti Island (Gulf of St. Lawrence)
    (
        name = "Anticosti Island",
        code = :anticosti,
        lons = [
            -63.60, -64.20, -64.30, -63.60, -62.80, -62.00, -61.60, -61.80,
            -62.30, -62.80, -63.60
        ],
        lats = [
            49.90,  49.85,  49.45,  49.10,  49.20,  49.35,  49.60,  49.90,
            49.95,  50.00,  49.90
        ]
    ),
    # 8. Magdalen Islands (Îles-de-la-Madeleine) — clustered archipelago
    (
        name = "Magdalen Islands",
        code = :magdalen_islands,
        lons = [-62.00, -61.60, -61.50, -61.80, -62.10, -62.40, -62.00],
        lats = [ 47.30,  47.20,  47.50,  47.65,  47.55,  47.35,  47.30]
    ),
    # 9. Sable Island (closed perimeter)
    (
        name = "Sable Island",
        code = :sable_island,
        lons = [-60.15, -59.90, -59.70, -59.90, -60.15],
        lats = [43.93, 43.95, 43.96, 43.91, 43.93]
    )
]

"""
    get_strata_buffered_envelope(
        polygons::AbstractVector{<:NamedTuple};
        buffer_km::Real = 100.0
    ) -> NamedTuple

Compute the collective bounding envelope across a set of administrative stratum polygons
(e.g., loaded CFAs) expanded by a user-defined buffer distance (default 100.0 km).
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

    buf_lon, buf_lat = GeoData.Data.expand_domain_with_buffer(raw_lon, raw_lat, buffer_km = buffer_km)
    dlon, dlat = GeoData.Data.buffer_distance_to_degrees(buffer_km, 0.5 * (raw_lat[1] + raw_lat[2]))

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