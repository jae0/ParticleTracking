"""
    voronoi_tessellation.jl

Multi-resolution spatial discretization via depth-stratified Voronoi tessellation
for regional marine seascape modeling, larval transport aggregation, and
high-resolution demographic connectivity networks.
"""

using Random
using LinearAlgebra
using Statistics
using NCDatasets

"""
    VoronoiUnit

Structured representation of an individual Voronoi spatial areal unit.

# Fields
- `id::Int`: 1-based unique identifier.
- `lon::Float64`: Centroid / generator longitude in degrees East.
- `lat::Float64`: Centroid / generator latitude in degrees North.
- `depth::Float64`: Bathymetric depth at centroid (negative meters).
- `stratum::Symbol`: Bathymetric stratum (`:core`, `:shallow`, or `:deep`).
- `area_km2::Float64`: Estimated spatial footprint in square kilometers.
- `slope::Float64`: Local topographic bathymetric slope (dimensionless, m/m).
"""
struct VoronoiUnit
    id       :: Int
    lon      :: Float64
    lat      :: Float64
    depth    :: Float64
    stratum  :: Symbol
    area_km2 :: Float64
    slope    :: Float64
end

VoronoiUnit(id::Int, lon::Real, lat::Real, depth::Real, stratum::Symbol, area_km2::Real;
            slope::Real = 0.0) =
    VoronoiUnit(id, Float64(lon), Float64(lat), Float64(depth), stratum,
                Float64(area_km2), Float64(slope))


"""
    VoronoiTessellation

Container holding a depth-stratified Voronoi partition of the marine seascape,
including fast spatial indexing tables for constant-time coordinate lookup.

# Fields
- `units::Vector{VoronoiUnit}`: Vector of all `N` discrete Voronoi units.
- `n_units::Int`: Total count of spatial units.
- `lon_range::Tuple{Float64, Float64}`: Domain longitude bounds.
- `lat_range::Tuple{Float64, Float64}`: Domain latitude bounds.
- `centroids_lon::Vector{Float64}`: Longitude coordinates of generators.
- `centroids_lat::Vector{Float64}`: Latitude coordinates of generators.
- `grid_bins::Dict{Tuple{Int, Int}, Vector{Int}}`: Spatial hash table.
- `bin_ddeg::Float64`: Spatial hash bucket grid resolution in degrees.
- `cos_lat_ref::Float64`: Metric latitude scaling factor.
"""
struct VoronoiTessellation
    units         :: Vector{VoronoiUnit}
    n_units       :: Int
    lon_range     :: Tuple{Float64, Float64}
    lat_range     :: Tuple{Float64, Float64}
    centroids_lon :: Vector{Float64}
    centroids_lat :: Vector{Float64}
    grid_bins     :: Dict{Tuple{Int, Int}, Vector{Int}}
    bin_ddeg      :: Float64
    cos_lat_ref   :: Float64
end

"""
    generate_depth_stratified_voronoi_units(
        bathymetry;
        lon_range::Tuple{Real, Real} = (-68.0, -57.0),
        lat_range::Tuple{Real, Real} = (42.0, 47.5),
        n_units::Int = 5000,
        core_depth_range::Tuple{Real, Real} = (-350.0, -50.0),
        core_prob::Real = 0.8,
        core_min_res_km::Real = 1.5,
        shallow_depth_range::Tuple{Real, Real} = (-50.0, 0.0),
        shallow_prob::Real = 0.1,
        shallow_min_res_km::Real = 5.0,
        deep_depth_range::Tuple{Real, Real} = (-6000.0, -350.0),
        deep_prob::Real = 0.1,
        deep_min_res_km::Real = 10.0,
        seed::Int = 42
    ) -> VoronoiTessellation

Generate a multi-resolution, depth-stratified Voronoi partition of the marine domain
by sampling generator points along bathymetric depth strata using spatial Poisson-disc
rejection filtering.

# Mathematical Formulation & Sampling Strategy
The regional domain is partitioned into three bathymetrically conditioned ecological strata:
1. **Core Nursery Habitat** (\$z \\in [-350, -50]\\text{ m}\$): Receives `core_prob`
   (default 0.80) of all generator units with fine minimum spacing `core_min_res_km`
   (default 1.5 km), yielding dense coverage of commercial snow crab grounds.
2. **Shallow Inshore** (\$z \\in [-50, 0]\\text{ m}\$): Receives `shallow_prob`
   (default 0.10) with intermediate spacing `shallow_min_res_km` (default 5.0 km).
3. **Deep Continental Slope & Abyss** (\$z < -350\\text{ m}\$): Receives `deep_prob`
   (default 0.10) with coarse spacing `deep_min_res_km` (default 10.0 km).

Generator points \$\\boldsymbol{p}_k = (\\lambda_k, \\phi_k)\$ define Voronoi cells:
```math
V_k = \\left\\{ \\boldsymbol{x} \\in \\Omega \\;\\middle|\\;
\\text{dist}(\\boldsymbol{x}, \\boldsymbol{p}_k) \\le \\text{dist}(\\boldsymbol{x}, \\boldsymbol{p}_j)
\\; \\forall j \\ne k \\right\\}
```
where geodesic metric distance is computed using standard spherical projection:
```math
\\text{dist}(\\boldsymbol{p}_1, \\boldsymbol{p}_2) = \\sqrt{
\\left[ (\\lambda_1 - \\lambda_2) \\cos\\phi_{\\text{ref}} \\cdot 111.13 \\right]^2
+ \\left[ (\\phi_1 - \\phi_2) \\cdot 111.13 \\right]^2 }
```

# Inputs
- `bathymetry`: Function `(lon, lat) -> depth` or NetCDF file path string.
- `lon_range, lat_range`: Geographic bounding box coordinates.
- `n_units`: Total target count of discrete Voronoi units (default 5000).
- `core_*, shallow_*, deep_*`: Depth thresholds, probabilities, and spatial resolution limits.
- `seed`: Random number generator seed.

# Outputs
- `VoronoiTessellation`: Container with all units and \$O(1)\$ spatial hash lookup table.
"""
function generate_depth_stratified_voronoi_units(
    bathymetry = :numerical_earth;
    lon_range::Tuple{Real, Real} = (-68.0, -57.0),
    lat_range::Tuple{Real, Real} = (42.0, 47.5),
    n_units::Int = 5000,
    core_depth_range::Tuple{Real, Real} = (-350.0, -50.0),
    core_prob::Union{Nothing, Real} = nothing,
    prob_core::Union{Nothing, Real} = nothing,
    core_min_res_km::Union{Nothing, Real} = nothing,
    min_res_core_km::Union{Nothing, Real} = nothing,
    shallow_depth_range::Tuple{Real, Real} = (-50.0, 0.0),
    shallow_prob::Union{Nothing, Real} = nothing,
    prob_shallow::Union{Nothing, Real} = nothing,
    shallow_min_res_km::Union{Nothing, Real} = nothing,
    min_res_shallow_km::Union{Nothing, Real} = nothing,
    deep_depth_range::Tuple{Real, Real} = (-6000.0, -350.0),
    deep_prob::Union{Nothing, Real} = nothing,
    prob_deep::Union{Nothing, Real} = nothing,
    deep_min_res_km::Union{Nothing, Real} = nothing,
    min_res_deep_km::Union{Nothing, Real} = nothing,
    slope_weighting::Bool = true,
    slope_factor::Float64 = 15.0,
    lloyd_iterations::Int = 1,
    seed::Int = 42
)::VoronoiTessellation
    c_prob = something(core_prob, prob_core, 0.8)
    s_prob = something(shallow_prob, prob_shallow, 0.1)
    d_prob = something(deep_prob, prob_deep, 0.1)
    c_res  = something(core_min_res_km, min_res_core_km, 1.5)
    s_res  = something(shallow_min_res_km, min_res_shallow_km, 5.0)
    d_res  = something(deep_min_res_km, min_res_deep_km, 10.0)

    rng = MersenneTwister(seed)
    ref_lat = 0.5 * (Float64(lat_range[1]) + Float64(lat_range[2]))
    cos_lat = cosd(ref_lat)
    km_per_deg_lat = 111.13
    km_per_deg_lon = 111.13 * cos_lat

    # Target unit counts per stratum
    target_core = round(Int, n_units * Float64(c_prob))
    target_shallow = round(Int, n_units * Float64(s_prob))
    target_deep = n_units - target_core - target_shallow

    # 1. Resolve bathymetry source through NumericalEarth.jl
    bathy_data = if bathymetry isa Function
        nothing
    elseif bathymetry isa NamedTuple && haskey(bathymetry, :elevation)
        bathymetry
    elseif bathymetry isa AbstractString && isfile(bathymetry)
        NumericalEarth.load_regional_bathymetry(
            filepath = bathymetry,
            lon_range = lon_range,
            lat_range = lat_range
        )
    else
        # Default :numerical_earth, :auto, or nothing
        NumericalEarth.load_regional_bathymetry(
            lon_range = lon_range,
            lat_range = lat_range
        )
    end

    # Load optional coastline polygons to strictly exclude land points
    coast_polys = NamedTuple[]
    coast_path = joinpath("inputs", "coastline.dat")
    if isfile(coast_path)
        try
            coast_polys = load_coastline_polygons("inputs")
        catch
        end
    end

    # 2. Extract discrete coordinate pool and calculate bathymetric slope ||∇H||
    pool_core = Tuple{Float64, Float64, Float64, Float64}[]
    pool_shallow = Tuple{Float64, Float64, Float64, Float64}[]
    pool_deep = Tuple{Float64, Float64, Float64, Float64}[]

    if !isnothing(bathy_data)
        raw_lons = bathy_data.lon
        raw_lats = bathy_data.lat
        raw_z = bathy_data.elevation
        nx_b = length(raw_lons)
        ny_b = length(raw_lats)

        # Compute 2D bathymetric slope ||∇H|| (m/m)
        slope_mat = zeros(Float64, nx_b, ny_b)
        for j in 1:ny_b
            j_prev = max(1, j - 1)
            j_next = min(ny_b, j + 1)
            dy_m = max(1.0, (raw_lats[j_next] - raw_lats[j_prev]) * km_per_deg_lat * 1000.0)
            cos_j = cosd(raw_lats[j])
            dx_scale = km_per_deg_lat * cos_j * 1000.0
            for i in 1:nx_b
                i_prev = max(1, i - 1)
                i_next = min(nx_b, i + 1)
                dx_m = max(1.0, (raw_lons[i_next] - raw_lons[i_prev]) * dx_scale)
                dz_dx = (raw_z[i_next, j] - raw_z[i_prev, j]) / dx_m
                dz_dy = (raw_z[i, j_next] - raw_z[i, j_prev]) / dy_m
                slope_mat[i, j] = sqrt(dz_dx^2 + dz_dy^2)
            end
        end

        for j in 1:ny_b
            y = raw_lats[j]
            (y < lat_range[1] || y > lat_range[2]) && continue
            for i in 1:nx_b
                x = raw_lons[i]
                (x < lon_range[1] || x > lon_range[2]) && continue
                d = raw_z[i, j]
                slp = slope_mat[i, j]
                # Negative elevation = marine water, strictly excluding emergent land
                if d < 0.0
                    if !isempty(coast_polys) && is_point_on_land(x, y; coastline = coast_polys)
                        continue
                    end
                    if d >= core_depth_range[1] && d <= core_depth_range[2]
                        push!(pool_core, (x, y, d, slp))
                    elseif d > shallow_depth_range[1] && d <= shallow_depth_range[2]
                        push!(pool_shallow, (x, y, d, slp))
                    elseif d < deep_depth_range[2]
                        push!(pool_deep, (x, y, d, slp))
                    end
                end
            end
        end
    else
        # Continuous function fallback
        bathy_fn = bathymetry isa Function ? bathymetry :
            ((x, y) -> -50.0 - 450.0 * (1.0 - clamp((y - lat_range[1]) /
                       max(1e-3, lat_range[2] - lat_range[1]), 0.0, 1.0))^1.5)
        nx_grid, ny_grid = 300, 200
        for y in range(lat_range[1], lat_range[2], length = ny_grid)
            for x in range(lon_range[1], lon_range[2], length = nx_grid)
                d = Float64(bathy_fn(x, y))
                if d < 0.0
                    if !isempty(coast_polys) && is_point_on_land(x, y; coastline = coast_polys)
                        continue
                    end
                    slp = 0.005
                    if d >= core_depth_range[1] && d <= core_depth_range[2]
                        push!(pool_core, (x, y, d, slp))
                    elseif d > shallow_depth_range[1] && d <= shallow_depth_range[2]
                        push!(pool_shallow, (x, y, d, slp))
                    elseif d < deep_depth_range[2]
                        push!(pool_deep, (x, y, d, slp))
                    end
                end
            end
        end
    end

    # Fallback coordinate synthesis if a stratum has insufficient points in synthetic datasets
    if isempty(pool_shallow) && !isempty(pool_core)
        for p in pool_core[1:min(length(pool_core), 2000)]
            push!(pool_shallow, (p[1], p[2], -30.0, p[4]))
        end
    end
    if isempty(pool_deep) && !isempty(pool_core)
        for p in pool_core[1:min(length(pool_core), 2000)]
            push!(pool_deep, (p[1], p[2], -500.0, p[4]))
        end
    end

    # 3. Stratified Poisson-disc rejection sampling function
    function sample_pool(pool, target_count, min_dist_km, stratum_sym)
        selected = VoronoiUnit[]
        isempty(pool) && return selected
        shuffled = shuffle(rng, pool)

        bin_deg = max(0.01, min_dist_km / km_per_deg_lat)
        bin_grid = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}()

        for (x, y, d, slp) in shuffled
            length(selected) >= target_count && break

            # Modulate minimum spacing by local slope if slope_weighting is enabled
            effective_min_dist = slope_weighting ?
                min_dist_km / (1.0 + slope_factor * slp) : min_dist_km

            bx = round(Int, x / bin_deg)
            by = round(Int, y / bin_deg)

            conflict = false
            for dx in -1:1, dy in -1:1
                pts = get(bin_grid, (bx + dx, by + dy), nothing)
                if !isnothing(pts)
                    for (px, py) in pts
                        dist = sqrt(((x - px) * km_per_deg_lon)^2 +
                                    ((y - py) * km_per_deg_lat)^2)
                        if dist < effective_min_dist
                            conflict = true
                            break
                        end
                    end
                end
                conflict && break
            end

            if !conflict
                u_id = length(selected) + 1
                push!(selected, VoronoiUnit(u_id, x, y, d, stratum_sym, 0.0, slp))
                if !haskey(bin_grid, (bx, by))
                    bin_grid[(bx, by)] = Tuple{Float64, Float64}[]
                end
                push!(bin_grid[(bx, by)], (x, y))
            end
        end

        # Fill deficit with relaxed spacing if min_dist density saturated
        if length(selected) < target_count
            cand_pts = setdiff(shuffled, [(u.lon, u.lat, u.depth, u.slope) for u in selected])
            for (x, y, d, slp) in cand_pts
                length(selected) >= target_count && break
                u_id = length(selected) + 1
                push!(selected, VoronoiUnit(u_id, x, y, d, stratum_sym, 0.0, slp))
            end
        end

        return selected
    end

    units_core = sample_pool(pool_core, target_core, Float64(c_res), :core)
    units_shallow = sample_pool(pool_shallow, target_shallow, Float64(s_res), :shallow)
    units_deep = sample_pool(pool_deep, target_deep, Float64(d_res), :deep)

    # 4. Combine units and assign sequential 1-based IDs
    combined_raw = vcat(units_core, units_shallow, units_deep)

    # Optional Centroidal Voronoi (Lloyd's Relaxation) to smooth cell compactness
    if lloyd_iterations > 0 && !isempty(combined_raw)
        bin_sz = 0.25 # Degrees spatial hash bucket
        for _ in 1:lloyd_iterations
            u_bins = Dict{Tuple{Int, Int}, Vector{Int}}()
            for (idx, u) in enumerate(combined_raw)
                bx = round(Int, u.lon / bin_sz)
                by = round(Int, u.lat / bin_sz)
                if !haskey(u_bins, (bx, by))
                    u_bins[(bx, by)] = Int[]
                end
                push!(u_bins[(bx, by)], idx)
            end

            sum_x = [u.lon for u in combined_raw]
            sum_y = [u.lat for u in combined_raw]
            cnt   = ones(Int, length(combined_raw))

            all_pool = vcat(pool_core, pool_shallow, pool_deep)
            stride_val = max(1, length(all_pool) ÷ 10000)
            for p_idx in 1:stride_val:length(all_pool)
                px, py, _, _ = all_pool[p_idx]
                bx = round(Int, px / bin_sz)
                by = round(Int, py / bin_sz)
                best_i = 0
                best_d = Inf
                for dx in -1:1, dy in -1:1
                    cands = get(u_bins, (bx + dx, by + dy), nothing)
                    isnothing(cands) && continue
                    for c_idx in cands
                        d = ((px - combined_raw[c_idx].lon) * km_per_deg_lon)^2 +
                            ((py - combined_raw[c_idx].lat) * km_per_deg_lat)^2
                        if d < best_d
                            best_d = d
                            best_i = c_idx
                        end
                    end
                end
                if best_i > 0
                    sum_x[best_i] += px
                    sum_y[best_i] += py
                    cnt[best_i] += 1
                end
            end

            for idx in 1:length(combined_raw)
                u = combined_raw[idx]
                nx = clamp(sum_x[idx] / cnt[idx], Float64(lon_range[1]), Float64(lon_range[2]))
                ny = clamp(sum_y[idx] / cnt[idx], Float64(lat_range[1]), Float64(lat_range[2]))
                combined_raw[idx] = VoronoiUnit(u.id, nx, ny, u.depth, u.stratum, u.area_km2, u.slope)
            end
        end
    end
    
    # Calculate domain marine area to apportion realistic unit footprints (km²)
    total_domain_area_km2 = (Float64(lon_range[2]) - Float64(lon_range[1])) * km_per_deg_lon *
                            (Float64(lat_range[2]) - Float64(lat_range[1])) * km_per_deg_lat * 0.75

    nominal_unit_area = total_domain_area_km2 / max(1, length(combined_raw))

    final_units = Vector{VoronoiUnit}(undef, length(combined_raw))
    centroids_lon = Vector{Float64}(undef, length(combined_raw))
    centroids_lat = Vector{Float64}(undef, length(combined_raw))

    for (k, u) in enumerate(combined_raw)
        # Approximate unit area based on stratum resolution
        footprint = if u.stratum == :core
            π * (Float64(c_res) / 2.0)^2 * 1.5
        elseif u.stratum == :shallow
            π * (Float64(s_res) / 2.0)^2 * 1.5
        else
            π * (Float64(d_res) / 2.0)^2 * 1.5
        end
        final_units[k] = VoronoiUnit(k, u.lon, u.lat, u.depth, u.stratum, footprint, u.slope)
        centroids_lon[k] = u.lon
        centroids_lat[k] = u.lat
    end


    # 4. Build fast spatial hash grid index for O(1) cell identification
    bin_ddeg = 0.15 # ~15 km spatial lookup bucket
    grid_bins = Dict{Tuple{Int, Int}, Vector{Int}}()

    for k in 1:length(final_units)
        bx = round(Int, centroids_lon[k] / bin_ddeg)
        by = round(Int, centroids_lat[k] / bin_ddeg)
        if !haskey(grid_bins, (bx, by))
            grid_bins[(bx, by)] = Int[]
        end
        push!(grid_bins[(bx, by)], k)
    end

    return VoronoiTessellation(
        final_units,
        length(final_units),
        (Float64(lon_range[1]), Float64(lon_range[2])),
        (Float64(lat_range[1]), Float64(lat_range[2])),
        centroids_lon,
        centroids_lat,
        grid_bins,
        bin_ddeg,
        cos_lat
    )
end

"""
    find_voronoi_cell(tess::VoronoiTessellation, lon::Real, lat::Real) -> Int

Locate the unique Voronoi areal unit index \$k \\in \\{1, \\dots, N\\}\$ containing
the target coordinates `(lon, lat)`.

Uses spatial hash buckets to bound the nearest-neighbor search, evaluating distances
in projected metric coordinates.
"""
@inline function find_voronoi_cell(tess::VoronoiTessellation, lon::Real, lat::Real)::Int
    x = Float64(lon)
    y = Float64(lat)
    bx = round(Int, x / tess.bin_ddeg)
    by = round(Int, y / tess.bin_ddeg)

    km_lon = 111.13 * tess.cos_lat_ref
    km_lat = 111.13

    best_id = 1
    best_dist2 = Inf

    # Search expanding bucket rings
    for r in 0:3
        for dx in -r:r, dy in -r:r
            if abs(dx) == r || abs(dy) == r
                candidates = get(tess.grid_bins, (bx + dx, by + dy), nothing)
                if !isnothing(candidates)
                    for id in candidates
                        cx = tess.centroids_lon[id]
                        cy = tess.centroids_lat[id]
                        d2 = ((x - cx) * km_lon)^2 + ((y - cy) * km_lat)^2
                        if d2 < best_dist2
                            best_dist2 = d2
                            best_id = id
                        end
                    end
                end
            end
        end
        # Terminate early when candidate is guaranteed closer than ring boundary
        if best_dist2 < ((r * tess.bin_ddeg * km_lon)^2)
            return best_id
        end
    end

    # Fallback to global scan if query coordinate is outside normal spatial index
    if isinf(best_dist2)
        for id in 1:tess.n_units
            cx = tess.centroids_lon[id]
            cy = tess.centroids_lat[id]
            d2 = ((x - cx) * km_lon)^2 + ((y - cy) * km_lat)^2
            if d2 < best_dist2
                best_dist2 = d2
                best_id = id
            end
        end
    end

    return best_id
end

"""
    find_voronoi_cells(
        tess::VoronoiTessellation,
        lons::AbstractVector,
        lats::AbstractVector
    ) -> Vector{Int}

Vectorized batch mapping of spatial coordinate arrays to Voronoi unit identifiers.
"""
function find_voronoi_cells(
    tess::VoronoiTessellation,
    lons::AbstractVector,
    lats::AbstractVector
)::Vector{Int}
    n = min(length(lons), length(lats))
    cell_ids = Vector{Int}(undef, n)
    for i in 1:n
        cell_ids[i] = find_voronoi_cell(tess, lons[i], lats[i])
    end
    return cell_ids
end

"""
    compute_tesselated_connectivity_matrix(
        trajectories,
        tess::VoronoiTessellation;
        settlement_only::Bool = true
    ) -> NamedTuple

Compute high-resolution demographic connectivity metrics and transition probabilities
across depth-stratified Voronoi areal units for a simulated larval cohort.

# Returns
- `NamedTuple`:
  - `tessellation::VoronoiTessellation`: Underlying spatial tessellation.
  - `count_matrix::Matrix{Int}`: \$N \\times N\$ discrete particle transition counts.
  - `probability_matrix::Matrix{Float64}`: Row-normalized transition probabilities \$P_{ij}\$.
  - `retention_vector::Vector{Float64}`: Diagonal self-recruitment fractions \$P_{ii}\$.
  - `settlement_counts::Vector{Int}`: Total settlement arrivals per unit.
  - `settlement_density::Vector{Float64}`: Larval recruits per square kilometer.
  - `strata_matrix::Matrix{Float64}`: \$3 \\times 3\$ connectivity across `:core`, `:shallow`, `:deep`.
  - `strata_labels::Vector{Symbol}`: `[:core, :shallow, :deep]`.
"""
function compute_tesselated_connectivity_matrix(
    trajectories,
    tess::VoronoiTessellation;
    settlement_only::Bool = true
)
    n_units = tess.n_units
    counts = zeros(Int, n_units, n_units)
    settlement_totals = zeros(Int, n_units)

    lons = trajectories.lons
    lats = trajectories.lats
    n_parts, n_times = size(lons)

    settlement_status = hasproperty(trajectories, :settlement_status) ?
        trajectories.settlement_status : fill(:unsettled, n_parts)

    alive = hasproperty(trajectories, :alive) ?
        trajectories.alive : fill(true, n_parts)

    # Check if settled particles exist; if not, use end-point dispersal
    has_settled = any(==( :settled_successful), settlement_status)
    use_settlement = settlement_only && has_settled

    for p in 1:n_parts
        # Release cell (t = 1)
        rel_lon = lons[p, 1]
        rel_lat = lats[p, 1]
        i = find_voronoi_cell(tess, rel_lon, rel_lat)

        # Final or settlement cell
        is_settled = settlement_status[p] == :settled_successful
        if use_settlement && !is_settled
            continue
        end

        dest_lon = lons[p, end]
        dest_lat = lats[p, end]
        j = find_voronoi_cell(tess, dest_lon, dest_lat)

        counts[i, j] += 1
        settlement_totals[j] += 1
    end

    # Row-normalize to transition probability matrix P_ij
    probs = zeros(Float64, n_units, n_units)
    for i in 1:n_units
        row_sum = sum(counts[i, :])
        if row_sum > 0
            probs[i, :] .= counts[i, :] ./ row_sum
        end
    end

    retention = diag(probs)

    densities = [
        settlement_totals[k] / max(0.1, tess.units[k].area_km2) for k in 1:n_units
    ]

    # Aggregate by bathymetric strata (:core, :shallow, :deep)
    strata_order = [:core, :shallow, :deep]
    strata_counts = zeros(Int, 3, 3)

    unit_strata = [u.stratum for u in tess.units]
    for i in 1:n_units
        si = findfirst(==(unit_strata[i]), strata_order)
        for j in 1:n_units
            counts[i, j] == 0 && continue
            sj = findfirst(==(unit_strata[j]), strata_order)
            if !isnothing(si) && !isnothing(sj)
                strata_counts[si, sj] += counts[i, j]
            end
        end
    end

    strata_probs = zeros(Float64, 3, 3)
    for si in 1:3
        r_sum = sum(strata_counts[si, :])
        if r_sum > 0
            strata_probs[si, :] .= strata_counts[si, :] ./ r_sum
        end
    end

    return (
        tessellation = tess,
        count_matrix = counts,
        probability_matrix = probs,
        retention_vector = retention,
        retention_indices = retention,
        settlement_counts = settlement_totals,
        settlement_density = densities,
        strata_matrix = strata_probs,
        macro_matrix = strata_probs,
        strata_labels = strata_order,
        strata_names = strata_order
    )
end

"""
    export_voronoi_geojson(
        tessellation::VoronoiTessellation,
        filepath::AbstractString;
        n_vertices::Int = 6
    ) -> String

Export a `VoronoiTessellation` as a standardized GeoJSON `FeatureCollection` polygon dataset
for GIS visualization (QGIS, ArcGIS, Python Geopandas) and spatial management analysis.

# Mathematical Formulation
Each Voronoi unit \$k\$ with generator \$\\boldsymbol{p}_k = (\\lambda_k, \\phi_k)\$ and footprint
area \$A_k\$ is represented as a closed regular polygon with radius:
```math
R_k = \\sqrt{\\frac{A_k}{\\pi}} \\quad (\\text{km})
```
Vertices at angles \$\\theta_v = 2\\pi v / n_v\$ for \$v \\in \\{0, \\dots, n_v\\}\$ are:
```math
\\lambda_v = \\lambda_k + \\frac{R_k \\cos\\theta_v}{111.13 \\cos\\phi_k}, \\quad
\\phi_v = \\phi_k + \\frac{R_k \\sin\\theta_v}{111.13}
```

# Inputs
- `tessellation::VoronoiTessellation`: Tessellation container.
- `filepath::AbstractString`: Destination path for `.geojson` file.
- `n_vertices::Int`: Number of vertices defining polygon boundary (default 6, hexagon).

# Outputs
- `String`: Filepath of the exported GeoJSON archive.
"""
function export_voronoi_geojson(
    tessellation::VoronoiTessellation,
    filepath::AbstractString;
    n_vertices::Int = 6
)::String
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        println(io, "{\n  \"type\": \"FeatureCollection\",\n  \"features\": [")
        n_u = length(tessellation.units)
        km_lon = 111.13 * tessellation.cos_lat_ref
        km_lat = 111.13
        for (idx, u) in enumerate(tessellation.units)
            radius_km = sqrt(max(0.01, u.area_km2) / π)
            dlon = radius_km / km_lon
            dlat = radius_km / km_lat

            poly_coords = String[]
            for v in 0:n_vertices
                ang = 2π * v / n_vertices
                px = round(u.lon + dlon * cos(ang), digits = 6)
                py = round(u.lat + dlat * sin(ang), digits = 6)
                push!(poly_coords, "[$px, $py]")
            end
            coords_str = join(poly_coords, ", ")

            comma = idx < n_u ? "," : ""
            println(io, "    {")
            println(io, "      \"type\": \"Feature\",")
            println(io, "      \"id\": $(u.id),")
            println(io, "      \"properties\": {")
            println(io, "        \"id\": $(u.id),")
            println(io, "        \"depth\": $(round(u.depth, digits = 1)),")
            println(io, "        \"stratum\": \"$(u.stratum)\",")
            println(io, "        \"area_km2\": $(round(u.area_km2, digits = 2)),")
            println(io, "        \"slope\": $(round(u.slope, digits = 5))")
            println(io, "      },")
            println(io, "      \"geometry\": {")
            println(io, "        \"type\": \"Polygon\",")
            println(io, "        \"coordinates\": [[$(coords_str)]]")
            println(io, "      }")
            println(io, "    }$(comma)")
        end
        println(io, "  ]\n}")
    end
    return filepath
end

