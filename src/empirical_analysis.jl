"""
    empirical_analysis.jl

Empirical movement analysis, turbulent dispersion estimation, recruitment metrics,
thermal exposure accounting, connectivity matrix derivation, and multi-layer
NetCDF/JLD2 export for snow crab (*Chionoecetes opilio*) Lagrangian particle tracking.

# Mathematical & Methodological Foundation
1. **Empirical Advection & Diffusivity Estimation**:
   Displacement increments \$\\Delta \\boldsymbol{x}_{i,n} = \\boldsymbol{x}_{i}(t_{n+1}) - \\boldsymbol{x}_{i}(t_n)\$
   are binned onto a spatial grid to estimate Eulerian-equivalent mean advection
   \$\\bar{\\boldsymbol{u}}_{\\text{emp}}\$ and turbulent eddy diffusivity \$D_{\\text{emp}}\$:
   ```math
   \\bar{\\boldsymbol{u}}_{\\text{emp}}(k,l) = \\frac{1}{N_{kl}} \\sum_{(i,n) \\in kl} \\frac{\\Delta \\boldsymbol{x}_{i,n}}{\\Delta t_n}
   ```
   ```math
   D_{\\text{emp}}(k,l) = \\frac{1}{4} \\left( \\sigma_u^2(k,l) + \\sigma_v^2(k,l) \\right) \\overline{\\Delta t}
   ```
2. **Gridded Recruitment & Retention**:
   Spatial accounting of natal release densities, settlement sink densities, and
   self-recruitment retention rates across nursery grounds.
3. **Thermal Exposure & Physiological Metrics**:
   Spatiotemporal tracking of degree-day accumulation (\$DD = \\int \\max(0, T - T_0) dt\$),
   ambient exposure temperatures, and thermal mortality risk.
4. **Demographic Connectivity Matrices**:
   Transition probability matrix \$P_{ij} = N_{i \\to j} / N_i\$ between management zones.

# References
- Okubo, A. (1971). Oceanic diffusion diagrams. *Deep Sea Research*, 18(8), 789-802.
  DOI: 10.1016/0011-7471(71)90046-5
- North, E. W., et al. (2009). Manual of recommended practices for modelling physical-biological
  interactions during fish early life. *ICES Cooperative Research Report*, No. 295, 111 pp.
- Kuhn, P. S., & Choi, J. S. (2011). Influence of temperature on embryo incubation and larval
  development in snow crab (*Chionoecetes opilio*). *Fisheries Research*, 107(1-3), 81-87.
  DOI: 10.1016/j.fishres.2010.10.011
"""

using NCDatasets
using JLD2
using Statistics
using LinearAlgebra

"""
    estimate_empirical_movement(
        trajectories::NamedTuple;
        lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
        lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30)
    )

Estimate gridded mean empirical advection velocity \$(\\bar{u}_{\\text{emp}}, \\bar{v}_{\\text{emp}})\$
and isotropic turbulent diffusivity \$D_{\\text{emp}}\$ from Lagrangian particle displacements.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`, containing `.lons`, `.lats`,
  `.depths`, `.times`.
- `lon_bins::AbstractVector{<:Real}`: Grid cell edges for longitude in degrees.
- `lat_bins::AbstractVector{<:Real}`: Grid cell edges for latitude in degrees.
- `min_samples::Int`: Minimum particle transit count required per cell to estimate
  velocities and diffusivity (default: 1; cells with fewer transits are assigned `NaN`).

# Outputs
- `NamedTuple`:
  - `lon_centers`: Longitudinal cell midpoints.
  - `lat_centers`: Latitudinal cell midpoints.
  - `u_mean`: Empirical zonal velocity field (\$\\text{m s}^{-1}\$).
  - `v_mean`: Empirical meridional velocity field (\$\\text{m s}^{-1}\$).
  - `speed_mean`: Empirical speed field (\$\\text{m s}^{-1}\$).
  - `diffusivity`: Empirical eddy diffusivity field (\$D_{\\text{emp}}\$, \$\\text{m}^2 \\text{s}^{-1}\$).
  - `sample_count`: Number of particle transit steps per cell.
"""
function estimate_empirical_movement(
    trajectories::NamedTuple;
    lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
    lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30),
    min_samples::Int = 1
)
    n_lon = length(lon_bins) - 1
    n_lat = length(lat_bins) - 1

    lon_centers = [0.5 * (lon_bins[i] + lon_bins[i + 1]) for i in 1:n_lon]
    lat_centers = [0.5 * (lat_bins[j] + lat_bins[j + 1]) for j in 1:n_lat]

    u_accum = zeros(Float64, n_lon, n_lat)
    v_accum = zeros(Float64, n_lon, n_lat)
    u_sq_accum = zeros(Float64, n_lon, n_lat)
    v_sq_accum = zeros(Float64, n_lon, n_lat)
    counts = zeros(Int, n_lon, n_lat)
    dt_accum = zeros(Float64, n_lon, n_lat)

    lons = trajectories.lons
    lats = trajectories.lats
    times = trajectories.times

    n_particles, n_times = size(lons)
    r_earth = 6.371e6

    for p in 1:n_particles
        for t in 1:(n_times - 1)
            x0, y0 = lons[p, t], lats[p, t]
            x1, y1 = lons[p, t + 1], lats[p, t + 1]
            dt = times[t + 1] - times[t]

            if dt <= 0.0 || isnan(x0) || isnan(x1) || isnan(y0) || isnan(y1)
                continue
            end

            # Binning index logic: Clamp x0, y0 to bins.
            # Using max(1, min(n, index)) ensures indices stay within [1, n_bins]
            i_lon = clamp(searchsortedlast(lon_bins, x0), 1, n_lon)
            j_lat = clamp(searchsortedlast(lat_bins, y0), 1, n_lat)
            
            # Metric spherical displacements
            mean_lat_rad = deg2rad(0.5 * (y0 + y1))
            dx_m = r_earth * cos(mean_lat_rad) * deg2rad(x1 - x0)
            dy_m = r_earth * deg2rad(y1 - y0)

            u_inst = dx_m / dt
            v_inst = dy_m / dt

            u_accum[i_lon, j_lat] += u_inst
            v_accum[i_lon, j_lat] += v_inst
            u_sq_accum[i_lon, j_lat] += u_inst^2
            v_sq_accum[i_lon, j_lat] += v_inst^2
            dt_accum[i_lon, j_lat] += dt
            counts[i_lon, j_lat] += 1
        end
    end

    u_mean = zeros(Float64, n_lon, n_lat)
    v_mean = zeros(Float64, n_lon, n_lat)
    speed_mean = zeros(Float64, n_lon, n_lat)
    diffusivity = zeros(Float64, n_lon, n_lat)

    for i in 1:n_lon, j in 1:n_lat
        n_c = counts[i, j]
        if n_c >= min_samples
            u_bar = u_accum[i, j] / n_c
            v_bar = v_accum[i, j] / n_c
            u_mean[i, j] = u_bar
            v_mean[i, j] = v_bar
            speed_mean[i, j] = hypot(u_bar, v_bar)

            var_u = max(0.0, (u_sq_accum[i, j] / n_c) - u_bar^2)
            var_v = max(0.0, (v_sq_accum[i, j] / n_c) - v_bar^2)
            mean_dt = dt_accum[i, j] / n_c

            # Taylor dispersion formulation D = 0.25 * (var_u + var_v) * mean_dt
            # Note: Ensure mean_dt satisfies CFL condition at grid resolution scale
            diffusivity[i, j] = 0.25 * (var_u + var_v) * mean_dt
        else
            u_mean[i, j] = NaN
            v_mean[i, j] = NaN
            speed_mean[i, j] = NaN
            diffusivity[i, j] = NaN
        end
    end

    return (
        lon_centers = lon_centers,
        lat_centers = lat_centers,
        lon_bins = collect(lon_bins),
        lat_bins = collect(lat_bins),
        u_mean = u_mean,
        v_mean = v_mean,
        speed_mean = speed_mean,
        diffusivity = diffusivity,
        sample_count = counts
    )
end

"""
    compute_gridded_recruitment_metrics(
        trajectories::NamedTuple;
        lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
        lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30)
    )

Compute gridded release density, final settlement sink density, self-retention density,
and settlement success fractions from particle trajectories.

# Inputs
- `trajectories::NamedTuple`: Particle tracking output.
- `lon_bins::AbstractVector{<:Real}`: Longitude bin edges.
- `lat_bins::AbstractVector{<:Real}`: Latitude bin edges.

# Outputs
- `NamedTuple`:
  - `release_density`: Count of particles initiated per grid cell.
  - `settlement_density`: Count of all settled particles per grid cell.
  - `successful_settlement_density`: Count of particles meeting benthic habitat criteria.
  - `retention_density`: Count of particles settled in their natal release cell.
  - `retention_rate`: Cell-level self-recruitment ratio.
  - `success_rate`: Ratio of successful recruits to total settling particles per cell.
  - `mean_settlement_age_days`: Mean drift duration in days for settling cohorts.
"""
function compute_gridded_recruitment_metrics(
    trajectories::NamedTuple;
    lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
    lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30)
)
    n_lon = length(lon_bins) - 1
    n_lat = length(lat_bins) - 1

    lon_centers = [0.5 * (lon_bins[i] + lon_bins[i + 1]) for i in 1:n_lon]
    lat_centers = [0.5 * (lat_bins[j] + lat_bins[j + 1]) for j in 1:n_lat]

    release_density = zeros(Int, n_lon, n_lat)
    settlement_density = zeros(Int, n_lon, n_lat)
    successful_settlement_density = zeros(Int, n_lon, n_lat)
    retention_density = zeros(Int, n_lon, n_lat)
    settlement_duration_accum = zeros(Float64, n_lon, n_lat)

    lons = trajectories.lons
    lats = trajectories.lats
    times = trajectories.times
    statuses = trajectories.settlement_status
    n_particles = size(lons, 1)
    ages = hasproperty(trajectories, :settlement_age) ?
           trajectories.settlement_age : fill(times[end], n_particles)

    for p in 1:n_particles
        lon_start, lat_start = lons[p, 1], lats[p, 1]
        lon_end, lat_end = lons[p, end], lats[p, end]
        stat = statuses[p]
        age_sec = ages[p]

        i_start = searchsortedlast(lon_bins, lon_start)
        j_start = searchsortedlast(lat_bins, lat_start)

        i_end = searchsortedlast(lon_bins, lon_end)
        j_end = searchsortedlast(lat_bins, lat_end)
        if 1 <= i_end <= n_lon && 1 <= j_end <= n_lat
            settlement_density[i_end, j_end] += 1
            settlement_duration_accum[i_end, j_end] += (age_sec / 86400.0)

            if stat == :settled_successful
                successful_settlement_density[i_end, j_end] += 1
            end

            if (i_start == i_end) && (j_start == j_end)
                retention_density[i_end, j_end] += 1
            end
        end
    end

    retention_rate = zeros(Float64, n_lon, n_lat)
    success_rate = zeros(Float64, n_lon, n_lat)
    mean_settlement_age_days = zeros(Float64, n_lon, n_lat)

    for i in 1:n_lon, j in 1:n_lat
        n_rel = release_density[i, j]
        n_set = settlement_density[i, j]
        n_succ = successful_settlement_density[i, j]

        retention_rate[i, j] = n_rel > 0 ? (retention_density[i, j] / n_rel) : 0.0
        success_rate[i, j] = n_set > 0 ? (n_succ / n_set) : 0.0
        mean_settlement_age_days[i, j] = n_set > 0 ?
            (settlement_duration_accum[i, j] / n_set) : NaN
    end

    return (
        lon_centers = lon_centers,
        lat_centers = lat_centers,
        lon_bins = collect(lon_bins),
        lat_bins = collect(lat_bins),
        release_density = release_density,
        settlement_density = settlement_density,
        successful_settlement_density = successful_settlement_density,
        retention_density = retention_density,
        retention_rate = retention_rate,
        success_rate = success_rate,
        mean_settlement_age_days = mean_settlement_age_days
    )
end

"""
    compute_gridded_thermal_metrics(
        trajectories::NamedTuple;
        lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
        lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30)
    )

Calculate gridded mean cumulative degree-days, mean in-situ thermal exposure
temperature, and cumulative thermal mortality fractions across the spatial domain.

# Inputs
- `trajectories::NamedTuple`: Particle tracking output containing `.temperatures`,
  `.degree_days`, `.survival_probability`.
- `lon_bins::AbstractVector{<:Real}`: Longitude bin edges.
- `lat_bins::AbstractVector{<:Real}`: Latitude bin edges.

# Outputs
- `NamedTuple`:
  - `mean_degree_days`: Mean degree-days accumulated by particles residing in cell.
  - `mean_exposure_temperature`: Average ambient temperature (\$^\\circ\\text{C}\$).
  - `max_exposure_temperature`: Maximum peak exposure temperature (\$^\\circ\\text{C}\$).
  - `mean_survival_probability`: Mean physiological survival probability in cell.
  - `thermal_mortality_fraction`: Fraction of mortality attributed to thermal stress.
"""
function compute_gridded_thermal_metrics(
    trajectories::NamedTuple;
    lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
    lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30)
)
    n_lon = length(lon_bins) - 1
    n_lat = length(lat_bins) - 1

    lon_centers = [0.5 * (lon_bins[i] + lon_bins[i + 1]) for i in 1:n_lon]
    lat_centers = [0.5 * (lat_bins[j] + lat_bins[j + 1]) for j in 1:n_lat]

    temp_accum = zeros(Float64, n_lon, n_lat)
    temp_max = fill(-Inf, n_lon, n_lat)
    dd_accum = zeros(Float64, n_lon, n_lat)
    surv_accum = zeros(Float64, n_lon, n_lat)
    counts = zeros(Int, n_lon, n_lat)

    lons = trajectories.lons
    lats = trajectories.lats
    n_particles, n_times = size(lons)

    temps = hasproperty(trajectories, :temperatures) ?
            trajectories.temperatures : fill(4.0, n_particles, n_times)
    dds = hasproperty(trajectories, :degree_days_timeseries) ?
          trajectories.degree_days_timeseries :
          (hasproperty(trajectories, :degree_days) && trajectories.degree_days isa AbstractMatrix ?
           trajectories.degree_days :
           repeat(hasproperty(trajectories, :degree_days) ? trajectories.degree_days : fill(40.0, n_particles), 1, n_times))
    survs = hasproperty(trajectories, :survival_probability) ?
            trajectories.survival_probability : fill(0.95, n_particles, n_times)

    for p in 1:n_particles
        for t in 1:n_times
            x = lons[p, t]
            y = lats[p, t]
            temp = temps[p, t]
            dd = dds[p, t]
            surv = survs[p, t]

            if isnan(x) || isnan(y) || isnan(temp)
                continue
            end

            i = searchsortedlast(lon_bins, x)
            j = searchsortedlast(lat_bins, y)

            if 1 <= i <= n_lon && 1 <= j <= n_lat
                temp_accum[i, j] += temp
                if temp > temp_max[i, j]
                    temp_max[i, j] = temp
                end
                dd_accum[i, j] += dd
                surv_accum[i, j] += surv
                counts[i, j] += 1
            end
        end
    end

    mean_degree_days = zeros(Float64, n_lon, n_lat)
    mean_exposure_temperature = zeros(Float64, n_lon, n_lat)
    max_exposure_temperature = zeros(Float64, n_lon, n_lat)
    mean_survival_probability = zeros(Float64, n_lon, n_lat)
    thermal_mortality_fraction = zeros(Float64, n_lon, n_lat)

    for i in 1:n_lon, j in 1:n_lat
        n_c = counts[i, j]
        if n_c > 0
            mean_degree_days[i, j] = dd_accum[i, j] / n_c
            mean_exposure_temperature[i, j] = temp_accum[i, j] / n_c
            max_exposure_temperature[i, j] = temp_max[i, j]
            surv_bar = surv_accum[i, j] / n_c
            mean_survival_probability[i, j] = surv_bar
            thermal_mortality_fraction[i, j] = max(0.0, 1.0 - surv_bar)
        else
            mean_degree_days[i, j] = NaN
            mean_exposure_temperature[i, j] = NaN
            max_exposure_temperature[i, j] = NaN
            mean_survival_probability[i, j] = NaN
            thermal_mortality_fraction[i, j] = NaN
        end
    end

    return (
        lon_centers = lon_centers,
        lat_centers = lat_centers,
        lon_bins = collect(lon_bins),
        lat_bins = collect(lat_bins),
        mean_degree_days = mean_degree_days,
        mean_exposure_temperature = mean_exposure_temperature,
        max_exposure_temperature = max_exposure_temperature,
        mean_survival_probability = mean_survival_probability,
        thermal_mortality_fraction = thermal_mortality_fraction,
        sample_count = counts
    )
end



"""
    closest_point_on_segment(px::Float64, py::Float64, ax::Float64, ay::Float64, bx::Float64, by::Float64) -> Tuple{Float64, Float64, Float64}

Compute the orthogonal/closest point on line segment AB to point P, returning `(qx, qy, dist_squared)`.
"""
function closest_point_on_segment(
    px::Float64, py::Float64,
    ax::Float64, ay::Float64,
    bx::Float64, by::Float64
)
    vx = bx - ax
    vy = by - ay
    len_sq = vx * vx + vy * vy
    if len_sq < 1e-12
        dx = px - ax
        dy = py - ay
        return (ax, ay, dx * dx + dy * dy)
    end
    wx = px - ax
    wy = py - ay
    t = clamp((wx * vx + wy * vy) / len_sq, 0.0, 1.0)
    qx = ax + t * vx
    qy = ay + t * vy
    dx = px - qx
    dy = py - qy
    return (qx, qy, dx * dx + dy * dy)
end

"""
    snap_point_to_coastline(
        x::Float64,
        y::Float64,
        coast_polys::AbstractVector{<:NamedTuple}
    ) -> Tuple{Float64, Float64}

Find the nearest point on regional coastline polygon boundaries for coordinate `(x, y)`.
If `(x, y)` lies within a specific landmass polygon, projection is strictly constrained
to that containing landmass boundary to avoid jumping across straits or open water.
"""
function snap_point_to_coastline(
    x::Float64,
    y::Float64,
    coast_polys::AbstractVector{<:NamedTuple}
)
    # Check if point falls inside a specific land polygon
    target_polys = NamedTuple[]
    for poly in coast_polys
        if point_in_polygon(x, y, poly.lons, poly.lats)
            push!(target_polys, poly)
        end
    end

    # If inside specific land polygon(s), constrain search to containing landmass
    candidate_polys = isempty(target_polys) ? coast_polys : target_polys

    best_qx = x
    best_qy = y
    best_dist_sq = Inf

    for poly in candidate_polys
        c_lons = poly.lons
        c_lats = poly.lats
        m = length(c_lons)
        for i in 1:(m - 1)
            ax, ay = Float64(c_lons[i]), Float64(c_lats[i])
            bx, by = Float64(c_lons[i + 1]), Float64(c_lats[i + 1])
            qx, qy, dsq = closest_point_on_segment(x, y, ax, ay, bx, by)
            if dsq < best_dist_sq
                best_dist_sq = dsq
                best_qx = qx
                best_qy = qy
            end
        end
    end

    return (best_qx, best_qy)
end

"""
    intersect_polygon_with_coastline(
        poly_lons::AbstractVector{<:Real},
        poly_lats::AbstractVector{<:Real};
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> Tuple{Vector{Float64}, Vector{Float64}}

Intersect and clip a management zone polygon with the regional coastline using topology-preserving
coastal snapping. Terrestrial interior vertices are projected to the nearest high-water coastline
segment of their containing landmass rather than deleted, preserving contiguous boundary ordering,
shared inter-strata borders, and preventing geometric shards, slivers, or artificial spatial gaps.

# Inputs
- `poly_lons`: Vector of polygon vertex longitudes.
- `poly_lats`: Vector of polygon vertex latitudes.
- `coastline`: Optional coastline land polygons list.

# Outputs
- `Tuple{Vector{Float64}, Vector{Float64}}`: `(clipped_lons, clipped_lats)` representing
  the marine-clipped polygon.
"""
function intersect_polygon_with_coastline(
    poly_lons::AbstractVector{<:Real},
    poly_lats::AbstractVector{<:Real};
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    n = length(poly_lons)
    if n < 3
        return (Float64.(poly_lons), Float64.(poly_lats))
    end

    coast = isnothing(coastline) ? load_coastline_polygons() : coastline

    out_lons = Float64[]
    out_lats = Float64[]

    for i in 1:n
        x = Float64(poly_lons[i])
        y = Float64(poly_lats[i])

        if is_point_on_land(x, y; coastline = coast)
            # Topology-preserving projection onto containing landmass coastline segment
            qx, qy = snap_point_to_coastline(x, y, coast)
            push!(out_lons, qx)
            push!(out_lats, qy)
        else
            push!(out_lons, x)
            push!(out_lats, y)
        end
    end

    # Remove consecutive duplicate points from collinear snaps
    clean_lons = Float64[]
    clean_lats = Float64[]
    for i in 1:length(out_lons)
        if isempty(clean_lons) ||
           !isapprox(clean_lons[end], out_lons[i], atol = 1e-6) ||
           !isapprox(clean_lats[end], out_lats[i], atol = 1e-6)
            push!(clean_lons, out_lons[i])
            push!(clean_lats, out_lats[i])
        end
    end

    # Ensure closed polygon ring
    if length(clean_lons) >= 3
        if clean_lons[1] != clean_lons[end] || clean_lats[1] != clean_lats[end]
            push!(clean_lons, clean_lons[1])
            push!(clean_lats, clean_lats[1])
        end
        return (clean_lons, clean_lats)
    else
        return (Float64.(poly_lons), Float64.(poly_lats))
    end
end

"""
    load_cfa_polygons(
        dir_path::AbstractString = "inputs";
        intersect_coastline::Bool = true
    ) -> Vector{NamedTuple}

Scan a directory for Crab Fishing Area (CFA) boundary polygon files matching `cfa*.dat`,
parse coordinate vectors (`lon,lat`), intersect with the regional coastline to exclude
terrestrial landmasses, and construct structured polygon definitions ready for
point-in-polygon classification and Leaflet/Makie visualizations.

# Recognized Stratum Conventions
- `cfa4x.dat` -> `"CFA 4X (Southwest NS)"`, color `"#F59E0B"` (Amber)
- `cfanorth.dat` -> `"CFA North (CFA 20-22)"`, color `"#3B82F6"` (Blue)
- `cfasouth.dat` -> `"CFA South (CFA 23-24)"`, color `"#10B981"` (Emerald)
- Custom files: Derived from filename title, with distinct palette assignment.

# Inputs
- `dir_path::AbstractString`: Directory path containing the polygon `.dat` files.
- `intersect_coastline::Bool`: Whether to ensure polygon boundaries are clipped strictly
  to marine waters excluding emergent land (default `true`).

# Outputs
- `Vector{NamedTuple}`: List of polygon objects containing `name`, `code`, `lons`, `lats`,
  `coordinates` (`[[lat, lon], ...]`), `lon` bounding range, `lat` bounding range, and `color`.
"""
function load_cfa_polygons(
    dir_path::AbstractString = "inputs";
    intersect_coastline::Bool = true
)
    if !isdir(dir_path)
        return NamedTuple[]
    end

    all_files = readdir(dir_path)
    cfa_files = filter(f -> occursin(r"^cfa.*\.dat$"i, f), all_files)
    sort!(cfa_files) # Deterministic ordering

    if isempty(cfa_files)
        return NamedTuple[]
    end

    palette = ["#3B82F6", "#10B981", "#F59E0B", "#8B5CF6", "#EC4899", "#06B6D4", "#F97316"]
    polygons = NamedTuple[]

    for (idx, filename) in enumerate(cfa_files)
        filepath = joinpath(dir_path, filename)
        stem = lowercase(replace(filename, r"\.dat$"i => ""))

        lines = readlines(filepath)
        lons = Float64[]
        lats = Float64[]

        for line in lines
            trimmed = strip(line)
            if isempty(trimmed) || startswith(trimmed, "#") || startswith(lowercase(trimmed), "lon")
                continue
            end
            parts = split(trimmed, ",")
            if length(parts) >= 2
                p_lon = tryparse(Float64, strip(parts[1]))
                p_lat = tryparse(Float64, strip(parts[2]))
                if !isnothing(p_lon) && !isnothing(p_lat)
                    push!(lons, p_lon)
                    push!(lats, p_lat)
                end
            end
        end

        if length(lons) < 3
            continue
        end

        # Optionally clip to marine domain
        if intersect_coastline
            lons, lats = intersect_polygon_with_coastline(lons, lats)
        end

        name, color, code = if stem == "cfa4x"
            ("CFA 4X (Southwest NS)", "#F59E0B", :cfa4x)
        elseif stem == "cfanorth" || stem == "cfa2022"
            ("CFA North (CFA 20-22)", "#3B82F6", :cfanorth)
        elseif stem == "cfasouth" || stem == "cfa2324"
            ("CFA South (CFA 23-24)", "#10B981", :cfasouth)
        else
            (uppercasefirst(stem), palette[mod1(idx, length(palette))], Symbol(stem))
        end

        coords = [[lats[i], lons[i]] for i in 1:length(lons)]

        push!(polygons, (
            name = name,
            code = code,
            filename = filename,
            lons = lons,
            lats = lats,
            coordinates = coords,
            lon = (minimum(lons), maximum(lons)),
            lat = (minimum(lats), maximum(lats)),
            color = color
        ))
    end

    return polygons
end

"""
    compute_empirical_connectivity(
        trajectories::NamedTuple;
        strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        strata_classifier::Union{Nothing, Function} = nothing,
        normalize_rows::Bool = true,
        empty_strata_mode::Symbol = :identity
    )

Derive the transition probability connectivity matrix \$P_{ij} = N_{i \\to j} / N_i\$
between spatial management zones (e.g. Crab Fishing Areas: CFA 20-22, CFA 23-24, CFA 4X).
Supports true polygon geometry (from `inputs/cfa*.dat`) using ray-casting point-in-polygon tests
or bounding box intervals.

# Inputs
- `trajectories::NamedTuple`: Particle tracking output.
- `strata_definitions::AbstractVector{<:NamedTuple}`: List of named zones with bounding boxes
  or explicit polygon vertex arrays `lons`, `lats`.
- `strata_classifier::Function`: Optional function `(lon, lat) -> Symbol/String/Int` returning
  stratum name or index.
- `normalize_rows::Bool`: Whether to normalize rows into conditional transition probabilities
  \$P_{ij} = N_{ij} / N_i\$ (default `true`).
- `empty_strata_mode::Symbol`: How to handle strata with zero particle releases (`:identity`, `:zero`, `:nan`).

# Outputs
- `NamedTuple`:
  - `matrix`: \$S \\times S\$ transition probability matrix \$P_{ij}\$ (or counts if `normalize_rows=false`).
  - `strata_names`: Vector of stratum names.
  - `self_retention`: Vector of self-recruitment rates (\$P_{ii}\$).
  - `export_fraction`: Vector of export fractions (\$1 - P_{ii}\$).
  - `counts_matrix`: Raw integer count of transitions \$N_{i \\to j}\$.
"""
function compute_empirical_connectivity(
    trajectories::NamedTuple;
    strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    strata_classifier::Union{Nothing, Function} = nothing,
    normalize_rows::Bool = true,
    empty_strata_mode::Symbol = :identity
)
    # 1. Resolve stratum definitions: check input, auto-detect inputs/cfa*.dat, or fallback
    defs = if !isnothing(strata_definitions)
        strata_definitions
    elseif isnothing(strata_classifier)
        auto_polys = load_cfa_polygons("inputs")
        if !isempty(auto_polys)
            vcat(auto_polys, [(name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))])
        else
            [
                (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
                (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
                (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
                (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
            ]
        end
    else
        nothing
    end

    # 2. Build classification function supporting both polygons and bounding boxes
    classify_fn = if !isnothing(strata_classifier)
        strata_classifier
    else
        (lon, lat) -> begin
            for (idx, d) in enumerate(defs)
                if hasproperty(d, :lons) && hasproperty(d, :lats) && length(d.lons) >= 3
                    if point_in_polygon(lon, lat, d.lons, d.lats)
                        return idx
                    end
                elseif hasproperty(d, :lon) && hasproperty(d, :lat)
                    if d.lon[1] <= lon <= d.lon[2] && d.lat[1] <= lat <= d.lat[2]
                        return idx
                    end
                end
            end
            return length(defs) # Default to last zone (e.g. Offshore)
        end
    end

    strata_names = if !isnothing(defs)
        [string(d.name) for d in defs]
    else
        ["Stratum 1", "Stratum 2", "Stratum 3", "Offshore"]
    end

    n_strata   = length(strata_names)
    # Weighted connectivity: each particle contributes its final survival
    # probability so P_ij reflects effective recruitment flux, not raw counts.
    has_surv_prob = hasproperty(trajectories, :survival_probability) &&
                    !isnothing(trajectories.survival_probability)
    counts          = zeros(Float64, n_strata, n_strata)
    counts_unweighted = zeros(Int, n_strata, n_strata)

    lons = trajectories.lons
    lats = trajectories.lats
    n_particles = size(lons, 1)

    for p in 1:n_particles
        lon_start, lat_start = lons[p, 1], lats[p, 1]
        lon_end, lat_end     = lons[p, end], lats[p, end]

        idx_src = classify_fn(lon_start, lat_start)
        idx_dst = classify_fn(lon_end, lat_end)

        if 1 <= idx_src <= n_strata && 1 <= idx_dst <= n_strata
            weight = has_surv_prob ?
                     Float64(trajectories.survival_probability[p, end]) : 1.0
            counts[idx_src, idx_dst]           += weight
            counts_unweighted[idx_src, idx_dst] += 1
        end
    end

    trans_prob = zeros(Float64, n_strata, n_strata)
    self_retention = zeros(Float64, n_strata)
    export_fraction = zeros(Float64, n_strata)

    for i in 1:n_strata
        row_sum = sum(counts[i, :])
        if row_sum > 0
            trans_prob[i, :] .= counts[i, :] ./ row_sum
            self_retention[i] = trans_prob[i, i]
            export_fraction[i] = 1.0 - self_retention[i]
        else
            if empty_strata_mode == :identity
                trans_prob[i, i] = 1.0
                self_retention[i] = 1.0
                export_fraction[i] = 0.0
            elseif empty_strata_mode == :nan
                trans_prob[i, :] .= NaN
                self_retention[i] = NaN
                export_fraction[i] = NaN
            else # :zero
                trans_prob[i, :] .= 0.0
                self_retention[i] = 0.0
                export_fraction[i] = 0.0
            end
        end
    end

    out_matrix = normalize_rows ? trans_prob : Float64.(counts)

    return (
        matrix = out_matrix,
        counts_matrix = counts,                   # survival-probability-weighted
        counts_unweighted = counts_unweighted,    # raw particle counts
        strata_names = strata_names,
        self_retention = self_retention,
        export_fraction = export_fraction
    )
end

"""
    export_larval_dispersal_netcdf(
        output_filename::AbstractString;
        trajectories::NamedTuple,
        model_hydrodynamics::Union{Nothing, NamedTuple} = nothing,
        lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
        lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30),
        strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    )

Export a comprehensive, CF-compliant NetCDF archive containing:
1. Gridded model velocity & turbulent diffusivity fields.
2. Gridded larval recruitment & retention metrics.
3. Gridded larval thermal exposure metrics.
4. Empirical velocity & diffusivity fields from tracked particles.
5. Macro-regional connectivity transition probability matrices.
6. Raw Lagrangian particle coordinates \$(lon, lat, depth, temp, DD)\$.

# Inputs
- `output_filename::AbstractString`: Destination path (e.g. `"outputs/larval_dispersal_analysis.nc"`).
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `model_hydrodynamics::NamedTuple`: Optional background model fields `(u, v, w, kappa_h, kappa_v)`.
- `lon_bins`, `lat_bins`: Grid bin partitions.
- `strata_definitions`: Regional zone definitions.

# Outputs
- `String`: Output file path.
"""
function export_larval_dispersal_netcdf(
    output_filename::AbstractString;
    trajectories::NamedTuple,
    model_hydrodynamics::Union{Nothing, NamedTuple} = nothing,
    lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
    lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30),
    strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    config::Union{Nothing, AbstractDict, AbstractString} = nothing
)::String
    out_dir = dirname(output_filename)
    if !isempty(out_dir) && !isdir(out_dir)
        mkpath(out_dir)
    end
    if isfile(output_filename)
        rm(output_filename, force = true)
    end

    # 1. Compute empirical movement, recruitment, thermal, and connectivity metrics
    emp_mov = estimate_empirical_movement(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    rec_met = compute_gridded_recruitment_metrics(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    therm_met = compute_gridded_thermal_metrics(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    conn_met = compute_empirical_connectivity(trajectories; strata_definitions = strata_definitions)

    lon_c = emp_mov.lon_centers
    lat_c = emp_mov.lat_centers
    n_lon = length(lon_c)
    n_lat = length(lat_c)
    n_particles, n_times = size(trajectories.lons)
    n_strata = length(conn_met.strata_names)

    # 2. Write multi-variable CF-compliant NetCDF dataset
    NCDataset(output_filename, "c") do ds
        # Define dimensions
        defDim(ds, "lon", n_lon)
        defDim(ds, "lat", n_lat)
        defDim(ds, "particle", n_particles)
        defDim(ds, "time", n_times)
        defDim(ds, "src_strata", n_strata)
        defDim(ds, "dst_strata", n_strata)

        # Coordinate variables
        v_lon = defVar(ds, "lon", Float64, ("lon",))
        v_lon.attrib["long_name"] = "Longitude"
        v_lon.attrib["units"] = "degrees_east"
        v_lon[:] = lon_c

        v_lat = defVar(ds, "lat", Float64, ("lat",))
        v_lat.attrib["long_name"] = "Latitude"
        v_lat.attrib["units"] = "degrees_north"
        v_lat[:] = lat_c

        v_time = defVar(ds, "time", Float64, ("time",))
        v_time.attrib["long_name"] = "Elapsed Simulation Time"
        v_time.attrib["units"] = "seconds"
        v_time[:] = trajectories.times

        # ── Group 1: Empirical Movement & Dispersion Fields ──────────────────
        v_u_emp = defVar(ds, "emp_u", Float64, ("lon", "lat"))
        v_u_emp.attrib["long_name"] = "Mean Empirical Zonal Velocity"
        v_u_emp.attrib["units"] = "m s-1"
        v_u_emp[:, :] = emp_mov.u_mean

        v_v_emp = defVar(ds, "emp_v", Float64, ("lon", "lat"))
        v_v_emp.attrib["long_name"] = "Mean Empirical Meridional Velocity"
        v_v_emp.attrib["units"] = "m s-1"
        v_v_emp[:, :] = emp_mov.v_mean

        v_diff_emp = defVar(ds, "emp_diffusivity", Float64, ("lon", "lat"))
        v_diff_emp.attrib["long_name"] = "Mean Empirical Turbulent Diffusivity"
        v_diff_emp.attrib["units"] = "m2 s-1"
        v_diff_emp[:, :] = emp_mov.diffusivity

        v_counts = defVar(ds, "particle_transit_count", Int32, ("lon", "lat"))
        v_counts.attrib["long_name"] = "Particle Transit Step Count"
        v_counts[:, :] = Int32.(emp_mov.sample_count)

        # ── Group 2: Gridded Recruitment & Retention Metrics ─────────────────
        v_rel = defVar(ds, "release_density", Int32, ("lon", "lat"))
        v_rel.attrib["long_name"] = "Natal Spawning Release Density"
        v_rel[:, :] = Int32.(rec_met.release_density)

        v_set = defVar(ds, "settlement_density", Int32, ("lon", "lat"))
        v_set.attrib["long_name"] = "Larval Settlement Density"
        v_set[:, :] = Int32.(rec_met.settlement_density)

        v_succ = defVar(ds, "successful_recruitment_density", Int32, ("lon", "lat"))
        v_succ.attrib["long_name"] = "Successful Benthic Recruitment Density"
        v_succ[:, :] = Int32.(rec_met.successful_settlement_density)

        v_ret = defVar(ds, "retention_rate", Float64, ("lon", "lat"))
        v_ret.attrib["long_name"] = "Cell Self-Retention Rate"
        v_ret.attrib["units"] = "fraction (0-1)"
        v_ret[:, :] = rec_met.retention_rate

        v_age = defVar(ds, "mean_settlement_age_days", Float64, ("lon", "lat"))
        v_age.attrib["long_name"] = "Mean Settlement Age of Recruits"
        v_age.attrib["units"] = "days"
        v_age[:, :] = rec_met.mean_settlement_age_days

        # ── Group 3: Gridded Thermal Exposure & Mortality ────────────────────
        v_dd = defVar(ds, "mean_degree_days", Float64, ("lon", "lat"))
        v_dd.attrib["long_name"] = "Mean In-Situ Accumulated Degree-Days"
        v_dd.attrib["units"] = "degC days"
        v_dd[:, :] = therm_met.mean_degree_days

        v_temp_exp = defVar(ds, "mean_exposure_temperature", Float64, ("lon", "lat"))
        v_temp_exp.attrib["long_name"] = "Mean Ambient Temperature Exposure"
        v_temp_exp.attrib["units"] = "degree_Celsius"
        v_temp_exp[:, :] = therm_met.mean_exposure_temperature

        v_temp_max = defVar(ds, "max_exposure_temperature", Float64, ("lon", "lat"))
        v_temp_max.attrib["long_name"] = "Maximum Ambient Temperature Exposure"
        v_temp_max.attrib["units"] = "degree_Celsius"
        v_temp_max[:, :] = therm_met.max_exposure_temperature

        v_mort = defVar(ds, "thermal_mortality_fraction", Float64, ("lon", "lat"))
        v_mort.attrib["long_name"] = "Thermal Stress Mortality Fraction"
        v_mort.attrib["units"] = "fraction (0-1)"
        v_mort[:, :] = therm_met.thermal_mortality_fraction

        # ── Group 4: Macro-Regional Connectivity Matrix ──────────────────────
        v_conn = defVar(ds, "connectivity_matrix", Float64, ("src_strata", "dst_strata"))
        v_conn.attrib["long_name"] = "Inter-Regional Transition Probability Matrix (P_ij)"
        v_conn.attrib["units"] = "probability fraction (0-1)"
        v_conn[:, :] = conn_met.matrix

        # ── Group 5: Raw Particle Trajectories ───────────────────────────────
        v_p_lon = defVar(ds, "particle_lon", Float64, ("particle", "time"))
        v_p_lon.attrib["long_name"] = "Lagrangian Particle Longitude"
        v_p_lon.attrib["units"] = "degrees_east"
        v_p_lon[:, :] = trajectories.lons

        v_p_lat = defVar(ds, "particle_lat", Float64, ("particle", "time"))
        v_p_lat.attrib["long_name"] = "Lagrangian Particle Latitude"
        v_p_lat.attrib["units"] = "degrees_north"
        v_p_lat[:, :] = trajectories.lats

        v_p_depth = defVar(ds, "particle_depth", Float64, ("particle", "time"))
        v_p_depth.attrib["long_name"] = "Lagrangian Particle Depth"
        v_p_depth.attrib["units"] = "meters"
        v_p_depth[:, :] = trajectories.depths

        v_p_temp = defVar(ds, "particle_temperature", Float64, ("particle", "time"))
        v_p_temp.attrib["long_name"] = "Ambient In-Situ Temperature"
        v_p_temp.attrib["units"] = "degree_Celsius"
        temp_matrix = hasproperty(trajectories, :temperatures) && !isnothing(trajectories.temperatures) ?
                      trajectories.temperatures : fill(4.0, n_particles, n_times)
        v_p_temp[:, :] = temp_matrix

        v_p_dd = defVar(ds, "particle_degree_days", Float64, ("particle", "time"))
        v_p_dd.attrib["long_name"] = "Cumulative Thermal Degree-Days"
        v_p_dd.attrib["units"] = "degC days"
        dd_matrix = hasproperty(trajectories, :degree_days_timeseries) && !isnothing(trajectories.degree_days_timeseries) ?
                    trajectories.degree_days_timeseries :
                    (hasproperty(trajectories, :degree_days) && trajectories.degree_days isa AbstractMatrix ?
                     trajectories.degree_days :
                     repeat(hasproperty(trajectories, :degree_days) ? trajectories.degree_days : fill(40.0, n_particles), 1, n_times))
        v_p_dd[:, :] = dd_matrix

        # Global Attributes
        ds.attrib["title"] = "Snow Crab Larval Dispersal, Empirical Movement, & Demographic Connectivity"
        ds.attrib["source"] = "Oceananigans.jl Regional Hydrodynamic Model + Lagrangian Particle Tracking"
        ds.attrib["institution"] = "Ocean Sciences & Fisheries Modeling Group"
        ds.attrib["conventions"] = "CF-1.8"
        ds.attrib["species"] = "Chionoecetes opilio (Snow Crab)"

        config_str = if !isnothing(config)
            if config isa AbstractString
                config
            else
                s_io = IOBuffer()
                TOML.print(s_io, config; sorted = true)
                String(take!(s_io))
            end
        else
            cfg = load_configuration()
            s_io = IOBuffer()
            TOML.print(s_io, cfg; sorted = true)
            String(take!(s_io))
        end
        ds.attrib["configuration"] = config_str
    end

    println("Exported comprehensive larval dispersal NetCDF to $(output_filename)")
    return output_filename
end

"""
    export_larval_dispersal_jld2(
        output_filename::AbstractString;
        trajectories::NamedTuple,
        lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
        lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30),
        strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        config::Union{Nothing, AbstractDict, AbstractString} = nothing
    )

Save full empirical movement, recruitment, thermal, connectivity diagnostics,
and simulation configuration metadata to JLD2.
"""
function export_larval_dispersal_jld2(
    output_filename::AbstractString;
    trajectories::NamedTuple,
    lon_bins::AbstractVector{<:Real} = range(-68.0, -57.0, length = 30),
    lat_bins::AbstractVector{<:Real} = range(42.0, 47.0, length = 30),
    strata_definitions::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    config::Union{Nothing, AbstractDict, AbstractString} = nothing
)::String
    out_dir = dirname(output_filename)
    if !isempty(out_dir) && !isdir(out_dir)
        mkpath(out_dir)
    end

    emp_mov = estimate_empirical_movement(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    rec_met = compute_gridded_recruitment_metrics(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    therm_met = compute_gridded_thermal_metrics(trajectories; lon_bins = lon_bins, lat_bins = lat_bins)
    conn_met = compute_empirical_connectivity(trajectories; strata_definitions = strata_definitions)

    cfg_dict = if !isnothing(config)
        if config isa AbstractDict
            config
        else
            TOML.parse(config)
        end
    else
        load_configuration()
    end

    JLD2.jldsave(output_filename;
        trajectories = trajectories,
        empirical_movement = emp_mov,
        recruitment_metrics = rec_met,
        thermal_metrics = therm_met,
        connectivity = conn_met,
        configuration = cfg_dict
    )

    println("Exported larval dispersal JLD2 bundle to $(output_filename)")
    return output_filename
end
