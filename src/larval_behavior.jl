"""
    larval_behavior.jl

Snow crab (*Chionoecetes opilio*) larval behavior, Diel Vertical Migration (DVM),
thermal degree-day ontogenetic molting, tidal current superposition, and
benthic nursery settlement suitability filtering.
"""

using Random

"""
    diel_vertical_migration_velocity(
        z::Real,
        t::Real;
        stage::Symbol = :zoea1,
        day_depth::Real = -50.0,
        night_depth::Real = -10.0,
        max_swim_speed::Real = 0.005,
        relaxation_depth::Real = 10.0
    )

Compute active vertical swimming velocity \$w_{\\text{swim}}(z, t)\$ for snow crab
larvae undergoing stage-dependent Diel Vertical Migration (DVM).

# Mathematical & Ecological Model
Snow crab (*Chionoecetes opilio*) larvae exhibit stage-dependent vertical
positioning (Incze et al., 1987; Lovrich et al., 1995; Sainte-Marie & Sainte-Marie, 1999).
Planktonic Zoea I and II stages occupy shallower depths at night and migrate to deeper,
colder layers during the day. Semiplanktonic Megalopa larvae actively seek benthic
nursery grounds prior to settlement.

The instantaneous target depth \$z_{\\text{target}}(t)\$ over a 24-hour diurnal
cycle (\$T = 86400\\text{ s}\$) is:
```math
z_{\\text{target}}(t) = \\frac{z_{\\text{day}} + z_{\\text{night}}}{2}
- \\frac{z_{\\text{day}} - z_{\\text{night}}}{2} \\cos\\left( \\frac{2\\pi t}{86400} \\right)
```
The active swimming velocity \$w_{\\text{swim}}(z, t)\$ is:
```math
w_{\\text{swim}}(z, t) = w_{\\max} \\tanh\\left( \\frac{z_{\\text{target}}(t) - z}{L_{\\text{relax}}} \\right)
```

# Inputs
- `z::Real`: Current vertical coordinate (depth in meters, negative downward).
- `t::Real`: Current simulation time in seconds.
- `stage::Symbol`: Larval development stage (`:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`).
- `day_depth::Real`: Daytime preferred depth in meters (default -50.0 m).
- `night_depth::Real`: Nighttime preferred depth in meters (default -10.0 m).
- `megalopa_day_depth::Real`: Megalopa daytime benthic seeking depth in meters (default -120.0 m).
- `megalopa_night_depth::Real`: Megalopa nighttime depth in meters (default -60.0 m).
- `zoea2_depth_factor::Real`: Deepening factor for Zoea II relative to Zoea I (default 1.2).
- `megalopa_swim_factor::Real`: Swimming speed multiplier for Megalopa (default 1.5).
- `max_swim_speed::Real`: Maximum swimming speed in \$m s^{-1}\$ (default 0.005 m/s ≈ 5 mm/s).
- `relaxation_depth::Real`: Vertical transition scale \$L_{\\text{relax}}\$ in meters.

# Outputs
- `Float64`: Vertical swimming velocity \$w_{\\text{swim}}\$ in \$m s^{-1}\$.

# References
- Incze, L. S., Armstrong, D. A., & Smith, S. L. (1987). Abundance of filter-feeding
  and pelagic stages of crab larvae in the southeastern Bering Sea.
  *Marine Biology*, 95(2), 195-200. DOI: 10.1007/BF00409006
- Lovrich, G. A., Sainte-Marie, B., & Smith, B. D. (1995). Depth distribution and
  seasonal movements of *Chionoecetes opilio* in Baie Sainte-Marguerite.
  *Canadian Journal of Fisheries and Aquatic Sciences*, 52(4), 903-913.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). Growth, developmental stages, and
  vertical distribution of snow crab larvae (*Chionoecetes opilio*).
  *Canadian Journal of Fisheries and Aquatic Sciences*, 56(11), 2181-2193.
"""
function diel_vertical_migration_velocity(
    z::Real,
    t::Real;
    stage::Symbol = :zoea1,
    day_depth::Real = -50.0,
    night_depth::Real = -10.0,
    megalopa_day_depth::Real = -120.0,
    megalopa_night_depth::Real = -60.0,
    zoea2_depth_factor::Real = 1.2,
    megalopa_swim_factor::Real = 1.5,
    max_swim_speed::Real = 0.005,
    relaxation_depth::Real = 10.0
)
    if stage == :instar1_settled || stage == :dead
        return 0.0
    end

    z_day, z_night, w_max = if stage == :zoea1
        (Float64(day_depth), Float64(night_depth), Float64(max_swim_speed))
    elseif stage == :zoea2
        # Zoea II deepens preferentially during the day (ontogenetic descent toward
        # the CIL boundary) while nighttime ascent is only slightly greater than Zoea I.
        # Day depth factor 1.3, night depth factor 1.05 (Sainte-Marie & Sainte-Marie 1999).
        (Float64(day_depth   * 1.3),
         Float64(night_depth * 1.05),
         Float64(max_swim_speed * zoea2_depth_factor))
    elseif stage == :megalopa
        (Float64(megalopa_day_depth), Float64(megalopa_night_depth),
         Float64(max_swim_speed * megalopa_swim_factor))
    else
        (Float64(day_depth), Float64(night_depth), Float64(max_swim_speed))
    end

    diurnal_phase = 2.0 * π * (mod(t, 86400.0)) / 86400.0
    z_mean = 0.5 * (z_day + z_night)
    # Amplitude is positive: z_day < z_night (deeper during day).
    # At phase=0 (noon): z_target = z_mean - z_amp = z_day   (deep) ✓
    # At phase=π (midnight): z_target = z_mean + z_amp = z_night (shallow) ✓
    z_amp    = 0.5 * (z_day - z_night)
    z_target = z_mean - z_amp * cos(diurnal_phase)

    dz = z_target - z
    w_swim = w_max * tanh(dz / relaxation_depth)
    return Float64(w_swim)
end

"""
    superpose_tidal_velocity(
        u_mean::Real,
        v_mean::Real,
        t::Real;
        u_amp::Real = 0.25,
        v_amp::Real = 0.12,
        period::Real = 44712.0,
        phase::Real = 0.0
    )

Superpose semi-diurnal \$M_2\$ tidal current oscillations onto mean hydrodynamic
background flow.

# Mathematical Formulation
```math
u_{\\text{total}}(t) = u_{\\text{mean}} + u_{\\text{amp}} \\cos\\left( \\frac{2\\pi t}{T_{M2}} + \\phi \\right)
```
```math
v_{\\text{total}}(t) = v_{\\text{mean}} + v_{\\text{amp}} \\sin\\left( \\frac{2\\pi t}{T_{M2}} + \\phi \\right)
```
where \$T_{M2} \\approx 12.42\\text{ hours} = 44712\\text{ s}\$.

# Inputs
- `u_mean, v_mean::Real`: Background Eulerian flow velocities in \$m s^{-1}\$.
- `t::Real`: Current time in seconds.
- `u_amp, v_amp::Real`: Semi-major and semi-minor tidal velocity amplitudes (m s⁻¹).
- `period::Real`: Tidal constituent period in seconds (default 44712.0 s for \$M_2\$).
- `phase::Real`: Initial phase offset in radians.

# Outputs
- `Tuple{Float64, Float64}`: Total velocities `(u_tot, v_tot)`.
"""
function superpose_tidal_velocity(
    u_mean::Real,
    v_mean::Real,
    t::Real;
    u_amp::Real = 0.25,
    v_amp::Real = 0.12,
    period::Real = 44712.0,
    phase::Real = 0.0
)
    omega = 2.0 * π / period
    u_tide = u_amp * cos(omega * t + phase)
    v_tide = v_amp * sin(omega * t + phase)
    return (Float64(u_mean + u_tide), Float64(v_mean + v_tide))
end

"""
    update_larval_stage(
        current_stage::Symbol,
        degree_days::Real;
        dd_zoea1_to_zoea2::Real = 150.0,
        dd_zoea2_to_megalopa::Real = 310.0,
        dd_megalopa_to_settle::Real = 510.0
    )

Determine the current ontogenetic development stage based on cumulative thermal degree-days.

# Mathematical Formulation
Degree-days (\$DD = \\int \\max(0, T - T_0) dt\$) accumulate in situ. Molting occurs
at stage-specific thermal thresholds calibrated from rearing experiments at 2°C
(Kuhn & Choi, 2011; Sainte-Marie & Sainte-Marie, 1999):
- \$0 \\le DD < 65\$: **Zoea I** (~32 days at T=2°C)
- \$65 \\le DD < 130\$: **Zoea II** (~65 days cumulative)
- \$130 \\le DD < 200\$: **Megalopa** (~100 days cumulative)
- \$DD \\ge 200\$: **Instar I (Settled)** (~130 days cumulative)

# Inputs
- `current_stage::Symbol`: Present developmental stage.
- `degree_days::Real`: Accumulated degree-days (°C · days).

# Outputs
- `Symbol`: Updated stage (`:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`).
"""
function update_larval_stage(
    current_stage::Symbol,
    degree_days::Real;
    dd_zoea1_to_zoea2::Real = 65.0,
    dd_zoea2_to_megalopa::Real = 130.0,
    dd_megalopa_to_settle::Real = 200.0
)
    if current_stage == :dead || current_stage == :instar1_settled
        return current_stage
    end

    if degree_days < dd_zoea1_to_zoea2
        return :zoea1
    elseif degree_days < dd_zoea2_to_megalopa
        return :zoea2
    elseif degree_days < dd_megalopa_to_settle
        return :megalopa
    else
        return :instar1_settled
    end
end

"""
    evaluate_settlement_suitability(
        bed_elevation::Real,
        bottom_temperature::Real;
        min_depth::Real = -250.0,
        max_depth::Real = -50.0,
        max_bottom_temp::Real = 6.0
    )

Evaluate whether a settling Megalopa/Instar I larva has arrived on suitable benthic
nursery grounds on the continental shelf.

# Criteria
1. **Depth Criterion**: Seafloor elevation must lie within shelf nursery limits
   (\$-250\\text{ m} \\le z_{\\text{bed}} \\le -50\\text{ m}\$).
2. **Thermal Habitat**: Bottom temperature must remain within the Cold Intermediate
   Layer (CIL) thermal window (\$T_{\\text{bottom}} \\le 6.0^\\circ\\text{C}\$).

# Inputs
- `bed_elevation::Real`: Seafloor elevation in meters (negative).
- `bottom_temperature::Real`: Bottom water temperature in °C.

# Outputs
- `NamedTuple`: `(suitable::Bool, reason::String)`
"""
function evaluate_settlement_suitability(
    bed_elevation::Real,
    bottom_temperature::Real;
    min_depth::Real = -250.0,
    max_depth::Real = -50.0,
    max_bottom_temp::Real = 6.0
)
    if bed_elevation > max_depth
        return (suitable = false, reason = "Too shallow / nearshore surf zone (z > $(max_depth)m)")
    elseif bed_elevation < min_depth
        return (suitable = false, reason = "Too deep / offshore abyssal slope (z < $(min_depth)m)")
    elseif bottom_temperature > max_bottom_temp
        return (suitable = false, reason = "Thermal stress / bottom temp > $(max_bottom_temp)°C")
    else
        return (suitable = true, reason = "Suitable cold-water shelf nursery ground")
    end
end

"""
    larval_transport_step(
        x::Real,
        y::Real,
        z::Real,
        u::Real,
        v::Real,
        w::Real,
        κ_h::Real,
        κ_v::Real,
        dt::Real;
        t::Real = 0.0,
        stage::Symbol = :zoea1,
        z_bottom::Real = -1000.0,
        z_surface::Real = 0.0,
        is_lat_lon::Bool = true,
        rng::AbstractRNG = Random.default_rng()
    )

Perform a single stochastic Lagrangian transport step combining hydrodynamic
advection, turbulent diffusion, and active larval vertical swimming.

# Inputs
- `x, y, z::Real`: Current coordinates of the larva (longitude/x, latitude/y, depth).
- `u, v, w::Real`: Fluid velocity components at particle position in \$m s^{-1}\$.
- `κ_h::Real`: Horizontal eddy diffusivity in \$m^2 s^{-1}\$.
- `κ_v::Real`: Vertical eddy diffusivity in \$m^2 s^{-1}\$.
- `dt::Real`: Time step in seconds.
- `t::Real`: Current simulation time in seconds.
- `stage::Symbol`: Larval stage (`:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`).
- `z_bottom::Real`: Local seabed elevation in meters.
- `z_surface::Real`: Sea surface elevation (default 0.0 m).
- `is_lat_lon::Bool`: If true, converts metric displacements into spherical degrees.
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `Tuple{Float64, Float64, Float64}`: Updated coordinates `(x_new, y_new, z_new)`.
"""
function larval_transport_step(
    x::Real,
    y::Real,
    z::Real,
    u::Real,
    v::Real,
    w::Real,
    κ_h::Real,
    κ_v::Real,
    dt::Real;
    t::Real = 0.0,
    stage::Symbol = :zoea1,
    z_bottom::Real = -1000.0,
    z_surface::Real = 0.0,
    is_lat_lon::Bool = true,
    rng::AbstractRNG = Random.default_rng()
)
    # If settled on bottom or dead, particle remains at seabed
    if stage == :instar1_settled
        return (Float64(x), Float64(y), Float64(z_bottom))
    end

    w_swim = diel_vertical_migration_velocity(z, t; stage = stage)

    ξ_x = randn(rng)
    ξ_y = randn(rng)
    ξ_z = randn(rng)

    dx_meters = u * dt + sqrt(max(0.0, 2.0 * κ_h * dt)) * ξ_x
    dy_meters = v * dt + sqrt(max(0.0, 2.0 * κ_h * dt)) * ξ_y
    dz_meters = (w + w_swim) * dt + sqrt(max(0.0, 2.0 * κ_v * dt)) * ξ_z

    if is_lat_lon
        r_earth = 6.371e6
        lat_rad = deg2rad(clamp(Float64(y), -89.9, 89.9))
        deg_per_meter_lat = 180.0 / (π * r_earth)
        deg_per_meter_lon = 180.0 / (π * r_earth * cos(lat_rad))

        x_new = x + dx_meters * deg_per_meter_lon
        y_new = y + dy_meters * deg_per_meter_lat
    else
        x_new = x + dx_meters
        y_new = y + dy_meters
    end
    z_new = z + dz_meters

    # Absorbing vertical boundaries: larvae cannot exit surface or seabed.
    # Elastic reflection at the bottom is unphysical — it would impart an upward
    # velocity impulse. Instead, clamp to the boundary and let settlement
    # evaluation handle benthic contact.
    z_new = clamp(z_new, z_bottom, z_surface)

    # Land rejection: polygon-based ray-casting via is_point_on_land.
    # Priority: (1) accept new position, (2) reflect displacement and retry,
    # (3) stick to previous position (absorbing coastline).
    if is_lat_lon && is_point_on_land(x_new, y_new)
        # Specular reflection: reverse the net horizontal displacement vector
        x_refl = x - dx_meters * deg_per_meter_lon
        y_refl = y - dy_meters * deg_per_meter_lat
        if is_point_on_land(x_refl, y_refl)
            # Both candidate positions are on land: remain at current location
            x_new = Float64(x)
            y_new = Float64(y)
        else
            x_new = x_refl
            y_new = y_refl
        end
    end

    return (Float64(x_new), Float64(y_new), Float64(z_new))
end

"""
    initialize_larval_particles(
        n_particles::Int;
        lon_range::Tuple{Real, Real} = (-64.0, -60.0),
        lat_range::Tuple{Real, Real} = (43.5, 45.5),
        depth_range::Tuple{Real, Real} = (-60.0, -20.0),
        min_seabed_depth::Real = 100.0,
        buffer_km::Real = 0.0,
        bathymetry = nothing,
        coastline = nothing,
        stratum::Union{Nothing, NamedTuple, Symbol, AbstractString} = nothing,
        stage::Symbol = :zoea1,
        max_attempts::Int = 100000,
        rng::AbstractRNG = Random.default_rng()
    )

Initialize an ensemble of snow crab larval particle states across a specified
spawning area, spatial stratum, or extended buffered marine envelope on the
continental shelf. Strict marine water and high-resolution coastline land exclusion
constraints are enforced so that larvae are placed strictly in active marine water
columns and never on emergent land.

# Spatial Buffering & Land Exclusion
1. **User-Defined Buffer Envelope**: When `buffer_km > 0.0` (e.g. default 100 km), the release
   domain is expanded outward by \$\\Delta\\phi = d_{\\text{buf}}/R_{\\text{earth}}\$ and
   \$\\Delta\\lambda = \\Delta\\phi / \\cos(\\phi_0)\$ to capture inflow, upstream advection,
   and transboundary exchange across administrative boundaries.
2. **Stratum & CFA Polygon Seeding**: When `stratum` is provided (e.g. `:cfa4x`, `:cfasouth`,
   `:cfanorth`, or a polygon `NamedTuple`), particles are sampled across the stratum or its
   buffered marine envelope.
3. **Strict Coastline Land Rejection**: For every candidate coordinate, `!is_point_on_land(lon, lat)`
   and \$z_{\\text{bed}} \\le -h_{\\text{min}}\$ are strictly enforced.

# Inputs
- `n_particles::Int`: Total number of particles to release.
- `lon_range::Tuple{Real, Real}`: Longitude release bounds (used if `stratum = nothing`).
- `lat_range::Tuple{Real, Real}`: Latitude release bounds (used if `stratum = nothing`).
- `depth_range::Tuple{Real, Real}`: Depth release bounds in meters (negative downward).
- `min_seabed_depth::Real`: Minimum water depth in meters (default 100.0 m, min 0.0 m).
- `buffer_km::Real`: Spatial buffer distance in km beyond stratum/domain bounds (default 0.0 km).
- `bathymetry`: Bathymetry source (`NamedTuple`, NetCDF path, elevation `Function`, or `nothing`).
- `coastline`: Optional coastline multi-polygon definitions.
- `stratum`: Optional target stratum (`NamedTuple` with `:lons, :lats`, or CFA symbol `:cfa4x`, etc.).
- `stage::Symbol`: Initial stage (`:zoea1`, `:zoea2`, `:megalopa`).
- `max_attempts::Int`: Maximum rejection sampling attempts per particle (default 100000).
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `NamedTuple`: Arrays `(lon, lat, depth, stage, id, degree_days, alive, settlement_status)`.
"""
function initialize_larval_particles(
    n_particles::Int;
    lon_range::Tuple{Real, Real} = (-64.0, -60.0),
    lat_range::Tuple{Real, Real} = (43.5, 45.5),
    depth_range::Tuple{Real, Real} = (-60.0, -20.0),
    min_seabed_depth::Real = 100.0,
    buffer_km::Real = 0.0,
    bathymetry::Union{Nothing, Function, NamedTuple, AbstractString} = nothing,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    stratum::Union{Nothing, NamedTuple, Symbol, AbstractString} = nothing,
    stage::Symbol = :zoea1,
    max_attempts::Int = 100000,
    rng::AbstractRNG = Random.default_rng()
)
    if n_particles <= 0
        error("Number of particles must be positive: $(n_particles)")
    end

    # Auto-resolve bathymetry if not explicitly provided
    resolved_bathy = if isnothing(bathymetry)
        def_path = "inputs/bathymetry_active.nc"
        isfile(def_path) ? def_path : nothing
    else
        bathymetry
    end

    # Auto-resolve stratum polygon if requested by name/symbol
    resolved_stratum_poly = if stratum isa NamedTuple && hasproperty(stratum, :lons) && hasproperty(stratum, :lats)
        stratum
    elseif stratum isa Union{Symbol, AbstractString}
        cfa_list = load_cfa_polygons("inputs")
        s_code = Symbol(lowercase(string(stratum)))
        found = findfirst(p -> p.code == s_code, cfa_list)
        if isnothing(found)
            error("Stratum '$(stratum)' not found in loaded CFA polygons.")
        end
        cfa_list[found]
    else
        nothing
    end

    # Determine base bounding box
    raw_lon_range, raw_lat_range = if !isnothing(resolved_stratum_poly)
        (extrema(resolved_stratum_poly.lons), extrema(resolved_stratum_poly.lats))
    else
        (lon_range, lat_range)
    end

    # Apply spatial buffer if buffer_km > 0
    eff_lon_range, eff_lat_range = if buffer_km > 0.0
        expand_domain_with_buffer(raw_lon_range, raw_lat_range, buffer_km = buffer_km)
    else
        (raw_lon_range, raw_lat_range)
    end

    lons = Vector{Float64}(undef, n_particles)
    lats = Vector{Float64}(undef, n_particles)
    depths = Vector{Float64}(undef, n_particles)

    h_min = max(0.0, Float64(min_seabed_depth))

    # Fast deterministic marine grid sampling if dataset is available and no strict unbuffered polygon
    use_grid_sampling = (resolved_bathy isa Union{NamedTuple, AbstractString}) &&
                        (isnothing(resolved_stratum_poly) || buffer_km > 0.0)

    if use_grid_sampling && isnothing(resolved_stratum_poly)
        sampled_lons, sampled_lats, sampled_zbed = sample_marine_coordinates(
            n_particles,
            resolved_bathy,
            lon_range = eff_lon_range,
            lat_range = eff_lat_range,
            min_seabed_depth = h_min,
            coastline = coastline,
            rng = rng
        )
        lons .= sampled_lons
        lats .= sampled_lats

        for p in 1:n_particles
            cand_z = depth_range[1] + rand(rng, Float64) * (depth_range[2] - depth_range[1])
            depths[p] = clamp(cand_z, sampled_zbed[p] + 2.0, -1.0)
        end
    else
        bathy_fn = if resolved_bathy isa Function
            resolved_bathy
        elseif resolved_bathy isa Union{NamedTuple, AbstractString}
            get_bathymetry_interpolator(resolved_bathy)
        else
            nothing
        end

        total_attempts = 0
        accepted = 0
        max_total_attempts = max_attempts * n_particles

        while accepted < n_particles
            total_attempts += 1
            if total_attempts > max_total_attempts
                target_desc = !isnothing(resolved_stratum_poly) ?
                              "stratum '$(resolved_stratum_poly.name)' (buffer = $(buffer_km) km)" :
                              "longitude $(eff_lon_range) and latitude $(eff_lat_range)"
                error(
                    "Unable to place $(n_particles) larvae with min_seabed_depth = $(min_seabed_depth) m " *
                    "within $(target_desc) after $(total_attempts) attempts. " *
                    "Please verify that the target area contains sufficient marine water depth."
                )
            end

            cand_lon = eff_lon_range[1] + rand(rng, Float64) * (eff_lon_range[2] - eff_lon_range[1])
            cand_lat = eff_lat_range[1] + rand(rng, Float64) * (eff_lat_range[2] - eff_lat_range[1])

            # 1. Stratum polygon enclosure check (if specified and buffer_km == 0)
            if !isnothing(resolved_stratum_poly) && buffer_km <= 0.0
                if !point_in_polygon(cand_lon, cand_lat, resolved_stratum_poly.lons, resolved_stratum_poly.lats)
                    continue
                end
            end

            # 2. Strict coastline land rejection
            if is_point_on_land(cand_lon, cand_lat; coastline = coastline)
                continue
            end

            # 3. Bathymetric elevation constraint
            z_bed = !isnothing(bathy_fn) ? bathy_fn(cand_lon, cand_lat) : -h_min - 10.0
            if z_bed > -h_min || z_bed >= 0.0
                continue
            end

            cand_z = depth_range[1] + rand(rng, Float64) * (depth_range[2] - depth_range[1])
            cand_z = clamp(cand_z, z_bed + 2.0, -1.0)

            accepted += 1
            lons[accepted] = cand_lon
            lats[accepted] = cand_lat
            depths[accepted] = cand_z
        end
    end

    stages = fill(stage, n_particles)
    ids = collect(1:n_particles)
    degree_days = zeros(Float64, n_particles)
    alive = fill(true, n_particles)
    settlement_status = fill(:pelagic, n_particles)

    return (
        lon = lons,
        lat = lats,
        depth = depths,
        stage = stages,
        id = ids,
        degree_days = degree_days,
        alive = alive,
        settlement_status = settlement_status
    )
end

"""
    track_larval_cohort(
        larvae::NamedTuple;
        velocity_fn::Function,
        temperature_fn::Union{Nothing, Function} = nothing,
        bathymetry_fn::Union{Nothing, Function} = nothing,
        total_duration::Real = 86400.0,
        dt::Real = 300.0,
        κ_h::Real = 10.0,
        κ_v::Real = 1e-4,
        is_lat_lon::Bool = true,
        enable_tides::Bool = false,
        enable_molting::Bool = true,
        rng::AbstractRNG = Random.default_rng()
    )

Simulate Lagrangian trajectories for an ensemble cohort of snow crab larvae,
integrating in-situ degree-day accumulation, stage molting, thermal mortality,
and benthic settlement suitability.

# Inputs
- `larvae::NamedTuple`: Initial particle ensemble from `initialize_larval_particles`.
- `velocity_fn::Function`: Function `(lon, lat, z, t) -> (u, v, w)` returning flow velocities.
- `temperature_fn::Union{Nothing, Function}`: Optional `(lon, lat, z, t) -> T_celsius`.
- `bathymetry_fn::Union{Nothing, Function}`: Optional `(lon, lat) -> z_bottom`.
- `total_duration::Real`: Total tracking duration in seconds.
- `dt::Real`: Lagrangian integration time step in seconds.
- `κ_h::Real`: Horizontal turbulent eddy diffusivity in \$m^2 s^{-1}\$.
- `κ_v::Real`: Vertical turbulent eddy diffusivity in \$m^2 s^{-1}\$.
- `is_lat_lon::Bool`: Coordinate format flag.
- `enable_tides::Bool`: Whether to superpose \$M_2\$ tidal currents.
- `enable_molting::Bool`: Whether to update stages via thermal degree-days.
- `default_temperature::Real`: Fallback temperature if `temperature_fn = nothing` (default 4.0 °C).
- `default_bottom_depth::Real`: Fallback seabed depth if `bathymetry_fn = nothing` (default -1000.0 m).
- `t_base::Real`: Baseline development threshold temperature \$T_0\$ for degree-days (default 0.0 °C).
- `dd_zoea1_to_zoea2::Real`: Degree-day threshold for Zoea I -> Zoea II molt (default 150.0).
- `dd_zoea2_to_megalopa::Real`: Degree-day threshold for Zoea II -> Megalopa molt (default 310.0).
- `dd_megalopa_to_settle::Real`: Degree-day threshold for Megalopa -> Instar I settlement (default 510.0).
- `mortality_base::Real`: Baseline daily natural mortality rate (default 0.02 day⁻¹).
- `mortality_thermal_threshold::Real`: Thermal stress mortality threshold in °C (default 10.0 °C).
- `mortality_thermal_sensitivity::Real`: Thermal mortality multiplier (default 0.015).
- `settlement_min_depth::Real`: Deepest boundary for nursery settlement habitat (default -250.0 m).
- `settlement_max_depth::Real`: Shallowest boundary for nursery settlement habitat (default -50.0 m).
- `settlement_max_temp::Real`: Maximum bottom temperature for suitable settlement (default 6.0 °C).
- `tidal_u_amp::Real`: Zonal tidal velocity amplitude in m/s (default 0.25 m/s).
- `tidal_v_amp::Real`: Meridional tidal velocity amplitude in m/s (default 0.12 m/s).
- `tidal_period::Real`: Tidal constituent period in seconds (default 44712.0 s).
- `tidal_phase::Real`: Tidal initial phase offset in radians (default 0.0).
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `NamedTuple`: Trajectory record containing `(lons, lats, depths, stages, degree_days, alive, settlement_status, times, ids)`.
"""
function track_larval_cohort(
    larvae::NamedTuple;
    velocity_fn::Function,
    temperature_fn::Union{Nothing, Function} = nothing,
    bathymetry_fn::Union{Nothing, Function} = nothing,
    total_duration::Real = 86400.0,
    dt::Real = 300.0,
    κ_h::Real = 10.0,
    κ_v::Real = 1e-4,
    is_lat_lon::Bool = true,
    enable_tides::Bool = false,
    enable_molting::Bool = true,
    default_temperature::Real = 4.0,
    default_bottom_depth::Real = -1000.0,
    t_base::Real = -1.5,
    dd_zoea1_to_zoea2::Real = 65.0,
    dd_zoea2_to_megalopa::Real = 130.0,
    dd_megalopa_to_settle::Real = 200.0,
    mortality_base::Real = 0.02,
    mortality_thermal_threshold::Real = 7.0,
    mortality_thermal_sensitivity::Real = 0.015,
    mortality_cold_threshold::Real = -1.5,
    mortality_cold_sensitivity::Real = 0.01,
    settlement_min_depth::Real = -250.0,
    settlement_max_depth::Real = -50.0,
    settlement_max_temp::Real = 6.0,
    tidal_u_amp::Real = 0.25,
    tidal_v_amp::Real = 0.12,
    tidal_period::Real = 44712.0,
    tidal_phase::Real = 0.0,
    min_survival_prob::Real = 0.01,
    rng::AbstractRNG = Random.default_rng()
)
    n_particles = length(larvae.lon)
    n_steps = Int(floor(total_duration / dt)) + 1
    times = range(0.0, total_duration, length = n_steps)
    dt_days = dt / 86400.0

    traj_lon = Matrix{Float64}(undef, n_particles, n_steps)
    traj_lat = Matrix{Float64}(undef, n_particles, n_steps)
    traj_depth = Matrix{Float64}(undef, n_particles, n_steps)
    traj_temp = Matrix{Float64}(undef, n_particles, n_steps)
    traj_dd = Matrix{Float64}(undef, n_particles, n_steps)
    traj_surv = Matrix{Float64}(undef, n_particles, n_steps)
    traj_stage = Matrix{Symbol}(undef, n_particles, n_steps)

    current_degree_days = hasproperty(larvae, :degree_days) ?
                          copy(larvae.degree_days) : zeros(Float64, n_particles)
    current_alive = hasproperty(larvae, :alive) ?
                    copy(larvae.alive) : fill(true, n_particles)
    current_settlement = hasproperty(larvae, :settlement_status) ?
                         copy(larvae.settlement_status) : fill(:pelagic, n_particles)
    current_settlement_age = fill(Float64(total_duration), n_particles)

    # Set initial states at t = 0
    traj_lon[:, 1] = copy(larvae.lon)
    traj_lat[:, 1] = copy(larvae.lat)
    traj_depth[:, 1] = copy(larvae.depth)
    traj_stage[:, 1] = copy(larvae.stage)
    traj_dd[:, 1] = copy(current_degree_days)
    traj_surv[:, 1] .= 1.0

    for p in 1:n_particles
        init_T = isnothing(temperature_fn) ? Float64(default_temperature) :
                 temperature_fn(traj_lon[p, 1], traj_lat[p, 1], traj_depth[p, 1], 0.0)
        traj_temp[p, 1] = init_T
    end

    for s in 1:(n_steps - 1)
        t_current = times[s]
        for p in 1:n_particles
            if !current_alive[p]
                traj_lon[p, s + 1] = traj_lon[p, s]
                traj_lat[p, s + 1] = traj_lat[p, s]
                traj_depth[p, s + 1] = traj_depth[p, s]
                traj_temp[p, s + 1] = traj_temp[p, s]
                traj_dd[p, s + 1] = traj_dd[p, s]
                traj_surv[p, s + 1] = 0.0
                traj_stage[p, s + 1] = :dead
                continue
            end

            cur_lon = traj_lon[p, s]
            cur_lat = traj_lat[p, s]
            cur_depth = traj_depth[p, s]
            cur_stage = traj_stage[p, s]

            # Query background velocity
            u, v, w = velocity_fn(cur_lon, cur_lat, cur_depth, t_current)
            if enable_tides
                u, v = superpose_tidal_velocity(
                    u, v, t_current,
                    u_amp = tidal_u_amp,
                    v_amp = tidal_v_amp,
                    period = tidal_period,
                    phase = tidal_phase
                )
            end

            # Query temperature and seabed depth
            cur_T = isnothing(temperature_fn) ? Float64(default_temperature) :
                    temperature_fn(cur_lon, cur_lat, cur_depth, t_current)
            z_bed = isnothing(bathymetry_fn) ? Float64(default_bottom_depth) :
                    bathymetry_fn(cur_lon, cur_lat)

            # Degree-day accumulation, survival rate, and molting
            mort_rate = larval_thermal_mortality_rate(
                cur_T,
                base_mortality       = mortality_base,
                thermal_threshold    = mortality_thermal_threshold,
                thermal_sensitivity  = mortality_thermal_sensitivity,
                cold_threshold       = mortality_cold_threshold,
                cold_sensitivity     = mortality_cold_sensitivity
            )
            traj_surv[p, s + 1] = max(0.0, traj_surv[p, s] * exp(-mort_rate * dt_days))
            traj_temp[p, s + 1] = cur_T

            if enable_molting && cur_stage != :instar1_settled
                current_degree_days[p] += max(0.0, cur_T - t_base) * dt_days
                new_stage = update_larval_stage(
                    cur_stage,
                    current_degree_days[p],
                    dd_zoea1_to_zoea2 = dd_zoea1_to_zoea2,
                    dd_zoea2_to_megalopa = dd_zoea2_to_megalopa,
                    dd_megalopa_to_settle = dd_megalopa_to_settle
                )

                # Check settlement if molted to Instar I
                if new_stage == :instar1_settled && cur_stage != :instar1_settled
                    suit = evaluate_settlement_suitability(
                        z_bed,
                        cur_T,
                        min_depth = settlement_min_depth,
                        max_depth = settlement_max_depth,
                        max_bottom_temp = settlement_max_temp
                    )
                    if suit.suitable
                        current_settlement[p] = :settled_successful
                        current_settlement_age[p] = t_current
                    else
                        # Unsuitable habitat: mark status but keep larva alive so
                        # it can continue drifting to potentially settle elsewhere.
                        # Real C. opilio megalopa can re-enter pelagic phase
                        # (Sainte-Marie & Sainte-Marie, 1999).
                        current_settlement[p] = :settled_unsuitable
                        current_settlement_age[p] = t_current
                        # Revert stage to megalopa so particle keeps moving
                        new_stage = :megalopa
                    end
                end
                cur_stage = new_stage
            end
            traj_dd[p, s + 1] = current_degree_days[p]

            # Step Lagrangian transport
            new_lon, new_lat, new_depth = larval_transport_step(
                cur_lon, cur_lat, cur_depth,
                u, v, w, κ_h, κ_v, dt,
                t = t_current, stage = cur_stage,
                z_bottom = z_bed,
                is_lat_lon = is_lat_lon, rng = rng
            )

            # Coastal shoreline boundary condition: absorbing boundary.
            # Larvae that step onto land remain at their previous marine position.
            # Note: larval_transport_step() already attempts specular reflection
            # inside itself; this outer check catches cases where the post-step
            # position still lands on the coast.
            if !isnothing(bathymetry_fn)
                z_bed_new = bathymetry_fn(new_lon, new_lat)
                if z_bed_new >= 0.0
                    new_lon = cur_lon
                    new_lat = cur_lat
                end
            end

            traj_lon[p, s + 1] = new_lon
            traj_lat[p, s + 1] = new_lat
            traj_depth[p, s + 1] = new_depth
            traj_stage[p, s + 1] = cur_stage

            # Mark as dead when cumulative survival probability falls below threshold
            if traj_surv[p, s + 1] < min_survival_prob
                current_alive[p] = false
            end
        end
    end

    return (
        lons = traj_lon,
        lats = traj_lat,
        depths = traj_depth,
        temperatures = traj_temp,
        degree_days = current_degree_days,
        degree_days_timeseries = traj_dd,
        survival_probability = traj_surv,
        stages = traj_stage,
        alive = current_alive,
        settlement_status = current_settlement,
        settlement_age = current_settlement_age,
        times = collect(times),
        ids = larvae.id
    )
end
