"""
    larval_behavior.jl

Snow crab (*Chionoecetes opilio*) larval behavior, Diel Vertical Migration (DVM),
thermal degree-day ontogenetic molting, tidal current superposition, and
benthic nursery settlement suitability filtering.
"""

using Random

"""
    larval_ascent_velocity(
        z::Real;
        z_surface_target::Real = -10.0,
        max_ascent_speed::Real = 0.010,
        relaxation_depth::Real = 10.0
    ) -> Float64

Compute active post-hatch vertical ascent velocity \$w_{\\text{ascent}}(z)\$ for
newly hatched snow crab larvae (*Chionoecetes opilio*) rising from benthic hatching
grounds on the ocean floor to the epipelagic surface mixed layer.

# Mathematical & Biological Formulation
Ovigerous female snow crabs dwell on the seabed (\$z_{\\text{bed}}\$) and release
larvae into near-bottom waters during the spring phytoplankton bloom (Lovrich et al.,
1995; Sainte-Marie & Sainte-Marie, 1999). Newly hatched Stage I zoeae exhibit strong
negative geotaxis (active upward swimming against gravity) and positive phototaxis
(Sulkin, 1984; Forward, 1988).

The active upward swimming velocity \$w_{\\text{ascent}}(z)\$ directed toward the
productive surface mixed layer (\$z_{\\text{target}}\$) is parameterized via a
hyperbolic tangent transition scale:
```math
w_{\\text{ascent}}(z) = w_{\\text{ascent,max}} \\tanh\\left(
    \\frac{z_{\\text{target}} - z}{L_{\\text{relax}}}
\\right)
```
where \$w_{\\text{ascent,max}}\$ is the maximum directed ascent speed (typically
\$0.008\\text{--}0.015\\text{ m s}^{-1} \\approx 8\\text{--}15\\text{ mm s}^{-1}\$),
\$z_{\\text{target}}\$ is the upper euphotic/mixed layer depth (default \$-10.0\\text{ m}\$),
and \$L_{\\text{relax}}\$ is the vertical relaxation scale (default \$10.0\\text{ m}\$).
Deep in the water column (\$z \\ll z_{\\text{target}}\$), \$\\tanh \\approx 1.0\$,
providing sustained upward swimming. As the larva nears the surface target,
the velocity smoothly decelerates to zero.

# Inputs
- `z::Real`: Current vertical coordinate in meters (negative downward).
- `z_surface_target::Real`: Target epipelagic depth in meters (default -10.0 m).
- `max_ascent_speed::Real`: Maximum vertical ascent speed in \$m s^{-1}\$ (default 0.010 m/s).
- `relaxation_depth::Real`: Vertical deceleration scale in meters (default 10.0 m).

# Outputs
- `Float64`: Vertical ascent velocity in \$m s^{-1}\$ (positive upward).

# References
- Forward, R. B. (1988). Diel vertical migration: zooplankton photobiology and
  behaviour. *Oceanography and Marine Biology: An Annual Review*, 26, 361-393.
- Lovrich, G. A., Sainte-Marie, B., & Smith, B. D. (1995). Depth distribution and
  seasonal movements of *Chionoecetes opilio* in Baie Sainte-Marguerite.
  *Canadian Journal of Fisheries and Aquatic Sciences*, 52(4), 903-913.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). Growth, developmental stages, and
  vertical distribution of snow crab larvae (*Chionoecetes opilio*).
  *Canadian Journal of Fisheries and Aquatic Sciences*, 56(11), 2181-2193.
- Sulkin, S. D. (1984). The behavioral basis of depth regulation in the planktonic
  larvae of brachyuran crabs. *Marine Ecology Progress Series*, 193-205.
"""
function larval_ascent_velocity(
    z::Real;
    z_surface_target::Real = -10.0,
    max_ascent_speed::Real = 0.010,
    relaxation_depth::Real = 10.0
)
    dz = Float64(z_surface_target) - Float64(z)
    return Float64(max_ascent_speed * tanh(dz / relaxation_depth))
end

"""
    diel_vertical_migration_velocity(
        z::Real,
        t::Real;
        stage::Symbol = :zoea1,
        day_depth::Real = -50.0,
        night_depth::Real = -10.0,
        zoea2_day_depth::Real = -55.0,
        zoea2_night_depth::Real = -8.0,
        megalopa_day_depth::Real = -120.0,
        megalopa_night_depth::Real = -60.0,
        zoea2_depth_factor::Real = 1.2,
        megalopa_swim_factor::Real = 1.5,
        max_swim_speed::Real = 0.005,
        relaxation_depth::Real = 10.0,
        turbidity_attenuation::Real = 1.0,
        cil_depth::Union{Nothing, Real} = nothing,
        initial_ascent::Bool = false,
        ascent_speed::Real = 0.010,
        surface_target::Real = -10.0
    )

Compute active vertical swimming velocity \$w_{\\text{swim}}(z, t)\$ for snow crab
larvae undergoing stage-dependent Diel Vertical Migration (DVM) or initial post-hatch
ascent from benthic spawning grounds.

# Mathematical & Ecological Model
Snow crab (*Chionoecetes opilio*) larvae exhibit stage-dependent vertical
positioning (Incze et al., 1987; Lovrich et al., 1995; Sainte-Marie & Sainte-Marie, 1999).
When `initial_ascent = true`, newly released larvae swim actively upward toward the
surface mixed layer at speed `ascent_speed` via `larval_ascent_velocity`.
Once established in the pelagic zone (`initial_ascent = false`), planktonic Zoea I and II
stages occupy shallower depths at night and migrate to deeper, colder layers during the day.
Semiplanktonic Megalopa larvae actively seek benthic nursery grounds prior to settlement.

The instantaneous target depth \$z_{\\text{target}}(t)\$ over a 24-hour diurnal
cycle (\$T = 86400\\text{ s}\$) is:
```math
z_{\\text{target}}(t) = \\frac{z_{\\text{day}} + z_{\\text{night}}}{2}
- \\frac{z_{\\text{night}} - z_{\\text{day}}}{2} \\cos\\left( \\frac{2\\pi t}{86400} \\right)
```
The active swimming velocity \$w_{\\text{swim}}(z, t)\$ is:
```math
w_{\\text{swim}}(z, t) = w_{\\max} \\tanh\\left( \\frac{z_{\\text{target}}(t) - z}{L_{\\text{relax}}} \\right)
```

# Inputs
- `z::Real`: Current vertical coordinate (depth in meters, negative downward).
- `t::Real`: Current simulation time in seconds.
- `stage::Symbol`: Larval development stage (`:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`).
- `day_depth::Real`: Zoea I daytime preferred depth in meters (default -50.0 m).
- `night_depth::Real`: Zoea I nighttime preferred depth in meters (default -10.0 m).
- `zoea2_day_depth::Real`: Zoea II daytime depth in meters (default -55.0 m).
- `zoea2_night_depth::Real`: Zoea II nighttime depth in meters (default -8.0 m).
- `megalopa_day_depth::Real`: Megalopa daytime benthic seeking depth in meters (default -120.0 m).
- `megalopa_night_depth::Real`: Megalopa nighttime depth in meters (default -60.0 m).
- `zoea2_depth_factor::Real`: Zoea II swim speed multiplier (default 1.2).
- `megalopa_swim_factor::Real`: Swimming speed multiplier for Megalopa (default 1.5).
- `max_swim_speed::Real`: Maximum swimming speed in \$m s^{-1}\$ (default 0.005 m/s ≈ 5 mm/s).
- `relaxation_depth::Real`: Vertical transition scale \$L_{\\text{relax}}\$ in meters (default 10.0 m).
- `turbidity_attenuation::Real`: Light cue attenuation factor in turbid coastal plumes (0.0 to 1.0, default 1.0).
- `cil_depth::Union{Nothing, Real}`: Optional upper depth of the Cold Intermediate Layer (thermocline).
- `initial_ascent::Bool`: Whether the larva is in the initial post-hatch ascent phase (default false).
- `ascent_speed::Real`: Upward swimming speed during ascent in \$m s^{-1}\$ (default 0.010 m/s).
- `surface_target::Real`: Target depth for ascent completion in meters (default -10.0 m).

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
    zoea2_day_depth::Real = -55.0,
    zoea2_night_depth::Real = -8.0,
    megalopa_day_depth::Real = -120.0,
    megalopa_night_depth::Real = -60.0,
    zoea2_depth_factor::Real = 1.2,
    megalopa_swim_factor::Real = 1.5,
    max_swim_speed::Real = 0.005,
    relaxation_depth::Real = 10.0,
    turbidity_attenuation::Real = 1.0,
    cil_depth::Union{Nothing, Real} = nothing,
    initial_ascent::Bool = false,
    ascent_speed::Real = 0.010,
    surface_target::Real = -10.0
)
    if stage == :instar1_settled || stage == :dead
        return 0.0
    end

    if initial_ascent
        return larval_ascent_velocity(
            z,
            z_surface_target = surface_target,
            max_ascent_speed = ascent_speed,
            relaxation_depth = relaxation_depth
        )
    end

    z_day, z_night, w_max = if stage == :zoea1
        (Float64(day_depth), Float64(night_depth), Float64(max_swim_speed))
    elseif stage == :zoea2
        # Non-linear ontogenetic niche: Zoea II deepens slightly during daytime
        # for feeding (-55 m), while nighttime ascent is slightly shallower (-8 m)
        # to avoid surface predators (Sainte-Marie & Sainte-Marie 1999).
        (Float64(zoea2_day_depth),
         Float64(zoea2_night_depth),
         Float64(max_swim_speed * zoea2_depth_factor))
    elseif stage == :megalopa
        (Float64(megalopa_day_depth), Float64(megalopa_night_depth),
         Float64(max_swim_speed * megalopa_swim_factor))
    else
        (Float64(day_depth), Float64(night_depth), Float64(max_swim_speed))
    end

    diurnal_phase = 2.0 * π * (mod(t, 86400.0)) / 86400.0
    z_mean = 0.5 * (z_day + z_night)
    # Positive amplitude modulated by turbidity/light attenuation (Incze et al., 1987)
    atten = clamp(Float64(turbidity_attenuation), 0.0, 1.0)
    z_amp = 0.5 * (z_night - z_day) * atten
    z_target = z_mean - z_amp * cos(diurnal_phase)

    # Ontogenetic niche partitioning relative to CIL boundary
    if !isnothing(cil_depth)
        z_cil = Float64(cil_depth)
        if stage == :zoea1
            # Zoea I stay within warm surface mixed layer above CIL
            z_target = max(z_target, z_cil)
        elseif stage == :zoea2
            # Zoea II target the pycnocline / CIL upper transition
            z_target = clamp(z_target, z_cil - 20.0, z_night)
        elseif stage == :megalopa
            # Megalopae seek deep CIL for pre-settlement benthic staging
            z_target = min(z_target, z_cil)
        end
    end

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
- `constituents::Union{Nothing, Vector{Symbol}}`: Optional vector of constituents (e.g. `[:M2, :S2]`).
- `u_amplitudes, v_amplitudes::Union{Nothing, Dict{Symbol, Float64}}`: Multi-constituent amplitudes.
- `phases::Union{Nothing, Dict{Symbol, Float64}}`: Multi-constituent phase offsets.

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
    phase::Real = 0.0,
    constituents::Union{Nothing, Vector{Symbol}} = nothing,
    u_amplitudes::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing,
    v_amplitudes::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing,
    phases::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing
)
    if !isnothing(constituents) && !isempty(constituents)
        u_amps = something(u_amplitudes, Dict(:M2 => Float64(u_amp)))
        v_amps = something(v_amplitudes, Dict(:M2 => Float64(v_amp)))
        phs    = something(phases, Dict(:M2 => Float64(phase)))
        u_tide, v_tide = tidal_velocity_vector(
            t,
            constituents = constituents,
            u_amplitudes = u_amps,
            v_amplitudes = v_amps,
            phases = phs
        )
        return (Float64(u_mean + u_tide), Float64(v_mean + v_tide))
    else
        omega = 2.0 * π / period
        u_tide = u_amp * cos(omega * t + phase)
        v_tide = v_amp * sin(omega * t + phase)
        return (Float64(u_mean + u_tide), Float64(v_mean + v_tide))
    end
end

"""
    bbl_velocity_factor(
        z::Real,
        z_bed::Real;
        h_bbl::Real = 10.0,
        z0::Real = 0.001
    ) -> Float64

Compute the logarithmic bottom boundary layer (BBL) velocity attenuation factor
\$f_{\\text{bbl}} \\in [0, 1]\$ based on the law of the wall:
```math
f_{\\text{bbl}}(z) = \\frac{\\ln(\\max(z_0, z - z_{\\text{bed}}) / z_0)}{\\ln(h_{\\text{bbl}} / z_0)}
```
where \$h_{\\text{bbl}}\$ is the BBL thickness (meters) and \$z_0\$ is bottom roughness length.

# Inputs
- `z::Real`: Particle vertical depth in meters (negative downward).
- `z_bed::Real`: Seabed elevation in meters.
- `h_bbl::Real`: Height of bottom boundary layer in meters (default 10.0 m).
- `z0::Real`: Bottom roughness length in meters (default 0.001 m).

# Outputs
- `Float64`: Attenuation factor in \$[0, 1]\$.
"""
function bbl_velocity_factor(
    z::Real,
    z_bed::Real;
    h_bbl::Real = 10.0,
    z0::Real = 0.001
)
    dz = Float64(z) - Float64(z_bed)
    if dz <= 0.0
        return 0.0
    elseif dz >= h_bbl
        return 1.0
    else
        return clamp(log(max(z0, dz) / z0) / log(h_bbl / z0), 0.0, 1.0)
    end
end

"""
    larval_passive_sinking_velocity(
        stage::Symbol;
        zoea1_sink::Real = -0.0005,
        zoea2_sink::Real = -0.0010,
        megalopa_sink::Real = -0.0025
    ) -> Float64

Compute the passive gravitational settling velocity \$w_{\\text{sink}}\$ in \$m s^{-1}\$
for snow crab larvae due to excess body density (\$\\Delta \\rho \\approx 15\\text{ kg/m}^3\$).

# Inputs
- `stage::Symbol`: Developmental stage (`:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`).
- `zoea1_sink::Real`: Zoea I passive sinking speed in \$m s^{-1}\$ (default -0.5 mm/s).
- `zoea2_sink::Real`: Zoea II passive sinking speed in \$m s^{-1}\$ (default -1.0 mm/s).
- `megalopa_sink::Real`: Megalopa passive sinking speed in \$m s^{-1}\$ (default -2.5 mm/s).

# Outputs
- `Float64`: Sinking velocity in \$m s^{-1}\$ (negative downward).

# References
- Sainte-Marie, G., & Sainte-Marie, B. (1999). *CJFAS*, 56(11), 2181-2193.
"""
function larval_passive_sinking_velocity(
    stage::Symbol;
    zoea1_sink::Real = -0.0005,
    zoea2_sink::Real = -0.0010,
    megalopa_sink::Real = -0.0025
)
    if stage == :zoea1
        return Float64(zoea1_sink)
    elseif stage == :zoea2
        return Float64(zoea2_sink)
    elseif stage == :megalopa
        return Float64(megalopa_sink)
    else
        return 0.0
    end
end

"""
    update_larval_stage(
        current_stage::Symbol,
        degree_days::Real;
        dd_zoea1_to_zoea2::Real = 65.0,
        dd_zoea2_to_megalopa::Real = 130.0,
        dd_megalopa_to_settle::Real = 200.0
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
        optimal_min_depth::Real = -180.0,
        optimal_max_depth::Real = -80.0,
        min_temp::Real = -1.0,
        optimal_min_temp::Real = 0.5,
        optimal_max_temp::Real = 3.5,
        max_bottom_temp::Real = 6.0,
        stochastic::Bool = false,
        rng::AbstractRNG = Random.default_rng()
    )

Evaluate whether a settling Megalopa/Instar I larva has arrived on suitable benthic
nursery grounds on the continental shelf using a continuous Habitat Suitability Index (HSI).

# Mathematical & Ecological Formulation
1. **Depth Suitability Index** \$S_z(z) \\in [0, 1]\$:
   - Optimal nursery depths: \$z \\in [-180, -80]\\text{ m}\$.
   - Tapers to 0 at nearshore surf zone (\$z > -50\\text{ m}\$) and offshore abyss (\$z < -250\\text{ m}\$).
2. **Thermal Suitability Index** \$S_T(T) \\in [0, 1]\$:
   - Optimal Cold Intermediate Layer (CIL) nursery: \$T \\in [0.5, 3.5]^\\circ\\text{C}\$.
   - Tapers to 0 at thermal stress threshold (\$T > 6.0^\\circ\\text{C}\$) and superchilling (\$T < -1.0^\\circ\\text{C}\$).
3. **Combined Habitat Suitability Index**:
   ```math
   \\text{HSI}(z, T) = S_z(z) \\times S_T(T) \\in [0, 1]
   ```
   In stochastic mode, settlement succeeds if \$r \\sim U(0, 1) \\le \\text{HSI}\$.
   In deterministic mode, settlement succeeds if \$\\text{HSI} > 0.0\$.

# Inputs
- `bed_elevation::Real`: Seafloor elevation in meters (negative).
- `bottom_temperature::Real`: Bottom water temperature in °C.
- `min_depth::Real`: Deepest boundary for nursery habitat (default -250.0 m).
- `max_depth::Real`: Shallowest boundary for nursery habitat (default -50.0 m).
- `optimal_min_depth::Real`: Deep limit of optimal nursery depth (default -180.0 m).
- `optimal_max_depth::Real`: Shallow limit of optimal nursery depth (default -80.0 m).
- `min_temp::Real`: Sub-zero superchilling limit in °C (default -1.0 °C).
- `optimal_min_temp::Real`: Minimum temperature of optimal CIL (default 0.5 °C).
- `optimal_max_temp::Real`: Maximum temperature of optimal CIL (default 3.5 °C).
- `max_bottom_temp::Real`: Maximum allowable bottom temperature in °C (default 6.0 °C).
- `stochastic::Bool`: If true, samples settlement probabilistically from HSI (default false).
- `rng::AbstractRNG`: Random number generator for stochastic settlement.

# Outputs
- `NamedTuple`: `(suitable::Bool, probability::Float64, hsi::Float64, reason::String)`

# References
- Dionne, M., Sainte-Marie, B., Bourget, E., & Gilbert, D. (2003). *MEPS*, 259, 117-128.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). *MEPS*, 182, 157-174.
"""
function evaluate_settlement_suitability(
    bed_elevation::Real,
    bottom_temperature::Real;
    min_depth::Real = -250.0,
    max_depth::Real = -50.0,
    optimal_min_depth::Real = -180.0,
    optimal_max_depth::Real = -80.0,
    min_temp::Real = -1.0,
    optimal_min_temp::Real = 0.5,
    optimal_max_temp::Real = 3.5,
    max_bottom_temp::Real = 6.0,
    stochastic::Bool = false,
    rng::AbstractRNG = Random.default_rng()
)
    z = Float64(bed_elevation)
    T = Float64(bottom_temperature)

    # Depth suitability S_z
    s_z = if z > max_depth || z < min_depth
        0.0
    elseif z >= optimal_min_depth && z <= optimal_max_depth
        1.0
    elseif z > optimal_max_depth
        (z - max_depth) / (optimal_max_depth - max_depth)
    else
        (z - min_depth) / (optimal_min_depth - min_depth)
    end

    # Thermal suitability S_T
    s_t = if T < min_temp || T > max_bottom_temp
        0.0
    elseif T >= optimal_min_temp && T <= optimal_max_temp
        1.0
    elseif T > optimal_max_temp
        (max_bottom_temp - T) / (max_bottom_temp - optimal_max_temp)
    else
        (T - min_temp) / (optimal_min_temp - min_temp)
    end

    hsi = clamp(s_z * s_t, 0.0, 1.0)

    suitable = if stochastic
        rand(rng) <= hsi
    else
        hsi > 0.0
    end

    reason = if hsi == 0.0
        if z > max_depth
            "Too shallow / nearshore surf zone (z > $(max_depth)m)"
        elseif z < min_depth
            "Too deep / offshore abyssal slope (z < $(min_depth)m)"
        elseif T > max_bottom_temp
            "Thermal stress / bottom temp > $(max_bottom_temp)°C"
        else
            "Sub-zero thermal stress (T < $(min_temp)°C)"
        end
    elseif suitable
        "Suitable cold-water shelf nursery ground (HSI = $(round(hsi, digits=3)))"
    else
        "Stochastic settlement search failed on marginal ground (HSI = $(round(hsi, digits=3)))"
    end

    return (suitable = suitable, probability = Float64(hsi), hsi = Float64(hsi), reason = reason)
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
advection, turbulent diffusion, active larval vertical swimming, and absorbing
shoreline boundary dynamics with alongshore tangential slip.

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
- `is_ascending::Bool`: Whether larva is currently in directed vertical ascent phase (default false).
- `ascent_speed::Real`: Upward swimming speed during ascent in \$m s^{-1}\$ (default 0.010 m/s).
- `surface_target::Real`: Target epipelagic depth for ascent in meters (default -10.0 m).
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `Tuple{Float64, Float64, Float64}`: Updated coordinates `(x_new, y_new, z_new)`.

# References
- North, E. W., et al. (2008). Simulating dispersal of marine organisms: validation and
  sensitivity. *Marine Ecology Progress Series*, 372, 283-294.
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
    coastline = nothing,
    enable_bbl::Bool = true,
    h_bbl::Real = 10.0,
    z0::Real = 0.001,
    enable_sinking::Bool = true,
    zoea1_sink::Real = -0.0005,
    zoea2_sink::Real = -0.0010,
    megalopa_sink::Real = -0.0025,
    κ_v_gradient::Real = 0.0,
    is_ascending::Bool = false,
    ascent_speed::Real = 0.010,
    surface_target::Real = -10.0,
    rng::AbstractRNG = Random.default_rng()
)
    # If settled on bottom or dead, particle remains at seabed
    if stage == :instar1_settled
        return (Float64(x), Float64(y), Float64(z_bottom))
    end

    # 1. Logarithmic Bottom Boundary Layer (BBL) velocity attenuation (law of the wall)
    f_bbl = enable_bbl ? bbl_velocity_factor(z, z_bottom, h_bbl = h_bbl, z0 = z0) : 1.0
    u_eff = u * f_bbl
    v_eff = v * f_bbl

    # 2. Vertical swimming and passive gravitational settling
    w_swim = diel_vertical_migration_velocity(
        z,
        t;
        stage = stage,
        initial_ascent = is_ascending,
        ascent_speed = ascent_speed,
        surface_target = surface_target
    )
    w_sink = enable_sinking ? larval_passive_sinking_velocity(
        stage,
        zoea1_sink = zoea1_sink,
        zoea2_sink = zoea2_sink,
        megalopa_sink = megalopa_sink
    ) : 0.0

    ξ_x = randn(rng)
    ξ_y = randn(rng)
    ξ_z = randn(rng)

    dx_meters = u_eff * dt + sqrt(max(0.0, 2.0 * κ_h * dt)) * ξ_x
    dy_meters = v_eff * dt + sqrt(max(0.0, 2.0 * κ_h * dt)) * ξ_y

    # 3. Visser (1997) diffusive pseudo-drift correction: dκ_v/dz term prevents artificial
    # accumulation of particles in low-diffusivity pycnoclines
    dz_meters = (w + w_swim + w_sink + κ_v_gradient) * dt + sqrt(max(0.0, 2.0 * κ_v * dt)) * ξ_z

    if is_lat_lon
        r_earth = 6.371e6
        # Protect near-pole singularity: clamp latitude and bound cos(lat) at 85° Web Mercator cutoff
        lat_rad = deg2rad(clamp(Float64(y), -85.0, 85.0))
        cos_lat = max(cos(lat_rad), cosd(85.0))
        deg_per_meter_lat = 180.0 / (π * r_earth)
        deg_per_meter_lon = 180.0 / (π * r_earth * cos_lat)

        x_raw = x + dx_meters * deg_per_meter_lon
        y_raw = y + dy_meters * deg_per_meter_lat

        # Periodic antimeridian wrapping to [-180, 180]
        x_new = mod(x_raw + 180.0, 360.0) - 180.0
        y_new = clamp(y_raw, -89.9, 89.9)
    else
        x_new = x + dx_meters
        y_new = y + dy_meters
    end
    z_new = z + dz_meters

    # Absorbing vertical boundaries: larvae cannot exit surface or penetrate seabed
    z_new = clamp(z_new, z_bottom, z_surface)

    # Land boundary condition: absorbing shoreline with alongshore tangential slip.
    # Larvae that step onto land test orthogonal alongshore displacement components.
    # If both components are blocked by land, the particle is absorbed at its previous
    # marine coordinate (x, y). Reversing displacement vectors offshore is unphysical.
    if is_lat_lon && is_point_on_land(x_new, y_new, coastline = coastline)
        # Attempt alongshore tangential slip
        x_zonal = x + dx_meters * deg_per_meter_lon
        y_merid = y + dy_meters * deg_per_meter_lat
        if !is_point_on_land(x_zonal, y, coastline = coastline)
            x_new = x_zonal
            y_new = Float64(y)
        elseif !is_point_on_land(x, y_merid, coastline = coastline)
            x_new = Float64(x)
            y_new = y_merid
        else
            # Both alongshore candidates on land: absorb at previous valid marine location
            x_new = Float64(x)
            y_new = Float64(y)
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
- `depth_range::Union{Nothing, Tuple{Real, Real}}`: Depth release bounds in meters (used if
  `release_depth_mode = :range`, default nothing).
- `release_depth_mode::Union{Nothing, Symbol}`: Vertical release mode: `:bottom` (default,
  near seabed), `:range` (within `depth_range`), or `:surface` (near surface mixed layer).
- `bottom_offset::Tuple{Real, Real}`: Vertical elevation range above seabed in meters when
  `release_depth_mode = :bottom` (default (0.5, 3.0) m).
- `ascent_target_depth::Real`: Target epipelagic depth in meters (default -10.0 m).
- `min_seabed_depth::Real`: Minimum water depth in meters (default 100.0 m, min 0.0 m).
- `buffer_km::Real`: Spatial buffer distance in km beyond stratum/domain bounds (default 0.0 km).
- `bathymetry`: Bathymetry source (`NamedTuple`, NetCDF path, elevation `Function`, or `nothing`).
- `bathymetry_fn`: Alias for `bathymetry`.
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
    depth_range::Union{Nothing, Tuple{Real, Real}} = nothing,
    release_depth_mode::Union{Nothing, Symbol} = nothing,
    bottom_offset::Tuple{Real, Real} = (0.5, 3.0),
    ascent_target_depth::Real = -10.0,
    min_seabed_depth::Real = 100.0,
    buffer_km::Real = 0.0,
    bathymetry::Union{Nothing, Function, NamedTuple, AbstractString} = nothing,
    bathymetry_fn::Union{Nothing, Function, NamedTuple, AbstractString} = nothing,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    stratum::Union{Nothing, NamedTuple, Symbol, AbstractString} = nothing,
    stage::Symbol = :zoea1,
    max_attempts::Int = 100000,
    rng::AbstractRNG = Random.default_rng()
)
    if n_particles <= 0
        error("Number of particles must be positive: $(n_particles)")
    end

    eff_mode = if !isnothing(release_depth_mode)
        release_depth_mode
    elseif !isnothing(depth_range)
        :range
    else
        :bottom
    end

    sample_initial_depth = function(bed_z::Float64)
        if eff_mode == :bottom
            Δz = bottom_offset[1] + rand(rng, Float64) * (bottom_offset[2] - bottom_offset[1])
            return clamp(bed_z + Δz, bed_z, -1.0)
        elseif eff_mode == :surface
            z_surf = ascent_target_depth + rand(rng, Float64) * 5.0
            return clamp(z_surf, bed_z, -1.0)
        else # :range
            drange = something(depth_range, (-60.0, -20.0))
            raw_z = drange[1] + rand(rng, Float64) * (drange[2] - drange[1])
            return clamp(raw_z, bed_z + 2.0, -1.0)
        end
    end

    # Auto-resolve bathymetry if not explicitly provided
    resolved_bathy = if !isnothing(bathymetry)
        bathymetry
    elseif !isnothing(bathymetry_fn)
        bathymetry_fn
    else
        def_path = "inputs/bathymetry_active.nc"
        isfile(def_path) ? def_path : nothing
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
            depths[p] = sample_initial_depth(sampled_zbed[p])
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

            cand_z = sample_initial_depth(z_bed)

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
- `settlement_stochastic::Bool`: Whether to evaluate settlement probabilistically via HSI (default true).
- `competence_window_days::Real`: Maximum pelagic search window for competent megalopae in days (default 15.0 days).
- `tidal_u_amp::Real`: Zonal tidal velocity amplitude in m/s (default 0.25 m/s).
- `tidal_v_amp::Real`: Meridional tidal velocity amplitude in m/s (default 0.12 m/s).
- `tidal_period::Real`: Tidal constituent period in seconds (default 44712.0 s).
- `tidal_phase::Real`: Tidal initial phase offset in radians (default 0.0).
- `min_survival_prob::Real`: Survival probability below which particle is marked dead (default 0.01).
- `enable_initial_ascent::Bool`: Whether newly released larvae perform directed vertical ascent
  from benthic release depths to the surface mixed layer (default true).
- `ascent_speed::Real`: Upward swimming speed during ascent in \$m s^{-1}\$ (default 0.010 m/s).
- `ascent_target_depth::Real`: Target epipelagic depth for ascent in meters (default -10.0 m).
- `max_ascent_duration::Real`: Maximum allowed duration for initial ascent in seconds (default 86400.0 s).
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `NamedTuple`: Trajectory record containing `(lons, lats, depths, stages, degree_days, alive, settlement_status, times, ids, ascent_duration)`.
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
    settlement_stochastic::Bool = true,
    competence_window_days::Real = 15.0,
    tidal_u_amp::Real = 0.25,
    tidal_v_amp::Real = 0.12,
    tidal_period::Real = 44712.0,
    tidal_phase::Real = 0.0,
    tidal_constituents::Union{Nothing, Vector{Symbol}} = nothing,
    tidal_u_amplitudes::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing,
    tidal_v_amplitudes::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing,
    tidal_phases::Union{Nothing, AbstractDict{Symbol, <:Real}} = nothing,
    enable_bbl::Bool = true,
    h_bbl::Real = 10.0,
    z0::Real = 0.001,
    enable_sinking::Bool = true,
    zoea1_sink::Real = -0.0005,
    zoea2_sink::Real = -0.0010,
    megalopa_sink::Real = -0.0025,
    κ_v_profile::Union{Nothing, Function} = nothing,
    min_survival_prob::Real = 0.01,
    coastline = nothing,
    enable_initial_ascent::Bool = true,
    ascent_speed::Real = 0.010,
    ascent_target_depth::Real = -10.0,
    max_ascent_duration::Real = 86400.0,
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
    competence_onset_step  = fill(-1, n_particles)
    t_bed_filtered = zeros(Float64, n_particles)
    stage_survival = fill(1.0, n_particles)
    ascent_complete = Vector{Bool}(undef, n_particles)
    ascent_duration = fill(0.0, n_particles)

    # Set initial states at t = 0
    traj_lon[:, 1] = copy(larvae.lon)
    traj_lat[:, 1] = copy(larvae.lat)
    traj_depth[:, 1] = copy(larvae.depth)
    traj_stage[:, 1] = copy(larvae.stage)
    traj_dd[:, 1] = copy(current_degree_days)
    traj_surv[:, 1] .= 1.0

    for p in 1:n_particles
        # Particles already at or above surface target have completed initial ascent
        ascent_complete[p] = !enable_initial_ascent || (traj_depth[p, 1] >= ascent_target_depth)
    end

    for p in 1:n_particles
        init_T = isnothing(temperature_fn) ? Float64(default_temperature) :
                 temperature_fn(traj_lon[p, 1], traj_lat[p, 1], traj_depth[p, 1], 0.0)
        traj_temp[p, 1] = init_T
        init_bed = isnothing(bathymetry_fn) ? Float64(default_bottom_depth) :
                   bathymetry_fn(traj_lon[p, 1], traj_lat[p, 1])
        t_bed_filtered[p] = isnothing(temperature_fn) ? Float64(default_temperature) :
                            temperature_fn(traj_lon[p, 1], traj_lat[p, 1], init_bed, 0.0)
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
                traj_stage[p, s + 1] = traj_stage[p, s]
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
                    phase = tidal_phase,
                    constituents = tidal_constituents,
                    u_amplitudes = tidal_u_amplitudes,
                    v_amplitudes = tidal_v_amplitudes,
                    phases = tidal_phases
                )
            end

            # Query temperature and seabed depth
            cur_T = isnothing(temperature_fn) ? Float64(default_temperature) :
                    temperature_fn(cur_lon, cur_lat, cur_depth, t_current)
            z_bed = isnothing(bathymetry_fn) ? Float64(default_bottom_depth) :
                    bathymetry_fn(cur_lon, cur_lat)
            cur_bed_T = isnothing(temperature_fn) ? Float64(default_temperature) :
                        temperature_fn(cur_lon, cur_lat, z_bed, t_current)

            # Filter high-frequency semi-diurnal M2 tidal temperature oscillations (tau ≈ 12.42 h = 44712 s)
            alpha_tidal = clamp(dt / 44712.0, 0.005, 1.0)
            t_bed_filtered[p] = (1.0 - alpha_tidal) * t_bed_filtered[p] + alpha_tidal * cur_bed_T

            # Stage-specific mortality scaling (Zoea I fragile, Megalopa robust)
            stage_mort_factor = cur_stage == :zoea1 ? 1.1 : (cur_stage == :zoea2 ? 1.0 : 0.8)
            mort_rate = larval_thermal_mortality_rate(
                cur_T,
                base_mortality       = mortality_base * stage_mort_factor,
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

                if new_stage != cur_stage
                    # Record survival fraction across stage molt
                    stage_survival[p] = traj_surv[p, s + 1]
                end

                # Settlement evaluation for competent larvae
                is_competent = (new_stage == :instar1_settled) ||
                               (cur_stage == :megalopa && current_degree_days[p] >= dd_megalopa_to_settle)

                if is_competent
                    if competence_onset_step[p] < 0
                        competence_onset_step[p] = s
                    end
                    days_competent = (s - competence_onset_step[p]) * dt_days

                    # Evaluate settlement suitability against tidally filtered benthic temperature
                    suit = evaluate_settlement_suitability(
                        z_bed,
                        t_bed_filtered[p],
                        min_depth = settlement_min_depth,
                        max_depth = settlement_max_depth,
                        max_bottom_temp = settlement_max_temp,
                        stochastic = settlement_stochastic,
                        rng = rng
                    )

                    if suit.suitable
                        current_settlement[p] = :settled_successful
                        current_settlement_age[p] = t_current
                        cur_stage = :instar1_settled
                    elseif days_competent > competence_window_days
                        # Competence search window expired without reaching nursery habitat
                        current_settlement[p] = :settled_unsuitable
                        current_settlement_age[p] = t_current
                        current_alive[p] = false
                        cur_stage = :megalopa
                    else
                        # Delay settlement: competent larva continues pelagic search phase
                        current_settlement[p] = :pelagic
                        cur_stage = :megalopa
                    end
                else
                    cur_stage = new_stage
                end
            end
            traj_dd[p, s + 1] = current_degree_days[p]

            # Depth-dependent vertical diffusivity and Visser (1997) gradient drift
            local_kv = if !isnothing(κ_v_profile)
                Float64(κ_v_profile(cur_depth))
            else
                Float64(κ_v)
            end
            kv_grad = if !isnothing(κ_v_profile)
                δz = 0.5
                (Float64(κ_v_profile(cur_depth + δz)) - Float64(κ_v_profile(cur_depth - δz))) / (2.0 * δz)
            else
                0.0
            end

            is_ascending = enable_initial_ascent && !ascent_complete[p]

            # Step Lagrangian transport with BBL shear, sinking, and diffusive drift
            new_lon, new_lat, new_depth = larval_transport_step(
                cur_lon, cur_lat, cur_depth,
                u, v, w, κ_h, local_kv, dt,
                t = t_current, stage = cur_stage,
                z_bottom = z_bed,
                is_lat_lon = is_lat_lon,
                coastline = coastline,
                enable_bbl = enable_bbl,
                h_bbl = h_bbl,
                z0 = z0,
                enable_sinking = enable_sinking,
                zoea1_sink = zoea1_sink,
                zoea2_sink = zoea2_sink,
                megalopa_sink = megalopa_sink,
                κ_v_gradient = kv_grad,
                is_ascending = is_ascending,
                ascent_speed = ascent_speed,
                surface_target = ascent_target_depth,
                rng = rng
            )

            # Check whether ascending larva reached the surface mixed layer
            if is_ascending
                if new_depth >= (ascent_target_depth - 0.5) || (t_current + dt) >= max_ascent_duration
                    ascent_complete[p] = true
                    ascent_duration[p] = t_current + dt
                end
            end

            # Coastal landmass boundary condition: absorbing boundary via polygon ray-casting and bathymetry.
            on_land = (is_lat_lon && is_point_on_land(new_lon, new_lat, coastline = coastline)) ||
                      (!isnothing(bathymetry_fn) && bathymetry_fn(new_lon, new_lat) >= 0.0)
            if on_land
                new_lon = cur_lon
                new_lat = cur_lat
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
        stage_survival = stage_survival,
        stages = traj_stage,
        alive = current_alive,
        settlement_status = current_settlement,
        settlement_age = current_settlement_age,
        ascent_duration = ascent_duration,
        times = collect(times),
        ids = larvae.id
    )
end

"""
    safe_mean(collection; default::Real = 0.0) -> Float64

Compute the arithmetic mean of a numerical collection safely. If the collection
is empty or contains only NaN values, returns `default` rather than throwing an
`ArgumentError` or runtime exception.

# Mathematical Formulation
```math
\\bar{x} = \\begin{cases}
\\frac{1}{N} \\sum_{i=1}^N x_i, & N > 0 \\\\
x_{\\text{default}}, & N = 0
\\end{cases}
```

# Inputs
- `collection`: Iterable or array of numeric values.
- `default::Real`: Fallback value returned when collection is empty (default 0.0).

# Outputs
- `Float64`: Mean value or default fallback.
"""
function safe_mean(collection; default::Real = 0.0)
    if isempty(collection)
        return Float64(default)
    end
    valid_vals = [Float64(x) for x in collection if !isnan(x)]
    if isempty(valid_vals)
        return Float64(default)
    end
    return mean(valid_vals)
end

"""
    canonicalize_trajectories(trajectories::NamedTuple) -> NamedTuple

Normalize and validate a Lagrangian trajectory dataset to guarantee uniform array
dimensions, consistent canonical Symbol types for stages and statuses, and complete
fields required for DuckDB storage, spatial analysis, and interactive visualization.

# Canonical Types & Array Specifications
- `lons::Matrix{Float64}`: 2D array of longitudes `(n_particles, n_times)`.
- `lats::Matrix{Float64}`: 2D array of latitudes `(n_particles, n_times)`.
- `depths::Matrix{Float64}`: 2D array of depths in meters `(n_particles, n_times)`.
- `times::Vector{Float64}`: 1D array of timestamps in seconds `(n_times,)`.
- `stages::Matrix{Symbol}`: 2D array of canonical stage symbols `(n_particles, n_times)`
  (e.g., `:zoea1`, `:zoea2`, `:megalopa`, `:instar1_settled`, `:dead`).
- `alive::Vector{Bool}`: 1D boolean array indicating final alive state `(n_particles,)`.
- `settlement_status::Vector{Symbol}`: 1D canonical symbols `(n_particles,)`
  (e.g., `:pelagic`, `:settled_successful`, `:settled_unsuitable`).
- `settlement_age::Vector{Float64}`: 1D array of settlement ages in seconds `(n_particles,)`.
- `degree_days::Vector{Float64}`: 1D cumulative degree-days per particle `(n_particles,)`.
- `degree_days_timeseries::Matrix{Float64}`: 2D degree-days timeseries `(n_particles, n_times)`.
- `temperatures::Matrix{Float64}`: 2D water temperatures in °C `(n_particles, n_times)`.
- `survival_probability::Matrix{Float64}`: 2D survival rates in `[0, 1]` `(n_particles, n_times)`.
- `ids::Vector{Int}`: 1D unique particle identifiers `(n_particles,)`.

# Inputs
- `trajectories::NamedTuple`: Raw or partially constructed trajectory record.

# Outputs
- `NamedTuple`: Fully normalized trajectory record with guaranteed schema.
"""
function canonicalize_trajectories(trajectories::NamedTuple)
    if !hasproperty(trajectories, :lons) || !hasproperty(trajectories, :lats)
        error("Trajectory structure must contain :lons and :lats coordinates.")
    end

    raw_lons = trajectories.lons
    raw_lats = trajectories.lats

    n_p, n_t = if ndims(raw_lons) == 2
        size(raw_lons)
    elseif ndims(raw_lons) == 1
        (length(raw_lons), 1)
    else
        error("Unsupported :lons dimensionality: $(ndims(raw_lons)). Expected 1D or 2D.")
    end

    lons = ndims(raw_lons) == 2 ? Float64.(raw_lons) : reshape(Float64.(raw_lons), n_p, 1)
    lats = ndims(raw_lats) == 2 ? Float64.(raw_lats) : reshape(Float64.(raw_lats), n_p, 1)

    depths = if hasproperty(trajectories, :depths)
        raw_z = trajectories.depths
        ndims(raw_z) == 2 ? Float64.(raw_z) : reshape(Float64.(raw_z), n_p, n_t)
    else
        fill(0.0, n_p, n_t)
    end

    times = if hasproperty(trajectories, :times)
        Float64.(collect(trajectories.times))
    else
        collect(range(0.0, step = 300.0, length = n_t))
    end

    ids = if hasproperty(trajectories, :ids)
        Int.(collect(trajectories.ids))
    elseif hasproperty(trajectories, :id)
        Int.(collect(trajectories.id))
    else
        collect(1:n_p)
    end

    alive = if hasproperty(trajectories, :alive)
        raw_a = trajectories.alive
        if length(raw_a) == n_p
            Bool.(collect(raw_a))
        else
            fill(true, n_p)
        end
    else
        fill(true, n_p)
    end

    # Normalize settlement_status to canonical Vector{Symbol}
    settle_syms = if hasproperty(trajectories, :settlement_status)
        raw_s = trajectories.settlement_status
        [Symbol(lowercase(string(s))) for s in raw_s]
    else
        fill(:pelagic, n_p)
    end

    settlement_age = if hasproperty(trajectories, :settlement_age)
        Float64.(collect(trajectories.settlement_age))
    else
        fill(times[end], n_p)
    end

    # Normalize stages to Matrix{Symbol}
    stages_mat = if hasproperty(trajectories, :stages)
        raw_st = trajectories.stages
        if ndims(raw_st) == 2 && size(raw_st) == (n_p, n_t)
            [Symbol(lowercase(string(raw_st[p, t]))) for p in 1:n_p, t in 1:n_t]
        elseif length(raw_st) == n_p
            repeat([Symbol(lowercase(string(s))) for s in raw_st], 1, n_t)
        else
            fill(:zoea1, n_p, n_t)
        end
    elseif hasproperty(trajectories, :stage)
        raw_st = trajectories.stage
        if length(raw_st) == n_p
            repeat([Symbol(lowercase(string(s))) for s in raw_st], 1, n_t)
        else
            fill(:zoea1, n_p, n_t)
        end
    else
        fill(:zoea1, n_p, n_t)
    end

    # Normalize degree_days to 1D vector and 2D timeseries
    dd_ts = if hasproperty(trajectories, :degree_days_timeseries) &&
               ndims(trajectories.degree_days_timeseries) == 2 &&
               size(trajectories.degree_days_timeseries) == (n_p, n_t)
        Float64.(trajectories.degree_days_timeseries)
    elseif hasproperty(trajectories, :degree_days)
        raw_dd = trajectories.degree_days
        if ndims(raw_dd) == 2 && size(raw_dd) == (n_p, n_t)
            Float64.(raw_dd)
        elseif length(raw_dd) == n_p
            repeat(reshape(Float64.(raw_dd), n_p, 1), 1, n_t)
        else
            fill(0.0, n_p, n_t)
        end
    else
        fill(0.0, n_p, n_t)
    end

    dd_vec = if hasproperty(trajectories, :degree_days) && length(trajectories.degree_days) == n_p
        Float64.(collect(trajectories.degree_days))
    else
        Float64.(dd_ts[:, end])
    end

    # Normalize temperatures to 2D matrix
    temps = if hasproperty(trajectories, :temperatures) &&
               ndims(trajectories.temperatures) == 2 &&
               size(trajectories.temperatures) == (n_p, n_t)
        Float64.(trajectories.temperatures)
    elseif hasproperty(trajectories, :temperature) && length(trajectories.temperature) == n_p
        repeat(reshape(Float64.(trajectories.temperature), n_p, 1), 1, n_t)
    else
        fill(4.0, n_p, n_t)
    end

    # Normalize survival_probability to 2D matrix
    survs = if hasproperty(trajectories, :survival_probability) &&
               ndims(trajectories.survival_probability) == 2 &&
               size(trajectories.survival_probability) == (n_p, n_t)
        Float64.(trajectories.survival_probability)
    elseif hasproperty(trajectories, :survival_probability) &&
           length(trajectories.survival_probability) == n_p
        repeat(reshape(Float64.(trajectories.survival_probability), n_p, 1), 1, n_t)
    else
        ones(Float64, n_p, n_t)
    end

    stage_survs = if hasproperty(trajectories, :stage_survival) &&
                     length(trajectories.stage_survival) == n_p
        Float64.(collect(trajectories.stage_survival))
    else
        Float64.(survs[:, end])
    end

    ascent_dur = if hasproperty(trajectories, :ascent_duration) &&
                    length(trajectories.ascent_duration) == n_p
        Float64.(collect(trajectories.ascent_duration))
    else
        zeros(Float64, n_p)
    end

    return (
        lons = lons,
        lats = lats,
        depths = depths,
        times = times,
        stages = stages_mat,
        alive = alive,
        settlement_status = settle_syms,
        settlement_age = settlement_age,
        ascent_duration = ascent_dur,
        degree_days = dd_vec,
        degree_days_timeseries = dd_ts,
        temperatures = temps,
        survival_probability = survs,
        stage_survival = stage_survs,
        ids = ids
    )
end

"""
    canonicalize_status(trajectories::NamedTuple) -> NamedTuple
    canonicalize_status!(trajectories::NamedTuple) -> NamedTuple

Alias for `canonicalize_trajectories` providing in-place or returning normalization of
status and stage symbols across Lagrangian trajectory records.
"""
canonicalize_status(trajectories::NamedTuple) = canonicalize_trajectories(trajectories)
canonicalize_status!(trajectories::NamedTuple) = canonicalize_trajectories(trajectories)

