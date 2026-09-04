"""
    tides.jl

Astronomical tidal constituents, tidal momentum body forcing, harmonic velocity
reconstruction, and Simpson-Hunter tidal mixing front diagnostics for regional
shelf modeling.
"""

using Oceananigans

"""
    get_tidal_frequency(constituent::Symbol)

Return the astronomical angular frequency \$\\omega\$ (\$\\text{rad s}^{-1}\$) for a
specified tidal constituent.

# Supported Constituents
- `:M2`: Principal lunar semidiurnal (\$T = 12.4206\\text{ h}, \\omega = 1.405189 \\times 10^{-4}\\text{ rad s}^{-1}\$)
- `:S2`: Principal solar semidiurnal (\$T = 12.0000\\text{ h}, \\omega = 1.454441 \\times 10^{-4}\\text{ rad s}^{-1}\$)
- `:N2`: Larger lunar elliptic semidiurnal (\$T = 12.6583\\text{ h}, \\omega = 1.378797 \\times 10^{-4}\\text{ rad s}^{-1}\$)
- `:K1`: Lunar diurnal (\$T = 23.9345\\text{ h}, \\omega = 7.292116 \\times 10^{-5}\\text{ rad s}^{-1}\$)
- `:O1`: Principal lunar diurnal (\$T = 25.8193\\text{ h}, \\omega = 6.759774 \\times 10^{-5}\\text{ rad s}^{-1}\$)

# References
- Pugh, D., & Woodworth, P. (2014). *Sea-Level Science: Understanding Tides,
  Surges, Tsunamis and Mean Sea-Level Changes*. Cambridge University Press.
  DOI: 10.1017/CBO9781139235778
"""
function get_tidal_frequency(constituent::Symbol)
    if constituent == :M2
        return 1.405189e-4
    elseif constituent == :S2
        return 1.454441e-4
    elseif constituent == :N2
        return 1.378797e-4
    elseif constituent == :K1
        return 7.292116e-5
    elseif constituent == :O1
        return 6.759774e-5
    else
        error(
            "Unknown tidal constituent '$(constituent)'. " *
            "Supported: :M2, :S2, :N2, :K1, :O1"
        )
    end
end

"""
    build_tidal_body_forcing(;
        constituents::Vector{Symbol} = [:M2],
        u_amplitudes::Dict{Symbol, Float64} = Dict(:M2 => 0.25),
        v_amplitudes::Dict{Symbol, Float64} = Dict(:M2 => 0.12),
        phases::Dict{Symbol, Float64} = Dict(:M2 => 0.0)
    )

Construct Oceananigans momentum forcing functions \$(F_u, F_v)\$ to drive barotropic
tidal oscillations across the model domain.

# Mathematical Formulation
The horizontal momentum body forcing terms \$F_u(t), F_v(t)\$ are:
```math
F_u(x, y, z, t) = \\sum_k U_k \\omega_k \\cos(\\omega_k t + \\phi_k)
```
```math
F_v(x, y, z, t) = \\sum_k V_k \\omega_k \\sin(\\omega_k t + \\phi_k)
```
where \$\\omega_k\$ is the astronomical tidal frequency (rad s⁻¹) and \$(U_k, V_k)\$ are
velocity amplitudes (m s⁻¹).

The factor \$\\omega_k\$ is intentional: \$F_u = U_k \\omega_k \\cos(\\omega_k t)\$ is the
time derivative \$\\partial_t[U_k \\sin(\\omega_k t)]\$, making it dimensionally an
acceleration (m s⁻²), the correct unit for an Oceananigans body force. The
amplitude \$U_k\$ therefore represents the target tidal velocity amplitude in m s⁻¹;
calibrate it against TPXO/FES tidal prediction or observed current meter ellipses.
For M₂ (\$\\omega = 1.405 \\times 10^{-4}\$ rad s⁻¹) and \$U = 0.25\$ m s⁻¹, the
resulting body force amplitude is \$\\approx 3.5 \\times 10^{-5}\$ m s⁻², consistent
with observed tidal acceleration magnitudes on the Scotian Shelf.

# Inputs
# Inputs
- `constituents::Vector{Symbol}`: List of constituents to include (e.g. `[:M2, :S2]`).
- `u_amplitudes::Dict{Symbol, Float64}`: Zonal tidal velocity amplitudes in \$m s^{-1}\$.
- `v_amplitudes::Dict{Symbol, Float64}`: Meridional tidal velocity amplitudes in \$m s^{-1}\$.
- `phases::Dict{Symbol, Float64}`: Tidal phase offsets in radians.
- `bottom_drag_linear::Real`: Linear bottom drag rate \$r_{\\text{drag}}\$ in \$s^{-1}\$.
  When \$r_{\\text{drag}} > 0\$, the acceleration amplitude is calibrated to
  \$U_k \\sqrt{\\omega_k^2 + r_{\\text{drag}}^2}\$ to compensate for frictional damping
  and maintain the target velocity amplitude \$U_k\$.

# Outputs
- `NamedTuple`: `(u = Fu, v = Fv)` forcing functions.

# References
- Egbert, G. D., & Erofeeva, S. Y. (2002). Efficient inverse modeling of barotropic
  ocean tides. *Journal of Atmospheric and Oceanic Technology*, 19(2), 183-204.
  DOI: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
"""
function build_tidal_body_forcing(;
    constituents::Vector{Symbol} = [:M2],
    u_amplitudes::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.25),
    v_amplitudes::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.12),
    phases::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.0),
    bottom_drag_linear::Real = 0.0,
    bottom_drag::Union{Nothing, Real} = nothing
)
    freqs = [get_tidal_frequency(c) for c in constituents]
    u_amps = [Float64(get(u_amplitudes, c, 0.0)) for c in constituents]
    v_amps = [Float64(get(v_amplitudes, c, 0.0)) for c in constituents]
    phs = [Float64(get(phases, c, 0.0)) for c in constituents]
    r_linear = !isnothing(bottom_drag) ? Float64(bottom_drag) : Float64(bottom_drag_linear)

    # Scaling factor: sqrt(ω² + r²) compensates for linear bottom drag damping
    # In the absence of drag (r = 0), this simplifies to ω.
    scales = [sqrt(ω^2 + r_linear^2) for ω in freqs]

    F_u(x, y, z, t) = begin
        val = 0.0
        for i in eachindex(constituents)
            val += u_amps[i] * scales[i] * cos(freqs[i] * t + phs[i])
        end
        val
    end

    F_v(x, y, z, t) = begin
        val = 0.0
        for i in eachindex(constituents)
            val += v_amps[i] * scales[i] * sin(freqs[i] * t + phs[i])
        end
        val
    end

    return (u = F_u, v = F_v)
end

"""
    tidal_velocity_vector(
        t::Real;
        constituents::Vector{Symbol} = [:M2],
        u_amplitudes::Dict{Symbol, Float64} = Dict(:M2 => 0.25),
        v_amplitudes::Dict{Symbol, Float64} = Dict(:M2 => 0.12),
        phases::Dict{Symbol, Float64} = Dict(:M2 => 0.0)
    )

Compute instantaneous horizontal tidal velocity vector \$(u_{\\text{tide}}, v_{\\text{tide}})\$
at time \$t\$ from multi-constituent harmonic synthesis.

# Inputs
- `t::Real`: Time in seconds.
- `constituents::Vector{Symbol}`: Active tidal constituents.
- `u_amplitudes::Dict`: Zonal velocity amplitudes.
- `v_amplitudes::Dict`: Meridional velocity amplitudes.
- `phases::Dict`: Phase offsets.

# Outputs
- `Tuple{Float64, Float64}`: `(u_tide, v_tide)` in \$m s^{-1}\$.
"""
function tidal_velocity_vector(
    t::Real;
    constituents::Vector{Symbol} = [:M2],
    u_amplitudes::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.25),
    v_amplitudes::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.12),
    phases::AbstractDict{Symbol, <:Real} = Dict(:M2 => 0.0)
)
    u_tot = 0.0
    v_tot = 0.0
    for c in constituents
        omega = get_tidal_frequency(c)
        u_a = Float64(get(u_amplitudes, c, 0.0))
        v_a = Float64(get(v_amplitudes, c, 0.0))
        phi = Float64(get(phases, c, 0.0))

        u_tot += u_a * cos(omega * t + phi)
        v_tot += v_a * sin(omega * t + phi)
    end
    return (Float64(u_tot), Float64(v_tot))
end

"""
    simpson_hunter_parameter(
        water_depth::Real,
        u_tidal_amplitude::Real;
        u_wind_speed::Real = 0.0,
        cd::Real = 2.5e-3,
        buoyancy_flux::Union{Nothing, Real} = nothing
    )

Compute the Simpson-Hunter tidal mixing front parameter \$\\chi = \\log_{10}(h / U^3)\$
or generalized shear-buoyancy index to identify where tidal and wind dissipation overcomes stratification.

# Mathematical Formulation
```math
\\chi = \\log_{10}\\left( \\frac{h}{U_{\\text{tide}}^3 + \\gamma U_{\\text{wind}}^3} \\right)
```
where \$h\$ is water column depth in meters, \$U_{\\text{tide}}\$ is tidal velocity amplitude,
and \$U_{\\text{wind}}\$ represents surface wind shear dissipation
(Garrett, Keeley & Greenberg 1978; Loder & Greenberg 1986).
- \$\\chi < 1.5\$: Well-mixed water column (e.g. shallow banks, Georges Bank, Bay of Fundy).
- \$\\chi > 2.0\$: Thermally stratified shelf waters.
- \$\\chi \\approx 1.5 - 2.0\$: Tidal mixing front (nursery retention zone).

# Inputs
- `water_depth::Real`: Water column thickness in meters (\$h > 0\$).
- `u_tidal_amplitude::Real`: Peak tidal current amplitude in \$m s^{-1}\$.
- `u_wind_speed::Real`: Optional 10-meter wind speed in \$m s^{-1}\$ (default 0.0).
- `cd::Real`: Bottom drag friction coefficient (default \$2.5 \\times 10^{-3}\$).
- `buoyancy_flux::Union{Nothing, Real}`: Optional surface buoyancy flux \$B\$ in \$m^2 s^{-3}\$.

# Outputs
- `Float64`: Simpson-Hunter parameter \$\\chi\$.

# References
- Simpson, J. H., & Hunter, J. R. (1974). Fronts in the Irish Sea. *Nature*, 250, 404-406.
- Garrett, C. J. R., Keeley, J. R., & Greenberg, D. A. (1978). *Continental Shelf Research*, 18(1), 17-33.
- Loder, J. W., & Greenberg, D. A. (1986). *Continental Shelf Research*, 5(6), 679-704.
"""
function simpson_hunter_parameter(
    water_depth::Real,
    u_tidal_amplitude::Real;
    u_wind_speed::Real = 0.0,
    cd::Real = 2.5e-3,
    buoyancy_flux::Union{Nothing, Real} = nothing
)
    h = max(1.0, abs(Float64(water_depth)))
    u_tide = max(1e-3, abs(Float64(u_tidal_amplitude)))
    u_wind = max(0.0, abs(Float64(u_wind_speed)))
    # Energy dissipation: ε ~ ρ * (U_tide³ + 0.05 * U_wind³)
    dissipation_u3 = u_tide^3 + 0.05 * u_wind^3
    if !isnothing(buoyancy_flux) && buoyancy_flux > 0.0
        return Float64(log10(buoyancy_flux * h / (cd * dissipation_u3)))
    else
        return Float64(log10(h / dissipation_u3))
    end
end
