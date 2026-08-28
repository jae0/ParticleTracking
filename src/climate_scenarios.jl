"""
    climate_scenarios.jl

Climate forcing scenarios (CMIP6 SSPs, downscaled regional projections, marine heatwaves)
and temperature-dependent larval biology (pelagic larval duration, thermal mortality,
and thermal DVM suppression) for snow crab (*Chionoecetes opilio*).
"""

using Oceananigans

"""
    get_climate_scenario_deltas(
        scenario::Symbol = :ssp245;
        year::Int = 2050
    )

Retrieve projected thermal, haline, and atmospheric wind stress anomaly deltas for the
Scotian Shelf and Northwest Atlantic domain under CMIP6 / IPCC AR6 climate scenarios.

# Mathematical Formulation
Anomalies \$\\Delta T_{\\text{surf}}, \\Delta T_{\\text{deep}}, \\Delta S_{\\text{surf}}\$
are scaled linearly relative to a 2015 baseline:
```math
\\Delta X(t) = \\Delta X_{2050} \\left( \\frac{t - 2015}{2050 - 2015} \\right)
```

# Inputs
- `scenario::Symbol`: Climate scenario identifier (`:historical`, `:ssp126`, `:ssp245`,
  `:ssp370`, `:ssp585`, `:marine_heatwave`).
- `year::Int`: Target projection year (e.g. 2050, 2100).
- `baseline_year::Real`: Baseline climatological reference year (default 2015.0).
- `horizon_year::Real`: Benchmark projection horizon year for deltas (default 2050.0).

# Outputs
- `NamedTuple`: `(ΔT_surface, ΔT_deep, ΔS_surface, Δwind_factor, description)`

# References
- Brickman, D., Wang, Z., & DeTracey, B. (2018). Variability and trends in the
  Scotian Shelf and Gulf of Maine region from a high-resolution regional ocean
  climate model. *Progress in Oceanography*, 164, 49-64.
  DOI: 10.1016/j.pocean.2018.04.004
- Saba, V. S., et al. (2016). Enhanced warming of the Northwest Atlantic Ocean
  under climate change. *Journal of Geophysical Research: Oceans*, 121(1), 118-132.
  DOI: 10.1002/2015JC011346
- O'Neill, B. C., et al. (2016). The Scenario Model Intercomparison Project
  (ScenarioMIP) for CMIP6. *Geoscientific Model Development*, 9(9), 3461-3482.
  DOI: 10.5194/gmd-9-3461-2016
"""
function get_climate_scenario_deltas(
    scenario::Symbol = :ssp245;
    year::Int = 2050,
    baseline_year::Real = 2015.0,
    horizon_year::Real = 2050.0
)
    span = horizon_year - baseline_year
    time_factor = span > 0.0 ? clamp((Float64(year) - Float64(baseline_year)) / span, 0.0, 5.0) : 1.0

    if scenario == :historical || scenario == :baseline
        return (
            ΔT_surface = 0.0,
            ΔT_deep = 0.0,
            ΔS_surface = 0.0,
            Δwind_factor = 1.0,
            description = "Historical / present-day climatological baseline"
        )
    elseif scenario == :ssp126
        return (
            ΔT_surface = 1.1 * time_factor,
            ΔT_deep = 0.5 * time_factor,
            ΔS_surface = -0.25 * time_factor,
            Δwind_factor = 1.0 + 0.05 * time_factor,
            description = "CMIP6 SSP1-2.6 (Low emissions / Paris target)"
        )
    elseif scenario == :ssp245
        return (
            ΔT_surface = 1.8 * time_factor,
            ΔT_deep = 0.9 * time_factor,
            ΔS_surface = -0.45 * time_factor,
            Δwind_factor = 1.0 + 0.10 * time_factor,
            description = "CMIP6 SSP2-4.5 (Intermediate / Middle of the road)"
        )
    elseif scenario == :ssp370
        return (
            ΔT_surface = 2.6 * time_factor,
            ΔT_deep = 1.4 * time_factor,
            ΔS_surface = -0.65 * time_factor,
            Δwind_factor = 1.0 + 0.15 * time_factor,
            description = "CMIP6 SSP3-7.0 (High emissions / Regional rivalry)"
        )
    elseif scenario == :ssp585
        return (
            ΔT_surface = 3.5 * time_factor,
            ΔT_deep = 1.9 * time_factor,
            ΔS_surface = -0.85 * time_factor,
            Δwind_factor = 1.0 + 0.20 * time_factor,
            description = "CMIP6 SSP5-8.5 (Very high emissions / Fossil-fueled)"
        )
    elseif scenario == :marine_heatwave || scenario == :mhw
        return (
            ΔT_surface = 3.5,
            ΔT_deep = 0.2,
            ΔS_surface = -0.1,
            Δwind_factor = 0.8,
            description = "Transient Category III/IV Marine Heatwave (MHW)"
        )
    else
        error(
            "Unknown climate scenario '$(scenario)'. " *
            "Available: :historical, :ssp126, :ssp245, :ssp370, :ssp585, :marine_heatwave"
        )
    end
end

"""
    apply_climate_scenario!(
        model;
        scenario::Symbol = :ssp245,
        year::Int = 2050,
        baseline_surface_T::Real = 15.0,
        baseline_grad_T::Real = 0.01,
        baseline_S::Real = 35.0,
        mixed_layer_depth::Real = 30.0
    )

Apply temperature and salinity climate anomalies to modify the initial stratification
of a hydrodynamic model.

# Mathematical Formulation
The vertical climate anomaly profiles \$\\Delta T(z)\$ and \$\\Delta S(z)\$ are:
```math
\\Delta T(z) = \\Delta T_{\\text{deep}} + (\\Delta T_{\\text{surface}} - \\Delta T_{\\text{deep}}) \\exp\\left( \\frac{z}{H_{\\text{mix}}} \\right)
```
```math
\\Delta S(z) = \\Delta S_{\\text{surface}} \\exp\\left( \\frac{z}{H_{\\text{mix}}} \\right)
```

# Inputs
- `model`: Configured Oceananigans `HydrostaticFreeSurfaceModel`.
- `scenario::Symbol`: Climate projection scenario.
- `year::Int`: Target year.
- `baseline_surface_T::Real`: Baseline sea surface temperature in °C.
- `baseline_grad_T::Real`: Baseline vertical temperature gradient in °C/m.
- `baseline_S::Real`: Baseline practical salinity in PSU.
- `mixed_layer_depth::Real`: Epipelagic transition depth \$H_{\\text{mix}}\$ in meters.

# Outputs
- `NamedTuple`: Applied climate delta values.

# References
- Brickman, D., Wang, Z., & DeTracey, B. (2018). *Progress in Oceanography*, 164, 49-64.
- Loder, J. W., van der Baaren, A., & Yashayaev, I. (2015). *Can. Tech. Rep. Hydrogr. Ocean Sci.*, 305, 142 pp.
"""
function apply_climate_scenario!(
    model;
    scenario::Symbol = :ssp245,
    year::Int = 2050,
    baseline_surface_T::Real = 15.0,
    baseline_grad_T::Real = 0.01,
    baseline_S::Real = 35.0,
    mixed_layer_depth::Real = 30.0
)
    deltas = get_climate_scenario_deltas(scenario, year = year)

    # Temperature profile under climate scenario
    climate_T(x, y, z) = begin
        base_t = baseline_surface_T + baseline_grad_T * z
        anom_t = deltas.ΔT_deep + (deltas.ΔT_surface - deltas.ΔT_deep) *
                 exp(z / mixed_layer_depth)
        base_t + anom_t
    end

    # Salinity profile under climate scenario
    climate_S(x, y, z) = begin
        base_s = baseline_S
        anom_s = deltas.ΔS_surface * exp(z / mixed_layer_depth)
        base_s + anom_s
    end

    set!(model, T = climate_T, S = climate_S)
    return deltas
end

"""
    temperature_dependent_pld(
        temperature_celsius::Real;
        a::Real = 135.0,
        b::Real = 0.75,
        t_ref::Real = 1.0
    )

Calculate total Pelagic Larval Duration (PLD in days) for snow crab (*Chionoecetes opilio*)
from egg hatch to benthic settlement as a function of ambient water temperature.

# Mathematical Formulation
Following empirical rearing models (Sainte-Marie & Sainte-Marie, 1999; Kuhn & Choi, 2011):
```math
\\text{PLD}(T) = a \\cdot (T + t_{\\text{ref}})^{-b}
```

# Inputs
- `temperature_celsius::Real`: Mean ambient water temperature in °C.
- `a::Real`: Empirical scaling constant (default 135.0).
- `b::Real`: Empirical power coefficient (default 0.75).
- `t_ref::Real`: Temperature offset parameter (default 1.0 °C).

# Outputs
- `Float64`: Estimated larval drift duration in days.

# References
- Kuhn, P. S., & Choi, J. S. (2011). Influence of temperature on embryo incubation
  and larval development in snow crab (*Chionoecetes opilio*).
  *Fisheries Research*, 107(1-3), 81-87. DOI: 10.1016/j.fishres.2010.10.011
- Sainte-Marie, G., & Sainte-Marie, B. (1999). Growth, developmental stages, and
  vertical distribution of snow crab larvae (*Chionoecetes opilio*).
  *Can. J. Fish. Aquat. Sci.*, 56(11), 2181-2193. DOI: 10.1139/f99-151
"""
function temperature_dependent_pld(
    temperature_celsius::Real;
    a::Real = 135.0,
    b::Real = 0.75,
    t_ref::Real = 1.0
)
    t_effective = max(0.0, Float64(temperature_celsius)) + t_ref
    pld_days = a * (t_effective)^(-b)
    return Float64(pld_days)
end

"""
    larval_thermal_mortality_rate(
        temperature_celsius::Real;
        base_mortality::Real = 0.02,
        thermal_threshold::Real = 10.0,
        thermal_sensitivity::Real = 0.015
    )

Compute the instantaneous daily mortality rate \$\\mu(T)\$ (day⁻¹) for snow crab larvae
accounting for thermal stress when temperatures exceed tolerance thresholds.

# Mathematical Formulation
```math
\\mu(T) = \\mu_{\\text{base}} + \\mu_{\\text{thermal}} \\max(0, T - T_{\\text{crit}})^2
```

# Inputs
- `temperature_celsius::Real`: Water temperature in °C.
- `base_mortality::Real`: Baseline natural daily mortality rate (default 0.02 day⁻¹).
- `thermal_threshold::Real`: Upper thermal stress threshold \$T_{\\text{crit}}\$ (default 10.0 °C).
- `thermal_sensitivity::Real`: Thermal penalty coefficient (default 0.015 day⁻¹ °C⁻²).

# Outputs
- `Float64`: Instantaneous daily mortality rate in \$\\text{day}^{-1}\$.

# References
- Epifanio, C. E., & Cohen, J. H. (2016). *Journal of Experimental Marine Biology
  and Ecology*, 482, 85-105.
"""
function larval_thermal_mortality_rate(
    temperature_celsius::Real;
    base_mortality::Real = 0.02,
    thermal_threshold::Real = 10.0,
    thermal_sensitivity::Real = 0.015
)
    excess_t = max(0.0, Float64(temperature_celsius) - thermal_threshold)
    mortality = base_mortality + thermal_sensitivity * (excess_t^2)
    return Float64(mortality)
end
