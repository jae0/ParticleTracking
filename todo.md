1.  eventually move movementanalysis here and use the temperature and salinity data to provide hsi information for snow crab

---

## Physical Correctness Assessment (2026-10-07)

### Overall: Good physical foundation with specific concerns

### 1. Hydrodynamic Model (`hydrodynamic_model.jl`) - GOOD
**Strengths:**
- TEOS-10 equation of state with ρ₀ = 1025 kg/m³ (correct for Boussinesq)
- Spherical Coriolis (`HydrostaticSphericalCoriolis`) - appropriate for 42-47.5°N where β-effect matters
- Split-explicit free surface - standard efficient choice
- Tidal forcing as **relaxation toward prescribed velocity** (correct: target velocity, not body force)
- Lateral sponge relaxation (Price & Aumont 2011) - proper formulation
- Tidal body forcing as relaxation toward tidal velocity - physically consistent

**Critical Concern - Surface Heat Flux:**
- T-dependent surface heat flux BLOCKS with StackOverflowError in `ContinuousBoundaryFunction`
- Falls back to constant flux with error if T-dependent function provided
- **Impact**: Cannot use realistic bulk heat flux (depends on SST). Known Oceananigans 0.111 limitation.
workaround:

Problem: Oceananigans 0.111 has a bug where ContinuousBoundaryFunction with field_dependencies = (:T,) causes StackOverflowError
Solution: Evaluate the T-dependent bulk flux function once at 10°C climatological SST to obtain a constant flux approximation
Warning: Clear message that this is an approximation; for time-varying bulk flux, user should provide pre-computed flux time series or upgrade Oceananigans when the bug is fixed

The proper long-term fix would require either:

Upgrading Oceananigans when the field_dependencies bug is fixed
Implementing a proper interior relaxation with a callback in the simulation loop (requires modifying ParticleTrackingRun.jl)
Providing pre-computed surface heat flux time series as input


to properly address this:

Verify if the bug is actually fixed in Oceananigans 0.113.5 by testing with a minimal ContinuousBoundaryFunction example
If fixed, update the code to use proper boundary conditions with field dependencies
Add version checking or conditional logic to maintain compatibility with older versions if needed
Remove the workaround comments and approximations
Until the bug is confirmed fixed and the code is updated, the current workarounds remain necessary to avoid GPU compilation errors.

 the Oceananigans bug were fixed in a newer version, the proper solution for surface heat flux would be to use ContinuousBoundaryFunction with appropriate field dependencies instead of the current workaround.

Here's what the proper implementation would look like:

# PROPER SOLUTION (when Oceananigans bug is fixed):
# Replace the current surface heat flux handling with:

if surface_heat_flux isa Function
    # Create a ContinuousBoundaryFunction that depends on temperature field
    # This allows the flux to vary with SST as needed
    T_flux_bc = ContinuousBoundaryFunction(
        (model, field, i, j, k, grid, clock, _) -> 
            -surface_heat_flux(
                grid, 
                model.velocities.u[i, j, k], 
                model.velocities.v[i, j, k], 
                model.velocities.w[i, j, k], 
                model.tracers.T[i, j, k], 
                model.tracers.S[i, j, k], 
                clock
            ) / rho0_cp,
        field_dependencies = (:T,)  # Key: depends on temperature field
    )
elseif surface_heat_flux isa AbstractMatrix
    FluxBoundaryCondition(-surface_heat_flux ./ rho0_cp)
else
    FluxBoundaryCondition(-Float64(surface_heat_flux) / rho0_cp)
end
This approach would:

Allow true temperature-dependent surface heat flux calculations
Maintain full coupling between SST and heat flux
Work correctly on both CPU and GPU (once the bug is fixed)
Eliminate the approximation of evaluating flux only at 10°C SST


---

### 2. Open Boundary / Sponge - INCOMPLETE (Documented)
Explicitly documented at lines 319-363:
1. **No momentum constraint** - T/S sponge only; velocity sponge uses constant `u_inflow/v_inflow`
2. **No time variability** - WOA23 climatology = static boundary, no tidal/synoptic/seasonal cycle
3. **No O₂ relaxation** if `enable_o2=true`
4. **No currents from WOA23** - needs GLORYS for real currents
**Assessment**: Honest documentation, but users must understand these are one-way nested boundaries.

---

### 3. Larval Behavior & Transport - STRONG
**Vertical Migration (DVM):** Stage-specific depths with hyperbolic tangent transitions; CIL awareness; turbidity attenuation
**Developmental Molting:** Mean-preserving lognormal thresholds (increment-based, avoids `max()` bias); per-larva quantile drawn once; `cv_molt` as pure dispersion knob
**Thermal Mortality:** Stage-specific scaling; thermal + cold stress; frailty model (lognormal) - correct unobserved heterogeneity
**Settlement:** HSI = S_z × S_T; Beta-distributed propensity; per-larva settlement quantile; competence window
**Transport Step:** Euler-Maruyama with Visser (1997) κ_v gradient drift; BBL logarithmic attenuation; passive sinking by stage; current speed clamp (3 m/s)

**Land Boundary Concern:** Only tries zonal then meridional slip; no true alongshore projection for oblique coastlines - could trap particles in corners.

---

### 4. Key Physical Inconsistencies / Risks

| Issue | Location | Severity |
|-------|----------|----------|
| Surface heat flux T-dependence broken | `hydrodynamic_model.jl:149-167` | HIGH |
| Sponge T/S only, no currents | `hydrodynamic_model.jl:319-363` | MEDIUM (documented) |
| Tidal temperature filter τ=12.42h hardcoded | `larval_behavior.jl:1652-1654` | MEDIUM |
| No Coriolis in larval transport | `larval_transport_step` | LOW |
| Land boundary: only zonal/meridional slip | `larval_behavior.jl:1143-1158` | MEDIUM |
| No Stokes drift / wave effects | Throughout | LOW |

---

### 5. Recommendations (Priority Order)

1. **Fix surface heat flux** - Use Oceananigans version where `ContinuousBoundaryFunction` with `field_dependencies` works, or implement as interior relaxation
2. **Add Coriolis to larval transport** - `f * k × u` term in `larval_transport_step`
3. **Improve land boundary** - Full tangential projection along coastline normal
4. **Consider Stokes drift** - For surface-trapped Zoea I, wave-driven Stokes drift matters
4. **Validate tidal filter** - `α = dt/44712` assumes pure M2; S2/fortnightly cycle gets filtered
5. **Add vertical diffusivity profile validation** - `κ_v_profile` should be positive-definite

---

### Bottom Line
Codebase is **physically sound in its core** (TEOS-10, proper tidal relaxation, correct molting statistics, Visser drift correction, frailty mortality). Main gaps are **Oceananigans version limitations** (surface flux) and **documented boundary incompleteness**. For Scotian Shelf snow crab application, **fit for purpose with documented caveats**.

 