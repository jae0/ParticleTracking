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

---

## Physical Correctness Assessment & Implementation Status (Updated 2026-10-07)

### Executive Summary
The core hydrodynamic solver and Lagrangian tracking pipelines have been verified and stabilized on CUDA GPU (`Oceananigans.jl` v0.111.0, CUDA.jl v6.0.0, Julia 1.13):
- **24-hour hydrodynamic validation completed successfully in 113.26 seconds on GPU**.
- Velocity fields remained fully bounded: final $\max |u| = 0.6391\text{ m/s}$, $\max |v| = 1.0292\text{ m/s}$ with adaptive CFL at $0.019$ (target $\le 0.2$), with zero numerical blowup and no artificial clamping.

---

### 1. Hydrodynamic Model (`hydrodynamic_model.jl`) - VERIFIED & RESOLVED
- **Surface Heat Flux Dynamic Coupling (RESOLVED)**:
  - The stack overflow previously seen with raw `ContinuousBoundaryFunction` was resolved by using `FluxBoundaryCondition(..., field_dependencies = :T)` in `src/model/hydrodynamic_model.jl:120-135`.
  - Kinematic heat flux $Q / (\rho_0 c_p)$ is computed with physical constants $\rho_0 = 1025.0\text{ kg/m}^3$ and $c_p = 3991.0\text{ J/(kg K)}$.
  - Verified to execute dynamically on both CPU and CUDA without blocking or lowering errors, allowing real SST-dependent heat exchange.
- **Empty / Non-Sponge Forcing (RESOLVED)**:
  - Replaced undefined `ZeroForcing()` placeholder with an idiomatic empty `NamedTuple()`.
- **Options Struct Synchronization (RESOLVED)**:
  - Added concrete `input_dir` and `output_dir` fields to `HydrodynamicOptions` across struct definitions and positional/keyword constructors.

---

### 2. Audit Clarifications & Physical Revisions
- **Clarification on Coriolis in Lagrangian Step (CORRECTED)**:
  - A prior audit recommended adding an explicit Coriolis acceleration ($f \cdot \mathbf{k} \times \mathbf{u}$) to `larval_transport_step`.
  - **Correction**: This was analytically incorrect. The Eulerian velocity field $\mathbf{u}$ produced by Oceananigans already includes the Coriolis acceleration. Adding Coriolis to Lagrangian kinematic advection ($d\mathbf{x}/dt = \mathbf{u}$) would double-count Coriolis forces and distort particle trajectories.
- **GeoData Multi-Format Stores (RESOLVED)**:
  - Updated `ParticleTrackingRun.jl` (Segments 1, 2, and 3) to uniformly discover both `.zarr` and `.nc` formats for bathymetry, surface winds, and WOA23 hydrography, matching GeoData's unified storage catalog.

---

### 3. Open Tasks & Remaining Priorities
1. **Open Boundary Dynamic Inflow Coupling**:
   - Replace constant inflow velocities ($u_{\text{inflow}} = -0.15\text{ m/s}$) with 3D Field targets (`CenterField(grid)`) populated from GLORYS reanalysis when dynamic velocity boundary conditions are enabled.
2. **Coastline Normal Projection for Larvae**:
   - Enhance land boundary collision handling in `larval_behavior.jl` with tangential projection along the local normal to prevent trapping in acute coastal embayments.
3. **Stokes Drift & Wave Coupling**:
   - Wire wave-driven Stokes drift into upper-column (0-20 m) Zoea I transport when wave fields are present.

 



 