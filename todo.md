# ParticleTracking — TODO

Continuation notes. Supersedes the earlier `REFACTOR_PLAN.md` (deleted 2026-09-26).

State at last update: **703/703 tests pass** (as of the last full run; the suite is now **parked**
by the standing constraint below). `configs/default.toml` runs green end to end on CPU.

`configs/opendata.toml` gets through Segments 1–2 against real data and now has WOA23 T/S/O₂
downloaded and verified; it is **parked** at Segment 3+ for the same reason. Earlier notes claimed
these configs "fail for lack of network access" — **that was wrong**. The network is partially
open, and the failures were code bugs: a `MethodError` before any request, wrong dataset ids, and
a reader that assumed 3-D data. See [section 5](#5-woa23-is-downloaded-and-working-and-the-network-blocker-was-never-the-blocker).

The GPU `InvalidIRError` has been **traced to a specific cause** and a verified replacement found
— see [section 2](#2-gpu-root-cause-found-field_dependencies-on-a-continuousforcing).

---

## Recently completed

| Area | Change | Where |
|---|---|---|
| Config | Driver's duplicate ~80-key TOML parser deleted; it now calls `configuration_to_options` and reads fields off the result. The TOML is interpreted in exactly one place. | `ParticleTrackingRun.jl:2218-2326` |
| Config | 8 keys the driver silently dropped now flow through: `boundaries.method`, `sponge_layer_width`, `sponge_timescale`, `u_inflow`, `v_inflow`, `hydrodynamics.max_dt_seconds`, `tides.s2_u_amp`, `s2_v_amp`. | `ParticleTrackingRun.jl` |
| Config | `options_to_configuration` no longer hardcodes `sponge_layer_width`/`sponge_timescale` and no longer writes `max_dt_seconds = 90.0`; all now read from `opts`, so `resolved_config.toml` is a faithful record. | `src/configuration.jl:873` |
| Config | New `hydrography_source` / `hydrography_month` options, read from `[data]`. | `src/configuration.jl` |
| Data | `woa23_regridded_tracers`, `read_woa_variable`, `build_woa_interpolator`: read WOA23 NetCDF and regrid onto the model grid (keyless NOAA NCEI). | `src/hydrodynamic_model.jl` |
| Data | New `configs/opendata.toml`: fully keyless run (ETOPO 2022 + Open-Meteo ERA5 winds + WOA23 + TPXO9). | `configs/opendata.toml` |
| Physics | Native `:ρ` and `:N2` diagnostics (TEOS-10, model grid), replacing post-hoc finite differences. | `src/simulation.jl` |
| Runtime | Flow interpolator reads the JLD2 sidecar; no longer falls back to an analytical jet. | `src/simulation.jl:1116` |

### Bugs found and fixed while doing the above

- **`configs/default.toml` documented `method = "none"` but ran with the sponge enabled.** The
  driver read neither `boundaries.method` nor the sponge parameters, so the constructor defaults
  (`:relaxation`, width 0.35) applied. This is why earlier notes referred to a "boundary/sponge not
  wired" state that was partly self-inflicted.
- **`fill_seawater_density!` scalar-indexed device arrays.** The TEOS-10 polynomial carries a
  coefficient vector, so `eos` is not `isbits` and cannot be captured in a GPU kernel; the fill is
  now done on the host and copied back.
- **`fill_seawater_density!` could throw `DomainError` and kill the run.** Uninitialised halo values
  of `T`/`S` pushed the polynomial outside its validity box, producing a huge negative discriminant
  for its internal `sqrt`. Inputs are now checked and clamped, with a reference-density fallback.
  This one also silently corrupted halo density feeding `∂z`.
- **JLD2 flow interpolator "missing sidecar" error was false.** The sidecar was always present; the
  `domain_lon` fallback was evaluated *eagerly* and threw even when the sidecar was about to be used.
- **EOS argument order was swapped** (`ρ(S, T, …)` instead of `ρ(T, S, …)`), giving ≈1004 kg m⁻³
  where shelf conditions must read ≈1027.
- **Density halos were never filled** because `axes` on a `Field` reports the *interior* while the
  data spans the halo. N² was wrong by ~1000× at the first/last interior cells.
- **`set!` has no NamedTuple method**; `set!(model, (T=…, S=…))` is a `MethodError`.
- **`interpolate` indexes by index coordinates**, not physical lon/lat/depth.
- **`fill_seawater_density!` read uninitialised halo memory** and could segfault the interpreter
  (`EXCEPTION_ACCESS_VIOLATION`). The buffer now starts filled with the EOS reference density and
  only the interior is evaluated, with the interior range derived from `axes` and `Nx/Ny/Nz`
  (`indices(field, k)` returns a `Colon` in 0.111 and is useless for this). This was the cause of
  the `default.toml --gpu` crash; that config now completes on GPU.
- **`apply_prescribed_surface_state!` was not merely a no-op — it could not run.** It referenced
  `model.boundary_conditions`, a field `HydrostaticFreeSurfaceModel` does not have, so it would
  have raised `FieldError`. The old body also passed all seven existing BC components straight
  back through `FieldBoundaryConditions`, and discarded the `sst0`/`sss0` it had just computed.
  Rewritten to validate shapes and write the top interior layer via `apply_to_surface_layer!`;
  the `Symbol`-free, version-stable route is documented, and `prescribed_ocean_component` was
  added for the supported `EarthSystemModel` + `PrescribedOcean` route. Verified: surface changes,
  interior preserved, 3-D input takes the first slice, bad shape raises.
- **`apply_climate_scenario!` replaced the baseline instead of perturbing it.** It called
  `set!(model, T = T₀ + γz + ΔT)` from hardcoded defaults (`T₀=15`, `γ=0.01`, `S=35`), so a
  scenario run silently discarded the stratification set immediately before it — including real
  WOA23/GLORYS hydrography — and its `climate_T(x, y, z)` ignored `x`/`y` entirely, erasing any
  cross-shelf or CIL structure. Now adds the surface-weighted anomaly to the model's actual
  state. Verified: ΔT grows towards the surface and floors at `ΔT_deep`, the anomaly is uniform
  *within* each layer, horizontal structure is bit-for-bit preserved, and the baseline salinity
  stays 33 rather than snapping to 35.
- **`LateralBoundaryRelaxation` took a `::Symbol` argument** in the GPU forcing path. A `Symbol` is
  not `isbits` and cannot cross a CUDA kernel ABI, so the enclosing `ContinuousForcing` was
  unlowerable. Added `sponge_relaxation_u` / `sponge_relaxation_v` with no selector and routed the
  forcing through them. This is a genuine latent GPU bug, but it is **not** the `snowcrab` cause —
  see section 2, where the trigger turned out to be `field_dependencies`.
- **`build_hydrodynamic_model` silently discarded a caller-supplied `forcing`.** With the sponge
  inactive it builds `composed_forcing = (; u = ZeroForcing(), v = ZeroForcing())` and then does
  `merge(forcing, composed_forcing)`, so the caller's `u`/`v` forcings are overwritten by zero
  forcing. A caller can hand in a perfectly good `ContinuousForcing` and get nothing, with no
  warning. Found while building a GPU reproducer that could not influence the model at all. Only
  the tide and sponge branches (which assign real forcings) survive the merge. Fix: only fill in
  keys that are actually absent, or merge the other way round.
- **A `ContinuousForcing` with `field_dependencies` does not lower on GPU** in the split-explicit
  free-surface model — the actual `snowcrab.toml --gpu` failure. Full isolation and the verified
  replacement are in section 2.

### Source layout

The package was reorganised into `src/{utils,data,model,config,biology,analysis,output}/`. Files
were **moved**, not rewritten, so all working code is preserved byte-for-byte; the package is a
single flat namespace, so these remain plain `include`s in dependency order rather than submodules.

---

## Open items

### Remaining tasks, in priority order

| # | Task | Status | Ref |
|---|---|---|---|
| 1 | GPU `InvalidIRError` — migrate the sponge to `Relaxation` | **DONE and verified**; migration is numerically exact, all 4 GPU combinations compile | [2](#2-gpu-is-fixed-the-sponge-now-uses-relaxation-and-all-combinations-compile) |
| 2 | Run the full test suite | not run since the biology + GPU + data changes | — |
| 3 | Complete `opendata.toml` / `snowcrab.toml` end to end | segments 1–2 pass; 3+ not reached since the WOA23 reader fix | [7](#7-network-dependent-configs) |
| 4 | Implement the GLORYS reader | **required**. Access proven anonymously; reader not written | [6](#6-glorys-support-is-required-and-access-is-solved) |
| 5 | Replace `heat_flux = 50.0` with real data | **required**. Bulk formula from cached wind needs no credentials; ERA5 needs a CDS key | [3](#3-atmospheric-forcing-is-a-placeholder-and-real-data-is-required) |
| 6 | Calibrate `cv_molt` / `cv_mortality` / `cv_settlement` off the `0.25` placeholder | functional form now lognormal / lognormal / Beta, verified; the **numbers** are still placeholders | [1](#1-larval-biology-stochasticity-is-wired-and-two-latent-bugs-are-fixed) |
| 7 | Confirm a full GPU production run reaches completion | compiles and constructs on GPU; never run end to end on GPU | [2](#2-gpu-is-fixed-the-sponge-now-uses-relaxation-and-all-combinations-compile) |
| 8 | Re-test Julia 1.12 | `julia +1.12` installed; earlier failure blamed on a missing network, which is now known to be wrong. Lower value now that the GPU bug is root-caused | [2](#2-gpu-is-fixed-the-sponge-now-uses-relaxation-and-all-combinations-compile) |
| 9 | Consider a horizontal nearest-wet fill for fully-masked WOA columns | optional refinement; currently the window mean | [5](#5-woa23-is-downloaded-and-working-and-the-network-blocker-was-never-the-blocker) |
| 10 | Remap the two deleted Copernicus dataset IDs | `GLOBAL_MULTIYEAR_PHY_001_033` and `GLOBAL_REANALYSIS_PHY_001_031` no longer exist | [6](#6-glorys-support-is-required-and-access-is-solved) |

**Settled decisions (2026-09-27)**

- GPU: proceed with any reasonable workaround — **done**, `Relaxation` migration is exact.
- Masked WOA cells: **NaN**, then repaired — done, see [5](#5-woa23-is-downloaded-and-working-and-the-network-blocker-was-never-the-blocker).
- GLORYS support: **required** — see [6](#6-glorys-support-is-required-and-access-is-solved).
- Atmospheric forcing: **use real data** — see [3](#3-atmospheric-forcing-is-a-placeholder-and-real-data-is-required).
- `work/baseline.jld2`: **deleted** (already absent; nothing over 50 MB remains in `work/`).
- Biology dispersion: **lognormal**, with a documented Beta exception for the bounded settlement
  index — see [1](#1-larval-biology-stochasticity-is-wired-and-two-latent-bugs-are-fixed).

### 1. Larval biology: stochasticity is wired, and two latent bugs are fixed

**Status: fixed and verified. Re-audited 2026-09-27.**

The earlier version of this section claimed growth and mortality were deterministic by design and
listed three open design decisions. That was stale — the mechanism already existed for all three
processes but was never switched on. Auditing the four layers found the wiring gap *and* two real
bugs in the stochasticity code, which were dormant only because the CV values defaulted to `0.0`.

**What was broken and is now fixed**

1. **The driver never forwarded the parameters.** Both `track_larval_cohort` call sites
   (`ParticleTrackingRun.jl:1125`, `:1809`) omitted `cv_molt`, `cv_mortality`, `cv_settlement` and
   `settlement_stochastic`, so `track_larval_cohort`'s `0.0` defaults applied. Worse,
   `resolved_config.toml` recorded `0.25` — the provenance file asserted a stochastic model that no
   run ever used. Both call sites now pass all four from `opts`.
2. **Molt thresholds were redrawn on every call.** `update_larval_stage` drew a fresh threshold
   triple from `rng` *inside* the function, and the time loop called it once per step. A larva's
   developmental rate was therefore resampled every 5 minutes: a stage could regress
   (`:zoea2` → `:zoea1`) and settlement competence could flicker on and off. Measured effect before
   the fix: the mean first-molt degree-day was **32.8** instead of the base **65**, because the
   scan effectively sampled a fresh threshold at every step. Fixed by adding
   [`draw_molt_thresholds`](@ref) and `stage_from_thresholds`, drawing each particle's triple
   **once** at cohort initialisation and passing it back via a new `thresholds` keyword. Measured
   after: mean **65.08**, all stage sequences monotone.
3. **Mortality was applied twice.** The mean-field decay `traj_surv *= exp(-μ·dt)` ran
   unconditionally, and the Bernoulli draw then removed survival a second time. The two modes are
   now mutually exclusive, and in stochastic mode `traj_surv` holds until death so the recorded
   value is an outcome rather than an expectation.

**Where the stochasticity lives** (`src/biology/larval_behavior.jl`):

| Process | Mechanism | Gate |
|---|---|---|
| Growth / molting | per-particle thresholds `Normal(dd, dd·cv_molt)`, strictly increasing, drawn once | `cv_molt > 0` |
| Mortality | per-particle `Normal(μ, μ·cv_mortality)` clamped at 0, then a Bernoulli death draw | `cv_mortality > 0` |
| Settlement | `p_accept = clamp(rand(Normal(hsi, hsi·cv_settlement)), 0, 1)` then a Bernoulli draw | `cv_settlement > 0` and `stochastic` |
| Turbulent transport | three `randn` diffusion increments | always |
| Initial positions / depths | `rand` | always |

`draw_molt_thresholds` and `stage_from_thresholds` are new and exported.

**Config layer** (`src/config/configuration.jl`) was already correct: `settlement_stochastic`,
`cv_molt`, `cv_mortality` and `cv_settlement` are struct fields (`:262-264`), read from TOML
`[biology]` with a `0.25` fallback (`:726-729`), and written to `resolved_config.toml`
(`:970-973`). The four keys were missing from all three TOML configs and have been added to
`[biology]` in `configs/default.toml`, `configs/snowcrab.toml` and `configs/opendata.toml`.

**Verification**

- Test suite **703/703 pass** (up from 661; 42 new tests in testset `6b`).
- `configs/default.toml` runs green end to end (140 s, 25 particles).
- Live cohort check (`work/verify_stochastic.jl`, 60 particles, 40 d):
  - `cv_molt = 0` → **1** distinct first-molt degree-day, identical across seeds (deterministic).
  - `cv_molt = 0.4` → **59** distinct values, sd ≈ 25.5, mean **64.54** vs base 65.
  - all stage sequences monotone non-decreasing (the regression is gone).
  - `cv_mortality = 0` → 1 unique final `traj_surv` (0.4685), 0 deaths — mean-field.
  - `cv_mortality = 0.5` → 2 unique values, 37/60 dead — independent per-particle outcomes.
- The `6b` testset ends with a source-level guard asserting every `track_larval_cohort(` call site
  in the driver also passes `cv_molt = opts.cv_molt`, so this gap cannot silently return.

**Dispersion is now lognormal (molt, mortality) and Beta (settlement) — implemented and verified.**

The `0.25` values remain placeholders, but the *functional form* is no longer a Normal, which
mattered because the Normal was actively wrong in two ways.

New helpers, all exported, all exactly identity at `cv = 0`:

- `lognormal_sigma(cv)` — `σ = sqrt(log1p(cv²))`.
- `draw_lognormal_mean(base, cv, rng)` — `base * exp(σ·z − σ²/2)`. The `−σ²/2` term makes it
  **mean preserving**, so `E[draw] = base` and `CV[draw] = cv` exactly. Without it the mean would
  be `base·exp(σ²/2)` and turning a `cv` up would also raise the expected rate.
- `draw_beta_index(base, cv, rng)` — for the settlement index; see below.

**`cv_molt` — dispersion is applied to the stage increments, not the cumulative thresholds.** The
base thresholds are decomposed into increments `(b₁, b₂−b₁, b₃−b₂)` = `(65, 65, 70)` degree-days,
each increment is drawn lognormally, and the thresholds are the running sums. This matters because
both obvious alternatives were measured and are wrong:

- Drawing the **cumulative thresholds** independently and forcing `t₂ > t₁`, `t₃ > t₂` with `max`
  looks equivalent but is not: the `max` guard propagates a slow-developer tail upward, inflating
  the mean of `t₃` to **217 degree-days instead of 200** at `cv = 0.5`, and **234 at `cv = 0.75`**.
  `cv_molt` would be a knob that silently *lengthens* the expected time to settlement.
- Drawing a Normal and clamping at zero piled an artificial spike on the floor: **2.4 %** of larvae
  got a Zoea I duration of about two hours at `cv = 0.5` (**9.5 %** at `cv = 0.75`).

With increments, ordering is automatic, nothing is clamped, and all three means are exact
(`E[tₖ] = bₖ`). The CV of `t₃` is naturally *below* `cv_molt` because variability averages out
across stages — correct behaviour, and a useful diagnostic.

**`cv_mortality` — lognormal on the individual rate.** Positive by construction, so no pile-up of
zero-rate individuals, and mean preserving so the expected mortality rate is unchanged.

**`cv_settlement` — a Beta, not a lognormal.** This one deviates from the instruction, deliberately.
The suitability index is a probability confined to `[0, 1]`, so a lognormal has the wrong support
outright. The code previously used `clamp(rand(Normal(hsi, cv·hsi)), 0, 1)`, which has two defects:
* It is **not mean preserving.** At the optimum (`hsi = 1`) about half the mass falls above 1 and
  is clipped to exactly 1 while the rest is clipped down, giving an expected acceptance of **0.90
  instead of 1.00** (measured). So raising `cv_settlement` *reduced* the settlement rate — a
  dispersion knob silently changing the mean.
* Support mismatch produces mass piling on the 0 and 1 boundaries.

A logit-normal was tried and rejected: its probability-space CV is a saturating, non-monotone
function of the log-odds sd (measured implied CV 0.25 / 0.416 / 0.541 for requested 0.25 / 0.5 /
0.75 at `p = 0.5`), so `cv` would not mean what it says. A **Beta** hits both the mean and the CV
exactly: `κ = (1−μ)/(μ·cv²) − 1`, `α = μκ`, `β = (1−μ)κ`. Measured implied CV 0.2499 / 0.5001 /
0.7503 against requested 0.25 / 0.5 / 0.75, means 0.4997–0.5002 at `p = 0.5`, and the old
suppression is gone (at `hsi = 1.0` the new expected acceptance is 1.00, versus 0.90 before).

One consequence to keep in mind when reading results: a variable on `[0, 1]` has
`var ≤ μ(1−μ)`, so a bounded index of mean `μ` **cannot** have CV above `sqrt((1−μ)/μ)` — 0.33 at
`μ = 0.9`, 1.0 at `μ = 0.5`. That is a property of the interval, not a modelling choice. `κ` is
floored so an over-large `cv` degrades to a near-Bernoulli draw rather than a negative
concentration, and the mean is still preserved.

**Verification** — `work/test_lognormal.jl`, 400 000 draws per case, all passing:

- mean preserving and CV exact at `cv` = 0.1, 0.25, 0.5, 0.75, 1.0 (e.g. `cv = 0.5` gives mean
  99.87 and empirical CV 0.4999 on a base of 100);
- strictly positive, with mean > median (right-skewed, as developmental times are);
- no boundary pile-up: the most frequent single value appears < 50 times out of 400 000, where a
  clamp spike would be thousands;
- all three cumulative thresholds strictly increasing, all three means matching 65 / 130 / 200;
- `cv = 0` returns the exact base thresholds; non-increasing bases throw a clear `ArgumentError`;
- `cv_molt = 0.4` cohort run: 59 distinct first-molt degree-days, mean 64.05 vs base 65, all stage
  sequences monotone; `cv_mortality = 0.5` gives 38/60 deaths with 2 distinct survival outcomes,
  while `cv = 0` gives the deterministic mean-field path unchanged.

**Still open.** The `0.25` values are still placeholders. Two useful bridges for whoever calibrates:
`CV(threshold) = CV(stage duration)` for a cohort at fixed temperature, since `DD = (T − T_base)·d`;
and the total-requirement CV is necessarily below the per-stage CV.

### 2. GPU is fixed, the sponge now uses `Relaxation`, and all combinations compile

**Status: fixed and verified on the GPU. 2026-09-27.**

> Superseded constraint: GPU work was originally parked pending sign-off because the migration
> changes how the sponge is applied. That approval was given, and the migration turned out to be
> **numerically exact**, so no result changes.

`snowcrab.toml --gpu` failed with `InvalidIRError` in
`Oceananigans.…HydrostaticFreeSurfaceModels.gpu_compute_hydrostatic_free_surface_Gu!`. Isolated
with `work/gpu_deps.jl`, which builds a tiny GPU grid and constructs the model four ways. Kernel
compilation happens during **construction**, so construction success is the signal:

| Case | `field_dependencies` | Result |
|---|---|---|
| A | `(:u, :v)` — instantaneous velocities | **FAIL** `InvalidIRError` |
| B | none (4-arg forcing) | **PASS** |
| C | `(:T,)` — a tracer | **FAIL** `InvalidIRError` |
| D | `Relaxation(rate, mask, target)` instead | **PASS** |

**The trigger was `field_dependencies` itself, not velocity specifically.** Any `ContinuousForcing`
that declares field dependencies produces invalid LLVM IR in the split-explicit free-surface kernel
on this stack (Oceananigans 0.111.0, CUDA.jl 6.0.0, Julia 1.13.1). The regularised forcing carries
`field_dependencies_indices = Tuple{Int64, Int64}` and interpolation operators (`identity4`,
`ℑxyᶠᶜᵃ` for A; `ℑxᶠᵃᵃ` for C); the no-dependency case carries `Tuple{}` for both and compiles.

Ruled out with evidence:

- **Our sponge code was innocent.** `isbitstype(LateralBoundaryRelaxation) = true`, and `code_llvm`
  on `compute_sponge_gamma` and the `sponge_relaxation_u`/`_v` closures yields valid,
  paren-balanced IR (3054 / 3263 bytes). `work/gpu_ir.jl`.
- **The tidal forcing was innocent.** `TidalBodyForcingU/V` inside a `ContinuousForcing` with no
  dependencies constructs fine — that is case B.
- **`configs/default.toml --gpu` passing was never evidence.** That config sets
  `boundaries.method = "none"` and `tides.constituents = []`, so it never builds a forcing with
  field dependencies.

**The fix, and why it is exact.** The sponge now uses Oceananigans' own primitive,
`Relaxation(; rate, mask, target)`, which contributes

    rate * mask(X) * (target(X, t) − field)

With `rate = 1/tau_relax` and `mask = compute_sponge_gamma(…)` that is algebraically the sponge
tendency `−gamma·(psi − ref)/tau_relax`. It reads the field it relaxes through `r.relaxed`, which
the materialiser rebinds to the model field and `Adapt`s per device, so **no dependency
declaration is needed at all** — which is precisely what the broken code path was.

`Relaxation` takes over the `u`/`v` forcing entry, so the tidal body forcing is kept as a separate
additive term via `MultipleForcings(Relaxation(…), ContinuousForcing(tide))`. The total is
therefore numerically identical to the previous construction,
`tide + (1/tau)·gamma·(ref − u)`. Three candidate constructions were tested on the GPU
(`work/gpu_forcing_combo.jl`) and all compile; this one was chosen because the alternatives were
not exact — folding the tide into the `target` instead would have made the tide a target offset
rather than a tendency, which is a different model.

**Verified through the real `build_hydrodynamic_model` on GPU** (`work/gpu_verify_sponge.jl`):

| Case | Result |
|---|---|
| no sponge, no tides (previously passing) | **PASS** |
| sponge only | **PASS** |
| sponge + tides (the `snowcrab` shape, previously `InvalidIRError`) | **PASS** |
| sponge + tides + caller-supplied `forcing` | **PASS**, caller forcing preserved |

So `use_gpu = true` is no longer known-broken. A full GPU production run has still not been
executed, so this is "compiles and constructs on GPU", not "ran to completion on GPU".

Reproduction scripts (all in `work/`, untracked):

- `work/gpu_deps.jl` — the A/B/C isolation that found the root cause.
- `work/gpu_relax.jl` — case D, `Relaxation` alone.
- `work/gpu_forcing_combo.jl` — the three forcing combinations, all compiling.
- `work/gpu_verify_sponge.jl` — the four-case verification through the real builder.
- `work/gpu_ir.jl` — proves the sponge function's own IR is valid.
- `work/gpu_isolate.jl` — end-to-end across sponge × tide combinations.

### 3. Atmospheric forcing is a placeholder, and real data is required

**Decision: use real data. Two routes, and the first needs no credentials.**

`atmo_craft = (stress_x = tau_x, stress_y = tau_y, heat_flux = 50.0)` — the wind stress is
computed from real cached winds, but the **heat flux is a hardcoded 50 W/m²**, which is not a
formula at all. It should be replaced in this order:

1. **Bulk formula driven by the wind we already have (no credentials).** `inputs/wind_active.nc`
   is real, cached data. A standard bulk formulation — latent + sensible + longwave, e.g. the
   standard `q_air`-based or Clarke & Hess coefficients — turns that wind plus SST into a
   space- and time-varying heat flux. This is the highest-value change because it is physically
   meaningful, needs no new data, and removes a constant that currently dominates the surface
   energy budget. ClimaOcean already provides the closures:
   `TwoColorRadiation`, `CoefficientBasedFluxes`, `SimilarityTheoryFluxes`, `BulkTemperature`,
   `ChlorophyllOptics`, `TabulatedAlbedo`.
2. **ERA5 surface fields via CDS (needs a key).** `ERA5PrescribedAtmosphere` /
   `JRA55PrescribedAtmosphere` give real surface temperature, humidity, cloud and radiation.
   `cds.climate.copernicus.eu` **is reachable**; this route needs `~/.cdsapirc` with a CDS API
   key, and is a user action. Note this is the **CDS**, a different service from Copernicus Marine,
   so `cdsapi` + `~/.cdsapirc` authenticate it and it must not emit the `copernicusmarine login`
   reminder.

Wind data itself is already keyless and working: `fetch_open_surface_winds` tries
Open-Meteo → NOAA PSL → CDS, and all three hosts are reachable.

**ClimaOcean cannot supply CMIP/SSP projections.** Grepping ClimaOcean and NumericalEarth for
`SSP|CMIP|IPCC|rcp85|ssp245` returns nothing. What they do provide is *reanalysis-era* surface
forcing — `ERA5PrescribedAtmosphere`, `JRA55PrescribedAtmosphere`, radiation closures
(`TwoColorRadiation`, `ChlorophyllOptics`, `TabulatedAlbedo`) and flux closures
(`SimilarityTheoryFluxes`, `CoefficientBasedFluxes`, `BulkTemperature`) — plus coupled scaffolding
(`EarthSystemModel`, `OceanOnlyModel`, `AtmosphereOceanModel`, `PrescribedOcean`). None of it is a
projection. So the coefficient approach in `climate_scenarios.jl` has no upstream replacement;
replacing it with ClimaOcean would *remove* projection support rather than add it. Real
projections need a CMIP difference field, which nothing here reads.

### 4. GLORYS reader is not written yet

Access is solved — see [section 6](#6-glorys-support-is-required-and-access-is-solved), where an
anonymous ARCO-Zarr chunk read is demonstrated working. What is missing is the Julia side:
`fetch_copernicus_physics_subset` downloads but has **no reader** for `thetao`/`so`, so the driver
errors after a successful download. `NumericalEarth.Column(grid, dataset; …)` is a struct
describing a *point location*, not a dataset reader, so that call never worked and should be
replaced by `ClimaOcean.Metadatum` + `DataWrangling.FieldRegridding`.

### 5. WOA23 is downloaded and working, and the "network blocker" was never the blocker

**Status: WOA23 obtained, read, interpolated, and masked-value handling fixed. Rewritten
2026-09-27.**

The earlier claim that "network-dependent configs cannot be verified here" was **wrong**, and so was
the stated cause. The network is partially open. What actually blocked WOA23 was three code bugs
plus a NOAA service outage, described below.

**Reachability measured 2026-09-27**

| Host | :443 | Used for |
|---|---|---|
| `www.ncei.noaa.gov` | reachable | main site 200, but **THREDDS is 503** |
| `www.ncei.noaa.gov/data/...` (static) | reachable | **WOA23 NetCDF — works** |
| `coastwatch.pfeg.noaa.gov` | reachable | ETOPO 2022 bathymetry via ERDDAP |
| `cds.climate.copernicus.eu` | reachable | CDS ERA5 winds |
| `copernicusmarine.climate.copernicus.eu` | **blocked** | Copernicus Marine API |
| `s3.waw3-1.cloudferro.com` | reachable | Copernicus ARCO/STAC metadata |
| `thredds.ucar.edu`, `tds.marine.rutgers.edu`, `data.nodc.noaa.gov` | reachable but no `ncei/woa` collection | 404 |

NCEI's THREDDS returned 503 for every path while the main host served pages normally — a genuine
NOAA-side service outage, not a local block. The OPeNDAP mirrors the code listed all 404 on
`ncei/woa`.

**Four real bugs found and fixed in the WOA23 path** (`src/data/open_data.jl` and
`src/model/hydrodynamic_model.jl`)

1. **`fetch_woa23_hydrography` threw unconditionally.** It always forwarded `include_o2` to
   `fetch_open_woa_climatology`, which had no such keyword, so *every* call raised a `MethodError`
   before any network access. This is why the path had "never worked" — it was blamed on the
   network. `fetch_open_woa_climatology` now accepts `include_o2` and the download loop honours it.
2. **Wrong revision tag and wrong directory layout.** The URL builder used `A5B7`; no `A5B7`
   directory exists. The real layout is
   `data/oceans/woa/WOA23/DATA/<var>/netcdf/A5B4/{0.25,1.00}/woa23_A5B4_<v><MM>_0{4,1}.nc`.
   Oxygen is published **only at 1.00° and 5.00°**, under a different subtree with an `all_` prefix
   (`oxygen/netcdf/all/1.00/woa23_all_o00_01.nc`). Rewritten with the verified paths, keeping the
   OPeNDAP subset URLs as a secondary option since they transfer MBs instead of ~1.7 GB.
3. **The fetcher's reader pulled the entire global array.** `ds[vname][:, :, :, 1]` on a
   1440 × 720 × 102 field is ~850 MB per variable as Float64, and `coalesce.` doubles it
   transiently. Now reads only the requested lon/lat window (coordinate axes are read in full
   first to resolve the indices), keeping the whole depth axis so no profile is truncated. A
   shelf-scale domain is a few percent of the global field.
4. **`read_woa_variable` — the *consumer-side* reader — crashed on real files.** It assumed a 3-D
   variable and did `permutedims(Array(var), (3, 2, 1))`, so the real 4-D
   `(lon, lat, depth, time=1)` field raised `ArgumentError: no valid permutation of dimensions`.
   It also only handled the old `-1e30` fill convention, not the `Missing` mask the current
   release uses. This surfaced the moment a real run reached Segment 3; the earlier "regridding
   validated against a synthetic WOA-format NetCDF" had used 3-D synthetic data, so the gap was
   invisible. Now accepts 3-D or 4-D, collapses a singleton time axis, maps `Missing` → `NaN`,
   and errors clearly if multiple time levels are present.

   Note there are **two independent WOA23 readers** — `make_woa_interpolator` in
   `src/data/open_data.jl` (download side) and `read_woa_variable` in
   `src/model/hydrodynamic_model.jl` (regrid side). They had drifted apart.

**Data obtained and verified** — cached in `inputs/` at the exact paths the fetcher checks, so
subsequent runs are cache hits:

| File | Size | Grid |
|---|---|---|
| `woa23_temperature_00_0.25deg.nc` | 884.0 MB | 1440 × 720 × 102, 0.25° |
| `woa23_salinity_00_0.25deg.nc` | 782.1 MB | 1440 × 720 × 102, 0.25° |
| `woa23_oxygen_00_0.25deg.nc` | 71.4 MB | 360 × 180 × 102, 1.00° |

Land is masked as `Missing` (not NaN) — 50.4% of cells — which is why the reader uses `coalesce`.
Interpolated values at (−62, 44) are physically correct for the Scotian Shelf: T 9.38 → 7.14 →
4.91 °C through the upper 50 m (the cold intermediate layer), S 31.48 → 34.00, O₂ 206 → 40.9
µmol/kg at −200 m, and genuine horizontal structure at −100 m (T = 6.78 at −62° vs 8.82 at −55°).

**Masked-value semantics — decided and implemented.** Masked WOA cells are now **`NaN`**, never
`0.0`, and the reason for a gap is respected when repairing it. This was the S = 0 hazard: a `0.0`
fill is not "absent", it is a physical claim (0 K, 0 PSU) that corrupts the equation of state.

`make_woa_interpolator` now reads masked cells as `NaN` and repairs in two stages, because the
two reasons a cell is masked want different fixes:

1. **Gaps inside a water column** (WOA masks the volume below the seabed) are filled by linear
   interpolation between the bracketing valid depths, carrying the nearest valid value down. Below
   the seabed the bottom water is a far better estimate than any global mean, and this is the case
   that actually occurs in the Scotian Shelf domain.
2. **Columns masked over their whole depth range** (land) fall back to the mean of valid values in
   the requested window, and the counts are printed so they are never silently invisible. These
   should not be sampled by the model, whose water column is built from bathymetry.

Verified by targeted test (`work/test_nan_fill.jl`): 25 344 queries across the domain returned
**0 NaN and 0 exact-zero** values; salinity 29.23 – 36.06 (the low end is genuine Gulf of St.
Lawrence dilution) and temperature −1.25 – 20.52 °C. Below the seabed the repair holds bottom
water rather than a mean — at (−62, 44) it gives T = 7.47 °C, S = 34.25 at −200 through −500 m —
and real horizontal structure survives: at −300 m, S ranges 34.01 → 34.96 and T 5.95 → 7.25 °C
across −70° to −57°.

A further refinement is available but not taken: stage 2 could fill fully-masked columns from the
nearest *horizontally adjacent* wet column rather than a window mean, which matters only if the
domain contains in-ocean cells that WOA masks. See task 8 in the list above.

### 5a. Bathymetry: the ERDDAP dataset ID was wrong

Found while verifying bathymetry reachability. `fetch_etopo2022_bathymetry` and
`fetch_open_bathymetry` both requested **`nceiEtopo2022`, which does not exist on the live ERDDAP
catalogue (404)**. Verified working ids, with the elevation variable each one uses:

| ERDDAP dataset id | elevation variable | status |
|---|---|---|
| `ETOPO_2022_v1_15s` | `z` | 200 — 15 arc-second, the intended product |
| `etopo180` | `altitude` | 200 — 1 arc-minute fallback |
| `srtm30plus` | `elevation` | 200 — added as a third fallback |

`ETOPO_2022_v1_15s` is now the default and primary, with an id → variable-name map so a caller
cannot pair an id with the wrong variable. A live subset request was verified to return a valid
NetCDF. `inputs/bathymetry_active.nc` (8.0 MB, 1449 × 725) was already cached, which is why this
never surfaced in a run — the cache check fires before the download.

### 6. GLORYS support is required, and access is solved

**Decision: GLORYS support is required. The hard part — getting the bytes — is now solved and
proven anonymously, with no Copernicus credentials.** What remains is the `thetao`/`so` reader in
Julia and wiring it to the regrid path.

**The `copernicusmarine` API host is blocked; the data host is not.** `copernicusmarine.climate.copernicus.eu`
fails at TCP, but the hosts that actually serve Copernicus Marine data are all reachable:
`s3.waw3-1.cloudferro.com` (200), `marine.copernicus.eu` (200),
`stac.marine.copernicus.eu/clients-config-v1` (200), and the STAC catalogue on the CloudFerro S3
(200). The installed `copernicusmarine` 2.4.1 client still cannot use them — its catalogue load
reports `0/2` and resolves every dataset to "not found", it ships no static catalogue to fall back
on, and `COPERNICUSMARINE_SERVICE_URL` does not change that. So the client is bypassed entirely.

**Access route, verified end to end (2026-09-27)**

1. `…/metadata/GLOBAL_MULTIYEAR_PHY_001_030/product.stac.json` → 200, and its `item` links name
   the concrete datasets, e.g. `cmems_mod_glo_phy_my_0.083deg_P1M-m_202311` (monthly) and
   `…_P1D-m_202311` (daily).
2. The dataset STAC exposes the ARCO assets directly:
   - `timeChunked` — `https://s3.waw3-1.cloudferro.com/mdl-arco-time-025/arco/<ID>/<DS>/timeChunked.zarr`
   - `downsampled4`, `geoChunked`, and a `native` bucket
3. **`.zmetadata` is anonymously readable (200)**, as are `.zgroup`, and every per-variable
   `.zarray` / `.zattrs`. Variables present: `thetao`, `so`, `uo`, `vo`, `usi`, `vsi`, `zos`,
   `bottomT`, `mlotst`, `siconc`, `sithick`, plus `time`, `latitude`, `longitude`, `elevation`.
   So the reanalysis supplies hydrography **and** currents/SSH, which is what an open boundary
   condition would need.
4. `thetao` / `so`: shape `(401, 50, 2041, 4320)` = (time, elevation, lat, lon), chunks
   `(1, 1, 512, 2048)`, `int16`, Blosc/LZ4, `fill_value = -32767`, with
   `scale_factor = 0.0007324442267417908`, `add_offset = 21.0`. Decoded °C = `raw * scale + offset`;
   the fill must be masked on the **raw** int16 value, because the fill decodes to −3.0 °C, which
   is a plausible-looking temperature rather than an obvious sentinel.
   `elevation` is ordered deepest-first, so the surface is the **last** index, not the first.
5. **A real chunk read was performed, anonymously** (`work/probe_glorys_chunks.py`): `.zmetadata`
   in 0.5 s, domain window resolved to lon idx 1308–1524 / lat idx 1440–1542, chunk
   `thetao/0.2.0.0` fetched in 0.6 s (2.1 MB compressed) and sliced to the exact (96, 216) domain
   window.

**Why not the `xarray`/`zarr` route.** `xr.open_zarr` fails against this store: zarr-python 3.x
probes the Zarr-v3 `zarr.json` key, and the bucket answers **403 Forbidden rather than 404**, so
store discovery raises `ClientResponseError` before reading anything. Reading `.zmetadata` and
fetching chunk keys directly with `requests` + `numcodecs` works and is fully deterministic.
`requests` and `aiohttp` were missing from `.venv` and are now installed (PyPI is reachable).
Available in `.venv`: `zarr 3.1.6`, `xarray 2026.7.0`, `fsspec 2026.9.0`, `numcodecs`,
`netCDF4 1.7.4`, `copernicusmarine 2.4.1`. `s3fs` is absent and not needed.

**Implementation plan**

1. Python helper in the venv: given lon/lat window, depth range, time window, and variable list,
   resolve chunk keys from `.zmetadata`, fetch, Blosc-decode, scale, and write a regional NetCDF
   into `inputs/` in the same layout the Julia readers already expect.
2. Julia reader for `thetao`/`so`, replacing the absent one in `fetch_copernicus_physics_subset`.
   Reuse the conventions now proven on WOA23: 4-D-with-singleton-time tolerance, `Missing`/fill →
   `NaN`, ascending axes, mask-before-decode for packed integers.
3. Regrid via `ClimaOcean.Metadatum` + `DataWrangling.FieldRegridding`. Note
   `NumericalEarth.Column(grid, dataset; …)` is a struct describing a *point location*, not a
   dataset reader, so that call never worked and should not be used.
4. Remap the two deleted IDs. `GLOBAL_MULTIYEAR_PHY_001_033` and `GLOBAL_REANALYSIS_PHY_001_031`
   no longer exist; `GLOBAL_MULTIYEAR_PHY_001_030` and `GLOBAL_ANALYSISFORECAST_PHY_001_024` both
   resolve. Credentials and per-dataset terms acceptance become unnecessary if the anonymous
   ARCO/Zarr route is used, which sidesteps the interactive `copernicusmarine login` entirely.
5. Cost: GLORYS12 is 1/12° and far larger than WOA23's 0.25°. Chunks are 2 MB each and
   `(1, 1, 512, 2048)`, so a regional monthly extract is tens of chunks — cheap. The full
   401-month global series is not, and should not be cached wholesale.

**Client setup that is already in place.** `copernicusmarine` **2.4.1** is installed in the project
`.venv`, and the Julia side resolves the CLI explicitly rather than relying on `PATH`:

- `project_python()` — `PARTICLETRACKING_PYTHON` → `PYTHON_EXECUTABLE` → `<root>/.venv` →
  `VIRTUAL_ENV` → `python` on `PATH`.
- `copernicusmarine_executable()` — `COPERNICUSMARINE_EXE` → the script directory beside the
  resolved interpreter → `PATH`. Verified to find `.venv\Scripts\copernicusmarine.exe`
  **with `.venv` stripped from `PATH`**.
- `copernicus_login_reminder()` — appended to every Copernicus download failure.

Still required, and not done:

1. `.\.venv\Scripts\copernicusmarine.exe login` — interactive; writes
   `~\.copernicusmarine\credentials.toml`. Check with `login --check-credentials-valid`.
2. **Per-dataset terms-and-conditions acceptance.** Valid credentials are necessary but not
   sufficient; each dataset must be accepted in the portal or the API returns nothing. Note that
   two of the four IDs this repo referenced no longer exist, so the T&C work has to be redone
   against whatever supersedes them.

`netCDF4` 1.7.4 was added to the venv for inspecting downloads. Note the CDS *wind* path
(`fetch_copernicus_surface_winds`) authenticates through `cdsapi` + `~/.cdsapirc`, a **different**
service from Copernicus Marine, so it deliberately does **not** emit the `copernicusmarine login`
reminder.

### 7. Network-dependent configs

Segments 1 and 2 of `opendata.toml` complete against real data: bathymetry read
(`inputs/bathymetry_active.nc`, 1449 × 725) and wind read (`inputs/wind_active.nc`, 0.9 MB),
grid and immersed boundary built (145 × 105 × 15, tanh-stretched, −3000 to 0 m). Segment 3 then
failed on the WOA23 read — bug 4 in [section 5](#5-woa23-is-downloaded-and-working-and-the-network-blocker-was-never-the-blocker),
now fixed and verified against the real files. **The run has not been repeated past that point**, so
segments 4–8 remain unverified against live data. `snowcrab.toml` is likewise unverified. Both are
`use_gpu = false`; the GPU blockers are separate and documented in
[section 2](#2-gpu-root-cause-found-field_dependencies-on-a-continuousforcing).

All non-Copernicus data hosts are reachable and their paths are now verified, so no further
network investigation should be needed for T/S/O2/bathymetry/wind.

### 8. Copernicus Marine client setup — superseded

`copernicusmarine` 2.4.1 is in `.venv` and the Julia side resolves it explicitly
(`project_python()`, `copernicusmarine_executable()`, `copernicus_login_reminder()`), but the
client cannot reach a usable catalogue and **is no longer on the critical path**: the anonymous
ARCO/Zarr route in [section 6](#6-glorys-support-is-required-and-access-is-solved) bypasses it,
and with it the interactive `copernicusmarine login` and the per-dataset terms-and-conditions
acceptance. Kept only so the resolution helpers are not mistaken for dead code.

`netCDF4` 1.7.4 is in the venv for inspecting downloads. The CDS *wind* path
(`fetch_copernicus_surface_winds`) authenticates through `cdsapi` + `~/.cdsapirc`, a **different**
service from Copernicus Marine, so it deliberately does **not** emit the `copernicusmarine login`
reminder.

### 9. Housekeeping

- `work/baseline.jld2` is **42 GB** and fully allocated (not sparse) on a 194 GB volume. Every other
  file in `work/` is 77 KB – 14 MB, so it looks like a runaway write, but it was left in place
  pending a decision. Delete it if it is not a deliberate reference.
- `inputs/` now holds ~1.7 GB of WOA23 NetCDF. Cached on purpose, but it is untracked — decide
  whether to `.gitignore` it or document how to re-fetch it.
- `outputs/` and `work/` are untracked and growing; consider a `.gitignore` entry.
- Do not edit these files with PowerShell `Set-Content` on a non-ASCII file. Reading with
  `Get-Content` and writing back re-encoded UTF-8 as Latin-1 and corrupted `ΔT` into `Î”T` and `°C`
  into `Â°C`, breaking the parse. It cost a full restore from git. Use the editor tooling, or
  `Get-Content -Encoding UTF8` / `Set-Content -Encoding UTF8` as a matched pair.

**Verified data sources, 2026-09-27 — what actually works from this machine**

| Data | Working source | Note |
|---|---|---|
| WOA23 T/S (0.25°) | `https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/{temperature,salinity}/netcdf/A5B4/0.25/woa23_A5B4_{t,s}<MM>_04.nc` | 884 / 782 MB, global |
| WOA23 O2 (1.00°) | `https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/1.00/woa23_all_o<MM>_01.nc` | 71 MB; **no 0.25° oxygen exists** |
| ETOPO 2022 | `https://coastwatch.pfeg.noaa.gov/erddap/griddap/ETOPO_2022_v1_15s.nc?z[…][…]` | variable is `z`, not `altitude` |
| Surface winds | Open-Meteo archive API → NOAA PSL → CDS | all three hosts reachable |
| CDS ERA5 | `cds.climate.copernicus.eu` | reachable; needs `~/.cdsapirc` |
| GLORYS STAC metadata | `s3.waw3-1.cloudferro.com/mdl-metadata/metadata/<id>/product.stac.json` | reachable, but **not** `copernicusmarine.climate.copernicus.eu` |
| NCEI THREDDS (OPeNDAP) | — | **503 outage**; all four listed mirrors 404 on `ncei/woa` |

Refetch script: `work/get_woa23.ps1` (plain `WebClient`; `[Net.Http.HttpClient]` is unavailable in
Windows PowerShell 5.1 without an explicit `Add-Type`).

Added to `.venv` on 2026-09-27: `requests 2.34.2` and `aiohttp 3.14.3`. fsspec's HTTP backend
requires both, and without them any fsspec-backed store over HTTPS fails with
`ImportError: HTTPFileSystem requires "requests" and "aiohttp" to be installed`. PyPI is reachable
from this machine, which is what disproved the earlier "no network" claim.

`copernicusmarine` **2.4.1 is installed** in the project `.venv`, and the Julia side now resolves
the CLI explicitly rather than relying on `PATH`:

- `project_python()` — `PARTICLETRACKING_PYTHON` → `PYTHON_EXECUTABLE` → `<root>/.venv` →
  `VIRTUAL_ENV` → `python` on `PATH`.
- `copernicusmarine_executable()` — `COPERNICUSMARINE_EXE` → the script directory beside the
  resolved interpreter → `PATH`. Verified to find
  `.venv\Scripts\copernicusmarine.exe` **with `.venv` stripped from `PATH`**.
- `copernicus_login_reminder()` — appended to every Copernicus download failure, quoting
  `.\.venv\Scripts\copernicusmarine.exe login`.

Still required, and not done:

1. `.\.venv\Scripts\copernicusmarine.exe login` — interactive; writes
   `~\.copernicusmarine\credentials.toml`. Check with `login --check-credentials-valid`.
2. **Per-dataset terms-and-conditions acceptance.** Valid credentials are necessary but not
   sufficient; each dataset must be accepted in the portal or the API returns nothing. Unverified
   which of `GLOBAL_MULTIYEAR_PHY_001_030/033`, `GLOBAL_REANALYSIS_PHY_001_031`,
   `GLOBAL_ANALYSISFORECAST_PHY_001_024` have been cleared.

`netCDF4` 1.7.4 was added to the venv for inspecting downloads. Note the CDS *wind* path
(`fetch_copernicus_surface_winds`) authenticates through `cdsapi` + `~/.cdsapirc`, a **different**
service from Copernicus Marine, so it deliberately does **not** emit the `copernicusmarine login`
reminder.

---

## Do not regress

- **One place reads the TOML.** `configuration_to_options` in `src/config/configuration.jl`. Do not
  reintroduce per-key `get(cfg, …)` reads in the driver; a key added to the loader is otherwise
  silently ignored. This exact failure mode hid 8 keys for a long time.
- **A config value that never reaches the call site is a silent lie in `resolved_config.toml`.**
  The biology CV parameters are the current example: loaded, defaulted to `0.25`, written to the
  provenance file — and never passed to `track_larval_cohort`, so every run used `0.0`. Adding a
  key to the loader and to `options_to_configuration` is **not** enough; trace the value all the
  way to the function that consumes it. When adding a knob, grep for the consumer in the same
  change.
- **A `cv_*` knob must control dispersion only, never the mean.** The previous Normal-then-clamp
  form broke this twice over: `clamp(Normal(hsi, cv·hsi), 0, 1)` pushed the expected settlement rate
  *below* `hsi` (0.90 instead of 1.00 at the optimum), and a `max(0, …)` on the mortality rate
  piled individuals onto a zero rate. Use `draw_lognormal_mean` (mean preserving, via the `−σ²/2`
  term) for positive rate-like quantities and `draw_beta_index` for bounded ones. If a dispersion
  knob ever changes an expected outcome, that is a bug.
- **Never clamp a physical draw to put it in range; choose a distribution with the right support.**
  Clamping puts a point mass on the boundary that is indistinguishable from real structure, and it
  is not mean preserving. A clamped draw gave 2.4 % of larvae a two-hour Zoea I stage.
- **Check the `LateralBoundaryRelaxation` construction against the new forcing path.** The sponge is
  now a `Relaxation`, so `sponge_relaxation_u` / `_v` are the reference formula only, not the code
  the model runs. Changing `compute_sponge_gamma` or the `ExponentialInflow` targets affects the
  live path through the `mask` and `target` arguments, not through those entry points.
- **Verify a remote dataset id and its variable name against the live catalogue before shipping a
  URL.** Three separate wrong assumptions were sitting in the data layer: `nceiEtopo2022` (404;
  the real id is `ETOPO_2022_v1_15s` and its variable is `z`, not `altitude`), WOA revision `A5B7`
  (no such directory; the real one is `A5B4`), and the WOA23 NetCDF layout itself. A one-line
  `Invoke-WebRequest` against the `.das`/`.dds` endpoint settles it in seconds.
- **Test data readers against real files, not just synthetic ones shaped like them.** The WOA23
  consumer reader passed against a synthetic 3-D WOA-format NetCDF and then failed on the real
  4-D file (`time` = 1). Synthetic fixtures have to copy the real rank and masking convention —
  current WOA23 masks land as `Missing`, not the older `-1e30` fill.
- A **cache hit hides a broken download path.** `inputs/bathymetry_active.nc` and the WOA23 files
  are checked before any request, so a wrong URL never fires while a cache is warm. Clear the
  cache (or point at a scratch `inputs/`) when testing acquisition specifically.
- **Do not assume a "no network" diagnosis.** Verify per host before concluding it. Every failure
  in this project attributed to the network turned out to be something else: an `A5B7` path that
  does not exist, a 404 dataset id, a 3-D assumption on a 4-D file, a THREDDS outage that left the
  static data server working fine, and missing Python packages that `pip` installs in seconds.
- **A 403 is not a 404.** A cloud bucket that answers 403 for absent keys breaks
  store/metadata discovery that would otherwise fall through to the next candidate. When a remote
  directory listing "does not exist", check which status it actually returns.
- **`axes` on a `Field` is the interior.** Iterate `axes(field.data, k)` whenever halos matter —
  `∂z` reads one column past the interior. Halo cells of a freshly built model hold
  *uninitialised* memory: reading it and feeding it to a polynomial can overflow to Inf/NaN and
  take the process down. `fill_seawater_density!` now initialises its buffer to the EOS reference
  density and evaluates only the interior (range derived from `axes` and `grid.Nx/Ny/Nz`).
- **A standalone `HydrostaticFreeSurfaceModel` has no `boundary_conditions` field.** Its fields are
  `architecture, grid, clock, advection, buoyancy, coriolis, free_surface, forcing, closure,
  particles, biogeochemistry, velocities, transport_velocities, tracers, pressure, closure_fields,
  timestepper, auxiliary_fields, vertical_coordinate, boundary_transport`. BCs live on the `Field`s
  and are fixed at construction, so post-hoc BC surgery raises `FieldError` and is not supported.
  Use `EarthSystemModel` + `PrescribedOcean` for prescribed surface state that evolves in time.
- **`indices(field, k)` returns a `Colon` in Oceananigans 0.111** — it is not a usable index range.
  Derive interior ranges from `axes(field.data, k)` and the grid's `Nx`/`Ny`/`Nz`.
- **Grid node accessors in 0.111 are `λᶠᵃᵃ`, `φᵃᶠᵃ`, `zᶜᶜᶜ`** (not `λᵃᶠᵃ`). The codebase already uses
  the correct ones; a wrong one fails at `getproperty`.
- **Forcing closures must stay `isbits`, and must not take a `Symbol`.** `isbitstype(Symbol)` is
  `false`, so a `Symbol`-dispatched forcing method is unlowerable into a GPU kernel. Use concrete
  structs with one method per component (`sponge_relaxation_u` / `sponge_relaxation_v`) rather
  than a runtime selector.
- **Declare `field_dependencies`** whenever a forcing needs model state, or Oceananigans will call
  it with 4 arguments instead of 6. **But do not use it on the GPU path** — any
  `ContinuousForcing` carrying `field_dependencies` produces invalid LLVM IR in
  `gpu_compute_hydrostatic_free_surface_Gu!` on this stack, whether the dependency is a velocity
  or a tracer. Use `Relaxation(; rate, mask, target)` for regional nudging; it reads the relaxed
  field through `r.relaxed` and needs no dependencies.
- **A `ContinuousForcing` is only the right tool for stateless-in-space tendencies.** Anything that
  needs to read the field it is forcing belongs in `Relaxation`. `field_dependencies` on a
  momentum forcing is also circular — the sponge reads the very `u`/`v` being updated.
- **A passing `configs/default.toml --gpu` run does not mean the GPU path works.** That config
  disables the sponge and the tides, so it never builds a forcing with field dependencies. When
  checking GPU health, reproduce the `snowcrab` forcing combination, not the `default` one.
- **Keep the honest-failure policy.** A run must never claim a dataset-forced state it did not
  read; raise instead, or make the synthetic choice explicit and recorded.
- **`_grid.jld2` sidecars** are the only supported way to read grid geometry back; deserialising
  `serialized/grid` yields an unusable `JLD2.ReconstructedStatic`.
- Geography comes only from TOML `[domain]` / `[boundaries]`; no hard-coded domain constants.
