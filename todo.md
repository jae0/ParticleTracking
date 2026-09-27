# ParticleTracking — TODO

Continuation notes. Supersedes the earlier `REFACTOR_PLAN.md` (deleted 2026-09-26).

State at last update: **661/661 tests pass**. `configs/default.toml` runs green end to end on
CPU. `configs/snowcrab.toml` and `configs/opendata.toml` request `hydrography_source = "woa23"`
and therefore fail on this machine for lack of network access (see [Open items](#open-items)).

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
  the same `InvalidIRError` persists after the fix.

### Source layout

The package was reorganised into `src/{utils,data,model,config,biology,analysis,output}/`. Files
were **moved**, not rewritten, so all working code is preserved byte-for-byte; the package is a
single flat namespace, so these remain plain `include`s in dependency order rather than submodules.

---

## Open items

### 1. Larval biology: growth and mortality are deterministic

**Settlement is stochastic. Growth and mortality are not.** This was investigated deliberately and
the asymmetry is real, not an oversight in one function.

| Process | Stochastic? | Evidence |
|---|---|---|
| Initial positions / depths | Yes | `rand` in `initialize_larval_particles` (`src/larval_behavior.jl:799-903`) |
| Turbulent transport | Yes | `randn` for the three diffusion increments (`:652-654`) |
| Settlement | Yes | `rand(rng) <= hsi` (`:527`) |
| **Growth / molting** | **No** | `update_larval_stage` is pure threshold logic on cumulative degree-days (`:409-427`) |
| **Mortality** | **No** | survival is an expected *fraction* `exp(-μ·dt)`; death fires when the fraction crosses `min_survival_prob` (`:1288-1289`) |

Consequences:

- Particles sharing a thermal history have **identical** degree-day accumulation and therefore molt
  in lockstep, so stage durations have no spread. Real cohorts show substantial individual
  variation in developmental rate; Sainte-Marie & Sainte-Marie (1999) report stage durations with
  real spread for *Chionoecetes opilio*.
- Mortality is applied as a mean, so a cohort's deaths are spread smoothly rather than drawn
  independently. Every particle in a cohort with similar exposure dies at the same step.
  `traj_surv` is reported as a fraction, so downstream summaries are expectations, not outcomes.

Decisions still needed before changing this (it alters results and needs calibration, not just code):

1. Whether to draw mortality per particle per step (`rand() < μ·dt`) and keep the fraction only for
   reporting, or keep the fraction and accept a mean-field cohort.
2. Where individual variation in the degree-day thresholds should enter — a per-particle draw at
   cohort initialisation, or a distribution over the threshold.
3. What published values to centre any new distributions on.

Also worth exposing: `settlement_stochastic` is **not** in any TOML config and is never passed by
the driver; it works only because `track_larval_cohort` defaults it to `true`. Add a
`[biology] settlement_stochastic` key and pass it explicitly so the choice is recorded in
`resolved_config.toml`.

### 2. GPU path is still broken — and it is NOT a Julia-version bug

`configs/default.toml --gpu` now **completes** (the `EXCEPTION_ACCESS_VIOLATION` was the
uninitialised-halo read in `fill_seawater_density!`, now fixed — see below). `snowcrab.toml --gpu`
still fails at kernel-compile time.

The environment was re-tested after a library upgrade (Oceananigans 0.111.0, ClimaOcean 0.10.0,
NumericalEarth 0.8.1, CUDA 6.0.0, Manifest refreshed 2026-09-27 06:02). **The failure is
unchanged**, so the lib upgrade did not fix it.

- **Julia is still 1.13.1** (`juliaup status` → `release 1.13.1+0.x64.w64.mingw32`). The only installed
  channel is 1.13.1, and Oceananigans prints an explicit warning that it *"is currently tested on
  Julia v1.12"*. **Installing a 1.12 channel and retrying is the single highest-value next test.**
- Remaining error: `InvalidIRError` in `gpu_compute_hydrostatic_free_surface_Gu!`, forcing tuple
  `MultipleForcings{2, Tuple{BarotropicPotentialForcing, ContinuousForcing{LateralBoundaryRelaxation,
  TidalBodyForcingU}}}`.

One real GPU bug was found and fixed while investigating (see `LateralBoundaryRelaxation` below),
but it was **not** the cause — the same `InvalidIRError` persists after the fix. Do not set
`use_gpu = true` in a config you intend to run.

### 3. Atmospheric forcing is still a placeholder

`atmo_craft = (stress_x = tau_x, stress_y = tau_y, heat_flux = 50.0)` — a constant, not a bulk
formula. `ERA5PrescribedAtmosphere` needs a download that has never worked here. Wind *is* keyless
via `fetch_open_surface_winds` (Open-Meteo → NOAA PSL → CDS), so the missing piece is the heat-flux
/ stress-from-wind coupling, not the wind data.

**ClimaOcean cannot supply CMIP/SSP projections.** Grepping ClimaOcean and NumericalEarth for
`SSP|CMIP|IPCC|rcp85|ssp245` returns nothing. What they do provide is *reanalysis-era* surface
forcing — `ERA5PrescribedAtmosphere`, `JRA55PrescribedAtmosphere`, radiation closures
(`TwoColorRadiation`, `ChlorophyllOptics`, `TabulatedAlbedo`) and flux closures
(`SimilarityTheoryFluxes`, `CoefficientBasedFluxes`, `BulkTemperature`) — plus coupled scaffolding
(`EarthSystemModel`, `OceanOnlyModel`, `AtmosphereOceanModel`, `PrescribedOcean`). None of it is a
projection. So the coefficient approach in `climate_scenarios.jl` has no upstream replacement;
replacing it with ClimaOcean would *remove* projection support rather than add it. Real
projections need a CMIP difference field, which nothing here reads.

### 4. GLORYS is not read at all

`NumericalEarth.Column(grid, dataset; …)` is a struct describing a point location, not a dataset
reader, so that call never worked. The supported path is `ClimaOcean.Metadatum` +
`DataWrangling.FieldRegridding`. Given WOA23 now covers the keyless case, decide whether GLORYS
support is still wanted at all before investing in it.

Separately, `fetch_copernicus_physics_subset` downloads fine but has **no reader** for the
`thetao`/`so` variables, so the driver errors after a successful download. That gap is still open
and blocks any end-to-end GLORYS run.

### 5. Network-dependent configs cannot be verified here

`snowcrab.toml` and `opendata.toml` both set `hydrography_source = "woa23"` and fail on this machine
because `www.ncei.noaa.gov` is unreachable. This is the intended honest-failure behaviour, and it
proves the wiring works. What is **not** yet verified is the full download → regrid → run path
against live data; the regridding itself was validated against a synthetic WOA-format NetCDF
(exact node reproduction, correct lon wrapping and depth reversal, NaN masking, out-of-range
clamping). Run both configs on a networked machine before trusting them.

### 6. Copernicus Marine setup (ready, not yet authenticated)

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

### 7. Housekeeping

- `work/baseline.jld2` is **42 GB** and fully allocated (not sparse) on a 194 GB volume. Every other
  file in `work/` is 77 KB – 14 MB, so it looks like a runaway write, but it was left in place
  pending a decision. Delete it if it is not a deliberate reference.
- `outputs/` and `work/` are untracked and growing; consider a `.gitignore` entry.

---

## Do not regress

- **One place reads the TOML.** `configuration_to_options` in `src/config/configuration.jl`. Do not
  reintroduce per-key `get(cfg, …)` reads in the driver; a key added to the loader is otherwise
  silently ignored. This exact failure mode hid 8 keys for a long time.
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
  it with 4 arguments instead of 6.
- **Keep the honest-failure policy.** A run must never claim a dataset-forced state it did not
  read; raise instead, or make the synthetic choice explicit and recorded.
- **`_grid.jld2` sidecars** are the only supported way to read grid geometry back; deserialising
  `serialized/grid` yields an unusable `JLD2.ReconstructedStatic`.
- Geography comes only from TOML `[domain]` / `[boundaries]`; no hard-coded domain constants.
