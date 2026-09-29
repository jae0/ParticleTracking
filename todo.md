# ParticleTracking — TODO

State at last update: **2026-09-29**. Supersedes the earlier `REFACTOR_PLAN.md` (deleted
2026-09-26).

**Where things stand.** The TOML file is now the single source of truth: struct field names
match the TOML keys, the driver reads nothing by hand, and there is one interpreter of the
config (`configuration_to_options`). All three `--all` pipelines complete end to end against
real data. Dead code from the refactor has been removed (ArgParse, the `data_mode` switch, the
`fetch_era5_atmospheric_forcing` fabricator, two unused wind wrappers).

**The single most important caveat.** The full-length hydrodynamics diverge — `max(|u|)` grows
to thousands of m/s while the median stays physical. Any result read off that field is
indicative, not quantitative. This is item 1 and it gates the rest.

---

## Open items

| # | Item | Status | Blocker |
|---|---|---|---|
| 1 | **Velocity divergence in the full-length run** | **open, highest priority.** `max(|u|)` climbs from 0 to ~5600 m/s over 730 days while the median stays at ~0.03 m/s. Not the CFL formula and not the time-step floor; both were ruled out with evidence. | needs a velocity time series; the probe harness has failed three times and produced nothing yet |
| 2 | **GLORYS download works; the credentials file format does not** | **route proven end to end 2026-09-29.** Client v2.4.1 authenticates against `data.marine.copernicus.eu` and a real 33 KB `thetao` NetCDF was retrieved and read. The blocker is that `~/.copernicusmarine/credentials` is **base64-encoded**, so the official client cannot parse it and reports "requires a username and password" | decode the file in `fetch_copernicus_physics_subset`, or rewrite it as plain INI |
| 3 | **HYCOM data path refused** | **code written, server side blocking.** `fetch_hycom_boundary` is implemented and exported; the catalogue, variables and grid are confirmed live, but every constrained data request returns HTTP 400. **Superseded by item 2** — GLORYS needs no further infrastructure. | only worth revisiting if account-free access becomes a requirement |
| 4 | **Benthic temperature dataset** | design settled, not built | awaiting the file from the user — see [Benthic temperature ingest](#benthic-temperature-ingest) |
| 5 | **Tidal amplitudes are too small by ~5–10x** | **open.** Bay of Fundy M2 head amplitude computes to 0.086 m/s; the published value is 0.5–1.0 m/s. The forcing is now dimensionally correct (`U_k·ω_k·cos(ω_k t)`, a tendency in m/s²) but the amplitude is not yet understood. | not diagnosed. Deliberately **not** rewired onto the sponge until verified — the earlier ~7000x error came from exactly this kind of rewiring |
| 6 | **Stochasticity numbers are placeholders** | functional form verified (lognormal / lognormal / Beta, all mean preserving); the `0.25` values are not calibrated | needs a cohort time series or field observations |
| 7 | **Copernicus dataset IDs to remap** | `GLOBAL_MULTIYEAR_PHY_001_033` and `GLOBAL_REANALYSIS_PHY_001_031` no longer exist | needs a Copernicus account to search the current catalogue |
| 8 | **Mixed-layer depth has no figure** | the model computes MLD; nothing plots it | low priority |
| 9 | **`min_dt_seconds` reachability** | hardcoded `2.0` removed; the driver's `opts` path still not confirmed end to end | needs a print at the call site |

**Settled decisions, not to reopen.** GPU work proceeds with any reasonable workaround (the
`Relaxation` migration is numerically exact, so nothing changes). Larval drift runs the full
length of the hydrodynamic model and stops early only when every larva is dead or settled.
Masked WOA cells become `NaN` and are then repaired. LTO is skipped (`allow_analytical_fallback`
is opt-in, because a silent fallback once produced recruitment numbers from a synthetic
current). Geography comes only from TOML; no hard-coded domain constants.

---

## 1. Velocity divergence — the main open problem

**Symptom.** Over a 730-day run the saved hydrodynamics diverges monotonically. `|u|` over wet
cells, in m/s:

| day | median | p90 | p99 | max |
|---|---|---|---|---|
| 73 | 0.032 | 7.1 | 50.8 | 562 |
| 219 | 0.035 | 22.2 | 160.9 | 1741 |
| 365 | 0.029 | 39.5 | 265.9 | 2979 |
| 511 | 0.023 | 59.5 | 382.4 | 4264 |
| 657 | 0.023 | 81.2 | 511.1 | **5600** |

The **median is physical**. The tail is four orders of magnitude above any real Scotian Shelf
current. The `divergence_velocity_limit` watchdog only fires at 20 m/s and merely warns, so the
run completes and reports numbers nobody should trust.

**Two hypotheses raised, both wrong — recorded so they are not re-proposed.**

1. *"The CFL formula caused it."* **Half right, and not causal.** `compute_advective_cfl`
   computed the vertical Courant number from the **mean** spacing (`Lz/Nz`) rather than the
   minimum, underestimating the constraint by ~20x on a stretched grid. But that function only
   feeds the **progress log line**; the step actually used is chosen by Oceananigans' own
   `TimeStepWizard`, which is correct. The fix made the diagnostic honest — the log read
   `CFL: 0.428` when the true vertical Courant was ~10, which is precisely why the failure was
   undiagnosable from the output — but it did not stop the instability.
2. *"`min_dt_seconds = 2.0` clamped the step."* **Wrong.** See item 9; the effective floor was
   0.1 s, so nothing was clamped.

**Still untested candidates:** the surface layer genuinely needing a ~1.1 s step with the initial
transient outpacing its dissipation; the sponge / open-boundary configuration; the surface
heat-flux term; the initial stratification; the resolved pycnocline. A climb from 0 to 20 m/s
over hours is consistent with several of these, and the watchdog cannot distinguish them.

**What is needed.** A velocity time series. The probe harness has failed three times (soft-scope
capture, `step!` not being in Oceananigans, and iterating a `Simulation` directly yielding
nothing). Nothing about stability should be assumed until it produces one.

**It may be fixed by item 4.** WOA23's annual-mean climatology smears the cold intermediate
layer, and an under- or over-resolved pycnocline is a classic source of spurious baroclinic
instability. A properly resolved CIL could be both more accurate and more stable — a
hypothesis to test against the velocity diagnostic, not a reassurance.

## 2. The GLORYS download works — the credentials file is the only problem

**Proven end to end, 2026-09-29.** An earlier note in this file said the boundary was blocked
for want of an account. **That was wrong**, and so was the stated cause. Three separate checks
were each wrong in turn:

1. `Test-Path "$HOME\.copernicusmarine\credentials"` returned `False` — a bad path test. The file
   **is** there.
2. `Get-Command copernicusmarine` and `import copernicusmarine` both failed — because they
   interrogated the **system** Python. The project uses `.venv\Scripts\python.exe`, which has
   `copernicusmarine` **2.4.1** and the `copernicusmarine.exe` console script. The code was
   right all along; the probe was not.
3. `copernicusmarine.climate.copernicus.eu` returns **NXDOMAIN** — true, and irrelevant. That
   host is retired. The live one is `data.marine.copernicus.eu`, which resolves and answers 200.

**What actually works, measured.** With credentials passed explicitly:

```
copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m \
  --variable thetao --minimum-longitude -58 --maximum-longitude -57 \
  --minimum-latitude 44 --maximum-latitude 44.5 \
  --start-datetime 2004-06-01 --end-datetime 2004-06-30 ...
  -> "Total size of the download: 22.20 KB" -> 33 KB .nc on disk
```

Read back with NCDatasets: `thetao` on a 50-level × 7 × 13 × 1 grid, and — worth noting for the
reader — **`thetao` is in `degrees_C`, not Kelvin.**

**The one real blocker.** `~/.copernicusmarine/credentials` is **base64-encoded**, not the plain
INI the client expects. Decoded it holds a normal `[credentials]` section with `username` and
`password`, and those values **are valid** — authentication succeeds when they are passed
explicitly. Left as-is, the client reports "requires a username and password" and aborts.

So this is a small, local, fixable defect rather than a missing account. Two options:

- Have `fetch_copernicus_physics_subset` detect a non-INI payload, base64-decode it, and parse
  the result. Robust to whoever wrote the file this way.
- Or rewrite the file as plain INI once. Simpler, but the encoding will come back if the same
  tool regenerates it.

The first is the better fix, and it is also the one that fails loudly: decoding a file that is
already plain INI must raise, not silently produce empty credentials.

**Also worth recording:** the `--x-min`/`--y-min` flags I first reached for do not exist. The
real ones are `--minimum-longitude`/`--maximum-longitude` and `--minimum-latitude`/
`--maximum-latitude`, and the time bounds are `--start-datetime`/`--end-datetime`, not
`--start-date`. A wrong flag produces a usage error naming the correct one, so this is
self-correcting.

**A caution on the climatology dataset.** The `_P1M-m` (monthly climatology) dataset spans
**2004 only**; asking for 2020 fails with a bounds error. The multi-year product
`GLOBAL_MULTIYEAR_PHY_001_030` has several datasets, and picking the right one is a decision,
not a default.

## 3. Tidal amplitudes are too small

`build_tidal_body_forcing` was returning a **velocity** where the model needs a **tendency** —
missing a factor of ω ≈ 1.4×10⁻⁴, so roughly 7000× too large. That is fixed: the forcing is
now `U_k·ω_k·cos(ω_k t)` in m/s², and `tidal_forcing_coefficients` accumulates
`s += -w * (re*sin(w*t) + im*cos(w*t))`.

The amplitudes are a separate problem. At the Bay of Fundy head the M2 amplitude computes to
**0.086 m/s** against a published **0.5–1.0 m/s** — too small by roughly an order of magnitude,
which is consistent with a decade of the wrong grid scale or a units slip somewhere upstream.

**The sponge has deliberately not been rewired onto the tidal forcing.** The previous ~7000×
error came from exactly that kind of rewiring, and until the amplitude is understood, connecting
the two would compound an unexplained discrepancy into a silent one.

## Benthic temperature ingest

**Blocked on the user obtaining access to the dataset.** The design is settled, so this is a
build task rather than a design task.

**The data.** `(lon, lat, depth, timestamp)`, multi-year, degrees Celsius, **no salinity**,
mostly *benthic* temperatures between 50 and 350 m. Some vertical profiles exist but are not
currently available.

**Why it is worth adding.** The depth range lands almost exactly on the settlement band
(`settlement_min_depth = -350.0`, `settlement_max_depth = -50.0`,
`settlement_max_temp = 6.0`). Settlement suitability is a function of tidally filtered
**near-bed** temperature and depth, so this dataset speaks directly to the term that gates
settlement — the most consequential input in the biology, and the one behind the 22/500 settled
figure.

**What it cannot do.** It is *not* a drop-in replacement for WOA23 as the initial stratification:
the model needs temperature at every level from 0 to 5000 m, and this is a benthic observation
confined to 50–350 m. It constrains the **shelf portion of the column**, nothing else.

**Design decisions already made.**

1. **Lower-shelf blend, not a replacement.** WOA23 retained for the surface and below 350 m so
   the column stays complete.
2. **Monthly climatology.** The run is climatological and repeating, so reduce the multi-year
   record to a 12-month climatology and drop it into the existing annual cycle.
3. **Keep WOA23 salinity.** None is supplied, but temperature inconsistent with salinity yields a
   density field consistent with neither — so the overlap must be *measured*, not assumed.
4. **Ingest-time validation.** Reject or flag non-finite, fill-value and out-of-range values.
   This class of bug has already bitten the project (`Missing` → `0.0` in the WOA23 reader), and
   a bad fill in a field feeding a 6 °C gate would do real damage.
5. **A blend weight, so three runs are comparable** — WOA23 only / user data only / blended — to
   show what actually changes.
6. **A startup overlap diagnostic** printing user value vs WOA23 at the same location as a
   difference map. Differences under roughly 0.5 °C mean use the data directly; large ones mean
   something is wrong with the file, the units, or the fill convention.
7. **2D now, 3D later.** The vertical profiles, when they arrive, slot into the same reader and
   upgrade the field without changing the interface.

**Needed to build it**

- File format (NetCDF variables and dimensions, or other) and the variable name.
- Depth convention (negative-down metres?) and whether the axis reaches the seabed or stops at 350 m.
- Fill/missing encoding: `NaN`, `-999`, `-1e30`, or empty.
- Timestamp units and the units attribute.

**Open questions for the user**

- Is the depth axis the *observation* depth (near-bed per cast) or a grid level? This decides
  whether the field is treated as benthic or as a level field.
- Is the multi-year record to be used as a climatology, or does a specific year matter for the
  climate-scenario comparison?

---

## Recently completed

Condensed. `git log` has the full history.

### 2026-09-29 — configuration is the single source of truth; dead code removed

- Struct fields renamed to match TOML keys; the driver's hand-parsed TOML reader deleted in
  favour of `configuration_to_options`. One interpreter of the config, not two.
- Colliding keys renamed: `[bathymetry] source` → `ingest_library`, `[atmosphere] source` →
  `forcing`. The exact bug this prevents has already happened once — `inshore_depth` and
  `shelf_slope` were set under `[data]` and read from `[bathymetry]`, and the two sections
  disagreed.
- Product ids live in `[data]`, mirroring the provenance registry one-to-one.
- `default.toml` switched from `"synthetic"` to real keyless products, so the baseline no
  longer runs on invented data.
- Deleted: ArgParse, the `data_mode` switch, `fetch_era5_atmospheric_forcing` (**a fabricator** —
  it generated a sinusoid from hardcoded constants `u_mean = 5.5`, `q_mean = 65.0` while its
  docstring claimed ERA5 reanalysis), `read_wind_stress` and `build_bulk_surface_flux` (neither
  had a runtime caller; the model builds its own flux).
- `fetch_hycom_boundary` added. See item 2 for its verification status.
- Accepted values for every `*_source` key are now documented in `configuration.jl` and in
  each TOML, in language meant for someone who is not a physical oceanographer.

### 2026-09-28 — visualization audit, biology, grid

- **Six computed quantities were never plotted** (`cohort_molt_fraction`, `degree_days_timeseries`,
  `temperatures`, `survival_probability`, `stage_survival`, `ascent_duration`): they existed only
  in the HTML payload, so the quantile/CDF stage model was unverifiable from output.
- **A figure fabricated its own data.** `plot_hydrodynamic_timeseries` never read the
  simulation — it took one static value and added a synthesised M2 cosine, and invented the
  time axis when only one snapshot existed. It now walks the real snapshots.
- Six new figures, including a time–depth diagram (the canonical Eulerian diagnostic, and the
  single most useful one missing) and a temperature–salinity scatter as a water-mass check.
- `--figures` / `--regenerate-figures`: redraws everything from stored results without
  re-simulating.
- Interactive HTML was 1,660 MB at production scale; now flat at ~19 MB via uniform decimation,
  with statistics still computed from full resolution.
- Larval biology reformulated to the quantile/CDF scheme: one `developmental_u ~ Uniform(0,1)`
  per larva, held for life and reused across all three transitions, so a fast developer stays
  fast within an individual.
- **Three silent-wrong-result bugs fixed:** `MoltCDF`'s default constructor was more specific
  than the `(base, cv)` method, so the lognormal got median `e^65` and a mean threshold of
  **1.8e28 degree-days** — nothing ever molted; settlement dispersion was applied twice; and
  `cohort_molt` column 1 was uninitialised memory.
- Vertical grid switched to an explicit CSV column (0, 10, 20, 40 … 5000 m): 10 cell centres in
  0–400 m and 4 in the cold intermediate layer, at 45% fewer cells than the stretched variant.
  The original 20-cell grid had its top cell spanning **0 to −512 m** — the entire shelf column
  was a single cell.
- Full 730-day run: exit 0, all 8 segments, **22/500 settled (4.4%)** — with the diverged tail
  present, so the count is more defensible than the trajectories.

### 2026-09-27 — verification pass, GPU fix, OOM fix

- GPU `InvalidIRError` root-caused: **any** `ContinuousForcing` carrying `field_dependencies`
  produces invalid LLVM IR in `gpu_compute_hydrostatic_free_surface_Gu!` on this stack, whether
  the dependency is a velocity or a tracer. Migrated the sponge to Oceananigans' own
  `Relaxation`, which needs no dependency declaration and is **numerically exact**.
- Segment 6 OOM: the flow interpolator materialised 378.7 GB of `Float64` for a 16.8 GB
  archive. Now selects *before* allocating, keeping ~22 GB. The OOM had also been **silently
  swallowed** — a `catch` fell back to an analytical jet, so the run reported recruitment numbers
  computed from a **synthetic current**. It now raises, with an explicit opt-in flag.
- `apply_climate_scenario!` replaced the baseline instead of perturbing it, discarding the
  stratification set immediately before it — including real WOA23 hydrography.
- EOS argument order was swapped (`ρ(S,T,…)` instead of `ρ(T,S,…)`), giving ≈1004 kg m⁻³ where
  shelf conditions must read ≈1027.
- Density halos were never filled, so N² was wrong by ~1000× at the first and last interior cells.
- `fill_seawater_density!` read uninitialised halo memory and could segfault the interpreter.
- `build_hydrodynamic_model` silently discarded a caller-supplied `forcing` by merging zero
  forcings over it.
- `ParticleTracking.stretched_tanh_z_faces` was called by the driver but **not defined at module
  scope**, so a fresh run could not build a grid at all.
- WOA23 obtained, read, interpolated; masked-value handling fixed. The long-standing "no network"
  diagnosis was **wrong** — what actually blocked it were a wrong revision directory (`A5B7` does
  not exist; the real one is `A5B4`), a 404 dataset id, and a 3-D assumption on a 4-D file.
- ETOPO: the dataset id was wrong (`nceiEtopo2022` 404s; the real one is `ETOPO_2022_v1_15s`)
  and the variable is `z`, not `altitude`.

---

## Do not regress

- **One place reads the TOML.** `configuration_to_options` in `src/config/configuration.jl`. Do not
  reintroduce per-key `get(cfg, …)` reads in the driver; a key added to the loader is otherwise
  silently ignored. This exact failure mode hid 8 keys for a long time.
- **A config value that never reaches the call site is a silent lie in `resolved_config.toml`.**
  The biology CV parameters were the example: loaded, defaulted to `0.25`, written to the
  provenance file — and never passed to `track_larval_cohort`, so every run used `0.0`. Adding a
  key to the loader and to `options_to_configuration` is **not** enough; trace the value all the
  way to the function that consumes it. When adding a knob, grep for the consumer in the same
  change.
- **A config key must be read from the section it is written in.** `inshore_depth` and
  `shelf_slope` were set under `[data]` and read from `[bathymetry]`. This shipped because both
  sections existed and neither complained.
- **A `cv_*` knob must control dispersion only, never the mean.** The previous Normal-then-clamp
  form broke this twice over: `clamp(Normal(hsi, cv·hsi), 0, 1)` pushed the expected settlement
  rate *below* `hsi` (0.90 instead of 1.00 at the optimum), and a `max(0, …)` on the mortality
  rate piled individuals onto a zero rate. Use `molt_cdf`, `draw_lognormal_mean` (mean preserving
  via the `−σ²/2` term), and `draw_beta_index`/`settlement_propensity` (log-odds scale, bounded,
  mean preserving). If a dispersion knob ever changes an expected outcome, that is a bug.
- **Before deleting anything, grep the top-level driver too — `ParticleTrackingRun.jl` and
  `*.jl` at the repo root.** It is not under `src/`, and it is 3273 lines. `read_wind_stress`
  and `build_bulk_surface_flux` were both deleted as dead code on the strength of a search that
  covered only `src/` and `test/`; the driver calls them at lines 397 and 634, so every run with
  atmospheric forcing was broken until they were restored. `src/` is not the codebase.
- **Check that a struct's parameterising constructor is actually the one being called.** `MoltCDF`
  had a default two-argument constructor more specific than any `(base, cv)` method, so the
  parameterisation was silently bypassed and the biology did nothing. Prefer a distinct factory
  function name (`molt_cdf(base, cv)`) over a same-arity constructor method.
- **Do not stack a new dispersion onto a function that already applies one.** The settlement
  propensity was added on top of `evaluate_settlement_suitability`'s internal `cv_settlement`,
  compounding the variance. When adding per-larva variation, first check what the callee already
  does and make the two modes mutually exclusive.
- **Never clamp a physical draw to put it in range; choose a distribution with the right support.**
  Clamping puts a point mass on the boundary that is indistinguishable from real structure, and
  it is not mean preserving. A clamped draw gave 2.4% of larvae a two-hour Zoea I stage.
- **Testset conventions** (`test/runtests.jl`): prefer one aggregate `@test all(...)` to a loop of
  many. A `for p in 1:n; @test …; end` block emits one test per iteration; testset 7 alone
  produced 305 of the old 732 this way, and a single failed assertion in a 260-test wall is hard
  to read. Only loop when the iterations are genuinely distinct cases. Numbers are append-only
  identifiers — do not renumber — and an assertion count is not a quality metric. There is
  deliberately no "skip the slow testsets" switch: an earlier attempt guarded only the first line
  of each body, so it printed "skipping" while still running everything.
- **Check the `LateralBoundaryRelaxation` construction against the new forcing path.** The sponge is
  now a `Relaxation`, so `sponge_relaxation_u` / `_v` are the reference formula only, not the code
  the model runs. Changing `compute_sponge_gamma` or the `ExponentialInflow` targets affects the
  live path through the `mask` and `target` arguments, not through those entry points.
- **Verify a remote dataset id and its variable name against the live catalogue before shipping a
  URL.** Three separate wrong assumptions were sitting in the data layer: `nceiEtopo2022` (404;
  the real id is `ETOPO_2022_v1_15s` and its variable is `z`, not `altitude`), WOA revision `A5B7`
  (no such directory; the real one is `A5B4`), and the WOA23 NetCDF layout itself. A one-line
  request against the `.das`/`.dds` endpoint settles it in seconds.
- **Verify a dataset can be *downloaded*, not just that its metadata loads.** HYCOM is the current
  example: the catalogue, the `.dds` and the `.html` all succeed, and every actual data request
  is refused. A metadata record is not a served file.
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
- **A 403 is not a 404.** A cloud bucket that answers 403 for absent keys breaks store/metadata
  discovery that would otherwise fall through to the next candidate.
- **`axes` on a `Field` is the interior.** Iterate `axes(field.data, k)` whenever halos matter —
  `∂z` reads one column past the interior. Halo cells of a freshly built model hold
  *uninitialised* memory: reading it and feeding it to a polynomial can overflow to Inf/NaN and
  take the process down. `fill_seawater_density!` now initialises its buffer to the EOS reference
  density and evaluates only the interior.
- **A standalone `HydrostaticFreeSurfaceModel` has no `boundary_conditions` field.** BCs live on
  the `Field`s and are fixed at construction, so post-hoc BC surgery raises `FieldError` and is
  not supported. Use `EarthSystemModel` + `PrescribedOcean` for prescribed surface state.
- **`indices(field, k)` returns a `Colon` in Oceananigans 0.111** — it is not a usable index
  range. Derive interior ranges from `axes(field.data, k)` and the grid's `Nx`/`Ny`/`Nz`.
- **Grid node accessors in 0.111 are `λᶠᵃᵃ`, `φᵃᶠᵃ`, `zᶜᶜᶜ`** (not `λᵃᶠᵃ`). A wrong one fails at
  `getproperty`.
- **Forcing closures must stay `isbits`, and must not take a `Symbol`.** `isbitstype(Symbol)` is
  `false`, so a `Symbol`-dispatched forcing method is unlowerable into a GPU kernel. Use concrete
  structs with one method per component (`sponge_relaxation_u` / `sponge_relaxation_v`) rather
  than a runtime selector.
- **Declare `field_dependencies`** whenever a forcing needs model state, or Oceananigans will call
  it with 4 arguments instead of 6. **But do not use it on the GPU path** — any
  `ContinuousForcing` carrying `field_dependencies` produces invalid LLVM IR in
  `gpu_compute_hydrostatic_free_surface_Gu!` on this stack. Use `Relaxation(; rate, mask, target)`
  for regional nudging; it reads the relaxed field through `r.relaxed` and needs no dependencies.
- **A `ContinuousForcing` is only the right tool for stateless-in-space tendencies.** Anything that
  needs to read the field it is forcing belongs in `Relaxation`. `field_dependencies` on a
  momentum forcing is also circular — the sponge reads the very `u`/`v` being updated.
- **A passing `configs/default.toml --gpu` run does not mean the GPU path works.** That config
  disables the sponge and the tides, so it never builds a forcing with field dependencies. When
  checking GPU health, reproduce the `snowcrab` forcing combination, not the `default` one.
- **Keep the honest-failure policy.** A run must never claim a dataset-forced state it did not
  read; raise instead, or make the synthetic choice explicit and recorded. The provenance
  registry's `keyless` flag means "retrievable with only a network connection" — do not set it
  true on the strength of a catalogue entry alone.
- **A test run that is killed leaves a root-level `ParticleTracking.toml` behind**, which the
  loader will then silently adopt for any later run that passes no `--config`. It is in
  `.gitignore`; if you see one, delete it.
- **`_grid.jld2` sidecars** are the only supported way to read grid geometry back; deserialising
  `serialized/grid` yields an unusable `JLD2.ReconstructedStatic`.
- Geography comes only from TOML `[domain]` / `[boundaries]`; no hard-coded domain constants.
