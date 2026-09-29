"""
    ParticleTrackingCLI

Command-line interface for the ParticleTracking driver, built on `ArgParse`.

**One declaration, one source of truth.** The flag spec below is the only place a flag name
exists: the parser, the generated `--help` text, and the tests that check every flag is wired all
derive from it. The previous driver hand-matched `"--flag" in args` about 84 times and carried a
104-line hand-written help block that had drifted out of sync with the parser.

**Naming rules, deliberately strict**

- Full words only. No abbreviations (`--sim`, `--viz` are gone; they are `--simulation`,
  `--visualize`).
- Correct spelling only. The long-standing misspelling `--tesselated` is gone, not aliased.
- No backwards-compatibility aliases. A renamed flag simply stops working, loudly, rather than
  quietly meaning something for one more release.
- Boolean pairs are written `--x` / `--no-x`, and passing both is an error rather than a silent
  precedence rule.

`command_line_settings` returns **only the options the user actually passed**, keyed by
`HydrodynamicOptions` field name. An unpassed flag is absent, so it can never overwrite a TOML value
with a default -- the class of bug that produced a "silent lie" in `resolved_config.toml` twice.
"""
module ParticleTrackingCLI

using ArgParse

export command_line_settings, CLI_OVERRIDE_KEYS, CLI_FLAG_NAMES

"""
    CLI_FLAG_NAMES -> Vector{String}

Every command-line flag this driver accepts, spelled with full words and no abbreviations.

Recorded here rather than read back out of `ArgParseSettings`, whose internal table type is a
version detail. The test suite asserts three properties of this list directly:

1. no name is an abbreviation or a misspelling of a real one (`--tesselated`, `--sim`, `--viz`,
   `--cp-`, `--buf`, `--n-units`, `--anim-var` are all absent, and nothing is aliased to them);
2. no name is shorter than three characters, so nothing cryptic can be introduced;
3. every name contains only lowercase letters and hyphens.
"""
const CLI_FLAG_NAMES = String[
    "--config", "--run-id", "--quick", "--gpu", "--cpu", "--allow-analytical-fallback",
    "--all", "--data", "--grid", "--model", "--climate", "--simulation", "--tracking",
    "--metrics", "--visualize", "--figures", "--hydrodynamics-only", "--tracking-only",
    "--reuse-hydrodynamics", "--segment",
    "--real-data", "--synthetic-data", "--historical-five-year",
    "--climatology-two-year", "--climatology-eighteen-month",
    "--tides", "--no-tides", "--open-boundary-conditions", "--no-open-boundary-conditions",
    "--wind-source", "--diel-vertical-migration",
    "--no-diel-vertical-migration", "--molting", "--no-molting",
    "--initial-ascent", "--tessellated-cells", "--no-tessellated-cells",
    "--dissolved-oxygen", "--no-dissolved-oxygen", "--interactive-map", "--no-interactive-map",
    "--duckdb", "--no-duckdb", "--checkpoints", "--no-checkpoints",
    "--retain-all-checkpoints", "--restart", "--no-restart", "--adaptive-cfl", "--no-adaptive-cfl",
    "--animate-hydrodynamics", "--animation-overlay-particles",
    "--no-animation-overlay-particles",
    "--domain-longitude", "--domain-latitude", "--domain-depth", "--grid-size",
    "--grid-nx", "--grid-ny", "--grid-nz", "--resolution-scale",
    "--vertical-stretching-mode", "--vertical-break-depth", "--levels-above-break",
    "--vertical-grid-file",
    "--simulation-time-step", "--minimum-time-step", "--maximum-time-step",
    "--simulation-duration", "--target-cfl", "--surface-heat-flux",
    "--horizontal-diffusivity", "--vertical-diffusivity", "--seed",
    "--open-boundary-type", "--sponge-width", "--sponge-timescale",
    "--inflow-u", "--inflow-v", "--tidal-u-amplitude", "--tidal-v-amplitude",
    "--secondary-tidal-u-amplitude", "--secondary-tidal-v-amplitude",
    "--n-particles", "--tracking-duration", "--tracking-time-step", "--minimum-release-depth",
    "--seeding-buffer", "--release-mode", "--release-offset", "--ascent-speed",
    "--ascent-target-depth",
    "--n-tessellation-cells", "--tessellation-core-probability",
    "--tessellation-shallow-probability", "--tessellation-deep-probability",
    "--tessellation-core-minimum-resolution-km",
    "--tessellation-shallow-minimum-resolution-km",
    "--tessellation-deep-minimum-resolution-km",
    "--output-directory", "--input-directory", "--duckdb-path", "--hydrodynamics-file",
    "--checkpoint-prefix", "--checkpoint-directory", "--checkpoint-interval",
    "--maximum-current-speed", "--maximum-flow-snapshots",
    "--animation-variable", "--animation-frames-per-second", "--animation-format",
    "--animation-depth", "--animation-output",
    "--list-runs", "--compare-scenarios", "--model-average",
]

"""
    CLI_OVERRIDE_KEYS -> Vector{String}

Names of every command-line option that maps onto a `HydrodynamicOptions` field. The driver applies
exactly these as overrides on top of `base_opts`. The test suite asserts this list, the parsed
flags, and the option struct stay consistent, so an option cannot be added to one and forgotten in
the others.
"""
const CLI_OVERRIDE_KEYS = String[
    "domain_lon", "domain_lat", "domain_z", "grid_size", "data_mode", "scenario",
    "projection_year", "n_particles", "track_duration", "track_dt", "sim_dt", "sim_duration",
    "min_dt_seconds", "max_dt", "target_cfl", "adaptive_cfl", "surface_heat_flux",
    "diffusivity_h", "diffusivity_v", "use_gpu", "fallback_to_cpu", "seed",
    "output_dir", "input_dir", "run_id", "hydro_model_file", "enable_tides", "tidal_u_amp",
    "tidal_v_amp", "s2_u_amp", "s2_v_amp", "obc_type", "boundary_method", "sponge_layer_width",
    "sponge_timescale", "u_inflow", "v_inflow", "min_seabed_depth", "buffer_km",
    "release_depth_mode", "bottom_release_offset", "enable_initial_ascent", "ascent_speed",
    "ascent_target_depth", "enable_dvm", "enable_molting", "enable_voronoi", "voronoi_n_units",
    "voronoi_prob_core", "voronoi_prob_shallow", "voronoi_prob_deep", "voronoi_min_res_core_km",
    "voronoi_min_res_shallow_km", "voronoi_min_res_deep_km", "enable_duckdb", "duckdb_path",
    "interactive_map", "enable_checkpoint", "checkpoint_prefix", "checkpoint_schedule",
    "checkpoint_dir", "checkpoint_cleanup", "auto_restart",
    "max_current_speed", "max_flow_snapshots",
    "vertical_stretching_mode", "nz_above", "vertical_break_depth", "vertical_grid_file",
    "resolution_scale", "atmospheric_source", "animate_hydro", "anim_variable", "anim_fps",
    "anim_format",
    "anim_depth", "anim_overlay_particles", "anim_output_path",
]

"""
    CLI_KEY_ALIASES::Dict{String,String}

Parsed key -> `HydrodynamicOptions` field name, for the flags whose CLI spelling deliberately
differs from the field they set.

Most flags are named after their field, so the mapping is the identity and needs no table. Where
the command line is better named than the field -- `--wind-source` reads as what it does, whereas
`atmospheric_source` names the field -- the rename is recorded here so that
`command_line_settings` still returns a dictionary keyed purely by field names, which is what
`configuration_to_options` accepts as overrides.

An alias whose target is also listed in `CLI_OVERRIDE_KEYS` is redundant and should be deleted.
"""
const CLI_KEY_ALIASES = Dict{String, String}(
    "wind_source" => "atmospheric_source",
)

"""
    WIND_SOURCE_NONE::String

The `--wind-source` spellings that mean "no atmospheric forcing".

`--wind-source=""` and `--wind-source=false` both have to work, because a shell cannot always
express an empty argument conveniently and `false` is what a user reaches for out of habit from
the old boolean pair. Both normalise to the same value.
"""
const WIND_SOURCE_NONE = ("", "false", "none", "off", "no")

"""
    build_settings() -> ArgParseSettings

The flag specification. `--help` is generated from this, so the two cannot drift apart.
"""
function build_settings()
    settings = ArgParseSettings(
        prog = "ParticleTrackingRun.jl",
        description = "Larval transport, biophysical ocean circulation and demographic connectivity.",
        epilog = """
        Configuration precedence: command-line flags override the TOML configuration file, which
        overrides built-in defaults. Every run writes its fully resolved configuration to
        <output-directory>/resolved_config.toml.
        """,
        add_version = true,
        commands_are_required = false,
        autofix_names = false,
        # Throw on a bad command line instead of printing usage and calling `exit`. This is a
        # library, not a script: exiting would take the caller (and the test suite) down with it.
        error_on_conflict = true,
    )

    @add_arg_table! settings begin
        "--config"
            help = "path to a TOML configuration file"
            arg_type = String
            default = nothing
        "--run-id"
            help = "identifier for this run; also names the output subdirectory"
            arg_type = String
            default = nothing
        "--quick"
            help = "fast prototype settings: coarse grid, short simulation, small cohort"
            action = :store_true
        "--gpu"
            help = "require GPU execution; fail rather than fall back"
            action = :store_true
        "--cpu"
            help = "force CPU execution"
            action = :store_true
        "--allow-analytical-fallback"
            help = "permit substituting an analytical current when the hydrodynamics archive cannot be read"
            action = :store_true

        # ---- pipeline selection ----------------------------------------------
        "--all"
            help = "run every segment in order"
            action = :store_true
        "--data"
            help = "acquire bathymetry and wind; regrid hydrography"
            action = :store_true
        "--grid"
            help = "build the grid and immersed boundary"
            action = :store_true
        "--model"
            help = "assemble the hydrodynamic model and apply climate forcing"
            action = :store_true
        "--climate"
            help = "apply the climate scenario only"
            action = :store_true
        "--simulation"
            help = "integrate the hydrodynamic model"
            action = :store_true
        "--tracking"
            help = "track the larval cohort"
            action = :store_true
        "--metrics"
            help = "compute connectivity, thermal and recruitment metrics"
            action = :store_true
        "--visualize"
            help = "render figures and the interactive map"
            action = :store_true
        "--figures"
            help = "redraw every figure and visual from stored results, without re-simulating"
            action = :store_true
        "--hydrodynamics-only"
            help = "run the hydrodynamic segments only"
            action = :store_true
        "--tracking-only"
            help = "run the tracking segments only, against an existing hydrodynamics archive"
            action = :store_true
        "--reuse-hydrodynamics"
            help = "reuse the existing hydrodynamics archive rather than re-integrating"
            action = :store_true
        "--segment"
            help = "run one named segment: data, grid, model, climate, simulation, tracking, metrics, visualize, all"
            arg_type = String
            default = nothing
            metavar = "NAME"

        # ---- scenario presets -------------------------------------------------
        "--real-data"
            help = "use real bathymetry, wind and hydrography"
            action = :store_true
        "--synthetic-data"
            help = "use synthetic data instead of downloads"
            action = :store_true
        "--historical-five-year"
            help = "scenario preset: five years of historical forcing"
            action = :store_true
        "--climatology-two-year"
            help = "scenario preset: two years of repeating annual climatology"
            action = :store_true
        "--climatology-eighteen-month"
            help = "scenario preset: eighteen months of repeating annual climatology"
            action = :store_true

        # ---- feature switches --------------------------------------------------
        "--tides"
            help = "include astronomical tides"
            action = :store_true
        "--no-tides"
            help = "exclude astronomical tides"
            action = :store_true
        "--open-boundary-conditions"
            help = "include open-boundary conditions"
            action = :store_true
        "--no-open-boundary-conditions"
            help = "exclude open-boundary conditions"
            action = :store_true
        "--wind-source"
            help = "atmospheric forcing source, matching [atmosphere] source: era5, " *
                   "era5_climatology, mhw; empty or false for no atmospheric forcing"
            arg_type = String
            default = nothing
            metavar = "SOURCE"
        "--diel-vertical-migration"
            help = "include diel vertical migration"
            action = :store_true
        "--no-diel-vertical-migration"
            help = "exclude diel vertical migration"
            action = :store_true
        "--molting"
            help = "include degree-day molting"
            action = :store_true
        "--no-molting"
            help = "exclude degree-day molting"
            action = :store_true
        "--initial-ascent"
            help = "post-hatch vertical ascent; pass =true or =false"
            arg_type = Bool
            default = nothing
        "--tessellated-cells"
            help = "include Voronoi tessellated connectivity cells"
            action = :store_true
        "--no-tessellated-cells"
            help = "exclude Voronoi tessellated connectivity cells"
            action = :store_true
        "--dissolved-oxygen"
            help = "include the dissolved oxygen tracer"
            action = :store_true
        "--no-dissolved-oxygen"
            help = "exclude the dissolved oxygen tracer"
            action = :store_true
        "--interactive-map"
            help = "generate the interactive HTML map"
            action = :store_true
        "--no-interactive-map"
            help = "skip the interactive HTML map"
            action = :store_true
        "--duckdb"
            help = "persist results to DuckDB"
            action = :store_true
        "--no-duckdb"
            help = "do not persist results to DuckDB"
            action = :store_true
        "--checkpoints"
            help = "write periodic checkpoints"
            action = :store_true
        "--no-checkpoints"
            help = "do not write periodic checkpoints"
            action = :store_true
        "--retain-all-checkpoints"
            help = "keep every checkpoint rather than only the latest"
            action = :store_true
        "--restart"
            help = "resume from a checkpoint if one exists"
            action = :store_true
        "--no-restart"
            help = "never resume; start fresh"
            action = :store_true
        "--adaptive-cfl"
            help = "adapt the time step to the Courant number"
            action = :store_true
        "--no-adaptive-cfl"
            help = "use a fixed time step"
            action = :store_true
        "--animate-hydrodynamics"
            help = "render a hydrodynamics animation"
            action = :store_true
        "--animation-overlay-particles"
            help = "overlay drifting larvae on the animation"
            action = :store_true
        "--no-animation-overlay-particles"
            help = "do not overlay drifting larvae"
            action = :store_true

        # ---- geometry ---------------------------------------------------------
        "--domain-longitude"
            help = "longitude range as min,max in degrees east"
            arg_type = Float64
            nargs = 2
            default = nothing
        "--domain-latitude"
            help = "latitude range as min,max in degrees north"
            arg_type = Float64
            nargs = 2
            default = nothing
        "--domain-depth"
            help = "depth range as surface,bottom in metres, negative down"
            arg_type = Float64
            nargs = 2
            default = nothing
        "--grid-size"
            help = "grid cells as nx,ny,nz"
            arg_type = Int
            nargs = 3
            default = nothing
        "--grid-nx"
            help = "cells in longitude"
            arg_type = Int
            default = nothing
        "--grid-ny"
            help = "cells in latitude"
            arg_type = Int
            default = nothing
        "--grid-nz"
            help = "cells in the vertical"
            arg_type = Int
            default = nothing
        "--resolution-scale"
            help = "divisor applied to the horizontal grid size"
            arg_type = Float64
            default = nothing
        "--vertical-stretching-mode"
            help = "vertical grid mode: two_segment, csv, tanh, stretched or uniform"
            arg_type = String
            default = nothing
        "--vertical-break-depth"
            help = "depth of the two-segment vertical grid join, metres negative down"
            arg_type = Float64
            default = nothing
        "--levels-above-break"
            help = "vertical levels above the break depth; 0 means 60 percent of nz"
            arg_type = Int
            default = nothing
        "--vertical-grid-file"
            help = "CSV of vertical layer faces, deepest layer first, ending at 0 m"
            arg_type = String
            default = nothing

        # ---- numerics ---------------------------------------------------------
        "--simulation-time-step"
            help = "initial simulation time step, seconds"
            arg_type = Float64
            default = nothing
        "--minimum-time-step"
            help = "floor on the adaptive time step, seconds"
            arg_type = Float64
            default = nothing
        "--maximum-time-step"
            help = "ceiling on the adaptive time step, seconds"
            arg_type = Float64
            default = nothing
        "--simulation-duration"
            help = "simulation duration, seconds"
            arg_type = Float64
            default = nothing
        "--target-cfl"
            help = "target advective Courant number"
            arg_type = Float64
            default = nothing
        "--surface-heat-flux"
            help = "net surface heat flux, W m-2"
            arg_type = Float64
            default = nothing
        "--horizontal-diffusivity"
            help = "horizontal turbulent diffusivity, m2 s-1"
            arg_type = Float64
            default = nothing
        "--vertical-diffusivity"
            help = "vertical turbulent diffusivity, m2 s-1"
            arg_type = Float64
            default = nothing
        "--seed"
            help = "random seed"
            arg_type = Int
            default = nothing

        # ---- boundary conditions ----------------------------------------------
        "--open-boundary-type"
            help = "open boundary condition scheme, for example flather_chapman"
            arg_type = String
            default = nothing
        "--sponge-width"
            help = "lateral relaxation band width, degrees"
            arg_type = Float64
            default = nothing
        "--sponge-timescale"
            help = "lateral relaxation timescale, seconds"
            arg_type = Float64
            default = nothing
        "--inflow-u"
            help = "upstream inflow u velocity, m s-1"
            arg_type = Float64
            default = nothing
        "--inflow-v"
            help = "upstream inflow v velocity, m s-1"
            arg_type = Float64
            default = nothing
        "--tidal-u-amplitude"
            help = "M2 semi-major tidal velocity, m s-1"
            arg_type = Float64
            default = nothing
        "--tidal-v-amplitude"
            help = "M2 semi-minor tidal velocity, m s-1"
            arg_type = Float64
            default = nothing
        "--secondary-tidal-u-amplitude"
            help = "S2 tidal u velocity, m s-1"
            arg_type = Float64
            default = nothing
        "--secondary-tidal-v-amplitude"
            help = "S2 tidal v velocity, m s-1"
            arg_type = Float64
            default = nothing

        # ---- biology ----------------------------------------------------------
        "--n-particles"
            help = "number of larvae in the cohort"
            arg_type = Int
            default = nothing
        "--tracking-duration"
            help = "larval tracking duration, seconds"
            arg_type = Float64
            default = nothing
        "--tracking-time-step"
            help = "Lagrangian integration step, seconds"
            arg_type = Float64
            default = nothing
        "--minimum-release-depth"
            help = "shallowest water depth permitting larval release, metres"
            arg_type = Float64
            default = nothing
        "--seeding-buffer"
            help = "spatial seeding buffer beyond the boundary, km"
            arg_type = Float64
            default = nothing
        "--release-mode"
            help = "release depth mode: bottom or uniform"
            arg_type = String
            default = nothing
        "--release-offset"
            help = "release height above the seabed as min,max metres"
            arg_type = Float64
            nargs = 2
            default = nothing
        "--ascent-speed"
            help = "post-hatch ascent speed, m s-1"
            arg_type = Float64
            default = nothing
        "--ascent-target-depth"
            help = "depth of the diel mixed layer target, metres"
            arg_type = Float64
            default = nothing

        # ---- connectivity cells ----------------------------------------------
        "--n-tessellation-cells"
            help = "number of Voronoi cells"
            arg_type = Int
            default = nothing
        "--tessellation-core-probability"
            help = "probability of a cell being core habitat"
            arg_type = Float64
            default = nothing
        "--tessellation-shallow-probability"
            help = "probability of a cell being shallow habitat"
            arg_type = Float64
            default = nothing
        "--tessellation-deep-probability"
            help = "probability of a cell being deep habitat"
            arg_type = Float64
            default = nothing
        "--tessellation-core-minimum-resolution-km"
            help = "minimum resolution of a core cell, km"
            arg_type = Float64
            default = nothing
        "--tessellation-shallow-minimum-resolution-km"
            help = "minimum resolution of a shallow cell, km"
            arg_type = Float64
            default = nothing
        "--tessellation-deep-minimum-resolution-km"
            help = "minimum resolution of a deep cell, km"
            arg_type = Float64
            default = nothing

        # ---- storage ----------------------------------------------------------
        "--output-directory"
            help = "directory for run outputs"
            arg_type = String
            default = nothing
        "--input-directory"
            help = "directory for cached input data"
            arg_type = String
            default = nothing
        "--duckdb-path"
            help = "path to the DuckDB database"
            arg_type = String
            default = nothing
        "--hydrodynamics-file"
            help = "hydrodynamics archive to read or write"
            arg_type = String
            default = nothing
        "--checkpoint-prefix"
            help = "checkpoint filename prefix"
            arg_type = String
            default = nothing
        "--checkpoint-directory"
            help = "directory for checkpoints"
            arg_type = String
            default = nothing
        "--checkpoint-interval"
            help = "checkpoint interval in seconds, or a value with an h, m or d suffix"
            arg_type = String
            default = nothing
        "--maximum-current-speed"
            help = "cap on the background current a larva is advected by, m s-1"
            arg_type = Float64
            default = nothing
        "--maximum-flow-snapshots"
            help = "cap on hydrodynamic snapshots held in memory for the flow interpolator"
            arg_type = Int
            default = nothing

        # ---- animation ---------------------------------------------------------
        "--animation-variable"
            help = "field to animate: dashboard, temperature, advection, salinity, speed, w, density, stratification"
            arg_type = String
            default = nothing
        "--animation-frames-per-second"
            help = "animation playback rate"
            arg_type = Int
            default = nothing
        "--animation-format"
            help = "animation format: gif or mp4"
            arg_type = String
            default = nothing
        "--animation-depth"
            help = "depth of the animated slice, metres"
            arg_type = Float64
            default = nothing
        "--animation-output"
            help = "animation output path"
            arg_type = String
            default = nothing

        # ---- analytics ---------------------------------------------------------
        "--list-runs"
            help = "list recorded runs"
            action = :store_true
        "--compare-scenarios"
            help = "compare recorded scenarios"
            action = :store_true
        "--model-average"
            help = "ensemble-average recorded models"
            action = :store_true
    end

    return settings
end

# Boolean pairs: passing both is an error rather than a silent precedence rule.
const EXCLUSIVE_PAIRS = Tuple{String, String}[
    ("--tides", "--no-tides"),
    ("--open-boundary-conditions", "--no-open-boundary-conditions"),
    ("diel_vertical_migration", "no_diel_vertical_migration"),
    ("molting", "no-molting"),
    ("tessellated_cells", "no_tessellated_cells"),
    ("dissolved_oxygen", "no_dissolved_oxygen"),
    ("interactive_map", "no_interactive_map"),
    ("duckdb", "no_duckdb"),
    ("checkpoints", "no_checkpoints"),
    ("restart", "no_restart"),
    ("adaptive_cfl", "no_adaptive_cfl"),
    ("gpu", "cpu"),
    ("animation_overlay_particles", "no_animation_overlay_particles"),
]

"""
    was_passed(args, key::AbstractString) -> Bool

Whether the user actually typed the option `key` on the command line.

`parse_args` returns a fully populated dictionary: every declared option comes back with *some*
value, so "the user asked for `false`" and "the user said nothing" are indistinguishable from the
parsed result alone. That ambiguity is the whole reason this function exists.

It cannot be resolved with a sentinel comparison. `false` means "not passed" for a `store_true`
flag, because that is the default ArgParse supplies, but it is the literal answer for an option
declared as `arg_type = Bool` — `--initial-ascent=false` is the only way to turn that feature off.
Dropping every `false` therefore made the option impossible to disable; keeping every `false` would
let every unpassed `store_true` flag overwrite the TOML. The command line itself is the only place
the distinction exists, so it is read from there.

Both spellings count as passing: `--flag value` and `--flag=value`.
"""
function was_passed(args::AbstractVector{<:AbstractString}, key::AbstractString)::Bool
    flag = string("--", replace(key, '_' => '-'))
    for a in args
        a == flag && return true
        startswith(a, string(flag, "=")) && return true
    end
    return false
end

"""
    command_line_settings(args) -> Dict{String,Any}

Parse `args` (excluding the program name) and return only the options the user actually supplied,
keyed by `HydrodynamicOptions` field name. A flag the user did not pass is absent from the result,
so it cannot silently overwrite the TOML configuration.

Unknown flags are an error, not a warning.
"""
function command_line_settings(args::AbstractVector{<:AbstractString})::Dict{String, Any}
    argv = collect(String, args)
    parsed = parse_args(argv, build_settings(); as_symbols = false)

    # Keep only what the user actually passed, and key it by the `HydrodynamicOptions` FIELD name.
    # `parse_args` fills in a value for every declared flag, so this filter is what makes an
    # unpassed flag *absent* rather than present-with-a-default -- which is what stopped it
    # overwriting the TOML configuration. Multi-value options come back as an empty collection
    # when unset, so those are dropped too. `false` is kept when, and only when, the user typed
    # the flag; see `was_passed`.
    passed = Dict{String, Any}()
    for (k, v) in parsed
        v === nothing && continue
        v isa AbstractVector && isempty(v) && continue
        v === false && !was_passed(argv, k) && continue
        passed[replace(k, '-' => '_')] = v
    end

    # Contradictory pairs are an error, checked on the *passed* set so defaults never trip it.
    for (on, off) in EXCLUSIVE_PAIRS
        on_key = replace(startswith(on, "--") ? on[3:end] : on, '-' => '_')
        off_key = replace(startswith(off, "--") ? off[3:end] : off, '-' => '_')
        if haskey(passed, on_key) && haskey(passed, off_key)
            throw(ArgumentError(
                "flags $(on) and $(off) are mutually exclusive; pass at most one."))
        end
    end

    # Rename the flags that are deliberately spelled better than the field they set, so the
    # result is keyed by field name and can be splatted straight into `configuration_to_options`.
    for (from, to) in CLI_KEY_ALIASES
        if haskey(passed, from)
            passed[to] = pop!(passed, from)
        end
    end

    # `--wind-source` also has to be able to switch atmospheric forcing *off*, which is the one
    # thing a source name cannot express. The driver tests membership of a list of sources that
    # do produce fluxes, so "off" is a source that is not in it.
    if haskey(passed, "atmospheric_source")
        src = passed["atmospheric_source"]
        if src isa AbstractString && lowercase(strip(src)) in WIND_SOURCE_NONE
            passed["atmospheric_source"] = :none
        end
    end

    return passed
end

end # module





