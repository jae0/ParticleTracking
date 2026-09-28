"""
    src/output/biology_diagnostics.jl

Cohort-level diagnostics for the larval biology that the main figure set did not cover.

The particle-tracking model computes several quantities that were previously only embedded in the
interactive HTML payload and never drawn as a figure, which made the stage model hard to check:

- `cohort_molt_fraction` — the CDF-based probability that a cohort has completed each transition,
  on the same 0-to-1 scale as the per-larva developmental quantiles;
- `degree_days_timeseries` — the cumulative thermal time each larva actually accumulated, which is
  the x-axis on which every stage transition is defined;
- `temperatures` — the temperature history the larvae were advected through;
- `survival_probability` and `stage_survival` — cohort survival and the fraction alive at each
  stage transition;
- the resolved `molt_schedule` — the three threshold CDFs themselves.

Together these let the stage model be verified from the output rather than inferred from it.
"""

using CairoMakie
using CairoMakie: Figure, Axis, heatmap!, lines!, scatter!, hlines!, vlines!, band!,
                  Legend, hist!, barplot!, text!
using Statistics: mean, var
using ParticleTracking: molt_probability, threshold, update_larval_stage

const _STAGE_ORDER = (:zoea1, :zoea2, :megalopa, :instar1_settled)
const _STAGE_COLORS = (:dodgerblue, :darkorange, :mediumseagreen, :crimson)

"""
    plot_molt_progression(trajs; output_path, title)

Diagnostic for the quantile/CDF stage model: how the cohort advances through Zoea I -> Zoea II ->
Megalopa -> Instar I over the simulated period.

# Panels
1. **Cohort transition probability** — `cohort_molt_fraction` against time, one curve per
   transition, on a 0-to-1 axis. This is the direct readout of the CDF the per-larva quantiles are
   compared against.
2. **Realised stage composition** — the fraction of the cohort in each stage over time, as stacked
   bands. Should track panel 1: a cohort fraction that has completed a transition should match the
   fraction of larvae at or beyond that stage.
3. **Developmental quantile spread** — molting degree-day against developmental quantile `u` for a
   sample of particles, with the resolved schedule's median and base threshold marked. Confirms the
   ordering (small `u` moults early) and the shape of the distribution.

# Inputs
- `trajs`: Trajectory NamedTuple from `track_larval_cohort`.
- `output_path`: PNG path, or `nothing` to skip saving.
"""
function plot_molt_progression(
    trajs::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/mole_progression.png",
    title::AbstractString = "Stage Progression: Cohort CDF and Realised Composition",
    n_sample::Int = 120
)
    times = [t / 86400.0 for t in trajs.times]
    n_t = length(times)
    n_p = size(trajs.lons, 1)
    cm = trajs.cohort_molt_fraction          # 3 x n_t
    fig = Figure(size = (1400, 450), fontsize = 12)
    labels = ["Zoea I → Zoea II", "Zoea II → Megalopa", "Megalopa → Instar I"]

    # Panel 1: cohort transition probability from the CDF
    ax1 = Axis(fig[1, 1], title = "Cohort transition probability (CDF)",
               xlabel = "Time (days)", ylabel = "Fraction of cohort completed",
               limits = (0, isempty(times) ? 1 : times[end], 0, 1))
    for k in 1:3
        size(cm, 1) >= k || break
        lines!(ax1, times, cm[k, :], color = _STAGE_COLORS[k], linewidth = 2.5, label = labels[k])
    end
    Legend(fig[1, 2], ax1; position = :rt)

    # Panel 2: realised stage composition, stacked
    ax2 = Axis(fig[1, 3], title = "Realised stage composition",
               xlabel = "Time (days)", ylabel = "Fraction of cohort", limits = (0, isempty(times) ? 1 : times[end], 0, 1))
    stage_frac = zeros(length(_STAGE_ORDER), n_t)
    for j in 1:n_t, s in _STAGE_ORDER
        stage_frac[findfirst(==(s), _STAGE_ORDER), j] = count(==(s), view(trajs.stages, :, j)) / n_p
    end
    lower = zeros(n_t)
    for i in eachindex(_STAGE_ORDER)
        band!(ax2, times, lower, lower .+ stage_frac[i, :], color = _STAGE_COLORS[i],
              alpha = 0.75, label = String(_STAGE_ORDER[i]))
        lower = lower .+ stage_frac[i, :]
    end
    Legend(fig[1, 4], ax2; position = :rt)

    # Panel 3: molting degree-day against developmental quantile
    has_sched = hasproperty(trajs, :molt_schedule) && trajs.molt_schedule !== nothing &&
                length(trajs.molt_schedule) == 3
    ax3 = Axis(fig[1, 5],
               title = has_sched ? "Molt threshold vs quantile (dashed: base 65/130/200)" :
                                   "Molt threshold vs quantile (no schedule recorded)",
               xlabel = "Developmental quantile u", ylabel = "Molt degree-day")
    if has_sched
        sched = trajs.molt_schedule
        us = range(0.02, 0.98, length = 60)
        for k in 1:3
            lines!(ax3, us, [threshold(sched[k], u) for u in us],
                   color = _STAGE_COLORS[k], linewidth = 2, label = labels[k])
        end
        for thr in (65.0, 130.0, 200.0)
            hlines!(ax3, [thr], color = (:grey, 0.5), linestyle = :dash)
        end
        Legend(fig[1, 6], ax3; position = :rt)
    end

    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

"""
    plot_degree_day_growth(trajs; output_path, title)

Cumulative thermal degree-day actually accumulated by each larva, against time, with the three
base thresholds marked.

This is the x-axis of the whole stage model — every transition is defined on it — so a run whose
curves stay flat has larvae that never develop, and one whose curves diverge is doing so because
larvae saw different temperatures. A 20-day cohort should reach ~65 DD and molt; a 730-day one
saturates the first two thresholds and is limited by the third.

# Inputs
- `trajs`: Trajectory NamedTuple.
- `output_path`: PNG path, or `nothing`.
"""
function plot_degree_day_growth(
    trajs::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/degree_day_growth.png",
    title::AbstractString = "Cumulative Degree-Day Accumulation",
    n_sample::Int = 200
)
    times = [t / 86400.0 for t in trajs.times]
    dds = trajs.degree_days_timeseries
    n_p = size(dds, 1)
    fig = Figure(size = (1200, 500), fontsize = 12)

    ax1 = Axis(fig[1, 1],
               title = "Cumulative degree-day, cohort mean (dashed: 65 / 130 / 200 DD)",
               xlabel = "Time (days)", ylabel = "Cumulative degree-day (°C·d)")
    idx = unique(round.(Int, range(1, n_p, length = min(n_sample, n_p))))
    for p in idx
        lines!(ax1, times, @view(dds[p, :]), color = (:steelblue, 0.25), linewidth = 0.8)
    end
    for (thr, lab) in ((65.0, "Z I→II"), (130.0, "II→M"), (200.0, "M→I"))
        hlines!(ax1, [thr], color = (:crimson, 0.55), linestyle = :dash, label = lab)
    end
    lines!(ax1, times, [mean(@view dds[:, j]) for j in axes(dds, 2)],
           color = :black, linewidth = 3, label = "cohort mean")
    axislegend(ax1; position = :rt)

    # Panel 2: the thermal forcing the larvae actually saw
    ax2 = Axis(fig[1, 2], title = "Temperature history", xlabel = "Time (days)",
               ylabel = "Temperature (°C)")
    temps = trajs.temperatures
    for p in idx
        lines!(ax2, times, @view(temps[p, :]), color = (:darkorange, 0.2), linewidth = 0.8)
    end
    lines!(ax2, times, [mean(@view temps[:, j]) for j in axes(temps, 2)],
           color = :black, linewidth = 3)

    # Panel 3: distribution of final accumulated degree-days
    ax3 = Axis(fig[1, 3], title = "Final degree-day per larva",
               xlabel = "Cumulative degree-day (°C·d)", ylabel = "Larvae")
    final_dd = [dds[p, end] for p in 1:n_p]
    hist!(ax3, final_dd; bins = min(40, max(10, n_p ÷ 4)), color = :steelblue)
    for thr in (65.0, 130.0, 200.0)
        vlines!(ax3, [thr], color = (:crimson, 0.5), linestyle = :dash)
    end

    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

"""
    plot_survival_curves(trajs; output_path, title)

Population-level survival diagnostics.

# Panels
1. **Survival probability** — `survival_probability` per larva over time, plus the cohort mean. In
   stochastic mode the curves are step functions (an individual dies once and stays dead); in
   mean-field mode they decay smoothly and every larva shares one curve.
2. **Settlement and mortality outcome** — how the cohort ended, as a single stacked bar. Distinguishes
   settled, still pelagic, and dead, which no other figure reports directly.
3. **Stage-transition survival** — `stage_survival`, the fraction of each cohort still alive at the
   moment it changed stage, as a bar per transition. A large drop at a transition indicates
   stage-specific mortality not captured by the background rate.

# Inputs
- `trajs`: Trajectory NamedTuple.
- `output_path`: PNG path, or `nothing`.
"""
function plot_survival_curves(
    trajs::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/survival_curves.png",
    title::AbstractString = "Survival, Settlement and Stage-Transition Mortality"
)
    times = [t / 86400.0 for t in trajs.times]
    surv = trajs.survival_probability
    n_p = size(surv, 1)
    fig = Figure(size = (1300, 450), fontsize = 12)

    # Panel 1: survival curves
    ax1 = Axis(fig[1, 1], title = "Survival probability", xlabel = "Time (days)",
               ylabel = "Survival probability", limits = (0, isempty(times) ? 1 : times[end], 0, 1))
    idx = unique(round.(Int, range(1, n_p, length = min(200, n_p))))
    for p in idx
        lines!(ax1, times, @view(surv[p, :]), color = (:seagreen, 0.15), linewidth = 0.8)
    end
    lines!(ax1, times, [mean(@view surv[:, j]) for j in axes(surv, 2)],
           color = :black, linewidth = 3, label = "cohort mean")
    Legend(fig[1, 2], ax1; position = :rt)

    # Panel 2: final disposition of the cohort. Makie's barplot does not accept string
    # categories in this version, so plot numeric positions and supply the labels as ticks.
    ax2 = Axis(fig[1, 3], title = "Final cohort disposition", ylabel = "Fraction of cohort")
    st = [trajs.settlement_status[p] for p in 1:n_p]
    settled = count(==(:settled_successful), st) / n_p
    dead = count(==(:dead), st) / n_p
    pelagic = max(0.0, 1.0 - settled - dead)
    vals = [settled, pelagic, dead]
    barplot!(ax2, 1:3, vals, color = [:crimson, :steelblue, :grey40], width = 0.6)
    ax2.xticks = (1:3, ["settled", "pelagic", "dead"])

    # Panel 3: survival at each stage transition
    ss = hasproperty(trajs, :stage_survival) ? collect(Float64, trajs.stage_survival) : Float64[]
    has_ss = !isempty(ss) && all(x -> !isnan(x), ss)
    ax3 = Axis(fig[1, 4],
               title = has_ss ? "Survival at stage transitions" : "No stage transitions recorded",
               ylabel = "Mean survival at transition", limits = (0, 4, 0, 1.05))
    has_ss && barplot!(ax3, 1:length(ss), ss, color = :darkorange, width = 0.6)

    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end
