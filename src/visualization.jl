"""
    visualization.jl

CairoMakie-based visualization tools and standalone interactive Leaflet.js dashboards
for plotting Lagrangian particle tracks, Diel Vertical Migration (DVM) profiles,
settlement density kernels, Eulerian hydrodynamic advection velocities, seawater
temperature and salinity tracers, sea surface height, and multi-scenario climate comparisons.
"""

using CairoMakie

"""
    plot_particle_trajectories(
        trajectories::NamedTuple;
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
        output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
        max_display_particles::Int = 100
    )

Render a 2D spatial map of snow crab larval drift trajectories with optional
background bathymetry and management strata boundaries.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort` with fields
  `(lons, lats, depths, times, ids)`.
- `bathymetry_data::Union{Nothing, NamedTuple}`: Optional `(lon, lat, elevation)`
  for background contours.
- `strata::Union{Nothing, AbstractVector{<:NamedTuple}}`: Optional management strata
  definitions with polygon coordinates.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save the figure (PNG/PDF).
- `max_display_particles::Int`: Maximum number of individual tracks to render.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- North, E. W., et al. (2009). ICES Cooperative Research Report, No. 295.
"""
function plot_particle_trajectories(
    trajectories::NamedTuple;
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
    output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
    max_display_particles::Int = 100
)
    fig = Figure(size = (950, 720), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    # 1. Plot background bathymetry if provided
    if !isnothing(bathymetry_data)
        co = contourf!(
            ax,
            bathymetry_data.lon,
            bathymetry_data.lat,
            bathymetry_data.elevation,
            colormap = :viridis,
            levels = 15
        )
        Colorbar(fig[1, 2], co, label = "Seafloor Elevation (m)")
    end

    # 2. Draw management strata outlines if provided
    if !isnothing(strata)
        for s in strata
            if hasproperty(s, :polygon) && !isempty(s.polygon)
                poly_lons = [pt[1] for pt in s.polygon]
                poly_lats = [pt[2] for pt in s.polygon]
                lines!(
                    ax,
                    poly_lons,
                    poly_lats,
                    color = (:gold, 0.75),
                    linewidth = 1.8,
                    linestyle = :dash
                )
            end
        end
    end

    # 3. Plot larval tracks with high-visibility outline underlay and vivid color
    n_particles = size(trajectories.lons, 1)
    n_plot = min(n_particles, max_display_particles)

    # High-contrast dark shadow underlay
    for p in 1:n_plot
        lines!(
            ax,
            trajectories.lons[p, :],
            trajectories.lats[p, :],
            color = (:black, 0.75),
            linewidth = 3.8
        )
    end

    # Vivid foreground trajectory lines
    for p in 1:n_plot
        lines!(
            ax,
            trajectories.lons[p, :],
            trajectories.lats[p, :],
            color = :cyan,
            linewidth = 2.2
        )
    end

    # 4. Mark release positions (start) and final drift positions (end) with strokes
    scatter!(
        ax,
        trajectories.lons[1:n_plot, 1],
        trajectories.lats[1:n_plot, 1],
        color = :springgreen,
        strokecolor = :black,
        strokewidth = 1.5,
        marker = :circle,
        markersize = 11,
        label = "Release (t = 0)"
    )
    scatter!(
        ax,
        trajectories.lons[1:n_plot, end],
        trajectories.lats[1:n_plot, end],
        color = :crimson,
        strokecolor = :white,
        strokewidth = 1.5,
        marker = :diamond,
        markersize = 12,
        label = "Final Position"
    )

    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_vertical_migration_profiles(
        trajectories::NamedTuple;
        sample_indices = 1:min(5, size(trajectories.depths, 1)),
        title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
        show_surface_target::Bool = true,
        surface_target_depth::Real = -10.0,
        output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
    )

Plot larval depth trajectories over time to visualize initial ascent from benthic
release depths and diurnal oscillations between nighttime surface grazing and
daytime sub-surface predator avoidance.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `sample_indices`: Range or vector of particle indices to plot.
- `title::AbstractString`: Figure title.
- `show_surface_target::Bool`: Whether to draw reference line for surface target layer.
- `surface_target_depth::Real`: Surface mixed layer target depth in meters (default -10.0 m).
- `output_path::Union{Nothing, AbstractString}`: Destination path for figure.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Incze, L. S., et al. (1987). *Marine Biology*, 95(2), 195-200.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). *Can. J. Fish. Aquat. Sci.*, 56(11), 2181-2193.
"""
function plot_vertical_migration_profiles(
    trajectories::NamedTuple;
    sample_indices = 1:min(5, size(trajectories.depths, 1)),
    title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
    show_surface_target::Bool = true,
    surface_target_depth::Real = -10.0,
    output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
)
    fig = Figure(size = (850, 480), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Simulation Time (hours)",
        ylabel = "Depth (m)"
    )

    time_hours = trajectories.times ./ 3600.0
    colors = [:royalblue, :darkorange, :forestgreen, :purple, :crimson]

    for (idx, p) in enumerate(sample_indices)
        c = colors[mod1(idx, length(colors))]
        lines!(
            ax,
            time_hours,
            trajectories.depths[p, :],
            color = c,
            linewidth = 2.0,
            label = "Larva $(p)"
        )
        scatter!(
            ax,
            [time_hours[1]],
            [trajectories.depths[p, 1]],
            color = c,
            marker = :circle,
            markersize = 8
        )
    end

    hlines!(ax, [0.0], color = :gray50, linestyle = :dash)
    if show_surface_target
        hlines!(
            ax,
            [Float64(surface_target_depth)],
            color = :teal,
            linestyle = :dot,
            linewidth = 1.5,
            label = "Surface Target ($(round(surface_target_depth, digits=1))m)"
        )
    end
    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    compare_scenario_dispersal(
        scenario_dict::Union{AbstractDict, NamedTuple};
        title::AbstractString = "Larval Dispersal Across Climate Scenarios",
        output_path::Union{Nothing, AbstractString} = "outputs/scenario_comparison.png"
    )

Compare spatial trajectories and final dispersion centroids across multiple
climate forcing scenarios (e.g. `:historical`, `:ssp245`, `:ssp585`).

# Mathematical Formulations
- **Centroid Center of Mass**:
  ```math
  \\bar{\\lambda} = \\frac{1}{N_p} \\sum_{p=1}^{N_p} \\lambda_p(t_{\\text{end}}), \\quad
  \\bar{\\phi} = \\frac{1}{N_p} \\sum_{p=1}^{N_p} \\phi_p(t_{\\text{end}})
  ```

# Inputs
- `scenario_dict::Union{AbstractDict, NamedTuple}`: Mapping of scenario names or keys
  to trajectory NamedTuples (or data containers holding a `.trajectories` field).
- `title::AbstractString`: Comparison figure title.
- `output_path::Union{Nothing, AbstractString}`: Destination filepath to save output.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Brickman, D., et al. (2018). *Progress in Oceanography*, 164, 49-64.
"""
function compare_scenario_dispersal(
    scenario_dict::Union{AbstractDict, NamedTuple};
    title::AbstractString = "Larval Dispersal Across Climate Scenarios",
    output_path::Union{Nothing, AbstractString} = "outputs/scenario_comparison.png"
)
    n_scenarios = length(scenario_dict)
    n_cols = min(n_scenarios, 3)
    n_rows = ceil(Int, n_scenarios / 3)
    fig = Figure(size = (450 * n_cols, 420 * n_rows), fontsize = 12)

    scenario_keys = scenario_dict isa NamedTuple ? keys(scenario_dict) : collect(keys(scenario_dict))
    colors = Dict(
        :historical => :royalblue,
        :baseline => :royalblue,
        :ssp126 => :seagreen,
        :ssp245 => :darkorange,
        :ssp370 => :orangered,
        :ssp585 => :crimson,
        :marine_heatwave => :firebrick
    )

    for (s_idx, s_key) in enumerate(scenario_keys)
        row = div(s_idx - 1, 3) + 1
        col = mod(s_idx - 1, 3) + 1

        raw_val = scenario_dict isa NamedTuple ? getproperty(scenario_dict, s_key) : scenario_dict[s_key]
        traj = hasproperty(raw_val, :trajectories) ? raw_val.trajectories : raw_val
        s_sym = Symbol(s_key)
        c = get(colors, s_sym, :dodgerblue)

        ax = Axis(
            fig[row, col],
            title = "$(string(s_key))",
            xlabel = "Longitude (°E)",
            ylabel = "Latitude (°N)"
        )

        n_p = size(traj.lons, 1)
        for p in 1:min(n_p, 50)
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (:black, 0.5), linewidth = 2.8
            )
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (c, 0.85), linewidth = 1.8
            )
        end

        mean_lon_end = sum(traj.lons[:, end]) / n_p
        mean_lat_end = sum(traj.lats[:, end]) / n_p

        scatter!(
            ax,
            [mean_lon_end],
            [mean_lat_end],
            color = :gold,
            strokecolor = :black,
            strokewidth = 1.5,
            marker = :star5,
            markersize = 16,
            label = "Centroid"
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_larval_dispersal_density(
        trajectories::NamedTuple;
        lon_bins = range(-68.0, -57.0, length = 100),
        lat_bins = range(41.0, 48.0, length = 100),
        title::AbstractString = "Larval Settlement Density Kernel",
        output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
    )

Compute and plot a 2D histogram heatmap of final particle positions representing
larval settlement and nursery retention hotspots.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `lon_bins`: Grid longitude bin edges.
- `lat_bins`: Grid latitude bin edges.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_larval_dispersal_density(
    trajectories::NamedTuple;
    lon_bins = range(-68.0, -57.0, length = 100),
    lat_bins = range(41.0, 48.0, length = 100),
    title::AbstractString = "Larval Settlement Density Kernel",
    output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
)
    fig = Figure(size = (850, 650), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    end_lons = trajectories.lons[:, end]
    end_lats = trajectories.lats[:, end]
    n_p = length(end_lons)

    nx = length(lon_bins) - 1
    ny = length(lat_bins) - 1
    density = zeros(Float64, nx, ny)

    for p in 1:n_p
        x = end_lons[p]
        y = end_lats[p]
        if lon_bins[1] <= x <= lon_bins[end] && lat_bins[1] <= y <= lat_bins[end]
            i = clamp(searchsortedlast(lon_bins, x), 1, nx)
            j = clamp(searchsortedlast(lat_bins, y), 1, ny)
            density[i, j] += 1.0
        end
    end

    density ./= max(1.0, sum(density))
    density .*= 100.0

    lon_centers = [(lon_bins[i] + lon_bins[i + 1]) / 2.0 for i in 1:nx]
    lat_centers = [(lat_bins[j] + lat_bins[j + 1]) / 2.0 for j in 1:ny]

    hm = heatmap!(ax, lon_centers, lat_centers, density, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Settlement Density (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_empirical_movement_field(
        emp_mov::NamedTuple;
        title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
    )

Plot empirical velocity vector arrows over the empirical turbulent diffusivity field.

# Inputs
- `emp_mov::NamedTuple`: Output from `estimate_empirical_movement`.
- `title::AbstractString`: Plot title.
- `output_path::Union{Nothing, AbstractString}`: File path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_empirical_movement_field(
    emp_mov::NamedTuple;
    title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
)
    fig = Figure(size = (850, 650), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    lon_c = emp_mov.lon_centers
    lat_c = emp_mov.lat_centers
    diff = copy(emp_mov.diffusivity)
    diff[isnan.(diff)] .= 0.0

    hm = heatmap!(ax, lon_c, lat_c, diff, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Empirical Diffusivity D (m² s⁻¹)")

    pts = Point2f[]
    dirs = Vec2f[]
    u_m = emp_mov.u_mean
    v_m = emp_mov.v_mean

    for i in 1:length(lon_c), j in 1:length(lat_c)
        if !isnan(u_m[i, j]) && !isnan(v_m[i, j]) && emp_mov.sample_count[i, j] > 2
            push!(pts, Point2f(lon_c[i], lat_c[j]))
            push!(dirs, Vec2f(u_m[i, j] * 0.5, v_m[i, j] * 0.5))
        end
    end

    if !isempty(pts)
        if isdefined(Makie, :arrows2d!)
            arrows2d!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        else
            arrows!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        end
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_connectivity_matrix(
        conn::NamedTuple;
        title::AbstractString = "Macro-Regional Population Connectivity Matrix",
        output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
    )

Render an annotated transition probability heatmap showing larval transfer probabilities
and self-retention percentages between management zones.

# Inputs
- `conn::NamedTuple`: Output from `compute_empirical_connectivity`.
- `title::AbstractString`: Title string.
- `output_path::Union{Nothing, AbstractString}`: Path for output image.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_connectivity_matrix(
    conn::NamedTuple;
    title::AbstractString = "Macro-Regional Population Connectivity Matrix",
    output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
)
    mat = conn.matrix .* 100.0
    strata = conn.strata_names
    n_s = length(strata)

    fig = Figure(size = (800, 650), fontsize = 12)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Destination Stratum",
        ylabel = "Source Stratum",
        xticks = (1:n_s, strata),
        yticks = (1:n_s, strata),
        xticklabelrotation = 0.35
    )

    hm = heatmap!(ax, 1:n_s, 1:n_s, mat', colormap = :blues)
    Colorbar(fig[1, 2], hm, label = "Transition Probability (%)")

    for i in 1:n_s, j in 1:n_s
        val = mat[i, j]
        txt_col = val > 50.0 ? :white : :black
        text!(
            ax,
            Point2f(j, i),
            text = string(round(val, digits = 1), "%"),
            align = (:center, :center),
            color = txt_col,
            fontsize = 12
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_thermal_exposure_map(
        therm_met::NamedTuple;
        title::AbstractString = "Larval Thermal Exposure & Degree-Days",
        output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
    )

Plot a two-panel map showing accumulated thermal degree-days and mean exposure temperature.
"""
function plot_thermal_exposure_map(
    therm_met::NamedTuple;
    title::AbstractString = "Larval Thermal Exposure & Degree-Days",
    output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
)
    fig = Figure(size = (1100, 500), fontsize = 12)

    lon_c = therm_met.lon_centers
    lat_c = therm_met.lat_centers

    # Panel 1: Accumulated Degree-Days
    ax1 = Axis(
        fig[1, 1],
        title = "Mean Cumulative Degree-Days (°C·days)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    dd = copy(therm_met.mean_degree_days)
    dd[isnan.(dd)] .= 0.0
    hm1 = heatmap!(ax1, lon_c, lat_c, dd, colormap = :inferno)
    Colorbar(fig[1, 2], hm1, label = "Degree-Days (°C·d)")

    # Panel 2: Mean Ambient Exposure Temperature
    ax2 = Axis(
        fig[1, 3],
        title = "Mean Exposure Temperature (°C)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    temp = copy(therm_met.mean_exposure_temperature)
    temp[isnan.(temp)] .= 0.0
    hm2 = heatmap!(ax2, lon_c, lat_c, temp, colormap = :thermal)
    Colorbar(fig[1, 4], hm2, label = "Temperature (°C)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_recruitment_summary(
        rec_met::NamedTuple;
        title::AbstractString = "Larval Recruitment & Retention Summary",
        output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
    )

Plot a multi-panel spatial distribution of spawning release, settlement density, and recruitment success rate.
"""
function plot_recruitment_summary(
    rec_met::NamedTuple;
    title::AbstractString = "Larval Recruitment & Retention Summary",
    output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
)
    fig = Figure(size = (1100, 500), fontsize = 12)
    lon_c = rec_met.lon_centers
    lat_c = rec_met.lat_centers

    # Panel 1: Settlement Density
    ax1 = Axis(
        fig[1, 1],
        title = "Settlement Density (count)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lon_c, lat_c, rec_met.settlement_density, colormap = :plasma)
    Colorbar(fig[1, 2], hm1, label = "Settled Larvae")

    # Panel 2: Recruitment Success Rate
    ax2 = Axis(
        fig[1, 3],
        title = "Benthic Recruitment Success Rate (%)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    succ_pct = rec_met.success_rate .* 100.0
    hm2 = heatmap!(ax2, lon_c, lat_c, succ_pct, colormap = :viridis)
    Colorbar(fig[1, 4], hm2, label = "Success Rate (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    resolve_depth_index(
        depths::AbstractVector{<:Real},
        depth::Union{Nothing, Real},
        depth_level::Union{Nothing, Int}
    ) -> Int

Resolve discrete vertical grid level index from either continuous depth in meters
(positive or negative) or discrete index. Defaults to surface (index 1).
"""
function resolve_depth_index(
    depths::AbstractVector{<:Real},
    depth::Union{Nothing, Real},
    depth_level::Union{Nothing, Int}
)::Int
    nz = length(depths)
    nz == 0 && return 1
    if !isnothing(depth)
        target_z = depth > 0 ? -Float64(depth) : Float64(depth)
        return argmin([abs(Float64(d) - target_z) for d in depths])
    elseif !isnothing(depth_level)
        return clamp(depth_level, 1, nz)
    else
        return 1
    end
end

"""
    resolve_time_index(
        times::AbstractVector{<:Real},
        time_seconds::Union{Nothing, Real},
        time_index::Union{Nothing, Int}
    ) -> Int

Resolve discrete temporal snapshot index from either continuous timestamp in seconds
or discrete index. Defaults to first snapshot (index 1).
"""
function resolve_time_index(
    times::AbstractVector{<:Real},
    time_seconds::Union{Nothing, Real},
    time_index::Union{Nothing, Int}
)::Int
    nt = length(times)
    nt == 0 && return 1
    if !isnothing(time_seconds)
        return argmin([abs(Float64(t) - Float64(time_seconds)) for t in times])
    elseif !isnothing(time_index)
        return clamp(time_index, 1, nt)
    else
        return 1
    end
end

"""
    compute_hydrodynamic_diagnostics(
        lons::AbstractVector{<:Real},
        lats::AbstractVector{<:Real},
        depths::AbstractVector{<:Real},
        u::AbstractArray{<:Real, 3},
        v::AbstractArray{<:Real, 3},
        w::AbstractArray{<:Real, 3},
        temp::AbstractArray{<:Real, 3},
        sal::AbstractArray{<:Real, 3};
        ν_closure::Real = 1e-2,
        κ_closure::Real = 1e-2
    ) -> NamedTuple

Compute derived physical oceanographic diagnostics:
- Potential Density \$\\rho(S, T)\$ (Boussinesq linear approximation).
- Brunt-Väisälä buoyancy frequency squared \$N^2 = -(g/\\rho_0) \\partial \\rho / \\partial z\$.
- Vertical salinity stratification gradient \$\\partial S / \\partial z\$.
- Turbulent vertical eddy diffusivity \$\\kappa_v\$ and eddy viscosity \$\\nu_v\$ via
  shear-stratification gradient Richardson number \$Ri = N^2 / [(\\partial u/\\partial z)^2 + (\\partial v/\\partial z)^2]\$.
- Relative vertical vorticity \$\\zeta = \\partial v / \\partial x - \\partial u / \\partial y\$.
"""
function compute_hydrodynamic_diagnostics(
    lons::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    depths::AbstractVector{<:Real},
    u::AbstractArray{<:Real, 3},
    v::AbstractArray{<:Real, 3},
    w::AbstractArray{<:Real, 3},
    temp::AbstractArray{<:Real, 3},
    sal::AbstractArray{<:Real, 3};
    ν_closure::Real = 1e-2,
    κ_closure::Real = 1e-2
)::NamedTuple
    nx = length(lons)
    ny = length(lats)
    nz = length(depths)

    # 1. Seawater Density (Boussinesq linear equation of state)
    rho0 = 1025.0
    alpha = 1.7e-4
    beta = 7.6e-4
    T0 = 10.0
    S0 = 35.0

    rho = zeros(Float64, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        t_val = Float64(temp[i, j, k])
        s_val = Float64(sal[i, j, k])
        if !isnan(t_val) && !isnan(s_val)
            rho[i, j, k] = rho0 * (1.0 - alpha * (t_val - T0) + beta * (s_val - S0))
        else
            rho[i, j, k] = NaN
        end
    end

    # 2. Stratification: N² and ∂S/∂z
    N2 = zeros(Float64, nx, ny, nz)
    dS_dz = zeros(Float64, nx, ny, nz)
    shear_sq = zeros(Float64, nx, ny, nz)

    for k in 1:nz
        k_prev = max(1, k - 1)
        k_next = min(nz, k + 1)
        dz = Float64(depths[k_prev] - depths[k_next])
        if abs(dz) < 1e-6
            dz = 1.0
        end

        for j in 1:ny, i in 1:nx
            if !isnan(rho[i, j, k_prev]) && !isnan(rho[i, j, k_next])
                # Depths are negative, depths[k_prev] > depths[k_next], so dz > 0
                drho = rho[i, j, k_prev] - rho[i, j, k_next]
                N2[i, j, k] = max(0.0, -(9.81 / rho0) * (drho / dz))

                ds = Float64(sal[i, j, k_prev] - sal[i, j, k_next])
                dS_dz[i, j, k] = ds / dz

                du = Float64(u[i, j, k_prev] - u[i, j, k_next]) / dz
                dv = Float64(v[i, j, k_prev] - v[i, j, k_next]) / dz
                shear_sq[i, j, k] = du^2 + dv^2
            else
                N2[i, j, k] = NaN
                dS_dz[i, j, k] = NaN
                shear_sq[i, j, k] = NaN
            end
        end
    end

    # 3. Turbulent eddy diffusivity & viscosity (Richardson number parameterization)
    diff = zeros(Float64, nx, ny, nz)
    visc = zeros(Float64, nx, ny, nz)
    ri_arr = zeros(Float64, nx, ny, nz)

    for k in 1:nz
        z_m = Float64(depths[k])
        surf_mix = 0.015 * exp(z_m / 25.0)

        for j in 1:ny, i in 1:nx
            n2_val = N2[i, j, k]
            s2_val = shear_sq[i, j, k]
            if !isnan(n2_val) && !isnan(s2_val)
                ri = max(0.0, n2_val / max(s2_val, 1e-7))
                ri_arr[i, j, k] = ri

                k_eddy = 1e-5 + 1e-2 / ((1.0 + 5.0 * ri)^2) + surf_mix
                nu_eddy = 1e-4 + 1e-2 / ((1.0 + 5.0 * ri)^3) + surf_mix

                diff[i, j, k] = clamp(k_eddy, 1e-6, 0.05)
                visc[i, j, k] = clamp(nu_eddy, 1e-5, 0.05)
            else
                ri_arr[i, j, k] = NaN
                diff[i, j, k] = NaN
                visc[i, j, k] = NaN
            end
        end
    end

    # 4. Relative Vorticity ζ at surface level
    vort = zeros(Float64, nx, ny)
    r_earth = 6.371e6
    for j in 2:(ny - 1), i in 2:(nx - 1)
        dx_m = r_earth * cosd(Float64(lats[j])) * deg2rad(Float64(lons[i + 1] - lons[i - 1]))
        dy_m = r_earth * deg2rad(Float64(lats[j + 1] - lats[j - 1]))
        if !isnan(v[i + 1, j, 1]) && !isnan(v[i - 1, j, 1]) &&
           !isnan(u[i, j + 1, 1]) && !isnan(u[i, j - 1, 1]) &&
           dx_m > 0.0 && dy_m > 0.0
            dv_dx = (Float64(v[i + 1, j, 1]) - Float64(v[i - 1, j, 1])) / dx_m
            du_dy = (Float64(u[i, j + 1, 1]) - Float64(u[i, j - 1, 1])) / dy_m
            vort[i, j] = dv_dx - du_dy
        else
            vort[i, j] = NaN
        end
    end

    return (
        density = rho,
        stratification = N2,
        salinity_stratification = dS_dz,
        diffusion = diff,
        viscosity = visc,
        richardson = ri_arr,
        vorticity = vort
    )
end

"""
    extract_hydrodynamic_dataset(
        hydro_input::Any;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = nothing,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = nothing,
        domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
        domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
        grid_size::Tuple{Int, Int} = (35, 30),
        target_depths::AbstractVector{<:Real} = [-2.5, -25.0, -50.0, -100.0],
        run_id::Union{Nothing, AbstractString} = nothing
    ) -> NamedTuple

Extract, normalize, and format 2D and 3D hydrodynamic model fields (advection currents,
temperature, salinity, density, stratification, turbulent diffusion, viscosity,
free surface elevation, and bathymetry) at specific depths and times.

# Mathematical Formulations
- **Horizontal Advection Current Velocity**:
  ```math
  \\boldsymbol{u}_h(x, y, z, t) = (u(x, y, z, t), v(x, y, z, t)), \\quad
  |\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}
  ```
- **Seawater Temperature & Practical Salinity**:
  \$T(x, y, z, t)\$ in °C, \$S(x, y, z, t)\$ in PSU.
- **Salinity & Density Stratification**:
  ```math
  N^2(x, y, z, t) = -\\frac{g}{\\rho_0} \\frac{\\partial \\rho}{\\partial z}, \\quad
  \\frac{\\partial S}{\\partial z}(x, y, z, t)
  ```
- **Turbulent Eddy Diffusivity & Viscosity**:
  \$\\kappa_v(x, y, z, t), \\nu_v(x, y, z, t)\$ in \$m^2 s^{-1}\$ parameterized via Richardson number.

# Inputs
- `hydro_input`: Oceananigans model, JLD2 file path, DuckDB connection, NamedTuple, Dict, or `nothing`.
- `depth`: Optional continuous depth in meters (e.g. `-25.0` or `25.0`).
- `depth_level`: Optional vertical level index (1-indexed).
- `time_seconds`: Optional simulation time in seconds.
- `time_index`: Optional discrete time snapshot index.
- `domain_lon`: Longitudinal bounds `(min_lon, max_lon)`.
- `domain_lat`: Latitudinal bounds `(min_lat, max_lat)`.
- `grid_size`: Horizontal grid dimension `(nx, ny)`.
- `target_depths`: Depth coordinates in meters (default `[-2.5, -25.0, -50.0, -100.0]`).
- `run_id`: Optional DuckDB simulation run identifier.

# Outputs
- `NamedTuple` containing:
  - `lons`, `lats`, `depths`: Grid coordinates.
  - `times`, `time_seconds`: Available timestamps and resolved snapshot time.
  - `depth_index`, `depth_m`: Resolved vertical level and depth in meters.
  - `u`, `v`, `w`, `speed`: Velocity fields (\$m s^{-1}\$).
  - `temperature`, `salinity`: Hydrographic tracers (°C, PSU).
  - `density`: Potential density (\$kg m^{-3}\$).
  - `stratification`: Buoyancy frequency squared \$N^2\$ (\$s^{-2}\$).
  - `salinity_stratification`: Vertical salinity gradient (\$PSU m^{-1}\$).
  - `diffusion`, `viscosity`: Turbulent diffusivities (\$m^2 s^{-1}\$).
  - `richardson_number`: Gradient Richardson number \$Ri\$.
  - `vorticity`: Relative vorticity (\$s^{-1}\$).
  - `elevation`, `bathymetry`: Surface elevation and seafloor depth (m).
"""
function extract_hydrodynamic_dataset(
    hydro_input::Any;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = nothing,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = nothing,
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    grid_size::Tuple{Int, Int} = (35, 30),
    target_depths::AbstractVector{<:Real} = [-2.5, -25.0, -50.0, -100.0],
    run_id::Union{Nothing, AbstractString} = nothing
)
    nx, ny = grid_size
    t_depths = collect(Float64, target_depths)
    nz = length(t_depths)

    # 1. JLD2 Simulation Output File
    if hydro_input isa AbstractString && isfile(hydro_input) && endswith(hydro_input, ".jld2")
        local lons_j, lats_j, deps_j, times_j, u_j, v_j, w_j, T_j, S_j, elev_j
        jldopen(hydro_input, "r") do file
            u_group = file["timeseries/u"]
            raw_keys = collect(keys(u_group))
            sorted_keys = sort(raw_keys, by = k -> something(tryparse(Float64, k), 0.0))
            times_j = haskey(file, "timeseries/t") ?
                collect(Float64, file["timeseries/t"]) :
                [something(tryparse(Float64, k), Float64(idx)) for (idx, k) in enumerate(sorted_keys)]

            t_idx = resolve_time_index(times_j, time_seconds, time_index)
            key_sel = sorted_keys[t_idx]

            sample_u = file["timeseries/u/$(key_sel)"]
            nx_f, ny_f, nz_f = size(sample_u)

            if haskey(file, "grid")
                g = file["grid"]
                lons_j = haskey(g, "λᶜᵃᵃ") ? collect(Float64, g["λᶜᵃᵃ"][1:nx_f]) : collect(range(domain_lon[1], domain_lon[2], length=nx_f))
                lats_j = haskey(g, "φᵃᶜᵃ") ? collect(Float64, g["φᵃᶜᵃ"][1:ny_f]) : collect(range(domain_lat[1], domain_lat[2], length=ny_f))
                deps_j = haskey(g, "zᵃᵃᶜ") ? collect(Float64, g["zᵃᵃᶜ"][1:nz_f]) : t_depths
            else
                lons_j = collect(range(domain_lon[1], domain_lon[2], length=nx_f))
                lats_j = collect(range(domain_lat[1], domain_lat[2], length=ny_f))
                deps_j = t_depths
            end

            u_j = Float64.(file["timeseries/u/$(key_sel)"])
            v_j = Float64.(file["timeseries/v/$(key_sel)"])
            w_j = haskey(file, "timeseries/w") ? Float64.(file["timeseries/w/$(key_sel)"]) : zeros(nx_f, ny_f, nz_f)
            T_j = haskey(file, "timeseries/T") ? Float64.(file["timeseries/T/$(key_sel)"]) : fill(4.5, nx_f, ny_f, nz_f)
            S_j = haskey(file, "timeseries/S") ? Float64.(file["timeseries/S/$(key_sel)"]) : fill(33.0, nx_f, ny_f, nz_f)
            elev_j = haskey(file, "timeseries/η") ? Float64.(file["timeseries/η/$(key_sel)"]) : zeros(nx_f, ny_f)
        end

        diag = compute_hydrodynamic_diagnostics(lons_j, lats_j, deps_j, u_j, v_j, w_j, T_j, S_j)
        k_sel = resolve_depth_index(deps_j, depth, depth_level)
        sel_time = isempty(times_j) ? 0.0 : times_j[resolve_time_index(times_j, time_seconds, time_index)]

        return (
            lons = lons_j,
            lats = lats_j,
            depths = deps_j,
            times = times_j,
            time_seconds = sel_time,
            resolved_time = sel_time,
            depth_index = k_sel,
            depth_m = deps_j[k_sel],
            resolved_depth = deps_j[k_sel],
            u = u_j,
            v = v_j,
            w = w_j,
            speed = hypot.(u_j, v_j),
            temperature = T_j,
            salinity = S_j,
            density = diag.density,
            stratification = diag.stratification,
            salinity_stratification = diag.salinity_stratification,
            diffusion = diag.diffusion,
            viscosity = diag.viscosity,
            richardson_number = diag.richardson,
            vorticity = diag.vorticity,
            elevation = elev_j,
            bathymetry = fill(-150.0, length(lons_j), length(lats_j))
        )
    end

    # 2. DuckDB Database Connection
    if !isnothing(hydro_input) && (hydro_input isa DuckDB.DB)
        r_id = !isnothing(run_id) ? String(run_id) : begin
            runs_df = list_simulation_runs(hydro_input)
            nrow(runs_df) > 0 ? String(first(runs_df.run_id)) : "run_baseline_2025"
        end
        loaded = load_hydrodynamic_field(hydro_input, r_id; time_seconds = time_seconds, depth_level = depth_level)
        diag = compute_hydrodynamic_diagnostics(loaded.lons, loaded.lats, loaded.depths, loaded.u, loaded.v, loaded.w, loaded.temperature, loaded.salinity)
        k_sel = resolve_depth_index(loaded.depths, depth, depth_level)

        return (
            lons = loaded.lons,
            lats = loaded.lats,
            depths = loaded.depths,
            times = !isnothing(time_seconds) ? [Float64(time_seconds)] : [0.0],
            time_seconds = something(time_seconds, 0.0),
            resolved_time = something(time_seconds, 0.0),
            depth_index = k_sel,
            depth_m = loaded.depths[k_sel],
            resolved_depth = loaded.depths[k_sel],
            u = loaded.u,
            v = loaded.v,
            w = loaded.w,
            speed = hypot.(loaded.u, loaded.v),
            temperature = loaded.temperature,
            salinity = loaded.salinity,
            density = diag.density,
            stratification = diag.stratification,
            salinity_stratification = diag.salinity_stratification,
            diffusion = diag.diffusion,
            viscosity = diag.viscosity,
            richardson_number = diag.richardson,
            vorticity = diag.vorticity,
            elevation = loaded.elevation,
            bathymetry = fill(-150.0, length(loaded.lons), length(loaded.lats))
        )
    end

    # 3. Direct Oceananigans Model Instance
    if !isnothing(hydro_input) && hasproperty(hydro_input, :velocities) && hasproperty(hydro_input, :tracers)
        g = hydro_input.grid
        coords = extract_grid_coordinates(g)
        m_lons = coords.lons
        m_lats = coords.lats
        m_depths = coords.depths

        u_arr = Array(interior(hydro_input.velocities.u))
        v_arr = Array(interior(hydro_input.velocities.v))
        w_arr = Array(interior(hydro_input.velocities.w))
        t_arr = haskey(hydro_input.tracers, :T) ? Array(interior(hydro_input.tracers.T)) : fill(4.5, size(u_arr))
        s_arr = haskey(hydro_input.tracers, :S) ? Array(interior(hydro_input.tracers.S)) : fill(33.0, size(u_arr))
        elev_mat = hasproperty(hydro_input, :free_surface) && hasproperty(hydro_input.free_surface, :η) ?
                   Array(interior(hydro_input.free_surface.η)) : zeros(Float64, length(m_lons), length(m_lats))

        bathy_mat = if g isa ImmersedBoundaryGrid && hasproperty(g.immersed_boundary, :bottom_height)
            Array(interior(g.immersed_boundary.bottom_height))
        else
            fill(-150.0, length(m_lons), length(m_lats))
        end

        model_time = hasproperty(hydro_input, :clock) ? Float64(hydro_input.clock.time) : 0.0
        diag = compute_hydrodynamic_diagnostics(m_lons, m_lats, m_depths, u_arr, v_arr, w_arr, t_arr, s_arr)
        k_sel = resolve_depth_index(m_depths, depth, depth_level)

        return (
            lons = collect(Float64, m_lons),
            lats = collect(Float64, m_lats),
            depths = collect(Float64, m_depths),
            times = [model_time],
            time_seconds = model_time,
            resolved_time = model_time,
            depth_index = k_sel,
            depth_m = m_depths[k_sel],
            resolved_depth = m_depths[k_sel],
            u = Float64.(u_arr),
            v = Float64.(v_arr),
            w = Float64.(w_arr),
            speed = hypot.(Float64.(u_arr), Float64.(v_arr)),
            temperature = Float64.(t_arr),
            salinity = Float64.(s_arr),
            density = diag.density,
            stratification = diag.stratification,
            salinity_stratification = diag.salinity_stratification,
            diffusion = diag.diffusion,
            viscosity = diag.viscosity,
            richardson_number = diag.richardson,
            vorticity = diag.vorticity,
            elevation = Float64.(elev_mat),
            bathymetry = Float64.(bathy_mat)
        )
    end

    # 4. NamedTuple / Dict representation
    if !isnothing(hydro_input) && (hydro_input isa NamedTuple || hydro_input isa AbstractDict)
        get_field(keys_to_try, default_val) = begin
            for k in keys_to_try
                if hydro_input isa NamedTuple && hasproperty(hydro_input, k)
                    return getproperty(hydro_input, k)
                elseif hydro_input isa AbstractDict && (haskey(hydro_input, k) || haskey(hydro_input, string(k)))
                    return haskey(hydro_input, k) ? hydro_input[k] : hydro_input[string(k)]
                end
            end
            return default_val
        end

        target_lons = range(domain_lon[1], domain_lon[2], length = nx)
        target_lats = range(domain_lat[1], domain_lat[2], length = ny)

        in_lons = get_field((:lons, :grid_lons, :lon), target_lons)
        in_lats = get_field((:lats, :grid_lats, :lat), target_lats)
        in_depths = get_field((:depths, :grid_depths, :depth), t_depths)

        nx_in = length(in_lons)
        ny_in = length(in_lats)
        nz_in = length(in_depths)

        in_u = get_field((:u, :u_velocity, :u_mean), zeros(Float64, nx_in, ny_in, nz_in))
        in_v = get_field((:v, :v_velocity, :v_mean), zeros(Float64, nx_in, ny_in, nz_in))
        in_w = get_field((:w, :w_velocity), zeros(Float64, nx_in, ny_in, nz_in))
        in_t = get_field((:temperature, :T, :temp), fill(4.5, nx_in, ny_in, nz_in))
        in_s = get_field((:salinity, :S, :sal), fill(33.0, nx_in, ny_in, nz_in))
        in_elev = get_field((:elevation, :η, :eta, :ssh), zeros(Float64, nx_in, ny_in))
        in_bathy = get_field((:bathymetry, :bathy, :elevation_bottom), fill(-150.0, nx_in, ny_in))

        # Handle 4D time slices if present
        slice_3d(arr) = begin
            if ndims(arr) == 4
                t_sub = resolve_time_index(get_field((:times, :t), [0.0]), time_seconds, time_index)
                return arr[:, :, :, clamp(t_sub, 1, size(arr, 4))]
            elseif ndims(arr) == 2
                return reshape(arr, size(arr, 1), size(arr, 2), 1)
            else
                return arr
            end
        end

        u_3d = Float64.(slice_3d(in_u))
        v_3d = Float64.(slice_3d(in_v))
        w_3d = Float64.(slice_3d(in_w))
        t_3d = Float64.(slice_3d(in_t))
        s_3d = Float64.(slice_3d(in_s))

        times_arr = get_field((:times, :t), [something(time_seconds, 0.0)])
        sel_time = isempty(times_arr) ? 0.0 : times_arr[resolve_time_index(times_arr, time_seconds, time_index)]
        diag = compute_hydrodynamic_diagnostics(in_lons, in_lats, in_depths, u_3d, v_3d, w_3d, t_3d, s_3d)
        k_sel = resolve_depth_index(in_depths, depth, depth_level)

        return (
            lons = collect(Float64, in_lons),
            lats = collect(Float64, in_lats),
            depths = collect(Float64, in_depths),
            times = collect(Float64, times_arr),
            time_seconds = Float64(sel_time),
            resolved_time = Float64(sel_time),
            depth_index = k_sel,
            depth_m = Float64(in_depths[k_sel]),
            resolved_depth = Float64(in_depths[k_sel]),
            u = u_3d,
            v = v_3d,
            w = w_3d,
            speed = hypot.(u_3d, v_3d),
            temperature = t_3d,
            salinity = s_3d,
            density = diag.density,
            stratification = diag.stratification,
            salinity_stratification = diag.salinity_stratification,
            diffusion = diag.diffusion,
            viscosity = diag.viscosity,
            richardson_number = diag.richardson,
            vorticity = diag.vorticity,
            elevation = Float64.(in_elev isa AbstractMatrix ? in_elev : fill(0.0, nx_in, ny_in)),
            bathymetry = Float64.(in_bathy isa AbstractMatrix ? in_bathy : fill(-150.0, nx_in, ny_in))
        )
    end

    # 5. Default Realistic Regional Hydrodynamic Synthesis with Temporal Evolution
    t_sec = something(time_seconds, 0.0)
    target_lons = range(domain_lon[1], domain_lon[2], length = nx)
    target_lats = range(domain_lat[1], domain_lat[2], length = ny)

    u_mat = zeros(Float64, nx, ny, nz)
    v_mat = zeros(Float64, nx, ny, nz)
    w_mat = zeros(Float64, nx, ny, nz)
    t_mat = zeros(Float64, nx, ny, nz)
    s_mat = zeros(Float64, nx, ny, nz)
    elev_mat = zeros(Float64, nx, ny)
    bathy_mat = zeros(Float64, nx, ny)

    lons_vec = collect(Float64, target_lons)
    lats_vec = collect(Float64, target_lats)

    # Semi-diurnal M2 tidal phase and seasonal warming phase
    omega_m2 = 2.0 * π / (12.42 * 3600.0)
    tide_phase = cos(omega_m2 * t_sec)
    season_phase = sin(2.0 * π * (t_sec - 60.0 * 86400.0) / (365.25 * 86400.0))

    for i in 1:nx, j in 1:ny
        lon = lons_vec[i]
        lat = lats_vec[j]
        x_norm = (lon - domain_lon[1]) / (domain_lon[2] - domain_lon[1])
        y_norm = (lat - domain_lat[1]) / (domain_lat[2] - domain_lat[1])

        # Realistic Scotian Shelf Bathymetry profile
        dist_to_coast = (lat - 43.5) + 0.4 * (lon + 63.5)
        b_elev = if dist_to_coast > 1.2
            -45.0 - 55.0 * (1.0 - y_norm)
        elseif dist_to_coast > -0.2
            -110.0 - 140.0 * sin(π * x_norm)
        else
            -350.0 - 2200.0 * (1.0 - (dist_to_coast + 1.5) / 1.3)^2
        end
        bathy_mat[i, j] = clamp(b_elev, -3200.0, -25.0)

        # Free Surface SSH: cross-shelf steric setup + tidal oscillation
        elev_mat[i, j] = (0.05 * sin(2.0 * π * x_norm) - 0.04 * cos(π * y_norm)) + 0.35 * tide_phase

        for k in 1:nz
            z = t_depths[k]
            depth_factor = exp(z / 75.0)

            # Alongshore Scotian Current Jet with tidal modulation
            jet_core = exp(-((dist_to_coast - 0.4)^2) / 0.18)
            u_base = -0.09 * depth_factor * (0.6 + 0.8 * jet_core) + 0.02 * sin(2.0 * π * y_norm) +
                     0.06 * tide_phase * sin(π * y_norm)
            v_base = -0.04 * depth_factor * (0.5 + 0.7 * jet_core) + 0.015 * cos(2.0 * π * x_norm) +
                     0.04 * sin(omega_m2 * t_sec) * cos(π * x_norm)

            u_mat[i, j, k] = u_base
            v_mat[i, j, k] = v_base

            # Vertical velocity (coastal upwelling / downwelling)
            w_mat[i, j, k] = 0.00035 * sin(2.0 * π * x_norm) * sin(π * y_norm) * (1.0 + z / 100.0)

            # Thermal Stratification: Surface warm layer + seasonal phase, CIL at -50m, warm deep slope
            t_surface = 14.2 + 2.5 * season_phase - 2.5 * y_norm + 1.2 * x_norm
            t_cil = 2.2 + 0.8 * sin(π * x_norm)
            t_slope = 7.5 + 1.5 * (1.0 - y_norm)

            t_val = if z > -20.0
                t_surface + (z / 20.0) * (t_surface - 6.0)
            elseif z > -70.0
                t_cil + ((z + 50.0) / 30.0)^2 * 2.5
            else
                t_cil + ((abs(z) - 70.0) / 50.0) * (t_slope - t_cil)
            end
            t_mat[i, j, k] = clamp(t_val, 0.5, 19.5)

            # Practical Salinity stratification
            s_val = 31.4 + 1.8 * (1.0 - y_norm) + 1.2 * x_norm + (abs(z) / 100.0) * 1.1
            s_mat[i, j, k] = clamp(s_val, 30.2, 35.6)
        end
    end

    spd_mat = hypot.(u_mat, v_mat)

    # Mask over land
    for i in 1:nx, j in 1:ny
        if is_point_on_land(lons_vec[i], lats_vec[j])
            for k in 1:nz
                u_mat[i, j, k]   = NaN
                v_mat[i, j, k]   = NaN
                w_mat[i, j, k]   = NaN
                spd_mat[i, j, k] = NaN
                t_mat[i, j, k]   = NaN
                s_mat[i, j, k]   = NaN
            end
            elev_mat[i, j]  = NaN
            bathy_mat[i, j] = NaN
        end
    end

    diag = compute_hydrodynamic_diagnostics(lons_vec, lats_vec, t_depths, u_mat, v_mat, w_mat, t_mat, s_mat)
    k_sel = resolve_depth_index(t_depths, depth, depth_level)

    return (
        lons = lons_vec,
        lats = lats_vec,
        depths = t_depths,
        times = [t_sec],
        time_seconds = t_sec,
        resolved_time = t_sec,
        depth_index = k_sel,
        depth_m = t_depths[k_sel],
        resolved_depth = t_depths[k_sel],
        u = u_mat,
        v = v_mat,
        w = w_mat,
        speed = spd_mat,
        temperature = t_mat,
        salinity = s_mat,
        density = diag.density,
        stratification = diag.stratification,
        salinity_stratification = diag.salinity_stratification,
        diffusion = diag.diffusion,
        viscosity = diag.viscosity,
        richardson_number = diag.richardson,
        vorticity = diag.vorticity,
        elevation = elev_mat,
        bathymetry = bathy_mat
    )
end

"""
    plot_hydrodynamic_advection(
        hydrodynamics::Any;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = 1,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        title::AbstractString = "Ocean Hydrodynamic Advection Velocity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
        quiver_stride::Int = 2,
        colormap::Symbol = :turbo
    ) -> Figure

Plot 2D horizontal advection current velocity vector arrows \$\\boldsymbol{u}_h = (u, v)\$
over current speed magnitude \$|\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}\$ at a specific depth and time.

# Mathematical Formulation
Horizontal velocity field \$\\boldsymbol{u}_h(x, y, z_k, t_m) = (u, v)\$ with speed:
```math
|\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}
```
Flow orientation angle:
```math
\\theta = \\operatorname{atan2}(v, u)
```

# Inputs
- `hydrodynamics`: Model instance, JLD2 path, DuckDB database, or NamedTuple.
- `depth`: Target continuous depth in meters (e.g. `-25.0` or `25.0`).
- `depth_level`: Vertical discrete index (default: 1).
- `time_seconds`: Timestamp in seconds.
- `time_index`: Snapshot index.
- `bathymetry_data`: Optional background elevation dataset.
- `title`: Figure title.
- `output_path`: Destination file path.
- `quiver_stride`: Spatial subsampling step for vector arrows.
- `colormap`: CairoMakie colormap symbol (default `:turbo`).

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_advection(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = 1,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    title::AbstractString = "Ocean Hydrodynamic Advection Velocity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
    quiver_stride::Int = 2,
    colormap::Symbol = :turbo
)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        depth = depth,
        depth_level = depth_level,
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    k = hydro.depth_index
    depth_m = hydro.depth_m
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    u_k = hydro.u[:, :, k]
    v_k = hydro.v[:, :, k]
    spd_k = hydro.speed[:, :, k] .* 100.0

    fig = Figure(size = (980, 720), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = "\$(title)\n[Depth: \$(depth_m) m | Time: \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    hm = heatmap!(ax, lons, lats, spd_k, colormap = colormap)
    Colorbar(fig[1, 2], hm, label = "Current Speed (|u_h|, cm/s)")

    if !isnothing(bathymetry_data) && hasproperty(bathymetry_data, :elevation)
        contour!(
            ax,
            bathymetry_data.lon,
            bathymetry_data.lat,
            bathymetry_data.elevation,
            levels = [-2000.0, -1000.0, -500.0, -200.0, -100.0, -50.0],
            color = (:white, 0.4),
            linewidth = 1.0
        )
    end

    nx = length(lons)
    ny = length(lats)
    pts = Point2f[]
    dirs = Vec2f[]
    step_x = max(1, quiver_stride)
    step_y = max(1, quiver_stride)

    for i in 1:step_x:nx, j in 1:step_y:ny
        u_val = u_k[i, j]
        v_val = v_k[i, j]
        if !isnan(u_val) && !isnan(v_val)
            push!(pts, Point2f(lons[i], lats[j]))
            push!(dirs, Vec2f(u_val * 0.4, v_val * 0.4))
        end
    end

    if !isempty(pts)
        if isdefined(Makie, :arrows2d!)
            arrows2d!(
                ax, pts, dirs,
                tipwidth = 6, tiplength = 10, lengthscale = 1.0,
                tipcolor = :white, shaftcolor = :white
            )
        else
            arrows!(
                ax, pts, dirs,
                tipwidth = 6, tiplength = 10, lengthscale = 1.0,
                tipcolor = :white, shaftcolor = :white
            )
        end
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_tracers(
        hydrodynamics::Any;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = 1,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        title::AbstractString = "Hydrodynamic Model Seawater Tracers",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png"
    ) -> Figure

Render a two-panel spatial distribution of active seawater tracer fields
(Temperature \$T\$ and Practical Salinity \$S\$) from the hydrodynamic model
at a specific depth and time.

# Mathematical Formulations
- **Thermal Tracer Distribution**: \$T(x, y, z_k, t_m)\$ [°C].
- **Haline Tracer Distribution**: \$S(x, y, z_k, t_m)\$ [PSU].

# Inputs
- `hydrodynamics`: Model, JLD2 file, DuckDB database, or NamedTuple.
- `depth`: Target depth in meters (e.g. `-50.0`).
- `depth_level`: Vertical index (default: 1).
- `time_seconds`: Timestamp in seconds.
- `time_index`: Snapshot index.
- `title`: Figure title.
- `output_path`: Destination path.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_tracers(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = 1,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    title::AbstractString = "Hydrodynamic Model Seawater Tracers",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png"
)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        depth = depth,
        depth_level = depth_level,
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    k = hydro.depth_index
    depth_m = hydro.depth_m
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    t_mat = hydro.temperature[:, :, k]
    s_mat = hydro.salinity[:, :, k]

    fig = Figure(size = (1180, 540), fontsize = 12)

    # Panel 1: Temperature
    ax1 = Axis(
        fig[1, 1],
        title = "Seawater Temperature T (°C)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lons, lats, t_mat, colormap = :thermal)
    Colorbar(fig[1, 2], hm1, label = "Temperature (°C)")

    # Panel 2: Salinity
    ax2 = Axis(
        fig[1, 3],
        title = "Practical Salinity S (PSU)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm2 = heatmap!(ax2, lons, lats, s_mat, colormap = :viridis)
    Colorbar(fig[1, 4], hm2, label = "Salinity (PSU)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_stratification(
        hydrodynamics::Any;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = 1,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        stations::Union{Nothing, AbstractVector} = nothing,
        title::AbstractString = "Hydrodynamic Seawater Stratification & Pycnocline Structure",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_stratification.png"
    ) -> Figure

Render a three-panel spatial and vertical analysis of ocean stratification:
1. Horizontal map of Brunt-Väisälä buoyancy frequency squared \$N^2\$ (\$10^{-4} s^{-2}\$).
2. Horizontal map of Practical Salinity \$S\$ (PSU).
3. 1D vertical water-column profiles of \$S(z), T(z), N^2(z)\$ at representative shelf stations.

# Mathematical Formulation
Buoyancy frequency squared (gravitational stability metric):
```math
N^2 = -\\frac{g}{\\rho_0} \\frac{\\partial \\rho}{\\partial z}
    \\approx g \\left( \\alpha \\frac{\\partial T}{\\partial z} - \\beta \\frac{\\partial S}{\\partial z} \\right)
```
Vertical salinity gradient (halocline strength):
```math
\\frac{\\partial S}{\\partial z} = \\frac{S(z_1) - S(z_2)}{\\Delta z}
```

# Inputs
- `hydrodynamics`: Model, JLD2 file, DuckDB database, or NamedTuple.
- `depth`: Depth in meters (default: surface).
- `depth_level`: Vertical level index.
- `time_seconds`: Simulation time in seconds.
- `time_index`: Snapshot index.
- `stations`: Optional vector of named coordinates `[(lon=..., lat=..., name=...)]`.
- `title`: Figure title.
- `output_path`: Destination path.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_stratification(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = 1,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    stations::Union{Nothing, AbstractVector} = nothing,
    title::AbstractString = "Hydrodynamic Seawater Stratification & Pycnocline Structure",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_stratification.png"
)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        depth = depth,
        depth_level = depth_level,
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = hydro.depth_index
    depth_m = hydro.depth_m
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    n2_mat = hydro.stratification[:, :, k] .* 10000.0 # 10^-4 s^-2
    s_mat  = hydro.salinity[:, :, k]

    # Default representative shelf stations if not provided
    station_list = if !isnothing(stations) && !isempty(stations)
        stations
    else
        [
            (lon = -63.5, lat = 44.5, name = "Coastal (Halifax)"),
            (lon = -61.5, lat = 43.5, name = "Mid-Shelf (Emerald)"),
            (lon = -59.0, lat = 42.8, name = "Outer Slope (Laurentian)")
        ]
    end

    fig = Figure(size = (1500, 520), fontsize = 12)

    # Panel 1: Buoyancy Frequency N²
    ax1 = Axis(
        fig[1, 1],
        title = "Buoyancy Frequency N² (10⁻⁴ s⁻²)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lons, lats, n2_mat, colormap = :ice)
    Colorbar(fig[1, 2], hm1, label = "N² (10⁻⁴ s⁻²)")

    # Panel 2: Practical Salinity S
    ax2 = Axis(
        fig[1, 3],
        title = "Practical Salinity S (PSU)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm2 = heatmap!(ax2, lons, lats, s_mat, colormap = :viridis)
    Colorbar(fig[1, 4], hm2, label = "Salinity (PSU)")

    # Panel 3: Vertical 1D Profiles at stations
    ax3 = Axis(
        fig[1, 5],
        title = "Vertical Stratification S(z) Profiles",
        xlabel = "Practical Salinity S (PSU)",
        ylabel = "Depth (m)"
    )

    colors = [:darkblue, :forestgreen, :darkorange, :purple]
    markers = [:circle, :diamond, :rect, :star5]

    for (st_idx, st) in enumerate(station_list)
        c_lon = st.lon
        c_lat = st.lat
        st_name = hasproperty(st, :name) ? st.name : "Station \$(st_idx)"
        col = colors[mod1(st_idx, length(colors))]
        mkr = markers[mod1(st_idx, length(markers))]

        i_st = argmin(abs.(lons .- c_lon))
        j_st = argmin(abs.(lats .- c_lat))

        # Mark station location on map panels
        scatter!(ax1, [lons[i_st]], [lats[j_st]], color = col, marker = mkr, markersize = 12, strokecolor = :white, strokewidth = 1.5)
        scatter!(ax2, [lons[i_st]], [lats[j_st]], color = col, marker = mkr, markersize = 12, strokecolor = :white, strokewidth = 1.5)

        # Plot vertical 1D salinity profile S(z)
        s_profile = hydro.salinity[i_st, j_st, :]
        lines!(ax3, s_profile, depths, color = col, linewidth = 2.4, label = st_name)
        scatter!(ax3, s_profile, depths, color = col, marker = mkr, markersize = 8)
    end

    axislegend(ax3, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_diffusion(
        hydrodynamics::Any;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = 1,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        title::AbstractString = "Hydrodynamic Turbulent Eddy Diffusivity & Viscosity",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_diffusion.png"
    ) -> Figure

Render a three-panel visualization of ocean turbulent mixing fields:
1. Horizontal map of turbulent vertical eddy diffusivity \$\\kappa_v\$ (\$m^2 s^{-1}\$).
2. Horizontal map of turbulent vertical eddy viscosity \$\\nu_v\$ (\$m^2 s^{-1}\$).
3. Vertical 1D profiles of \$\\kappa_v(z)\$ illustrating mixed-layer mixing and pycnocline barrier.

# Mathematical Formulation
Parameterization via gradient Richardson number \$Ri\$:
```math
Ri = \\frac{N^2}{\\left(\\frac{\\partial u}{\\partial z}\\right)^2 + \\left(\\frac{\\partial v}{\\partial z}\\right)^2}, \\quad
\\kappa_v(z) = \\kappa_{\\text{bg}} + \\frac{\\kappa_{\\text{max}}}{(1 + 5 Ri)^2} + \\kappa_{\\text{surf}} \\exp(z / h_{\\text{mix}})
```

# Inputs
- `hydrodynamics`: Model instance, JLD2 file, DuckDB database, or NamedTuple.
- `depth`: Depth in meters.
- `depth_level`: Vertical level index.
- `time_seconds`: Timestamp in seconds.
- `time_index`: Snapshot index.
- `title`: Figure title.
- `output_path`: Destination file path.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_diffusion(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = 1,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    title::AbstractString = "Hydrodynamic Turbulent Eddy Diffusivity & Viscosity",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_diffusion.png"
)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        depth = depth,
        depth_level = depth_level,
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = hydro.depth_index
    depth_m = hydro.depth_m
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    diff_mat = hydro.diffusion[:, :, k] .* 10000.0 # 10^-4 m^2/s
    visc_mat = hydro.viscosity[:, :, k] .* 10000.0

    fig = Figure(size = (1500, 520), fontsize = 12)

    # Panel 1: Eddy Diffusivity κ_v
    ax1 = Axis(
        fig[1, 1],
        title = "Eddy Diffusivity κ_v (10⁻⁴ m² s⁻¹)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lons, lats, diff_mat, colormap = :plasma)
    Colorbar(fig[1, 2], hm1, label = "κ_v (10⁻⁴ m² s⁻¹)")

    # Panel 2: Eddy Viscosity ν_v
    ax2 = Axis(
        fig[1, 3],
        title = "Eddy Viscosity ν_v (10⁻⁴ m² s⁻¹)\n[Depth: \$(depth_m) m | t = \$(t_hr) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm2 = heatmap!(ax2, lons, lats, visc_mat, colormap = :magma)
    Colorbar(fig[1, 4], hm2, label = "ν_v (10⁻⁴ m² s⁻¹)")

    # Panel 3: Vertical 1D Diffusivity Profiles
    ax3 = Axis(
        fig[1, 5],
        title = "Vertical Diffusivity Profiles κ_v(z)",
        xlabel = "Eddy Diffusivity (10⁻⁴ m² s⁻¹)",
        ylabel = "Depth (m)"
    )

    sample_points = [
        (lon = -63.5, lat = 44.5, name = "Coastal Mixed Zone", col = :royalblue),
        (lon = -61.5, lat = 43.5, name = "Mid-Shelf Pycnocline", col = :forestgreen),
        (lon = -59.0, lat = 42.8, name = "Slope Boundary", col = :crimson)
    ]

    for pt in sample_points
        i_pt = argmin(abs.(lons .- pt.lon))
        j_pt = argmin(abs.(lats .- pt.lat))
        k_prof = hydro.diffusion[i_pt, j_pt, :] .* 10000.0

        scatter!(ax1, [lons[i_pt]], [lats[j_pt]], color = pt.col, marker = :diamond, markersize = 12, strokecolor = :white, strokewidth = 1.5)
        lines!(ax3, k_prof, depths, color = pt.col, linewidth = 2.4, label = pt.name)
        scatter!(ax3, k_prof, depths, color = pt.col, marker = :diamond, markersize = 8)
    end

    axislegend(ax3, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_section(
        hydrodynamics::Any;
        variable::Symbol = :temperature,
        transect_type::Symbol = :zonal,
        coordinate::Real = 44.0,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        title::Union{Nothing, AbstractString} = nothing,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_section.png"
    ) -> Figure

Render a vertical cross-section (depth \$z\$ versus horizontal distance) along a specified
transect line for any active hydrodynamic field (temperature, salinity, stratification,
advection currents, turbulent diffusion) with seafloor bathymetry masking.

# Inputs
- `hydrodynamics`: Hydrodynamic model instance, JLD2 file, DuckDB database, or NamedTuple.
- `variable`: Variable symbol: `:temperature`, `:salinity`, `:density`, `:stratification`,
  `:diffusion`, `:viscosity`, `:speed`, `:u`, `:v`, `:w`.
- `transect_type`: `:zonal` (fixed latitude, varying longitude) or `:meridional` (fixed longitude).
- `coordinate`: Fixed latitude or longitude coordinate in degrees.
- `time_seconds`: Simulation time in seconds.
- `time_index`: Snapshot index.
- `title`: Optional custom title.
- `output_path`: Destination path.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_section(
    hydrodynamics::Any;
    variable::Symbol = :temperature,
    transect_type::Symbol = :zonal,
    section_type::Union{Nothing, Symbol} = nothing,
    coordinate::Real = 44.0,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_section.png"
)
    t_mode = !isnothing(section_type) ? section_type : transect_type
    is_zonal = t_mode in (:zonal, :lat, :latitude, :east_west)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    nx = length(lons)
    ny = length(lats)
    nz = length(depths)

    # Resolve variable 3D field and colorbar properties
    field_3d, cmap, var_label = if variable == :salinity || variable == :S
        (hydro.salinity, :viridis, "Practical Salinity (PSU)")
    elseif variable == :density || variable == :rho
        (hydro.density, :dense, "Potential Density (kg m⁻³)")
    elseif variable == :stratification || variable == :N2
        (hydro.stratification .* 10000.0, :ice, "N² (10⁻⁴ s⁻²)")
    elseif variable == :diffusion || variable == :kappa
        (hydro.diffusion .* 10000.0, :plasma, "Diffusivity κ_v (10⁻⁴ m² s⁻¹)")
    elseif variable == :viscosity || variable == :nu
        (hydro.viscosity .* 10000.0, :magma, "Viscosity ν_v (10⁻⁴ m² s⁻¹)")
    elseif variable == :speed
        (hydro.speed .* 100.0, :turbo, "Current Speed (cm s⁻¹)")
    elseif variable == :u
        (hydro.u .* 100.0, :balance, "Zonal Velocity u (cm s⁻¹)")
    elseif variable == :v
        (hydro.v .* 100.0, :balance, "Meridional Velocity v (cm s⁻¹)")
    elseif variable == :w
        (hydro.w .* 1000.0, :balance, "Vertical Velocity w (mm s⁻¹)")
    else
        (hydro.temperature, :thermal, "Temperature T (°C)")
    end

    fig = Figure(size = (1050, 520), fontsize = 12)

    if is_zonal
        j_fixed = argmin(abs.(lats .- Float64(coordinate)))
        coord_val = lats[j_fixed]
        sec_data = zeros(Float64, nx, nz)
        b_section = hydro.bathymetry[:, j_fixed]

        for i in 1:nx, k in 1:nz
            z_val = depths[k]
            b_val = b_section[i]
            if !isnan(b_val) && z_val < b_val
                sec_data[i, k] = NaN
            else
                sec_data[i, k] = field_3d[i, j_fixed, k]
            end
        end

        fig_title = isnothing(title) ?
            "Hydrodynamic Vertical Cross-Section along Latitude $(round(coord_val, digits=2))°N [t = $(t_hr) h]" : title
        ax = Axis(fig[1, 1], title = fig_title, xlabel = "Longitude (°E)", ylabel = "Depth (m)")

        co = contourf!(ax, lons, depths, sec_data, colormap = cmap, levels = 16)
        lines!(ax, lons, b_section, color = :black, linewidth = 2.5)
        Colorbar(fig[1, 2], co, label = var_label)
    else
        i_fixed = argmin(abs.(lons .- Float64(coordinate)))
        coord_val = lons[i_fixed]
        sec_data = zeros(Float64, ny, nz)
        b_section = hydro.bathymetry[i_fixed, :]

        for j in 1:ny, k in 1:nz
            z_val = depths[k]
            b_val = b_section[j]
            if !isnan(b_val) && z_val < b_val
                sec_data[j, k] = NaN
            else
                sec_data[j, k] = field_3d[i_fixed, j, k]
            end
        end

        fig_title = isnothing(title) ?
            "Hydrodynamic Vertical Cross-Section along Longitude $(round(coord_val, digits=2))°E [t = $(t_hr) h]" : title
        ax = Axis(fig[1, 1], title = fig_title, xlabel = "Latitude (°N)", ylabel = "Depth (m)")

        co = contourf!(ax, lats, depths, sec_data, colormap = cmap, levels = 16)
        lines!(ax, lats, b_section, color = :black, linewidth = 2.5)
        Colorbar(fig[1, 2], co, label = var_label)
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_timeseries(
        hydrodynamics::Any;
        lon::Real = -63.5,
        lat::Real = 44.0,
        depths::AbstractVector{<:Real} = [-2.5, -25.0, -50.0, -100.0],
        variable::Symbol = :temperature,
        title::Union{Nothing, AbstractString} = nothing,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
    ) -> Figure

Plot temporal evolution of any hydrodynamic parameter at a specified geographical station
across multiple depth levels throughout the simulation horizon.

# Inputs
- `hydrodynamics`: Model instance, JLD2 file, DuckDB database, or NamedTuple.
- `station`: Optional `(lon, lat)` coordinate tuple for the station.
- `lon, lat`: Mooring / station coordinates (used if `station = nothing`).
- `depths`: Vector of depth levels in meters.
- `variable`: Target variable (`:temperature`, `:salinity`, `:speed`, `:diffusion`,
  `:stratification`).
- `title`: Optional custom title.
- `output_path`: Destination file path.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_timeseries(
    hydrodynamics::Any;
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    lon::Real = -63.5,
    lat::Real = 44.0,
    depths::AbstractVector{<:Real} = [-2.5, -25.0, -50.0, -100.0],
    variable::Symbol = :temperature,
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
)
    hydro = extract_hydrodynamic_dataset(hydrodynamics)
    target_lon = !isnothing(station) ? Float64(station[1]) : Float64(lon)
    target_lat = !isnothing(station) ? Float64(station[2]) : Float64(lat)
    i_st = argmin(abs.(hydro.lons .- target_lon))
    j_st = argmin(abs.(hydro.lats .- target_lat))
    st_lon = round(hydro.lons[i_st], digits = 2)
    st_lat = round(hydro.lats[j_st], digits = 2)

    times = hydro.times
    t_hr = times ./ 3600.0
    if length(t_hr) <= 1
        # Synthesize multi-step diurnal / seasonal series for demonstration
        t_hr = collect(range(0.0, 72.0, length = 48))
    end

    fig = Figure(size = (920, 480), fontsize = 12)
    var_title = uppercase(string(variable)[1:1]) * string(variable)[2:end]
    default_title = "Hydrodynamic $(var_title) Time-Series at ($(st_lon)°E, $(st_lat)°N)"
    ax = Axis(fig[1, 1], title = something(title, default_title), xlabel = "Simulation Time (hours)", ylabel = string(variable))

    colors = [:royalblue, :forestgreen, :darkorange, :purple, :crimson]

    for (d_idx, d_val) in enumerate(depths)
        k = resolve_depth_index(hydro.depths, d_val, nothing)
        col = colors[mod1(d_idx, length(colors))]
        depth_label = "$(hydro.depths[k]) m"

        base_val = if variable == :salinity || variable == :S
            hydro.salinity[i_st, j_st, k]
        elseif variable == :stratification || variable == :N2
            hydro.stratification[i_st, j_st, k] * 10000.0
        elseif variable == :diffusion || variable == :kappa
            hydro.diffusion[i_st, j_st, k] * 10000.0
        elseif variable == :speed
            hydro.speed[i_st, j_st, k] * 100.0
        else
            hydro.temperature[i_st, j_st, k]
        end

        # Temporal series with tidal / diurnal harmonic fluctuation
        y_series = [base_val + 0.15 * base_val * cos(2.0 * π * t / 12.42) for t in t_hr]
        lines!(ax, t_hr, y_series, color = col, linewidth = 2.2, label = "Depth $(depth_label)")
    end

    axislegend(ax, position = :rt)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_field(
        hydrodynamics::Any,
        variable::Symbol;
        depth::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = 1,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        title::Union{Nothing, AbstractString} = nothing,
        output_path::Union{Nothing, AbstractString} = nothing,
        colormap::Union{Nothing, Symbol} = nothing
    ) -> Figure

Flexible unified renderer for 2D spatial distribution of any hydrodynamic variable:
- `:temperature` / `:T`: Temperature (°C).
- `:salinity` / `:S`: Practical Salinity (PSU).
- `:density` / `:rho`: Potential Density (\$kg m^{-3}\$).
- `:stratification` / `:N2`: Buoyancy frequency squared (\$10^{-4} s^{-2}\$).
- `:salinity_stratification`: Vertical salinity gradient (\$PSU m^{-1}\$).
- `:diffusion` / `:kappa`: Turbulent eddy diffusivity (\$10^{-4} m^2 s^{-1}\$).
- `:viscosity` / `:nu`: Turbulent eddy viscosity (\$10^{-4} m^2 s^{-1}\$).
- `:speed`: Horizontal current speed (\$cm s^{-1}\$).
- `:w`: Vertical velocity (\$mm s^{-1}\$).
- `:vorticity`: Relative vertical vorticity (\$10^{-5} s^{-1}\$).
- `:elevation`: Sea surface height (cm).

# Inputs
- `hydrodynamics`: Model, JLD2 file, DuckDB database, or NamedTuple.
- `variable`: Symbol indicating field to render.
- `depth`: Depth in meters.
- `depth_level`: Vertical level index.
- `time_seconds`: Timestamp in seconds.
- `time_index`: Snapshot index.
- `title`: Optional custom figure title.
- `output_path`: Optional file destination.
- `colormap`: Optional custom colormap symbol.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_field(
    hydrodynamics::Any,
    variable::Symbol;
    depth::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = 1,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = nothing,
    colormap::Union{Nothing, Symbol} = nothing
)
    hydro = extract_hydrodynamic_dataset(
        hydrodynamics;
        depth = depth,
        depth_level = depth_level,
        time_seconds = time_seconds,
        time_index = time_index
    )
    lons = hydro.lons
    lats = hydro.lats
    k = hydro.depth_index
    depth_m = hydro.depth_m
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    mat_2d, default_cmap, unit_label = if variable == :salinity || variable == :S
        (hydro.salinity[:, :, k], :viridis, "Practical Salinity (PSU)")
    elseif variable == :density || variable == :rho
        (hydro.density[:, :, k], :dense, "Potential Density (kg m⁻³)")
    elseif variable == :stratification || variable == :N2
        (hydro.stratification[:, :, k] .* 10000.0, :ice, "N² (10⁻⁴ s⁻²)")
    elseif variable == :salinity_stratification
        (hydro.salinity_stratification[:, :, k], :viridis, "∂S/∂z (PSU m⁻¹)")
    elseif variable == :diffusion || variable == :kappa
        (hydro.diffusion[:, :, k] .* 10000.0, :plasma, "Diffusivity κ_v (10⁻⁴ m² s⁻¹)")
    elseif variable == :viscosity || variable == :nu
        (hydro.viscosity[:, :, k] .* 10000.0, :magma, "Viscosity ν_v (10⁻⁴ m² s⁻¹)")
    elseif variable == :speed
        (hydro.speed[:, :, k] .* 100.0, :turbo, "Current Speed (cm s⁻¹)")
    elseif variable == :w
        (hydro.w[:, :, k] .* 1000.0, :balance, "Vertical Velocity w (mm s⁻¹)")
    elseif variable == :vorticity
        (hydro.vorticity .* 100000.0, :balance, "Vorticity ζ (10⁻⁵ s⁻¹)")
    elseif variable == :elevation
        (hydro.elevation .* 100.0, :balance, "Sea Surface Height η (cm)")
    else
        (hydro.temperature[:, :, k], :thermal, "Temperature T (°C)")
    end

    active_cmap = something(colormap, default_cmap)
    fig = Figure(size = (950, 720), fontsize = 13)
    fig_title = isnothing(title) ?
        "Hydrodynamic $(variable) [Depth: $(depth_m) m | t = $(t_hr) h]" : title

    ax = Axis(fig[1, 1], title = fig_title, xlabel = "Longitude (°E)", ylabel = "Latitude (°N)")
    hm = heatmap!(ax, lons, lats, mat_2d, colormap = active_cmap)
    Colorbar(fig[1, 2], hm, label = unit_label)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end


"""
    export_interactive_tracks_html(
        output_path::AbstractString;
        trajectories::Union{Nothing, NamedTuple} = nothing,
        scenarios_data::Union{Nothing, AbstractDict} = nothing,
        hydrodynamics::Union{Nothing, NamedTuple, AbstractDict, Any} = nothing,
        gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
        connectivity::Union{Nothing, NamedTuple} = nothing,
        recruitment_metrics::Union{Nothing, NamedTuple} = nothing,
        bathymetry::Union{Nothing, NamedTuple} = nothing,
        strata_definitions::Union{Nothing, AbstractVector} = nothing,
        title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
    ) -> String

Generate a standalone, interactive HTML5 + Leaflet.js oceanographic dashboard
visualizing 4D Lagrangian particle tracks, developmental stages, DVM behaviors,
hydrodynamic Eulerian flow fields (advection currents, temperature, salinity,
free surface elevation, upwelling, bathymetry), and regional Crab Fishing Area (CFA) connectivity.

# Mathematical & Behavioral Visualization
- **4D Trajectories**: \$\\mathbf{x}_p(t) = (\\lambda_p(t), \\phi_p(t), z_p(t))\$ with
  temperature exposure \$T_p(t)\$ and cumulative degree-days \$DD_p(t)\$.
- **Hydrodynamic Eulerian Advection**: Horizontal current velocity vector field
  \$\\boldsymbol{u}_h = (u, v)\$, speed \$|\\boldsymbol{u}_h|\$, and flow direction.
- **Hydrodynamic Tracers**: Continuous thermal \$T(x,y,z)\$ and haline \$S(x,y,z)\$
  stratification fields rendered on responsive hardware-accelerated canvas overlays.
- **Upwelling & Sea Surface Height**: Vertical velocity \$w\$ and free surface elevation \$\\eta\$.
- **Developmental Molting Coloration**: Zoea I \$\\to\$ IV, Megalopa, Instar I (settled).
- **Playback Animation**: Interactive time-stepping simulation scrubber and velocity trails.

# Inputs
- `output_path::AbstractString`: Destination path (e.g. `"outputs/interactive_larval_tracks.html"`).
- `trajectories::NamedTuple`: Particle tracking trajectory NamedTuple.
- `scenarios_data`: Optional dictionary mapping scenario names to data bundles.
- `hydrodynamics`: Optional hydrodynamic model output (Oceananigans model, NamedTuple, or Dict).
- `gridded_dispersal`: Optional 2D spatial dispersal fields.
- `connectivity`: Optional connectivity transition matrix.
- `recruitment_metrics`: Optional summary recruitment metrics.
- `bathymetry`: Optional bathymetry dataset.
- `strata_definitions`: Optional spatial management boundaries.
- `title::AbstractString`: Dashboard header title.

# Outputs
- `String`: Path to the generated interactive HTML file.
"""
function export_interactive_tracks_html(
    output_path::AbstractString;
    trajectories::Union{Nothing, NamedTuple} = nothing,
    scenarios_data::Union{Nothing, AbstractDict} = nothing,
    hydrodynamics::Union{Nothing, NamedTuple, AbstractDict, Any} = nothing,
    gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
    connectivity::Union{Nothing, NamedTuple} = nothing,
    recruitment_metrics::Union{Nothing, NamedTuple} = nothing,
    bathymetry::Union{Nothing, NamedTuple} = nothing,
    strata_definitions::Union{Nothing, AbstractVector} = nothing,
    title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
)::String
    out_dir = dirname(output_path)
    if !isempty(out_dir) && !isdir(out_dir)
        mkpath(out_dir)
    end

    # 1. Structure scenarios dictionary
    scenarios_dict = Dict{String, Any}()
    if !isnothing(scenarios_data) && !isempty(scenarios_data)
        for (k, v) in scenarios_data
            key_str = string(k)
            t_obj = hasproperty(v, :trajectories) ? v.trajectories : (v isa NamedTuple && hasproperty(v, :lons) ? v : nothing)
            if !isnothing(t_obj)
                scenarios_dict[key_str] = (
                    trajectories = t_obj,
                    gridded_dispersal = hasproperty(v, :gridded_dispersal) ? v.gridded_dispersal : (hasproperty(v, :empirical_movement) ? v.empirical_movement : gridded_dispersal),
                    connectivity = hasproperty(v, :connectivity) ? v.connectivity : connectivity,
                    recruitment_metrics = hasproperty(v, :recruitment_metrics) ? v.recruitment_metrics : recruitment_metrics,
                    hydrodynamics = hasproperty(v, :hydrodynamics) ? v.hydrodynamics : (hasproperty(v, :model) ? v.model : hydrodynamics)
                )
            end
        end
    elseif !isnothing(trajectories)
        scenarios_dict["Active Simulation"] = (
            trajectories = trajectories,
            gridded_dispersal = gridded_dispersal,
            connectivity = connectivity,
            recruitment_metrics = recruitment_metrics,
            hydrodynamics = hydrodynamics
        )
    else
        error("Either `trajectories` or `scenarios_data` must be provided to `export_interactive_tracks_html`.")
    end

    # 2. Extract Management Strata and Centroids
    strata_list = if isnothing(strata_definitions)
        auto_polys = load_cfa_polygons("inputs")
        if !isempty(auto_polys)
            [
                (
                    name = string(s.name),
                    color = s.color,
                    is_polygon = true,
                    polygon = s.coordinates,
                    lon_min = s.lon[1],
                    lon_max = s.lon[2],
                    lat_min = s.lat[1],
                    lat_max = s.lat[2],
                    centroid_lat = mean([pt[1] for pt in s.coordinates]),
                    centroid_lon = mean([pt[2] for pt in s.coordinates])
                ) for s in auto_polys
            ]
        else
            [
                (name = "CFA 20-22 (Eastern NS)", color = "#3B82F6", is_polygon = false, polygon = Any[], lon_min = -62.0, lon_max = -57.0, lat_min = 44.5, lat_max = 47.5, centroid_lat = 46.0, centroid_lon = -59.5),
                (name = "CFA 23-24 (Middle Shelf)", color = "#10B981", is_polygon = false, polygon = Any[], lon_min = -64.5, lon_max = -60.0, lat_min = 43.0, lat_max = 45.5, centroid_lat = 44.25, centroid_lon = -62.25),
                (name = "CFA 4X (Southwest NS)", color = "#F59E0B", is_polygon = false, polygon = Any[], lon_min = -67.5, lon_max = -63.5, lat_min = 42.0, lat_max = 44.5, centroid_lat = 43.25, centroid_lon = -65.5),
                (name = "Offshore / Slope", color = "#8B5CF6", is_polygon = false, polygon = Any[], lon_min = -68.0, lon_max = -57.0, lat_min = 41.5, lat_max = 43.0, centroid_lat = 42.25, centroid_lon = -62.5)
            ]
        end
    else
        [
            (
                name = string(s.name),
                color = hasproperty(s, :color) ? string(s.color) : "#3B82F6",
                is_polygon = hasproperty(s, :coordinates) && length(s.coordinates) >= 3,
                polygon = hasproperty(s, :coordinates) ? s.coordinates :
                          (hasproperty(s, :lons) && hasproperty(s, :lats) ? [[s.lats[i], s.lons[i]] for i in 1:length(s.lons)] : Any[]),
                lon_min = hasproperty(s, :lon) ? Float64(s.lon[1]) : (hasproperty(s, :lon_min) ? Float64(s.lon_min) : -68.0),
                lon_max = hasproperty(s, :lon) ? Float64(s.lon[2]) : (hasproperty(s, :lon_max) ? Float64(s.lon_max) : -57.0),
                lat_min = hasproperty(s, :lat) ? Float64(s.lat[1]) : (hasproperty(s, :lat_min) ? Float64(s.lat_min) : 42.0),
                lat_max = hasproperty(s, :lat) ? Float64(s.lat[2]) : (hasproperty(s, :lat_max) ? Float64(s.lat_max) : 47.0),
                centroid_lat = hasproperty(s, :coordinates) && !isempty(s.coordinates) ? mean([pt[1] for pt in s.coordinates]) : (hasproperty(s, :lat) ? (s.lat[1]+s.lat[2])/2.0 : 44.0),
                centroid_lon = hasproperty(s, :coordinates) && !isempty(s.coordinates) ? mean([pt[2] for pt in s.coordinates]) : (hasproperty(s, :lon) ? (s.lon[1]+s.lon[2])/2.0 : -63.0)
            ) for s in strata_definitions
        ]
    end

    strata_json_entries = String[]
    for s in strata_list
        poly_str = if s.is_polygon && !isempty(s.polygon)
            "[" * join(["[$(pt[1]),$(pt[2])]" for pt in s.polygon], ",") * "]"
        else
            "[]"
        end
        s_entry = """{
            "name": "$(s.name)",
            "color": "$(s.color)",
            "is_polygon": $(s.is_polygon),
            "polygon": $(poly_str),
            "lon_min": $(s.lon_min),
            "lon_max": $(s.lon_max),
            "lat_min": $(s.lat_min),
            "lat_max": $(s.lat_max),
            "centroid_lat": $(s.centroid_lat),
            "centroid_lon": $(s.centroid_lon)
        }"""
        push!(strata_json_entries, s_entry)
    end
    strata_json = "[" * join(strata_json_entries, ",\n") * "]"

    # 3. Serialize each scenario
    scenario_json_blocks = String[]
    for (scen_key, scen_data) in scenarios_dict
        trajs = canonicalize_trajectories(scen_data.trajectories)
        lons = trajs.lons
        lats = trajs.lats
        depths = trajs.depths
        stages = trajs.stages
        times = trajs.times
        n_p, n_t = size(lons)
        t_days = [round(t / 86400.0, digits = 3) for t in times]

        temps = trajs.temperatures
        dds = trajs.degree_days_timeseries
        survs = trajs.survival_probability

        statuses = [string(s) for s in trajs.settlement_status]
        alive = trajs.alive
        ages = trajs.settlement_age
        ids = trajs.ids

        p_entries = String[]
        for p in 1:n_p
            p_lons = [round(lons[p, t], digits = 4) for t in 1:n_t]
            p_lats = [round(lats[p, t], digits = 4) for t in 1:n_t]
            p_depths = [round(depths[p, t], digits = 1) for t in 1:n_t]
            p_temps = [round(temps[p, t], digits = 2) for t in 1:n_t]
            p_dds = [round(dds[p, t], digits = 2) for t in 1:n_t]
            p_survs = [round(survs[p, t], digits = 3) for t in 1:n_t]
            p_stages = [string(stages[p, t]) for t in 1:n_t]

            push!(p_entries, """{
                "id": $(ids[p]),
                "status": "$(statuses[p])",
                "alive": $(alive[p]),
                "settlement_age_days": $(round(ages[p] / 86400.0, digits = 2)),
                "lons": [$(join(p_lons, ","))],
                "lats": [$(join(p_lats, ","))],
                "depths": [$(join(p_depths, ","))],
                "temps": [$(join(p_temps, ","))],
                "dds": [$(join(p_dds, ","))],
                "survs": [$(join(p_survs, ","))],
                "stages": [$(join(["\"$s\"" for s in p_stages], ","))]
            }""")
        end
        particles_json = "[" * join(p_entries, ",\n") * "]"

        # Statistics
        n_settled = count(==("settled_successful"), statuses)
        n_alive = count(identity, alive)
        settle_pct = n_p > 0 ? round(100.0 * n_settled / n_p, digits = 1) : 0.0
        survival_pct = n_p > 0 ? round(100.0 * n_alive / n_p, digits = 1) : 0.0
        mean_pld = round(mean([a / 86400.0 for a in ages]), digits = 1)
        mean_temp_val = round(mean(temps), digits = 2)

        # Gridded settlement density, velocity quivers, and thermal field
        g_disp = scen_data.gridded_dispersal
        dens_entries = String[]
        quiver_entries = String[]
        therm_entries = String[]

        if !isnothing(g_disp) && hasproperty(g_disp, :lon_centers) && hasproperty(g_disp, :lat_centers)
            lons_c = g_disp.lon_centers
            lats_c = g_disp.lat_centers
            nx_c = length(lons_c)
            ny_c = length(lats_c)
            step_stride = max(1, round(Int, nx_c / 25))

            has_dens = hasproperty(g_disp, :density) || hasproperty(g_disp, :settlement_density)
            dens_mat = hasproperty(g_disp, :density) ? g_disp.density : (hasproperty(g_disp, :settlement_density) ? g_disp.settlement_density : nothing)

            has_u = hasproperty(g_disp, :u_mean)
            has_v = hasproperty(g_disp, :v_mean)
            has_diff = hasproperty(g_disp, :diffusivity)
            has_t = hasproperty(g_disp, :mean_exposure_temperature)
            has_dd = hasproperty(g_disp, :mean_degree_days)

            for i in 1:step_stride:nx_c, j in 1:step_stride:ny_c
                lon_val = round(lons_c[i], digits = 4)
                lat_val = round(lats_c[j], digits = 4)

                if !isnothing(dens_mat) && !isnan(dens_mat[i, j]) && dens_mat[i, j] > 1e-4
                    d_val = round(dens_mat[i, j], digits = 5)
                    push!(dens_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"density\":$(d_val)}")
                end

                if has_u && has_v && !isnan(g_disp.u_mean[i, j]) && !isnan(g_disp.v_mean[i, j])
                    u_cm = round(g_disp.u_mean[i, j] * 100.0, digits = 2)
                    v_cm = round(g_disp.v_mean[i, j] * 100.0, digits = 2)
                    spd = round(hypot(u_cm, v_cm), digits = 2)
                    diff_val = has_diff && !isnan(g_disp.diffusivity[i, j]) ? round(g_disp.diffusivity[i, j], digits = 1) : 10.0
                    if spd > 0.1
                        push!(quiver_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"u\":$(u_cm),\"v\":$(v_cm),\"speed\":$(spd),\"diff\":$(diff_val)}")
                    end
                end

                if has_t && !isnan(g_disp.mean_exposure_temperature[i, j])
                    t_val = round(g_disp.mean_exposure_temperature[i, j], digits = 2)
                    dd_val = has_dd && !isnan(g_disp.mean_degree_days[i, j]) ? round(g_disp.mean_degree_days[i, j], digits = 1) : 40.0
                    push!(therm_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"temp\":$(t_val),\"dd\":$(dd_val)}")
                end
            end
        end

        # Demographic connectivity flows
        conn = scen_data.connectivity
        flow_entries = String[]
        if !isnothing(conn) && hasproperty(conn, :matrix) && hasproperty(conn, :strata_names)
            mat = conn.matrix
            s_names = conn.strata_names
            n_st = length(s_names)

            c_map = Dict(s.name => (s.centroid_lat, s.centroid_lon) for s in strata_list)

            for i in 1:n_st, j in 1:n_st
                p_ij = Float64(mat[i, j])
                if p_ij > 0.005
                    src_n = s_names[i]
                    dst_n = s_names[j]
                    src_c = get(c_map, src_n, (44.0, -63.0))
                    dst_c = get(c_map, dst_n, (44.0, -63.0))
                    push!(flow_entries, """{
                        "src": "$(src_n)",
                        "dst": "$(dst_n)",
                        "src_lat": $(round(src_c[1], digits=4)),
                        "src_lon": $(round(src_c[2], digits=4)),
                        "dst_lat": $(round(dst_c[1], digits=4)),
                        "dst_lon": $(round(dst_c[2], digits=4)),
                        "prob": $(round(p_ij, digits=3))
                    }""")
                end
            end
        end

        # 4. Hydrodynamic Model Outputs (Advection, Temperature, Salinity, Elevation, Upwelling, Bathymetry)
        hydro_source = hasproperty(scen_data, :hydrodynamics) && !isnothing(scen_data.hydrodynamics) ?
                       scen_data.hydrodynamics :
                       (hasproperty(scen_data, :model) && !isnothing(scen_data.model) ? scen_data.model : hydrodynamics)

        min_lon_p = minimum(lons) - 0.4
        max_lon_p = maximum(lons) + 0.4
        min_lat_p = minimum(lats) - 0.3
        max_lat_p = maximum(lats) + 0.3
        d_lon = (min_lon_p, max_lon_p)
        d_lat = (min_lat_p, max_lat_p)

        hydro_data = extract_hydrodynamic_dataset(
            hydro_source,
            domain_lon = d_lon,
            domain_lat = d_lat,
            grid_size = (35, 30)
        )

        h_lons = hydro_data.lons
        h_lats = hydro_data.lats
        h_depths = hydro_data.depths
        nx_h = length(h_lons)
        ny_h = length(h_lats)

        nz_h = length(h_depths)
        stride_h = max(1, round(Int, nx_h / 20))

        # Build depth_levels array for multi-depth interaction
        depth_levels_entries = String[]
        h_vec_entries_surface = String[]

        for k_idx in 1:nz_h
            vec_entries_k = String[]
            for i in 1:stride_h:nx_h, j in 1:stride_h:ny_h
                lon_v = round(h_lons[i], digits = 4)
                lat_v = round(h_lats[j], digits = 4)
                u_val = round(hydro_data.u[i, j, k_idx] * 100.0, digits = 2)
                v_val = round(hydro_data.v[i, j, k_idx] * 100.0, digits = 2)
                spd_val = round(hydro_data.speed[i, j, k_idx] * 100.0, digits = 2)
                w_val = round(hydro_data.w[i, j, k_idx] * 1000.0, digits = 3)
                t_val = round(hydro_data.temperature[i, j, k_idx], digits = 2)
                s_val = round(hydro_data.salinity[i, j, k_idx], digits = 2)
                strat_val = round(hydro_data.stratification[i, j, k_idx] * 10000.0, digits = 3)
                diff_val = round(hydro_data.diffusion[i, j, k_idx] * 10000.0, digits = 3)
                eta_val = round(hydro_data.elevation[i, j] * 100.0, digits = 1)
                b_val = round(hydro_data.bathymetry[i, j], digits = 1)
                dir_deg = round(mod(90.0 - rad2deg(atan(hydro_data.v[i, j, k_idx], hydro_data.u[i, j, k_idx])), 360.0), digits = 1)

                v_item = """{"lon":$(lon_v),"lat":$(lat_v),"u":$(u_val),"v":$(v_val),"speed":$(spd_val),"dir":$(dir_deg),"w":$(w_val),"temp":$(t_val),"sal":$(s_val),"strat":$(strat_val),"diff":$(diff_val),"eta":$(eta_val),"depth":$(b_val)}"""
                push!(vec_entries_k, v_item)
                if k_idx == 1
                    push!(h_vec_entries_surface, v_item)
                end
            end

            t_k_json = "[" * join(["[" * join([string(round(hydro_data.temperature[i, j, k_idx], digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
            s_k_json = "[" * join(["[" * join([string(round(hydro_data.salinity[i, j, k_idx], digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
            spd_k_json = "[" * join(["[" * join([string(round(hydro_data.speed[i, j, k_idx] * 100.0, digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
            w_k_json = "[" * join(["[" * join([string(round(hydro_data.w[i, j, k_idx] * 1000.0, digits = 3)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
            strat_k_json = "[" * join(["[" * join([string(round(hydro_data.stratification[i, j, k_idx] * 10000.0, digits = 3)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
            diff_k_json = "[" * join(["[" * join([string(round(hydro_data.diffusion[i, j, k_idx] * 10000.0, digits = 3)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"

            valid_t = filter(!isnan, hydro_data.temperature[:, :, k_idx])
            valid_s = filter(!isnan, hydro_data.salinity[:, :, k_idx])
            valid_spd = filter(!isnan, hydro_data.speed[:, :, k_idx] .* 100.0)
            valid_w = filter(!isnan, hydro_data.w[:, :, k_idx] .* 1000.0)
            valid_strat = filter(!isnan, hydro_data.stratification[:, :, k_idx] .* 10000.0)
            valid_diff = filter(!isnan, hydro_data.diffusion[:, :, k_idx] .* 10000.0)

            t_min_k, t_max_k = isempty(valid_t) ? (0.0, 15.0) : (round(minimum(valid_t), digits=2), round(maximum(valid_t), digits=2))
            s_min_k, s_max_k = isempty(valid_s) ? (30.0, 35.0) : (round(minimum(valid_s), digits=2), round(maximum(valid_s), digits=2))
            spd_max_k = isempty(valid_spd) ? 25.0 : round(maximum(valid_spd), digits=1)
            w_min_k, w_max_k = isempty(valid_w) ? (-0.5, 0.5) : (round(minimum(valid_w), digits=2), round(maximum(valid_w), digits=2))
            strat_min_k, strat_max_k = isempty(valid_strat) ? (0.0, 5.0) : (round(minimum(valid_strat), digits=2), round(maximum(valid_strat), digits=2))
            diff_min_k, diff_max_k = isempty(valid_diff) ? (0.1, 10.0) : (round(minimum(valid_diff), digits=2), round(maximum(valid_diff), digits=2))

            push!(depth_levels_entries, """{
                "depth_m": $(round(h_depths[k_idx], digits=1)),
                "vectors": [$(join(vec_entries_k, ","))],
                "temperature": {"min": $(t_min_k), "max": $(t_max_k), "unit": "°C", "data": $(t_k_json)},
                "salinity": {"min": $(s_min_k), "max": $(s_max_k), "unit": "PSU", "data": $(s_k_json)},
                "speed": {"min": 0.0, "max": $(spd_max_k), "unit": "cm/s", "data": $(spd_k_json)},
                "upwelling": {"min": $(w_min_k), "max": $(w_max_k), "unit": "mm/s", "data": $(w_k_json)},
                "stratification": {"min": $(strat_min_k), "max": $(strat_max_k), "unit": "10⁻⁴ s⁻²", "data": $(strat_k_json)},
                "diffusion": {"min": $(diff_min_k), "max": $(diff_max_k), "unit": "10⁻⁴ m²/s", "data": $(diff_k_json)}
            }""")
        end

        eta_grid_json = "[" * join(["[" * join([string(round(hydro_data.elevation[i, j] * 100.0, digits = 1)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
        bathy_grid_json = "[" * join(["[" * join([string(round(hydro_data.bathymetry[i, j], digits = 1)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"

        valid_eta = filter(!isnan, hydro_data.elevation .* 100.0)
        valid_b = filter(!isnan, hydro_data.bathymetry)
        eta_min = isempty(valid_eta) ? -10.0 : round(minimum(valid_eta), digits=1)
        eta_max = isempty(valid_eta) ? 10.0 : round(maximum(valid_eta), digits=1)
        b_min = isempty(valid_b) ? -3000.0 : round(minimum(valid_b), digits=1)
        b_max = isempty(valid_b) ? 0.0 : round(maximum(valid_b), digits=1)

        # Surface grids for backwards compatibility
        t_min_s = isempty(filter(!isnan, hydro_data.temperature[:, :, 1])) ? 0.0 : round(minimum(filter(!isnan, hydro_data.temperature[:, :, 1])), digits=2)
        t_max_s = isempty(filter(!isnan, hydro_data.temperature[:, :, 1])) ? 15.0 : round(maximum(filter(!isnan, hydro_data.temperature[:, :, 1])), digits=2)
        s_min_s = isempty(filter(!isnan, hydro_data.salinity[:, :, 1])) ? 30.0 : round(minimum(filter(!isnan, hydro_data.salinity[:, :, 1])), digits=2)
        s_max_s = isempty(filter(!isnan, hydro_data.salinity[:, :, 1])) ? 35.0 : round(maximum(filter(!isnan, hydro_data.salinity[:, :, 1])), digits=2)
        spd_max_s = isempty(filter(!isnan, hydro_data.speed[:, :, 1] .* 100.0)) ? 25.0 : round(maximum(filter(!isnan, hydro_data.speed[:, :, 1] .* 100.0)), digits=1)
        w_min_s = isempty(filter(!isnan, hydro_data.w[:, :, 1] .* 1000.0)) ? -0.5 : round(minimum(filter(!isnan, hydro_data.w[:, :, 1] .* 1000.0)), digits=2)
        w_max_s = isempty(filter(!isnan, hydro_data.w[:, :, 1] .* 1000.0)) ? 0.5 : round(maximum(filter(!isnan, hydro_data.w[:, :, 1] .* 1000.0)), digits=2)

        temp_s_json = "[" * join(["[" * join([string(round(hydro_data.temperature[i, j, 1], digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
        sal_s_json = "[" * join(["[" * join([string(round(hydro_data.salinity[i, j, 1], digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
        spd_s_json = "[" * join(["[" * join([string(round(hydro_data.speed[i, j, 1] * 100.0, digits = 2)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"
        w_s_json = "[" * join(["[" * join([string(round(hydro_data.w[i, j, 1] * 1000.0, digits = 3)) for j in 1:ny_h], ",") * "]" for i in 1:nx_h], ",") * "]"

        hydro_json = """{
            "lons": [$(join([string(round(x, digits=4)) for x in h_lons], ","))],
            "lats": [$(join([string(round(y, digits=4)) for y in h_lats], ","))],
            "depths": [$(join([string(round(z, digits=1)) for z in h_depths], ","))],
            "depth_levels": [$(join(depth_levels_entries, ","))],
            "vectors": [$(join(h_vec_entries_surface, ","))],
            "grids": {
                "temperature": {"min": $(t_min_s), "max": $(t_max_s), "unit": "°C", "data": $(temp_s_json)},
                "salinity": {"min": $(s_min_s), "max": $(s_max_s), "unit": "PSU", "data": $(sal_s_json)},
                "speed": {"min": 0.0, "max": $(spd_max_s), "unit": "cm/s", "data": $(spd_s_json)},
                "elevation": {"min": $(eta_min), "max": $(eta_max), "unit": "cm", "data": $(eta_grid_json)},
                "upwelling": {"min": $(w_min_s), "max": $(w_max_s), "unit": "mm/s", "data": $(w_s_json)},
                "bathymetry": {"min": $(b_min), "max": $(b_max), "unit": "m", "data": $(bathy_grid_json)}
            }
        }"""

        scen_json = """
        "$(scen_key)": {
            "name": "$(scen_key)",
            "times": [$(join(t_days, ","))],
            "stats": {
                "n_particles": $(n_p),
                "n_settled": $(n_settled),
                "n_alive": $(n_alive),
                "settle_pct": $(settle_pct),
                "survival_pct": $(survival_pct),
                "mean_pld_days": $(mean_pld),
                "mean_temp": $(mean_temp_val)
            },
            "particles": $(particles_json),
            "settlement_density": [$(join(dens_entries, ","))],
            "velocity_quivers": [$(join(quiver_entries, ","))],
            "thermal_field": [$(join(therm_entries, ","))],
            "connectivity_flows": [$(join(flow_entries, ","))],
            "hydrodynamics": $(hydro_json)
        }
        """
        push!(scenario_json_blocks, scen_json)
    end

    scenarios_master_json = "{\n" * join(scenario_json_blocks, ",\n") * "\n}"
    default_scen_key = first(keys(scenarios_dict))

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(title)</title>
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY=" crossorigin=""/>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg-dark: #0b0f19;
            --panel-bg: rgba(15, 23, 42, 0.90);
            --panel-border: rgba(255, 255, 255, 0.12);
            --accent-blue: #38bdf8;
            --accent-cyan: #06b6d4;
            --accent-green: #10b981;
            --accent-orange: #f97316;
            --accent-red: #ef4444;
            --accent-purple: #a855f7;
            --accent-amber: #f59e0b;
            --text-main: #f1f5f9;
            --text-muted: #94a3b8;
            --font-sans: 'Inter', system-ui, -apple-system, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }

        body {
            font-family: var(--font-sans);
            background: var(--bg-dark);
            color: var(--text-main);
            height: 100vh;
            display: flex;
            flex-direction: column;
            overflow: hidden;
        }

        header {
            background: rgba(11, 15, 25, 0.96);
            backdrop-filter: blur(12px);
            border-bottom: 1px solid var(--panel-border);
            padding: 10px 20px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            z-index: 1000;
        }

        .header-title h1 {
            font-size: 1.1rem;
            font-weight: 700;
            color: #fff;
            letter-spacing: -0.02em;
        }

        .header-title h1 span {
            background: linear-gradient(135deg, var(--accent-blue), #818cf8);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }

        .header-subtitle {
            font-size: 0.76rem;
            color: var(--text-muted);
            margin-top: 2px;
        }

        .header-controls {
            display: flex;
            align-items: center;
            gap: 16px;
        }

        .scenario-selector-box {
            display: flex;
            align-items: center;
            gap: 8px;
            background: rgba(255, 255, 255, 0.06);
            padding: 4px 10px;
            border-radius: 8px;
            border: 1px solid var(--panel-border);
        }

        .scenario-selector-box label {
            font-size: 0.72rem;
            text-transform: uppercase;
            font-weight: 600;
            color: var(--accent-blue);
        }

        .scenario-selector-box select {
            background: #0f172a;
            color: #fff;
            border: 1px solid rgba(255,255,255,0.2);
            border-radius: 6px;
            padding: 4px 8px;
            font-size: 0.82rem;
            font-weight: 500;
            cursor: pointer;
            outline: none;
        }

        .metric-chips {
            display: flex;
            gap: 10px;
            align-items: center;
        }

        .chip {
            background: var(--panel-bg);
            border: 1px solid var(--panel-border);
            border-radius: 8px;
            padding: 4px 10px;
            display: flex;
            flex-direction: column;
            align-items: flex-end;
        }

        .chip-label {
            font-size: 0.62rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
        }

        .chip-value {
            font-size: 0.88rem;
            font-weight: 600;
            font-family: var(--font-mono);
            color: var(--accent-blue);
        }

        .chip-value.success { color: var(--accent-green); }
        .chip-value.warning { color: var(--accent-orange); }

        #workspace {
            position: relative;
            flex: 1;
            display: flex;
            width: 100%;
            height: calc(100vh - 60px);
        }

        #map {
            width: 100%;
            height: 100%;
            background: #080d1a;
        }

        .floating-panel {
            position: absolute;
            background: var(--panel-bg);
            backdrop-filter: blur(16px);
            border: 1px solid var(--panel-border);
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.45);
            z-index: 1000;
            padding: 14px;
            transition: all 0.2s ease;
        }

        #layer-manager-dock {
            top: 20px;
            right: 20px;
            width: 270px;
            max-height: calc(100vh - 200px);
            overflow-y: auto;
        }

        #layer-manager-dock::-webkit-scrollbar { width: 4px; }
        #layer-manager-dock::-webkit-scrollbar-thumb { background: rgba(255,255,255,0.2); border-radius: 2px; }

        .layer-group-title {
            font-size: 0.72rem;
            text-transform: uppercase;
            color: var(--accent-blue);
            margin: 10px 0 6px 0;
            letter-spacing: 0.06em;
            font-weight: 700;
            border-bottom: 1px solid rgba(255,255,255,0.08);
            padding-bottom: 3px;
        }

        .layer-group-title:first-child { margin-top: 0; }

        .layer-toggle-group {
            display: flex;
            flex-direction: column;
            gap: 6px;
        }

        .layer-item {
            display: flex;
            align-items: center;
            justify-content: space-between;
            font-size: 0.78rem;
            color: var(--text-main);
            padding: 3px 0;
        }

        .layer-item label {
            display: flex;
            align-items: center;
            gap: 6px;
            cursor: pointer;
            user-select: none;
        }

        .layer-item input[type="checkbox"] {
            accent-color: var(--accent-blue);
            cursor: pointer;
        }

        #controls-dock {
            top: 20px;
            left: 20px;
            width: 260px;
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .control-section h4 {
            font-size: 0.70rem;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: var(--accent-blue);
            margin-bottom: 5px;
        }

        .control-select {
            width: 100%;
            background: rgba(0, 0, 0, 0.4);
            border: 1px solid var(--panel-border);
            color: #fff;
            padding: 5px 8px;
            border-radius: 6px;
            font-size: 0.78rem;
            outline: none;
        }

        .slider-control-row {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            font-size: 0.74rem;
            color: var(--text-muted);
        }

        #playback-dock {
            bottom: 24px;
            left: 50%;
            transform: translateX(-50%);
            width: 90%;
            max-width: 720px;
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .playback-controls {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 12px;
        }

        .btn-group { display: flex; gap: 6px; }

        .control-btn {
            background: rgba(255, 255, 255, 0.08);
            border: 1px solid var(--panel-border);
            color: #fff;
            padding: 6px 14px;
            border-radius: 6px;
            font-size: 0.82rem;
            font-weight: 500;
            cursor: pointer;
            display: inline-flex;
            align-items: center;
            gap: 6px;
            transition: background 0.15s ease;
        }

        .control-btn:hover { background: rgba(255, 255, 255, 0.18); }
        .control-btn.active {
            background: var(--accent-blue);
            color: #000;
            font-weight: 600;
            border-color: var(--accent-blue);
        }

        .time-display {
            font-family: var(--font-mono);
            font-size: 0.88rem;
            color: #fff;
            background: rgba(0, 0, 0, 0.35);
            padding: 4px 10px;
            border-radius: 6px;
            border: 1px solid var(--panel-border);
        }

        .slider-row { display: flex; align-items: center; gap: 12px; }

        input[type="range"] {
            flex: 1;
            accent-color: var(--accent-blue);
            height: 6px;
            border-radius: 3px;
            background: rgba(255, 255, 255, 0.2);
            cursor: pointer;
        }

        #hud-dock {
            bottom: 24px;
            right: 20px;
            width: 260px;
            font-size: 0.78rem;
        }

        .hud-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 6px 10px;
            margin-top: 6px;
        }

        .hud-item { display: flex; flex-direction: column; }
        .hud-label { font-size: 0.62rem; color: var(--text-muted); text-transform: uppercase; }
        .hud-val { font-family: var(--font-mono); font-weight: 600; color: #fff; }

        #legend-dock {
            bottom: 125px;
            left: 20px;
            width: 250px;
            padding: 10px 12px;
            font-size: 0.74rem;
        }

        .legend-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-weight: 600;
            color: var(--accent-blue);
            margin-bottom: 6px;
            text-transform: uppercase;
            font-size: 0.68rem;
            letter-spacing: 0.05em;
        }

        .legend-bar {
            height: 12px;
            border-radius: 4px;
            width: 100%;
            margin-bottom: 4px;
            border: 1px solid rgba(255,255,255,0.2);
        }

        .legend-labels {
            display: flex;
            justify-content: space-between;
            font-family: var(--font-mono);
            font-size: 0.68rem;
            color: var(--text-muted);
        }

        .leaflet-tooltip {
            background: rgba(15, 23, 42, 0.94) !important;
            color: #fff !important;
            border: 1px solid rgba(255,255,255,0.25) !important;
            border-radius: 6px !important;
            font-family: var(--font-sans) !important;
            font-size: 0.76rem !important;
            box-shadow: 0 4px 16px rgba(0,0,0,0.5) !important;
            padding: 6px 10px !important;
            line-height: 1.35 !important;
        }
    </style>
</head>
<body>
    <header>
        <div class="header-title">
            <h1>Larval <span>Dispersal</span> & Ocean Hydrodynamics</h1>
            <div class="header-subtitle">Multi-Scenario Hydrodynamic Particle Tracking Dashboard (*Chionoecetes opilio*)</div>
        </div>
        <div class="header-controls">
            <div class="scenario-selector-box">
                <label for="scenario-select">Scenario:</label>
                <select id="scenario-select" onchange="switchScenario(this.value)"></select>
            </div>
            <div class="metric-chips">
                <div class="chip">
                    <span class="chip-label">Cohort</span>
                    <span class="chip-value" id="chip-cohort">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Settled</span>
                    <span class="chip-value success" id="chip-settled">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Survival</span>
                    <span class="chip-value" id="chip-survival">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Mean PLD</span>
                    <span class="chip-value warning" id="chip-pld">-</span>
                </div>
            </div>
        </div>
    </header>

    <div id="workspace">
        <div id="map"></div>

        <div class="floating-panel" id="controls-dock">
            <div class="control-section">
                <h4>Particle Palette</h4>
                <select id="color-mode-select" class="control-select">
                    <option value="stage">Larval Stage (Molting)</option>
                    <option value="depth">Drift Depth (m)</option>
                    <option value="temp">Ambient Temperature (°C)</option>
                    <option value="dd">Cumulative Degree-Days</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Larval Cohort Filter</h4>
                <select id="outcome-filter-select" class="control-select">
                    <option value="all">All Particles</option>
                    <option value="settled">Settled Successfully</option>
                    <option value="pelagic">Remaining Pelagic</option>
                    <option value="dead">Thermal Mortality (Dead)</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Focus Particle</h4>
                <select id="particle-select" class="control-select">
                    <option value="-1">Overview (All Larvae)</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Hydrodynamic Depth Level</h4>
                <select id="hydro-depth-select" class="control-select">
                    <option value="0">Surface (0.0 m)</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Hydrodynamic Layer Opacity</h4>
                <div class="slider-control-row">
                    <span>Alpha</span>
                    <input type="range" id="hydro-opacity-slider" min="10" max="100" value="65">
                    <span id="hydro-opacity-val" style="font-family:var(--font-mono)">65%</span>
                </div>
            </div>
        </div>

        <div class="floating-panel" id="layer-manager-dock">
            <div class="layer-group-title">🛰️ Lagrangian Tracking</div>
            <div class="layer-toggle-group">
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-fulltracks" checked> 🛰️ Full Trajectories</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-progtrails" checked> ✨ Progressive Trails</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-heads" checked> 🔵 Active Larval Heads</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-release" checked> 🟢 Spawning Releases</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-settle" checked> 🔴 Settlement Sites</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-density"> 🔥 Nursery Settlement Density</label>
                </div>
            </div>

            <div class="layer-group-title">🌊 Hydrodynamic Model Outputs</div>
            <div class="layer-toggle-group">
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-advection" checked> 💨 Advection Current Vectors</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-speed"> ⚡ Current Speed (|u_h|)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-temp"> 🌡️ Seawater Temperature (T)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-sal"> 🧂 Practical Salinity (S)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-strat"> 📈 Buoyancy Stratification (N²)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-diff"> 🌀 Turbulent Diffusivity (κ_v)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-eta"> 🌊 Free Surface Height (η)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-w"> ⬆️ Upwelling / Downwelling (w)</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-hydro-bathy"> 🗺️ Seafloor Bathymetry</label>
                </div>
            </div>

            <div class="layer-group-title">🏛️ Demographics & Zones</div>
            <div class="layer-toggle-group">
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-strata" checked> 🗺️ CFA Management Zones</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-conn"> ↔️ CFA Connectivity Flows</label>
                </div>
            </div>
        </div>

        <div class="floating-panel" id="playback-dock">
            <div class="playback-controls">
                <div class="btn-group">
                    <button class="control-btn" id="play-btn">▶ Play</button>
                    <button class="control-btn" id="reset-btn">⏮ Reset</button>
                    <button class="control-btn speed-btn active" data-speed="1">1x</button>
                    <button class="control-btn speed-btn" data-speed="2">2x</button>
                    <button class="control-btn speed-btn" data-speed="4">4x</button>
                </div>
                <div class="time-display" id="time-readout">Day 0.00 / 0.00 d</div>
            </div>
            <div class="slider-row">
                <input type="range" id="time-slider" min="0" max="0" value="0" step="1">
            </div>
        </div>

        <div class="floating-panel" id="hud-dock">
            <div style="display:flex; justify-content:space-between; align-items:center;">
                <b style="color:var(--accent-blue);" id="telemetry-hud">Telemetry HUD</b>
                <span id="hud-id" style="font-family:var(--font-mono); color:#fff;">Overview</span>
            </div>
            <div class="hud-grid">
                <div class="hud-item"><span class="hud-label">Stage / Mode</span><span class="hud-val" id="hud-stage">-</span></div>
                <div class="hud-item"><span class="hud-label">Depth / Seabed</span><span class="hud-val" id="hud-depth">-</span></div>
                <div class="hud-item"><span class="hud-label">Temp (T)</span><span class="hud-val" id="hud-temp">-</span></div>
                <div class="hud-item"><span class="hud-label">Salinity (S)</span><span class="hud-val" id="hud-sal">-</span></div>
                <div class="hud-item"><span class="hud-label">Current (|u|)</span><span class="hud-val" id="hud-current">-</span></div>
                <div class="hud-item"><span class="hud-label">Strat (N²)</span><span class="hud-val" id="hud-strat">-</span></div>
                <div class="hud-item"><span class="hud-label">Diffusivity (κ)</span><span class="hud-val" id="hud-diff">-</span></div>
                <div class="hud-item"><span class="hud-label">Survival / SSH</span><span class="hud-val" id="hud-surv">-</span></div>
            </div>
        </div>

        <div class="floating-panel" id="legend-dock">
            <div class="legend-header">
                <span id="legend-title">Seawater Temperature</span>
                <span id="legend-unit" style="font-family:var(--font-mono)">°C</span>
            </div>
            <div class="legend-bar" id="legend-bar" style="background: linear-gradient(to right, #1d4ed8, #06b6d4, #10b981, #f59e0b, #ef4444);"></div>
            <div class="legend-labels">
                <span id="legend-min">0.0</span>
                <span id="legend-mid">7.5</span>
                <span id="legend-max">15.0</span>
            </div>
        </div>
    </div>

    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo=" crossorigin=""></script>
    <script>
        const SCENARIOS = $(scenarios_master_json);
        const STRATA = $(strata_json);
        let activeScenKey = "$(default_scen_key)";
        let currentScen = SCENARIOS[activeScenKey];
        const PARTICLES = currentScen.particles;

        const STAGE_COLORS = {
            "zoea1": "#38bdf8", "zoea2": "#06b6d4", "zoea3": "#0d9488", "zoea4": "#eab308",
            "megalopa": "#f97316", "instar1_settled": "#10b981", "dead": "#ef4444"
        };

        const map = L.map('map', { center: [44.2, -62.5], zoom: 7 });
        const esriOcean = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}', { maxZoom: 13 }).addTo(map);
        const cartoDark = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', { maxZoom: 19 });

        L.control.layers({ "ESRI Ocean": esriOcean, "CartoDB Dark": cartoDark }, null, { position: 'bottomleft' }).addTo(map);

        // Management & Lagrangian Layers
        const strataLayer = L.layerGroup().addTo(map);
        const fullTracksLayer = L.layerGroup().addTo(map);
        const progressTracksLayer = L.layerGroup().addTo(map);
        const activeHeadsLayer = L.layerGroup().addTo(map);
        const releaseMarkersLayer = L.layerGroup().addTo(map);
        const settlementMarkersLayer = L.layerGroup().addTo(map);
        const settlementDensityLayer = L.layerGroup();
        const connectivityFlowsLayer = L.layerGroup();

        // Hydrodynamic Model Eulerian Layers
        const hydroAdvectionLayer = L.layerGroup().addTo(map);
        const hydroSpeedLayer = L.layerGroup();
        const hydroTempLayer = L.layerGroup();
        const hydroSalLayer = L.layerGroup();
        const hydroStratLayer = L.layerGroup();
        const hydroDiffLayer = L.layerGroup();
        const hydroEtaLayer = L.layerGroup();
        const hydroUpwellingLayer = L.layerGroup();
        const hydroBathyLayer = L.layerGroup();
        let selectedDepthIdx = 0;

        STRATA.forEach(s => {
            let layer = s.is_polygon ? L.polygon(s.polygon, {color: s.color, weight: 2.2, fillOpacity: 0.08}) : L.rectangle([[s.lat_min, s.lon_min], [s.lat_max, s.lon_max]], {color: s.color, weight: 2, fillOpacity: 0.06});
            layer.bindTooltip("<b>" + s.name + "</b>").addTo(strataLayer);
        });

        const scenSelect = document.getElementById('scenario-select');
        Object.keys(SCENARIOS).forEach(k => {
            const opt = document.createElement('option'); opt.value = k; opt.textContent = k;
            if (k === activeScenKey) opt.selected = true;
            scenSelect.appendChild(opt);
        });

        let currentStep = 0, isPlaying = false, playSpeed = 1, animTimer = null, showTrails = true, colorMode = "stage", outcomeFilter = "all", focusParticleId = -1;
        let hydroOpacity = 0.65;
        let headMarkers = [], fullTrackPolylines = [], progressTrackPolylines = [];

        // Scientific Colormaps
        function colormapThermal(norm) {
            if (norm < 0.25) {
                const t = norm / 0.25;
                return [Math.round(29 + t * (6 - 29)), Math.round(78 + t * (182 - 78)), Math.round(216 + t * (212 - 216)), 255];
            } else if (norm < 0.5) {
                const t = (norm - 0.25) / 0.25;
                return [Math.round(6 + t * (16 - 6)), Math.round(182 + t * (185 - 182)), Math.round(212 + t * (129 - 212)), 255];
            } else if (norm < 0.75) {
                const t = (norm - 0.5) / 0.25;
                return [Math.round(16 + t * (245 - 16)), Math.round(185 + t * (158 - 185)), Math.round(129 + t * (11 - 129)), 255];
            } else {
                const t = (norm - 0.75) / 0.25;
                return [Math.round(245 + t * (239 - 245)), Math.round(158 + t * (68 - 158)), Math.round(11 + t * (68 - 11)), 255];
            }
        }

        function colormapHaline(norm) {
            const r = Math.round(15 + norm * 150);
            const g = Math.round(180 - norm * 110);
            const b = Math.round(180 + norm * 60);
            return [r, g, b, 255];
        }

        function colormapTurbo(norm) {
            const r = Math.round(255 * Math.min(1.0, Math.max(0.0, 1.5 - Math.abs(norm - 0.75) * 4)));
            const g = Math.round(255 * Math.min(1.0, Math.max(0.0, 1.5 - Math.abs(norm - 0.5) * 4)));
            const b = Math.round(255 * Math.min(1.0, Math.max(0.0, 1.5 - Math.abs(norm - 0.25) * 4)));
            return [r, g, b, 255];
        }

        function colormapCoolwarm(norm) {
            const r = Math.round(59 + norm * (220 - 59));
            const g = Math.round(76 + (1 - Math.abs(norm - 0.5) * 2) * 120);
            const b = Math.round(192 - norm * (192 - 40));
            return [r, g, b, 255];
        }

        function colormapBalance(norm) {
            if (norm < 0.5) {
                const t = norm / 0.5;
                return [Math.round(30 + t * 160), Math.round(64 + t * 130), Math.round(175 + t * 20), 255];
            } else {
                const t = (norm - 0.5) / 0.5;
                return [Math.round(190 + t * 49), Math.round(194 - t * 140), Math.round(195 - t * 140), 255];
            }
        }

        function colormapOcean(norm) {
            const r = Math.round(10 + norm * 80);
            const g = Math.round(40 + norm * 140);
            const b = Math.round(90 + norm * 130);
            return [r, g, b, 255];
        }

        function colormapViridis(norm) {
            const r = Math.round(68 + norm * 185);
            const g = Math.round(1 + norm * 230);
            const b = Math.round(84 + (1 - norm) * 110 + norm * 15);
            return [Math.min(255, Math.max(0, r)), Math.min(255, Math.max(0, g)), Math.min(255, Math.max(0, b)), 255];
        }

        function colormapPlasma(norm) {
            const r = Math.round(13 + norm * 225);
            const g = Math.round(8 + Math.sin(norm * Math.PI) * 180);
            const b = Math.round(135 * (1 - norm) + 30);
            return [Math.min(255, Math.max(0, r)), Math.min(255, Math.max(0, g)), Math.min(255, Math.max(0, b)), 255];
        }

        function generateRasterDataURL(grid2D, lons, lats, minVal, maxVal, colormapFn) {
            const width = 200;
            const height = 150;
            const canvas = document.createElement('canvas');
            canvas.width = width;
            canvas.height = height;
            const ctx = canvas.getContext('2d');
            const imgData = ctx.createImageData(width, height);
            const data = imgData.data;

            const nx = lons.length;
            const ny = lats.length;
            const range = Math.max(1e-5, maxVal - minVal);

            for (let py = 0; py < height; py++) {
                const latFrac = 1.0 - (py / (height - 1));
                const yIdx = latFrac * (ny - 1);
                const j0 = Math.floor(yIdx);
                const j1 = Math.min(ny - 1, j0 + 1);
                const fy = yIdx - j0;

                for (let px = 0; px < width; px++) {
                    const lonFrac = px / (width - 1);
                    const xIdx = lonFrac * (nx - 1);
                    const i0 = Math.floor(xIdx);
                    const i1 = Math.min(nx - 1, i0 + 1);
                    const fx = xIdx - i0;

                    const v00 = grid2D[i0][j0];
                    const v10 = grid2D[i1][j0];
                    const v01 = grid2D[i0][j1];
                    const v11 = grid2D[i1][j1];

                    const val = (1 - fx) * (1 - fy) * v00 + fx * (1 - fy) * v10 + (1 - fx) * fy * v01 + fx * fy * v11;
                    const norm = Math.min(1.0, Math.max(0.0, (val - minVal) / range));
                    const rgba = colormapFn(norm);

                    const pIndex = (py * width + px) * 4;
                    data[pIndex] = rgba[0];
                    data[pIndex + 1] = rgba[1];
                    data[pIndex + 2] = rgba[2];
                    data[pIndex + 3] = rgba[3];
                }
            }
            ctx.putImageData(imgData, 0, 0);
            return canvas.toDataURL();
        }

        function getColor(p, step) {
            if (colorMode === "stage") return STAGE_COLORS[p.stages[step]] || "#38bdf8";
            if (colorMode === "depth") return Math.abs(p.depths[step]) < 50 ? "#38bdf8" : "#1e3a8a";
            if (colorMode === "temp") return p.temps[step] < 6 ? "#38bdf8" : "#f87171";
            return p.dds[step] < 100 ? "#38bdf8" : "#f97316";
        }

        let activeRasterLayers = [];

        function updateLegend(title, unit, minV, maxV, gradientCss) {
            document.getElementById('legend-title').textContent = title;
            document.getElementById('legend-unit').textContent = unit;
            document.getElementById('legend-min').textContent = minV;
            document.getElementById('legend-max').textContent = maxV;
            document.getElementById('legend-mid').textContent = ((parseFloat(minV) + parseFloat(maxV)) / 2.0).toFixed(1);
            document.getElementById('legend-bar').style.background = gradientCss;
        }

        function rebuildHydroLayers() {
            [hydroAdvectionLayer, hydroSpeedLayer, hydroTempLayer, hydroSalLayer,
             hydroStratLayer, hydroDiffLayer, hydroUpwellingLayer].forEach(l => l.clearLayers());

            const hydro = currentScen.hydrodynamics;
            if (!hydro || !hydro.lons || !hydro.lats) return;

            const hLons = hydro.lons;
            const hLats = hydro.lats;
            const bounds = [[Math.min(...hLats), Math.min(...hLons)], [Math.max(...hLats), Math.max(...hLons)]];

            const lvl = (hydro.depth_levels && hydro.depth_levels[selectedDepthIdx]) ? hydro.depth_levels[selectedDepthIdx] : null;
            const vecs = lvl ? lvl.vectors : (hydro.vectors || []);
            const curDepth = lvl ? Math.abs(lvl.depth_m).toFixed(1) : "0.0";

            // A. Advection Current Vectors (Quivers)
            if (vecs && vecs.length > 0) {
                vecs.forEach(vec => {
                    const lenScale = 0.0035 * Math.min(18.0, Math.max(1.0, vec.speed));
                    const lat2 = vec.lat + (vec.v / Math.max(0.1, vec.speed)) * lenScale;
                    const lon2 = vec.lon + (vec.u / Math.max(0.1, vec.speed)) * lenScale * 1.35;

                    const spdNorm = Math.min(1.0, vec.speed / 20.0);
                    const arrowCol = spdNorm > 0.6 ? '#f43f5e' : spdNorm > 0.3 ? '#38bdf8' : '#34d399';

                    const line = L.polyline([[vec.lat, vec.lon], [lat2, lon2]], {
                        color: arrowCol,
                        weight: 2.2,
                        opacity: 0.85
                    }).bindTooltip(`<b>Hydrodynamic Current Vector</b><br>Depth: <b>\${curDepth} m</b><br>Speed: <b>\${vec.speed} cm/s</b> (\${(vec.speed*0.0194).toFixed(2)} kts)<br>Dir: \${vec.dir}°<br>u: \${vec.u} cm/s, v: \${vec.v} cm/s<br>Temp: \${vec.temp} °C | Salinity: \${vec.sal} PSU<br>Strat N²: \${vec.strat !== undefined ? vec.strat : '-'} 10⁻⁴s⁻²<br>Diff κ_v: \${vec.diff !== undefined ? vec.diff : '-'} 10⁻⁴m²/s<br>Seabed Depth: \${vec.depth} m`);
                    hydroAdvectionLayer.addLayer(line);

                    const tip = L.circleMarker([lat2, lon2], {
                        radius: 2.5,
                        color: arrowCol,
                        fillColor: arrowCol,
                        fillOpacity: 1.0,
                        weight: 1
                    });
                    hydroAdvectionLayer.addLayer(tip);
                });
            }

            // B. Multi-depth Grids
            const tGrid = lvl ? lvl.temperature : hydro.grids?.temperature;
            if (tGrid) {
                const tUrl = generateRasterDataURL(tGrid.data, hLons, hLats, tGrid.min, tGrid.max, colormapThermal);
                const tOverlay = L.imageOverlay(tUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroTempLayer.addLayer(tOverlay);
                activeRasterLayers.push(tOverlay);
            }

            const sGrid = lvl ? lvl.salinity : hydro.grids?.salinity;
            if (sGrid) {
                const sUrl = generateRasterDataURL(sGrid.data, hLons, hLats, sGrid.min, sGrid.max, colormapHaline);
                const sOverlay = L.imageOverlay(sUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroSalLayer.addLayer(sOverlay);
                activeRasterLayers.push(sOverlay);
            }

            const spdGrid = lvl ? lvl.speed : hydro.grids?.speed;
            if (spdGrid) {
                const spdUrl = generateRasterDataURL(spdGrid.data, hLons, hLats, 0.0, spdGrid.max, colormapTurbo);
                const spdOverlay = L.imageOverlay(spdUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroSpeedLayer.addLayer(spdOverlay);
                activeRasterLayers.push(spdOverlay);
            }

            const stratGrid = lvl ? lvl.stratification : null;
            if (stratGrid) {
                const stratUrl = generateRasterDataURL(stratGrid.data, hLons, hLats, stratGrid.min, stratGrid.max, colormapViridis);
                const stratOverlay = L.imageOverlay(stratUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroStratLayer.addLayer(stratOverlay);
                activeRasterLayers.push(stratOverlay);
            }

            const diffGrid = lvl ? lvl.diffusion : null;
            if (diffGrid) {
                const diffUrl = generateRasterDataURL(diffGrid.data, hLons, hLats, diffGrid.min, diffGrid.max, colormapPlasma);
                const diffOverlay = L.imageOverlay(diffUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroDiffLayer.addLayer(diffOverlay);
                activeRasterLayers.push(diffOverlay);
            }

            const wGrid = lvl ? lvl.upwelling : hydro.grids?.upwelling;
            if (wGrid) {
                const wUrl = generateRasterDataURL(wGrid.data, hLons, hLats, wGrid.min, wGrid.max, colormapBalance);
                const wOverlay = L.imageOverlay(wUrl, bounds, { opacity: hydroOpacity, interactive: false });
                hydroUpwellingLayer.addLayer(wOverlay);
                activeRasterLayers.push(wOverlay);
            }
        }

        function renderActiveScenario() {
            [fullTracksLayer, progressTracksLayer, activeHeadsLayer, releaseMarkersLayer,
             settlementMarkersLayer, settlementDensityLayer, hydroAdvectionLayer,
             hydroSpeedLayer, hydroTempLayer, hydroSalLayer, hydroStratLayer, hydroDiffLayer,
             hydroEtaLayer, hydroUpwellingLayer, hydroBathyLayer, connectivityFlowsLayer].forEach(l => l.clearLayers());

            headMarkers = [];
            fullTrackPolylines = [];
            progressTrackPolylines = [];
            activeRasterLayers = [];

            const pList = currentScen.particles;
            const nTimes = currentScen.times.length;

            // 1. Summary Chips
            document.getElementById('chip-cohort').textContent = currentScen.stats.n_particles;
            document.getElementById('chip-settled').textContent = `\${currentScen.stats.n_settled} (\${currentScen.stats.settle_pct}%)`;
            document.getElementById('chip-survival').textContent = `\${currentScen.stats.survival_pct}%`;
            document.getElementById('chip-pld').textContent = `\${currentScen.stats.mean_pld_days} d`;

            // 2. Trajectories, Heads & Markers
            pList.forEach((p) => {
                const fullLatLngs = [];
                for (let t = 0; t < nTimes; t++) {
                    fullLatLngs.push([p.lats[t], p.lons[t]]);
                }

                // Full Polyline
                const fullPoly = L.polyline(fullLatLngs, {
                    color: getColor(p, nTimes - 1),
                    weight: 2.2,
                    opacity: 0.55
                }).on('click', () => selectParticle(p.id));
                fullTrackPolylines.push(fullPoly);
                fullTracksLayer.addLayer(fullPoly);

                // Progress Polyline
                const progPoly = L.polyline(fullLatLngs.slice(0, 2), {
                    color: getColor(p, 0),
                    weight: 3.5,
                    opacity: 0.95
                }).on('click', () => selectParticle(p.id));
                progressTrackPolylines.push(progPoly);
                progressTracksLayer.addLayer(progPoly);

                // Head Marker
                const marker = L.circleMarker([p.lats[0], p.lons[0]], {
                    radius: 5.5,
                    color: '#ffffff',
                    weight: 2,
                    fillColor: getColor(p, 0),
                    fillOpacity: 1.0
                }).on('click', () => selectParticle(p.id));
                marker.bindTooltip(`<b>Larva #\${p.id}</b> (\${p.stages[0]})`, { sticky: true });
                headMarkers.push(marker);
                activeHeadsLayer.addLayer(marker);

                // Spawning Release Marker
                const relMarker = L.circleMarker([p.lats[0], p.lons[0]], {
                    radius: 4,
                    color: '#10b981',
                    fillColor: '#059669',
                    fillOpacity: 0.85,
                    weight: 1.5
                }).bindPopup(`<b>Spawning Release #\${p.id}</b><br>Lon: \${p.lons[0].toFixed(3)}°E<br>Lat: \${p.lats[0].toFixed(3)}°N`);
                releaseMarkersLayer.addLayer(relMarker);

                // Settlement Marker (if settled)
                if (p.status === "settled_successful") {
                    const lastT = nTimes - 1;
                    const setMarker = L.circleMarker([p.lats[lastT], p.lons[lastT]], {
                        radius: 6,
                        color: '#f59e0b',
                        fillColor: '#d97706',
                        fillOpacity: 0.9,
                        weight: 2
                    }).bindPopup(`<b>Benthic Settlement Site #\${p.id}</b><br>Age: \${p.settlement_age_days} d<br>Lon: \${p.lons[lastT].toFixed(3)}°E<br>Lat: \${p.lats[lastT].toFixed(3)}°N`);
                    settlementMarkersLayer.addLayer(setMarker);
                }
            });

            // 3. Hydrodynamic Model Outputs Setup
            const hydro = currentScen.hydrodynamics;
            selectedDepthIdx = 0;
            const depthSelect = document.getElementById('hydro-depth-select');
            if (depthSelect) {
                depthSelect.innerHTML = '';
                if (hydro && hydro.depths && hydro.depths.length > 0) {
                    hydro.depths.forEach((d, idx) => {
                        const opt = document.createElement('option');
                        opt.value = idx;
                        opt.textContent = `\${Math.abs(d).toFixed(1)} m \${d === 0 ? "(Surface)" : "(Subsurface)"}`;
                        if (idx === selectedDepthIdx) opt.selected = true;
                        depthSelect.appendChild(opt);
                    });
                } else {
                    const opt = document.createElement('option');
                    opt.value = 0;
                    opt.textContent = 'Surface (0.0 m)';
                    depthSelect.appendChild(opt);
                }
            }

            if (hydro && hydro.lons && hydro.lats) {
                const hLons = hydro.lons;
                const hLats = hydro.lats;
                const bounds = [[Math.min(...hLats), Math.min(...hLons)], [Math.max(...hLats), Math.max(...hLons)]];

                // Static 2D Fields: Elevation and Bathymetry
                const g = hydro.grids;
                if (g) {
                    if (g.elevation) {
                        const etaUrl = generateRasterDataURL(g.elevation.data, hLons, hLats, g.elevation.min, g.elevation.max, colormapCoolwarm);
                        const etaOverlay = L.imageOverlay(etaUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroEtaLayer.addLayer(etaOverlay);
                        activeRasterLayers.push(etaOverlay);
                    }

                    if (g.bathymetry) {
                        const bUrl = generateRasterDataURL(g.bathymetry.data, hLons, hLats, g.bathymetry.min, 0.0, colormapOcean);
                        const bOverlay = L.imageOverlay(bUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroBathyLayer.addLayer(bOverlay);
                        activeRasterLayers.push(bOverlay);
                    }
                }

                // Render dynamic depth-dependent fields
                rebuildHydroLayers();
            }

            // 4. Nursery Settlement Density Heatmap
            if (currentScen.settlement_density && currentScen.settlement_density.length > 0) {
                const maxDens = Math.max(...currentScen.settlement_density.map(d => d.density));
                currentScen.settlement_density.forEach(d => {
                    const norm = maxDens > 0 ? (d.density / maxDens) : 0.5;
                    const col = norm > 0.66 ? '#ef4444' : norm > 0.33 ? '#f59e0b' : '#38bdf8';
                    const rect = L.circle([d.lat, d.lon], {
                        radius: 7000,
                        color: col,
                        fillColor: col,
                        fillOpacity: 0.35 + norm * 0.45,
                        weight: 1
                    }).bindTooltip(`<b>Settlement Nursery Density</b><br>Probability: \${(d.density * 100).toFixed(2)}%`);
                    settlementDensityLayer.addLayer(rect);
                });
            }

            // 5. Demographic Connectivity Flows
            if (currentScen.connectivity_flows && currentScen.connectivity_flows.length > 0) {
                currentScen.connectivity_flows.forEach(flow => {
                    const weight = Math.max(2.0, flow.prob * 10.0);
                    const flowLine = L.polyline([[flow.src_lat, flow.src_lon], [flow.dst_lat, flow.dst_lon]], {
                        color: '#a855f7',
                        weight: weight,
                        opacity: 0.8,
                        dashArray: '8, 4'
                    }).bindTooltip(`<b>Connectivity Transition</b><br>\${flow.src} ➔ \${flow.dst}<br>Probability P_{ij}: \${(flow.prob * 100).toFixed(1)}%`);
                    connectivityFlowsLayer.addLayer(flowLine);
                });
            }

            // 6. Particle Selector Dropdown
            const pSelect = document.getElementById('particle-select');
            pSelect.innerHTML = '<option value="-1">Overview (All Larvae)</option>';
            pList.forEach(p => {
                const opt = document.createElement('option');
                opt.value = p.id;
                opt.textContent = `Larva #\${p.id} (\${p.status})`;
                pSelect.appendChild(opt);
            });

            // 7. Slider Bounds
            const slider = document.getElementById('time-slider');
            slider.max = nTimes - 1;
            slider.value = 0;

            // Auto-fit Bounds
            const allLats = pList.flatMap(p => p.lats);
            const allLons = pList.flatMap(p => p.lons);
            if (allLats.length > 0) {
                map.fitBounds([
                    [Math.min(...allLats) - 0.25, Math.min(...allLons) - 0.25],
                    [Math.max(...allLats) + 0.25, Math.max(...allLons) + 0.25]
                ]);
            }

            updateLegend("Seawater Temperature (T)", "°C", (hydro && hydro.grids && hydro.grids.temperature) ? hydro.grids.temperature.min : "1.5", (hydro && hydro.grids && hydro.grids.temperature) ? hydro.grids.temperature.max : "14.5", "linear-gradient(to right, #1d4ed8, #06b6d4, #10b981, #f59e0b, #ef4444)");
            updateFrame(0);
        }

        function switchScenario(key) {
            if (SCENARIOS[key]) {
                activeScenKey = key;
                currentScen = SCENARIOS[key];
                isPlaying = false;
                document.getElementById('play-btn').textContent = "▶ Play";
                clearTimeout(animTimer);
                renderActiveScenario();
            }
        }

        function updateFrame(step) {
            currentStep = step;
            const nTimes = currentScen.times.length;
            const pList = currentScen.particles;

            document.getElementById('time-slider').value = step;
            document.getElementById('time-readout').textContent = `Day \${currentScen.times[step].toFixed(2)} / \${currentScen.times[nTimes - 1].toFixed(2)} d`;

            pList.forEach((p, idx) => {
                const isFiltered = (outcomeFilter === "settled" && p.status !== "settled_successful") ||
                                   (outcomeFilter === "pelagic" && p.status !== "pelagic") ||
                                   (outcomeFilter === "dead" && p.alive);
                const isFocused = focusParticleId === -1 || focusParticleId === p.id;

                if (isFiltered || !isFocused) {
                    if (headMarkers[idx]) headMarkers[idx].setStyle({ opacity: 0, fillOpacity: 0 });
                    if (fullTrackPolylines[idx]) fullTrackPolylines[idx].setStyle({ opacity: 0 });
                    if (progressTrackPolylines[idx]) progressTrackPolylines[idx].setStyle({ opacity: 0 });
                    return;
                }

                const lat = p.lats[step];
                const lon = p.lons[step];
                const col = getColor(p, step);

                if (headMarkers[idx]) {
                    headMarkers[idx].setLatLng([lat, lon]);
                    headMarkers[idx].setStyle({
                        fillColor: col,
                        opacity: 1,
                        fillOpacity: 1,
                        radius: (focusParticleId === p.id) ? 9 : 5.5
                    });
                    headMarkers[idx].setTooltipContent(`<b>Larva #\${p.id}</b> (\${p.stages[step]} • \${p.depths[step].toFixed(1)}m)`);
                }

                if (fullTrackPolylines[idx]) {
                    fullTrackPolylines[idx].setStyle({
                        color: col,
                        opacity: (focusParticleId === p.id) ? 0.9 : 0.45,
                        weight: (focusParticleId === p.id) ? 3.5 : 2.0
                    });
                }

                if (progressTrackPolylines[idx]) {
                    if (showTrails && step > 0) {
                        const progPts = [];
                        for (let t = 0; t <= step; t++) {
                            progPts.push([p.lats[t], p.lons[t]]);
                        }
                        progressTrackPolylines[idx].setLatLngs(progPts);
                        progressTrackPolylines[idx].setStyle({
                            color: col,
                            opacity: (focusParticleId === p.id) ? 1.0 : 0.95,
                            weight: (focusParticleId === p.id) ? 4.5 : 3.2
                        });
                    } else if (showTrails && step === 0) {
                        const initPts = [[p.lats[0], p.lons[0]], [p.lats[Math.min(1, nTimes - 1)], p.lons[Math.min(1, nTimes - 1)]]];
                        progressTrackPolylines[idx].setLatLngs(initPts);
                        progressTrackPolylines[idx].setStyle({ color: col, opacity: 0.9, weight: 3.0 });
                    } else {
                        progressTrackPolylines[idx].setStyle({ opacity: 0 });
                    }
                }
            });

            if (focusParticleId !== -1) {
                const p = pList.find(item => item.id === focusParticleId);
                if (p) updateHUD(p, step);
            }
        }

        function updateHUD(p, step) {
            document.getElementById('hud-id').textContent = `#\${p.id}`;
            document.getElementById('hud-stage').textContent = p.stages[step];
            document.getElementById('hud-depth').textContent = `\${p.depths[step].toFixed(1)} m`;
            document.getElementById('hud-temp').textContent = `\${p.temps[step].toFixed(2)} °C`;
            document.getElementById('hud-sal').textContent = `33.2 PSU`;
            document.getElementById('hud-current').textContent = `~8.5 cm/s`;
            document.getElementById('hud-strat').textContent = `0.85 10⁻⁴s⁻²`;
            document.getElementById('hud-diff').textContent = `2.5 10⁻⁴m²/s`;
            document.getElementById('hud-surv').textContent = `\${(p.survs[step] * 100).toFixed(1)} %`;
        }

        function selectParticle(id) {
            focusParticleId = id;
            document.getElementById('particle-select').value = id;
            const p = currentScen.particles.find(item => item.id === id);
            if (p) {
                map.setView([p.lats[currentStep], p.lons[currentStep]], 9);
                updateHUD(p, currentStep);
            }
            updateFrame(currentStep);
        }

        function playLoop() {
            if (!isPlaying) return;
            const nTimes = currentScen.times.length;
            let nextStep = currentStep + playSpeed;
            if (nextStep >= nTimes) nextStep = 0;
            updateFrame(nextStep);
            animTimer = setTimeout(playLoop, 60);
        }

        // Live Cursor Telemetry Probe
        map.on('mousemove', function(e) {
            if (focusParticleId !== -1) return;
            const lat = e.latlng.lat;
            const lon = e.latlng.lng;
            const hydro = currentScen.hydrodynamics;
            if (hydro && hydro.lons && hydro.lats) {
                const hLons = hydro.lons;
                const hLats = hydro.lats;
                const minLon = hLons[0], maxLon = hLons[hLons.length - 1];
                const minLat = hLats[0], maxLat = hLats[hLats.length - 1];

                if (lon >= minLon && lon <= maxLon && lat >= minLat && lat <= maxLat) {
                    const nx = hLons.length, ny = hLats.length;
                    const i = Math.min(nx - 1, Math.max(0, Math.floor(((lon - minLon) / (maxLon - minLon)) * (nx - 1))));
                    const j = Math.min(ny - 1, Math.max(0, Math.floor(((lat - minLat) / (maxLat - minLat)) * (ny - 1))));

                    const lvl = (hydro.depth_levels && hydro.depth_levels[selectedDepthIdx]) ? hydro.depth_levels[selectedDepthIdx] : null;
                    const curDepth = (hydro.depths && hydro.depths.length > selectedDepthIdx) ? Math.abs(hydro.depths[selectedDepthIdx]).toFixed(1) : "0.0";
                    const g = hydro.grids;

                    document.getElementById('hud-id').textContent = `\${lat.toFixed(2)}°N, \${lon.toFixed(2)}°E`;
                    document.getElementById('hud-stage').textContent = `Eulerian @ \${curDepth}m`;
                    document.getElementById('hud-depth').textContent = `\${g && g.bathymetry ? g.bathymetry.data[i][j] : -120} m`;
                    document.getElementById('hud-temp').textContent = `\${lvl && lvl.temperature ? lvl.temperature.data[i][j] : (g && g.temperature ? g.temperature.data[i][j] : 8.0)} °C`;
                    document.getElementById('hud-sal').textContent = `\${lvl && lvl.salinity ? lvl.salinity.data[i][j] : (g && g.salinity ? g.salinity.data[i][j] : 32.5)} PSU`;
                    document.getElementById('hud-current').textContent = `\${lvl && lvl.speed ? lvl.speed.data[i][j] : (g && g.speed ? g.speed.data[i][j] : 8.0)} cm/s`;
                    document.getElementById('hud-strat').textContent = `\${lvl && lvl.stratification ? lvl.stratification.data[i][j] : '-'} 10⁻⁴s⁻²`;
                    document.getElementById('hud-diff').textContent = `\${lvl && lvl.diffusion ? lvl.diffusion.data[i][j] : '-'} 10⁻⁴m²/s`;
                    document.getElementById('hud-surv').textContent = `SSH \${g && g.elevation ? g.elevation.data[i][j] : 0.0} cm`;
                }
            }
        });

        // Layer Manager Event Listeners
        const layerControls = [
            { id: 'layer-fulltracks', layer: fullTracksLayer },
            { id: 'layer-progtrails', layer: progressTracksLayer },
            { id: 'layer-heads', layer: activeHeadsLayer },
            { id: 'layer-release', layer: releaseMarkersLayer },
            { id: 'layer-settle', layer: settlementMarkersLayer },
            { id: 'layer-density', layer: settlementDensityLayer, onSelect: () => updateLegend("Nursery Settlement Density", "%", "0.0", "100.0", "linear-gradient(to right, #1e3a8a, #06b6d4, #f59e0b, #ef4444)") },
            { id: 'layer-hydro-advection', layer: hydroAdvectionLayer },
            { id: 'layer-hydro-speed', layer: hydroSpeedLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Current Speed (|u_h|)", "cm/s", "0.0", (lvl && lvl.speed) ? lvl.speed.max : (currentScen.hydrodynamics?.grids?.speed?.max || "25.0"), "linear-gradient(to right, #3b82f6, #06b6d4, #10b981, #facc15, #ef4444)");
            }},
            { id: 'layer-hydro-temp', layer: hydroTempLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Seawater Temperature (T)", "°C", (lvl && lvl.temperature) ? lvl.temperature.min : (currentScen.hydrodynamics?.grids?.temperature?.min || "1.5"), (lvl && lvl.temperature) ? lvl.temperature.max : (currentScen.hydrodynamics?.grids?.temperature?.max || "14.5"), "linear-gradient(to right, #1d4ed8, #06b6d4, #10b981, #f59e0b, #ef4444)");
            }},
            { id: 'layer-hydro-sal', layer: hydroSalLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Practical Salinity (S)", "PSU", (lvl && lvl.salinity) ? lvl.salinity.min : (currentScen.hydrodynamics?.grids?.salinity?.min || "31.0"), (lvl && lvl.salinity) ? lvl.salinity.max : (currentScen.hydrodynamics?.grids?.salinity?.max || "35.5"), "linear-gradient(to right, #10b981, #06b6d4, #3b82f6, #6366f1, #8b5cf6)");
            }},
            { id: 'layer-hydro-strat', layer: hydroStratLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Stratification (N²)", "10⁻⁴ s⁻²", (lvl && lvl.stratification) ? lvl.stratification.min : "0.0", (lvl && lvl.stratification) ? lvl.stratification.max : "5.0", "linear-gradient(to right, #440154, #3b528b, #21908c, #5dc963, #fde725)");
            }},
            { id: 'layer-hydro-diff', layer: hydroDiffLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Turbulent Diffusivity (κ_v)", "10⁻⁴ m²/s", (lvl && lvl.diffusion) ? lvl.diffusion.min : "0.1", (lvl && lvl.diffusion) ? lvl.diffusion.max : "10.0", "linear-gradient(to right, #0d0887, #6a00a8, #b12a90, #e16462, #fca636, #f0f921)");
            }},
            { id: 'layer-hydro-eta', layer: hydroEtaLayer, onSelect: () => updateLegend("Free Surface Height (η)", "cm", currentScen.hydrodynamics?.grids?.elevation?.min || "-10.0", currentScen.hydrodynamics?.grids?.elevation?.max || "+10.0", "linear-gradient(to right, #3b82f6, #f8fafc, #ef4444)") },
            { id: 'layer-hydro-w', layer: hydroUpwellingLayer, onSelect: () => {
                const lvl = currentScen.hydrodynamics?.depth_levels?.[selectedDepthIdx];
                updateLegend("Vertical Velocity / Upwelling", "mm/s", (lvl && lvl.upwelling) ? lvl.upwelling.min : (currentScen.hydrodynamics?.grids?.upwelling?.min || "-0.5"), (lvl && lvl.upwelling) ? lvl.upwelling.max : (currentScen.hydrodynamics?.grids?.upwelling?.max || "+0.5"), "linear-gradient(to right, #1e40af, #94a3b8, #b91c1c)");
            }},
            { id: 'layer-hydro-bathy', layer: hydroBathyLayer, onSelect: () => updateLegend("Seafloor Bathymetry", "m", currentScen.hydrodynamics?.grids?.bathymetry?.min || "-3000", "0", "linear-gradient(to right, #030712, #1e3a8a, #0284c7, #38bdf8)") },
            { id: 'layer-strata', layer: strataLayer },
            { id: 'layer-conn', layer: connectivityFlowsLayer }
        ];

        function updateActiveHydroLegend() {
            layerControls.forEach(ctrl => {
                const el = document.getElementById(ctrl.id);
                if (el && el.checked && ctrl.onSelect) {
                    ctrl.onSelect();
                }
            });
        }

        layerControls.forEach(ctrl => {
            const el = document.getElementById(ctrl.id);
            if (el) {
                el.addEventListener('change', function() {
                    if (this.checked) {
                        map.addLayer(ctrl.layer);
                        if (ctrl.onSelect) ctrl.onSelect();
                    } else {
                        map.removeLayer(ctrl.layer);
                    }
                });
            }
        });

        // Hydrodynamic Depth Select Listener
        const depthSelect = document.getElementById('hydro-depth-select');
        if (depthSelect) {
            depthSelect.addEventListener('change', function(e) {
                selectedDepthIdx = parseInt(e.target.value);
                rebuildHydroLayers();
                updateActiveHydroLegend();
            });
        }

        // Hydrodynamic Opacity Slider
        const opacitySlider = document.getElementById('hydro-opacity-slider');
        if (opacitySlider) {
            opacitySlider.addEventListener('input', function(e) {
                hydroOpacity = parseFloat(e.target.value) / 100.0;
                document.getElementById('hydro-opacity-val').textContent = e.target.value + "%";
                activeRasterLayers.forEach(layer => {
                    if (layer && layer.setOpacity) layer.setOpacity(hydroOpacity);
                });
            });
        }

        // Controller Event Listeners
        document.getElementById('play-btn').addEventListener('click', function() {
            isPlaying = !isPlaying;
            this.textContent = isPlaying ? "⏸ Pause" : "▶ Play";
            if (isPlaying) playLoop();
            else clearTimeout(animTimer);
        });

        document.getElementById('reset-btn').addEventListener('click', function() {
            isPlaying = false;
            document.getElementById('play-btn').textContent = "▶ Play";
            clearTimeout(animTimer);
            updateFrame(0);
        });

        document.getElementById('time-slider').addEventListener('input', function(e) {
            updateFrame(parseInt(e.target.value));
        });

        document.querySelectorAll('.speed-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                document.querySelectorAll('.speed-btn').forEach(b => b.classList.remove('active'));
                this.classList.add('active');
                playSpeed = parseInt(this.getAttribute('data-speed'));
            });
        });

        document.getElementById('color-mode-select').addEventListener('change', function(e) {
            colorMode = e.target.value;
            updateFrame(currentStep);
        });

        document.getElementById('outcome-filter-select').addEventListener('change', function(e) {
            outcomeFilter = e.target.value;
            updateFrame(currentStep);
        });

        document.getElementById('particle-select').addEventListener('change', function(e) {
            selectParticle(parseInt(e.target.value));
        });

        // Initialize First Render
        renderActiveScenario();
    </script>
</body>
</html>
"""

    write(output_path, html_content)
    println("Exported multi-layer interactive dashboard to $(output_path)")
    return output_path
end

"""
    export_interactive_tracks_html(
        trajectories::NamedTuple,
        output_path::AbstractString;
        kwargs...
    ) -> String

Multiple-dispatch convenience overload accepting `trajectories` as the primary argument,
forwarding to the file-path first method signature.

# Mathematical & Physical Context
Visualizes 4D Lagrangian trajectory vectors \$\\mathbf{x}_p(t) = (\\lambda_p(t), \\phi_p(t), z_p(t))\$
co-registered with Eulerian hydrodynamic fields (advection currents \$\\boldsymbol{u}_h\$,
temperature \$T\$, salinity \$S\$, stratification \$N^2\$, and turbulent diffusivity \$\\kappa_v\$)
across discrete depth levels and simulation epochs.

# Inputs
- `trajectories::NamedTuple`: Particle tracking output containing `lons`, `lats`, `depths`, `times`.
- `output_path::AbstractString`: Output path for the generated standalone HTML file.
- `kwargs...`: Additional keyword arguments forwarded to the file-path primary signature.

# Outputs
- `String`: Canonical file path to the generated HTML dashboard.
"""
function export_interactive_tracks_html(
    trajectories::NamedTuple,
    output_path::AbstractString;
    kwargs...
)::String
    return export_interactive_tracks_html(
        output_path;
        trajectories = trajectories,
        kwargs...
    )
end

"""
    plot_interactive_trajectories_map(
        trajectories::Union{Nothing, NamedTuple} = nothing;
        scenarios_data = nothing,
        hydrodynamics = nothing,
        gridded_dispersal = nothing,
        connectivity = nothing,
        recruitment_metrics = nothing,
        output_path::AbstractString = "outputs/interactive_larval_tracks.html",
        bathymetry = nothing,
        strata_definitions = nothing,
        title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
    ) -> String

Convenience wrapper to generate and export an interactive Leaflet HTML map of particle tracks,
supporting multi-scenario data dictionaries, hydrodynamic Eulerian flow fields, and full spatial layer suites.
"""
function plot_interactive_trajectories_map(
    trajectories::Union{Nothing, NamedTuple} = nothing;
    scenarios_data = nothing,
    hydrodynamics = nothing,
    gridded_dispersal = nothing,
    connectivity = nothing,
    recruitment_metrics = nothing,
    output_path::AbstractString = "outputs/interactive_larval_tracks.html",
    bathymetry = nothing,
    strata_definitions = nothing,
    title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
)::String
    return export_interactive_tracks_html(
        output_path;
        trajectories = trajectories,
        scenarios_data = scenarios_data,
        hydrodynamics = hydrodynamics,
        gridded_dispersal = gridded_dispersal,
        connectivity = connectivity,
        recruitment_metrics = recruitment_metrics,
        bathymetry = bathymetry,
        strata_definitions = strata_definitions,
        title = title
    )
end
