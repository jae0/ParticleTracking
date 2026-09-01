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
        output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
    )

Plot larval depth trajectories over time to visualize diurnal oscillations
between nighttime surface grazing and daytime sub-surface predator avoidance.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `sample_indices`: Range or vector of particle indices to plot.
- `title::AbstractString`: Figure title.
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
    end

    hlines!(ax, [0.0], color = :gray50, linestyle = :dash)
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
        i = searchsortedlast(lon_bins, x)
        j = searchsortedlast(lat_bins, y)
        if 1 <= i <= nx && 1 <= j <= ny
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
    extract_hydrodynamic_dataset(
        hydro_input::Any;
        domain_lon::Tuple{<:Real, <:Real} = (-68.0, -57.0),
        domain_lat::Tuple{<:Real, <:Real} = (42.0, 47.0),
        grid_size::Tuple{Int, Int} = (35, 30)
    ) -> NamedTuple

Extract, normalize, and format 2D and 3D hydrodynamic model fields (current velocity
vectors, speed, temperature, salinity, vertical velocity, free-surface elevation,
and bathymetry) for visualization and Leaflet layer serialization.

# Mathematical Formulations
- **Horizontal Advection Current Velocity**:
  ```math
  \\boldsymbol{u}_h(x, y, z) = (u(x, y, z), v(x, y, z))
  ```
  ```math
  |\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}
  ```
  ```math
  \\theta_{\\text{dir}} = \\text{mod}\\left(90^\\circ - \\text{atan2}(v, u) \\cdot \\frac{180^\\circ}{\\pi}, 360^\\circ\\right)
  ```
- **Seawater Temperature & Thermal Stratification**:
  \$T(x, y, z)\$ in °C with vertical thermocline gradient.
- **Practical Salinity**:
  \$S(x, y, z)\$ in PSU.
- **Vertical Velocity (Upwelling/Downwelling)**:
  \$w(x, y, z)\$ in \$\\text{m s}^{-1}\$.
- **Free-Surface Height**:
  \$\\eta(x, y)\$ in meters.
- **Seafloor Bathymetric Elevation**:
  \$H(x, y)\$ in meters.

# Inputs
- `hydro_input`: Oceananigans `HydrostaticFreeSurfaceModel`, `NamedTuple`, `AbstractDict`,
  DuckDB loaded dataset, or `nothing`.
- `domain_lon`: Longitudinal range `(min_lon, max_lon)`.
- `domain_lat`: Latitudinal range `(min_lat, max_lat)`.
- `grid_size`: Target 2D grid resolution `(nx, ny)` for raster serialization.

# Outputs
- `NamedTuple`:
  - `lons::Vector{Float64}`: Longitude coordinates.
  - `lats::Vector{Float64}`: Latitude coordinates.
  - `depths::Vector{Float64}`: Depth level coordinates.
  - `u::Array{Float64, 3}`: Zonal velocity (\$m s^{-1}\$).
  - `v::Array{Float64, 3}`: Meridional velocity (\$m s^{-1}\$).
  - `w::Array{Float64, 3}`: Vertical velocity (\$m s^{-1}\$).
  - `speed::Array{Float64, 3}`: Current speed (\$m s^{-1}\$).
  - `temperature::Array{Float64, 3}`: Seawater temperature (°C).
  - `salinity::Array{Float64, 3}`: Practical salinity (PSU).
  - `elevation::Matrix{Float64}`: Free-surface elevation \$\\eta\$ (m).
  - `bathymetry::Matrix{Float64}`: Seafloor elevation (m).
"""
function extract_hydrodynamic_dataset(
    hydro_input::Any;
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    grid_size::Tuple{Int, Int} = (35, 30)
)
    nx, ny = grid_size
    target_lons = range(domain_lon[1], domain_lon[2], length = nx)
    target_lats = range(domain_lat[1], domain_lat[2], length = ny)
    target_depths = [-2.5, -25.0, -50.0, -100.0]
    nz = length(target_depths)

    # 1. Direct Oceananigans Model Instance
    if !isnothing(hydro_input) && hasproperty(hydro_input, :velocities) && hasproperty(hydro_input, :tracers)
        g = hydro_input.grid
        under_g = g isa ImmersedBoundaryGrid ? g.underlying_grid : g
        m_lons = collect(under_g.λᶠᵃᵃ[1:under_g.Nx])
        m_lats = collect(under_g.φᵃᶠᵃ[1:under_g.Ny])
        m_depths = collect(under_g.zᵃᵃᶜ[1:under_g.Nz])

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

        spd_arr = hypot.(u_arr, v_arr)
        return (
            lons = collect(Float64, m_lons),
            lats = collect(Float64, m_lats),
            depths = collect(Float64, m_depths),
            u = Float64.(u_arr),
            v = Float64.(v_arr),
            w = Float64.(w_arr),
            speed = Float64.(spd_arr),
            temperature = Float64.(t_arr),
            salinity = Float64.(s_arr),
            elevation = Float64.(elev_mat),
            bathymetry = Float64.(bathy_mat)
        )
    end

    # 2. NamedTuple / Dict representation
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

        in_lons = get_field((:lons, :grid_lons, :lon), target_lons)
        in_lats = get_field((:lats, :grid_lats, :lat), target_lats)
        in_depths = get_field((:depths, :grid_depths, :depth), target_depths)

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

        to_3d(arr) = arr isa AbstractMatrix ? reshape(arr, size(arr, 1), size(arr, 2), 1) : arr
        u_3d = to_3d(in_u)
        v_3d = to_3d(in_v)
        w_3d = to_3d(in_w)
        t_3d = to_3d(in_t)
        s_3d = to_3d(in_s)
        spd_3d = hypot.(u_3d, v_3d)

        return (
            lons = collect(Float64, in_lons),
            lats = collect(Float64, in_lats),
            depths = collect(Float64, in_depths),
            u = Float64.(u_3d),
            v = Float64.(v_3d),
            w = Float64.(w_3d),
            speed = Float64.(spd_3d),
            temperature = Float64.(t_3d),
            salinity = Float64.(s_3d),
            elevation = Float64.(in_elev isa AbstractMatrix ? in_elev : fill(0.0, nx_in, ny_in)),
            bathymetry = Float64.(in_bathy isa AbstractMatrix ? in_bathy : fill(-150.0, nx_in, ny_in))
        )
    end

    # 3. Default Realistic Scotian Shelf Regional Hydrodynamic Synthesis
    u_mat = zeros(Float64, nx, ny, nz)
    v_mat = zeros(Float64, nx, ny, nz)
    w_mat = zeros(Float64, nx, ny, nz)
    t_mat = zeros(Float64, nx, ny, nz)
    s_mat = zeros(Float64, nx, ny, nz)
    elev_mat = zeros(Float64, nx, ny)
    bathy_mat = zeros(Float64, nx, ny)

    lons_vec = collect(Float64, target_lons)
    lats_vec = collect(Float64, target_lats)

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

        # Free Surface SSH: cross-shelf steric setup
        elev_mat[i, j] = 0.05 * sin(2.0 * π * x_norm) - 0.04 * cos(π * y_norm)

        for k in 1:nz
            z = target_depths[k]
            depth_factor = exp(z / 75.0)

            # Alongshore Scotian Current Jet
            jet_core = exp(-((dist_to_coast - 0.4)^2) / 0.18)
            u_base = -0.09 * depth_factor * (0.6 + 0.8 * jet_core) + 0.02 * sin(2.0 * π * y_norm)
            v_base = -0.04 * depth_factor * (0.5 + 0.7 * jet_core) + 0.015 * cos(2.0 * π * x_norm)

            u_mat[i, j, k] = u_base
            v_mat[i, j, k] = v_base

            # Coastal upwelling / downwelling vertical velocity
            w_mat[i, j, k] = 0.00035 * sin(2.0 * π * x_norm) * sin(π * y_norm) * (1.0 + z / 100.0)

            # Thermal Stratification: Surface warm layer, Cold Intermediate Layer (CIL) at -50m
            t_surface = 14.2 - 2.5 * y_norm + 1.2 * x_norm
            t_cil = 2.2 + 0.8 * sin(π * x_norm)
            t_slope = 7.5 + 1.5 * (1.0 - y_norm)

            t_val = if z > -20.0
                t_surface + (z / 20.0) * (t_surface - 6.0)
            elseif z > -70.0
                t_cil + ((z + 50.0) / 30.0)^2 * 2.5
            else
                t_cil + ((abs(z) - 70.0) / 50.0) * (t_slope - t_cil)
            end
            t_mat[i, j, k] = clamp(t_val, 0.5, 18.0)

            # Salinity: Fresh coastal water, Mid-shelf, Deep slope water
            s_val = 31.4 + 1.8 * (1.0 - y_norm) + 1.2 * x_norm + (abs(z) / 100.0) * 1.1
            s_mat[i, j, k] = clamp(s_val, 30.2, 35.6)
        end
    end

    spd_mat = hypot.(u_mat, v_mat)

    # Mask all ocean fields to NaN over land so quivers and heatmaps do not
    # render over coastal landmasses. is_point_on_land uses polygon ray-casting.
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

    return (
        lons = lons_vec,
        lats = lats_vec,
        depths = target_depths,
        u = u_mat,
        v = v_mat,
        w = w_mat,
        speed = spd_mat,
        temperature = t_mat,
        salinity = s_mat,
        elevation = elev_mat,
        bathymetry = bathy_mat
    )
end

"""
    plot_hydrodynamic_advection(
        hydrodynamics::Any;
        depth_level::Int = 1,
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        title::AbstractString = "Ocean Hydrodynamic Advection Velocity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
        quiver_stride::Int = 2
    ) -> Figure

Plot 2D horizontal advection current velocity vector arrows \$\\boldsymbol{u}_h = (u, v)\$
over current speed magnitude \$|\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}\$ or background
bathymetry contours.

# Mathematical Formulation
Horizontal flow field \$\\boldsymbol{u}_h(x, y, z_k) = (u, v)\$ with speed:
```math
|\\boldsymbol{u}_h| = \\sqrt{u^2 + v^2}
```
Directional orientation:
```math
\\theta = \\operatorname{atan2}(v, u)
```

# Inputs
- `hydrodynamics`: Oceananigans model instance, NamedTuple, or Dict containing `(lons, lats, u, v)`.
- `depth_level::Int`: Vertical level index to plot (default: 1 for surface).
- `bathymetry_data::Union{Nothing, NamedTuple}`: Optional `(lon, lat, elevation)` background.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save figure (PNG/PDF).
- `quiver_stride::Int`: Spatial stride for subsampling vector arrows.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Marshall, J., et al. (1997). A finite-volume, incompressible Navier Stokes model for
  studies of the ocean on parallel computers. *J. Geophys. Res. Oceans*, 102(C3), 5753-5766.
"""
function plot_hydrodynamic_advection(
    hydrodynamics::Any;
    depth_level::Int = 1,
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    title::AbstractString = "Ocean Hydrodynamic Advection Velocity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
    quiver_stride::Int = 2
)
    hydro = extract_hydrodynamic_dataset(hydrodynamics)
    lons = hydro.lons
    lats = hydro.lats
    nz = length(hydro.depths)
    k = clamp(depth_level, 1, nz)
    depth_m = hydro.depths[k]

    u_k = hydro.u[:, :, k]
    v_k = hydro.v[:, :, k]
    spd_k = hydro.speed[:, :, k] .* 100.0

    fig = Figure(size = (950, 720), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = "\$(title) [Depth: \$(depth_m) m]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    # 1. Background speed heatmap
    hm = heatmap!(ax, lons, lats, spd_k, colormap = :turbo)
    Colorbar(fig[1, 2], hm, label = "Current Speed (|u_h|, cm/s)")

    # 2. Optional bathymetry contour lines
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

    # 3. Vector arrows
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
        depth_level::Int = 1,
        title::AbstractString = "Hydrodynamic Model Seawater Tracers",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png"
    ) -> Figure

Render a two-panel spatial distribution of active seawater tracer fields
(Temperature \$T\$ and Practical Salinity \$S\$) from the ocean hydrodynamic model.

# Mathematical Formulations
- **Thermal Tracer Distribution**:
  ```math
  T(x, y, z_k) \\quad [^\\circ\\text{C}]
  ```
- **Haline Tracer Distribution**:
  ```math
  S(x, y, z_k) \\quad [\\text{PSU}]
  ```

# Inputs
- `hydrodynamics`: Oceananigans model instance, NamedTuple, or Dict containing `(lons, lats, temperature, salinity)`.
- `depth_level::Int`: Vertical level index to plot (default: 1 for surface).
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Destination path for figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_hydrodynamic_tracers(
    hydrodynamics::Any;
    depth_level::Int = 1,
    title::AbstractString = "Hydrodynamic Model Seawater Tracers",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png"
)
    hydro = extract_hydrodynamic_dataset(hydrodynamics)
    lons = hydro.lons
    lats = hydro.lats
    nz = length(hydro.depths)
    k = clamp(depth_level, 1, nz)
    depth_m = hydro.depths[k]

    t_mat = hydro.temperature[:, :, k]
    s_mat = hydro.salinity[:, :, k]

    fig = Figure(size = (1150, 520), fontsize = 12)

    # Panel 1: Temperature
    ax1 = Axis(
        fig[1, 1],
        title = "Seawater Temperature T (°C) [Depth: \$(depth_m) m]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lons, lats, t_mat, colormap = :thermal)
    Colorbar(fig[1, 2], hm1, label = "Temperature (°C)")

    # Panel 2: Salinity
    ax2 = Axis(
        fig[1, 3],
        title = "Practical Salinity S (PSU) [Depth: \$(depth_m) m]",
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
        trajs = scen_data.trajectories
        lons = trajs.lons
        lats = trajs.lats
        depths = trajs.depths
        stages = trajs.stages
        times = trajs.times
        n_p, n_t = size(lons)
        t_days = [round(t / 86400.0, digits = 3) for t in times]

        temps = hasproperty(trajs, :temperatures) && !isnothing(trajs.temperatures) ?
                trajs.temperatures : fill(4.0, n_p, n_t)
        dds = hasproperty(trajs, :degree_days_timeseries) && !isnothing(trajs.degree_days_timeseries) ?
              trajs.degree_days_timeseries :
              (hasproperty(trajs, :degree_days) && trajs.degree_days isa AbstractMatrix ?
               trajs.degree_days :
               repeat(hasproperty(trajs, :degree_days) ? trajs.degree_days : fill(40.0, n_p), 1, n_t))
        survs = hasproperty(trajs, :survival_probability) && !isnothing(trajs.survival_probability) ?
                trajs.survival_probability : fill(0.95, n_p, n_t)

        statuses = [string(s) for s in trajs.settlement_status]
        alive = collect(trajs.alive)
        ages = hasproperty(trajs, :settlement_age) ? collect(trajs.settlement_age) : fill(times[end], n_p)
        ids = collect(trajs.ids)

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

        # Sampled vector current quivers
        h_vec_entries = String[]
        stride_h = max(1, round(Int, nx_h / 20))
        for i in 1:stride_h:nx_h, j in 1:stride_h:ny_h
            lon_v = round(h_lons[i], digits = 4)
            lat_v = round(h_lats[j], digits = 4)
            u_val = round(hydro_data.u[i, j, 1] * 100.0, digits = 2)
            v_val = round(hydro_data.v[i, j, 1] * 100.0, digits = 2)
            spd_val = round(hydro_data.speed[i, j, 1] * 100.0, digits = 2)
            w_val = round(hydro_data.w[i, j, 1] * 1000.0, digits = 3)
            t_val = round(hydro_data.temperature[i, j, 1], digits = 2)
            s_val = round(hydro_data.salinity[i, j, 1], digits = 2)
            eta_val = round(hydro_data.elevation[i, j] * 100.0, digits = 1)
            b_val = round(hydro_data.bathymetry[i, j], digits = 1)
            dir_deg = round(mod(90.0 - rad2deg(atan(hydro_data.v[i, j, 1], hydro_data.u[i, j, 1])), 360.0), digits = 1)

            push!(h_vec_entries, """{"lon":$(lon_v),"lat":$(lat_v),"u":$(u_val),"v":$(v_val),"speed":$(spd_val),"dir":$(dir_deg),"w":$(w_val),"temp":$(t_val),"sal":$(s_val),"eta":$(eta_val),"depth":$(b_val)}""")
        end

        # 2D scalar matrices for smooth bilinear raster interpolation
        temp_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.temperature[i, j, 1], digits = 2)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        sal_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.salinity[i, j, 1], digits = 2)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        spd_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.speed[i, j, 1] * 100.0, digits = 2)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        eta_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.elevation[i, j] * 100.0, digits = 1)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        w_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.w[i, j, 1] * 1000.0, digits = 3)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        bathy_grid_json = "[" * join([
            "[" * join([string(round(hydro_data.bathymetry[i, j], digits = 1)) for j in 1:ny_h], ",") * "]"
            for i in 1:nx_h
        ], ",") * "]"

        t_min, t_max = round(minimum(hydro_data.temperature[:, :, 1]), digits=2), round(maximum(hydro_data.temperature[:, :, 1]), digits=2)
        s_min, s_max = round(minimum(hydro_data.salinity[:, :, 1]), digits=2), round(maximum(hydro_data.salinity[:, :, 1]), digits=2)
        spd_max = round(maximum(hydro_data.speed[:, :, 1] * 100.0), digits=1)
        eta_min, eta_max = round(minimum(hydro_data.elevation * 100.0), digits=1), round(maximum(hydro_data.elevation * 100.0), digits=1)
        w_min, w_max = round(minimum(hydro_data.w[:, :, 1] * 1000.0), digits=2), round(maximum(hydro_data.w[:, :, 1] * 1000.0), digits=2)
        b_min, b_max = round(minimum(hydro_data.bathymetry), digits=1), round(maximum(hydro_data.bathymetry), digits=1)

        hydro_json = """{
            "lons": [$(join([string(round(x, digits=4)) for x in h_lons], ","))],
            "lats": [$(join([string(round(y, digits=4)) for y in h_lats], ","))],
            "depths": [$(join([string(round(z, digits=1)) for z in h_depths], ","))],
            "vectors": [$(join(h_vec_entries, ","))],
            "grids": {
                "temperature": {"min": $(t_min), "max": $(t_max), "unit": "°C", "data": $(temp_grid_json)},
                "salinity": {"min": $(s_min), "max": $(s_max), "unit": "PSU", "data": $(sal_grid_json)},
                "speed": {"min": 0.0, "max": $(spd_max), "unit": "cm/s", "data": $(spd_grid_json)},
                "elevation": {"min": $(eta_min), "max": $(eta_max), "unit": "cm", "data": $(eta_grid_json)},
                "upwelling": {"min": $(w_min), "max": $(w_max), "unit": "mm/s", "data": $(w_grid_json)},
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
        const hydroEtaLayer = L.layerGroup();
        const hydroUpwellingLayer = L.layerGroup();
        const hydroBathyLayer = L.layerGroup();

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

        function renderActiveScenario() {
            [fullTracksLayer, progressTracksLayer, activeHeadsLayer, releaseMarkersLayer,
             settlementMarkersLayer, settlementDensityLayer, hydroAdvectionLayer,
             hydroSpeedLayer, hydroTempLayer, hydroSalLayer, hydroEtaLayer,
             hydroUpwellingLayer, hydroBathyLayer, connectivityFlowsLayer].forEach(l => l.clearLayers());

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

            // 3. Hydrodynamic Model Outputs
            const hydro = currentScen.hydrodynamics;
            if (hydro && hydro.lons && hydro.lats) {
                const hLons = hydro.lons;
                const hLats = hydro.lats;
                const bounds = [[Math.min(...hLats), Math.min(...hLons)], [Math.max(...hLats), Math.max(...hLons)]];

                // A. Advection Current Vectors (Quivers)
                if (hydro.vectors && hydro.vectors.length > 0) {
                    hydro.vectors.forEach(vec => {
                        const lenScale = 0.0035 * Math.min(18.0, Math.max(1.0, vec.speed));
                        const lat2 = vec.lat + (vec.v / Math.max(0.1, vec.speed)) * lenScale;
                        const lon2 = vec.lon + (vec.u / Math.max(0.1, vec.speed)) * lenScale * 1.35;

                        const spdNorm = Math.min(1.0, vec.speed / 20.0);
                        const arrowCol = spdNorm > 0.6 ? '#f43f5e' : spdNorm > 0.3 ? '#38bdf8' : '#34d399';

                        const line = L.polyline([[vec.lat, vec.lon], [lat2, lon2]], {
                            color: arrowCol,
                            weight: 2.2,
                            opacity: 0.85
                        }).bindTooltip(`<b>Hydrodynamic Current Vector</b><br>Speed: <b>\${vec.speed} cm/s</b> (\${(vec.speed*0.0194).toFixed(2)} kts)<br>Dir: \${vec.dir}°<br>u: \${vec.u} cm/s, v: \${vec.v} cm/s<br>Temp: \${vec.temp} °C | Salinity: \${vec.sal} PSU<br>Seabed Depth: \${vec.depth} m`);
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

                // B. Continuous Raster Overlays (Temperature, Salinity, Speed, Elevation, Upwelling, Bathymetry)
                const g = hydro.grids;
                if (g) {
                    if (g.temperature) {
                        const tUrl = generateRasterDataURL(g.temperature.data, hLons, hLats, g.temperature.min, g.temperature.max, colormapThermal);
                        const tOverlay = L.imageOverlay(tUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroTempLayer.addLayer(tOverlay);
                        activeRasterLayers.push(tOverlay);
                    }

                    if (g.salinity) {
                        const sUrl = generateRasterDataURL(g.salinity.data, hLons, hLats, g.salinity.min, g.salinity.max, colormapHaline);
                        const sOverlay = L.imageOverlay(sUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroSalLayer.addLayer(sOverlay);
                        activeRasterLayers.push(sOverlay);
                    }

                    if (g.speed) {
                        const spdUrl = generateRasterDataURL(g.speed.data, hLons, hLats, 0.0, g.speed.max, colormapTurbo);
                        const spdOverlay = L.imageOverlay(spdUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroSpeedLayer.addLayer(spdOverlay);
                        activeRasterLayers.push(spdOverlay);
                    }

                    if (g.elevation) {
                        const etaUrl = generateRasterDataURL(g.elevation.data, hLons, hLats, g.elevation.min, g.elevation.max, colormapCoolwarm);
                        const etaOverlay = L.imageOverlay(etaUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroEtaLayer.addLayer(etaOverlay);
                        activeRasterLayers.push(etaOverlay);
                    }

                    if (g.upwelling) {
                        const wUrl = generateRasterDataURL(g.upwelling.data, hLons, hLats, g.upwelling.min, g.upwelling.max, colormapBalance);
                        const wOverlay = L.imageOverlay(wUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroUpwellingLayer.addLayer(wOverlay);
                        activeRasterLayers.push(wOverlay);
                    }

                    if (g.bathymetry) {
                        const bUrl = generateRasterDataURL(g.bathymetry.data, hLons, hLats, g.bathymetry.min, 0.0, colormapOcean);
                        const bOverlay = L.imageOverlay(bUrl, bounds, { opacity: hydroOpacity, interactive: false });
                        hydroBathyLayer.addLayer(bOverlay);
                        activeRasterLayers.push(bOverlay);
                    }
                }
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
            if (hydro && hydro.lons && hydro.lats && hydro.grids) {
                const hLons = hydro.lons;
                const hLats = hydro.lats;
                const minLon = hLons[0], maxLon = hLons[hLons.length - 1];
                const minLat = hLats[0], maxLat = hLats[hLats.length - 1];

                if (lon >= minLon && lon <= maxLon && lat >= minLat && lat <= maxLat) {
                    const nx = hLons.length, ny = hLats.length;
                    const i = Math.min(nx - 1, Math.max(0, Math.floor(((lon - minLon) / (maxLon - minLon)) * (nx - 1))));
                    const j = Math.min(ny - 1, Math.max(0, Math.floor(((lat - minLat) / (maxLat - minLat)) * (ny - 1))));

                    const g = hydro.grids;
                    document.getElementById('hud-id').textContent = `\${lat.toFixed(2)}°N, \${lon.toFixed(2)}°E`;
                    document.getElementById('hud-stage').textContent = "Eulerian Field";
                    document.getElementById('hud-depth').textContent = `\${g.bathymetry ? g.bathymetry.data[i][j] : -120} m`;
                    document.getElementById('hud-temp').textContent = `\${g.temperature ? g.temperature.data[i][j] : 8.0} °C`;
                    document.getElementById('hud-sal').textContent = `\${g.salinity ? g.salinity.data[i][j] : 32.5} PSU`;
                    document.getElementById('hud-current').textContent = `\${g.speed ? g.speed.data[i][j] : 8.0} cm/s`;
                    document.getElementById('hud-surv').textContent = `SSH \${g.elevation ? g.elevation.data[i][j] : 0.0} cm`;
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
            { id: 'layer-hydro-speed', layer: hydroSpeedLayer, onSelect: () => updateLegend("Current Speed (|u_h|)", "cm/s", "0.0", currentScen.hydrodynamics?.grids?.speed?.max || "25.0", "linear-gradient(to right, #3b82f6, #06b6d4, #10b981, #facc15, #ef4444)") },
            { id: 'layer-hydro-temp', layer: hydroTempLayer, onSelect: () => updateLegend("Seawater Temperature (T)", "°C", currentScen.hydrodynamics?.grids?.temperature?.min || "1.5", currentScen.hydrodynamics?.grids?.temperature?.max || "14.5", "linear-gradient(to right, #1d4ed8, #06b6d4, #10b981, #f59e0b, #ef4444)") },
            { id: 'layer-hydro-sal', layer: hydroSalLayer, onSelect: () => updateLegend("Practical Salinity (S)", "PSU", currentScen.hydrodynamics?.grids?.salinity?.min || "31.0", currentScen.hydrodynamics?.grids?.salinity?.max || "35.5", "linear-gradient(to right, #10b981, #06b6d4, #3b82f6, #6366f1, #8b5cf6)") },
            { id: 'layer-hydro-eta', layer: hydroEtaLayer, onSelect: () => updateLegend("Free Surface Height (η)", "cm", currentScen.hydrodynamics?.grids?.elevation?.min || "-10.0", currentScen.hydrodynamics?.grids?.elevation?.max || "+10.0", "linear-gradient(to right, #3b82f6, #f8fafc, #ef4444)") },
            { id: 'layer-hydro-w', layer: hydroUpwellingLayer, onSelect: () => updateLegend("Vertical Velocity / Upwelling", "mm/s", currentScen.hydrodynamics?.grids?.upwelling?.min || "-0.5", currentScen.hydrodynamics?.grids?.upwelling?.max || "+0.5", "linear-gradient(to right, #1e40af, #94a3b8, #b91c1c)") },
            { id: 'layer-hydro-bathy', layer: hydroBathyLayer, onSelect: () => updateLegend("Seafloor Bathymetry", "m", currentScen.hydrodynamics?.grids?.bathymetry?.min || "-3000", "0", "linear-gradient(to right, #030712, #1e3a8a, #0284c7, #38bdf8)") },
            { id: 'layer-strata', layer: strataLayer },
            { id: 'layer-conn', layer: connectivityFlowsLayer }
        ];

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
