# ============================================================================
# Additional Visualization Functions
# ============================================================================

"""
    plot_voronoi_tessellation(
        voronoi_data::NamedTuple;
        bathymetry::Union{Nothing, NamedTuple} = nothing,
        trajectories::Union{Nothing, NamedTuple} = nothing,
        strata::Union{Nothing, AbstractVector} = nothing,
        output_path::Union{Nothing, AbstractString} = "outputs/voronoi_tessellation.png",
        title::Union{Nothing, AbstractString} = nothing
    ) -> Figure

Visualize depth-stratified Voronoi tessellation with settlement density and connectivity.

# Inputs
- `voronoi_data`: NamedTuple from Voronoi tessellation (units, centroids, depths, strata, areas, settlement)
- `bathymetry`: Optional bathymetry NamedTuple for background
- `trajectories`: Optional Lagrangian trajectories for overlay
- `strata`: Optional management strata definitions (CFA polygons)
- `output_path`: Destination file path
- `title`: Optional custom title

# Outputs
- `Figure`: CairoMakie figure with tessellation polygons, depth coloring, and settlement metrics
"""
function plot_voronoi_tessellation(
    voronoi_data::NamedTuple;
    bathymetry::Union{Nothing, NamedTuple} = nothing,
    trajectories::Union{Nothing, NamedTuple} = nothing,
    strata::Union{Nothing, AbstractVector} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/voronoi_tessellation.png",
    title::Union{Nothing, AbstractString} = nothing
)::Figure
    units = voronoi_data.units
    n_units = length(units.centroids_lon)
    
    fig = Figure(size = (1400, 1000), fontsize = 11)
    Label(fig[0, 1:2], isnothing(title) ? "Depth-Stratified Voronoi Tessellation" : title, fontsize = 16, font = :bold)
    
    # Panel 1: Voronoi cells colored by stratum
    ax1 = Axis(fig[1, 1], title = "Voronoi Units by Depth Stratum",
               xlabel = "Longitude (°E)", ylabel = "Latitude (°N)", aspect = DataAspect())
    
    stratum_colors = Dict("core" => :steelblue, "shallow" => :gold, "deep" => :purple)
    for i in 1:n_units
        if isnothing(units.polygons[i]) || isempty(units.polygons[i])
            continue
        end
        poly = units.polygons[i]
        stratum = units.stratum[i]
        col = get(stratum_colors, stratum, :gray)
        poly!(ax1, poly, color = (col, 0.5), strokecolor = :white, strokewidth = 0.5)
    end
    
    if !isnothing(bathymetry) && hasproperty(bathymetry, :bathymetry)
        lons = bathymetry.lons
        lats = bathymetry.lats
        bathy = bathymetry.bathymetry
        contour!(ax1, lons, lats, bathy, levels = [-50, -100, -200, -500, -1000],
                 color = :gray, linewidth = 0.5, linestyle = :dash)
    end
    
    # Panel 2: Settlement density by Voronoi unit
    ax2 = Axis(fig[1, 2], title = "Settlement Density per Voronoi Unit",
               xlabel = "Longitude (°E)", ylabel = "Latitude (°N)", aspect = DataAspect())
    settlement = units.settlement_count
    if maximum(settlement) > 0
        norm_settle = settlement ./ maximum(settlement)
        for i in 1:n_units
            if isnothing(units.polygons[i]) || isempty(units.polygons[i])
                continue
            end
            poly = units.polygons[i]
            intensity = norm_settle[i]
            col = RGBf(intensity, 0.0, 1.0 - intensity)
            poly!(ax2, poly, color = (col, 0.7), strokecolor = :white, strokewidth = 0.3)
            if settlement[i] > 0
                text!(ax2, units.centroids_lon[i], units.centroids_lat[i],
                      text = "$(settlement[i])", fontsize = 8, color = :black)
            end
        end
    else
        scatter!(ax2, units.centroids_lon, units.centroids_lat,
                 color = units.depth, colormap = :deep, markersize = 8)
    end
    
    # Panel 3: Depth distribution of Voronoi centroids
    ax3 = Axis(fig[2, 1], title = "Voronoi Centroid Depth Distribution",
               xlabel = "Depth (m)", ylabel = "Count")
    hist!(ax3, -units.depth, bins = 20, color = :steelblue, strokewidth = 0.5)
    
    # Panel 4: Stratum area summary
    ax4 = Axis(fig[2, 2], title = "Area by Depth Stratum",
               xlabel = "Stratum", ylabel = "Total Area (km²)")
    stratum_names = ["shallow", "core", "deep"]
    stratum_areas = Float64[]
    for s in stratum_names
        mask = units.stratum .== s
        push!(stratum_areas, sum(units.area_km2[mask]) / 1e6)
    end
    barplot!(ax4, stratum_names, stratum_areas, color = [:gold, :steelblue, :purple])
    
    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end
    
    return fig
end

"""
    plot_particle_fate_summary(
        trajectories::NamedTuple;
        output_path::Union{Nothing, AbstractString} = "outputs/particle_fate_summary.png",
        title::Union{Nothing, AbstractString} = nothing
    ) -> Figure

Generate a comprehensive particle fate summary dashboard.

Panels:
1. Stage transition timeline (horizontal bar chart per particle)
2. Final status distribution (pie chart)
3. Settlement depth histogram
4. Degree-day accumulation vs. settlement success
5. Thermal exposure (degree-days) distribution
6. Mortality causes breakdown

# Inputs
- `trajectories`: Lagrangian tracking output from `track_larval_cohort`
- `output_path`: Destination file path
- `title`: Optional custom title

# Outputs
- `Figure`: CairoMakie figure with 6 panels summarizing particle fates
"""
function plot_particle_fate_summary(
    trajectories::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/particle_fate_summary.png",
    title::Union{Nothing, AbstractString} = nothing
)::Figure
    n_p, n_t = size(trajectories.lons)
    stages = trajectories.stages
    statuses = trajectories.settlement_status
    degree_days = trajectories.degree_days
    degree_days_ts = trajectories.degree_days_timeseries
    temperatures = trajectories.temperatures
    alive = trajectories.alive
    settlement_age = trajectories.settlement_age
    ids = trajectories.ids
    
    fig = Figure(size = (1600, 1200), fontsize = 11)
    Label(fig[0, 1:2], isnothing(title) ? "Particle Fate Summary Dashboard" : title, fontsize = 16, font = :bold)
    
    # Panel 1: Stage transition timeline
    ax1 = Axis(fig[1, 1], title = "Stage Progression Timeline",
               xlabel = "Time (days)", ylabel = "Particle ID")
    stage_order = [:zoea1, :zoea2, :megalopa, :instar1]
    stage_y = Dict(s => i for (i, s) in enumerate(stage_order))
    for p in 1:n_p
        current_stage = stages[p, 1]
        stage_start = 1
        for t in 2:n_t
            if stages[p, t] != current_stage
                s_y = stage_y[current_stage]
                days_start = trajectories.times[stage_start] / 86400.0
                days_end = trajectories.times[t-1] / 86400.0
                barplot!(ax1, [days_start, days_end], [s_y], [s_y],
                         color = :steelblue, direction = :x)
                current_stage = stages[p, t]
                stage_start = t
            end
        end
        s_y = stage_y[current_stage]
        days_start = trajectories.times[stage_start] / 86400.0
        days_end = trajectories.times[end] / 86400.0
        barplot!(ax1, [days_start, days_end], [s_y], [s_y],
                 color = :steelblue, direction = :x)
    end
    
    ax1.yticks = (collect(values(stage_y)), collect(keys(stage_order)))
    
    # Panel 2: Final status distribution (pie)
    ax2 = Axis(fig[1, 2], title = "Final Status Distribution", aspect = DataAspect())
    status_counts = Dict()
    for s in statuses
        status_counts[s] = get(status_counts, s, 0) + 1
    end
    labels = collect(keys(status_counts))
    values = collect(values(status_counts))
    pie!(ax2, values, label = labels, colors = [:red, :orange, :green, :purple, :gray, :blue])
    
    # Panel 3: Settlement depth histogram
    ax3 = Axis(fig[2, 1], title = "Settlement Depth Distribution",
               xlabel = "Depth (m)", ylabel = "Count")
    settled_depths = Float64[]
    for p in 1:n_p
        if statuses[p] == :settled && trajectories.depths[p, end] < 0
            push!(settled_depths, trajectories.depths[p, end])
        end
    end
    if !isempty(settled_depths)
        hist!(ax3, -settled_depths, bins = 20, color = :steelblue, strokewidth = 0.5)
    end
    
    # Panel 4: Degree-days vs Settlement Success
    ax4 = Axis(fig[2, 2], title = "Cumulative Degree-Days vs Settlement",
               xlabel = "Degree-Days (°C·days)", ylabel = "Settlement Success (0/1)")
    dd_final = [degree_days[p, end] for p in 1:n_p]
    settled_binary = [statuses[p] == :settled ? 1.0 : 0.0 for p in 1:n_p]
    scatter!(ax4, dd_final, settled_binary, color = :steelblue, markersize = 8, alpha = 0.6)
    
    # Panel 5: Thermal exposure distribution
    ax5 = Axis(fig[3, 1], title = "Thermal Exposure (Degree-Days) Distribution",
               xlabel = "Degree-Days (°C·days)", ylabel = "Count")
    hist!(ax5, dd_final, bins = 20, color = :thermal, strokewidth = 0.5)
    
    # Panel 6: Mortality causes
    ax6 = Axis(fig[3, 2], title = "Mortality & Fate Breakdown",
               xlabel = "Category", ylabel = "Count")
    thermal_mort = sum([any(t > 15.0 for t in temperatures[p, :]) for p in 1:n_p])
    cold_mort = sum([any(t < 0.0 for t in temperatures[p, :]) for p in 1:n_p])
    not_settled = sum(statuses .!= :settled)
    settled_count = sum(statuses .== :settled)
    categories = ["Settled", "Not Settled", "Thermal Stress", "Cold Stress"]
    counts = [settled_count, not_settled, thermal_mort, cold_mort]
    barplot!(ax6, categories, counts, color = [:green, :red, :orange, :blue])
    
    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end
    
    return fig
end