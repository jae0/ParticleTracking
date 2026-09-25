"""
    gl_visualization.jl

Professional 3D hydrodynamic and Lagrangian visualizations using GLMakie (OpenGL backend).
Provides high-quality, interactive 3D figures with GPU-accelerated rendering.
"""

module GLVisualization

using GLMakie
using GLMakie.Makie
using GeometryBasics
using Colors
using ParticleTracking
import ParticleTracking: extract_hydrodynamic_dataset

"""
    plot_3d_hydrodynamic_field(
        hydrodynamics::Any;
        variable::Symbol = :temperature,
        time_seconds::Union{Nothing, Real} = nothing,
        time_index::Union{Nothing, Int} = 1,
        colormap::Symbol = :thermal,
        output_path::Union{Nothing, AbstractString} = "outputs/3d_hydrodynamic_field.png",
        title::Union{Nothing, AbstractString} = nothing,
        resolution::Tuple{Int, Int} = (1920, 1080),
        show_bathymetry::Bool = true,
        show_surface::Bool = false,
        show_particles::Union{Nothing, NamedTuple} = nothing,
        isosurface_levels::Union{Nothing, AbstractVector} = nothing,
        transparency::Float64 = 0.3
    ) -> Figure

Generate a professional 3D visualization of hydrodynamic fields using GLMakie.

Features:
- Volume rendering with customizable transparency
- Isosurface extraction for key thresholds
- Bathymetry as 3D terrain mesh
- Optional particle overlay with depth coloration
- High-resolution output suitable for publications

# Inputs
- `hydrodynamics`: JLD2 file, DuckDB, model instance, or NamedTuple from `extract_hydrodynamic_dataset`
- `variable`: Variable to visualize (`:temperature`, `:salinity`, `:speed`, `:vorticity`, `:stratification`, `:diffusion`, `:density`)
- `time_seconds`: Simulation time in seconds (nearest snapshot)
- `time_index`: Snapshot index (1-based, overrides `time_seconds`)
- `colormap`: Makie colormap symbol (default: `:thermal` for temperature)
- `output_path`: Destination file path
- `title`: Optional custom title
- `resolution`: Output resolution (default: 1920×1080)
- `show_bathymetry`: Render seafloor as 3D mesh
- `show_surface`: Render sea surface elevation
- `show_particles`: Optional Lagrangian trajectories for overlay
- `isosurface_levels`: Values for isosurface extraction (e.g., [0.0, 2.0, 4.0] for temperature)
- `transparency`: Volume rendering alpha (0.0-1.0)

# Outputs
- `Figure`: GLMakie figure with 3D scene
"""
function plot_3d_hydrodynamic_field(
    hydrodynamics::Any;
    variable::Symbol = :temperature,
    time_seconds::Union{Nothing, Real} = nothing,
    time_index::Union{Nothing, Int} = 1,
    colormap::Symbol = :thermal,
    output_path::Union{Nothing, AbstractString} = "outputs/3d_hydrodynamic_field.png",
    title::Union{Nothing, AbstractString} = nothing,
    resolution::Tuple{Int, Int} = (1920, 1080),
    show_bathymetry::Bool = true,
    show_surface::Bool = false,
    show_particles::Union{Nothing, NamedTuple} = nothing,
    isosurface_levels::Union{Nothing, AbstractVector} = nothing,
    transparency::Float64 = 0.3
)::Figure
    hydro = extract_hydrodynamic_dataset(hydrodynamics; time_seconds = time_seconds, time_index = time_index)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    bathy = hydro.bathymetry
    t_hr = round(hydro.time_seconds / 3600.0, digits = 1)

    nx, ny, nz = length(lons), length(lats), length(depths)

    # Resolve variable 3D field
    field_3d, var_label = if variable == :salinity || variable == :S
        (hydro.salinity, "Practical Salinity (PSU)")
    elseif variable == :density || variable == :rho
        (hydro.density, "Potential Density (kg/m³)")
    elseif variable == :stratification || variable == :N2
        (hydro.stratification .* 10000.0, "N² (10⁻⁴ s⁻²)")
    elseif variable == :diffusion || variable == :kappa
        (hydro.diffusion .* 10000.0, "Diffusivity κ_v (10⁻⁴ m²/s)")
    elseif variable == :viscosity || variable == :nu
        (hydro.viscosity .* 10000.0, "Viscosity ν_v (10⁻⁴ m²/s)")
    elseif variable == :speed
        (hydro.speed .* 100.0, "Current Speed (cm/s)")
    elseif variable == :vorticity || variable == :zeta
        (hydro.vorticity .* 100000.0, "Relative Vorticity ζ (10⁻⁵ s⁻¹)")
    elseif variable == :elevation || variable == :eta
        (hydro.elevation .* 100.0, "Sea Surface Elevation (cm)")
    else
        (hydro.temperature, "Temperature T (°C)")
    end

    fig_title = isnothing(title) ? "3D Hydrodynamic Field: $variable [t = $(t_hr) h]" : title

    # Create figure with GLMakie backend
    fig = Figure(size = resolution, fontsize = 12)
    ax = Axis3(fig[1, 1], title = fig_title,
               xlabel = "Longitude (°E)", ylabel = "Latitude (°N)", zlabel = "Depth (m)",
               aspect = :data)

    # Create 3D coordinate arrays for mesh
    # We'll use a regular grid for the volume
    lon_grid = repeat(reshape(lons, nx, 1, 1), 1, ny, nz)
    lat_grid = repeat(reshape(lats, 1, ny, 1), nx, 1, nz)
    depth_grid = repeat(reshape(depths, 1, 1, nz), nx, ny, 1)

    # Mask values below bathymetry
    masked_field = copy(field_3d)
    for i in 1:nx, j in 1:ny
        b = bathy[i, j]
        if !isnan(b)
            for k in 1:nz
                if depths[k] < b
                    masked_field[i, j, k] = NaN
                end
            end
        end
    end

    # Volume rendering
    valid_mask = .!isnan.(masked_field)
    if any(valid_mask)
        vol = volume!(ax, lon_grid, lat_grid, depth_grid, masked_field,
                      algorithm = :mip,  # maximum intensity projection
                      colormap = colormap,
                      colorrange = (minimum(masked_field[valid_mask]), maximum(masked_field[valid_mask])),
                      transparency = transparency)
    end

    # Isosurfaces
    if !isnothing(isosurface_levels)
        for level in isosurface_levels
            if level >= minimum(masked_field[valid_mask]) && level <= maximum(masked_field[valid_mask])
                iso = isosurface!(ax, lon_grid, lat_grid, depth_grid, masked_field, level,
                                 color = (:white, 0.5), show_axis = false)
            end
        end
    end

    # Bathymetry as 3D terrain
    if show_bathymetry
        # Create bathymetry mesh
        bathy_mesh_x = repeat(reshape(lons, nx, 1), 1, ny)
        bathy_mesh_y = repeat(reshape(lats, 1, ny), nx, 1)
        bathy_mesh_z = bathy

        surface!(ax, bathy_mesh_x, bathy_mesh_y, bathy_mesh_z,
                 color = :gray90, transparency = 0.1, shading = NoShading)
    end

    # Sea surface
    if show_surface
        surf_x = repeat(reshape(lons, nx, 1), 1, ny)
        surf_y = repeat(reshape(lats, 1, ny), nx, 1)
        surf_z = hydro.elevation
        surface!(ax, surf_x, surf_y, surf_z,
                 color = :lightblue, transparency = 0.3, shading = NoShading)
    end

    # Particle overlay
    if !isnothing(show_particles)
        n_p, n_t = size(show_particles.lons)
        t_idx = time_index === nothing ? n_t : time_index
        t_idx = clamp(t_idx, 1, n_t)
        
        for p in 1:n_p
            if !isnothing(show_particles.alive) && !show_particles.alive[p]
                continue
            end
            scatter!(ax,
                     show_particles.lons[p, t_idx],
                     show_particles.lats[p, t_idx],
                     show_particles.depths[p, t_idx],
                     color = show_particles.temperatures[p, t_idx],
                     colormap = :thermal,
                     markersize = 15,
                     marker = :sphere)
        end
    end

    Colorbar(fig[1, 2], label = var_label, height = Relative(0.8))

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end


"""
    plot_3d_particle_trajectories(
        trajectories::NamedTuple;
        hydrodynamics::Union{Nothing, NamedTuple} = nothing,
        output_path::Union{Nothing, AbstractString} = "outputs/3d_trajectories.png",
        title::Union{Nothing, AbstractString} = nothing,
        resolution::Tuple{Int, Int} = (1920, 1080),
        show_bathymetry::Bool = true,
        max_display_particles::Int = 100,
        tail_length::Int = 50
    ) -> Figure

Generate professional 3D Lagrangian particle trajectory visualization with GLMakie.

Features:
- 3D trajectory tubes with color gradients
- Bathymetry terrain mesh
- Stage-based color coding
- Animated playback capability

# Inputs
- `trajectories`: Lagrangian tracking output from `track_larval_cohort`
- `hydrodynamics`: Optional hydrodynamic data for bathymetry background
- `output_path`: Destination file path
- `title`: Optional custom title
- `resolution`: Output resolution
- `show_bathymetry`: Render seafloor terrain
- `max_display_particles`: Maximum particles to display
- `tail_length`: Number of time steps to show as trail

# Outputs
- `Figure`: GLMakie figure with 3D trajectories
"""
function plot_3d_particle_trajectories(
    trajectories::NamedTuple;
    hydrodynamics::Union{Nothing, NamedTuple} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/3d_trajectories.png",
    title::Union{Nothing, AbstractString} = nothing,
    resolution::Tuple{Int, Int} = (1920, 1080),
    show_bathymetry::Bool = true,
    max_display_particles::Int = 100,
    tail_length::Int = 50
)::Figure
    n_p, n_t = size(trajectories.lons)
    n_plot = min(n_p, max_display_particles)

    fig_title = isnothing(title) ? "3D Lagrangian Particle Trajectories" : title

    fig = Figure(size = resolution, fontsize = 12)
    ax = Axis3(fig[1, 1], title = fig_title,
               xlabel = "Longitude (°E)", ylabel = "Latitude (°N)", zlabel = "Depth (m)")

    # Bathymetry background
    if show_bathymetry && !isnothing(hydrodynamics)
        hydro = extract_hydrodynamic_dataset(hydrodynamics)
        bathy = hydro.bathymetry
        lons = hydro.lons
        lats = hydro.lats
        
        bathy_mesh_x = repeat(reshape(lons, length(lons), 1), 1, length(lats))
        bathy_mesh_y = repeat(reshape(lats, 1, length(lats)), length(lons), 1)
        surface!(ax, bathy_mesh_x, bathy_mesh_y, bathy,
                 color = :gray80, transparency = 0.2, shading = NoShading)
    end

    # Stage color mapping
    stage_colors = Dict(
        :zoea1 => RGBf(0.4, 0.7, 1.0),
        :zoea2 => RGBf(0.2, 0.5, 1.0),
        :megalopa => RGBf(1.0, 0.6, 0.0),
        :instar1 => RGBf(0.2, 0.8, 0.2),
        :dead => RGBf(1.0, 0.0, 0.0),
        :settled => RGBf(0.6, 0.2, 0.8)
    )

    # Plot trajectories with tails
    for p in 1:n_plot
        if !isnothing(trajectories.alive) && !trajectories.alive[p]
            continue
        end
        
        stage_final = trajectories.stages[p, end]
        col = get(stage_colors, stage_final, RGBf(0.5, 0.5, 0.5))

        # Full trajectory with fading tail
        for t in max(1, n_t - tail_length + 1):n_t
            alpha = 0.2 + 0.8 * (t - max(1, n_t - tail_length + 1)) / tail_length
            if t > 1
                lines!(ax,
                       trajectories.lons[p, t-1:t],
                       trajectories.lats[p, t-1:t],
                       trajectories.depths[p, t-1:t],
                       color = (col, alpha),
                       linewidth = 2 + 3 * (t - max(1, n_t - tail_length + 1)) / tail_length)
            end
        end

        # Final position marker
        scatter!(ax,
                 trajectories.lons[p, end],
                 trajectories.lats[p, end],
                 trajectories.depths[p, end],
                 color = col, markersize = 20, marker = :sphere,
                 strokecolor = :white, strokewidth = 2)
    end

    # Start positions
    start_lons = trajectories.lons[1:n_plot, 1]
    start_lats = trajectories.lats[1:n_plot, 1]
    start_depths = trajectories.depths[1:n_plot, 1]
    scatter!(ax, start_lons, start_lats, start_depths,
             color = :springgreen, markersize = 25, marker = :circle,
             strokecolor = :black, strokewidth = 2)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end


"""
    plot_3d_connectivity(
        connectivity::NamedTuple;
        strata_coords::Union{Nothing, AbstractDict} = nothing,
        output_path::Union{Nothing, AbstractString} = "outputs/3d_connectivity.png",
        title::Union{Nothing, AbstractString} = nothing,
        resolution::Tuple{Int, Int} = (1920, 1080)
    ) -> Figure

Generate 3D connectivity visualization between management strata.

# Inputs
- `connectivity`: Output from `compute_empirical_connectivity`
- `strata_coords`: Dict mapping stratum names to (lon, lat, depth) centroids
- `output_path`: Destination file path
- `title`: Optional custom title
- `resolution`: Output resolution
"""
function plot_3d_connectivity(
    connectivity::NamedTuple;
    strata_coords::Union{Nothing, AbstractDict} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/3d_connectivity.png",
    title::Union{Nothing, AbstractString} = nothing,
    resolution::Tuple{Int, Int} = (1920, 1080)
)::Figure
    mat = connectivity.matrix
    strata = connectivity.strata_names
    n_s = length(strata)

    # Default centroid positions if not provided
    if isnothing(strata_coords)
        # Arrange in circle for visualization
        angles = range(0, 2π, length = n_s + 1)[1:n_s]
        strata_coords = Dict()
        for (i, s) in enumerate(strata)
            strata_coords[s] = (cos(angles[i]) * 5, sin(angles[i]) * 5, -100.0)
        end
    end

    fig_title = isnothing(title) ? "3D Population Connectivity Network" : title

    fig = Figure(size = resolution, fontsize = 12)
    ax = Axis3(fig[1, 1], title = fig_title,
               xlabel = "X", ylabel = "Y", zlabel = "Depth (m)")

    # Stratum nodes
    for (i, s) in enumerate(strata)
        coord = strata_coords[s]
        size_val = 50 + 100 * sum(mat[i, :])  # size by total connectivity
        scatter!(ax, [coord[1]], [coord[2]], [coord[3]],
                 color = :steelblue, markersize = size_val,
                 marker = :sphere, strokecolor = :white, strokewidth = 2)
        text!(ax, coord[1], coord[2], coord[3] + 10, text = s, align = (:center, :center), fontsize = 12)
    end

    # Connectivity arcs
    for i in 1:n_s, j in 1:n_s
        if i != j && mat[i, j] > 0.01
            c1 = strata_coords[strata[i]]
            c2 = strata_coords[strata[j]]
            # Arc path
            n_arc = 20
            arc_x = [c1[1] + (c2[1] - c1[1]) * t for t in range(0, 1, length = n_arc)]
            arc_y = [c1[2] + (c2[2] - c1[2]) * t for t in range(0, 1, length = n_arc)]
            arc_z = [c1[3] + (c2[3] - c1[3]) * t + 50 * sin(π * t) for t in range(0, 1, length = n_arc)]
            
            lines!(ax, arc_x, arc_y, arc_z,
                   color = (:steelblue, mat[i, j]), linewidth = 2 + 5 * mat[i, j])
        end
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end


export plot_3d_hydrodynamic_field, plot_3d_particle_trajectories, plot_3d_connectivity

end # module GLVisualization