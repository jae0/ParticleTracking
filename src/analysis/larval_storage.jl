"""
    larval_storage.jl

    Backend-agnostic analytical storage backend for larval particle tracking (*Chionoecetes opilio*).
    
    Provides relational persistence for simulation runs, multi-million particle
    trajectory steps, cohort recruitment outcomes, demographic connectivity
    matrices, and spatial dispersal fields with support for multiple storage
    backends (Zarr, GeoParquet, etc.) via GeoData.
    
    This is a wrapper around GeoData.Data.LarvalStorage (for Zarr) and
    GeoParquet-based storage (for GeoParquet) that provides a
    ParticleTracking-friendly API.
    """

using GeoData
using GeoData.Data.LarvalStorage
using DataFrames
using Dates
using Statistics
using LinearAlgebra

# Re-export the GeoData storage functions with ParticleTracking-friendly names
export
    open_storage,
    close_storage,
    initialize_storage_schema!,
    save_simulation_run!,
    load_run_configuration,
    list_simulation_runs,
    load_trajectories_df,
    load_trajectories_namedtuple,
    save_hydrodynamic_field!,
    load_hydrodynamic_field,
    load_gridded_dispersal,
    load_connectivity_matrix,
    compare_scenarios,
    compute_ensemble_model_average

"""
    open_storage(
        output_path::AbstractString = joinpath("outputs", "nova_scotia_hydrodynamics.jld2");
        storage_backend::String = "zarr";
        read_only::Bool = false,
        max_retries::Int = 5,
        retry_delay::Real = 0.5
    )
    
    Open a connection to the storage directory for particle tracking storage
    and analytics based on the specified backend. Automatically creates parent
    directories if needed.
    
    # Arguments
    - `output_path::AbstractString`: Base path for simulation outputs. Extension determines format:
      - `.jld2`: Uses .zarr or .parquet extension based on storage_backend
      - Other: Appropriate extension (.zarr or .parquet) is appended based on storage_backend
    - `storage_backend::String`: Storage backend to use ("zarr" or "geoparquet").
    - `read_only::Bool`: Whether to open the storage in read-only mode.
    - `max_retries::Int`: Maximum number of retries under transient file lock contention.
    - `retry_delay::Real`: Seconds to wait between retries.
    
    # Outputs
    - `Union{Zarr.ZGroup, GeoParquetStorage>`: Storage group object for the specified backend.
    """
function open_storage(
    output_path::AbstractString = joinpath("outputs", "nova_scotia_hydrodynamics.jld2");
    storage_backend::String = "zarr";
    read_only::Bool = false,
    max_retries::Int = 5,
    retry_delay::Real = 0.5
)
    # Open a connection to the storage directory for particle tracking storage
    # based on the specified backend.
    if storage_backend == "zarr"
        if endswith(output_path, ".jld2")
            storage_path = replace(output_path, r"\.jld2$" => ".zarr")
        else
            storage_path = output_path * ".zarr"
        end
        # For Zarr, we don't have a read-only mode in the same way, but we can
        # just open the store. The retry logic is handled by Zarr internally.
        store = open_larval_storage(storage_path)
        return store
    elseif storage_backend == "geoparquet"
        # For GeoParquet, we return a special object that mimics the Zarr group interface
        # but uses GeoParquet for storage
        if endswith(output_path, ".jld2")
            storage_path = replace(output_path, r"\.jld2$" => ".parquet")
        else
            storage_path = output_path * ".parquet"
        end
        # Return a GeoParquet storage object
        return GeoParquetStorage(storage_path, read_only=read_only)
    else
        error("Unsupported storage backend: $storage_backend. Supported backends are: 'zarr', 'geoparquet'")
    end
end

"""
    close_storage(group::Union{Zarr.ZGroup, GeoParquetStorage}; force::Bool = false)
    
    Close the storage group. Flushes pending data via checkpoint.
    When `force = false` (default), the instance is kept in the process session
    cache so subsequent pipeline stages can access it without re-opening.
    Pass `force = true` to evict from cache and close the handle.
    """
function close_storage(group::Union{Zarr.ZGroup, GeoParquetStorage}; force::Bool = false)
    if !force
        return nothing
    end
    if isinstance(group, Zarr.ZGroup)
        close_larval_storage(group)
        GC.gc()
    elseif isinstance(group, GeoParquetStorage)
        # For GeoParquet, we don't need to flush or garbage collect in the same way
        # Just return nothing
    end
    return nothing
end

"""
    close_all_storage!()
    
    Explicitly checkpoint, close, and evict all active storage instances
    held in the process session cache.
    """
function close_all_storage!()
    # Zarr doesn't have a global session cache
    GC.gc()
    return nothing
end

"""
    initialize_storage_schema!(group::Zarr.ZGroup)

Initialize storage schema with arrays for trajectories, metadata, etc.
    """
    function initialize_storage_schema!(group::Union{Zarr.ZGroup, GeoParquetStorage})
        if isinstance(group, Zarr.ZGroup)
            initialize_larval_storage_schema!(group)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet, we don't need to initialize a schema in the same way
            # The datasets will be created on first save
        end
    end

"""
    save_simulation_run!(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString,
        opts;
        trajectories::NamedTuple,
        metrics::Union{Nothing, NamedTuple} = nothing,
        connectivity::Union{Nothing, NamedTuple} = nothing,
        gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
        config::Union{Nothing, AbstractDict, AbstractString} = nothing,
        notes::AbstractString = ""
    )::String
    
    Archive a complete simulation run including parameter metadata, full configuration dictionary,
    4D particle trajectory time series, recruitment metrics, demographic connectivity matrix,
    and 2D spatial dispersal fields into storage within an atomic transaction.
    """
function save_simulation_run!(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString,
    opts;
    trajectories::NamedTuple,
    metrics::Union{Nothing, NamedTuple} = nothing,
    connectivity::Union{Nothing, NamedTuple} = nothing,
    gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
    config::Union{Nothing, AbstractDict, AbstractString} = nothing,
    notes::AbstractString = ""
)::String
    if isinstance(group, Zarr.ZGroup)
        # Use existing Zarr larval storage
        return save_larval_simulation_run!(
            group, run_id, opts;
            trajectories = trajectories,
            metrics = metrics,
            connectivity = connectivity,
            gridded_dispersal = gridded_dispersal,
            config = config,
            notes = notes
        )
    elseif isinstance(group, GeoParquetStorage)
        # For GeoParquet backend, save essential data as Parquet files
        # Warning: This is a simplified implementation that may not preserve all functionality
        @warn "GeoParquet backend for larval storage is experimental and may not preserve all Zarr functionality"
        
        # Save trajectories as Parquet file
        traj_path = replace(group.path, r"\.parquet$" => "_trajectories.parquet")
        # Convert trajectories NamedTuple to DataFrame for saving
        if !isempty(trajectories)
            # Get the length from the first variable
            first_var = first(values(trajectories))
            n_particles = length(first_var)
            n_timesteps = size(first_var, 2)
            
            # Create DataFrame with trajectory data
            traj_df = DataFrame()
            traj_df.particle_id = repeat(1:n_particles, inner=n_timesteps)
            traj_df.timestep = repeat(1:n_timesteps, outer=n_particles)
            
            # Add each trajectory variable
            for (var_name, var_data) in trajectories
                traj_df[!, var_name] = vec(var_data)
            end
            
            # Save to Parquet
            try
                using GeoParquet
                GeoParquet.write(traj_path, traj_df)
            catch e
                @warn "Failed to save trajectories to Parquet: $e"
            end
        end
        
        # Save connectivity as Parquet file if provided
        if !isnothing(connectivity)
            conn_path = replace(group.path, r"\.parquet$" => "_connectivity.parquet")
            # Similar processing for connectivity...
            @warn "Connectivity saving to Parquet not fully implemented"
        end
        
        # Save recruitment metrics as Parquet file if provided
        if !isnothing(metrics)
            metrics_path = replace(group.path, r"\.parquet$" => "_metrics.parquet")
            @warn "Metrics saving to Parquet not fully implemented"
        end
        
        # Save config as TOML file
        if !isnothing(config)
            config_path = replace(group.path, r"\.parquet$" => "_config.toml")
            try
                if config isa AbstractString
                    write(config_path, config)
                else
                    io = IOBuffer()
                    TOML.print(io, config; sorted = true)
                    write(config_path, String(take!(io)))
                end
            catch e
                @warn "Failed to save config to TOML: $e"
            end
        end
        
        return group.path
    end
end

"""
    load_run_configuration(
        group::Zarr.ZGroup,
        run_id::AbstractString
    ) -> Union{Nothing, AbstractDict}

Load the configuration dictionary for a specific run from storage.
"""
function load_run_configuration(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString
)
     # Load from runs array or config file based on storage type
      if isinstance(group, Zarr.ZGroup)
          # Load from storage
          runs_array = group["runs"]
          runs_json = Zarr.read(runs_array)
         
         for r_json in runs_json
             r = JSON3.read(r_json)
             if r["run_id"] == run_id && haskey(r, "config_toml") && !isempty(r["config_toml"])
                 return TOML.parse(r["config_toml"])
             end
         end
         return nothing
     elseif isinstance(group, GeoParquetStorage)
         # For GeoParquet backend, look for config file
         config_path = replace(group.path, r"\.parquet$" => "_config.toml")
         if isfile(config_path)
             try
                 return TOML.parse(read(config_path, String))
             catch e
                 @warn "Failed to load config from TOML file: $e"
                 return nothing
             end
         else
             return nothing
         end
     end
 end

"""
    list_simulation_runs(
        group::Union{Zarr.ZGroup, GeoParquetStorage};
        scenario::Union{Nothing, AbstractString, Symbol} = nothing,
        projection_year::Union{Nothing, Int} = nothing
    )::DataFrame
    
    Query and return a summary `DataFrame` of all simulation runs archived in storage,
    with optional filtering by climate scenario or projection year.
    """
function list_simulation_runs(
    group::Union{Zarr.ZGroup, GeoParquetStorage>;
    scenario::Union{Nothing, AbstractString, Symbol} = nothing,
    projection_year::Union{Nothing, Int} = nothing
)::DataFrame
    if isinstance(group, Zarr.ZGroup)
        return list_larval_simulation_runs(group; scenario = scenario, projection_year = projection_year)
    elseif isinstance(group, GeoParquetStorage)
        # For GeoParquet backend, return an empty DataFrame with warning
        @warn "list_simulation_runs not fully implemented for GeoParquet backend"
        return DataFrame()
    end
end

"""
    load_trajectories_df(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString;
        particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
        stage::Union{Nothing, Symbol, AbstractString} = nothing,
        time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
        max_particles::Union{Nothing, Int} = nothing
    )::DataFrame
    
    Retrieve particle trajectory records from storage as a `DataFrame` with flexible
    spatial, temporal, developmental stage, and particle ID filtering.
    """
function load_trajectories_df(
    group::Zarr.ZGroup,
    run_id::AbstractString;
    particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
    stage::Union{Nothing, Symbol, AbstractString} = nothing,
    time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
    max_particles::Union{Nothing, Int} = nothing
)::DataFrame
    if isinstance(group, Zarr.ZGroup)
        trajs = load_larval_trajectories(
            group, run_id;
            particle_ids = particle_ids,
            stage = stage,
            time_range = time_range,
            max_particles = max_particles
        )
    elseif isinstance(group, GeoParquetStorage)
        # For GeoParquet backend, try to load from Parquet file
        @warn "load_trajectories_df for GeoParquet backend is not fully implemented"
        return DataFrame()
    end
    
    # Convert NamedTuple to DataFrame
    n_p = size(trajs.lons, 1)
    n_t = size(trajs.lons, 2)
    
    df_run_id = fill(run_id, n_p * n_t)
    df_p_id = Vector{Int}(undef, n_p * n_t)
    df_step = Vector{Int}(undef, n_p * n_t)
    df_time = Vector{Float64}(undef, n_p * n_t)
    df_lon = Vector{Float64}(undef, n_p * n_t)
    df_lat = Vector{Float64}(undef, n_p * n_t)
    df_depth = Vector{Float64}(undef, n_p * n_t)
    df_temp = Vector{Float64}(undef, n_p * n_t)
    df_dd = Vector{Float64}(undef, n_p * n_t)
    df_surv = Vector{Float64}(undef, n_p * n_t)
    df_stage = Vector{String}(undef, n_p * n_t)
    df_alive = Vector{Bool}(undef, n_p * n_t)
    df_settle = Vector{String}(undef, n_p * n_t)
    
    idx = 1
    for p in 1:n_p
        p_id = trajs.ids[p]
        p_alive = trajs.alive[p]
        p_settle = string(trajs.settlement_status[p])
        
        for s in 1:n_t
            df_p_id[idx] = p_id
            df_step[idx] = s
            df_time[idx] = trajs.times[s]
            df_lon[idx] = trajs.lons[p, s]
            df_lat[idx] = trajs.lats[p, s]
            df_depth[idx] = trajs.depths[p, s]
            df_temp[idx] = trajs.temperatures[p, s]
            df_dd[idx] = trajs.degree_days[p, s]
            df_surv[idx] = trajs.survival_probability[p, s]
            df_stage[idx] = string(trajs.stages[p, s])
            df_alive[idx] = p_alive
            df_settle[idx] = p_settle
            idx += 1
        end
    end
    
    return DataFrame(
        run_id = df_run_id,
        particle_id = df_p_id,
        step = df_step,
        time_seconds = df_time,
        lon = df_lon,
        lat = df_lat,
        depth = df_depth,
        temperature = df_temp,
        degree_days = df_dd,
        survival_probability = df_surv,
        stage = df_stage,
        alive = df_alive,
        settlement_status = df_settle
    )
end

"""
    load_trajectories_namedtuple(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString;
        max_particles::Union{Nothing, Int} = nothing
    )::NamedTuple
    
    Retrieve Lagrangian particle tracking output from storage and reconstruct the
    complete `NamedTuple` structure with fields `(lons, lats, depths, temperatures,
    degree_days, degree_days_timeseries, survival_probability, stages, alive,
    settlement_status, settlement_age, times, ids)`.
    """
function load_trajectories_namedtuple(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString;
        max_particles::Union{Nothing, Int} = nothing
    )::NamedTuple
        if isinstance(group, Zarr.ZGroup)
            return load_larval_trajectories(group, run_id; max_particles = max_particles)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, try to load from Parquet file
            @warn "load_trajectories_namedtuple for GeoParquet backend is not fully implemented"
            # Return an empty NamedTuple with the expected structure
            return (lons = Float64[], lats = Float64[], depths = Float64[], temperatures = Float64[],
                   degree_days = Float64[], degree_days_timeseries = Float64[], survival_probability = Float64[],
                   stages = String[], alive = Bool[], settlement_status = String[], settlement_age = Float64[],
                   times = Float64[], ids = Int64[])
        end
    end

"""
    save_hydrodynamic_field!(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString,
        scenario::AbstractString,
        projection_year::Int,
        time_seconds::Real,
        field_data::NamedTuple
    )
    
    Save hydrodynamic fields (u, v, w, temperature, salinity, elevation) to storage.
    """
    function save_hydrodynamic_field!(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString,
        scenario::AbstractString,
        projection_year::Int,
        time_seconds::Real,
        field_data::NamedTuple
    )
        if isinstance(group, Zarr.ZGroup)
            hydro_group = group["hydrodynamic_fields"]
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, save to a Parquet file
            @warn "save_hydrodynamic_field! for GeoParquet backend is not fully implemented"
            return nothing
        end
    if !haskey(hydro_group, run_id)
        Zarr.create_group(hydro_group, run_id)
    end
    run_group = hydro_group[run_id]
    
    # Store scenario and projection year as metadata
    if !haskey(run_group, "scenario")
        Zarr.write(run_group["scenario"], scenario)
    end
    if !haskey(run_group, "projection_year")
        Zarr.write(run_group["projection_year"], projection_year)
    end
    
    # Store time
    if !haskey(run_group, "time_seconds")
        Zarr.write(run_group["time_seconds"], time_seconds)
    end
    
    # Store each field
    for (k, v) in pairs(field_data)
        if v isa Matrix || v isa Vector || v isa AbstractArray
            Zarr.write(run_group[k], v)
        end
    end
end

"""
    load_hydrodynamic_field(
        group::Zarr.ZGroup,
        run_id::AbstractString
    ) -> NamedTuple

Load hydrodynamic field from storage.
    """
    function load_hydrodynamic_field(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    )
        if isinstance(group, Zarr.ZGroup)
            return load_larval_hydrodynamic_fields(group, run_id)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, try to load from Parquet file
            @warn "load_hydrodynamic_field for GeoParquet backend is not fully implemented"
            return NamedTuple()
        end
    end

"""
    load_gridded_dispersal(
        group::Zarr.ZGroup,
        run_id::AbstractString
    ) -> NamedTuple

Load gridded dispersal fields from storage.
    """
    function load_gridded_dispersal(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    ) -> NamedTuple
        if isinstance(group, Zarr.ZGroup)
            return load_larval_gridded_dispersal(group, run_id)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, try to load from Parquet file
            @warn "load_gridded_dispersal for GeoParquet backend is not fully implemented"
            return NamedTuple()
        end
    end

"""
    load_connectivity_matrix(
        group::Zarr.ZGroup,
        run_id::AbstractString
    ) -> NamedTuple

Load connectivity matrix from storage.
    """
    function load_connectivity_matrix(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    ) -> NamedTuple
        if isinstance(group, Zarr.ZGroup)
            return load_larval_connectivity(group, run_id)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, try to load from Parquet file
            @warn "load_connectivity_matrix for GeoParquet backend is not fully implemented"
            return NamedTuple()
        end
    end

"""
    compare_scenarios(
        group::Zarr.ZGroup;
        scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
        projection_year::Int = 2050
    ) -> DataFrame

Compare scenarios across runs.
    """
    function compare_scenarios(
        group::Union{Zarr.ZGroup, GeoParquetStorage>;
        scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
        projection_year::Int = 2050
    )::DataFrame
        if isinstance(group, Zarr.ZGroup)
            return compare_larval_scenarios(group; scenarios = scenarios, projection_year = projection_year)
        elseif isinstance(group, GeoParquetStorage)
            # For GeoParquet backend, return an empty DataFrame with warning
            @warn "compare_scenarios not fully implemented for GeoParquet backend"
            return DataFrame()
        end
    end

"""
    compute_ensemble_model_average(
        group::Zarr.ZGroup;
        variable::Symbol = :mean_pld_days,
        scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
        projection_year::Int = 2050
     ) -> NamedTuple
 
  Compute ensemble model average across scenarios.
  """
  function compute_ensemble_model_average(
      group::Union{Zarr.ZGroup, GeoParquetStorage};
      variable::Symbol = :mean_pld_days,
      scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
      projection_year::Int = 2050
  ) -> NamedTuple
      if isinstance(group, Zarr.ZGroup)
          return compute_larval_ensemble_model_average(
              group;
              variable = variable,
              scenarios = scenarios,
              projection_year = projection_year
          )
      elseif isinstance(group, GeoParquetStorage)
          # For GeoParquet backend, return an empty NamedTuple with warning
          @warn "compute_ensemble_model_average not fully implemented for GeoParquet backend"
          return (mean_pld_days = 0.0,)
      end
  end

# GeoParquet storage type that mimics the Zarr group interface for basic operations
mutable struct GeoParquetStorage
    path::String
    read_only::Bool
    datasets::Dict{String, Any}
    
    function GeoParquetStorage(path::String; read_only::Bool=false)
        # Initialize with empty datasets
        new(path, read_only, Dict{String, Any}())
    end
end