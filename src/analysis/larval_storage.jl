"""
    larval_storage.jl

Zarr and GeoParquet analytical storage backends for larval particle tracking
(*Chionoecetes opilio*).

Provides relational persistence for simulation runs, multi-million particle
trajectory steps, cohort recruitment outcomes, demographic connectivity
matrices, and spatial dispersal fields with:
- Lagrangian particle trajectories stored in GeoParquet (.parquet) with Point geometries
- Eulerian hydrodynamic fields stored in chunked Zarr (.zarr) arrays
- Multi-scenario benchmarking and ensemble model averaging
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Zarr
using JSON3
using Dates
using Statistics
using LinearAlgebra
using DataFrames
using TOML
using GeoData: GeoInterface, GeoParquet

"""
    GeoParquetStorage

Container representation for tabular GeoParquet simulation outputs and metrics.
"""
mutable struct GeoParquetStorage
    path::String
    read_only::Bool
    datasets::Dict{String, Any}

    function GeoParquetStorage(path::String; read_only::Bool = false)
        new(path, read_only, Dict{String, Any}())
    end
end

"""
    _storage_base_dir(group::GeoParquetStorage) -> String

Extract the base directory containing the parquet files.
"""
function _storage_base_dir(group::GeoParquetStorage)
    p = group.path
    if isdir(p)
        return p
    elseif endswith(p, ".parquet")
        d = dirname(p)
        return isempty(d) ? "." : d
    else
        return p
    end
end

export
    GeoParquetStorage,
    open_storage,
    close_storage,
    close_all_storage!,
    initialize_storage_schema!,
    save_simulation_run!,
    load_run_configuration,
    list_simulation_runs,
    load_trajectories_df,
    load_trajectories_namedtuple,
    load_all_scenario_trajectories,
    save_hydrodynamic_field!,
    load_hydrodynamic_field,
    load_gridded_dispersal,
    load_connectivity_matrix,
    compare_scenarios,
    compute_ensemble_model_average,
    geopublish_lakehouse!

# ============================================================================
# Core Storage Lifecycle
# ============================================================================

"""
    open_storage(
        output_path::AbstractString = joinpath("outputs", "larval_trajectories.parquet");
        storage_backend::Union{Nothing, String} = nothing,
        read_only::Bool = false,
        max_retries::Int = 5,
        retry_delay::Real = 0.5
    ) -> Union{Zarr.ZGroup, GeoParquetStorage}

Open a connection to the storage location for particle tracking storage
and analytics based on the specified backend or file extension.

# Arguments
- `output_path`: Target filesystem path or directory.
- `storage_backend`: Backend driver (`"geoparquet"` or `"zarr"`).
- `read_only`: Whether to open in read-only mode.
- `max_retries`: Maximum attempts for lock acquisition.
- `retry_delay`: Delay between retry attempts.

# Returns
`Union{Zarr.ZGroup, GeoParquetStorage}` handle.
"""
function open_storage(
    output_path::AbstractString = joinpath("outputs", "larval_trajectories.parquet");
    storage_backend::Union{Nothing, String} = nothing,
    read_only::Bool = false,
    max_retries::Int = 5,
    retry_delay::Real = 0.5
)
    backend = if !isnothing(storage_backend)
        lowercase(storage_backend)
    elseif endswith(output_path, ".parquet")
        "geoparquet"
    elseif endswith(output_path, ".zarr")
        "zarr"
    else
        "geoparquet"
    end

    if backend == "zarr"
        storage_path = endswith(output_path, ".zarr") ? String(output_path) :
            (endswith(output_path, ".jld2") ?
             replace(output_path, r"\.jld2$" => ".zarr") : output_path * ".zarr")
        st = GeoData.open_geostorage(storage_path; backend = :zarr, read_only = read_only)
        return st.handle
    elseif backend == "geoparquet"
        storage_path = endswith(output_path, ".parquet") ? String(output_path) :
            (endswith(output_path, ".jld2") ?
             replace(output_path, r"\.jld2$" => ".parquet") :
             (endswith(output_path, ".zarr") ?
              replace(output_path, r"\.zarr$" => ".parquet") : output_path * ".parquet"))
        return GeoParquetStorage(storage_path; read_only = read_only)
    else
        error("Unsupported storage backend: $backend. Supported: 'geoparquet', 'zarr'")
    end
end

"""
    close_storage(group::Union{Zarr.ZGroup, GeoParquetStorage}; force::Bool = false)

Close the storage group and flush pending data.
"""
function close_storage(group::Union{Zarr.ZGroup, GeoParquetStorage}; force::Bool = false)
    GC.gc()
    return nothing
end

"""
    close_all_storage!()

Explicitly checkpoint, close, and run garbage collection.
"""
function close_all_storage!()
    GC.gc()
    return nothing
end

"""
    initialize_storage_schema!(group::Union{Zarr.ZGroup, GeoParquetStorage})

Initialize schema hierarchy for runs, trajectories, and gridded fields.
"""
function initialize_storage_schema!(group::Union{Zarr.ZGroup, GeoParquetStorage})
    if group isa Zarr.ZGroup
        for gname in ("metadata", "trajectories", "gridded_fields",
                      "hydrodynamic_fields", "metrics", "connectivity", "gridded_dispersal")
            if !haskey(group, gname)
                Zarr.zgroup(group, gname)
            end
        end
        if !haskey(group, "runs")
            arr = Zarr.zcreate(
                String, group, "runs", 0;
                chunks = (1000,),
                compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
            )
        end
    elseif group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        mkpath(base_dir)
    end
    return nothing
end

# ============================================================================
# Simulation Run Archival & Querying
# ============================================================================

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

Persist a complete simulation run including parameters, spatiotemporal trajectories,
recruitment metrics, and demographic connectivity matrices.
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
    canonical_trajs = canonicalize_trajectories(trajectories)

    scenario = string(hasproperty(opts, :scenario) ? opts.scenario : :baseline)
    proj_year = Int(hasproperty(opts, :projection_year) ? opts.projection_year : 2050)
    n_parts = Int(hasproperty(opts, :n_particles) ? opts.n_particles : size(canonical_trajs.lons, 1))
    duration = Float64(hasproperty(opts, :track_duration) ? opts.track_duration : (canonical_trajs.times[end] - canonical_trajs.times[1]))
    dt_val = Float64(hasproperty(opts, :track_dt) ? opts.track_dt : (length(canonical_trajs.times) > 1 ? canonical_trajs.times[2] - canonical_trajs.times[1] : 300.0))
    tides = Bool(hasproperty(opts, :enable_tides) ? opts.enable_tides : false)
    dvm = Bool(hasproperty(opts, :enable_dvm) ? opts.enable_dvm : true)
    molt = Bool(hasproperty(opts, :enable_molting) ? opts.enable_molting : true)
    diff_h = Float64(hasproperty(opts, :diffusivity_h) ? opts.diffusivity_h : 10.0)
    diff_v = Float64(hasproperty(opts, :diffusivity_v) ? opts.diffusivity_v : 1e-4)
    min_depth = Float64(hasproperty(opts, :min_seabed_depth) ? opts.min_seabed_depth : 100.0)
    seed_val = Int(hasproperty(opts, :seed) ? opts.seed : 42)
    created_at = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    rec_rate = (!isnothing(metrics) && hasproperty(metrics, :settlement_success_rate)) ?
        Float64(metrics.settlement_success_rate) : 0.0
    mean_pld = (!isnothing(metrics) && hasproperty(metrics, :mean_pld_days)) ?
        Float64(metrics.mean_pld_days) : 0.0
    mean_temp = (!isnothing(metrics) && hasproperty(metrics, :mean_exposure_temperature)) ?
        Float64(metrics.mean_exposure_temperature) : 0.0
    mean_disp = (!isnothing(metrics) && hasproperty(metrics, :mean_dispersal_distance_km)) ?
        Float64(metrics.mean_dispersal_distance_km) : 0.0
    mean_dd = (!isnothing(metrics) && hasproperty(metrics, :mean_degree_days)) ?
        Float64(metrics.mean_degree_days) : 0.0
    mort_therm = (!isnothing(metrics) && hasproperty(metrics, :thermal_mortality_rate)) ?
        Float64(metrics.thermal_mortality_rate) : 0.0

    config_toml_str = if !isnothing(config)
        if config isa AbstractString
            String(config)
        else
            sprint(TOML.print, config; sorted = true)
        end
    else
        ""
    end

    if group isa Zarr.ZGroup
        initialize_storage_schema!(group)

        run_metadata = Dict(
            "run_id" => run_id,
            "scenario" => scenario,
            "projection_year" => proj_year,
            "created_at" => created_at,
            "n_particles" => n_parts,
            "duration_seconds" => duration,
            "dt_seconds" => dt_val,
            "enable_tides" => tides,
            "enable_dvm" => dvm,
            "enable_molting" => molt,
            "diffusivity_h" => diff_h,
            "diffusivity_v" => diff_v,
            "min_seabed_depth" => min_depth,
            "seed" => seed_val,
            "settlement_success_rate" => rec_rate,
            "mean_pld_days" => mean_pld,
            "mean_exposure_temperature" => mean_temp,
            "mean_dispersal_distance_km" => mean_disp,
            "mean_degree_days" => mean_dd,
            "mean_temperature_celsius" => mean_temp,
            "mean_dispersal_km" => mean_disp,
            "mean_settlement_success" => rec_rate,
            "mean_thermal_mortality_rate" => mort_therm,
            "config_toml" => config_toml_str,
            "notes" => notes
        )

        runs_arr = group["runs"]
        existing_runs = try
            runs_arr[:]
        catch
            String[]
        end
        push!(existing_runs, JSON3.write(run_metadata))
        Zarr.resize!(runs_arr, length(existing_runs))
        runs_arr[:] = existing_runs

        traj_group = haskey(group, "trajectories") ? group["trajectories"] : Zarr.zgroup(group, "trajectories")
        run_traj = haskey(traj_group, run_id) ? traj_group[run_id] : Zarr.zgroup(traj_group, run_id)

        st_strs = [string(s) for s in canonical_trajs.stages]
        set_strs = [string(s) for s in canonical_trajs.settlement_status]
        al_vec = Bool.(canonical_trajs.alive)

        for (name, data) in (
            "lons" => canonical_trajs.lons,
            "lats" => canonical_trajs.lats,
            "depths" => canonical_trajs.depths,
            "temperatures" => canonical_trajs.temperatures,
            "degree_days" => canonical_trajs.degree_days_timeseries,
            "degree_days_timeseries" => canonical_trajs.degree_days_timeseries,
            "survival_probability" => canonical_trajs.survival_probability,
            "stages" => st_strs,
            "alive" => al_vec,
            "settlement_status" => set_strs,
            "settlement_age" => get(canonical_trajs, :settlement_age, nothing),
            "ids" => canonical_trajs.ids,
            "times" => canonical_trajs.times
        )
            if !isnothing(data)
                if haskey(run_traj, name)
                    run_traj[name][:] = data
                else
                    arr = Zarr.zcreate(
                        eltype(data), run_traj, name, size(data)...;
                        chunks = size(data),
                        compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
                    )
                    arr[:] = data
                end
            end
        end

        if !isnothing(metrics)
            metrics_group = haskey(group, "metrics") ? group["metrics"] : Zarr.zgroup(group, "metrics")
            rg = haskey(metrics_group, run_id) ? metrics_group[run_id] : Zarr.zgroup(metrics_group, run_id)
            for (k, v) in pairs(metrics)
                if v isa Number
                    if haskey(rg, string(k))
                        rg[string(k)][:] = [Float64(v)]
                    else
                        arr = Zarr.zcreate(Float64, rg, string(k), 1; chunks = (1,))
                        arr[:] = [Float64(v)]
                    end
                end
            end
        end

        if !isnothing(connectivity)
            conn_group = haskey(group, "connectivity") ? group["connectivity"] : Zarr.zgroup(group, "connectivity")
            rg = haskey(conn_group, run_id) ? conn_group[run_id] : Zarr.zgroup(conn_group, run_id)
            mat = Float64.(connectivity.matrix)
            if haskey(rg, "matrix")
                rg["matrix"][:] = mat
            else
                arr = Zarr.zcreate(Float64, rg, "matrix", size(mat)...; chunks = size(mat))
                arr[:] = mat
            end
            strata_strs = [string(s) for s in connectivity.strata_names]
            if haskey(rg, "strata_names")
                rg["strata_names"][:] = strata_strs
            else
                arr = Zarr.zcreate(String, rg, "strata_names", length(strata_strs); chunks = (length(strata_strs),))
                arr[:] = strata_strs
            end
        end

        if !isnothing(gridded_dispersal)
            grid_group = haskey(group, "gridded_dispersal") ? group["gridded_dispersal"] : Zarr.zgroup(group, "gridded_dispersal")
            rg = haskey(grid_group, run_id) ? grid_group[run_id] : Zarr.zgroup(grid_group, run_id)
            for (k, v) in pairs(gridded_dispersal)
                if v isa AbstractArray
                    if haskey(rg, string(k))
                        rg[string(k)][:] = v
                    else
                        arr = Zarr.zcreate(eltype(v), rg, string(k), size(v)...; chunks = size(v))
                        arr[:] = v
                    end
                end
            end
        end

        return String(run_id)
    elseif group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        mkpath(base_dir)

        n_p = size(canonical_trajs.lons, 1)
        n_t = size(canonical_trajs.lons, 2)
        n_total = n_p * n_t

        p_ids = canonical_trajs.ids
        times = canonical_trajs.times
        lons = canonical_trajs.lons
        lats = canonical_trajs.lats
        depths = canonical_trajs.depths
        temps = canonical_trajs.temperatures
        dds = canonical_trajs.degree_days_timeseries
        stages = canonical_trajs.stages
        alive = canonical_trajs.alive
        settle = canonical_trajs.settlement_status

        col_run_id = fill(String(run_id), n_total)
        col_pid = Vector{Int}(undef, n_total)
        col_step = Vector{Int}(undef, n_total)
        col_time = Vector{Float64}(undef, n_total)
        col_lon = Vector{Float64}(undef, n_total)
        col_lat = Vector{Float64}(undef, n_total)
        col_depth = Vector{Float64}(undef, n_total)
        col_temp = Vector{Float64}(undef, n_total)
        col_dd = Vector{Float64}(undef, n_total)
        col_stage = Vector{String}(undef, n_total)
        col_alive = Vector{Bool}(undef, n_total)
        col_settle = Vector{String}(undef, n_total)

        idx = 1
        for p in 1:n_p
            pid = p_ids[p]
            p_alive = alive isa AbstractMatrix ? alive[p, 1] : alive[p]
            p_settle = settle isa AbstractMatrix ? string(settle[p, end]) : string(settle[p])
            for t in 1:n_t
                col_pid[idx] = pid
                col_step[idx] = t - 1
                col_time[idx] = times[t]
                ln = Float64(lons[p, t])
                lt = Float64(lats[p, t])
                col_lon[idx] = ln
                col_lat[idx] = lt
                col_depth[idx] = Float64(depths[p, t])
                col_temp[idx] = Float64(temps[p, t])
                col_dd[idx] = Float64(dds[p, t])
                col_stage[idx] = stages isa AbstractMatrix ? string(stages[p, t]) : string(stages[p])
                col_alive[idx] = alive isa AbstractMatrix ? Bool(alive[p, t]) : Bool(p_alive)
                col_settle[idx] = settle isa AbstractMatrix ? string(settle[p, t]) : string(p_settle)
                idx += 1
            end
        end

        col_geom = [GeoInterface.Point(col_lon[i], col_lat[i]) for i in 1:n_total]

        df_trajs = DataFrame(
            run_id = col_run_id,
            particle_id = col_pid,
            step = col_step,
            time_seconds = col_time,
            lon = col_lon,
            lat = col_lat,
            longitude = col_lon,
            latitude = col_lat,
            depth = col_depth,
            temperature = col_temp,
            degree_days = col_dd,
            survival_probability = fill(1.0, n_total),
            stage = col_stage,
            alive = col_alive,
            settlement_status = col_settle,
            geometry = col_geom
        )

        traj_path = joinpath(base_dir, "trajectories_$(run_id).parquet")
        GeoData.GeoParquet.write(traj_path, df_trajs, (:geometry,))

        if endswith(group.path, ".parquet")
            GeoData.GeoParquet.write(group.path, df_trajs, (:geometry,))
        end

        centroid_lon = isempty(col_lon) ? 0.0 : sum(col_lon) / length(col_lon)
        centroid_lat = isempty(col_lat) ? 0.0 : sum(col_lat) / length(col_lat)
        run_row = DataFrame(
            run_id = [String(run_id)],
            scenario = [scenario],
            projection_year = [proj_year],
            n_particles = [Int(n_parts)],
            duration_days = [Float64(duration / 86400.0)],
            track_dt_seconds = [Float64(dt_val)],
            settlement_success_rate = [rec_rate],
            mean_pld_days = [mean_pld],
            mean_exposure_temperature = [mean_temp],
            mean_dispersal_distance_km = [mean_disp],
            mean_degree_days = [mean_dd],
            mean_temperature_celsius = [mean_temp],
            mean_dispersal_km = [mean_disp],
            mean_settlement_success = [rec_rate],
            mean_thermal_mortality_rate = [mort_therm],
            n_runs = [1],
            created_at = [created_at],
            config_toml = [config_toml_str],
            notes = [String(notes)],
            geometry = [GeoInterface.Point(centroid_lon, centroid_lat)]
        )

        if isfile(runs_path)
            try
                existing_runs = copy(GeoData.GeoParquet.read(runs_path))
                filter!(r -> r.run_id != run_id, existing_runs)
                append!(existing_runs, run_row)
                GC.gc()
                GeoData.GeoParquet.write(runs_path, existing_runs, (:geometry,))
            catch
                GeoData.GeoParquet.write(runs_path, run_row, (:geometry,))
            end
        else
            GeoData.GeoParquet.write(runs_path, run_row, (:geometry,))
        end

        if !isnothing(connectivity) && hasproperty(connectivity, :matrix)
            conn_path = joinpath(base_dir, "connectivity_$(run_id).parquet")
            strata = hasproperty(connectivity, :strata_names) ?
                [string(s) for s in connectivity.strata_names] :
                ["stratum_$i" for i in 1:size(connectivity.matrix, 1)]
            n_s = length(strata)
            mat = connectivity.matrix
            c_src = String[]
            c_dst = String[]
            c_prob = Float64[]
            for i in 1:n_s, j in 1:n_s
                push!(c_src, strata[i])
                push!(c_dst, strata[j])
                push!(c_prob, Float64(mat[i, j]))
            end
            conn_geom = [GeoInterface.Point(centroid_lon, centroid_lat) for _ in 1:length(c_src)]
            conn_df = DataFrame(source = c_src, target = c_dst, probability = c_prob, geometry = conn_geom)
            GeoData.GeoParquet.write(conn_path, conn_df, (:geometry,))
        end

        if !isnothing(gridded_dispersal)
            grid_path = joinpath(base_dir, "gridded_dispersal_$(run_id).parquet")
            lons_c = Float64.(gridded_dispersal.lon_centers)
            lats_c = Float64.(gridded_dispersal.lat_centers)
            nx, ny = length(lons_c), length(lats_c)
            g_lons = Vector{Float64}(undef, nx * ny)
            g_lats = Vector{Float64}(undef, nx * ny)
            g_u = Vector{Float64}(undef, nx * ny)
            g_v = Vector{Float64}(undef, nx * ny)
            g_diff = Vector{Float64}(undef, nx * ny)
            g_dens = Vector{Float64}(undef, nx * ny)

            g_idx = 1
            for i in 1:nx, j in 1:ny
                ln = lons_c[i]
                lt = lats_c[j]
                g_lons[g_idx] = ln
                g_lats[g_idx] = lt
                g_u[g_idx] = Float64(gridded_dispersal.u_mean[i, j])
                g_v[g_idx] = Float64(gridded_dispersal.v_mean[i, j])
                g_diff[g_idx] = Float64(gridded_dispersal.diffusivity[i, j])
                g_dens[g_idx] = Float64(gridded_dispersal.settlement_density[i, j])
                g_idx += 1
            end
            g_geom = [GeoInterface.Point(g_lons[k], g_lats[k]) for k in 1:(nx * ny)]
            grid_df = DataFrame(
                lon = g_lons,
                lat = g_lats,
                empirical_u = g_u,
                empirical_v = g_v,
                empirical_diffusivity = g_diff,
                settlement_density = g_dens,
                geometry = g_geom
            )
            GeoData.GeoParquet.write(grid_path, grid_df, (:geometry,))
        end

        return String(run_id)
    end
end

"""
    list_simulation_runs(
        group::Union{Zarr.ZGroup, GeoParquetStorage};
        scenario::Union{Nothing, AbstractString, Symbol} = nothing,
        projection_year::Union{Nothing, Int} = nothing
    )::DataFrame

Query and return a summary `DataFrame` of all simulation runs in storage.
"""
function list_simulation_runs(
    group::Union{Zarr.ZGroup, GeoParquetStorage};
    scenario::Union{Nothing, AbstractString, Symbol} = nothing,
    projection_year::Union{Nothing, Int} = nothing
)::DataFrame
    if group isa Zarr.ZGroup
        if !haskey(group, "runs")
            return DataFrame()
        end
        runs_arr = group["runs"]
        runs_json = try
            runs_arr[:]
        catch
            String[]
        end
        runs = DataFrame()
        for r_json in runs_json
            r = JSON3.read(r_json)
            if !isnothing(scenario) && r["scenario"] != string(scenario)
                continue
            end
            if !isnothing(projection_year) && r["projection_year"] != projection_year
                continue
            end
            push!(runs, r)
        end
        isempty(runs) || sort!(runs, :created_at, rev = true)
        return runs
    elseif group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        runs_path = joinpath(base_dir, "runs.parquet")
        if !isfile(runs_path)
            return DataFrame()
        end
        runs = GeoData.GeoParquet.read(runs_path)
        if !isnothing(scenario)
            sc_str = string(scenario)
            filter!(r -> r.scenario == sc_str, runs)
        end
        if !isnothing(projection_year)
            filter!(r -> r.projection_year == projection_year, runs)
        end
        if "created_at" in names(runs) && nrow(runs) > 1
            sort!(runs, :created_at, rev = true)
        end
        return runs
    else
        return DataFrame()
    end
end

"""
    load_run_configuration(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    ) -> Union{Nothing, AbstractDict}

Load configuration dictionary for a given run from storage.
"""
function load_run_configuration(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString
)
    runs = list_simulation_runs(group)
    for row in eachrow(runs)
        if row.run_id == run_id && hasproperty(row, :config_toml) && !isempty(row.config_toml)
            try
                return TOML.parse(row.config_toml)
            catch
                return nothing
            end
        end
    end
    return nothing
end

"""
    load_trajectories_namedtuple(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString;
        max_particles::Union{Nothing, Int} = nothing
    )::NamedTuple

Retrieve Lagrangian particle trajectories reconstructed as a `NamedTuple`.
"""
function load_trajectories_namedtuple(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString;
    max_particles::Union{Nothing, Int} = nothing
)::NamedTuple
    if group isa Zarr.ZGroup
        traj_group = group["trajectories"]
        tg = haskey(traj_group, run_id) ? traj_group[run_id] : traj_group
        lons = tg["lons"][:]
        lats = tg["lats"][:]
        depths = tg["depths"][:]
        temperatures = tg["temperatures"][:]
        degree_days = haskey(tg, "degree_days") ? tg["degree_days"][:] : zeros(size(lons))
        degree_days_ts = haskey(tg, "degree_days_timeseries") ? tg["degree_days_timeseries"][:] : degree_days
        survival_prob = haskey(tg, "survival_probability") ? tg["survival_probability"][:] : ones(size(lons))
        stages = tg["stages"][:]
        alive = tg["alive"][:]
        settle_status = tg["settlement_status"][:]
        ids = tg["ids"][:]
        times = tg["times"][:]

        n_p = size(lons, 1)
        p_ids = collect(1:n_p)
        if !isnothing(max_particles) && max_particles > 0
            p_ids = p_ids[1:min(max_particles, n_p)]
        end

        al_sub = alive isa AbstractMatrix ? alive[p_ids, :] : alive[p_ids]
        set_sub = settle_status isa AbstractMatrix ? settle_status[p_ids, :] : settle_status[p_ids]

        return (
            lons = lons[p_ids, :],
            lats = lats[p_ids, :],
            depths = depths[p_ids, :],
            temperatures = temperatures[p_ids, :],
            degree_days = degree_days[p_ids, :],
            degree_days_timeseries = degree_days_ts[p_ids, :],
            survival_probability = survival_prob[p_ids, :],
            stages = stages[p_ids, :],
            alive = al_sub,
            settlement_status = set_sub,
            settlement_age = haskey(tg, "settlement_age") ? tg["settlement_age"][:][p_ids] : nothing,
            times = times,
            ids = ids[p_ids]
        )
    elseif group isa GeoParquetStorage
        df = load_trajectories_df(group, run_id; max_particles = max_particles)
        if nrow(df) == 0
            return (
                lons = Matrix{Float64}(undef, 0, 0),
                lats = Matrix{Float64}(undef, 0, 0),
                depths = Matrix{Float64}(undef, 0, 0),
                temperatures = Matrix{Float64}(undef, 0, 0),
                degree_days = Matrix{Float64}(undef, 0, 0),
                degree_days_timeseries = Matrix{Float64}(undef, 0, 0),
                survival_probability = Matrix{Float64}(undef, 0, 0),
                stages = Matrix{String}(undef, 0, 0),
                alive = Bool[],
                settlement_status = String[],
                settlement_age = nothing,
                times = Float64[],
                ids = Int[]
            )
        end

        p_ids = sort(unique(df.particle_id))
        times = sort(unique(df.time_seconds))
        n_p = length(p_ids)
        n_t = length(times)

        pid_map = Dict(pid => i for (i, pid) in enumerate(p_ids))
        time_map = Dict(t => j for (j, t) in enumerate(times))

        lons = fill(NaN, n_p, n_t)
        lats = fill(NaN, n_p, n_t)
        depths = fill(NaN, n_p, n_t)
        temps = fill(NaN, n_p, n_t)
        dds = fill(NaN, n_p, n_t)
        stages = fill("", n_p, n_t)
        alive_vec = fill(true, n_p)
        settle_vec = fill("pelagic", n_p)

        for row in eachrow(df)
            i = pid_map[row.particle_id]
            j = time_map[row.time_seconds]
            lons[i, j] = row.longitude
            lats[i, j] = row.latitude
            depths[i, j] = row.depth
            temps[i, j] = row.temperature
            dds[i, j] = row.degree_days
            stages[i, j] = row.stage
            alive_vec[i] = row.alive
            settle_vec[i] = row.settlement_status
        end

        return (
            lons = lons,
            lats = lats,
            depths = depths,
            temperatures = temps,
            degree_days = dds,
            degree_days_timeseries = dds,
            survival_probability = fill(1.0, n_p, n_t),
            stages = stages,
            alive = alive_vec,
            settlement_status = settle_vec,
            settlement_age = nothing,
            times = times,
            ids = p_ids
        )
    else
        return (
            lons = Float64[], lats = Float64[], depths = Float64[],
            temperatures = Float64[], degree_days = Float64[],
            degree_days_timeseries = Float64[], survival_probability = Float64[],
            stages = String[], alive = Bool[], settlement_status = String[],
            settlement_age = Float64[], times = Float64[], ids = Int64[]
        )
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

Retrieve particle trajectory records from storage as a tidy `DataFrame`.
"""
function load_trajectories_df(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString;
    particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
    stage::Union{Nothing, Symbol, AbstractString} = nothing,
    time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
    max_particles::Union{Nothing, Int} = nothing
)::DataFrame
    if group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        traj_path = joinpath(base_dir, "trajectories_$(run_id).parquet")
        target_file = isfile(traj_path) ? traj_path :
            (isfile(group.path) ? group.path : "")
        if !isempty(target_file) && isfile(target_file)
            df = copy(GeoData.GeoParquet.read(target_file))
            if "run_id" in names(df)
                filter!(row -> row.run_id == run_id, df)
            end
            if !isnothing(particle_ids)
                filter!(row -> row.particle_id in particle_ids, df)
            end
            if !isnothing(stage)
                st_str = string(stage)
                filter!(row -> row.stage == st_str, df)
            end
            if !isnothing(time_range)
                filter!(row -> time_range[1] <= row.time_seconds <= time_range[2], df)
            end
            if !isnothing(max_particles) && max_particles > 0
                unique_pids = unique(df.particle_id)
                if length(unique_pids) > max_particles
                    keep_pids = Set(unique_pids[1:max_particles])
                    filter!(row -> row.particle_id in keep_pids, df)
                end
            end
            return df
        end
    end

    trajs = load_trajectories_namedtuple(group, run_id; max_particles = max_particles)
    n_p = size(trajs.lons, 1)
    n_t = size(trajs.lons, 2)
    n_total = n_p * n_t

    df_run_id = fill(run_id, n_total)
    df_p_id = Vector{Int}(undef, n_total)
    df_step = Vector{Int}(undef, n_total)
    df_time = Vector{Float64}(undef, n_total)
    df_lon = Vector{Float64}(undef, n_total)
    df_lat = Vector{Float64}(undef, n_total)
    df_depth = Vector{Float64}(undef, n_total)
    df_temp = Vector{Float64}(undef, n_total)
    df_dd = Vector{Float64}(undef, n_total)
    df_surv = Vector{Float64}(undef, n_total)
    df_stage = Vector{String}(undef, n_total)
    df_alive = Vector{Bool}(undef, n_total)
    df_settle = Vector{String}(undef, n_total)
    df_geom = Vector{GeoInterface.Point{false, false, Float64}}(undef, n_total)

    idx = 1
    for p in 1:n_p
        p_id = trajs.ids[p]
        p_alive = trajs.alive isa AbstractMatrix ? trajs.alive[p, 1] : trajs.alive[p]
        p_settle = trajs.settlement_status isa AbstractMatrix ?
            string(trajs.settlement_status[p, end]) : string(trajs.settlement_status[p])

        for s in 1:n_t
            df_p_id[idx] = p_id
            df_step[idx] = s - 1
            df_time[idx] = trajs.times[s]
            ln = trajs.lons[p, s]
            lt = trajs.lats[p, s]
            df_lon[idx] = ln
            df_lat[idx] = lt
            df_depth[idx] = trajs.depths[p, s]
            df_temp[idx] = trajs.temperatures[p, s]
            df_dd[idx] = trajs.degree_days[p, s]
            df_surv[idx] = trajs.survival_probability[p, s]
            df_stage[idx] = string(trajs.stages[p, s])
            df_alive[idx] = p_alive
            df_settle[idx] = p_settle
            df_geom[idx] = GeoInterface.Point(ln, lt)
            idx += 1
        end
    end

    df = DataFrame(
        run_id = df_run_id,
        particle_id = df_p_id,
        step = df_step,
        time_seconds = df_time,
        lon = df_lon,
        lat = df_lat,
        longitude = df_lon,
        latitude = df_lat,
        depth = df_depth,
        temperature = df_temp,
        degree_days = df_dd,
        survival_probability = df_surv,
        stage = df_stage,
        alive = df_alive,
        settlement_status = df_settle,
        geometry = df_geom
    )

    if !isnothing(particle_ids)
        filter!(row -> row.particle_id in particle_ids, df)
    end
    if !isnothing(stage)
        st_str = string(stage)
        filter!(row -> row.stage == st_str, df)
    end
    if !isnothing(time_range)
        filter!(row -> time_range[1] <= row.time_seconds <= time_range[2], df)
    end

    return df
end

"""
    load_all_scenario_trajectories(
        group::Union{Zarr.ZGroup, GeoParquetStorage};
        projection_year::Union{Nothing, Int} = nothing
    ) -> Dict{Symbol, NamedTuple}

Load Lagrangian trajectories for all scenarios present in storage.
"""
function load_all_scenario_trajectories(
    group::Union{Zarr.ZGroup, GeoParquetStorage};
    projection_year::Union{Nothing, Int} = nothing
)
    runs_df = list_simulation_runs(group; projection_year = projection_year)
    if nrow(runs_df) == 0
        return Dict{Symbol, NamedTuple}()
    end

    scenario_dict = Dict{Symbol, NamedTuple}()
    for scen_name in unique(runs_df.scenario)
        subset_df = filter(r -> r.scenario == scen_name, runs_df)
        latest_run_id = string(subset_df.run_id[1])
        try
            scenario_dict[Symbol(scen_name)] = load_trajectories_namedtuple(group, latest_run_id)
        catch err
            @warn "Failed to load trajectories for scenario $(scen_name): $(err)"
        end
    end

    return scenario_dict
end

"""
    save_hydrodynamic_field!(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString,
        opts;
        grid_lons::Union{Nothing, AbstractVector} = nothing,
        grid_lats::Union{Nothing, AbstractVector} = nothing,
        grid_depths::Union{Nothing, AbstractVector} = nothing,
        u = nothing,
        v = nothing,
        w = nothing,
        temperature = nothing,
        salinity = nothing,
        elevation = nothing,
        time_seconds::Real = 0.0,
        field_data::Union{Nothing, NamedTuple} = nothing
    )

Save Eulerian hydrodynamic fields (velocity components, tracers, surface elevation)
to storage.
"""
function save_hydrodynamic_field!(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString,
    opts;
    grid_lons::Union{Nothing, AbstractVector} = nothing,
    grid_lats::Union{Nothing, AbstractVector} = nothing,
    grid_depths::Union{Nothing, AbstractVector} = nothing,
    u = nothing,
    v = nothing,
    w = nothing,
    temperature = nothing,
    salinity = nothing,
    elevation = nothing,
    time_seconds::Real = 0.0,
    field_data::Union{Nothing, NamedTuple} = nothing
)
    if group isa Zarr.ZGroup
        hydro_group = haskey(group, "hydrodynamic_fields") ? group["hydrodynamic_fields"] : Zarr.zgroup(group, "hydrodynamic_fields")
        run_group = haskey(hydro_group, run_id) ? hydro_group[run_id] : Zarr.zgroup(hydro_group, run_id)

        for (k, val) in (
            :lons => grid_lons, :lats => grid_lats, :depths => grid_depths,
            :u => u, :v => v, :w => w, :temperature => temperature,
            :salinity => salinity, :elevation => elevation
        )
            if !isnothing(val)
                name = string(k)
                if haskey(run_group, name)
                    run_group[name][:] = val
                else
                    arr = Zarr.zcreate(eltype(val), run_group, name, size(val)...; chunks = size(val))
                    arr[:] = val
                end
            end
        end

        if !isnothing(field_data)
            for (k, v) in pairs(field_data)
                if v isa AbstractArray
                    name = string(k)
                    if haskey(run_group, name)
                        run_group[name][:] = v
                    else
                        arr = Zarr.zcreate(eltype(v), run_group, name, size(v)...; chunks = size(v))
                        arr[:] = v
                    end
                end
            end
        end
    elseif group isa GeoParquetStorage
        group.datasets["hydro_$(run_id)"] = (
            lons = grid_lons, lats = grid_lats, depths = grid_depths,
            u = u, v = v, w = w,
            temperature = temperature, salinity = salinity, elevation = elevation
        )
    end
    return nothing
end

"""
    load_hydrodynamic_field(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString;
        time_seconds::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = nothing
    ) -> NamedTuple

Load Eulerian hydrodynamic field arrays from storage.
"""
function load_hydrodynamic_field(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString;
    time_seconds::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = nothing
)::NamedTuple
    if group isa Zarr.ZGroup && haskey(group, "hydrodynamic_fields")
        hg = group["hydrodynamic_fields"]
        if haskey(hg, run_id)
            rg = hg[run_id]
            lons = haskey(rg, "lons") ? rg["lons"][:] : Float64[]
            lats = haskey(rg, "lats") ? rg["lats"][:] : Float64[]
            depths = haskey(rg, "depths") ? rg["depths"][:] : Float64[]
            u_arr = haskey(rg, "u") ? rg["u"][:] : zeros(Float64, length(lons), length(lats), length(depths))
            v_arr = haskey(rg, "v") ? rg["v"][:] : zeros(Float64, length(lons), length(lats), length(depths))
            w_arr = haskey(rg, "w") ? rg["w"][:] : zeros(Float64, length(lons), length(lats), length(depths))
            t_arr = haskey(rg, "temperature") ? rg["temperature"][:] : fill(4.0, length(lons), length(lats), length(depths))
            s_arr = haskey(rg, "salinity") ? rg["salinity"][:] : fill(32.5, length(lons), length(lats), length(depths))
            elev_arr = haskey(rg, "elevation") ? rg["elevation"][:] : zeros(Float64, length(lons), length(lats))
            return (
                lons = lons, lats = lats, depths = depths,
                u = u_arr, v = v_arr, w = w_arr,
                temperature = t_arr, salinity = s_arr, elevation = elev_arr
            )
        end
    elseif group isa GeoParquetStorage
        if haskey(group.datasets, "hydro_$(run_id)")
            return group.datasets["hydro_$(run_id)"]
        end
    end
    return (
        lons = Float64[], lats = Float64[], depths = Float64[],
        u = Array{Float64, 3}(undef, 0, 0, 0),
        v = Array{Float64, 3}(undef, 0, 0, 0),
        w = Array{Float64, 3}(undef, 0, 0, 0),
        temperature = Array{Float64, 3}(undef, 0, 0, 0),
        salinity = Array{Float64, 3}(undef, 0, 0, 0),
        elevation = Matrix{Float64}(undef, 0, 0)
    )
end

"""
    load_gridded_dispersal(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    ) -> NamedTuple

Load gridded dispersal and settlement density summary fields.
"""
function load_gridded_dispersal(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString
)::NamedTuple
    if group isa Zarr.ZGroup && haskey(group, "gridded_dispersal")
        gg = group["gridded_dispersal"]
        if haskey(gg, run_id)
            rg = gg[run_id]
            return (
                lon_centers = haskey(rg, "lon_centers") ? rg["lon_centers"][:] : Float64[],
                lat_centers = haskey(rg, "lat_centers") ? rg["lat_centers"][:] : Float64[],
                u_mean = haskey(rg, "u_mean") ? rg["u_mean"][:] : Matrix{Float64}(undef, 0, 0),
                v_mean = haskey(rg, "v_mean") ? rg["v_mean"][:] : Matrix{Float64}(undef, 0, 0),
                diffusivity = haskey(rg, "diffusivity") ? rg["diffusivity"][:] : Matrix{Float64}(undef, 0, 0),
                settlement_density = haskey(rg, "settlement_density") ? rg["settlement_density"][:] : Matrix{Float64}(undef, 0, 0),
                mean_exposure_temperature = haskey(rg, "mean_exposure_temperature") ? rg["mean_exposure_temperature"][:] : Matrix{Float64}(undef, 0, 0),
                mean_degree_days = haskey(rg, "mean_degree_days") ? rg["mean_degree_days"][:] : Matrix{Float64}(undef, 0, 0),
                sample_count = haskey(rg, "sample_count") ? rg["sample_count"][:] : Matrix{Int}(undef, 0, 0)
            )
        end
    elseif group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        g_path = joinpath(base_dir, "gridded_dispersal_$(run_id).parquet")
        if isfile(g_path)
            df = GeoData.GeoParquet.read(g_path)
            lons = sort(unique(df.lon))
            lats = sort(unique(df.lat))
            nx, ny = length(lons), length(lats)
            u_grid = fill(NaN, nx, ny)
            v_grid = fill(NaN, nx, ny)
            diff_grid = fill(NaN, nx, ny)
            dens_grid = fill(NaN, nx, ny)
            lon_to_i = Dict(lon => i for (i, lon) in enumerate(lons))
            lat_to_j = Dict(lat => j for (j, lat) in enumerate(lats))
            for row in eachrow(df)
                i = get(lon_to_i, row.lon, 0)
                j = get(lat_to_j, row.lat, 0)
                if i > 0 && j > 0
                    u_grid[i, j] = Float64(row.empirical_u)
                    v_grid[i, j] = Float64(row.empirical_v)
                    diff_grid[i, j] = Float64(row.empirical_diffusivity)
                    dens_grid[i, j] = Float64(row.settlement_density)
                end
            end
            return (
                lon_centers = lons,
                lat_centers = lats,
                u_mean = u_grid,
                v_mean = v_grid,
                diffusivity = diff_grid,
                settlement_density = dens_grid,
                mean_exposure_temperature = fill(NaN, nx, ny),
                mean_degree_days = fill(NaN, nx, ny),
                sample_count = zeros(Int, nx, ny)
            )
        end
    end
    return (
        lon_centers = Float64[],
        lat_centers = Float64[],
        u_mean = Matrix{Float64}(undef, 0, 0),
        v_mean = Matrix{Float64}(undef, 0, 0),
        diffusivity = Matrix{Float64}(undef, 0, 0),
        settlement_density = Matrix{Float64}(undef, 0, 0),
        mean_exposure_temperature = Matrix{Float64}(undef, 0, 0),
        mean_degree_days = Matrix{Float64}(undef, 0, 0),
        sample_count = Matrix{Int}(undef, 0, 0)
    )
end

"""
    load_connectivity_matrix(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        run_id::AbstractString
    ) -> NamedTuple

Retrieve transition probability connectivity matrix \$P_{ij}\$ and stratum names.
"""
function load_connectivity_matrix(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    run_id::AbstractString
)::NamedTuple
    if group isa Zarr.ZGroup && haskey(group, "connectivity")
        cg = group["connectivity"]
        if haskey(cg, run_id)
            rg = cg[run_id]
            return (
                matrix = rg["matrix"][:],
                strata_names = [string(s) for s in rg["strata_names"][:]],
                counts_unweighted = haskey(rg, "counts_unweighted") ? rg["counts_unweighted"][:] : nothing,
                counts_matrix = haskey(rg, "counts_matrix") ? rg["counts_matrix"][:] : nothing
            )
        end
    elseif group isa GeoParquetStorage
        base_dir = _storage_base_dir(group)
        conn_path = joinpath(base_dir, "connectivity_$(run_id).parquet")
        if isfile(conn_path)
            df = GeoData.GeoParquet.read(conn_path)
            sources = unique(df.source)
            targets = unique(df.target)
            strata = unique(vcat(sources, targets))
            n_s = length(strata)
            s_map = Dict(s => i for (i, s) in enumerate(strata))
            mat = zeros(Float64, n_s, n_s)
            for row in eachrow(df)
                i = s_map[row.source]
                j = s_map[row.target]
                mat[i, j] = Float64(row.probability)
            end
            return (
                matrix = mat,
                strata_names = strata,
                counts_unweighted = nothing,
                counts_matrix = nothing
            )
        end
    end
    return (
        matrix = Matrix{Float64}(undef, 0, 0),
        strata_names = String[],
        counts_unweighted = nothing,
        counts_matrix = nothing
    )
end

"""
    compare_scenarios(
        group::Union{Zarr.ZGroup, GeoParquetStorage};
        scenario_names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
        projection_years::Union{Nothing, AbstractVector{Int}} = nothing,
        scenarios::Union{Nothing, AbstractVector{Symbol}} = nothing,
        projection_year::Union{Nothing, Int} = nothing
    )::DataFrame

Compare aggregated recruitment and transport metrics across climate projection pathways.
"""
function compare_scenarios(
    group::Union{Zarr.ZGroup, GeoParquetStorage};
    scenario_names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    projection_years::Union{Nothing, AbstractVector{Int}} = nothing,
    scenarios::Union{Nothing, AbstractVector{Symbol}} = nothing,
    projection_year::Union{Nothing, Int} = nothing
)::DataFrame
    runs_df = list_simulation_runs(group)
    if nrow(runs_df) == 0
        return DataFrame()
    end

    scen_filter = if !isnothing(scenario_names)
        [string(s) for s in scenario_names]
    elseif !isnothing(scenarios)
        [string(s) for s in scenarios]
    else
        nothing
    end
    if !isnothing(scen_filter)
        filter!(r -> r.scenario in scen_filter, runs_df)
    end

    year_filter = if !isnothing(projection_years)
        projection_years
    elseif !isnothing(projection_year)
        [projection_year]
    else
        nothing
    end
    if !isnothing(year_filter)
        filter!(r -> r.projection_year in year_filter, runs_df)
    end

    if nrow(runs_df) == 0
        return DataFrame()
    end

    gdf = groupby(runs_df, [:scenario, :projection_year])
    out = combine(gdf) do sdf
        rec = hasproperty(sdf, :settlement_success_rate) ? sdf.settlement_success_rate :
              (hasproperty(sdf, :mean_settlement_success) ? sdf.mean_settlement_success : [0.0])
        pld = hasproperty(sdf, :mean_pld_days) ? sdf.mean_pld_days : [0.0]
        dd = hasproperty(sdf, :mean_degree_days) ? sdf.mean_degree_days : [0.0]
        temp = hasproperty(sdf, :mean_exposure_temperature) ? sdf.mean_exposure_temperature :
               (hasproperty(sdf, :mean_temperature_celsius) ? sdf.mean_temperature_celsius : [0.0])
        disp = hasproperty(sdf, :mean_dispersal_distance_km) ? sdf.mean_dispersal_distance_km :
               (hasproperty(sdf, :mean_dispersal_km) ? sdf.mean_dispersal_km : [0.0])
        mort = hasproperty(sdf, :thermal_mortality_rate) ? sdf.thermal_mortality_rate :
               (hasproperty(sdf, :mean_thermal_mortality_rate) ? sdf.mean_thermal_mortality_rate : [0.0])

        (
            n_runs = nrow(sdf),
            mean_settlement_success = mean(rec),
            std_settlement_success = length(rec) > 1 ? std(rec) : 0.0,
            mean_pld_days = mean(pld),
            std_pld_days = length(pld) > 1 ? std(pld) : 0.0,
            mean_degree_days = mean(dd),
            mean_temperature_celsius = mean(temp),
            mean_dispersal_km = mean(disp),
            mean_thermal_mortality_rate = mean(mort)
        )
    end
    sort!(out, [:scenario, :projection_year])
    return out
end

"""
    compute_ensemble_model_average(
        group::Union{Zarr.ZGroup, GeoParquetStorage},
        scenario_names::AbstractVector{<:AbstractString};
        weights::Union{Nothing, AbstractVector{<:Real}} = nothing
    ) -> NamedTuple

Compute weighted multi-model ensemble mean and variance across climate scenarios.
"""
function compute_ensemble_model_average(
    group::Union{Zarr.ZGroup, GeoParquetStorage},
    scenario_names::AbstractVector{<:AbstractString};
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing
)
    n_scen = length(scenario_names)
    if n_scen == 0
        error("At least one scenario name must be provided for ensemble averaging.")
    end

    w_vec = if isnothing(weights)
        fill(1.0 / n_scen, n_scen)
    else
        if length(weights) != n_scen
            error("Length of weights ($(length(weights))) must match scenario_names ($(n_scen)).")
        end
        w_sum = sum(weights)
        if w_sum <= 0.0
            error("Sum of weights must be positive.")
        end
        Float64.(weights) ./ w_sum
    end

    matrices = Matrix{Float64}[]
    rec_rates = Float64[]
    pld_vals = Float64[]
    mort_vals = Float64[]
    strata_names_ref = String[]

    for (s_idx, scen) in enumerate(scenario_names)
        runs_df = list_simulation_runs(group; scenario = scen)
        if nrow(runs_df) == 0
            error("No runs found in storage for scenario: '$(scen)'.")
        end
        latest_run_id = string(runs_df.run_id[1])
        conn = load_connectivity_matrix(group, latest_run_id)

        if s_idx == 1
            strata_names_ref = [string(s) for s in conn.strata_names]
        end

        push!(matrices, conn.matrix)
        rec = hasproperty(runs_df, :settlement_success_rate) ? Float64(runs_df.settlement_success_rate[1]) :
              (hasproperty(runs_df, :mean_settlement_success) ? Float64(runs_df.mean_settlement_success[1]) : 0.0)
        pld = hasproperty(runs_df, :mean_pld_days) ? Float64(runs_df.mean_pld_days[1]) : 0.0
        temp = hasproperty(runs_df, :mean_exposure_temperature) ? Float64(runs_df.mean_exposure_temperature[1]) :
               (hasproperty(runs_df, :mean_temperature_celsius) ? Float64(runs_df.mean_temperature_celsius[1]) : 0.0)

        push!(rec_rates, rec)
        push!(pld_vals, pld)
        push!(mort_vals, temp)
    end

    n_strata = length(strata_names_ref)
    mean_conn = zeros(Float64, n_strata, n_strata)
    var_conn = zeros(Float64, n_strata, n_strata)

    for m in 1:n_scen
        mean_conn .+= w_vec[m] .* matrices[m]
    end

    for m in 1:n_scen
        var_conn .+= w_vec[m] .* ((matrices[m] .- mean_conn) .^ 2)
    end
    std_conn = sqrt.(max.(0.0, var_conn))

    ens_rec = sum(w_vec .* rec_rates)
    ens_pld = sum(w_vec .* pld_vals)
    ens_mort = sum(w_vec .* mort_vals)

    return (
        mean_connectivity = mean_conn,
        std_connectivity = std_conn,
        strata_names = strata_names_ref,
        mean_recruitment_rate = ens_rec,
        mean_pld_days = ens_pld,
        mean_thermal_exposure = ens_mort,
        scenarios = scenario_names,
        weights = w_vec
    )
end

function compute_ensemble_model_average(
    group::Union{Zarr.ZGroup, GeoParquetStorage};
    variable::Symbol = :mean_pld_days,
    scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
    projection_year::Int = 2050
)
    scen_strs = [string(s) for s in scenarios]
    return compute_ensemble_model_average(group, scen_strs)
end

"""
    geopublish_lakehouse!(catalog::GeoData.CatalogModule.GeoDataCatalog,
                          storage_path::AbstractString;
                          key::Symbol,
                          name::String,
                          scenario::String = "",
                          projection_year::Int = 2024,
                          derived_from::Vector{Symbol} = Symbol[])

Publish particle tracking simulation outputs (e.g. dispersal kernels, connectivity)
directly into the GeoData Lakehouse repository under the `:derived` (Gold) tier.
"""
function geopublish_lakehouse!(catalog::GeoData.CatalogModule.GeoDataCatalog,
                               storage_path::AbstractString;
                               key::Symbol,
                               name::String,
                               scenario::String = "",
                               projection_year::Int = 2024,
                               derived_from::Vector{Symbol} = Symbol[])
    return GeoData.geopublish_dataset!(
        catalog,
        storage_path;
        key = key,
        name = name,
        tier = :derived,
        producer = "ParticleTracking.jl",
        derived_from = derived_from,
        format = endswith(storage_path, ".parquet") ? :geoparquet : :zarr,
        notes = "Scenario: $(scenario), Year: $(projection_year)"
    )
end