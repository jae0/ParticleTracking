"""
    storage_duckdb.jl

DuckDB analytical storage backend for regional hydrodynamic modeling and
Lagrangian larval particle tracking (*Chionoecetes opilio*).

Provides relational persistence for simulation runs, multi-million particle
trajectory steps, cohort recruitment outcomes, demographic connectivity
matrices, and spatial dispersal fields with fast SQL analytics, multi-scenario
benchmarking, and Bayesian/ensemble model averaging.
"""

using DuckDB
using DataFrames
using DBInterface
using Dates
using Statistics
using LinearAlgebra

# Global session cache mapping canonical database paths to active DuckDB instances
const _DUCKDB_SESSION_CACHE = Dict{String, DuckDB.DB}()

const _DUCKDB_CACHE_LOCK = ReentrantLock()




"""
    open_duckdb_storage(
        db_path::AbstractString = joinpath("outputs", "particle_tracking.duckdb");
        read_only::Bool = false,
        max_retries::Int = 5,
        retry_delay::Real = 0.5
    )::DuckDB.DB

Open a connection to a DuckDB database file for particle tracking storage
and analytics. Automatically creates parent directories if needed.

# Inputs
- `db_path::AbstractString`: Path to the DuckDB database file.
- `read_only::Bool`: Whether to open the database in read-only mode.

# Outputs
- `DuckDB.DB`: Active DuckDB database connection object.

Open a connection to an embedded DuckDB database, maintaining a process-wide session
cache so multiple pipeline stages (e.g. Segments 5, 6, 7, 8) share the same underlying
database instance without triggering Windows multi-instance file lock conflicts.

# Inputs
- `db_path::AbstractString`: Path to the DuckDB database file.
- `read_only::Bool`: Whether to open the database in read-only mode.
- `max_retries::Int`: Maximum number of retries under transient file lock contention.
- `retry_delay::Real`: Seconds to wait between retries.

# Outputs
- `DuckDB.DB`: Active DuckDB database connection object.

# References
- Raasveldt, M., & Mühleisen, H. (2019). DuckDB: an embeddable analytical
  database. *Proceedings of the 2019 International Conference on Management
  of Data*, 1981-1984. DOI: 10.1145/3299869.3320212
"""
function open_duckdb_storage(
    db_path::AbstractString = joinpath("outputs", "particle_tracking.duckdb");
    read_only::Bool = false,
    max_retries::Int = 5,
    retry_delay::Real = 0.5
)::DuckDB.DB

    canon_path = abspath(db_path)
    mkpath(dirname(canon_path))

    lock(_DUCKDB_CACHE_LOCK) do
        # 1. Reuse existing active database instance within this process if available
        if haskey(_DUCKDB_SESSION_CACHE, canon_path)
            cached_db = _DUCKDB_SESSION_CACHE[canon_path]
            if isfile(canon_path)
                try
                    DBInterface.execute(cached_db, "SELECT 1;")
                    return cached_db
                catch
                end
            end
            delete!(_DUCKDB_SESSION_CACHE, canon_path)
        end

        # 2. Open new database instance with retry loop for transient locks
        db = nothing
        for attempt in 1:max_retries
            try
                db = DuckDB.DB(canon_path; readonly = read_only)
                break
            catch err
                err_msg = sprint(showerror, err)
                is_locked = occursin("used by another process", err_msg) ||
                            occursin("Cannot open file", err_msg)
                if is_locked && attempt < max_retries
                    sleep(Float64(retry_delay))
                else
                    rethrow(err)
                end
            end
        end

        if !isnothing(db)
            _DUCKDB_SESSION_CACHE[canon_path] = db
            if !read_only
                initialize_duckdb_schema!(db)
            end
        end
        return db
    end
end

"""
    close_duckdb_storage(db::DuckDB.DB; force::Bool = false)

Close the DuckDB database connection. Flushes pending write-ahead log (WAL) data
via `CHECKPOINT;`. When `force = false` (default), the instance is kept in the
process session cache so subsequent pipeline stages can access it without Windows
file lock conflicts. Pass `force = true` to evict from cache and close the handle.

# Inputs
- `db::DuckDB.DB`: Database connection to close.
- `force::Bool`: Whether to forcefully close and evict from cache.
"""
function close_duckdb_storage(db::DuckDB.DB; force::Bool = false)
    try
        DBInterface.execute(db, "CHECKPOINT;")
    catch
    end

    if !force
        return nothing
    end

    lock(_DUCKDB_CACHE_LOCK) do
        for (path, cached_db) in _DUCKDB_SESSION_CACHE
            if cached_db === db
                delete!(_DUCKDB_SESSION_CACHE, path)
                break
            end
        end
    end

    try
        DuckDB.close(db)
    catch
    end
    GC.gc()
    return nothing
end

"""
    close_all_duckdb_storage!()

Explicitly checkpoint, close, and evict all active DuckDB database instances
held in the process session cache.
"""
function close_all_duckdb_storage!()
    lock(_DUCKDB_CACHE_LOCK) do
        for (path, db) in _DUCKDB_SESSION_CACHE
            try
                DBInterface.execute(db, "CHECKPOINT;")
                DuckDB.close(db)
            catch
            end
        end
        empty!(_DUCKDB_SESSION_CACHE)
    end
    GC.gc()
    return nothing
end

"""
    initialize_duckdb_schema!(db::DuckDB.DB)

Initialize relational tables and indices for simulation metadata, particle
trajectories, demographic connectivity matrices, and recruitment metrics.

# Schema Architecture
- `simulation_runs`: Run metadata, climate scenarios, and physical parameters.
- `particle_trajectories`: High-resolution columnar spatiotemporal particle states.
- `recruitment_metrics`: Cohort summary statistics and thermal mortality totals.
- `connectivity_transitions`: Directed demographic transition probabilities across strata.
- `gridded_dispersal_summary`: Spatial 2D dispersal, density, and empirical velocity fields.

# Inputs
- `db::DuckDB.DB`: Database connection.
"""
function initialize_duckdb_schema!(db::DuckDB.DB)
    # 1. Simulation runs metadata table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS simulation_runs (
            run_id VARCHAR PRIMARY KEY,
            scenario VARCHAR NOT NULL,
            projection_year INTEGER NOT NULL,
            created_at TIMESTAMP NOT NULL,
            n_particles INTEGER NOT NULL,
            duration_seconds DOUBLE NOT NULL,
            dt_seconds DOUBLE NOT NULL,
            enable_tides BOOLEAN NOT NULL,
            enable_dvm BOOLEAN NOT NULL,
            enable_molting BOOLEAN NOT NULL,
            diffusivity_h DOUBLE NOT NULL,
            diffusivity_v DOUBLE NOT NULL,
            min_seabed_depth DOUBLE NOT NULL,
            seed INTEGER NOT NULL,
            config_toml VARCHAR,
            notes VARCHAR
        );
    """)

    # Schema migration: ensure config_toml exists in older database files
    try
        DBInterface.execute(db, "ALTER TABLE simulation_runs ADD COLUMN IF NOT EXISTS config_toml VARCHAR;")
    catch
    end

    # 2. Particle trajectories time series table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS particle_trajectories (
            run_id VARCHAR NOT NULL,
            particle_id INTEGER NOT NULL,
            step INTEGER NOT NULL,
            time_seconds DOUBLE NOT NULL,
            lon DOUBLE NOT NULL,
            lat DOUBLE NOT NULL,
            depth DOUBLE NOT NULL,
            temperature DOUBLE NOT NULL,
            degree_days DOUBLE NOT NULL,
            survival_probability DOUBLE NOT NULL,
            stage VARCHAR NOT NULL,
            alive BOOLEAN NOT NULL,
            settlement_status VARCHAR NOT NULL,
            PRIMARY KEY (run_id, particle_id, step)
        );
    """)

    # 3. Cohort recruitment metrics table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS recruitment_metrics (
            run_id VARCHAR PRIMARY KEY,
            total_released INTEGER NOT NULL,
            total_settled_successful INTEGER NOT NULL,
            total_settled_unsuitable INTEGER NOT NULL,
            total_dead_thermal INTEGER NOT NULL,
            total_pelagic_remaining INTEGER NOT NULL,
            settlement_success_rate DOUBLE NOT NULL,
            mean_pld_days DOUBLE,
            mean_degree_days DOUBLE NOT NULL,
            mean_exposure_temperature DOUBLE NOT NULL,
            mean_dispersal_distance_km DOUBLE NOT NULL
        );
    """)

    # 4. Connectivity transitions table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS connectivity_transitions (
            run_id VARCHAR NOT NULL,
            source_stratum VARCHAR NOT NULL,
            destination_stratum VARCHAR NOT NULL,
            particle_count INTEGER NOT NULL,
            transition_probability DOUBLE NOT NULL,
            PRIMARY KEY (run_id, source_stratum, destination_stratum)
        );
    """)

    # 5. Gridded dispersal spatial fields summary table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS gridded_dispersal_summary (
            run_id VARCHAR NOT NULL,
            lon DOUBLE NOT NULL,
            lat DOUBLE NOT NULL,
            empirical_u DOUBLE,
            empirical_v DOUBLE,
            empirical_diffusivity DOUBLE,
            settlement_density DOUBLE,
            mean_temp DOUBLE,
            mean_degree_days DOUBLE,
            sample_count INTEGER NOT NULL,
            PRIMARY KEY (run_id, lon, lat)
        );
    """)

    # 6. Hydrodynamic fields summary and snapshots table
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS hydrodynamic_fields (
            run_id VARCHAR NOT NULL,
            scenario VARCHAR NOT NULL,
            projection_year INTEGER NOT NULL,
            time_seconds DOUBLE NOT NULL,
            grid_x INTEGER NOT NULL,
            grid_y INTEGER NOT NULL,
            depth_level INTEGER NOT NULL,
            lon DOUBLE NOT NULL,
            lat DOUBLE NOT NULL,
            depth DOUBLE NOT NULL,
            u DOUBLE,
            v DOUBLE,
            w DOUBLE,
            temperature DOUBLE,
            salinity DOUBLE,
            elevation DOUBLE,
            PRIMARY KEY (run_id, time_seconds, grid_x, grid_y, depth_level)
        );
    """)

    # Create analytical indexes for fast queries
    DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS idx_traj_run_time ON particle_trajectories(run_id, time_seconds);")
    DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS idx_traj_stage ON particle_trajectories(stage);")
    DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS idx_conn_run ON connectivity_transitions(run_id);")
    DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS idx_hydro_run_time ON hydrodynamic_fields(run_id, time_seconds);")
    DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS idx_hydro_scen_yr ON hydrodynamic_fields(scenario, projection_year);")

    return nothing
end

"""
    save_simulation_run!(
        db::DuckDB.DB,
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
and 2D spatial dispersal fields into DuckDB within an atomic transaction.

# Inputs
- `db::DuckDB.DB`: DuckDB database connection.
- `run_id::AbstractString`: Unique identifier for the simulation run.
- `opts`: Runtime options struct or named tuple.
- `trajectories::NamedTuple`: Particle tracking output trajectory NamedTuple.
- `metrics::Union{Nothing, NamedTuple}`: Optional cohort summary metrics.
- `connectivity::Union{Nothing, NamedTuple}`: Optional connectivity matrix.
- `gridded_dispersal::Union{Nothing, NamedTuple}`: Optional gridded fields.
- `config::Union{Nothing, AbstractDict, AbstractString}`: Optional configuration dict or TOML string.
- `notes::AbstractString`: User-provided description or experiment notes.

# Outputs
- `String`: The registered `run_id`.
"""
function save_simulation_run!(
    db::DuckDB.DB,
    run_id::AbstractString,
    opts;
    trajectories::NamedTuple,
    metrics::Union{Nothing, NamedTuple} = nothing,
    connectivity::Union{Nothing, NamedTuple} = nothing,
    gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
    config::Union{Nothing, AbstractDict, AbstractString} = nothing,
    notes::AbstractString = ""
)::String

    # Canonicalize trajectories immediately to guarantee uniform 2D/1D shapes, Symbols, and complete fields
    canonical_trajs = canonicalize_trajectories(trajectories)

    # Extract metadata properties safely
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
    created_at = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    # Serialize configuration dictionary to TOML string
    config_toml_str = if !isnothing(config)
        if config isa AbstractString
            config
        else
            s_io = IOBuffer()
            TOML.print(s_io, config; sorted = true)
            String(take!(s_io))
        end
    elseif isdefined(Main, :options_to_configuration) && opts isa HydrodynamicOptions
        s_io = IOBuffer()
        TOML.print(s_io, options_to_configuration(opts); sorted = true)
        String(take!(s_io))
    else
        ""
    end

    # Delete existing records for run_id if re-saving (parameterized to avoid injection)
    for table in ("simulation_runs", "particle_trajectories",
                  "recruitment_metrics", "connectivity_transitions",
                  "gridded_dispersal_summary")
        DBInterface.execute(
            db,
            "DELETE FROM $(table) WHERE run_id = ?;",
            [String(run_id)]
        )
    end

    # 1. Insert simulation run metadata
    stmt_meta = DBInterface.prepare(db, """
        INSERT INTO simulation_runs (
            run_id, scenario, projection_year, created_at, n_particles,
            duration_seconds, dt_seconds, enable_tides, enable_dvm,
            enable_molting, diffusivity_h, diffusivity_v, min_seabed_depth,
            seed, config_toml, notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """)
    DBInterface.execute(
        stmt_meta,
        [
            run_id, scenario, proj_year, created_at, n_parts,
            duration, dt_val, tides, dvm, molt, diff_h, diff_v,
            min_depth, seed_val, config_toml_str, notes
        ]
    )

    # 2. Insert particle trajectories via DataFrame batch append for high speed
    n_p, n_t = size(canonical_trajs.lons)
    total_rows = n_p * n_t

    df_run_id = fill(String(run_id), total_rows)
    df_p_id = Vector{Int}(undef, total_rows)
    df_step = Vector{Int}(undef, total_rows)
    df_time = Vector{Float64}(undef, total_rows)
    df_lon = Vector{Float64}(undef, total_rows)
    df_lat = Vector{Float64}(undef, total_rows)
    df_depth = Vector{Float64}(undef, total_rows)
    df_temp = Vector{Float64}(undef, total_rows)
    df_dd = Vector{Float64}(undef, total_rows)
    df_surv = Vector{Float64}(undef, total_rows)
    df_stage = Vector{String}(undef, total_rows)
    df_alive = Vector{Bool}(undef, total_rows)
    df_settle = Vector{String}(undef, total_rows)

    idx = 1
    for p in 1:n_p
        p_id     = canonical_trajs.ids[p]
        p_alive  = canonical_trajs.alive[p]
        p_settle = string(canonical_trajs.settlement_status[p])

        for s in 1:n_t
            df_p_id[idx]   = p_id
            df_step[idx]   = s
            df_time[idx]   = canonical_trajs.times[s]
            df_lon[idx]    = canonical_trajs.lons[p, s]
            df_lat[idx]    = canonical_trajs.lats[p, s]
            df_depth[idx]  = canonical_trajs.depths[p, s]
            df_temp[idx]   = canonical_trajs.temperatures[p, s]
            df_dd[idx]     = canonical_trajs.degree_days_timeseries[p, s]
            df_surv[idx]   = canonical_trajs.survival_probability[p, s]
            df_stage[idx]  = string(canonical_trajs.stages[p, s])
            df_alive[idx]  = p_alive
            df_settle[idx] = p_settle
            idx += 1
        end
    end

    traj_df = DataFrame(
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

    DuckDB.register_data_frame(db, traj_df, "temp_trajectories_view")
    DBInterface.execute(db, "INSERT INTO particle_trajectories SELECT * FROM temp_trajectories_view;")
    DuckDB.unregister_data_frame(db, "temp_trajectories_view")

    # 3. Insert cohort recruitment metrics
    r_earth = 6371.0 # Earth radius in km
    displacements = [
        r_earth * acos(clamp(
            sind(canonical_trajs.lats[p, 1]) * sind(canonical_trajs.lats[p, end]) +
            cosd(canonical_trajs.lats[p, 1]) * cosd(canonical_trajs.lats[p, end]) *
            cosd(canonical_trajs.lons[p, end] - canonical_trajs.lons[p, 1]),
            -1.0, 1.0
        ))
        for p in 1:n_p
    ]
    mean_disp = safe_mean(displacements; default = 0.0)

    tot_settled_succ   = count(==( :settled_successful ), canonical_trajs.settlement_status)
    tot_settled_unsuit = count(==( :settled_unsuitable ), canonical_trajs.settlement_status)
    tot_dead           = count(==( :dead ), canonical_trajs.stages[:, end])
    tot_pelagic        = max(0, n_p - tot_settled_succ - tot_settled_unsuit - tot_dead)
    succ_rate          = n_p > 0 ? (tot_settled_succ / n_p) : 0.0

    # Safe mean PLD: guard against empty settlement set to prevent ArgumentError
    ages_valid = hasproperty(canonical_trajs, :settlement_age) ?
                 [age for age in canonical_trajs.settlement_age if age < duration] : Float64[]
    mean_pld   = isempty(ages_valid) ? (duration / 86400.0) : (safe_mean(ages_valid) / 86400.0)

    mean_dd = safe_mean(canonical_trajs.degree_days; default = 0.0)
    mean_t  = safe_mean(canonical_trajs.temperatures; default = 4.0)

    stmt_rec = DBInterface.prepare(db, """
        INSERT INTO recruitment_metrics (
            run_id, total_released, total_settled_successful, total_settled_unsuitable,
            total_dead_thermal, total_pelagic_remaining, settlement_success_rate,
            mean_pld_days, mean_degree_days, mean_exposure_temperature, mean_dispersal_distance_km
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """)
    DBInterface.execute(
        stmt_rec,
        [
            run_id, n_p, tot_settled_succ, tot_settled_unsuit,
            tot_dead, tot_pelagic, succ_rate,
            mean_pld, mean_dd, mean_t, mean_disp
        ]
    )

    # 4. Insert connectivity transitions
    if !isnothing(connectivity)
        stmt_conn = DBInterface.prepare(db, """
            INSERT INTO connectivity_transitions (
                run_id, source_stratum, destination_stratum, particle_count, transition_probability
            ) VALUES (?, ?, ?, ?, ?);
        """)
        s_names = connectivity.strata_names
        n_strata = length(s_names)
        mat = connectivity.matrix
        c_mat = if hasproperty(connectivity, :counts_unweighted)
            # Prefer raw integer counts for the particle_count DB column
            connectivity.counts_unweighted
        elseif hasproperty(connectivity, :counts_matrix)
            # Fallback: round float-weighted counts to nearest integer
            round.(Int, connectivity.counts_matrix)
        else
            round.(Int, mat .* n_p)
        end

        for i in 1:n_strata, j in 1:n_strata
            src_name = s_names[i]
            dst_name = s_names[j]
            p_count  = Int(c_mat[i, j])
            t_prob = Float64(mat[i, j])
            DBInterface.execute(stmt_conn, [run_id, src_name, dst_name, p_count, t_prob])
        end
    end

    # 5. Insert gridded dispersal fields if provided
    if !isnothing(gridded_dispersal)
        lon_c = gridded_dispersal.lon_centers
        lat_c = gridded_dispersal.lat_centers
        n_lx = length(lon_c)
        n_ly = length(lat_c)
        grid_rows = n_lx * n_ly

        g_run = fill(String(run_id), grid_rows)
        g_lon = Vector{Float64}(undef, grid_rows)
        g_lat = Vector{Float64}(undef, grid_rows)
        g_u = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_v = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_diff = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_dens = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_temp = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_dd = Vector{Union{Float64, Missing}}(undef, grid_rows)
        g_count = Vector{Int}(undef, grid_rows)

        has_u = hasproperty(gridded_dispersal, :u_mean)
        has_v = hasproperty(gridded_dispersal, :v_mean)
        has_d = hasproperty(gridded_dispersal, :diffusivity)
        has_dens = hasproperty(gridded_dispersal, :density)
        has_gtemp = hasproperty(gridded_dispersal, :mean_exposure_temperature)
        has_gdd = hasproperty(gridded_dispersal, :mean_degree_days)
        has_cnt = hasproperty(gridded_dispersal, :sample_count)

        g_idx = 1
        for i in 1:n_lx, j in 1:n_ly
            g_lon[g_idx] = Float64(lon_c[i])
            g_lat[g_idx] = Float64(lat_c[j])
            g_u[g_idx] = (has_u && !isnan(gridded_dispersal.u_mean[i, j])) ? Float64(gridded_dispersal.u_mean[i, j]) : missing
            g_v[g_idx] = (has_v && !isnan(gridded_dispersal.v_mean[i, j])) ? Float64(gridded_dispersal.v_mean[i, j]) : missing
            g_diff[g_idx] = (has_d && !isnan(gridded_dispersal.diffusivity[i, j])) ? Float64(gridded_dispersal.diffusivity[i, j]) : missing
            g_dens[g_idx] = (has_dens && !isnan(gridded_dispersal.density[i, j])) ? Float64(gridded_dispersal.density[i, j]) : missing
            g_temp[g_idx] = (has_gtemp && !isnan(gridded_dispersal.mean_exposure_temperature[i, j])) ? Float64(gridded_dispersal.mean_exposure_temperature[i, j]) : missing
            g_dd[g_idx] = (has_gdd && !isnan(gridded_dispersal.mean_degree_days[i, j])) ? Float64(gridded_dispersal.mean_degree_days[i, j]) : missing
            g_count[g_idx] = has_cnt ? Int(gridded_dispersal.sample_count[i, j]) : 0
            g_idx += 1
        end

        grid_df = DataFrame(
            run_id = g_run,
            lon = g_lon,
            lat = g_lat,
            empirical_u = g_u,
            empirical_v = g_v,
            empirical_diffusivity = g_diff,
            settlement_density = g_dens,
            mean_temp = g_temp,
            mean_degree_days = g_dd,
            sample_count = g_count
        )

        DuckDB.register_data_frame(db, grid_df, "temp_grid_view")
        DBInterface.execute(db, "INSERT INTO gridded_dispersal_summary SELECT * FROM temp_grid_view;")
        DuckDB.unregister_data_frame(db, "temp_grid_view")
    end

    return String(run_id)
end

"""
    list_simulation_runs(
        db::DuckDB.DB;
        scenario::Union{Nothing, AbstractString, Symbol} = nothing,
        projection_year::Union{Nothing, Int} = nothing
    )::DataFrame

Query and return a summary `DataFrame` of all simulation runs archived in DuckDB,
with optional filtering by climate scenario or projection year.

# Inputs
- `db::DuckDB.DB`: Database connection.
- `scenario`: Optional climate scenario filter (`:ssp245`, `:ssp585`, etc.).
- `projection_year`: Optional projection year filter (e.g. 2050).

# Outputs
- `DataFrame`: Table of archived runs joined with recruitment outcomes.
"""
function list_simulation_runs(
    db::DuckDB.DB;
    scenario::Union{Nothing, AbstractString, Symbol} = nothing,
    projection_year::Union{Nothing, Int} = nothing
)::DataFrame

    query = """
        SELECT
            r.run_id,
            r.scenario,
            r.projection_year,
            r.created_at,
            r.n_particles,
            r.duration_seconds / 86400.0 AS duration_days,
            r.enable_tides,
            r.enable_dvm,
            r.enable_molting,
            m.settlement_success_rate,
            m.mean_pld_days,
            m.mean_degree_days,
            m.mean_exposure_temperature,
            m.mean_dispersal_distance_km,
            r.notes
        FROM simulation_runs r
        LEFT JOIN recruitment_metrics m ON r.run_id = m.run_id
        WHERE 1=1
    """
    if !isnothing(scenario)
        query *= " AND r.scenario = '$(string(scenario))'"
    end
    if !isnothing(projection_year)
        query *= " AND r.projection_year = $(projection_year)"
    end
    query *= " ORDER BY r.created_at DESC;"

    return DataFrame(DBInterface.execute(db, query))
end

"""
    load_trajectories_df(
        db::DuckDB.DB,
        run_id::AbstractString;
        particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
        stage::Union{Nothing, Symbol, AbstractString} = nothing,
        time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
        max_particles::Union{Nothing, Int} = nothing
    )::DataFrame

Retrieve particle trajectory records from DuckDB as a `DataFrame` with flexible
spatial, temporal, developmental stage, and particle ID filtering.

# Inputs
- `db::DuckDB.DB`: Database connection.
- `run_id::AbstractString`: Target simulation run ID.
- `particle_ids`: Optional list of particle IDs.
- `stage`: Optional stage filter (`:zoea1`, `:megalopa`, `:instar1_settled`).
- `time_range`: Optional `(min_t, max_t)` in seconds.
- `max_particles`: Maximum number of particles to load.

# Outputs
- `DataFrame`: Columnar particle trajectory records.
"""
function load_trajectories_df(
    db::DuckDB.DB,
    run_id::AbstractString;
    particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
    stage::Union{Nothing, Symbol, AbstractString} = nothing,
    time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
    max_particles::Union{Nothing, Int} = nothing
)::DataFrame

    query = """
        SELECT
            run_id, particle_id, step, time_seconds,
            lon, lat, depth, temperature, degree_days,
            survival_probability, stage, alive, settlement_status
        FROM particle_trajectories
        WHERE run_id = ?
    """
    params = Any[String(run_id)]

    if !isnothing(particle_ids) && !isempty(particle_ids)
        placeholders = join(fill("?", length(particle_ids)), ", ")
        query *= " AND particle_id IN ($(placeholders))"
        append!(params, particle_ids)
    elseif !isnothing(max_particles) && max_particles > 0
        query *= " AND particle_id <= ?"
        push!(params, max_particles)
    end
    if !isnothing(stage)
        query *= " AND stage = ?"
        push!(params, string(stage))
    end
    if !isnothing(time_range)
        query *= " AND time_seconds BETWEEN ? AND ?"
        push!(params, Float64(time_range[1]))
        push!(params, Float64(time_range[2]))
    end
    query *= " ORDER BY particle_id, step;"

    stmt = DBInterface.prepare(db, query)
    return DataFrame(DBInterface.execute(stmt, params))
end

"""
    load_trajectories_namedtuple(
        db::DuckDB.DB,
        run_id::AbstractString;
        max_particles::Union{Nothing, Int} = nothing
    )::NamedTuple

Retrieve Lagrangian particle tracking output from DuckDB and reconstruct the
complete `NamedTuple` structure with fields `(lons, lats, depths, temperatures,
degree_days, degree_days_timeseries, survival_probability, stages, alive,
settlement_status, settlement_age, times, ids)`.

# Inputs
- `db::DuckDB.DB`: DuckDB database connection.
- `run_id::AbstractString`: Target simulation run ID.
- `max_particles::Union{Nothing, Int}`: Optional limit on the number of particles.

# Outputs
- `NamedTuple`: Trajectory cohort object identical to `track_larval_cohort` return value.
"""
function load_trajectories_namedtuple(
    db::DuckDB.DB,
    run_id::AbstractString;
    max_particles::Union{Nothing, Int} = nothing
)::NamedTuple

    query = """
        SELECT
            particle_id, step, time_seconds,
            lon, lat, depth, temperature, degree_days,
            survival_probability, stage, alive, settlement_status
        FROM particle_trajectories
        WHERE run_id = ?
    """
    params = Any[String(run_id)]
    if !isnothing(max_particles) && max_particles > 0
        query *= " AND particle_id <= ?"
        push!(params, max_particles)
    end
    query *= " ORDER BY particle_id, step;"

    stmt = DBInterface.prepare(db, query)
    df = DataFrame(DBInterface.execute(stmt, params))
    if nrow(df) == 0
        error("No trajectory records found in DuckDB for run_id: '$(run_id)'")
    end

    p_ids = sort(unique(df.particle_id))
    n_p = length(p_ids)
    time_vals = sort(unique(df.time_seconds))
    n_t = length(time_vals)

    lons = zeros(Float64, n_p, n_t)
    lats = zeros(Float64, n_p, n_t)
    depths = zeros(Float64, n_p, n_t)
    temps = zeros(Float64, n_p, n_t)
    dds_ts = zeros(Float64, n_p, n_t)
    survs = zeros(Float64, n_p, n_t)
    stages = fill(:zoea1, n_p, n_t)
    alive_vec = fill(true, n_p)
    settle_vec = fill(:pelagic, n_p)
    settle_age = fill(time_vals[end], n_p)
    degree_days = zeros(Float64, n_p)

    id_to_idx = Dict(id => idx for (idx, id) in enumerate(p_ids))

    for row in eachrow(df)
        p = id_to_idx[row.particle_id]
        s = row.step
        if 1 <= s <= n_t
            lons[p, s] = Float64(row.lon)
            lats[p, s] = Float64(row.lat)
            depths[p, s] = Float64(row.depth)
            temps[p, s] = Float64(row.temperature)
            dds_ts[p, s] = Float64(row.degree_days)
            survs[p, s] = Float64(row.survival_probability)
            stages[p, s] = Symbol(row.stage)
        end
        if s == n_t
            alive_vec[p] = Bool(row.alive)
            settle_vec[p] = Symbol(row.settlement_status)
            degree_days[p] = Float64(row.degree_days)
        end
    end

    return (
        lons = lons,
        lats = lats,
        depths = depths,
        temperatures = temps,
        degree_days = degree_days,
        degree_days_timeseries = dds_ts,
        survival_probability = survs,
        stages = stages,
        alive = alive_vec,
        settlement_status = settle_vec,
        settlement_age = settle_age,
        times = time_vals,
        ids = p_ids
    )
end

"""
    load_all_scenario_trajectories(
        db::DuckDB.DB;
        projection_year::Union{Nothing, Int} = nothing
    )::Dict{Symbol, NamedTuple}

Load trajectories for all distinct climate scenarios present in DuckDB.

# Inputs
- `db::DuckDB.DB`: DuckDB database connection.
- `projection_year::Union{Nothing, Int}`: Optional projection year filter.

# Outputs
- `Dict{Symbol, NamedTuple}`: Mapping scenario names to trajectory NamedTuples.
"""
function load_all_scenario_trajectories(
    db::DuckDB.DB;
    projection_year::Union{Nothing, Int} = nothing
)::Dict{Symbol, NamedTuple}

    runs_df = list_simulation_runs(db; projection_year = projection_year)
    if nrow(runs_df) == 0
        return Dict{Symbol, NamedTuple}()
    end

    scenario_dict = Dict{Symbol, NamedTuple}()
    for scen_name in unique(runs_df.scenario)
        subset_df = filter(r -> r.scenario == scen_name, runs_df)
        latest_run_id = string(subset_df.run_id[1])
        try
            scenario_dict[Symbol(scen_name)] = load_trajectories_namedtuple(db, latest_run_id)
        catch err
            @warn "Failed to load trajectories for scenario $(scen_name): $(err)"
        end
    end

    return scenario_dict
end

"""
    load_connectivity_matrix(
        db::DuckDB.DB,
        run_id::AbstractString
    )::NamedTuple

Retrieve the transition probability connectivity matrix \$P_{ij}\$ and stratum
names for a designated run from DuckDB.

# Inputs
- `db::DuckDB.DB`: Database connection.
- `run_id::AbstractString`: Run identifier.

# Outputs
- `NamedTuple`: `(matrix::Matrix{Float64}, counts_matrix::Matrix{Int}, strata_names::Vector{String})`
"""
function load_connectivity_matrix(
    db::DuckDB.DB,
    run_id::AbstractString
)::NamedTuple

    query = """
        SELECT source_stratum, destination_stratum, particle_count, transition_probability
        FROM connectivity_transitions
        WHERE run_id = '$(run_id)'
        ORDER BY source_stratum, destination_stratum;
    """
    df = DataFrame(DBInterface.execute(db, query))
    if nrow(df) == 0
        error("No connectivity transitions found for run_id: $(run_id)")
    end

    strata_names = sort(unique(df.source_stratum))
    n_strata = length(strata_names)
    prob_mat = zeros(Float64, n_strata, n_strata)
    count_mat = zeros(Int, n_strata, n_strata)

    name_to_idx = Dict(name => idx for (idx, name) in enumerate(strata_names))

    for row in eachrow(df)
        i = get(name_to_idx, row.source_stratum, 0)
        j = get(name_to_idx, row.destination_stratum, 0)
        if i > 0 && j > 0
            prob_mat[i, j] = Float64(row.transition_probability)
            count_mat[i, j] = Int(row.particle_count)
        end
    end

    return (
        matrix = prob_mat,
        counts_matrix = count_mat,
        strata_names = strata_names
    )
end

"""
    compare_scenarios(
        db::DuckDB.DB;
        scenario_names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
        projection_years::Union{Nothing, AbstractVector{Int}} = nothing
    )::DataFrame

Perform comparative aggregation across archived simulation runs, evaluating
changes in recruitment success, pelagic larval duration (PLD), thermal mortality,
and dispersal displacement under climate scenarios.

# Mathematical Formulation
Computes scenario-level means and standard deviations:
```math
\\bar{Y}_s = \\frac{1}{K_s} \\sum_{k=1}^{K_s} Y_{s, k}, \\quad
\\sigma(Y_s) = \\sqrt{\\frac{1}{K_s - 1} \\sum_{k=1}^{K_s} (Y_{s, k} - \\bar{Y}_s)^2}
```

# Inputs
- `db::DuckDB.DB`: Database connection.
- `scenario_names`: Optional list of scenario names to compare.
- `projection_years`: Optional list of projection years.

# Outputs
- `DataFrame`: Aggregated comparative table of scenario outcomes.
"""
function compare_scenarios(
    db::DuckDB.DB;
    scenario_names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    projection_years::Union{Nothing, AbstractVector{Int}} = nothing
)::DataFrame

    query = """
        SELECT
            r.scenario,
            r.projection_year,
            COUNT(r.run_id) AS n_runs,
            AVG(m.settlement_success_rate) AS mean_settlement_success,
            STDDEV_SAMP(m.settlement_success_rate) AS std_settlement_success,
            AVG(m.mean_pld_days) AS mean_pld_days,
            STDDEV_SAMP(m.mean_pld_days) AS std_pld_days,
            AVG(m.mean_degree_days) AS mean_degree_days,
            AVG(m.mean_exposure_temperature) AS mean_temperature_celsius,
            AVG(m.mean_dispersal_distance_km) AS mean_dispersal_km,
            AVG(m.total_dead_thermal * 1.0 / r.n_particles) AS mean_thermal_mortality_rate
        FROM simulation_runs r
        JOIN recruitment_metrics m ON r.run_id = m.run_id
        WHERE 1=1
    """
    if !isnothing(scenario_names) && !isempty(scenario_names)
        scen_list = join(["'$(s)'" for s in scenario_names], ", ")
        query *= " AND r.scenario IN ($(scen_list))"
    end
    if !isnothing(projection_years) && !isempty(projection_years)
        yr_list = join(projection_years, ", ")
        query *= " AND r.projection_year IN ($(yr_list))"
    end
    query *= """
        GROUP BY r.scenario, r.projection_year
        ORDER BY r.scenario, r.projection_year;
    """

    return DataFrame(DBInterface.execute(db, query))
end

"""
    compute_ensemble_model_average(
        db::DuckDB.DB,
        scenario_names::AbstractVector{<:AbstractString};
        weights::Union{Nothing, AbstractVector{<:Real}} = nothing
    )::NamedTuple

Compute weighted ensemble model-averaged connectivity matrices and demographic
metrics across multiple climate projections or multi-realization ensembles.

# Mathematical Formulation
For \$M\$ models/scenarios with weights \$w_m\$ (\$\\sum_m w_m = 1\$):
```math
\\bar{P}_{ij} = \\sum_{m=1}^M w_m P_{ij}^{(m)}
```
```math
\\sigma^2(P_{ij}) = \\sum_{m=1}^M w_m \\left( P_{ij}^{(m)} - \\bar{P}_{ij} \\right)^2
```

# Inputs
- `db::DuckDB.DB`: Database connection.
- `scenario_names`: List of scenarios to average (e.g. `["ssp126", "ssp245", "ssp585"]`).
- `weights`: Optional weighting vector \$w_m\$. If `nothing`, equal weighting is applied.

# Outputs
- `NamedTuple`:
  - `mean_connectivity`: Ensemble mean transition probability matrix \$\\bar{P}_{ij}\$.
  - `std_connectivity`: Ensemble standard deviation matrix \$\\sigma(P_{ij})\$.
  - `strata_names`: Names of spatial management strata.
  - `mean_recruitment_rate`: Weighted mean settlement success rate.
  - `mean_pld_days`: Weighted mean pelagic larval duration.
  - `mean_thermal_exposure`: Weighted mean thermal exposure temperature.
  - `scenarios`: Included scenarios.
  - `weights`: Normalized weights applied.

# References
- Hoeting, J. A., et al. (1999). Bayesian model averaging: a tutorial.
  *Statistical Science*, 14(4), 382-401. DOI: 10.1214/ss/1009212519
"""
function compute_ensemble_model_average(
    db::DuckDB.DB,
    scenario_names::AbstractVector{<:AbstractString};
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing
)::NamedTuple

    n_scen = length(scenario_names)
    if n_scen == 0
        error("At least one scenario name must be provided for ensemble averaging.")
    end

    # Normalize weights
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

    # Gather connectivity matrices and metrics for each scenario
    matrices = Matrix{Float64}[]
    rec_rates = Float64[]
    pld_vals = Float64[]
    mort_vals = Float64[]
    strata_names_ref = String[]

    for (s_idx, scen) in enumerate(scenario_names)
        runs_df = list_simulation_runs(db, scenario = scen)
        if nrow(runs_df) == 0
            error("No runs found in DuckDB for scenario: '$(scen)'.")
        end
        # Use latest run for each scenario
        latest_run_id = string(runs_df.run_id[1])
        conn = load_connectivity_matrix(db, latest_run_id)

        if s_idx == 1
            strata_names_ref = conn.strata_names
        end

        push!(matrices, conn.matrix)
        push!(rec_rates, Float64(runs_df.settlement_success_rate[1]))
        push!(pld_vals, Float64(runs_df.mean_pld_days[1]))
        push!(mort_vals, Float64(runs_df.mean_exposure_temperature[1]))
    end

    n_strata = length(strata_names_ref)
    mean_conn = zeros(Float64, n_strata, n_strata)
    var_conn = zeros(Float64, n_strata, n_strata)

    # Compute weighted mean
    for m in 1:n_scen
        mean_conn .+= w_vec[m] .* matrices[m]
    end

    # Compute weighted variance
    for m in 1:n_scen
        var_conn .+= w_vec[m] .* ((matrices[m] .- mean_conn) .^ 2)
    end
    std_conn = sqrt.(max.(0.0, var_conn))

    # Compute weighted ensemble scalar metrics
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

"""
    save_hydrodynamic_field!(
        db::DuckDB.DB,
        run_id::AbstractString,
        opts;
        grid_lons::AbstractVector{<:Real},
        grid_lats::AbstractVector{<:Real},
        grid_depths::AbstractVector{<:Real},
        u::AbstractArray{<:Real, 3},
        v::AbstractArray{<:Real, 3},
        w::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
        temperature::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
        salinity::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
        elevation::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        time_seconds::Real = 0.0
    )::String

Archive a 3D Eulerian hydrodynamic velocity and tracer snapshot into DuckDB.

# Inputs
- `db::DuckDB.DB`: DuckDB database connection.
- `run_id::AbstractString`: Simulation run identifier.
- `opts`: Runtime options struct or named tuple.
- `grid_lons`: 1D array of cell center longitudes.
- `grid_lats`: 1D array of cell center latitudes.
- `grid_depths`: 1D array of cell center depths.
- `u`: 3D zonal velocity array u(x, y, z) (m/s).
- `v`: 3D meridional velocity array v(x, y, z) (m/s).
- `w`: Optional 3D vertical velocity array w(x, y, z) (m/s).
- `temperature`: Optional 3D potential temperature array (°C).
- `salinity`: Optional 3D practical salinity array (PSU).
- `elevation`: Optional 2D seafloor bathymetric elevation (m).
- `time_seconds`: Timestamp of snapshot in simulation seconds.

# Outputs
- `String`: Registered `run_id`.
"""
function save_hydrodynamic_field!(
    db::DuckDB.DB,
    run_id::AbstractString,
    opts;
    grid_lons::AbstractVector{<:Real},
    grid_lats::AbstractVector{<:Real},
    grid_depths::AbstractVector{<:Real},
    u::AbstractArray{<:Real, 3},
    v::AbstractArray{<:Real, 3},
    w::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
    temperature::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
    salinity::Union{Nothing, AbstractArray{<:Real, 3}} = nothing,
    elevation::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    time_seconds::Real = 0.0
)::String

    scenario = string(hasproperty(opts, :scenario) ? opts.scenario : :baseline)
    proj_year = Int(hasproperty(opts, :projection_year) ? opts.projection_year : 2050)

    nx = length(grid_lons)
    ny = length(grid_lats)
    nz = length(grid_depths)
    total_cells = nx * ny * nz

    df_run = fill(String(run_id), total_cells)
    df_scen = fill(scenario, total_cells)
    df_yr = fill(proj_year, total_cells)
    df_time = fill(Float64(time_seconds), total_cells)
    df_gx = Vector{Int}(undef, total_cells)
    df_gy = Vector{Int}(undef, total_cells)
    df_gz = Vector{Int}(undef, total_cells)
    df_lon = Vector{Float64}(undef, total_cells)
    df_lat = Vector{Float64}(undef, total_cells)
    df_depth = Vector{Float64}(undef, total_cells)
    df_u = Vector{Union{Float64, Missing}}(undef, total_cells)
    df_v = Vector{Union{Float64, Missing}}(undef, total_cells)
    df_w = Vector{Union{Float64, Missing}}(undef, total_cells)
    df_temp = Vector{Union{Float64, Missing}}(undef, total_cells)
    df_sal = Vector{Union{Float64, Missing}}(undef, total_cells)
    df_elev = Vector{Union{Float64, Missing}}(undef, total_cells)

    idx = 1
    for i in 1:nx, j in 1:ny, k in 1:nz
        df_gx[idx] = i
        df_gy[idx] = j
        df_gz[idx] = k
        df_lon[idx] = Float64(grid_lons[i])
        df_lat[idx] = Float64(grid_lats[j])
        df_depth[idx] = Float64(grid_depths[k])
        df_u[idx] = !isnan(u[i, j, k]) ? Float64(u[i, j, k]) : missing
        df_v[idx] = !isnan(v[i, j, k]) ? Float64(v[i, j, k]) : missing
        df_w[idx] = (!isnothing(w) && !isnan(w[i, j, k])) ? Float64(w[i, j, k]) : missing
        df_temp[idx] = (!isnothing(temperature) && !isnan(temperature[i, j, k])) ? Float64(temperature[i, j, k]) : missing
        df_sal[idx] = (!isnothing(salinity) && !isnan(salinity[i, j, k])) ? Float64(salinity[i, j, k]) : missing
        df_elev[idx] = (!isnothing(elevation) && !isnan(elevation[i, j])) ? Float64(elevation[i, j]) : missing
        idx += 1
    end

    df_hydro = DataFrame(
        run_id = df_run,
        scenario = df_scen,
        projection_year = df_yr,
        time_seconds = df_time,
        grid_x = df_gx,
        grid_y = df_gy,
        depth_level = df_gz,
        lon = df_lon,
        lat = df_lat,
        depth = df_depth,
        u = df_u,
        v = df_v,
        w = df_w,
        temperature = df_temp,
        salinity = df_sal,
        elevation = df_elev
    )

    DBInterface.execute(db, "DELETE FROM hydrodynamic_fields WHERE run_id = '$(run_id)' AND time_seconds = $(time_seconds);")
    DuckDB.register_data_frame(db, df_hydro, "temp_hydro_view")
    DBInterface.execute(db, "INSERT INTO hydrodynamic_fields SELECT * FROM temp_hydro_view;")
    DuckDB.unregister_data_frame(db, "temp_hydro_view")

    return String(run_id)
end

"""
    load_hydrodynamic_field(
        db::DuckDB.DB,
        run_id::AbstractString;
        time_seconds::Union{Nothing, Real} = nothing,
        depth_level::Union{Nothing, Int} = nothing
    )::NamedTuple

Retrieve 3D/2D hydrodynamic velocity and tracer field from DuckDB as structured arrays.

# Inputs
- `db::DuckDB.DB`: DuckDB database connection.
- `run_id::AbstractString`: Target simulation run ID.
- `time_seconds`: Optional timestamp in seconds.
- `depth_level`: Optional vertical level index.

# Outputs
- `NamedTuple`: `(lons, lats, depths, u, v, w, temperature, salinity, elevation)`.
"""
function load_hydrodynamic_field(
    db::DuckDB.DB,
    run_id::AbstractString;
    time_seconds::Union{Nothing, Real} = nothing,
    depth_level::Union{Nothing, Int} = nothing
)::NamedTuple

    query = """
        SELECT
            grid_x, grid_y, depth_level, lon, lat, depth,
            u, v, w, temperature, salinity, elevation
        FROM hydrodynamic_fields
        WHERE run_id = '$(run_id)'
    """
    if !isnothing(time_seconds)
        query *= " AND time_seconds = $(time_seconds)"
    end
    if !isnothing(depth_level)
        query *= " AND depth_level = $(depth_level)"
    end
    query *= " ORDER BY grid_x, grid_y, depth_level;"

    df = DataFrame(DBInterface.execute(db, query))
    if nrow(df) == 0
        error("No hydrodynamic field records found in DuckDB for run_id: '$(run_id)'")
    end

    xs = sort(unique(df.grid_x))
    ys = sort(unique(df.grid_y))
    zs = sort(unique(df.depth_level))
    nx, ny, nz = length(xs), length(ys), length(zs)

    lons = zeros(Float64, nx)
    lats = zeros(Float64, ny)
    depths = zeros(Float64, nz)
    u_mat = zeros(Float64, nx, ny, nz)
    v_mat = zeros(Float64, nx, ny, nz)
    w_mat = zeros(Float64, nx, ny, nz)
    t_mat = zeros(Float64, nx, ny, nz)
    s_mat = zeros(Float64, nx, ny, nz)
    elev_mat = zeros(Float64, nx, ny)

    for row in eachrow(df)
        i, j, k = row.grid_x, row.grid_y, row.depth_level
        if 1 <= i <= nx && 1 <= j <= ny && 1 <= k <= nz
            lons[i] = Float64(row.lon)
            lats[j] = Float64(row.lat)
            depths[k] = Float64(row.depth)
            u_mat[i, j, k] = ismissing(row.u) ? 0.0 : Float64(row.u)
            v_mat[i, j, k] = ismissing(row.v) ? 0.0 : Float64(row.v)
            w_mat[i, j, k] = ismissing(row.w) ? 0.0 : Float64(row.w)
            t_mat[i, j, k] = ismissing(row.temperature) ? 4.0 : Float64(row.temperature)
            s_mat[i, j, k] = ismissing(row.salinity) ? 32.5 : Float64(row.salinity)
            elev_mat[i, j] = ismissing(row.elevation) ? 0.0 : Float64(row.elevation)
        end
    end

    return (
        lons = lons,
        lats = lats,
        depths = depths,
        u = u_mat,
        v = v_mat,
        w = w_mat,
        temperature = t_mat,
        salinity = s_mat,
        elevation = elev_mat
    )
end

"""
    load_gridded_dispersal(
        db::DuckDB.DB,
        run_id::AbstractString
    )::NamedTuple

Retrieve 2D gridded dispersal, settlement density, empirical advection, and
turbulent diffusivity fields for a run from DuckDB.

# Inputs
- `db::DuckDB.DB`: DuckDB connection.
- `run_id::AbstractString`: Run identifier.

# Outputs
- `NamedTuple`: `(lon_centers, lat_centers, u_mean, v_mean, diffusivity,
  settlement_density, mean_exposure_temperature, mean_degree_days, sample_count)`.
"""
function load_gridded_dispersal(
    db::DuckDB.DB,
    run_id::AbstractString
)::NamedTuple

    query = """
        SELECT
            lon, lat, empirical_u, empirical_v, empirical_diffusivity,
            settlement_density, mean_temp, mean_degree_days, sample_count
        FROM gridded_dispersal_summary
        WHERE run_id = '$(run_id)'
        ORDER BY lon, lat;
    """
    df = DataFrame(DBInterface.execute(db, query))
    if nrow(df) == 0
        error("No gridded dispersal summary found in DuckDB for run_id: '$(run_id)'")
    end

    lons = sort(unique(df.lon))
    lats = sort(unique(df.lat))
    nx = length(lons)
    ny = length(lats)

    u_grid = fill(NaN, nx, ny)
    v_grid = fill(NaN, nx, ny)
    diff_grid = fill(NaN, nx, ny)
    dens_grid = fill(NaN, nx, ny)
    temp_grid = fill(NaN, nx, ny)
    dd_grid = fill(NaN, nx, ny)
    count_grid = zeros(Int, nx, ny)

    lon_to_i = Dict(lon => i for (i, lon) in enumerate(lons))
    lat_to_j = Dict(lat => j for (j, lat) in enumerate(lats))

    for row in eachrow(df)
        i = get(lon_to_i, row.lon, 0)
        j = get(lat_to_j, row.lat, 0)
        if i > 0 && j > 0
            u_grid[i, j] = ismissing(row.empirical_u) ? NaN : Float64(row.empirical_u)
            v_grid[i, j] = ismissing(row.empirical_v) ? NaN : Float64(row.empirical_v)
            diff_grid[i, j] = ismissing(row.empirical_diffusivity) ? NaN : Float64(row.empirical_diffusivity)
            dens_grid[i, j] = ismissing(row.settlement_density) ? NaN : Float64(row.settlement_density)
            temp_grid[i, j] = ismissing(row.mean_temp) ? NaN : Float64(row.mean_temp)
            dd_grid[i, j] = ismissing(row.mean_degree_days) ? NaN : Float64(row.mean_degree_days)
            count_grid[i, j] = Int(row.sample_count)
        end
    end

    return (
        lon_centers = lons,
        lat_centers = lats,
        u_mean = u_grid,
        v_mean = v_grid,
        diffusivity = diff_grid,
        settlement_density = dens_grid,
        mean_exposure_temperature = temp_grid,
        mean_degree_days = dd_grid,
        sample_count = count_grid
    )
end

"""
    export_duckdb_to_parquet(
        db::DuckDB.DB,
        output_dir::AbstractString = joinpath("outputs", "parquet")
    )::Vector{String}

Export all relational tables from DuckDB into high-performance Apache Parquet
files for seamless interoperability with Python, R, and BSTM modeling frameworks.

# Inputs
- `db::DuckDB.DB`: Active DuckDB connection.
- `output_dir::AbstractString`: Destination directory for Parquet files.

# Outputs
- `Vector{String}`: File paths to the generated Parquet files.
"""
function export_duckdb_to_parquet(
    db::DuckDB.DB,
    output_dir::AbstractString = joinpath("outputs", "parquet")
)::Vector{String}

    mkpath(output_dir)
    tables = [
        "simulation_runs",
        "particle_trajectories",
        "recruitment_metrics",
        "connectivity_transitions",
        "gridded_dispersal_summary",
        "hydrodynamic_fields"
    ]
    exported_paths = String[]

    for tbl in tables
        out_file = joinpath(output_dir, "$(tbl).parquet")
        # DuckDB native COPY TO parquet
        try
            DBInterface.execute(db, "COPY $(tbl) TO '$(replace(out_file, "\\" => "/"))' (FORMAT PARQUET, COMPRESSION ZSTD);")
            push!(exported_paths, out_file)
        catch err
            @warn "Failed to export table $(tbl) to Parquet: $(err)"
        end
    end

    return exported_paths
end

"""
    load_run_configuration(db::DuckDB.DB, run_id::AbstractString)::Dict{String, Any}

Extract and parse the archived TOML configuration dictionary associated with a simulation run.

# Inputs
- `db::DuckDB.DB`: DuckDB connection.
- `run_id::AbstractString`: Simulation run identifier.

# Outputs
- `Dict{String, Any}`: Parsed configuration dictionary.
"""
function load_run_configuration(db::DuckDB.DB, run_id::AbstractString)::Dict{String, Any}

    stmt = DBInterface.prepare(db, "SELECT config_toml FROM simulation_runs WHERE run_id = ?;")
    cursor = DBInterface.execute(stmt, [run_id])
    df = DataFrame(cursor)
    if nrow(df) == 0 || ismissing(df.config_toml[1]) || isempty(df.config_toml[1])
        return Dict{String, Any}()
    end
    try
        return TOML.parse(df.config_toml[1])
    catch err
        @warn "Failed to parse archived configuration for run $(run_id): $(err)"
        return Dict{String, Any}()
    end

end
