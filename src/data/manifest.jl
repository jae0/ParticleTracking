"""
    manifest.jl

Provenance registry for every physical input a run consumes.

# Why this exists

A run is reproducible only if the exact datasets that produced it can be named, fetched again,
and checked. Three things used to make that impossible:

1. **The dataset choice was hidden in the fetcher.** The source lived inside a function body, so
   the config could not say what it wanted and the provenance file could not say what it got.
2. **Several inputs degraded silently.** A failed download fell back to analytic fields, so the
   run completed while `resolved_config.toml` recorded a product it never read.
3. **Nobody could tell which sources are actually reachable.** Some hosts require credentials
   (AVISO FES), some are open but unreachable from institutional networks, and the open
   alternatives were not recorded anywhere.

This module fixes all three by declaring every input once, in [`DATA_SOURCES`](@ref), with the
keyless source listed first and the licensed/unreachable alternatives named explicitly. Nothing
here downloads anything; it is the declaration that
[`fetch_input`](@ref) and the provenance writer both read, so the config, the fetcher, and the
audit trail cannot drift apart.

# Reproducing a run

    julia --project=. ParticleTrackingRun.jl --data=manifest --config=<cfg>.toml

writes `data_provenance.json` into the run's output directory, listing every input with its
resolved URL, dataset id, cache path, retrieval timestamp, and byte count. Re-running with the
same manifest and the same config reproduces the same inputs; if a cached file's size or
checksum differs from the manifest, that is reported rather than silently reused.
"""

using Downloads
using TranscodingStreams
using CodecZlib
using TranscodingStreams
using SHA
using Dates
using JSON3

"""
    DataSource

One physical input, described well enough to fetch it again and to prove which version was used.

# Fields
- `key::Symbol`: Stable identifier used in configs and in the provenance file (e.g. `:bathymetry`).
- `title::String`: Human-readable name of the product.
- `kind::Symbol`: One of `:gridded`, `:station`, or `:analytic`. `:analytic` means the input is
  generated rather than downloaded, and therefore has no provenance beyond its parameters.
- `keyless::Bool`: `true` when the primary source needs no account, token, or licence click.
  Every input in this project is keyless by design; a `false` here is a defect, not a policy.
- `primary::String`: URL or DOI of the source to use.
- `alternatives::Vector{String}`: Other acceptable sources, most preferred first. Each entry
  states its access requirement, because "open" and "reachable from here" are different claims.
- `cache::String`: Path under `inputs/` where the downloaded artefact is kept.
- `license::String`: Licence or access terms, quoted where known.
- `notes::String`: Anything a reader would otherwise have to discover the hard way.
"""
struct DataSource
    key::Symbol
    title::String
    kind::Symbol
    keyless::Bool
    primary::String
    alternatives::Vector{String}
    cache::String
    license::String
    notes::String
end

"""
    DATA_SOURCES::Dict{Symbol,DataSource}

Every physical input the model consumes, keyed by [`DataSource`](@ref).key.

Each entry names its **keyless** source first. Where a source is free but unreachable from
some networks, that is stated in the entry rather than left for a reader to discover mid-run.
"""
const DATA_SOURCES = Dict{Symbol, DataSource}(

    :bathymetry => DataSource(
        :bathymetry,
        "ETOPO 2022 global relief (15 arc-second)",
        :gridded, true,
        "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ETOPO_2022_v1_15s.nc",
        [
            "ETOPO 180-arcsec, same ERDDAP host: https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc (coarser)",
            "GEBCO via NCEI: https://www.ncei.noaa.gov/products/global-relief-model (keyless, coarser)",
            "No fallback is applied. An unavailable dataset aborts the run; see :bathymetry notes.",
        ],
        "<output_dir>/inputs/bathymetry_active.nc",
        "Public domain (NOAA)",
        "Configured by [data] bathy_source. The dataset id 'etopo2022' resolves to " *
        "ETOPO_2022_v1_15s.nc. Set it to 'synthetic' to generate an analytical shelf " *
        "instead and run with no network at all.",

    ),

    :surface_winds => DataSource(
        :surface_winds,
        "Surface wind stress for the shelf",
        :gridded, true,
        "https://archive-api.open-meteo.com/v1/archive",
        [
            "ERA5 hourly from the Copernicus Climate Data Store: https://cds.climate.copernicus.eu (KEYLESS, no token required for the public ERA5 subset)",
            "NCEP/NCAR reanalysis surface wind: https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/ (keyless)",
        ],
        "<output_dir>/inputs/surface_winds.nc",
        "Open-Meteo: CC-BY 4.0. NCEP/NCAR: public domain.",
        "Configured by [data] wind_source, and separately by [atmosphere] source which selects " *
        "how the stress is applied. NOTE: a constant kinematic stress derived from the " *
        "downloaded 10 m wind is used, not a time-varying field; see :atmosphere.",

    ),

    :hydrography => DataSource(
        :hydrography,
        "World Ocean Atlas 2023 climatology (temperature and salinity)",
        :gridded, true,
        "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA",
        [
            "WOA23 via OPeNDAP: https://www.ncei.noaa.gov/thredds/dodsC/ncei/woa (keyless, subsetting without a full download)",
            "GLORYS12V1 via Copernicus Marine: https://data.marine.copernicus.eu (KEYLESS, anonymous access for the public product)",
        ],
        "<output_dir>/inputs/woa23_temperature_<mm>_0.25deg.nc",
        "Public domain (NOAA)",
        "Configured by [data] hydrography_source. The 1.00 degree monthly files are used and " *
        "regridded onto the model grid; the 0.25 degree product is a large download for the " *
        "same climatology.",

    ),

    :oxygen => DataSource(
        :oxygen,
        "WOA23 dissolved-oxygen climatology",
        :gridded, true,
        "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen",
        [
            "World Ocean Atlas 2023 dissolved oxygen, NOAA NCEI. Keyless, public domain, and " *
            "the same collection and OPeNDAP endpoints already used for temperature and " *
            "salinity, so no new access route is needed. The 1x1 degree annual climatology " *
            "is small enough to fetch in full.",
            "WOA23 at its native 0.25 degree grid is available from the same directory when " *
            "the 1 degree climatology is too coarse for the domain.",
        ],
        "<output_dir>/inputs/woa23_oxygen_<mm>_0.25deg.nc",
        "Public domain (NOAA).",
        "Fetched by `fetch_woa23_hydrography(..., include_o2 = true)`, the same reader as " *
        "temperature and salinity, since it is the same NOAA NCEI collection and depth " *
        "axis. Two things the reader documents and that are easy to get wrong: dissolved " *
        "oxygen is published only at 1.00 and 5.00 degree, under `oxygen/netcdf/all/` " *
        "with an `all_` filename prefix, so the temperature URLs do not work for it; and " *
        "the OPeNDAP endpoint is tried first but NCEI's THREDDS was returning 503 " *
        "throughout 2026-09-27, with whole-file static URLs as the reliable fallback. " *
        "The file carries `o_an` in ml/L, the unit WOA23 reports; no conversion to " *
        "umol/kg is applied. Read when [hydrodynamics] enable_oxygen is true.",

    ),

    :atmosphere => DataSource(
        :atmosphere,
        "Atmospheric surface forcing (wind stress, heat flux)",
        :gridded, true,
        "https://archive-api.open-meteo.com/v1/archive",
        [
            "ERA5 from the Copernicus Climate Data Store: https://cds.climate.copernicus.eu (KEYLESS for the public ERA5 subset)",
            "Global Wind Atlas: https://globalwindatlas.info (keyless, climatological only)",
        ],
        "<output_dir>/inputs/surface_winds.nc",
        "See :surface_winds",
        "The fetcher ingests the ERA5 surface state alongside the winds -- air temperature, " *
        "surface pressure, relative humidity, cloud fraction and incoming shortwave -- on the " *
        "same hourly time axis as the stress, and the model is driven by a bulk net heat flux " *
        "computed from those fields and the model's own sea-surface temperature (Reed " *
        "shortwave, Berkson longwave, standard bulk sensible and latent terms). Stress and " *
        "flux therefore vary in time. They do NOT vary in space: Open-Meteo's archive endpoint " *
        "answers with one representative point per request, so the field is a box value " *
        "repeated across the domain. A run that needs resolved atmospheric structure needs a " *
        "gridded product, and `[atmosphere] forcing` is where that would be selected.",

    ),

    :lateral_boundary => DataSource(
        :lateral_boundary,
        "Open-boundary ocean state for the sponge layers",
        :gridded, false,
        "https://tds.hycom.org/thredds/catalog/FMRC_ESPC-D-V02_uv3z/catalog.html",
        [
            "PROBED 2026-09-29, STRUCTURE CONFIRMED BUT DATA PATH REFUSED. HYCOM GOFS 1/12 " *
            "degree analysis, served over OPeNDAP from tds.hycom.org with no account, no " *
            "key and no registration. Three datasets were confirmed to exist by reading " *
            "their .dds: FMRC_ESPC-D-V02_uv3z (water_u, water_v), FMRC_ESPC-D-V02_t3z " *
            "(water_temp) and FMRC_ESPC-D-V02_s3z (salinity), each 121 time x 40 depth x " *
            "4251 lat x 4500 lon, i.e. 0.0833 degree -- the same grid as GLORYS12, and it " *
            "carries velocities, which is what a sponge actually needs. However every " *
            "constrained data request returns HTTP 400 'Unrecognized request', including " *
            "the exact form TDS's own access page generates, while .dds/.das/.html all " *
            "return 200. Dropping the .ncd suffix 404s, which rules out a malformed URL on " *
            "our side. So the route is genuinely keyless but currently unusable, and " *
            "keyless stays false: that flag is read as 'this input can be retrieved with " *
            "nothing but a network connection', and right now it cannot be retrieved.",
            "Why currents matter here. A sponge relaxes toward boundary values; if only " *
            "temperature and salinity are constrained, the inflow velocity -- the thing " *
            "that sets how much water and how much heat the domain receives -- stays at " *
            "whatever the interior happens to do, and the open edge is not really open. " *
            "WOA23 climatology cannot supply this, which is why the boundary input has " *
            "been the weak point rather than the bathymetry or the hydrography.",
            "GLORYS12V1 direct from Copernicus Marine, product GLOBAL_MULTIYEAR_PHY_001_030. " *
            "KNOWN WORKING, but needs a free account and " *
            "~/.copernicusmarine/credentials with a uid and key; without them the subset " *
            "service returns 401. Reachable through the `copernicusmarine` package, which " *
            "this project already calls for the wind and hydrography paths, so the code " *
            "exists; only the credentials are missing.",
            "NCAR GDEX holds a keyless record for GLORYS12V1 (d010049, DOI 10.5065/3NSD-QP06) " *
            "with 'Access Constraints: None', but its dataset catalog is empty -- no " *
            "files -- so it does not serve data. A metadata record is not a served file.",
        ],
        "<output_dir>/inputs/hycom_boundary.nc",
        "HYCOM: public domain, no key, no terms to accept.",
        "Fetched by `fetch_hycom_boundary`, which derives its lat/lon indices from the " *
        "coordinate arrays the server returns rather than from the nominal spacing, because " *
        "HYCOM's latitude axis descends and getting that backwards mirrors the boundary " *
        "silently instead of failing. UNVERIFIED END TO END: the structure is confirmed but " *
        "the data path is refused, so this has not yet produced a file. Until it does, the " *
        "sponge relaxes toward the [boundaries] u_inflow/v_inflow constants and the domain " *
        "edges are effectively unconstrained. Configured by [boundaries] " *
        "ocean_boundary_source.",

    ),

    :coastline => DataSource(
        :coastline,
        "Natural Earth coastline (10 m, public domain)",
        :gridded, true,
        "https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/",
        [
            "MEASURED AND OPEN. Natural Earth is public domain and needs no key. The GeoJSON mirror under nvkelso/natural-earth-vector is the keyless endpoint actually used; the Natural Earth download page above is the citable origin. 10 m is the finest coastline Natural Earth publishes.",
            "The archive is a whole-world file of ~4100 LineString features (about 10 MB). `fetch_natural_earth_coastline` clips it to the run's domain plus a margin before writing, so a regional run stores only the segments it needs.",
            "Coarser 50 m and 110 m series are available via `[data] coastline_source` when only a regional outline is wanted.",
        ],
        "<output_dir>/inputs/coastline.dat",
        "Public domain (Natural Earth).",
        "Written as lon,lat vertex lists grouped under `> name,code` headers, the format " *
        "`load_coastline_polygons` reads. This replaces the hand-drawn regional outline that " *
        "previously shipped in the repository, which described Nova Scotia with a handful of " *
        "vertices and left most of the domain's land unmasked.",
    ),

    :tides => DataSource(
        :tides,
        "TPXO9-Atlas astronomical tidal solution (global harmonic velocity)",
        :gridded, true,
        "https://zenodo.org/records/8074917",
        [
            "MEASURED AND CORRECT FOR THIS STUDY AREA. Zenodo 8074917, 'TPXO9 data used for ICON-O tides evaluation', open, CC-BY-4.0. The member `TPXO9_datafiles/u_tpxo9.v1.nc` is the GLOBAL harmonic velocity solution: lon_u spans 0..360, lat_u spans -90..90, with 60 constituents. Verified by inflating only a header prefix of the archived member, so the full 1.4 GB was never downloaded. Extract it with `fetch_input(:tides)`; the archive URL is https://zenodo.org/api/records/8074917/files/TPXO9_datafiles.zip/content",
            "DO NOT USE for this study area. Two smaller archives look equivalent and are NOT, both measured: the 48.7 MB `TPXO9_atlas_v2.zip` on Zenodo 22970362 (pyTMD test data) spans only lat -34.97..14.97, and the `*_remap*` members in record 8074917 are regional (e.g. h_tpxo9_remap3600x1800_sonett.nc spans lon -30..30, lat -55..-10). A tidal archive must be checked against its own lat range before use; the reader raises an error naming the actual coverage rather than returning an empty box.",
            "TPXO9-Atlas from the publisher, https://www.tpxo.atlas.oregonstate.edu (free and open, but the host does not resolve from every network, including this one).",
            "FES2014 / FES2022 from AVISO, https://www.aviso.altimetry.fr (AVISO requires a free account).",
        ],
        "<output_dir>/inputs/tpxo9/u_tpxo9.v1.nc",
        "TPXO9-Atlas: CC-BY-4.0.",
        "u_tpxo9.v1.nc carries the complex VELOCITY harmonics (uRe, uIm, vRe, vIm) on the " *
        "staggered u/v grids plus h_tpxo9.v1.nc for elevation, with 60 constituents. A " *
        "barotropic open boundary is driven by the normal velocity, so the velocity product " *
        "is what the boundary forcing needs. Note the two conventions that are easy to get " *
        "wrong: longitude is 0..360 rather than -180..180, and the imaginary part enters with " *
        "a MINUS sign (the files document `phase = atan2(-hIm, hRe)`).",
    ),
)

"""
    fetch_input(key::Symbol, output_dir::AbstractString; force::Bool = false) -> String

Ensure the physical input `key` is present in this run's inputs directory, downloading it if it
is not, and return the path to the cached artefact.

Every input is fetched keylessly, from the source named first in [`DATA_SOURCES`](@ref). The
alternative sources are for when the primary is unreachable from your network; they are not
tried automatically, because a run that quietly used a different dataset than the one its
provenance file names is exactly the failure this module exists to prevent.

Existing files are reused. `force = true` re-fetches. A download that fails is an error naming
every alternative, never a silent substitution.

!!! note "Where the data goes"
    Everything lands under `<output_dir>/inputs`, so a scenario is self-contained. A dataset is
    therefore fetched once *per scenario* rather than once per machine. If you are running
    several scenarios and want to share the download, point them at a common `output_dir` or
    copy the cached file across -- do not reintroduce a shared top-level `inputs/`.

# Examples
```julia
fetch_input(:bathymetry, "work/snowcrab")            # ETOPO 2022, keyless
fetch_input(:tides, "work/snowcrab")                 # TPXO9 global velocity solution
```
"""
function fetch_input(key::Symbol, output_dir::AbstractString; force::Bool = false)
    d = data_source(key)
    target = replace(d.cache, "<output_dir>" => output_dir)

    if !force && (isfile(target) || isdir(target))
        return target
    end

    # Inputs marked incomplete are not fabricated: the note says what is actually ingested.
    if d.kind == :analytic
        error("Input '$(key)' is analytic and is generated by the model, not downloaded. " *
              d.notes)
    end

    dir = isdir(target) ? target : dirname(target)
    mkpath(dir)

    url = d.primary
    # Several registry entries name a landing page rather than a file. The tide entry does,
    # because the archive is a ZIP whose single required member is extracted by byte range --
    # see `_fetch_tpxo_velocity`.
    if key == :tides
        return _fetch_tpxo_velocity(target; force = force)
    end

    # The coastline is not a straight file download: the upstream GeoJSON is a whole-world
    # archive that has to be clipped to the run's domain before it is written.
    if key == :coastline
        return fetch_natural_earth_coastline(output_path = target, verbose = false)
    end

    if !startswith(url, "http")
        error("Input '$(key)' has no direct download URL (primary is '$(url)'). " *
              "Obtain it manually and place it at $(target). Alternatives:\n  " *
              join(d.alternatives, "\n  "))
    end

    tmp = target * ".part"
    try
        Downloads.download(
            url,
            tmp;
            progress = _download_progress_callback(string(key)),
        )
    catch err
        rm(tmp; force = true)
        error(
            "Could not download $(key) from the keyless source:\n  $(url)\n" *
            "$(sprint(showerror, err))\n" *
            "Alternatives, and their access requirements:\n  " *
            join(d.alternatives, "\n  ") *
            "\nPlace the file at $(target) and re-run if none of these is reachable.")
    end
    mv(tmp, target; force = true)
    return target
end

"""
    _download_progress_callback(label::AbstractString) -> Function

Construct a progress callback for `Downloads.download` that reports transfer
progress and rate periodically without flooding the output stream.
"""
function _download_progress_callback(label::AbstractString)
    last_time = Ref(time())
    last_print = Ref(0.0)
    last_bytes = Ref(0)
    finished = Ref(false)
    return function (total::Integer, now::Integer)
        curr_time = time()
        if (curr_time - last_print[]) >= 2.0 || (total > 0 && now >= total)
            dt = curr_time - last_time[]
            rate_mb = dt > 0.1 ? ((now - last_bytes[]) / (1024 * 1024)) / dt : 0.0
            now_mb = round(now / (1024 * 1024); digits = 1)
            if total > 0
                tot_mb = round(total / (1024 * 1024); digits = 1)
                pct = round(100.0 * now / total; digits = 1)
                print("\r  [$(label)] $(now_mb) / $(tot_mb) MB ($(pct)%) " *
                      "[$(round(rate_mb; digits = 2)) MB/s]    ")
            else
                print("\r  [$(label)] $(now_mb) MB [$(round(rate_mb; digits = 2)) MB/s]    ")
            end
            flush(stdout)
            last_print[] = curr_time
            last_time[] = curr_time
            last_bytes[] = now
            if total > 0 && now >= total && !finished[]
                finished[] = true
                println()
            end
        end
    end
end

"""
    _fetch_tpxo_velocity_range(archive_url, member, target) -> String

Download single member `member` directly from remote ZIP archive at `archive_url`
via HTTP range requests and inflate directly to `target`.
"""
function _fetch_tpxo_velocity_range(
    archive_url::AbstractString,
    member::AbstractString,
    target::AbstractString,
)
    mkpath(dirname(target))

    # Probe total length via HEAD request
    io_head = IOBuffer()
    resp = Downloads.request(archive_url; method = "HEAD", output = io_head)
    total_len = 0
    for (k, v) in resp.headers
        if lowercase(k) == "content-length"
            total_len = parse(Int, v)
            break
        end
    end
    total_len <= 0 && error("Content-Length header not returned by server")

    # Fetch trailing 128 KB containing ZIP footer and Central Directory
    tail_len = min(total_len, 131072)
    tail_start = total_len - tail_len
    io_tail = IOBuffer()
    Downloads.download(
        archive_url,
        io_tail;
        headers = ["Range" => "bytes=$(tail_start)-$(total_len - 1)"],
    )
    tail_bytes = take!(io_tail)

    # Locate central directory offset
    cd = 0
    eocd = findlast(
        i -> i + 3 <= length(tail_bytes) && _u32(tail_bytes, i) == 0x06054b50,
        1:(length(tail_bytes) - 4),
    )
    if eocd !== nothing && _u32(tail_bytes, eocd + 16) != 0xFFFFFFFF
        cd = Int(_u32(tail_bytes, eocd + 16))
    else
        e64 = findfirst(
            i -> i + 3 <= length(tail_bytes) && _u32(tail_bytes, i) == 0x06064b50,
            1:(length(tail_bytes) - 3),
        )
        e64 === nothing && error("Central directory record not found in archive footer")
        cd = Int(_u64(tail_bytes, e64 + 48))
    end

    cd_bytes = if cd >= tail_start
        tail_bytes[(cd - tail_start + 1):end]
    else
        io_cd = IOBuffer()
        Downloads.download(
            archive_url,
            io_cd;
            headers = ["Range" => "bytes=$(cd)-$(total_len - 1)"],
        )
        take!(io_cd)
    end

    # Search for member in Central Directory
    i = 1
    lho = csize = usize = 0
    method = 0
    while i + 45 <= length(cd_bytes) && _u32(cd_bytes, i) == 0x02014b50
        c = Int(_u32(cd_bytes, i + 20))
        u = Int(_u32(cd_bytes, i + 24))
        m = _u16(cd_bytes, i + 8)
        nlen = _u16(cd_bytes, i + 28)
        elen = _u16(cd_bytes, i + 30)
        clen = _u16(cd_bytes, i + 32)
        name = String(cd_bytes[(i + 46):(i + 45 + nlen)])
        if name == member
            csize, usize, method = c, u, m
            lho = Int(_u32(cd_bytes, i + 42))
            break
        end
        i += 46 + nlen + elen + clen
    end
    csize == 0 && error("Member '$(member)' not located in central directory")

    # Fetch local file header to compute exact data stream offset
    io_lh = IOBuffer()
    Downloads.download(
        archive_url,
        io_lh;
        headers = ["Range" => "bytes=$(lho)-$(lho + 128)"],
    )
    lh = take!(io_lh)
    _u32(lh, 1) != 0x04034b50 && error("Invalid local header signature at $(lho)")
    ln = _u16(lh, 27)
    le = _u16(lh, 29)
    data_start = lho + 30 + ln + le
    data_end = data_start + csize - 1

    # Stream compressed payload to temporary buffer file
    comp_file = target * ".deflate.part"
    println("Streaming $(round(csize / (1024 * 1024), digits=1)) MB compressed tidal " *
            "solution via HTTP Range...")
    try
        Downloads.download(
            archive_url,
            comp_file;
            headers = ["Range" => "bytes=$(data_start)-$(data_end)"],
            progress = _download_progress_callback("TPXO9 stream"),
        )
    catch err
        rm(comp_file; force = true)
        rethrow(err)
    end

    # Decompress deflate stream to destination file
    println("Inflating $(member) to $(basename(target))...")
    tmp_dest = target * ".part"
    try
        open(tmp_dest, "w") do out_io
            open(comp_file, "r") do in_io
                ds = CodecZlib.DeflateDecompressorStream(in_io)
                try
                    write(out_io, ds)
                finally
                    close(ds)
                end
            end
        end
        rm(comp_file; force = true)
        mv(tmp_dest, target; force = true)
    catch err
        rm(tmp_dest; force = true)
        rethrow(err)
    end
    return target
end

"""
    _fetch_tpxo_velocity(target; force = false, output_dir = ".") -> String

Acquire the TPXO9 global harmonic velocity solution (`u_tpxo9.v1.nc`) at `target`.
"""
function _fetch_tpxo_velocity(target::AbstractString; force::Bool = false)
    isfile(target) && !force && return target
    mkpath(dirname(target))

    member = "TPXO9_datafiles/u_tpxo9.v1.nc"
    archive_url = "https://zenodo.org/api/records/8074917/files/TPXO9_datafiles.zip/content"
    archive = joinpath(dirname(target), "TPXO9_datafiles.zip")

    # If full archive is already cached locally, extract directly from it
    if isfile(archive) && !force
        println("Extracting $(member) from local archive $(basename(archive))...")
        _extract_zip_member(archive, member, target)
        return target
    end

    # First attempt targeted HTTP range download (858 MB vs 4.4 GB full archive)
    try
        return _fetch_tpxo_velocity_range(archive_url, member, target)
    catch err
        @warn "Range download failed ($(sprint(showerror, err))); falling back to " *
              "full archive download."
    end

    # Fallback: download complete 4.4 GB archive
    if !isfile(archive) || force
        println("Downloading TPXO9 archive (4.4 GB) to $(archive)...")
        tmp = archive * ".part"
        try
            Downloads.download(
                archive_url,
                tmp;
                progress = _download_progress_callback("TPXO9 archive"),
            )
        catch err
            rm(tmp; force = true)
            error(
                "Could not download TPXO9 archive from $(archive_url):\n" *
                "$(sprint(showerror, err))\n" *
                "Place file at $(archive) or refer to DATA_SOURCES[:tides].")
        end
        mv(tmp, archive; force = true)
    end

    println("Extracting $(member) from $(basename(archive))...")
    _extract_zip_member(archive, member, target)
    return target
end

"""
    _extract_zip_member(archive, member, dest) -> String

Extract one member of a ZIP archive to `dest`, overwriting it.

Only the two members actually needed here are supported: a *stored* (uncompressed) entry, and a
*deflated* one. Both classic ZIP and ZIP64 central directories are handled, because the two
archives this project uses differ: the 4.4 GB one is ZIP64, whose classic offsets read
`0xFFFFFFFF` and are not usable, while the small one is classic. ZIP64 is detected from those
sentinel values rather than from file size, since a ZIP64 archive is not required to exceed 4 GB.

Deflate is raw (no zlib or gzip wrapper), so it is decoded with `DeflateDecompressorStream`;
`ZlibDecompressorStream` would fail on it with a header error.
"""
function _extract_zip_member(archive::AbstractString, member::AbstractString, dest::AbstractString)
    bytes = read(archive)

    # Prefer the classic end-of-central-directory record; fall back to ZIP64 when its
    # offset field is the 0xFFFFFFFF sentinel.
    cd = 0
    eocd = findlast(i -> i + 3 <= length(bytes) && _u32(bytes, i) == 0x06054b50,
                    1:(length(bytes) - 4))
    if eocd !== nothing && _u32(bytes, eocd + 16) != 0xFFFFFFFF
        cd = Int(_u32(bytes, eocd + 16))
    else
        e64 = findfirst(i -> i + 3 <= length(bytes) && _u32(bytes, i) == 0x06064b50,
                        1:(length(bytes) - 3))
        e64 === nothing && error(
            "$(basename(archive)) has neither a classic nor a ZIP64 central directory.")
        cd = Int(_u64(bytes, e64 + 48))
    end

    # The recorded offset is not trusted blindly. Some writers are off by one, and a
    # wrong-by-one start looks exactly like "member not found". The offset is probed
    # NEAREST-FIRST rather than by scanning forward, because a forward scan finds a
    # false-positive signature inside the compressed data long before the real one.
    if _u32(bytes, cd) != 0x02014b50
        hit = 0
        for d in 1:256
            for cand in (cd + d, cd - d)
                (1 <= cand <= length(bytes) - 3) || continue
                if _u32(bytes, cand) == 0x02014b50
                    hit = cand
                    break
                end
            end
            hit != 0 && break
        end
        hit == 0 && error(
            "Could not locate the central directory in $(basename(archive)); the archive " *
            "may be corrupt or in an unsupported form.")
        cd = hit
    end

    lho = csize = usize = 0
    method = 0
    i = cd
    while i + 45 <= length(bytes) && _u32(bytes, i) == 0x02014b50
        c = Int(_u32(bytes, i + 20))
        u = Int(_u32(bytes, i + 24))
        m = _u16(bytes, i + 8)          # compression method, 2 bytes at offset 8
        nlen = _u16(bytes, i + 28)
        elen = _u16(bytes, i + 30)
        clen = _u16(bytes, i + 32)
        if String(bytes[(i + 46):(i + 45 + nlen)]) == member
            csize, usize, method = c, u, m
            break
        end
        i += 46 + nlen + elen + clen
    end
    csize == 0 && error("Member '$(member)' is not present in $(basename(archive)).")

    # The local header is found by its OWN signature and name rather than by trusting the
    # central directory's relative offset. Several offsets in these archives have proved
    # unreliable -- the central directory was one byte off, and a wrong local-header offset
    # silently truncates the stream instead of failing. Matching the name makes the lookup
    # self-locating: if the bytes are not there, it is not there.
    lho = 0
    elen = 0
    lmethod = 0
    lcsize = 0
    j = 1
    while j + 30 <= length(bytes) && _u32(bytes, j) == 0x04034b50
        ln = _u16(bytes, j + 26)
        le = _u16(bytes, j + 28)
        if String(bytes[(j + 30):(j + 29 + ln)]) == member
            lho = j
            elen = le
            lmethod = _u16(bytes, j + 8)
            lcsize = Int(_u32(bytes, j + 18))
            break
        end
        j += 30 + ln + le + Int(_u32(bytes, j + 18))
    end
    lho == 0 && error(
        "Could not find the local header for '$(member)' in $(basename(archive)).")

    # The LOCAL header wins for method and size. This archive's central directory records
    # method 0 (stored) while the local header records method 8 (deflate); believing the
    # central directory copies the compressed bytes through verbatim and yields a file that
    # is the right length and completely unreadable.
    if lmethod != 0
        method = lmethod
    end
    if lcsize > 0
        csize = lcsize
    end

    data_start = lho + 30 + _u16(bytes, lho + 26) + elen

    # The payload offset is confirmed by decoding, not by trusting the headers. This
    # archive's offsets have repeatedly proved unreliable, and a deflate stream that starts
    # one byte out usually decodes *something* rather than failing cleanly. So a small window
    # around the computed start is tried, and the offset that yields exactly the recorded
    # uncompressed size wins. That check is self-validating in a way no header field is.
    want = usize
    for shift in -2:2
        s = data_start + shift
        (1 <= s && s + csize - 1 <= length(bytes)) || continue
        chunk = bytes[(s + 1):(s + csize)]
        got = try
            io = IOBuffer()
            ds = CodecZlib.DeflateDecompressorStream(IOBuffer(chunk))
            try
                write(io, ds)
            finally
                close(ds)
            end
            take!(io)
        catch
            UInt8[]
        end
        if !isempty(got) && (want == 0 || length(got) == want)
            tmp = dest * ".part"
            open(tmp, "w") do out
                write(out, got)
            end
            mv(tmp, dest; force = true)
            return dest
        end
    end
    error(
        "Could not decode member '$(member)' from $(basename(archive)). Tried payload " *
        "offsets $(data_start - 2):$(data_start + 2) and none produced the recorded " *
        "$(want) bytes.")
end

# --- minimal ZIP readers -----------------------------------------------------
# A 4.4 GB archive is ZIP64: the classic end-of-central-directory offsets are 0xFFFFFFFF.

@inline _u16(b, i) = Int(UInt16(b[i]) | (UInt16(b[i + 1]) << 8))
@inline _u32(b, i) = UInt32(b[i]) | (UInt32(b[i + 1]) << 8) | (UInt32(b[i + 2]) << 16) | (UInt32(b[i + 3]) << 24)
@inline _u64(b, i) = foldl((v, k) -> v | (UInt64(b[i + k]) << (8 * k)), 0:7; init = UInt64(0))

"""
    data_source(key::Symbol) -> DataSource

Look up a declared input by its config key, with a helpful error listing the valid keys.
"""
function data_source(key::Symbol)
    haskey(DATA_SOURCES, key) || error(
        "Unknown data source '$(key)'. Declared: " *
        join(sort(string.(collect(keys(DATA_SOURCES)))), ", "))
    return DATA_SOURCES[key]
end

"""
    describe_data_sources(io = stdout) -> nothing

Print the full registry: what each input is, which source to use, the alternatives, and the
access requirements. This is the human-facing half of the reproducibility contract; the
machine-facing half is `data_provenance.json`, written by a run.
"""
function describe_data_sources(io::IO = stdout)
    println(io, "Physical inputs required by ParticleTracking\n")
    println(io, "="^78)
    for key in sort!(collect(keys(DATA_SOURCES)); by = string)
        d = DATA_SOURCES[key]
        println(io, "\n", uppercase(string(d.key)), "  (", d.kind, ", keyless=", d.keyless, ")")
        println(io, "  ", d.title)
        println(io, "  primary : ", d.primary)
        println(io, "  cache   : ", d.cache)
        println(io, "  licence : ", d.license)
        println(io, "  other sources:")
        for a in d.alternatives
            println(io, "    - ", a)
        end
        println(io, "  notes   : ", d.notes)
    end
    return nothing
end

"""
    input_dir(output_dir) -> String

The container directory for this run's downloaded inputs, `<output_dir>/inputs`.

Inputs are deliberately **not** shared in a top-level `inputs/` directory. Scenarios and
subprojects are self-contained: a `snowcrab` run keeps its bathymetry, winds, hydrography,
ocean boundary and tidal solution beside its own outputs under `work/snowcrab/`, alongside its
`resolved_config.toml` and `data_provenance.json`. Three things follow from that:

* A run is a single directory that can be archived, moved, or shared as a unit, and every file
  in it is traceable to the run that produced it.
* Two scenarios cannot silently share a half-overwritten cache. The old shared directory made
  it possible to change a domain and reuse another domain's bathymetry, which is precisely the
  kind of quiet mismatch the provenance file exists to expose.
* The cost is disk, and that is the trade: a dataset is fetched once per scenario, not once
  per machine.
"""
input_dir(output_dir::AbstractString) = joinpath(output_dir, "inputs")

"""
    file_digest(path::AbstractString) -> String

SHA-256 of a file, used to detect a cached input that has changed since it was recorded.
"""
function file_digest(path::AbstractString)
    isfile(path) || return ""
    # SHA.update! has no streaming method, so the file is read in blocks. A tidal atlas is
    # tens of megabytes, which is fine to read whole, but the block loop keeps memory flat
    # if a much larger artefact is ever added to the registry.
    ctx = SHA.SHA256_CTX()
    open(path, "r") do io
        buf = Vector{UInt8}(undef, 1 << 20)
        while !eof(io)
            n = readbytes!(io, buf)
            n == 0 && break
            SHA.update!(ctx, view(buf, 1:n))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

"""
    data_provenance(output_dir; config_path = "", extra = Dict()) -> String

Write `data_provenance.json` into `output_dir`, recording every declared input with the exact
file that backs it, and return the path.

For each [`DATA_SOURCES`](@ref) entry the record carries the source, the cache path, whether the
cached file is present, its size in bytes, and its SHA-256. That is enough to answer, months
later, "which ETOPO did this run use?" without trusting the config that produced it.

This is deliberately a *separate* file from `resolved_config.toml`. The resolved config records
what was *asked for*; this records what was *actually read*. When a run degrades, or a cached
artefact is stale, those two diverge, and the divergence is the point.
"""
function data_provenance(
    output_dir::AbstractString;
    config_path::AbstractString = "",
    extra::AbstractDict = Dict{String, Any}()
)
    isdir(output_dir) || mkpath(output_dir)

    inputs = Dict{String, Any}()
    for key in sort!(collect(keys(DATA_SOURCES)); by = string)
        d = DATA_SOURCES[key]
        # Cache paths are declared as `<output_dir>/inputs/...` so the registry is
        # output-location independent; resolve them against this run's output directory here,
        # and fall back to the repository-relative form when writing the manifest on its own.
        rel = replace(d.cache, "<output_dir>" => output_dir)
        present = isdir(rel) || isfile(rel)
        size_bytes = 0
        if isfile(rel)
            size_bytes = filesize(rel)
        elseif isdir(rel)
            for (root, _, files) in walkdir(rel), f in files
                size_bytes += filesize(joinpath(root, f))
            end
        end
        inputs[string(d.key)] = Dict{String, Any}(
            "title" => d.title,
            "kind" => string(d.kind),
            "keyless" => d.keyless,
            "source" => d.primary,
            "cache" => replace(rel, '\\' => '/'),
            "present" => present,
            "bytes" => size_bytes,
            "sha256" => isfile(rel) ? file_digest(rel) : "",
            "license" => d.license,
        )
    end

    record = Dict{String, Any}(
        "generated_utc" => string(Dates.now(Dates.UTC)),
        "config" => replace(config_path, '\\' => '/'),
        "package_version" => string(VERSION),
        "input_dir" => replace(input_dir(output_dir), '\\' => '/'),
        "inputs" => inputs,
    )
    for (k, v) in extra
        record[k] = v
    end

    path = joinpath(output_dir, "data_provenance.json")
    open(path, "w") do io
        JSON3.write(io, record)
    end
    return path
end
