"""
    open_data.jl

ParticleTracking-specific data acquisition functions not covered by GeoData.

This module contains:
- Configuration-dependent domain resolution
- Source dispatchers that use ParticleTracking configuration

All generic geospatial data operations (bathymetry, winds, WOA23, coastlines, regridding)
are now in GeoData.Data.
Copernicus Marine and synthetic data functions are now in GeoData.Data.Copernicus.
"""

using Downloads
using TOML

"""
    resolve_domain_bounds(lon_range, lat_range; embedding = false) -> NamedTuple

Resolve geographic bounds for a data-retrieval function, reading them from the active TOML
configuration when not supplied explicitly.

No geographic domain is compiled into this package: the `[domain]` section defines the study
region and the `[boundaries]` section defines the larger embedding region used for boundary
and forcing downloads. Passing `nothing` for either range therefore defers to the
configuration file, which is the single source of truth.
"""
function resolve_domain_bounds(lon_range, lat_range; embedding::Bool = false)
    rngs = if isnothing(lon_range) || isnothing(lat_range)
        cfg = try
            load_configuration(find_default_config_path())
        catch err
            error(
                "No lon_range/lat_range supplied and the active configuration could not be " *
                "read ($(err)). Supply the bounds explicitly, or fix the TOML file."
            )
        end
        embedding ? embedding_domain_ranges(cfg) : study_domain_ranges(cfg)
    else
        (lon = (Float64(lon_range[1]), Float64(lon_range[2])),
         lat = (Float64(lat_range[1]), Float64(lat_range[2])))
    end
    return rngs
end

"""
    fetch_surface_winds(source; lon_range, lat_range, time_iso, output_path, verbose) -> String

Fetch surface winds from exactly one named source, and fail loudly if it is not one this
package can actually retrieve. Currently only `open_meteo` is implemented.
"""
function fetch_surface_winds(source::AbstractString; lon_range, lat_range, time_iso,
                            output_path, verbose::Bool = true)
    s = replace(lowercase(strip(source)), "-" => "_")
    if s in ("open_meteo", "era5", "open", "openmeteo")
        return GeoData.Data.fetch_open_meteo_winds(
            lon_range = lon_range, lat_range = lat_range,
            time_iso = time_iso, output_path = output_path,
            verbose = verbose
        )
    end
    error(_unknown_source_message("wind_source", s,
        ["open_meteo" => "ERA5 via the keyless Open-Meteo archive API (only implemented source)"]))
end

"""
    fetch_bathymetry(source; lon_range, lat_range, output_path, verbose) -> String

Fetch bathymetry from exactly one named source. Currently only `etopo2022` is implemented.
"""
function fetch_bathymetry(source::AbstractString; lon_range, lat_range, output_path,
                           verbose::Bool = true)
    s = replace(lowercase(strip(source)), "-" => "_")
    if s in ("etopo2022", "etopo", "open")
        return GeoData.Data.fetch_erddap_bathymetry(
            lon_range = lon_range, lat_range = lat_range,
            output_path = output_path, verbose = verbose
        )
    end
    error(_unknown_source_message("bathy_source", s,
        ["etopo2022" => "ETOPO 2022 15-arcsec global relief via the keyless NOAA ERDDAP"]))
end

"""
    _unknown_source_message(key, given, choices) -> String

Build the error for an unrecognised `[data] *_source` value.
"""
function _unknown_source_message(key::AbstractString, given::AbstractString,
                                choices::Vector{Pair{String, String}})
    io = IOBuffer()
    println(io, "[data] $key = \"$given\" is not a source this project can retrieve.")
    println(io, "Available sources:")
    for (name, note) in choices
        println(io, "  $name - ", note)
    end
    println(io, "Set [data] $key in the TOML and re-run. Nothing was downloaded, and no ",
            "substitute was written in its place.")
    String(take!(io))
end

end # module