"""
    manifest.jl

Provenance registry for every physical input a run consumes.

This module re-exports the GeoData manifest framework. All data source
declarations and specialized fetch functions are now in GeoData.Data.GeoDataManifest.
"""
module Manifest

using ..GeoData.Data.GeoDataManifest

# Re-export all public functions from GeoDataManifest
export DATA_SOURCES, DataSource, fetch_input, input_dir, file_digest, data_provenance, data_source, describe_data_sources

end # module Manifest