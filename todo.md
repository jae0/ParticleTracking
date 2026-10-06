



create data handling system: make it simple and abstract away underlying storage details, agnostic using the following as defaults but able to add others easily:
 
 - replace files operations with Zarr.jl + YAXArrays.jl + GeoParquet
 - ingest netcdf, hdf, geojson, shp, etc and save locally as zarr /yaxarrays /  and sf polygons as GeoParquet 

 - save outputs as parquet, zarr /yaxarrays instead of duckdb/parquet/jld2 (unless for specific use cases)
 - make the data storage and access system simple and easily used by other projects, with minimal dependencies .. and a common accessible file location 

 - ensure that the data is compressed efficiently and can be used by R or python as well (? backend = :nczarr, # Forces xarray/netcdf compatibility)

- abstract data query / slice / index / join / save, etc and make the same code usable across all backends ... make it as simple and generic as possible ... use functions that read write / save / load / query / slice / index / join / etc and use these functions everywhere instead of directly accessing the data ... use sensible default backends ... and allow users to specify backends when needed 


eventually move movementanalysis here and use the temperature and salinity data to provide hsi information for snow crab

