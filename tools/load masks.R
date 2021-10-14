# set study area to northern mid- to high-latitudes
e <- extent(-180, 180, 30, 70)

# load QA (quality assessment) and number of season (NOS) masks
file <- c("./data/Masks_0.5deg.nc") # World
masks0.5 <- brick(file)
QAmask0.5 <- masks0.5[[1]]
NOSmask0.5 <- masks0.5[[2]]

# mask for continents
file <- "./data/ne_10m_geography_regions_polys/ne_10m_geography_regions_polys.shp"
layer <- ogrListLayers(file)
continent <- readOGR(file, layer = layer)
Asmask <- rasterize(x = continent[continent@data$region == "Asia", ], y = raster("./data/MAT.nc"), field = 1)
Eumask <- rasterize(x = continent[continent@data$region == "Europe", ], y = raster("./data/MAT.nc"), field = 1)
NAmemask <- rasterize(x = continent[continent@data$region == "North America", ], y = raster("./data/MAT.nc"), field = 1)
