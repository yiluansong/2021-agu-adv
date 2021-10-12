# set study area to northern mid- to high-latitudes
e <- extent(-180, 180, 30, 70)

# load QA (quality assessment) and number of season (NOS) masks
file <- c("./data/Masks_0.5deg.nc")#World
masks0.5 <- brick(file)
QAmask0.5<-masks0.5[[1]]
NOSmask0.5<-masks0.5[[2]]

# mask for continents 
file<-"./data/ne_10m_geography_regions_polys/ne_10m_geography_regions_polys.shp"
layer <- ogrListLayers(file)
continent <- readOGR(file, layer=layer)
continentras<-rasterize(continent,raster("./data/MAT.nc"), field="region")
Asmask<-Eumask<-NAmemask<-continentras

Asmask[Asmask!=3]<-NA
Asmask[!is.na(Asmask)]<-1

Eumask[Eumask!=4]<-NA
Eumask[!is.na(Eumask)]<-1

NAmemask[NAmemask!=5]<-NA
NAmemask[!is.na(NAmemask)]<-1
