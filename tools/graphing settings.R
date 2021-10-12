# levelplot layout
my.settings <- list(
  layout.heights = list(
    bottom.padding = -1,
    top.padding = -1,
    left.padding = -1,
    right.padding = -1,
    key.sub.padding = -1,
    axis.xlab.padding = 0.5,
    axis.ylab.padding = -1,
    key.axis.padding = -1,
    main.key.padding = -1
  ),
  aspect.fill = TRUE,
  par.main.text = list(just = "left", x = grid::unit(10, "mm"), y = grid::unit(10, "mm"))
)

# label legend at selected values
labelpts <- function(layer,vector) {
  labelpts<-rep(NA,length(vector))
  for (i in 1:length(vector)) {
    labelpts[i]<-mean(na.omit(values(layer)) < vector[i])
  }
  return (labelpts)
}

# retrieve jpeg schematic diagrams
get_schematic<-function(season) {
  schematic<-vector(mode="list", length=4)
  if (season=="year") {
    # LOS vs MAR
    schematic[[1]]<-as.raster(readJPEG("./figures/schematics/schematic1.jpg", native = FALSE))
    schematic[[2]]<-as.raster(readJPEG("./figures/schematics/schematic2.jpg", native = FALSE))
    schematic[[3]]<-as.raster(readJPEG("./figures/schematics/schematic3.jpg", native = FALSE))
    schematic[[4]]<-as.raster(readJPEG("./figures/schematics/schematic4.jpg", native = FALSE))
  }
  if(season=="spring") {
    #SOS vs MST
    schematic[[1]]<-as.raster(readJPEG("./figures/schematics/schematic5.jpg", native = FALSE))
    schematic[[2]]<-as.raster(readJPEG("./figures/schematics/schematic6.jpg", native = FALSE))
    schematic[[3]]<-as.raster(readJPEG("./figures/schematics/schematic7.jpg", native = FALSE))
    schematic[[4]]<-as.raster(readJPEG("./figures/schematics/schematic8.jpg", native = FALSE))
  }
  if (season=="fall") {
    #EOS vs MFT
    schematic[[1]]<-as.raster(readJPEG("./figures/schematics/schematic9.jpg", native = FALSE))
    schematic[[2]]<-as.raster(readJPEG("./figures/schematics/schematic10.jpg", native = FALSE))
    schematic[[3]]<-as.raster(readJPEG("./figures/schematics/schematic11.jpg", native = FALSE))
    schematic[[4]]<-as.raster(readJPEG("./figures/schematics/schematic12.jpg", native = FALSE))
  }
  
  if (season=="year_prcp") {
    #LOS vs MAP
    schematic[[1]]<-as.raster(readJPEG("./figures/schematics/schematic17.jpg", native = FALSE))
    schematic[[2]]<-as.raster(readJPEG("./figures/schematics/schematic18.jpg", native = FALSE))
    schematic[[3]]<-as.raster(readJPEG("./figures/schematics/schematic19.jpg", native = FALSE))
    schematic[[4]]<-as.raster(readJPEG("./figures/schematics/schematic20.jpg", native = FALSE))
  }
  return (schematic)
}

# landmass layer
file<-"./data/ne_10m_land/ne_10m_land.shp"
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/
layer <- ogrListLayers(file)
land <- readOGR(file, layer=layer)
