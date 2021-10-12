# Code for 2.1.2 (leading to Figure 2 and supplementary figures)

source("./tools/packages.R")
source("./tools/load masks.R")
source("./tools/graphing settings.R")

##### Calculate climate-phenology metrics (Figure 2, S9, S10, S12)
season_list<-c ("year", #  (for MAT vs GSL Figure 2)
                "spring", # (for MST vs SOS Figure S9)
                "fall", # (for MFT vs EOS Figure S10)
                "year_prcp") # (for MAP vs GSL Figure S12)

for (season in season_list) {
  ### Vc
  if (season=="year") {
    file <- c("./data/MAT.nc") # mean annual temperature (for Figure 2)
  }
  if (season=="spring") {
    file <- c("./data/Tspring.nc") # spring temperature (for Figure S9)
  }
  if (season=="fall") {
    file <- c("./data/Tfall.nc") # fall temperature (for Figure S10)
  }
  if (season=="year_prcp") {
    file <- c("./data/MAP.nc") # mean annual precipitation (for Figure S12)
  }
  
  C <- brick(file)
  C <- C[[81:114]]
  
  # temporal gradient
  vt <- tempTrend(C, th = 10)
  # spatial gradient
  vg <- spatGrad(C, th = -Inf, projected = FALSE)
  # velocity
  gv <- gVoCC(vt, vg)
  
  gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
  
  # decompose
  Vc <- abs(gv[[1]])
  AngC <- gv[[2]]
  VcX <- Vc * sin(AngC * pi / 180)
  VcY <- Vc * cos(AngC * pi / 180)
  
  # summarize
  quantile(na.omit(values(Vc)), c(0.025, 0.5, 0.975))
  sum(na.omit(values(VcY)) > 0) / length(na.omit(values(VcY)))
  sum(na.omit(values(VcY)) < 0) / length(na.omit(values(VcY)))
  
  ### Vp
  if (season=="year"|season=="year_prcp") {
    file <- c("./data/LOS_World_0.5deg_QA.nc") # length of season  (for Figure 2 and Figure S12)
  }
  if (season=="spring") {
    file <- c("./data/SOS_World_0.5deg_QA.nc") # start of season  (for Figure S9)
  }
  if (season=="fall") {
    file <- c("./data/EOS_World_0.5deg_QA.nc") # end of season (for Figure S10)
  }
  
  P <- brick(file)
  
  # temporal gradient
  vt <- tempTrend(P, th = 10)
  # spatial gradient
  vg <- spatGrad(P, th = -Inf, projected = FALSE)
  # velocity
  gv <- gVoCC(vt, vg)
  
  gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
  
  # decompose
  Vp <- abs(gv[[1]])
  AngP <- gv[[2]]
  VpX <- Vp * sin(AngP * pi / 180)
  VpY <- Vp * cos(AngP * pi / 180)
  
  # summarize
  quantile(na.omit(values(Vp)), c(0.025, 0.5, 0.975))
  sum(na.omit(values(VpY)) > 0) / length(na.omit(values(VpY)))
  sum(na.omit(values(VpY)) < 0) / length(na.omit(values(VpY)))
  
  ### Compare Vp to Vc
  projection <- (VcX * VpX + VcY * VpY) / Vc
  anglediff <- abs(((AngC - AngP) + 180) %% 360 - 180)
  mismatch <- projection-Vc
  
  quantile(na.omit(values(anglediff)), c(0.025, 0.5, 0.975))
  sum(na.omit(values(anglediff)) > 90) / length(na.omit(values(anglediff)))
  
  quantile(na.omit(values(mismatch)), c(0.025, 0.5, 0.975))
  sum(na.omit(values(mismatch)) > 0) / length(na.omit(values(mismatch)))
  sum(na.omit(values(mismatch)) < 0) / length(na.omit(values(mismatch)))
  
  ##### Plot
  ### Maps of four metrics
  cutpts <- quantile(na.omit(values(stack(Vc,Vp))), probs = seq(0, 1, length.out = 21))
  fieldC <- stack(Vc / Vc, AngC)
  names(fieldC) <- c("slope", "aspect")
  bg_Vc <- levelplot(Vc, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = rev(colorspace::heat_hcl(20)), main = "", 
                     xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 1, 2, 5, 10), at = labelpts(stack(Vc,Vp), c(0, 1, 2, 5, 10)))))
  arrow_Vc <- vectorplot(fieldC, isField = TRUE, unit = "degrees", col.arrows = "black", narrows = 500, aspX = 4, aspY = 4, lwd.arrows = 1, alpha = 1, reverse = FALSE, colorkey = FALSE, region = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(draw = FALSE))
  p_C <- bg_Vc + arrow_Vc
  # p_C + latticeExtra::layer(sp.polygons(land))
  
  fieldP <- stack(Vp / Vp, AngP)
  names(fieldP) <- c("slope", "aspect")
  bg_Vp <- levelplot(Vp, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = rev(colorspace::heat_hcl(20)), main = "", 
                     xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 1, 2, 5, 10), at = labelpts(stack(Vc,Vp), c(0, 1, 2, 5, 10)))))
  arrow_Vp <- vectorplot(fieldP, isField = TRUE, unit = "degrees", col.arrows = "black", narrows = 500, aspX = 4, aspY = 4, lwd.arrows = 1, alpha = 1, reverse = FALSE, colorkey = FALSE, region = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(draw = FALSE))
  p_P <- bg_Vp + arrow_Vp
  # p_P + latticeExtra::layer(sp.polygons(land))
  
  if (season=="year") {
    col_pal_diverge <- c(colorspace::sequential_hcl(13, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(7, h = 330, l = c(40, 90), power = 1))) # whole year (for Figure 2)
  }
  if (season=="spring") {
    col_pal_diverge <- c(colorspace::sequential_hcl(12, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(8, h = 330, l = c(40, 90), power = 1))) #spring  (for Figure S9)
  }
  if (season=="fall") {
    col_pal_diverge <- c(colorspace::sequential_hcl(13, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(7, h = 330, l = c(40, 90), power = 1))) #fall (for Figure S10)
  }
  if (season=="year_prcp") {
    col_pal_diverge <- c(colorspace::sequential_hcl(10, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(10, h = 330, l = c(40, 90), power = 1))) #prcp  (for Figure S12)
  }
  cutpts <- quantile(na.omit(values(anglediff)), probs = seq(0, 1, length.out = 21))
  p_anglediff <- levelplot(anglediff, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = col_pal_diverge, main = "", 
                           xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 45, 90, 135, 180), at = labelpts(anglediff, c(0, 45, 90, 135, 180)))))
  # p_anglediff + latticeExtra::layer(sp.polygons(land))
  
  if (season=="year") {
    col_pal_diverge <- rev(c(colorspace::sequential_hcl(9, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(11, h = 330, l = c(40, 90), power = 1)))) # whole year (for Figure 2)
  }
  if (season=="spring") {
    col_pal_diverge <- rev(c(colorspace::sequential_hcl(6, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(14, h = 330, l = c(40, 90), power = 1)))) #spring (for Figure S9)
  }
  if (season=="fall") {
    col_pal_diverge <- rev(c(colorspace::sequential_hcl(7, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(13, h = 330, l = c(40, 90), power = 1)))) #fall (for Figure S10)
  }
  if (season=="year_prcp") {
    col_pal_diverge <- rev(c(colorspace::sequential_hcl(6, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(14, h = 330, l = c(40, 90), power = 1)))) #prcp (for Figure S12)
  }
  
  cutpts <- quantile(na.omit(values(mismatch)), probs = seq(0, 1, length.out = 21))
  p_mismatch <- levelplot(mismatch, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = col_pal_diverge, main = "", 
                          # xlab = list(label = expression("Spatial lag between two velocities (km yr"^-1 * ")"), cex = 0.8), 
                          xlab = list(cex = 0),ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(-10, -5, -2, -1, 0, 1, 2, 5, 10), at = labelpts(mismatch, c(-10, -5, -2, -1, 0, 1, 2, 5, 10)))))
  # p_mismatch + latticeExtra::layer(sp.polygons(land))
  
  if (season == "year") {
    ### Compare whole-year statistics with spring and fall statistics
    df_season<-data.frame(lower=rep(NA, 3), median=rep(NA, 3), upper=rep(NA, 3))
    df_season_list<-vector(mode="list")
    for (i in c("Vc", "Vp", "anglediff", "mismatch")) {
      df_season_list[[i]]<-df_season
    }
    for (i in 1:3) {
      if (i==1) {
        file <- c("./data/MAT.nc")
        C <- brick(file)
        C <- C[[81:114]]
        
        file <- c("./data/LOS_World_0.5deg_QA.nc")
        P <- brick(file)
      }
      if (i==2) {
        file <- c("./data/Tspring.nc")
        C <- brick(file)
        C <- C[[81:114]]
        
        file <- c("./data/SOS_World_0.5deg_QA.nc")
        P <- brick(file)
      }
      if (i==3) {
        file <- c("./data/Tfall.nc")
        C <- brick(file)
        C <- C[[81:114]]
        
        file <- c("./data/EOS_World_0.5deg_QA.nc")
        P <- brick(file)
      }
      
      vt <- tempTrend(C, th = 10)
      vg <- spatGrad(C, th = -Inf, projected = FALSE)
      gv <- gVoCC(vt, vg)
      gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
      Vc <- abs(gv[[1]])
      AngC <- gv[[2]]
      gradangleC<-vg[[2]]
      VcX <- Vc * sin(AngC * pi / 180)
      VcY <- Vc * cos(AngC * pi / 180)
      
      df_season_list[["Vc"]][i,]<-quantile(na.omit(values(Vc)), c(0.025, 0.5, 0.975))
      
      vt <- tempTrend(P, th = 10)
      vg <- spatGrad(P, th = -Inf, projected = FALSE)
      gv <- gVoCC(vt, vg)
      gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
      Vp <- abs(gv[[1]])
      AngP <- gv[[2]]
      gradangleP<-vg[[2]]
      VpX <- Vp * sin(AngP * pi / 180)
      VpY <- Vp * cos(AngP * pi / 180)
      
      df_season_list[["Vp"]][i,]<-quantile(na.omit(values(Vp)), c(0.025, 0.5, 0.975))
      
      projection <- (VcX * VpX + VcY * VpY) / Vc
      anglediff <- abs(((AngC - AngP) + 180) %% 360 - 180)
      mismatch <- projection-Vc
      
      df_season_list[["anglediff"]][i,]<-quantile(na.omit(values(anglediff)), c(0.025, 0.5, 0.975))
      
      df_season_list[["mismatch"]][i,]<-quantile(na.omit(values(mismatch)), c(0.025, 0.5, 0.975))
    }
    
    p_compare<-vector(mode="list")
    x_label_list<-list(c(expression(bolditalic(v)["MAT"]),
                         expression(bolditalic(v)["MST"]),
                         expression(bolditalic(v)["MFT"])),
                       c(expression(bolditalic(v)["GSL"]),
                         expression(bolditalic(v)["SOS"]),
                         expression(bolditalic(v)["EOS"])),
                       c(expression(italic(θ)["GSL,MAT"]),
                         expression(italic(θ)["SOS,MST"]),
                         expression(italic(θ)["EOS,MFT"])),
                       c(expression(italic(δ)["GSL,MAT"]),
                         expression(italic(δ)["SOS,MST"]),
                         expression(italic(δ)["EOS,MFT"])))
    y_label_list <- c(expression(bolditalic(v)* " (km yr"^-1 * ")"), 
                      expression(bolditalic(v)* " (km yr"^-1 * ")"),
                      expression(italic(θ)* " (°)"),
                      expression(italic(δ)* " (km yr"^-1 * ")")
    )
    
    for ( i in 1:4) {
      p_compare[[i]]<-
        ggplot(df_season_list[[i]])+
        geom_point(aes(y=median, x=1:3), col="red", cex=2)+
        geom_errorbar(aes(x=1:3, ymin=lower, ymax=upper), width=0)+
        geom_text(aes(x=1:3, y=median, label=round(median,2)), nudge_x = 0.3, cex=3.5)+
        geom_text(aes(x=1:3, y=lower, label=round(lower,2)), nudge_x = 0.3, cex=3.5, alpha=0.75)+
        geom_text(aes(x=1:3, y=upper, label=round(upper,2)), nudge_x = 0.3, cex=3.5, alpha=0.75)+
        scale_x_continuous(breaks=1:3, labels = x_label_list[[i]], limits = c(1,3.5))+
        xlab("")+
        ylab(y_label_list[i])+
        theme_classic()
    }
    
      cairo_pdf("./figures/Figure 2.pdf", width = 10, height = 9)  
      schematic<-get_schematic(season)
    grid.arrange(annotate_figure(p_C + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[1]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
                 annotate_figure(p_P + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[2]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
                 annotate_figure(p_anglediff + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[3]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
                 annotate_figure(p_mismatch + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[4]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
                 p_compare[[1]],
                 p_compare[[2]],
                 p_compare[[3]],
                 p_compare[[4]],
                 widths = c(5, 1.5),
                 heights = c(1, 1, 1, 1),
                 layout_matrix = rbind(
                   c(1, 5),
                   c(2, 6),
                   c(3, 7),
                   c(4, 8)
                 )
    )
    dev.off()
  }
  if (season!="year") {
    if (season=="spring") {
      cairo_pdf("./figures/Figure S9.pdf", width = 8, height = 9)
    }
    if (season=="fall") {
      cairo_pdf("./figures/Figure S10.pdf", width = 8, height = 9)
    }
    if (season=="year_prcp") {
      cairo_pdf("./figures/Figure S12.pdf", width = 8, height = 9)
    }
  }
  schematic<-get_schematic(season)
  grid.arrange(annotate_figure(p_C + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[1]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
               annotate_figure(p_P + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[2]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
               annotate_figure(p_anglediff + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[3]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
               annotate_figure(p_mismatch + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[4]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
               ncol=1
  )
  dev.off()
}

##### MAT and GSL correlation (Figure S1)
file <- c("./data/MAT.nc")
C <- brick(file)
C <- C[[81:114]]

file <- c("./data/LOS_World_0.5deg_QA.nc")
P <- brick(file)

df_C<-as.data.frame(C, xy=T) 
colnames(df_C)<-c("x","y",1:34)
df_C<-df_C %>% 
  gather(key="year", value="C", -x, -y)

df_P<-as.data.frame(P, xy=T)
colnames(df_P)<-c("x","y",1:34)
df_P<-df_P %>% 
  gather(key="year", value="P", -x, -y)

df_combined<-left_join(df_C, df_P, by=c("x","y","year")) %>% 
  drop_na()

# regression
summary(lm(data=df_combined, P~C+I(C^2)))

cairo_pdf("./figures/Figure S1.pdf")
ggplot(df_combined)+
  stat_binhex(aes(x=C, y=P),bins=100)+
  # geom_point(aes(x=C, y=P), alpha=0.2)+
  geom_smooth(aes(x=C, y=P),method ="lm" ,formula = y ~ x + I(x^2), col="red")+
  xlab("Mean annual temperature (°C)")+
  ylab("Growing season length (day)")+
  theme_classic()+
  ylim(c(min(df_combined$P), max(df_combined$P)))+
  scale_fill_viridis_c()
dev.off()

##### Choice of spatial scale (Table S1)
### 0.5 deg
file <- c("./data/MAT.nc")
C <- brick(file)
C <- C[[81:114]]
C<-crop(C,e)*NOSmask0.5*QAmask0.5

vt <- tempTrend(C, th = 10)
vg <- spatGrad(C, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 1.026765
median(values(vg[[1]]),na.rm=T)*111.325/2*3
# [1] 1.937107

### 1 deg
C<-aggregate(C, fact=2)
vt <- tempTrend(C, th = 10)
vg <- spatGrad(C, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 1.074569
median(values(vg[[1]]),na.rm=T)*111.325*3
# [1] 3.250164

### 2 deg
C<-aggregate(C, fact=2)
vt <- tempTrend(C, th = 10)
vg <- spatGrad(C, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 1.088881
median(values(vg[[1]]),na.rm=T)*111.325*2*3
# [1] 5.53448

### 0.5 deg
file <- c("./Data/Phenology/LOS_World_0.5deg_QA.nc")
P <- brick(file)
P<-crop(P,e)*NOSmask0.5*QAmask0.5

vt <- tempTrend(P, th = 10)
vg <- spatGrad(P, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 30.59171
median(values(vg[[1]]),na.rm=T)*111.325/2*3
# [1] 28.8513

### 1 deg
P<-aggregate(P, fact=2)
vt <- tempTrend(P, th = 10)
vg <- spatGrad(P, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 30.47709
median(values(vg[[1]]),na.rm=T)*111.325*3
# [1] 37.1471

### 2 deg
P<-aggregate(P, fact=2)
vt <- tempTrend(P, th = 10)
vg <- spatGrad(P, th = -Inf, projected = FALSE)

median(values(vt[[1]]),na.rm=T)*33
# [1] 30.39347
median(values(vg[[1]]),na.rm=T)*111.325*2*3
# [1] 49.39891

##### Accounting for possible improper space-for-time substitution (Figure S2)
### Calculate climate-phenology metrics
file <- c("./data/MAT.nc")
C <- brick(file)
C <- C[[81:114]]
vt <- tempTrend(C, th = 10)
vg <- spatGrad(C, th = -Inf, projected = FALSE)
gv <- gVoCC(vt, vg)
gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
Vc <- abs(gv[[1]])
AngC <- gv[[2]]
gradangleC<-vg[[2]]
VcX <- Vc * sin(AngC * pi / 180)
VcY <- Vc * cos(AngC * pi / 180)

file <- c("./data/LOS_World_0.5deg_QA.nc")
P <- brick(file)
vt <- tempTrend(P, th = 10)
vg <- spatGrad(P, th = -Inf, projected = FALSE)
gv <- gVoCC(vt, vg)
gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
Vp <- abs(gv[[1]])
AngP <- gv[[2]]
gradangleP<-vg[[2]]
VpX <- Vp * sin(AngP * pi / 180)
VpY <- Vp * cos(AngP * pi / 180)

projection <- (VcX * VpX + VcY * VpY) / Vc
anglediff <- abs(((AngC - AngP) + 180) %% 360 - 180)
mismatch <- projection-Vc

### Use directions of spatial gradients to create mask
gradanglediff <- abs(((gradangleC - gradangleP) + 180) %% 360 - 180)
gradanglediff<-crop(gradanglediff, e)
sum(na.omit(values(gradanglediff)) > 90) / length(na.omit(values(gradanglediff)))
sum(na.omit(values(gradanglediff)) < 90) / length(na.omit(values(gradanglediff)))
anglemask <- gradanglediff
anglemask[anglemask > 90] <- NA
anglemask[!is.na(anglemask)] <- 1

Vc<-Vc*anglemask
Vp<-Vp*anglemask
anglediff<-anglediff*anglemask
mismatch<-mismatch*anglemask

### Summarize
quantile(na.omit(values(Vc)), c(0.025, 0.5, 0.975))
quantile(na.omit(values(Vp)), c(0.025, 0.5, 0.975))
quantile(na.omit(values(anglediff)), c(0.025, 0.5, 0.975))
quantile(na.omit(values(mismatch)), c(0.025, 0.5, 0.975))

### Plot
cutpts <- quantile(na.omit(values(stack(Vc,Vp))), probs = seq(0, 1, length.out = 21))
fieldC <- stack(Vc / Vc, AngC)
names(fieldC) <- c("slope", "aspect")
bg_Vc <- levelplot(Vc, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = rev(colorspace::heat_hcl(20)), main = "", 
                   xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 1, 2, 5, 10), at = labelpts(stack(Vc,Vp), c(0, 1, 2, 5, 10)))))
arrow_Vc <- vectorplot(fieldC, isField = TRUE, unit = "degrees", col.arrows = "black", narrows = 500, aspX = 4, aspY = 4, lwd.arrows = 1, alpha = 1, reverse = FALSE, colorkey = FALSE, region = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(draw = FALSE))
p_C <- bg_Vc + arrow_Vc

fieldP <- stack(Vp / Vp, AngP)
names(fieldP) <- c("slope", "aspect")
bg_Vp <- levelplot(Vp, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = rev(colorspace::heat_hcl(20)), main = "", 
                   xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 1, 2, 5, 10), at = labelpts(stack(Vc,Vp), c(0, 1, 2, 5, 10)))))
arrow_Vp <- vectorplot(fieldP, isField = TRUE, unit = "degrees", col.arrows = "black", narrows = 500, aspX = 4, aspY = 4, lwd.arrows = 1, alpha = 1, reverse = FALSE, colorkey = FALSE, region = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(draw = FALSE))
p_P <- bg_Vp + arrow_Vp

col_pal_diverge <- c(colorspace::sequential_hcl(17, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(3, h = 330, l = c(40, 90), power = 1))) #anglemask
cutpts <- quantile(na.omit(values(anglediff)), probs = seq(0, 1, length.out = 21))
p_anglediff <- levelplot(anglediff, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = col_pal_diverge, main = "", 
                         xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(0, 45, 90, 135, 180), at = labelpts(anglediff, c(0, 45, 90, 135, 180)))))
# p_anglediff + latticeExtra::layer(sp.polygons(land))

col_pal_diverge <- rev(c(colorspace::sequential_hcl(12, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(8, h = 330, l = c(40, 90), power = 1)))) #anglemask
cutpts <- quantile(na.omit(values(mismatch)), probs = seq(0, 1, length.out = 21))
p_mismatch <- levelplot(mismatch, par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = col_pal_diverge, main = "", 
                        xlab = list(cex = 0),ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.05), labels = list(labels = c(-10, -5, -2, -1, 0, 1, 2, 5, 10), at = labelpts(mismatch, c(-10, -5, -2, -1, 0, 1, 2, 5, 10)))))
# p_mismatch + latticeExtra::layer(sp.polygons(land))

cairo_pdf("./figures/Figure S2.pdf", width = 8, height = 9)
grid.arrange(annotate_figure(p_C + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[1]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_P + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[2]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_anglediff + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[3]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_mismatch + latticeExtra::layer(sp.polygons(land))+ latticeExtra::layer(grid.raster(schematic[[4]],x=0.025,y=0.3,height=0.6,just="left",interpolate=F),under=FALSE), fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
            ncol=1
)
dev.off()

# continue to code 3 for spatial regression (Table S2)

### Summary statistics for different study areas (Figure S8)
file <- c("./data/MAT.nc")
C <- brick(file)
C <- C[[81:114]]
vt <- tempTrend(C, th = 10)
vg <- spatGrad(C, th = -Inf, projected = FALSE)
gv <- gVoCC(vt, vg)
gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
Vc <- abs(gv[[1]])
AngC <- gv[[2]]
VcX <- Vc * sin(AngC * pi / 180)
VcY <- Vc * cos(AngC * pi / 180)

file <- c("./data/LOS_World_0.5deg_QA.nc")
P <- brick(file)
vt <- tempTrend(P, th = 10)
vg <- spatGrad(P, th = -Inf, projected = FALSE)
gv <- gVoCC(vt, vg)
gv <- crop(gv, e) * NOSmask0.5 * QAmask0.5
Vp <- abs(gv[[1]])
AngP <- gv[[2]]
VpX <- Vp * sin(AngP * pi / 180)
VpY <- Vp * cos(AngP * pi / 180)

projection <- (VcX * VpX + VcY * VpY) / Vc
anglediff <- abs(((AngC - AngP) + 180) %% 360 - 180)
mismatch <- projection-Vc

df_area<-data.frame(lower=rep(NA, 6), median=rep(NA, 6), upper=rep(NA, 6))
df_area_list<-vector(mode="list")
for (i in c("Vc", "Vp", "anglediff", "mismatch")) {
  df_area_list[[i]]<-df_area
}
crop_layer<-function (i, var) {
  if (var=="Vc") {layer<-Vc}
  if (var=="Vp") {layer<-Vp}
  if (var=="anglediff") {layer<-anglediff}
  if (var=="mismatch") {layer<-mismatch}
  if(i==1) {temp<-layer}
  if (i==2) {temp<-crop(layer, extent(-180, 180, 50, 70))}
  if (i==3) {temp<-crop(layer, extent(-180, 180, 30, 50))}
  if (i==4) {temp<-layer * Eumask}
  if (i==5) {temp<-layer * NAmemask}
  if (i==6) {temp<-layer * Asmask}
  return(temp)
}
for (i in 1:6) {
  for (var in c("Vc", "Vp", "anglediff", "mismatch")) {
    temp<-crop_layer(i, var)
    df_area_list[[var]][i,]<-quantile(na.omit(values(temp)), c(0.025, 0.5, 0.975))
  }
}

label_list<-c("Total Area Studied", "50°–70° N", "30°–50° N", "Europe", "North America", "Asia")
y_label_list <- c(expression(bolditalic(v)["MAT"]* " (km yr"^-1 * ")"), 
                  expression(bolditalic(v)["GSL"]* " (km yr"^-1 * ")"),
                  expression(italic(θ)["GSL,MAT"]* " (°)"),
                  expression(italic(δ)["GSL,MAT"]* " (km yr"^-1 * ")")
)
p_compare<-vector(mode="list")
for ( i in 1:4) {
  p_compare[[i]]<-
    ggplot(df_area_list[[i]])+
    geom_point(aes(y=median, x=6:1), col="red", cex=2)+
    geom_errorbar(aes(x=6:1, ymin=lower, ymax=upper), width=0)+
    geom_text(aes(x=6:1, y=median, label=round(median,2)), nudge_x = 0.3, cex=3.5)+
    geom_text(aes(x=6:1, y=lower, label=round(lower,2)), nudge_x = 0.3, cex=3.5, alpha=0.75)+
    geom_text(aes(x=6:1, y=upper, label=round(upper,2)), nudge_x = 0.3, cex=3.5, alpha=0.75)+
    scale_x_continuous(breaks=6:1, labels = (label_list), limits = c(1,6.5))+
    xlab("")+
    ylab(y_label_list[i])+
    theme_classic()+
    coord_flip()
}

cairo_pdf("./figures/Figure S8.pdf", width = 6, height = 8)
grid.arrange(annotate_figure(p_compare[[1]], fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_compare[[2]], fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_compare[[3]], fig.lab = "c", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_compare[[4]], fig.lab = "d", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             ncol=1
)
dev.off()
