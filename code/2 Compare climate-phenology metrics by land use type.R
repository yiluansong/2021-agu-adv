# Code for  the first part of 2.1.3 (leading to Figure 3 and supplementary figures)

source("./tools/packages.R")
source("./tools/load masks.R")
source("./tools/graphing settings.R")

##### Compare climate-phenology metrics by land use type (Figure 3)
### Get climate-phenology metrics
file <- c("./data/MAT.nc") # mean annual temperature
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

file <- c("./data/LOS_World_0.5deg_QA.nc") # length of season
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

### Get land use type
file <- c("./data/anthro2_a2000.nc")
cover <- raster(file)
cover <- crop(cover,e)
cover[cover<30]<-1 #settlements
cover[cover>30&cover<40]<-2 #croplands
cover[cover>40&cover<50]<-3 #rangelands
cover[cover>50&cover<60]<-4 #semilatural
cover[cover>60]<-5 #wild

QAmask_1_6 <- resample(QAmask0.5, cover, method = "ngb")
NOSmask_1_6 <- resample(NOSmask0.5, cover, method = "ngb")

Vc_1_6 <- resample(Vc, cover, method = "ngb")
Vp_1_6 <- resample(Vp, cover, method = "ngb")
projection_1_6 <- resample(projection, cover, method = "ngb")
mismatch_1_6 <- resample(mismatch, cover, method = "ngb")
anglediff_1_6 <- resample(anglediff, cover, method = "ngb")

combined <- brick(Vc_1_6, Vp_1_6, projection_1_6, mismatch_1_6, anglediff_1_6, cover)
names(combined) <- c("Vc", "Vp", "Proj", "Mis", "Angdiff", "Cover")
combineddf <- as.data.frame(combined)
combineddf <- na.omit(combineddf)

cover_name <- c("Settlements", "Croplands", "Rangelands", "Semi-natural", "Wildlands")
Cover <- c(1,2,3,4,5)
cover_key <- data.frame(Cover, cover_name)

combineddf <- merge(x = combineddf, y = cover_key, by = "Cover", all = TRUE)

### Summarize
combineddf %>%
  group_by(cover_name) %>%
  summarize(
    lower = quantile(Vc, probs = 0.025),
    median = quantile(Vc, probs = 0.5),
    upper = quantile(Vc, probs = 0.975)
  ) %>%
  arrange(desc(median))

combineddf %>%
  group_by(cover_name) %>%
  dplyr::summarize(
    lower = quantile(Vp, probs = 0.025),
    median = quantile(Vp, probs = 0.5),
    upper = quantile(Vp, probs = 0.975)
  ) %>%
  arrange(desc(median))

combineddf %>%
  group_by(cover_name) %>%
  dplyr::summarize(
    lower = quantile(Angdiff, probs = 0.025),
    median = quantile(Angdiff, probs = 0.5),
    upper = quantile(Angdiff, probs = 0.975)
  ) %>%
  arrange(desc(median))

combineddf %>%
  group_by(cover_name) %>%
  dplyr::summarize(
    lower = quantile(Mis, probs = 0.025),
    median = quantile(Mis, probs = 0.5),
    upper = quantile(Mis, probs = 0.975)
  ) %>%
  arrange(desc(median))

# percentage of settlements
(combineddf %>% filter(cover_name=="Settlements") %>% nrow())/(combineddf  %>% nrow())

# test difference in distribution
anova_res<-aov(lm( Vc ~ cover_name, data=combineddf ))
tukey_res <- TukeyHSD(x=anova_res , conf.level=0.95)
tukey_res

anova_res<-aov(lm( Vp ~ cover_name, data=combineddf ))
tukey_res <- TukeyHSD(x=anova_res , conf.level=0.95)
tukey_res

anova_res<-aov(lm( Angdiff ~ cover_name, data=combineddf ))
tukey_res <- TukeyHSD(x=anova_res , conf.level=0.95)
tukey_res

anova_res<-aov(lm(Mis ~ cover_name, data=combineddf ))
tukey_res <- TukeyHSD(x=anova_res , conf.level=0.95)
tukey_res

### Plot
cover<-resample(cover * QAmask_1_6 * NOSmask_1_6,disaggregate(Vc,fact=6),method="ngb")

p_cover <-
  levelplot(cover,
            par.settings = my.settings,
            # main="Land Cover Type",
            at = 0.5:5.5,
            col.regions = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9","#2c7bb6"), #http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=5
            colorkey = list(space = "bottom", width = 0.8, height = 1, at = 0.5:5.5,
                            col = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9","#2c7bb6"),
                            labels = list(cex = 1, at = rev(1:5), labels = rev(cover_name))
            ),
            margin = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = c(30,50,70)), x = list(at = seq(-180, 180, 60)))
  )
# p_cover + latticeExtra::layer(sp.polygons(land))

plot_group <- function(variable, group) {
  if (group == "Settlements") color<-"#d7191c"
  if (group == "Croplands") color<-"#fdae61"
  if (group == "Rangelands") color<-"#ffffbf"
  if (group == "Semi-natural") color<-"#abd9e9"
  if (group == "Wildlands") color<-"#2c7bb6"
  if (variable == "Vc") {
    xrange <- c(quantile(combineddf$Vc, probs = 0.025), quantile(combineddf$Vc, probs = 0.975))
    yrange <- c(0, 0.4)
    ylabel <- c(0, 4)
    yrange2 <- ylabel / (quantile(combineddf$Vc, probs = 0.975) - quantile(combineddf$Vc, probs = 0.025))
    median <- quantile(combineddf[combineddf$cover_name == group, ]$Vc, 0.5)
    lower <- quantile(combineddf[combineddf$cover_name == group, ]$Vc, 0.025)
    upper <- quantile(combineddf[combineddf$cover_name == group, ]$Vc, 0.975)
    p <- ggplot(data = combineddf[combineddf$cover_name == group & combineddf$Vc >= quantile(combineddf$Vc, probs = 0.025) & combineddf$Vc <= quantile(combineddf$Vc, probs = 0.975), ], aes(x = Vc, y = stat(density)))
  }
  if (variable == "Vp") {
    xrange <- c(quantile(combineddf$Vp, probs = 0.025), quantile(combineddf$Vp, probs = 0.975))
    yrange <- c(0, 0.3)
    ylabel <- c(0, 8)
    yrange2 <- ylabel / (quantile(combineddf$Vp, probs = 0.975) - quantile(combineddf$Vp, probs = 0.025))
    median <- quantile(combineddf[combineddf$cover_name == group, ]$Vp, 0.5)
    lower <- quantile(combineddf[combineddf$cover_name == group, ]$Vp, 0.025)
    upper <- quantile(combineddf[combineddf$cover_name == group, ]$Vp, 0.975)
    p <- ggplot(data = combineddf[combineddf$cover_name == group & combineddf$Vp >= quantile(combineddf$Vp, probs = 0.025) & combineddf$Vp <= quantile(combineddf$Vp, probs = 0.975), ], aes(x = Vp, y = stat(density)))
  }
  if (variable == "Angdiff") {
    xrange <- c(0, 180)
    yrange <- c(0, 0.018)
    ylabel <- c(0, 2)
    yrange2 <- ylabel / (180 - 0)
    median <- quantile(combineddf[combineddf$cover_name == group, ]$Angdiff, 0.5)
    lower <- quantile(combineddf[combineddf$cover_name == group, ]$Angdiff, 0.025)
    upper <- quantile(combineddf[combineddf$cover_name == group, ]$Angdiff, 0.975)
    p <- ggplot(data = combineddf[combineddf$cover_name == group, ], aes(x = Angdiff, y = stat(density)))
  }
  if (variable == "Mis") {
    xrange <- c(quantile(combineddf$Mis, probs = 0.025), quantile(combineddf$Mis, probs = 0.975))
    yrange <- c(0, 0.18)
    ylabel <- c(0, 6)
    yrange2 <- ylabel / (quantile(combineddf$Mis, probs = 0.975) - quantile(combineddf$Mis, probs = 0.025))
    median <- quantile(combineddf[combineddf$cover_name == group, ]$Mis, 0.5)
    lower <- quantile(combineddf[combineddf$cover_name == group, ]$Mis, 0.025)
    upper <- quantile(combineddf[combineddf$cover_name == group, ]$Mis, 0.975)
    p <- ggplot(data = combineddf[combineddf$cover_name == group & combineddf$Mis >= quantile(combineddf$Mis, probs = 0.025) & combineddf$Mis <= quantile(combineddf$Mis, probs = 0.975), ], aes(x = Mis, y = stat(density)))
  }
  summ<-paste0(round(median, 2), " (", round(lower, 2), ", ", round(upper, 2), ")" )
  
  p +
    geom_histogram(bins = 100, fill = color) +
    xlim(xrange[1], xrange[2]) +
    scale_y_continuous(limits = c(yrange[1], yrange[2]), breaks = seq(yrange2[1], yrange2[2], length.out = 3), labels = format(seq(ylabel[1], ylabel[2], length.out = 3), nsmall = 0)) +
    # ylim(yrange[1],yrange[2])+
    xlab("") +
    ylab("") +
    annotate(geom = "text", x = xrange[2] - (xrange[2] - xrange[1]) * 0.3, y = yrange[2] * 0.7, label = summ, size = 3) +
    geom_vline(xintercept = median) +
    scale_fill_discrete(guide = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_blank())
}

p_Vc1 <- plot_group("Vc", "Settlements")
p_Vc2 <- plot_group("Vc", "Croplands")
p_Vc3 <- plot_group("Vc", "Rangelands")
p_Vc4 <- plot_group("Vc", "Semi-natural")
p_Vc5 <- plot_group("Vc", "Wildlands")

p_Vp1 <- plot_group("Vp", "Settlements")
p_Vp2 <- plot_group("Vp", "Croplands")
p_Vp3 <- plot_group("Vp", "Rangelands")
p_Vp4 <- plot_group("Vp", "Semi-natural")
p_Vp5 <- plot_group("Vp", "Wildlands")

p_Angdiff1 <- plot_group("Angdiff", "Settlements")
p_Angdiff2 <- plot_group("Angdiff", "Croplands")
p_Angdiff3 <- plot_group("Angdiff", "Rangelands")
p_Angdiff4 <- plot_group("Angdiff", "Semi-natural")
p_Angdiff5 <- plot_group("Angdiff", "Wildlands")

p_Mis1 <- plot_group("Mis", "Settlements")
p_Mis2 <- plot_group("Mis", "Croplands")
p_Mis3 <- plot_group("Mis", "Rangelands")
p_Mis4 <- plot_group("Mis", "Semi-natural")
p_Mis5 <- plot_group("Mis", "Wildlands")

p_Vc_axis <-
  ggplot(combineddf, aes(x = Vc, y = stat(density))) +
  geom_histogram(bins = 100, fill = "white") +
  ylab("") +
  xlab(expression(bolditalic(v)[MAT]* " (km yr"^-1 * ")")) +
  xlim(quantile(combineddf$Vc, probs = 0.025), quantile(combineddf$Vc, probs = 0.975)) +
  scale_y_continuous(limits = c(0, 0), breaks = 0, labels = format(0, nsmall = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(color = "white")) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank())+
  theme(axis.title.x = element_text(size = 10))

p_Vp_axis <-
  ggplot(combineddf, aes(x = Vp, y = stat(density))) +
  geom_histogram(bins = 100, fill = "white") +
  ylab("") +
  xlab(expression(bolditalic(v)[GSL]* " (km yr"^-1 * ")")) +
  xlim(quantile(combineddf$Vp, probs = 0.025), quantile(combineddf$Vp, probs = 0.975)) +
  scale_y_continuous(limits = c(0, 0), breaks = 0, labels = format(0, nsmall = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(color = "white")) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank())+
  theme(axis.title.x = element_text(size = 10))

p_Angdiff_axis <-
  ggplot(combineddf, aes(x = Angdiff, y = stat(density))) +
  geom_histogram(bins = 100, fill = "white") +
  ylab("") +
  xlab(expression(italic(θ)["GSL, MAT"]* " (°)")) +
  xlim(0, 180) +
  scale_y_continuous(limits = c(0, 0), breaks = 0, labels = format(0, nsmall = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(color = "white")) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank())+
  theme(axis.title.x = element_text(size = 10))

p_Mis_axis <-
  ggplot(combineddf, aes(x = Mis, y = stat(density))) +
  geom_histogram(bins = 100, fill = "white") +
  ylab("") +
  xlab(expression("lagging "*italic(δ)["GSL, MAT"]* " (km yr"^-1 * ")"*" leading")) +
  xlim(quantile(combineddf$Mis, probs = 0.025), quantile(combineddf$Mis, probs = 0.975)) +
  scale_y_continuous(limits = c(0, 0), breaks = 0, labels = format(0, nsmall = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(color = "white")) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank())+
  theme(axis.title.x = element_text(size = 10))

p_y_axis <-
  ggplot(combineddf, aes(x = Vc, y = stat(density))) +
  geom_histogram(bins = 100,fill="white") +
  ylab("Percentage of pixels") +
  xlab("") +
  xlim(quantile(combineddf$Vc, probs = 0.025), quantile(combineddf$Vc, probs = 0.975)) +
  ylim(0, 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.1)))

p_void <- ggplot() + theme_void()

cairo_pdf("./figures/Figure 3.pdf", width = 10, height = 7)
grid.arrange(
  grobs = list(
    annotate_figure(p_cover + latticeExtra::layer(sp.polygons(land)), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
    # p_Vc1, p_Vp1, p_Angdiff1, p_Mis1,
    p_Vc2, p_Vp2, p_Angdiff2, p_Mis2,
    p_Vc3, p_Vp3, p_Angdiff3, p_Mis3,
    p_Vc4, p_Vp4, p_Angdiff4, p_Mis4,
    p_Vc5, p_Vp5, p_Angdiff5, p_Mis5,
    p_Vc_axis, p_Vp_axis, p_Angdiff_axis, p_Mis_axis,
    p_y_axis,
    annotate_figure(p_void, fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold")
  ),
  widths = c(0.09, 1, 1, 1, 1),
  heights = c(1.3,  0.5, 0.5, 0.5, 0.5, 0.2),
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1),
    # c(27, 2, 3, 4, 5),
    c(27, 6, 7, 8, 9),
    c(26, 10, 11, 12, 13),
    c(26, 14, 15, 16, 17),
    c(26, 18, 19, 20, 21),
    c(NA, 22, 23, 24, 25)
  )
)
dev.off()

##### Potential limitations on the estimation of spatial gradients (Figure S3, S4)
file <- c("./data/MAT.nc") # mean annual temperature
C <- brick(file)
C <- C[[81:114]]
vg <- spatGrad(C, th = -Inf, projected = FALSE)
grad_C<-data.frame(gradient=as.data.frame(vg$Grad)$Grad, variable="MAT") %>% drop_na()

file <- c("./data/LOS_World_0.5deg_QA.nc") # length of season
P <- brick(file)
vg <- spatGrad(P, th = -Inf, projected = FALSE)
grad_P<-data.frame(gradient=as.data.frame(vg$Grad)$Grad, variable="GSL") %>% drop_na()

grad_df<-bind_rows(grad_C, grad_P)

p_gradient_cp<-
  ggplot(grad_df)+
  geom_density(aes(x=gradient, fill=variable), alpha=0.5)+
  xlim(0,1)+
  xlab(expression("spatial gradient of focal variable (day km"^-1 * ")")) +
  theme_classic()

cairo_pdf("./figures/Figure S3.pdf", width = 6, height = 6)
print(p_gradient_cp)
dev.off()

vg <- spatGrad(P, th = -Inf, projected = FALSE)
vg_1_6 <- resample(vg, cover, method = "ngb")
combined <- brick(vg_1_6$Grad, cover)
names(combined) <- c("gradient", "Cover")
combineddf <- as.data.frame(combined)
combineddf <- na.omit(combineddf)

cover_name <- c("Settlements", "Croplands", "Rangelands", "Semi-natural", "Wildlands")
Cover <- c(1,2,3,4,5)
cover_key <- data.frame(Cover, cover_name)

combineddf <- merge(x = combineddf, y = cover_key, by = "Cover", all = TRUE)

p_gradient_hist<-
  ggplot(combineddf)+
  geom_density(aes(x=gradient, fill=cover_name), alpha=0.5)+
  xlim(0,1)+
  xlab(expression("spatial gradient of growing season length (day km"^-1 * ")")) +
  guides(fill=guide_legend("land cover type"))+
  theme_classic()

cairo_pdf("./figures/Figure S4.pdf", width = 6, height = 6)
print(p_gradient_hist)
dev.off()

##### Plant functional type (Figure S11)
file <- c("./data/MODIS land cover.tif")
# https://lpdaac.usgs.gov/products/mcd12c1v006/
# https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pdf
pft <- raster(file)
pft<-crop(pft, e)
pft<-aggregate(pft,fact=10, fun=modal)

pft[pft==0]<-NA
pft[pft==255]<-NA
pft_smooth <- focal(pft, w=matrix(1, 5, 5), fun=modal,na.rm=TRUE,pad=TRUE)
pft_polygon<-rasterToPolygons(pft_smooth, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
area<-gUnaryUnion(pft_polygon)

pft_legend<-c("Evergreen Needleleaf Forest", "Evergreen Broadleaf Forest", "Deciduous Needleleaf Forest", "Deciduous Broadleaf Forest", "Mixed Forest", "Closed Shrubland", "Open Shrubland", "Woody Savannas", "Savannas", "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built-Up", "Cropland/Natural Vegetation", "Snow and Ice", "Barren or Sparsely Vegetated")

cols<-c("#00678A", "#FFFFFF", "#FFFFFF", "#56641A", "#C0AFFB", "#FFFFFF", "#5ECCAB", "#FFFFFF", "#FFFFFF", "#E6A176", "#FFFFFF", "#984464", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")

p_pft <-
  levelplot(pft_smooth,
            par.settings = my.settings,
            at = 0.5:16.5,
            col.regions = cols,
            colorkey=NULL,
            # colorkey = list(space = "bottom", width = 0.8, height = 1, at = 0.5:16.5,
            #                 col = cols,
            #                 labels = list(cex = 1, at = rev(1:16), labels = rev(pft_legend))
            # ),
            margin = FALSE, xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = c(30,50,70)), x = list(at = seq(-180, 180, 60)))
  )
# p_pft + latticeExtra::layer(sp.polygons(land))

combined<-brick(mismatch, pft_smooth)
names(combined)<-c("mismatch","pft")
combineddf<-as.data.frame(combined)
combineddf<-na.omit(combineddf)

pft_name<-c("Water","Evergreen Needleleaf Forest", "Evergreen Broadleaf Forest", "Deciduous Needleleaf Forest", "Deciduous Broadleaf Forest", "Mixed Forest", "Closed Shrubland", "Open Shrubland", "Woody Savannas", "Savannas", "Grasslands", "Permanent Wetlands", "Croplands", "Urban and Built-Up", "Cropland/Natural Vegetation", "Snow and Ice", "Barren or Sparsely Vegetated","Unclassified")
pft_index<-c(seq(0,16,1), 255)
pft_key<-data.frame(pft_index,pft_name)
pft_key$pft_name<-as.character(pft_key$pft_name)
for(i in 1:NROW(combineddf)) {
  combineddf$pft_name[i]<-pft_key[pft_key$pft_index==combineddf$pft[i],2]
}

cols <- c("Croplands" = "#984464", "Deciduous Broadleaf Forest" = "#56641A", "Evergreen Needleleaf Forest" = "#00678A", "Grasslands" = "#E6A176", "Mixed Forest" = "#C0AFFB","Open Shrubland" = "#5ECCAB")

p_pft_hist<-
  ggplot(combineddf %>% filter(pft_name %in% names(cols)))+
  geom_density(aes(x=mismatch, fill=pft_name), alpha=0.75)+
  scale_fill_manual(values=cols)+
  xlim(quantile(combineddf$mismatch, probs = 0.025),quantile(combineddf$mismatch, probs = 0.975))+
  xlab(expression("lagging      "*italic(δ)["GSL, MAT"]* " (km yr"^-1 * ")"*"      leading")) +
  guides(fill=guide_legend("plant functional type"))+
  theme_classic()

cairo_pdf("./figures/Figure S11.pdf", width = 10, height = 5)
grid.arrange(annotate_figure(p_pft + latticeExtra::layer(sp.polygons(land)), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             annotate_figure(p_pft_hist, fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
             ncol=1
)
dev.off()
