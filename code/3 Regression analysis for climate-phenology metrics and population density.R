# Code for  the second part of 2.1.3 (leading to Figure 4 and supplementary figures)

source("./code/setup/packages.R")
source("./code/setup/masks.R")
source("./code/setup/graphics.R")

##### Regression analysis for climate-phenology metrics and population density (Figure 4, Table S3, S4, S5, Figure S6, S7)
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
mismatch <- projection - Vc

### Get human population
file <- c("./data/gpw_v4_population_density_rev11_30_min.nc")
Human <- brick(file)[[1]]
Human <- crop(Human, e) * NOSmask0.5 * QAmask0.5

cutpts <- quantile(na.omit(values(Human)), probs = seq(0, 1, length.out = 6))
p_human <- levelplot(Human,
  par.settings = my.settings, at = cutpts, margin = FALSE, col.regions = rev(colorspace::heat_hcl(5)), main = "",
  xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8, at = seq(0, 1, by = 0.2), labels = list(labels = c(0, 0.1, 1, 10, 100), at = labelpts(Human, c(0, 0.1, 1, 10, 100))))
)
# p_human + latticeExtra::layer(sp.polygons(land))

# check human population density in each land use type
file <- c("./data/anthro2_a2000.nc")
cover <- raster(file)
cover <- crop(cover, e)
cover[cover < 30] <- 1 # settlements
cover[cover > 30 & cover < 40] <- 2 # croplands
cover[cover > 40 & cover < 50] <- 3 # rangelands
cover[cover > 50 & cover < 60] <- 4 # semilatural
cover[cover > 60] <- 5 # wild
cover_name <- c("Settlements", "Croplands", "Rangelands", "Semi-natural", "Wildlands")
Cover <- c(1, 2, 3, 4, 5)
cover_key <- data.frame(Cover, cover_name)

Human_1_6 <- resample(Human, cover, method = "ngb")
human_cover_brick <- brick(Human_1_6, cover)
human_cover_df <- as.data.frame(human_cover_brick)
colnames(human_cover_df) <- c("Human", "Cover")
human_cover_df <- merge(x = human_cover_df, y = cover_key, by = "Cover", all = TRUE)
human_cover_df$cover_name <- fct_relevel(human_cover_df$cover_name, c("Wildlands", "Seminatural", "Rangelands", "Croplands", "Settlements"))
human_cover_df$Human <- human_cover_df$Human + 1
human_cover_df <- na.omit(human_cover_df)

summary(res.aov <- aov(Human ~ cover_name, data = human_cover_df))

middle <- aggregate(human_cover_df$Human, list(human_cover_df$cover_name), median)[, 2]
ymin <- aggregate(human_cover_df$Human, list(human_cover_df$cover_name), function(x) {
  quantile(x, 0.025)
})[, 2]
ymax <- aggregate(human_cover_df$Human, list(human_cover_df$cover_name), function(x) {
  quantile(x, 0.975)
})[, 2]
lower <- aggregate(human_cover_df$Human, list(human_cover_df$cover_name), function(x) {
  quantile(x, 0.25)
})[, 2]
upper <- aggregate(human_cover_df$Human, list(human_cover_df$cover_name), function(x) {
  quantile(x, 0.75)
})[, 2]

df_boxplot <- data.frame(bin = unique(human_cover_df$cover_name), middle=log(middle), ymin=log(ymin), ymax=log(ymax), lower=log(lower), upper=log(upper))
p_human_cover <-
  ggplot(df_boxplot, aes(bin, middle)) +
  geom_point(pch=NA) +
  geom_segment(aes(x = seq(1, 5, 1), y = lower, xend = seq(1, 5, 1), yend = upper), colour = "darkblue", size = 5) +
  geom_segment(aes(x = seq(1, 5, 1), y = ymin, xend = seq(1, 5, 1), yend = ymax), colour = "black") +
  geom_segment(aes(x = seq(0.75, 4.75, 1), y = middle, xend = seq(1.25, 5.25, 1), yend = middle), colour = "black", size = 1) +
  # scale_x_discrete(labels=c("Wildlands", "Seminatural","Rangelands","Croplands","Settlements"))+
  scale_y_continuous(breaks=log(c(1,10,100,1000)), labels=c(1,10,100,1000))+
  xlab("Land use type") +
  ylab(expression("Population density (persons km"^-2 * ")")) +
  # ylab(expression(bolditalic(v)[MAT]* " (km yr"^-1 * ")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_human_cover

### Prepare dataframe for regression
projection_df <- as.data.frame(projection, xy = TRUE)
latitude_df <- data.frame(projection_df$x, projection_df$y, projection_df$y)
latitude <- rasterFromXYZ(latitude_df, crs = proj4string(projection))

Vc_df <- as.data.frame(Vc, xy = TRUE)
Vp_df <- as.data.frame(Vp, xy = TRUE)
projection_df <- as.data.frame(projection, xy = TRUE)
mismatch_df <- as.data.frame(mismatch, xy = TRUE)
anglediff_df <- as.data.frame(anglediff, xy = TRUE)
human_df <- as.data.frame(Human, xy = TRUE)
latitude_df <- as.data.frame(latitude, xy = TRUE)

df <- data.frame(human_df, Vc_df[, 3], Vp_df[, 3], projection_df[, 3], mismatch_df[, 3], anglediff_df[, 3], latitude_df[, 3])
colnames(df) <- c("x", "y", "human", "Vc", "Vp", "projection", "mismatch", "anglediff", "latitude")
df$Vcnew <- log(df$Vc)
df$Vpnew <- log(df$Vp)
df$anglediffnew <- log(df$anglediff / 180 / (1 - df$anglediff / 180)) # logistic
df$population.density <- log(df$human + 1)
df_valid <- na.omit(df)

### Plot
df_valid$bin <- cut(df_valid$human, cutpts, include.lowest = T)

middle <- aggregate(df_valid$Vc, list(df_valid$bin), median)[, 2]
ymin <- aggregate(df_valid$Vc, list(df_valid$bin), function(x) {
  quantile(x, 0.025)
})[, 2]
ymax <- aggregate(df_valid$Vc, list(df_valid$bin), function(x) {
  quantile(x, 0.975)
})[, 2]
lower <- aggregate(df_valid$Vc, list(df_valid$bin), function(x) {
  quantile(x, 0.25)
})[, 2]
upper <- aggregate(df_valid$Vc, list(df_valid$bin), function(x) {
  quantile(x, 0.75)
})[, 2]

df_boxplot <- data.frame(bin = unique(df_valid$bin), middle = log(middle), ymin = log(ymin), ymax = log(ymax), lower = log(lower), upper = log(upper))
p_Vc_human <-
  ggplot(df_boxplot, aes(bin, middle)) +
  geom_point() +
  geom_segment(aes(x = seq(1, 5, 1), y = lower, xend = seq(1, 5, 1), yend = upper), colour = rev(colorspace::heat_hcl(5)), size = 5) +
  geom_segment(aes(x = seq(1, 5, 1), y = ymin, xend = seq(1, 5, 1), yend = ymax), colour = "black") +
  geom_segment(aes(x = seq(0.75, 4.75, 1), y = middle, xend = seq(1.25, 5.25, 1), yend = middle), colour = "black", size = 1) +
  scale_x_discrete(label = c("[0,0.023]", "(0.023,0.439]", "(0.439,2.94]", "(2.94,25.5]", "(25.5,8010]")) +
  scale_y_continuous(breaks = log(c(1, 2, 5, 10)), labels = c(1, 2, 5, 10)) +
  xlab("") +
  ylab(expression(bolditalic(v)[MAT] * " (km yr"^-1 * ")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 8))

middle <- aggregate(df_valid$Vp, list(df_valid$bin), median)[, 2]
ymin <- aggregate(df_valid$Vp, list(df_valid$bin), function(x) {
  quantile(x, 0.025)
})[, 2]
ymax <- aggregate(df_valid$Vp, list(df_valid$bin), function(x) {
  quantile(x, 0.975)
})[, 2]
lower <- aggregate(df_valid$Vp, list(df_valid$bin), function(x) {
  quantile(x, 0.25)
})[, 2]
upper <- aggregate(df_valid$Vp, list(df_valid$bin), function(x) {
  quantile(x, 0.75)
})[, 2]

df_boxplot <- data.frame(bin = unique(df_valid$bin), middle = log(middle), ymin = log(ymin), ymax = log(ymax), lower = log(lower), upper = log(upper))
p_Vp_human <-
  ggplot(df_boxplot, aes(bin, middle)) +
  geom_point() +
  geom_segment(aes(x = seq(1, 5, 1), y = lower, xend = seq(1, 5, 1), yend = upper), colour = rev(colorspace::heat_hcl(5)), size = 5) +
  geom_segment(aes(x = seq(1, 5, 1), y = ymin, xend = seq(1, 5, 1), yend = ymax), colour = "black") +
  geom_segment(aes(x = seq(0.75, 4.75, 1), y = middle, xend = seq(1.25, 5.25, 1), yend = middle), colour = "black", size = 1) +
  scale_x_discrete(label = c("[0,0.023]", "(0.023,0.439]", "(0.439,2.94]", "(2.94,25.5]", "(25.5,8010]")) +
  scale_y_continuous(breaks = log(c(1, 2, 5, 10)), labels = c(1, 2, 5, 10)) +
  xlab("") +
  ylab(expression(bolditalic(v)[GSL] * " (km yr"^-1 * ")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 8))

middle <- aggregate(df_valid$anglediff, list(df_valid$bin), median)[, 2]
ymin <- aggregate(df_valid$anglediff, list(df_valid$bin), function(x) {
  quantile(x, 0.025)
})[, 2]
ymax <- aggregate(df_valid$anglediff, list(df_valid$bin), function(x) {
  quantile(x, 0.975)
})[, 2]
lower <- aggregate(df_valid$anglediff, list(df_valid$bin), function(x) {
  quantile(x, 0.25)
})[, 2]
upper <- aggregate(df_valid$anglediff, list(df_valid$bin), function(x) {
  quantile(x, 0.75)
})[, 2]

df_boxplot <- data.frame(bin = unique(df_valid$bin), middle = log(middle / (180 - middle)), ymin = log(ymin / (180 - ymin)), ymax = log(ymax / (180 - ymax)), lower = log(lower / (180 - lower)), upper = log(upper / (180 - upper)))
p_angdiff_human <-
  ggplot(df_boxplot, aes(bin, middle)) +
  geom_point() +
  geom_segment(aes(x = seq(1, 5, 1), y = lower, xend = seq(1, 5, 1), yend = upper), colour = rev(colorspace::heat_hcl(5)), size = 5) +
  geom_segment(aes(x = seq(1, 5, 1), y = ymin, xend = seq(1, 5, 1), yend = ymax), colour = "black") +
  geom_segment(aes(x = seq(0.75, 4.75, 1), y = middle, xend = seq(1.25, 5.25, 1), yend = middle), colour = "black", size = 1) +
  scale_x_discrete(label = c("[0,0.023]", "(0.023,0.439]", "(0.439,2.94]", "(2.94,25.5]", "(25.5,8010]")) +
  scale_y_continuous(breaks = log(c(10, 45, 90, 135, 170) / (180 - c(10, 45, 90, 135, 170))), labels = c(10, 45, 90, 135, 170)) +
  xlab("") +
  ylab(expression(italic(θ)["GSL, MAT"] * " (°)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 8))

middle <- aggregate(df_valid$mismatch, list(df_valid$bin), median)[, 2]
ymin <- aggregate(df_valid$mismatch, list(df_valid$bin), function(x) {
  quantile(x, 0.025)
})[, 2]
ymax <- aggregate(df_valid$mismatch, list(df_valid$bin), function(x) {
  quantile(x, 0.975)
})[, 2]
lower <- aggregate(df_valid$mismatch, list(df_valid$bin), function(x) {
  quantile(x, 0.25)
})[, 2]
upper <- aggregate(df_valid$mismatch, list(df_valid$bin), function(x) {
  quantile(x, 0.75)
})[, 2]

df_boxplot <- data.frame(bin = unique(df_valid$bin), middle, ymin, ymax, lower, upper)
p_mismatch_human <-
  ggplot(df_boxplot, aes(bin, middle)) +
  geom_point() +
  geom_segment(aes(x = seq(1, 5, 1), y = lower, xend = seq(1, 5, 1), yend = upper), colour = rev(colorspace::heat_hcl(5)), size = 5) +
  geom_segment(aes(x = seq(1, 5, 1), y = ymin, xend = seq(1, 5, 1), yend = ymax), colour = "black") +
  geom_segment(aes(x = seq(0.75, 4.75, 1), y = middle, xend = seq(1.25, 5.25, 1), yend = middle), colour = "black", size = 1) +
  ylim(-30, 25) +
  scale_x_discrete(label = c("[0,0.023]", "(0.023,0.439]", "(0.439,2.94]", "(2.94,25.5]", "(25.5,8010]")) +
  xlab("") +
  ylab(expression("lagging     " * italic(δ)["GSL, MAT"] * " (km yr"^-1 * ")" * "     leading")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 8))

p_x_axis <-
  ggplot(df_boxplot, aes(bin, middle), fill = "white") +
  xlab(expression("Population density (persons km"^-2 * ")")) +
  ylab("") +
  ylim(0, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

cairo_pdf("./figures/Figure 4.pdf", width = 10, height = 6)
grid.arrange(
  grobs = list(
    annotate_figure(p_human + latticeExtra::layer(sp.polygons(land)), fig.lab = "a", fig.lab.pos = "top.left", fig.lab.face = "bold"),
    annotate_figure(p_Vc_human, fig.lab = "b", fig.lab.pos = "top.left", fig.lab.face = "bold"),
    # p_Vc_human,
    p_Vp_human,
    p_angdiff_human,
    p_mismatch_human,
    p_x_axis
  ),
  widths = c(1, 1, 1, 1),
  heights = c(1.6, 2, 0.2),
  layout_matrix = rbind(
    c(1, 1, 1, 1),
    c(2, 3, 4, 5),
    c(6, 6, 6, 6)
  )
)
dev.off()

### Reproject to aeqd
compile_brick <- rasterFromXYZ(df, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
compile_brick <- aggregate(compile_brick, fact = 10, fun = stats::median) # stats:: is needed
ae_proj <- "+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
compile_points <- rasterToPoints(compile_brick, spatial = T)
compile_points_ae <- spTransform(compile_points, ae_proj)
df_ae <- as.data.frame(compile_points_ae)
df_sub <- na.omit(df_ae)
dim(df_sub)

### Fit nonspatial and spatial model
fit.spLM <- function(data, response) {
  # Prepare data frame
  if (response == "Vcnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vcnew)
  if (response == "Vpnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vpnew)
  if (response == "anglediffnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$anglediffnew)
  if (response == "mismatch") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$mismatch)

  coords <- as.matrix(data[, c("x", "y")])

  # Preparing to run Gaussian spatial model; get lognormal model phi (spatial decay parameter) estimate
  z <- resid(lm(response ~ population.density + latitude, spdata)) # residuals of non-spatial model
  spdata$z <- z
  coordinates(spdata) <- spdata[, c("x", "y")]
  test.vgm <- variogram(z ~ 1, data = spdata)
  test.fit <- fit.variogram(test.vgm, model = vgm("Exp"))
  phi <- test.fit[2, ]$range # extract phi value from fitted variogram

  # Set priors: loose priors on beta and residual error variance (tausq), and on spatial variance parameter (sigma.sq), but very tight on Phi (spatial decay parameter).
  priors <- list("beta.Norm" = list(rep(0, 3), diag(100, 3)), "phi.Unif" = c(-log(0.05) / (phi * 100), -log(0.05) / (phi / 100)), "sigma.sq.IG" = c(2, 2), "tau.sq.IG" = c(2, 0.1)) # shape and scale for IG
  # Set starting and tuning values
  starting <- list("phi" = -log(0.05) / phi, "sigma.sq" = 50, "tau.sq" = 1)
  tuning <- list("phi" = (log(0.05) / (phi * 100) - log(0.05) / (phi / 100)) / 10, "sigma.sq" = 0.01, "tau.sq" = 0.01)

  # Knots for Gaussian models
  # knots = kmeans(coords, 20,iter.max=100)$centers

  # Run Gaussian spatial model
  splm <- spLM(response ~ population.density + latitude, data = spdata, coords = coords, cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 10000, n.report = 1000)
  # knots=knots,

  return(splm)
}

### Parameter inference
splm.recover <- function(model) {
  # Recover parameter samples, ignoring the first half of the run as burn-in.
  splm <- spRecover(model, get.beta = TRUE, get.w = TRUE, start = 5001, n.report = 1000)
  beta.hat <- splm$p.beta.recover.samples
  theta.hat <- splm$p.theta.recover.samples

  beta.theta.hat <- cbind(beta.hat, theta.hat)
  # Coefficient summaries
  median <- apply(beta.theta.hat, 2, median)
  lower <- apply(beta.theta.hat, 2, quantile, probs = 0.025)
  upper <- apply(beta.theta.hat, 2, quantile, probs = 0.975)

  coef.summary <- data.frame(median, lower, upper)

  return(list(splm, coef.summary))
}

nonspatial <- function(data, response) {
  # Prepare data frame
  if (response == "Vcnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vcnew)
  if (response == "Vpnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vpnew)
  if (response == "anglediffnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$anglediffnew)
  if (response == "mismatch") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$mismatch)

  nonsplm <- bayesLMConjugate(response ~ population.density + latitude, data = spdata, n.samples = 10000, beta.prior.mean = rep(0, 3), beta.prior.precision = diag(1 / 100, 3), prior.shape = 2, prior.rate = 1 / 0.1)
  nonsplm$p.beta.tauSq.samples <- window(nonsplm$p.beta.tauSq.samples, start = 5001) # burn
  nonsplm$p.beta.tauSq.samples <- mcmc(data = nonsplm$p.beta.tauSq.samples, start = 1, end = 5000, thin = 1) # rename mcmc iterations to 1-5000
  nonsplm$n.samples <- 5000
  class(nonsplm) <- "bayesLMRef"

  beta.hat <- nonsplm$p.beta.tauSq.samples
  median <- apply(beta.hat, 2, median)
  lower <- apply(beta.hat, 2, quantile, probs = 0.025)
  upper <- apply(beta.hat, 2, quantile, probs = 0.975)

  coef.summary <- data.frame(median, lower, upper)

  return(list(nonsplm, coef.summary))
}

### Model validation
predicted.observed <- function(data, response, nonsplm, splm) {
  # Prepare data frame
  if (response == "Vcnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vcnew)
  if (response == "Vpnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$Vpnew)
  if (response == "anglediffnew") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$anglediffnew)
  if (response == "mismatch") spdata <- data.frame(x = data$x, y = data$y, population.density = data$population.density, latitude = data$latitude, response = data$mismatch)

  # nonspatial
  beta0_ppd <- as.numeric(nonsplm$p.beta.tauSq.samples[, 1])
  beta1_ppd <- as.numeric(nonsplm$p.beta.tauSq.samples[, 2])
  beta2_ppd <- as.numeric(nonsplm$p.beta.tauSq.samples[, 3])
  tau_sq_ppd <- as.numeric(nonsplm$p.beta.tauSq.samples[, 4])

  nonsp_median <- nonsp_lower <- nonsp_upper <- rep(NA, dim(data)[1])
  for (i in 1:dim(data)[1]) {
    y_ppd <- rep(NA, 1000)
    for (j in 1:1000) {
      y_ppd[j] <- rnorm(1, (sample(beta0_ppd, 1) + sample(beta1_ppd, 1) * spdata$population.density[i] + sample(beta2_ppd, 1) * spdata$latitude[i]), sqrt(sample(tau_sq_ppd, 1))) # generate ppd for y
    }
    nonsp_median[i] <- quantile(y_ppd, 0.5)
    nonsp_lower[i] <- quantile(y_ppd, 0.025)
    nonsp_upper[i] <- quantile(y_ppd, 0.975)
  }

  nonsp_df <- data.frame(observed = spdata$response, predicted = nonsp_median, lower = nonsp_lower, upper = nonsp_upper)

  # spatial
  beta0_ppd <- as.numeric(splm$p.beta.recover.samples[, 1])
  beta1_ppd <- as.numeric(splm$p.beta.recover.samples[, 2])
  beta2_ppd <- as.numeric(splm$p.beta.recover.samples[, 3])
  # sigma_sq_ppd<-as.numeric(splm$p.theta.recover.samples[,1])
  tau_sq_ppd <- as.numeric(splm$p.theta.recover.samples[, 2])
  # phi_sq_ppd<-as.numeric(splm$p.theta.recover.samples[,3])

  sp_median <- sp_lower <- sp_upper <- rep(NA, dim(data)[1])
  for (i in 1:dim(data)[1]) {
    w_ppd <- as.numeric(splm$p.w.recover.samples[i, ])
    y_ppd <- rep(NA, 1000)
    for (j in 1:1000) {
      y_ppd[j] <- rnorm(1, (sample(beta0_ppd, 1) + sample(beta1_ppd, 1) * spdata$population.density[i] + sample(beta2_ppd, 1) * spdata$latitude[i] + sample(w_ppd, 1)), sqrt(sample(tau_sq_ppd, 1))) # generate ppd for y
    }
    sp_median[i] <- quantile(y_ppd, 0.5)
    sp_lower[i] <- quantile(y_ppd, 0.025)
    sp_upper[i] <- quantile(y_ppd, 0.975)
  }

  sp_df <- data.frame(observed = spdata$response, predicted = sp_median, lower = sp_lower, upper = sp_upper)

  minval <- min(nonsp_df$lower, sp_df$lower, spdata$response, na.rm = T)
  maxval <- max(nonsp_df$upper, sp_df$upper, spdata$response, na.rm = T)

  p1 <-
    ggplot(nonsp_df, aes(x = observed, y = predicted)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, lwd = 0.2, alpha = 0.2) +
    ylim(minval, maxval) +
    xlim(minval, maxval) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Nonspatial model") +
    coord_fixed() +
    theme_classic()

  p2 <-
    ggplot(sp_df, aes(x = observed, y = predicted)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, lwd = 0.2, alpha = 0.2) +
    ylim(minval, maxval) +
    xlim(minval, maxval) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Spatial model") +
    coord_fixed() +
    theme_classic()

  # residuals
  combined <- data.frame(x = spdata$x, y = spdata$y, observed = spdata$response, nonsp.predicted = nonsp_df$predicted, sp.predicted = sp_df$predicted)
  combined$nonsp.resid <- combined$nonsp.predicted - combined$observed
  combined$sp.resid <- combined$sp.predicted - combined$observed

  coordinates(combined) <- c("x", "y")
  proj4string(combined) <- ae_proj
  combined_lonlat <- spTransform(combined, proj4string(C))
  gridded(combined_lonlat) <- TRUE
  combined_brick <- brick(combined_lonlat)

  col_pal_diverge <- c(colorspace::sequential_hcl(10, h = 255, l = c(40, 90), power = 1), rev(colorspace::sequential_hcl(10, h = 330, l = c(40, 90), power = 1)))
  p3 <- levelplot(combined_brick[[4]], par.settings = my.settings, margin = FALSE, col.regions = col_pal_diverge, main = "Nonspatial model", xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8))

  p4 <- levelplot(combined_brick[[5]], par.settings = my.settings, margin = FALSE, col.regions = col_pal_diverge, main = "Spatial model", xlab = list(cex = 0), ylab = list(cex = 0), scales = list(y = list(at = seq(30, 70, 20)), x = list(at = seq(-180, 180, 60))), colorkey = list(space = "bottom", width = 0.8, height = 0.8))

  # Variogram
  nonsp.var <- variogram(nonsp.resid ~ 1, data = combined)
  nonsp.fit <- fit.variogram(nonsp.var, model = vgm("Exp"))
  nonsp.line <- variogramLine(vgm(psill = nonsp.fit[2, 2], model = "Exp", range = nonsp.fit[2, 3], nugget = 0), nonsp.var[15, 2])

  sp.var <- variogram(sp.resid ~ 1, data = combined)
  sp.fit <- fit.variogram(sp.var, model = vgm("Exp"))
  sp.line <- variogramLine(vgm(psill = sp.fit[2, 2], model = "Exp", range = sp.fit[2, 3], nugget = 0), sp.var[15, 2])

  nonsp.var$model <- "Nonspatial model"
  nonsp.line$model <- "Nonspatial model"
  sp.var$model <- "Spatial model"
  sp.line$model <- "Spatial model"
  var.df <- rbind(nonsp.var, sp.var)
  line.df <- rbind(nonsp.line, sp.line)

  # maxval<-max(nonsp.var$gamma,sp.var$gamma)

  p5 <-
    ggplot(var.df, aes(dist, gamma)) +
    geom_point(aes(col = model)) +
    geom_line(data = line.df, aes(dist, gamma, col = model)) +
    xlab("Distance") +
    ylab("Semivariogram") +
    # ylim(0,maxval)+
    theme_classic() +
    theme(legend.position = "top")

  return(list(p1, p2, p3 + latticeExtra::layer(sp.polygons(land)), p4 + latticeExtra::layer(sp.polygons(land)), p5))
}

### Summarize models for four climate-phenology metrics and make diagnostic plots
make_diag_plot <- function(response) {
  if (response == "Vcnew") label <- "a"
  if (response == "Vpnew") label <- "b"
  if (response == "anglediffnew") label <- "c"
  if (response == "mismatch") label <- "d"

  nonsplm <- nonspatial(df_sub, response) # nonspatial model
  splm <- fit.spLM(df_sub, response) # spatial model
  splm <- splm.recover(splm) # recover spatial model

  nonsplm_mcmc <- ggs(nonsplm[[1]]$p.beta.tauSq.samples)
  splm_mcmc1 <- ggs(splm[[1]]$p.beta.recover.samples)
  splm_mcmc2 <- ggs(splm[[1]]$p.theta.recover.samples)

  p1 <- ggs_traceplot(nonsplm_mcmc) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # nonspatial model traceplots
  p2 <- ggs_traceplot(splm_mcmc1) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # spatial model traceplots 1
  p3 <- ggs_traceplot(splm_mcmc2) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # spatial model traceplots 2

  p4 <- ggs_density(nonsplm_mcmc) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # nonspatial model density plots
  p5 <- ggs_density(splm_mcmc1) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # spatial model density plots 1
  p6 <- ggs_density(splm_mcmc2) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0)) # spatial model density plots 2

  p7 <- grid.arrange(p1, p4, ncol = 2, nrow = 1)
  p8 <- grid.arrange(p2, p5, ncol = 2, nrow = 1)
  p9 <- grid.arrange(p3, p6, ncol = 2, nrow = 1)

  # MCMC
  file <- paste("./figures/regression/", response, "_1.pdf", sep = "")
  pdf(file, width = 10, height = 10)

  grid.arrange(
    grobs = list(
      annotate_figure(p7, top = text_grob("Nonspatial model"), fig.lab = label, fig.lab.pos = "top.left", fig.lab.face = "bold"),
      annotate_figure(p8, top = text_grob("Spatial model")),
      p9
    ),
    widths = c(1, 1),
    heights = c(1, 1),
    layout_matrix = rbind(
      c(1, 2),
      c(1, 3)
    )
  )
  dev.off()

  p1.1 <- ggAcf(nonsplm[[1]]$p.beta.tauSq.samples[, 1]) + theme_classic() + ggtitle("(Intercept)") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial acf
  p1.2 <- ggAcf(nonsplm[[1]]$p.beta.tauSq.samples[, 2]) + theme_classic() + ggtitle("population.density") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial acf
  p1.3 <- ggAcf(nonsplm[[1]]$p.beta.tauSq.samples[, 3]) + theme_classic() + ggtitle("latitude") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial acf
  p1.4 <- ggAcf(nonsplm[[1]]$p.beta.tauSq.samples[, 4]) + theme_classic() + ggtitle("tau.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial acf

  p2.1 <- ggPacf(nonsplm[[1]]$p.beta.tauSq.samples[, 1]) + theme_classic() + ggtitle("(Intercept)") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial pacf
  p2.2 <- ggPacf(nonsplm[[1]]$p.beta.tauSq.samples[, 2]) + theme_classic() + ggtitle("population.density") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial pacf
  p2.3 <- ggPacf(nonsplm[[1]]$p.beta.tauSq.samples[, 3]) + theme_classic() + ggtitle("latitude") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial pacf
  p2.4 <- ggPacf(nonsplm[[1]]$p.beta.tauSq.samples[, 4]) + theme_classic() + ggtitle("tau.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # nonspatial pacf

  p3.1 <- ggAcf(splm[[1]]$p.beta.recover.samples[, 1]) + theme_classic() + ggtitle("(Intercept)") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf
  p3.2 <- ggAcf(splm[[1]]$p.beta.recover.samples[, 2]) + theme_classic() + ggtitle("population.density") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf
  p3.3 <- ggAcf(splm[[1]]$p.beta.recover.samples[, 3]) + theme_classic() + ggtitle("latitude") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf
  p3.4 <- ggAcf(splm[[1]]$p.theta.recover.samples[, 1]) + theme_classic() + ggtitle("sigma.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf
  p3.5 <- ggAcf(splm[[1]]$p.theta.recover.samples[, 2]) + theme_classic() + ggtitle("tau.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf
  p3.6 <- ggAcf(splm[[1]]$p.theta.recover.samples[, 3]) + theme_classic() + ggtitle("phi") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial acf

  p4.1 <- ggPacf(splm[[1]]$p.beta.recover.samples[, 1]) + theme_classic() + ggtitle("(Intercept)") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf
  p4.2 <- ggPacf(splm[[1]]$p.beta.recover.samples[, 2]) + theme_classic() + ggtitle("population.density") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf
  p4.3 <- ggPacf(splm[[1]]$p.beta.recover.samples[, 3]) + theme_classic() + ggtitle("latitude") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf
  p4.4 <- ggPacf(splm[[1]]$p.theta.recover.samples[, 1]) + theme_classic() + ggtitle("sigma.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf
  p4.5 <- ggPacf(splm[[1]]$p.theta.recover.samples[, 2]) + theme_classic() + ggtitle("tau.sq") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf
  p4.6 <- ggPacf(splm[[1]]$p.theta.recover.samples[, 3]) + theme_classic() + ggtitle("phi") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) # spatial pacf

  p7 <- grid.arrange(p1.1, p2.1, p1.2, p2.2, p1.3, p2.2, p1.4, p2.4, nrow = 4, ncol = 2)
  p8 <- grid.arrange(p3.1, p4.1, p3.2, p4.2, p3.3, p4.2, p3.4, p4.4, p3.5, p4.5, p3.6, p4.6, nrow = 6, ncol = 2)

  # ACF and PACF
  file <- paste("./figures/regression/", response, "_2.pdf", sep = "")
  pdf(file, width = 10, height = 10)

  grid.arrange(
    grobs = list(
      annotate_figure(p7, top = text_grob("Nonspatial model"), fig.lab = label, fig.lab.pos = "top.left", fig.lab.face = "bold"),
      annotate_figure(p8, top = text_grob("Spatial model"))
    ),
    ncol = 2, nrow = 1
  )
  dev.off()

  # summary statistics (Table S3, S5)
  file <- paste("./figures/regression/", response, "_3.pdf", sep = "")
  pdf(file, width = 10, height = 6)

  nonspcoef <- nonsplm[[2]] # nonspatial model coef summary
  spcoef <- splm[[2]] # spatial model coef summary
  grid.arrange(text_grob("Nonspatial", just = "left"), tableGrob(nonspcoef), text_grob("Spatial", just = "left"), tableGrob(spcoef), ncol = 1, nrow = 4, heights = c(1, 4, 1, 6))

  dev.off()

  # diagnostics (Table S4)
  file <- paste("./figures/regression/", response, "_4.pdf", sep = "")
  pdf(file, width = 10, height = 3)

  nonspdiag <- spDiag(nonsplm[[1]]) # nonspatial model diagnostics
  nonspdiag <- rbind(nonspdiag$DIC, nonspdiag$GP, nonspdiag$GRS)
  rownames(nonspdiag)[8] <- "GRS"
  spdiag <- spDiag(splm[[1]]) # spatial model diagnostics
  spdiag <- rbind(spdiag$DIC, spdiag$GP, spdiag$GRS)
  rownames(spdiag)[8] <- "GRS"
  diag <- cbind(nonspdiag, spdiag)
  colnames(diag) <- c("Nonspatial", "Spatial")
  grid.table(diag)

  dev.off()

  # predictive error
  compare <- predicted.observed(df_sub, response, nonsplm[[1]], splm[[1]])
  file <- paste("./figures/regression/", response, "_5.pdf", sep = "")
  pdf(file, width = 10, height = 5)

  grid.arrange(annotate_figure(compare[[1]], fig.lab = label, fig.lab.pos = "top.left", fig.lab.face = "bold"), compare[[2]], ncol = 2, nrow = 1) # predicted vs. observed

  dev.off()

  # map of residuals (Figure S6)
  file <- paste("./figures/regression/", response, "_6.pdf", sep = "")
  pdf(file, width = 10, height = 6)

  grid.arrange(annotate_figure(compare[[3]], fig.lab = label, fig.lab.pos = "top.left", fig.lab.face = "bold"), compare[[4]], ncol = 1, nrow = 2) # residual maps

  dev.off()

  # Semivariograms (Figure S7)
  file <- paste("./figures/regression/", response, "_7.pdf", sep = "")
  pdf(file, width = 5, height = 3)

  grid.arrange(annotate_figure(compare[[5]], fig.lab = label, fig.lab.pos = "top.left", fig.lab.face = "bold"), ncol = 1, nrow = 1) # residual variograms

  dev.off()
}

dir.create("./figures/regression/")
make_diag_plot("Vcnew")
make_diag_plot("Vpnew")
make_diag_plot("anglediffnew")
make_diag_plot("mismatch")

##### MAT and elevation correlation (Figure S5)
elevation <- elevatr::get_elev_raster(data.frame(long = c(-180, 0, -180), lat = c(30, 50, 70)), prj = "+proj=longlat +datum=WGS84 +no_defs", z = 1)
elevation <- crop(elevation, e)
elev_re <- resample(elevation, C[[1]])
df_elev <- as.data.frame(elev_re, xy = T)
colnames(df_elev) <- c("x", "y", "elev")
df_Vc <- as.data.frame(Vc, xy = T) %>%
  rename(Vc = "layer")
df_mnC <- as.data.frame(calc(C, mean), xy = T) %>%
  rename(mnC = "layer")
df_combined <- left_join(df_elev, df_mnC, by = c("x", "y")) %>%
  drop_na() %>%
  filter(y < 70 & y > 30)

summary(lm(data = df_combined, mnC ~ elev))

cairo_pdf("./figures/Figure S5.pdf")
ggplot(df_combined) +
  stat_binhex(aes(x = elev, y = mnC), bins = 100) +
  xlab("Elevation (m)") +
  ylab("Mean annual temperature (°C)") +
  xlim(0, max(df_combined$elev)) +
  theme_classic() +
  scale_fill_viridis_c()
dev.off()
