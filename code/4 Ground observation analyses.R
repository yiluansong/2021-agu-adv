# Code for 2.2 (leading to Figure 5)
source("./tools/packages.R")

##### Ground observation analyses (Figure 5)
### Get human population
file <- c("./data/gpw_v4_population_density_rev11_30_min.nc")
Human <- brick(file)[[1]]

### Get NPN data
path_npn <- "./data/NPN/"
dir.create(path_npn, recursive = T)
rnpn::npn_phenophases() %>%
  filter(pheno_class_id == 1) %>%
  arrange(phenophase_name)

species_list <- rnpn::npn_species() %>%
  filter(kingdom == "Plantae") %>%
  rowwise() %>%
  filter(sum(species_type$Species_Type == "Calibration") > 0)
# https://www.usanpn.org/files/articles/developing_a_plant_profile.pdf

# run only once to download data
# for (i in 1:nrow(species_list)) {
#   npn_df <- rnpn::npn_download_site_phenometrics(
#     request_source = "YS",
#     years = as.character(1980:2021),
#     species_ids = species_list$species_id[i],
#     pheno_class_ids = 1,
#     climate_data = T
#   )
#   write_csv(npn_df, paste0(path_npn, species_list$species_id[i], ".csv"))
#   print(i)
# }

### Calculate mismatch
data_list <- vector(mode = "list")
for (i in 1:nrow(species_list)) {
  npn_df <- read_csv(paste0(path_npn, species_ids = species_list$species_id[i], ".csv"))

  if (nrow(npn_df) > 0) {
    # doy~temp regression
    data <- npn_df %>%
      dplyr::select(species_id,
        genus,
        species,
        common_name,
        lat = latitude,
        lon = longitude,
        doy = mean_first_yes_doy,
        tmax_spring,
        tmin_spring,
        gdd = mean_gdd
      ) %>%
      filter(
        doy != -9999,
        gdd != -9999,
        tmax_spring != -9999,
        tmin_spring != -9999
      )

    summary(model <- lm(data$doy ~ data$gdd + data$tmax_spring + data$tmin_spring))

    data$predict <- model$fitted.values
    data$mis <- model$residuals
    data$abs_mis <- abs(model$residuals)

    # retrieve human population density
    npn_coords_sp <- SpatialPoints(
      coords = data[, c("lon", "lat")],
      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs ")
    )

    data$pop <- raster::extract(Human, npn_coords_sp) + 1

    data_list[[i]] <- data
  }
}

data_all <- bind_rows(data_list) %>%
  group_by(species_id) %>%
  filter(n() > 100) %>%
  ungroup() %>%
  mutate(log_pop = log(pop + 1))

### Summarize
quantile(data_all$abs_mis, c(0.025, 0.5, 0.975))

species_list_common_name <- unique(data_all$common_name)
coef_all <- newdata_all <- vector(mode = "list")
for (i in 1:length(species_list_common_name)) {
  data_sp <- data_all %>% filter(common_name == species_list_common_name[i])
  model_sp <- lm(abs(mis) ~ log_pop + lat, data = data_sp)
  model_sp_summ <- summary(model_sp)
  coef <- data.frame(
    common_name = species_list_common_name[i],
    estimate = model_sp_summ$coefficients[2, 1],
    p.value = model_sp_summ$coefficients[2, 4],
    n = nrow(data_sp)
  ) %>%
    mutate(p.value = round(p.value * 10000) / 10000)
  coef_all[[i]] <- coef

  newdata_sp <- data.frame(log_pop = seq(min(data_sp$log_pop), max(data_sp$log_pop), length.out = 100), lat = 50, common_name = species_list_common_name[i])
  newdata_sp <- cbind(newdata_sp, predict.lm(model_sp, newdata = newdata_sp, interval = "confidence"))
  newdata_all[[i]] <- newdata_sp
}

coef_all <- bind_rows(coef_all)
newdata_all <- bind_rows(newdata_all)
coef_all
# write_csv(coef_all,  paste0(path_npn,"coef summary.csv"))

model <- lm(abs_mis ~ log_pop + lat, data = data_all)
summary(model)

### Plot
png(
  filename = "./figures/Figure 5_map.png",
  width = 4800, height = 4800, bg = "transparent", res = 600
)
us_map <- map_data("world") %>% filter(region == "USA")
p_map <- ggplot(data_all) +
  geom_polygon(data = us_map %>% filter(lat > -150), aes(x = long, y = lat, group = group), color = "gray", fill = "gray") +
  coord_map(projection = "albers", parameters = c(20, 50)) +
  geom_point(aes(x = lon, y = lat, col = common_name), alpha = 0.25) +
  xlim(-170, -70) +
  theme_void() +
  theme(legend.position = "none")
print(p_map)
dev.off()

cairo_pdf("./figures/Figure 5_trend.pdf", width = 6, height = 6)
newdata <- data.frame(log_pop = seq(min(data_all$log_pop), max(data_all$log_pop), length.out = 100), lat = 50)
newdata <- cbind(newdata, predict.lm(model, newdata = newdata, interval = "confidence"))
p_trend <- ggplot() +
  geom_point(data = data_all, aes(x = log_pop, y = (abs_mis)^(1 / 2), col = common_name), alpha = 0.1) +
  geom_line(data = newdata, aes(x = log_pop, y = (fit)^(1 / 2)), lwd = 1.5, alpha = 0.5) +
  geom_ribbon(data = newdata, aes(x = log_pop, ymin = (lwr)^(1 / 2), ymax = (upr)^(1 / 2)), alpha = 0.25) +
  geom_line(data = newdata_all, aes(x = log_pop, y = (fit)^(1 / 2), col = common_name, group = common_name), alpha = 0.75) +
  scale_y_continuous(
    breaks = sqrt(c(0, 4, 16, 36, 64)),
    labels = c(0, 4, 16, 36, 64)
  ) +
  scale_x_continuous(
    breaks = log(c(0, 10, 100, 1000) + 1),
    labels = c(0, 10, 100, 1000)
  ) +
  ylab("Absolute predictive error (day)") +
  xlab(expression("Population density (persons km"^-2 * ")")) +
  theme_classic() +
  theme(legend.position = "right") +
  guides(col = guide_legend(title = "Species"))
print(p_trend)
dev.off()
