###############################################################################
# LOAD LIBRARIES
###############################################################################
###############################################################################
# BAYESIAN SPATIOTEMPORAL CLIMATE MODELLING – SOUTH AFRICA
###############################################################################

library(INLA)
library(sf)
library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(terra)
library(scales)
library(lattice)
library(latticeExtra)
library(latex2exp)

set.seed(123)

###############################################################################
# 1. DATA PREPARATION
###############################################################################

clim <- read.csv("C:/Users/HP/Desktop/Research 2025/rainfall.temp SA data.csv")

clim <- clim %>%
  filter(
    !is.na(rainfall), rainfall > 0,
    !is.na(temp.max), !is.na(temp.min),
    !is.na(Longitude), !is.na(Latitude)
  ) %>%
  mutate(
    log_rain = log(rainfall),
    time = (year - min(year)) * 12 + Months,
    tmax = as.numeric(scale(temp.max)),
    tmin = as.numeric(scale(temp.min)),
    t_interact = tmax * tmin
  )

coords <- as.matrix(clim[, c("Longitude", "Latitude")])

###############################################################################
# 2. SPDE MESH (STABLE SETTINGS)
###############################################################################

mesh <- inla.mesh.2d(
  loc = coords,
  max.edge = c(2, 8),
  cutoff = 0.5
)

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(300, 0.5),
  prior.sigma = c(1, 0.5)
)

A <- inla.spde.make.A(mesh, loc = coords)
sidx <- inla.spde.make.index("spatial", spde$n.spde)

###############################################################################
# 3. INLA STACK
###############################################################################

stack <- inla.stack(
  data = list(y = clim$log_rain),
  A = list(A, 1),
  effects = list(
    spatial = sidx,
    data.frame(
      intercept = 1,
      time = clim$time,
      tmax = clim$tmax,
      tmin = clim$tmin,
      t_interact = clim$t_interact
    )
  ),
  tag = "est"
)

###############################################################################
# 4. COMPETING BAYESIAN MODELS
###############################################################################

# M0: Baseline (no space, no time)
f_m0 <- y ~ 1 + tmax + tmin

# M1: Temporal trend only
f_m1 <- y ~ 1 + tmax + tmin +
  f(time, model = "rw1", scale.model = TRUE)

# M2: Spatial only
f_m2 <- y ~ 1 + tmax + tmin +
  f(spatial, model = spde)

# M3 (PROPOSED): Spatiotemporal + interaction
f_proposed <- y ~ 1 + tmax + tmin + t_interact +
  f(time, model = "rw1", scale.model = TRUE) +
  f(spatial, model = spde)

###############################################################################
###############################################################################
# INLA CONTROL SETTINGS (STABLE & COMPARABLE)
###############################################################################

ctrl_pred <- list(
  A = inla.stack.A(stack),
  compute = TRUE,
  link = 1
)

ctrl_comp <- list(
  dic = TRUE,
  waic = TRUE,
  cpo = TRUE,
  mlik = TRUE,
  config = TRUE
)

ctrl_inla <- list(
  strategy = "simplified.laplace",
  int.strategy = "eb",
  diagonal = 1e-5,
  tolerance = 1e-6
)

###############################################################################
# FIT ALL MODELS
###############################################################################

models <- list(
  M0_Baseline = f_m0,
  M1_Temporal = f_m1,
  M2_Spatial = f_m2,
  M3_Proposed = f_proposed
)

fits <- list()

for (m in names(models)) {
  cat("Fitting", m, "...\n")
  fits[[m]] <- inla(
    formula = models[[m]],
    data = inla.stack.data(stack),
    family = "gaussian",
    control.predictor = ctrl_pred,
    control.compute = ctrl_comp,
    control.inla = ctrl_inla,
    verbose = FALSE
  )
}

###############################################################################
# MODEL COMPARISON TABLE
###############################################################################

model_comparison <- data.frame(
  Model = names(fits),
  DIC = sapply(fits, function(x) x$dic$dic),
  WAIC = sapply(fits, function(x) x$waic$waic),
  Marginal_LogLik = sapply(fits, function(x) x$mlik[1, 1]),
  Mean_CPO = sapply(fits, function(x) mean(-log(x$cpo$cpo), na.rm = TRUE))
)

model_comparison <- model_comparison %>%
  mutate(
    Delta_DIC = DIC - min(DIC),
    Delta_WAIC = WAIC - min(WAIC),
    WAIC_Rank = rank(WAIC)
  ) %>%
  arrange(WAIC_Rank)

###############################################################################
# OUTPUT RESULTS
###############################################################################

cat("\n", strrep("=", 65), "\n")
cat("BAYESIAN MODEL COMPARISON (INLA)\n")
cat(strrep("=", 65), "\n")
print(model_comparison)

###############################################################################
fit_proposed <- fits$M3_Proposed


###############################################################################
# OUT-OF-SAMPLE CROSS VALIDATION (K-FOLD) – CORRECTED VERSION
###############################################################################

set.seed(123)

K <- 5
n <- nrow(clim)

fold_id <- sample(rep(1:K, length.out = n))

cv_results <- data.frame(
  Fold = integer(),
  RMSE_no_interaction = numeric(),
  RMSE_interaction = numeric(),
  MAE_no_interaction = numeric(),
  MAE_interaction = numeric()
)

###############################################################################
# MODEL FORMULAS
###############################################################################

f_no_interact <- y ~ 1 + tmax + tmin +
  f(time, model="rw1", scale.model=TRUE) +
  f(spatial, model=spde)

f_with_interact <- y ~ 1 + tmax + tmin + t_interact +
  f(time, model="rw1", scale.model=TRUE) +
  f(spatial, model=spde)

###############################################################################
# CROSS VALIDATION LOOP
###############################################################################

for(k in 1:K){
  
  cat("Running fold:", k, "\n")
  
  train_index <- which(fold_id != k)
  test_index  <- which(fold_id == k)
  
  train_data <- clim[train_index,]
  test_data  <- clim[test_index,]
  
  coords_train <- as.matrix(train_data[,c("Longitude","Latitude")])
  
  A_train <- inla.spde.make.A(mesh, loc = coords_train)
  
  stack_train <- inla.stack(
    data=list(y=train_data$log_rain),
    A=list(A_train,1),
    effects=list(
      spatial=sidx,
      data.frame(
        intercept=1,
        time=train_data$time,
        tmax=train_data$tmax,
        tmin=train_data$tmin,
        t_interact=train_data$t_interact
      )
    ),
    tag="train"
  )
  
  ctrl_pred_cv <- list(
    A=inla.stack.A(stack_train),
    compute=TRUE
  )
  
  ###########################################################################
  # FIT MODELS
  ###########################################################################
  
  fit_no_int <- inla(
    f_no_interact,
    data=inla.stack.data(stack_train),
    family="gaussian",
    control.predictor=ctrl_pred_cv,
    control.inla=ctrl_inla,
    verbose=FALSE
  )
  
  fit_int <- inla(
    f_with_interact,
    data=inla.stack.data(stack_train),
    family="gaussian",
    control.predictor=ctrl_pred_cv,
    control.inla=ctrl_inla,
    verbose=FALSE
  )
  
  ###########################################################################
  # EXTRACT FIXED EFFECTS
  ###########################################################################
  
  beta_no_int <- fit_no_int$summary.fixed
  beta_int <- fit_int$summary.fixed
  
  b0_no <- beta_no_int["(Intercept)","mean"]
  b1_no <- beta_no_int["tmax","mean"]
  b2_no <- beta_no_int["tmin","mean"]
  
  b0_int <- beta_int["(Intercept)","mean"]
  b1_int <- beta_int["tmax","mean"]
  b2_int <- beta_int["tmin","mean"]
  b3_int <- beta_int["t_interact","mean"]
  
  ###########################################################################
  # PREDICTIONS
  ###########################################################################
  
  pred_no_int <- b0_no +
    b1_no * test_data$tmax +
    b2_no * test_data$tmin
  
  pred_int <- b0_int +
    b1_int * test_data$tmax +
    b2_int * test_data$tmin +
    b3_int * test_data$t_interact
  
  ###########################################################################
  # COMPUTE PREDICTION ERRORS
  ###########################################################################
  
  obs <- test_data$log_rain
  
  rmse_no_int <- sqrt(mean((obs - pred_no_int)^2, na.rm=TRUE))
  rmse_int <- sqrt(mean((obs - pred_int)^2, na.rm=TRUE))
  
  mae_no_int <- mean(abs(obs - pred_no_int), na.rm=TRUE)
  mae_int <- mean(abs(obs - pred_int), na.rm=TRUE)
  
  cv_results <- rbind(
    cv_results,
    data.frame(
      Fold=k,
      RMSE_no_interaction=rmse_no_int,
      RMSE_interaction=rmse_int,
      MAE_no_interaction=mae_no_int,
      MAE_interaction=mae_int
    )
  )
}

###############################################################################
# SUMMARY OF CROSS VALIDATION RESULTS
###############################################################################

cv_summary <- cv_results %>%
  summarise(
    Mean_RMSE_no_interaction = mean(RMSE_no_interaction, na.rm=TRUE),
    Mean_RMSE_interaction = mean(RMSE_interaction, na.rm=TRUE),
    Mean_MAE_no_interaction = mean(MAE_no_interaction, na.rm=TRUE),
    Mean_MAE_interaction = mean(MAE_interaction, na.rm=TRUE)
  )

cat("\n",strrep("=",60),"\n")
cat("OUT-OF-SAMPLE CROSS VALIDATION RESULTS\n")
cat(strrep("=",60),"\n")

print(cv_results)
print(cv_summary)

###############################################################################
# READ SOUTH AFRICA SHAPEFILE
###############################################################################

sa_shape <- st_read(
  "C:/Users/HP/Desktop/SA SHP/zaf_admin2.shp",
  quiet = TRUE
) |> 
  st_transform(4326)

###############################################################################
# PROJECT SPATIAL FIELD TO A GRID THAT MATCHES SA EXTENT
###############################################################################

# Get SA bounding box (ensures full coverage)
bbox <- st_bbox(sa_shape)

proj <- inla.mesh.projector(
  mesh,
  xlim = c(bbox["xmin"], bbox["xmax"]),
  ylim = c(bbox["ymin"], bbox["ymax"]),
  dims = c(300, 300)   # higher resolution for smooth coverage
)

field_base <- inla.mesh.project(
  proj,
  fit_proposed$summary.random$spatial$mean
)

###############################################################################
# CREATE GRID DATAFRAME
###############################################################################

grid_df <- expand.grid(
  lon = proj$x,
  lat = proj$y
)

grid_df$rain_spatial <- as.vector(field_base)
grid_df$interaction_spatial <- as.vector(field_base)
###############################################################################
# CONVERT TO RASTER & MASK TO SA
###############################################################################

r <- rast(grid_df, type = "xyz", crs = "EPSG:4326")

r_masked <- mask(r, vect(sa_shape))

grid_df_masked <- as.data.frame(r_masked, xy = TRUE)
names(grid_df_masked) <- c("lon", "lat", "rain_spatial")

###############################################################################
# COMPUTE FULL COLOUR RANGE
###############################################################################

rain_range <- range(grid_df_masked$rain_spatial, na.rm = TRUE)

###############################################################################
# PLOT — FULL SA COVERAGE
###############################################################################
#dat<- read.csv("trend.csv", h=T)
#head(dat)

#dim(dat)

# Aggregate yearly national averages
trend_data <- clim %>%
  filter(year >= 1980, year <= 2022) %>%
  group_by(year) %>%  # Group by year only
  summarise(
    rain = mean(rainfall, na.rm = TRUE),
    tmax = mean(temp.max, na.rm = TRUE),
    tmin = mean(temp.min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year)


###############################################################################
# Custom X-axis for spacing
###############################################################################
###############################################################################
# Custom X-axis for spacing
###############################################################################

# Get the actual year values from your data
years <- unique(trend_data$year)
year_breaks <- seq(1980, 2022, by = 2)
x_scale <- list(
  rot = 90,                    # rotate labels for readability
  cex = 0.8,                   # reduce label size slightly
  alternating = 1,             # avoid label crowding
  at = seq(1, length(years), by = 2),  # Show every 2nd year position
  labels = years[seq(1, length(years), by = 2)]  # Corresponding year labels
)
#############################################################################

# Get the actual year values from your data
years <- unique(trend_data$year)

dat.graf1 <- with(trend_data, data.frame(year=year, Temp=rep(c("Maximum","Minimum"), each=nrow(trend_data)),
                                  Y=c(tmax, tmin)))
head(dat.graf1)

graf1 <- xyplot(Y ~ factor(year), groups=Temp, data=dat.graf1, type="l",
                ylab=TeX("Temperature ($^{\\circ}C$)"), col=c("red", "blue"), 
                ylim=c(0, 30), lwd=2, scales = list(x = list(at = which(unique(dat.graf1$year) %in% year_breaks), labels = year_breaks,
                                                              rot = 45))) # Rotate labels if needed

graf1
graf2 <- barchart(rain ~ factor(year), data=trend_data, ylab="Rainfall (mm)", ylim=c(0,80), col="#ECFF7A",
                  xlab="Year") 

graf2

pdf("maxmim.pdf", width=16, height=12)
update(doubleYScale(graf2, graf1, add.ylab2=TRUE, use.style=FALSE, 
                    text=c("Maximum","Minimum"), columns=2),
       par.settings = simpleTheme(col.line = c('red','blue')))
dev.off()

##############################################################################
p_rain <- ggplot() +
  geom_raster(
    data = grid_df_masked,
    aes(lon, lat, fill = rain_spatial)
  ) +
  geom_sf(
    data = sa_shape,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,          # strongest = darkest
    limits = rain_range,
    oob = squish,
    name = "Rainfall Effect\n(Stronger → Weaker)",
    na.value = NA
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(x = "Longitude", y = "Latitude")

p_rain

########################################################################
###############################################################################
# RAINFALL EFFECT AND TEMP–RAIN INTERACTION
###############################################################################

#library(terra)
# First, create the interaction_spatial column in grid_df
grid_df$interaction_spatial <- as.vector(field_base)

# Convert projected grid to raster
r <- rast(grid_df, type = "xyz", crs = "EPSG:4326")

# Mask to SA boundary
sa_vect <- vect(sa_shape)
r_masked <- mask(r, sa_vect)

# Convert back to dataframe
grid_df_masked <- as.data.frame(r_masked, xy = TRUE)
names(grid_df_masked) <- c("lon", "lat", "rain_spatial", "interaction_spatial")

# Remove NA values
grid_df_masked <- na.omit(grid_df_masked)

# Compute ranges
rain_range <- range(grid_df_masked$rain_spatial, na.rm = TRUE)
inter_range <- range(grid_df_masked$interaction_spatial, na.rm = TRUE)

###############################################################################
# RAINFALL EFFECT MAP
###############################################################################

p_rain <- ggplot() +
  geom_raster(
    data = grid_df_masked,
    aes(x = lon, y = lat, fill = rain_spatial)
  ) +
  geom_sf(
    data = sa_shape,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    limits = rain_range,
    oob = scales::squish,
    name = "Rainfall Effect\n(Stronger → Weaker)",
    na.value = "transparent"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme

###############################################################################
# INTERACTION EFFECT MAP
###############################################################################

p_inter <- ggplot() +
  geom_raster(
    data = grid_df_masked,
    aes(x = lon, y = lat, fill = interaction_spatial)
  ) +
  geom_sf(
    data = sa_shape,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    option = "D",
    direction = -1,
    limits = inter_range,
    oob = scales::squish,
    name = "Temp–Rain Interaction\n(Stronger → Weaker)",
    na.value = "transparent"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme

###############################################################################
# DISPLAY
###############################################################################

p_rain
p_inter
###############################################################################

###############################################################################
# POINT TEMPERATURE MAPS
###############################################################################

p_tmax <- ggplot() +
  geom_sf(data = sa_shape, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = clim, aes(Longitude, Latitude, color = temp.max), size = 1) +
  scale_color_viridis_c(
    option = "A",
    direction = -1,   # Darker at top
    name = "Max Temp (°C)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme


p_tmin <- ggplot() +
  geom_sf(data = sa_shape, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = clim, aes(Longitude, Latitude, color = temp.min), size = 1) +
  scale_color_viridis_c(
    option = "B",
    direction = -1,   # Darker at top
    name = "Min Temp (°C)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme


print(p_tmax)
print(p_tmin)


###############################################################################
#SEASONAL FACETED MAPS 
###############################################################################
clim <- clim %>%
  mutate(
    season = case_when(
      Months %in% c(12, 1, 2) ~ "WIN",
      Months %in% c(3, 4, 5)  ~ "SPG",
      Months %in% c(6, 7, 8)  ~ "SUM",
      TRUE                   ~ "AUT"
    )
  )

seasonal <- clim %>%
  group_by(season, Longitude, Latitude) %>%
  summarise(
    rain = mean(rainfall, na.rm = TRUE),
    tmax = mean(temp.max, na.rm = TRUE),
    tmin = mean(temp.min, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(seasonal) +
  geom_sf(data = sa_shape, inherit.aes = FALSE,
          fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_point(aes(Longitude, Latitude, color = rain), size = 1.5) +
  scale_color_viridis_c(
    option = "C",
    direction = -1,   # Darker at top, lighter at bottom
    name = "Seasonal Rainfall"
  ) +
  facet_wrap(~season) +
  coord_sf(expand = FALSE) +
  scale_x_continuous(
    breaks = seq(
      floor(min(seasonal$Longitude)),
      ceiling(max(seasonal$Longitude)),
      by = 4
    ),
    labels = function(x) paste0(x, "°E")
  ) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme


###############################################################################
# CLIMATE CHANGE HOTSPOT DETECTION (SPATIAL SLOPES)
###############################################################################

hotspots <- clim %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    trend = coef(lm(temp.max ~ year))[2],
    .groups = "drop"
  )

ggplot(hotspots) +
  geom_sf(data = sa_shape,
          fill = "grey95",
          color = "grey40",
          linewidth = 0.3) +
  geom_point(aes(Longitude, Latitude, color = trend), size = 1) +
  scale_color_viridis_c(
    option = "E",
    direction = -1,   # ← Makes top darker, bottom lighter
    name = "Temp Trend (°C/year)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  gis_theme


###############################################################################
# 5. PIT HISTOGRAM – M3 DIAGNOSTIC
###############################################################################

m3 <- fits$M3_Proposed

# Extract PIT values
pit_vals <- m3$cpo$pit
pit_vals <- pit_vals[!is.na(pit_vals)]

# Plot PIT histogram
ggplot(data.frame(pit = pit_vals), aes(x = pit)) +
  geom_histogram(bins = 20, color = "black", fill = "skyblue") +
  geom_hline(yintercept = length(pit_vals)/20, linetype = "dashed", color = "red") +
  labs(
    #title = "PIT Histogram – M3 Spatiotemporal Model",
    x = "PIT values",
    y = "Frequency"
  ) +
  theme_minimal()



ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", color = "grey60") +
  geom_sf(data = clim_sf, aes(color = residual_std), size = 1.8) +
  scale_color_viridis(
    option = "C",
    direction = -1,
    name = "Std Residual"
  ) +
  labs(
    #title = "Standardized Spatial Residuals – M3 Spatiotemporal Model"
  ) +
  theme_minimal()


###############################################################################
# 1. EXTRACT FITTED VALUES CORRESPONDING TO OBSERVED DATA
###############################################################################

# Get indices for estimation stack
index_est <- inla.stack.index(stack, tag = "est")$data

# Extract fitted values for observations only
fitted_vals <- fit_proposed$summary.fitted.values$mean[index_est]

# Observed response
y_obs <- clim$log_rain

# Residuals
residuals <- y_obs - fitted_vals

# Attach residuals to dataset
clim$residuals <- residuals

###############################################################################
# 2. CREATE SPATIAL NEIGHBOUR STRUCTURE
###############################################################################

coords <- as.matrix(clim[, c("Longitude", "Latitude")])

# 5-nearest neighbour structure (robust for irregular stations)
knn_nb <- knearneigh(coords, k = 5)

nb <- knn2nb(knn_nb)

# Spatial weights
lw <- nb2listw(nb, style = "W")

###############################################################################
# 3. GLOBAL MORAN'S I TEST
###############################################################################

moran_test <- moran.test(clim$residuals, lw)

cat("\n====================================================\n")
cat("GLOBAL MORAN'S I TEST FOR MODEL RESIDUALS\n")
cat("====================================================\n")

print(moran_test)

###############################################################################
# 4. LOCAL MORAN'S I (LISA)
###############################################################################

local_moran <- localmoran(clim$residuals, lw)

clim$local_I <- local_moran[,1]
clim$local_p <- local_moran[,5]

###############################################################################
# 5. CONVERT DATA TO SPATIAL OBJECT
###############################################################################

clim_sf <- st_as_sf(
  clim,
  coords = c("Longitude","Latitude"),
  crs = 4326
)

###############################################################################
# 6. LOCAL MORAN'S I MAP
###############################################################################

ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_sf(data = clim_sf,
          aes(color = local_I),
          size = 2) +
  scale_color_viridis(option = "plasma", direction = -1) +
  labs(
    #title = "Local Moran's I Spatial Autocorrelation of Residuals",
    #subtitle = "Bayesian Spatiotemporal Rainfall Model – South Africa",
    color = "Local Moran's I"
  ) +
  theme_minimal()


###############################################################################
# 7. SIGNIFICANT SPATIAL CLUSTERS
###############################################################################

clim_sf$cluster <- ifelse(clim_sf$local_p < 0.05, "Significant", "Not Significant")

ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_sf(data = clim_sf,
          aes(color = cluster),
          size = 2) +
  scale_color_manual(values = c("red","grey40")) +
  labs(
    #title = "Significant Residual Spatial Clusters",
    #subtitle = "Local Moran's I (p < 0.05)"
  ) +
  theme_minimal()


