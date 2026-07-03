###############################################################################
# LOAD LIBRARIES
###############################################################################
###############################################################################
# BAYESIAN SPATIOTEMPORAL CLIMATE MODELLING – SOUTH AFRICA
###############################################################################

library(ggpubr)
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
library(fmesher)

inla.setOption(num.threads = "6:2")


set.seed(123)

load("res.RData")

###############################################################################
# 1. DATA PREPARATION
###############################################################################

#clim.orig <- read.csv("rainfall.temp SA data.csv")
#dim(clim.orig)
#head(clim.orig)
#
#pontos.tab <- aggregate(cbind(Latitude, Longitude) ~ district + Province, data=clim.orig, FUN=mean)
#head(pontos.tab)
#
#write.csv(pontos.tab, "estation.csv")
#


clim <- read.csv("rainfall.temp SA data.csv")


head(clim)

length(with(clim, which(is.na(rainfall))))

dim(clim)


clim <- clim %>%
  filter(
    !is.na(rainfall),
    !is.na(temp.max), !is.na(temp.min),
    !is.na(Longitude), !is.na(Latitude)
  ) %>%
  mutate(
    year_idx  = year - min(year) + 1L,   # índice anual: 1 a 43
    month_idx = as.integer(Months),      # índice sazonal: 1 a 12
    tmax      = as.numeric(scale(temp.max)),
    tmin      = as.numeric(scale(temp.min))
  )

head(clim)

coords <- as.matrix(clim[, c("Longitude", "Latitude")])

head(clim[, c("Longitude", "Latitude")])

head(coords)

###############################################################################
# 2. SPDE MESH (STABLE SETTINGS)
###############################################################################

#mesh <- inla.mesh.2d(
#  loc = coords,
#  max.edge = c(2, 8),
#  cutoff = 0.5
#)

amplitude_x <- diff(range(coords[,1])) # Longitude
amplitude_y <- diff(range(coords[,2])) # Latitude
maior_dimensao <- max(amplitude_x, amplitude_y)

max.edge.interno <- maior_dimensao * 0.05 
max.edge.interno

max.edge.externo <- max.edge.interno * 3
max.edge.externo

cutoff_ideal <- max.edge.interno * 0.1
cutoff_ideal

mesh <- fm_mesh_2d_inla(
  loc = coords, 
  max.edge = c(max.edge.interno, max.edge.externo),
  cutoff = cutoff_ideal
)

mesh
mesh$n

pdf("mesh.pdf", height=10, width=10)
plot(mesh)
points(coords, col = "red", pch = 19, cex = 1.5)
axis(1, at = seq(13, 36, by = 5), labels = paste0(seq(13, 36, by = 5), "°E"))
axis(2, at = seq(-38, -20, by = 4), labels = paste0(abs(seq(-38, -20, by = 4)), "°S"))

dev.off()

spde <- inla.spde2.pcmatern(
  mesh        = mesh,
  alpha       = 2,
  prior.range = c(11, 0.05),  # P(range < 11°) = 0.05
  prior.sigma = c(1,  0.05)   # P(sigma > 1)   = 0.05
)

n_years <- max(clim$year_idx)   # 43
n_years

# Para M3 (Tipo IV): A e sidx com grupo anual
A    <- inla.spde.make.A(mesh, loc = coords, group = clim$year_idx)
sidx <- inla.spde.make.index("spatial.field",
                              n.spde  = spde$n.spde,
                              n.group = n_years)

# Para M0, M1, M2: A e sidx sem grupo
A_m2    <- inla.spde.make.A(mesh, loc = coords)
sidx_m2 <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

###############################################################################
# 3. INLA STACKS
###############################################################################

# Stack principal — M3 (Tipo IV): campo espacial com grupo anual
stack <- inla.stack(
  data    = list(y = clim$rainfall),
  A       = list(A, 1),
  effects = list(
    sidx,
    data.frame(
      intercept = 1,
      year_idx  = clim$year_idx,
      month_idx = clim$month_idx,
      tmax      = clim$tmax,
      tmin      = clim$tmin
    )
  ),
  tag = "est"
)

# Stack M2 — M0, M1, M2: sem grupo temporal no campo espacial
stack_m2 <- inla.stack(
  data    = list(y = clim$rainfall),
  A       = list(A_m2, 1),
  effects = list(
    sidx_m2,
    data.frame(
      intercept = 1,
      year_idx  = clim$year_idx,
      month_idx = clim$month_idx,
      tmax      = clim$tmax,
      tmin      = clim$tmin
    )
  ),
  tag = "est_m2"
)

###############################################################################
# 4. COMPETING BAYESIAN MODELS
###############################################################################

# M0: Baseline — apenas efeitos fixos
f_m0 <- y ~ -1 + intercept + tmax + tmin

# M1: Temporal — RW1 interanual + RW2 sazonal cíclico
f_m1 <- y ~ -1 + intercept + tmax + tmin +
  f(year_idx,  model = "rw1", scale.model = TRUE) +
  f(month_idx, model = "rw2", cyclic = TRUE, scale.model = TRUE)

# M2: Espacial — campo SPDE + componentes temporais (sem interação)
f_m2 <- y ~ -1 + intercept + tmax + tmin +
  f(year_idx,      model = "rw1", scale.model = TRUE) +
  f(month_idx,     model = "rw2", cyclic = TRUE, scale.model = TRUE) +
  f(spatial.field, model = spde)

# M3 (PROPOSTO): Tipo IV — Q_delta = Q_SPDE ⊗ Q_RW1
f_m3 <- y ~ -1 + intercept + tmax + tmin +
  f(year_idx,      model = "rw1", scale.model = TRUE) +
  f(month_idx,     model = "rw2", cyclic = TRUE, scale.model = TRUE) +
  f(spatial.field,
    model         = spde,
    group         = spatial.field.group,
    control.group = list(model = "rw1"))

###############################################################################
###############################################################################
# INLA CONTROL SETTINGS (STABLE & COMPARABLE)
###############################################################################

ctrl_inla_cv <- list(
  strategy     = "simplified.laplace",
  int.strategy = "eb",
  diagonal     = 1e-4,
  tolerance    = 1e-3
)
 
ctrl_comp <- list(
  dic     = TRUE,
  waic    = TRUE,
  cpo     = TRUE,
  mlik    = TRUE,
  config  = TRUE
)

ctrl_inla <- list(
  strategy     = "simplified.laplace",
  int.strategy = "eb",
  diagonal     = 1e-5,
  tolerance    = 1e-6
)

ctrl_pred_m3 <- list(
  A       = inla.stack.A(stack),
  compute = TRUE,
  link    = 1
)

ctrl_pred_m2 <- list(
  A       = inla.stack.A(stack_m2),
  compute = TRUE,
  link    = 1
)

###############################################################################
# FIT ALL MODELS
###############################################################################

fit_m0 <- inla(
  formula           = f_m0,
  data              = inla.stack.data(stack_m2),
  family            = "tweedie",
  control.predictor = list(A = inla.stack.A(stack_m2), compute = TRUE), 
  control.compute   = ctrl_comp,
  control.inla      = ctrl_inla,
  verbose           = FALSE
)



fit_m1 <- inla(
  formula           = f_m1,
  data              = inla.stack.data(stack_m2),
  family            = "tweedie",
  control.predictor = ctrl_pred_m2,
  control.compute   = ctrl_comp,
  control.inla      = ctrl_inla,
  verbose           = FALSE
)

fit_m2 <- inla(
  formula           = f_m2,
  data              = inla.stack.data(stack_m2),
  family            = "tweedie",
  control.predictor = ctrl_pred_m2,
  control.compute   = ctrl_comp,
  control.inla      = ctrl_inla,
  verbose           = FALSE
)


fit_m3 <- inla(
  formula           = f_m3,
  data              = inla.stack.data(stack),
  family            = "tweedie",
  control.predictor = ctrl_pred_m3,
  control.compute   = ctrl_comp,
  control.inla      = ctrl_inla,
  verbose           = FALSE
)



# save.image("res.RData")

# Iniciando a partir daqui. !!!!

fits <- list(
  M0_Baseline = fit_m0,
  M1_Temporal = fit_m1,
  M2_Spatial  = fit_m2,
  M3_Proposed = fit_m3
)

fit_proposed <- fit_m3

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
# CROSS-VALIDATION SETUP
# Folds independentes — correr um de cada vez
# Após cada fold: save.image("cv_fold_k.RData")
###############################################################################

set.seed(123)
K       <- 5
n       <- nrow(clim)
fold_id <- sample(rep(1:K, length.out = n))

# Função auxiliar partilhada pelos 5 folds
cv_run_fold <- function(k, clim, fold_id, mesh, spde, n_years,
                        f_m2, f_m3, ctrl_inla) {

  cat("\n", strrep("=", 60), "\n")
  cat("CV FOLD:", k, "\n")
  cat(strrep("=", 60), "\n")

  train_idx    <- which(fold_id != k)
  test_idx     <- which(fold_id == k)
  train        <- clim[train_idx, ]
  test         <- clim[test_idx,  ]
  coords_train <- as.matrix(train[, c("Longitude", "Latitude")])
  coords_test  <- as.matrix(test[,  c("Longitude", "Latitude")])

  # Verificar anos representados no treino
  anos_faltam <- setdiff(1:n_years, unique(train$year_idx))
  if (length(anos_faltam) > 0)
    cat("  Aviso: anos sem obs no treino:", anos_faltam, "\n")

  # ── M2: sem grupo ────────────────────────────────────────────────────────
  cat("  Fitting M2...\n")
  sidx_cv_m2 <- inla.spde.make.index("spatial.field", spde$n.spde)
  A_train_m2 <- inla.spde.make.A(mesh, loc = coords_train)
  A_test_m2  <- inla.spde.make.A(mesh, loc = coords_test)

  stk_cv_m2 <- inla.stack(
    inla.stack(
      data    = list(y = train$rainfall),
      A       = list(A_train_m2, 1),
      effects = list(sidx_cv_m2,
                     data.frame(intercept = 1,
                                year_idx  = train$year_idx,
                                month_idx = train$month_idx,
                                tmax      = train$tmax,
                                tmin      = train$tmin)),
      tag = "train"
    ),
    inla.stack(
      data    = list(y = NA),
      A       = list(A_test_m2, 1),
      effects = list(sidx_cv_m2,
                     data.frame(intercept = 1,
                                year_idx  = test$year_idx,
                                month_idx = test$month_idx,
                                tmax      = test$tmax,
                                tmin      = test$tmin)),
      tag = "test"
    )
  )

  fit_cv_m2 <- tryCatch(
    inla(f_m2,
         data              = inla.stack.data(stk_cv_m2),
         family            = "tweedie",
         control.predictor = list(A       = inla.stack.A(stk_cv_m2),
                                  compute = TRUE, link = 1),
         control.compute   = list(dic = FALSE, waic = FALSE,
                                  cpo = FALSE, config = TRUE),
         control.inla      = ctrl_inla,
         verbose           = FALSE),
    error = function(e) {
      cat("  M2 falhou no fold", k, ":", conditionMessage(e), "\n")
      NULL
    }
  )

  pred_m2 <- if (!is.null(fit_cv_m2)) {
    idx <- inla.stack.index(stk_cv_m2, tag = "test")$data
    fit_cv_m2$summary.fitted.values$mean[idx]
  } else rep(NA, nrow(test))

  # ── M3: com grupo anual (Tipo IV) ────────────────────────────────────────
  cat("  Fitting M3 (Type IV)...\n")
  sidx_cv_m3 <- inla.spde.make.index("spatial.field",
                                      n.spde  = spde$n.spde,
                                      n.group = n_years)
  A_train_m3 <- inla.spde.make.A(mesh, loc = coords_train,
                                  group = train$year_idx)
  A_test_m3  <- inla.spde.make.A(mesh, loc = coords_test,
                                  group = test$year_idx)

  stk_cv_m3 <- inla.stack(
    inla.stack(
      data    = list(y = train$rainfall),
      A       = list(A_train_m3, 1),
      effects = list(sidx_cv_m3,
                     data.frame(intercept = 1,
                                year_idx  = train$year_idx,
                                month_idx = train$month_idx,
                                tmax      = train$tmax,
                                tmin      = train$tmin)),
      tag = "train"
    ),
    inla.stack(
      data    = list(y = NA),
      A       = list(A_test_m3, 1),
      effects = list(sidx_cv_m3,
                     data.frame(intercept = 1,
                                year_idx  = test$year_idx,
                                month_idx = test$month_idx,
                                tmax      = test$tmax,
                                tmin      = test$tmin)),
      tag = "test"
    )
  )

  fit_cv_m3 <- tryCatch(
    inla(f_m3,
         data              = inla.stack.data(stk_cv_m3),
         family            = "tweedie",
         control.predictor = list(A       = inla.stack.A(stk_cv_m3),
                                  compute = TRUE, link = 1),
         control.compute   = list(dic = FALSE, waic = FALSE,
                                  cpo = FALSE, config = TRUE),
         control.inla      = ctrl_inla,
         verbose           = FALSE),
    error = function(e) {
      cat("  M3 falhou no fold", k, ":", conditionMessage(e), "\n")
      NULL
    }
  )

  pred_m3 <- if (!is.null(fit_cv_m3)) {
    idx <- inla.stack.index(stk_cv_m3, tag = "test")$data
    fit_cv_m3$summary.fitted.values$mean[idx]
  } else rep(NA, nrow(test))

  # ── Erros na escala original (mm) ────────────────────────────────────────
  obs <- test$rainfall

  result <- data.frame(
    Fold    = k,
    RMSE_M2 = sqrt(mean((obs - pred_m2)^2, na.rm = TRUE)),
    RMSE_M3 = sqrt(mean((obs - pred_m3)^2, na.rm = TRUE)),
    MAE_M2  = mean(abs(obs - pred_m2), na.rm = TRUE),
    MAE_M3  = mean(abs(obs - pred_m3), na.rm = TRUE)
  )

  cat("  RMSE M2:", round(result$RMSE_M2, 4),
      "| RMSE M3:", round(result$RMSE_M3, 4), "\n")
  cat("  MAE  M2:", round(result$MAE_M2,  4),
      "| MAE  M3:", round(result$MAE_M3,  4), "\n")

  result
}

###############################################################################
# FOLD 1
###############################################################################

cv_result_fold1 <- cv_run_fold(1, clim, fold_id, mesh, spde, n_years,
                                f_m2, f_m3, ctrl_inla_cv)
save.image("cv_fold1.RData")
 
###############################################################################
# FOLD 2
###############################################################################
 
cv_result_fold2 <- cv_run_fold(2, clim, fold_id, mesh, spde, n_years,
                                f_m2, f_m3, ctrl_inla_cv)
save.image("cv_fold2.RData")

 
###############################################################################
# FOLD 3
###############################################################################
 
cv_result_fold3 <- cv_run_fold(3, clim, fold_id, mesh, spde, n_years,
                                f_m2, f_m3, ctrl_inla_cv)

cv_result_fold3

save.image("cv_fold3.RData")

 
###############################################################################
# FOLD 4
###############################################################################
 
cv_result_fold4 <- cv_run_fold(4, clim, fold_id, mesh, spde, n_years,
                                f_m2, f_m3, ctrl_inla_cv)
cv_result_fold4
save.image("cv_fold4.RData")

 
###############################################################################
# FOLD 5
###############################################################################
 
cv_result_fold5 <- cv_run_fold(5, clim, fold_id, mesh, spde, n_years,
                                f_m2, f_m3, ctrl_inla_cv)

cv_result_fold5

save.image("cv_fold5.RData")
 
###############################################################################
# CONSOLIDAR RESULTADOS
# Correr apenas após todos os folds concluídos
###############################################################################
 
cv_results <- rbind(
  cv_result_fold1,
  cv_result_fold2,
  cv_result_fold3,
  cv_result_fold4,
  cv_result_fold5
)
 
cv_summary <- cv_results %>%
  summarise(
    Mean_RMSE_M2 = mean(RMSE_M2, na.rm = TRUE),
    Mean_RMSE_M3 = mean(RMSE_M3, na.rm = TRUE),
    Mean_MAE_M2  = mean(MAE_M2,  na.rm = TRUE),
    Mean_MAE_M3  = mean(MAE_M3,  na.rm = TRUE)
  )
 
cat("\n", strrep("=", 60), "\n")
cat("CROSS-VALIDATION: M2 vs M3 (Tweedie, escala mm)\n")
cat(strrep("=", 60), "\n")
print(cv_results)
print(cv_summary)
 
save.image("cv_final.RData")
 
###############################################################################
# READ SOUTH AFRICA SHAPEFILE
###############################################################################
 
sa_shape <- st_read(paste(getwd(),"/zaf_admin_boundaries/zaf_admin2.shp",sep="")) |> 
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
  fit_proposed$summary.random$spatial.field$mean[1:spde$n.spde]
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

aggregate(Y ~ Temp, data=dat.graf1, FUN=function(x) c("Mean"=mean(x), "Min"=min(x), "Max"=max(x)))
 
graf1 <- xyplot(Y ~ factor(year), groups=Temp, data=dat.graf1, type="l",
								col=c("red", "blue"), 
                ylim=c(0, 30), lwd=2, scales = list(x = list(at = which(unique(dat.graf1$year) %in% year_breaks), labels = year_breaks,
                                                              rot = 45))) 
 
graf1

head(trend_data)
with(trend_data, c("Mean"=mean(rain), "Min"=min(rain), "Max"=max(rain)))

graf2 <- barchart(rain ~ factor(year), data=trend_data, ylab=list("Rainfall (mm)", cex=2), ylim=c(0,80), col="#ECFF7A",
                  xlab=list("Year", cex=2),
									scales = list(x = list(at = which(unique(dat.graf1$year) %in% year_breaks), labels = year_breaks,
                                                              rot = 45, cex=2))) 
 
graf2

pdf("maxmim.pdf", width=16, height=12)
update(doubleYScale(graf2, graf1, add.ylab2 = TRUE, use.style = FALSE,        
                    text = c("Maximum", "Minimum"), cex = 2, columns = 2),
       par.settings = simpleTheme(col.line = c('red', 'blue')),
       ylab.right = list(label = expression(Temperature ~ (degree * C)), cex = 2))
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

ggsave(
  filename = "RainfallEffect.pdf", 
  plot = p_rain,
  width = 8,    
  height = 6,  
  dpi = 300   
)
 
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
  labs(x = "Longitude", y = "Latitude") 
  gis_theme()

p_rain

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
  labs(x = "Longitude", y = "Latitude") 
  
p_inter

ggsave(
  filename = "TempINT.pdf", 
  plot = p_inter,
  width = 8,    
  height = 6,  
  dpi = 300   
)
 

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
  geom_point(data = clim, aes(Longitude, Latitude, color = temp.max), size = 3) +
  scale_color_viridis_c(
    option = "inferno",
    direction = -1,   # Darker at top
    name = "Max Temp (°C)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
	theme_minimal()

p_tmax
 
p_tmin <- ggplot() +
  geom_sf(data = sa_shape, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = clim, aes(Longitude, Latitude, color = temp.min), size = 3) +
  scale_color_viridis_c(
    option = "mako",
    direction = 1,   # Darker at top
    name = "Min Temp (°C)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
	theme_minimal()
 
p_tmin
  
ggsave(
  filename = "Maximum_Temp1.pdf", 
	plot = p_tmax,
  width = 8,    
  height = 6,  
  dpi = 300   
)
 
 ggsave(
  filename = "Minimum_Temp1.pdf", 
	plot = p_tmin,
  width = 8,    
  height = 6,  
  dpi = 300   
)
 
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
  geom_point(aes(Longitude, Latitude, color = rain), size = 3) +
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
  theme_minimal()
 
 
ggsave(
  filename = ".pdf", 
  width = 8,    
  height = 6,  
  dpi = 300   
)
 

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
  geom_point(aes(Longitude, Latitude, color = trend), size = 3) +
  scale_color_viridis_c(
    option = "E",
    direction = -1,   # ← Makes top darker, bottom lighter
    name = "Temp Trend (°C/year)"
  ) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()
 
  
ggsave(
  filename = ".pdf", 
  width = 8,    
  height = 6,  
  dpi = 300   
)
 
###############################################################################
# 5. PIT HISTOGRAM – M3 DIAGNOSTIC
###############################################################################
 
m3 <- fit_m3
 
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
 
 ggsave(
  filename = "PIT_plot.pdf", 
  width = 8,    
  height = 6,  
  dpi = 300   
)
 
###############################################################################
# 1. EXTRACT FITTED VALUES CORRESPONDING TO OBSERVED DATA
###############################################################################


# Get indices for estimation stack
index_est <- inla.stack.index(stack, tag = "est")$data
 
# Extract fitted values for observations only
fitted_vals <- fit_proposed$summary.fitted.values$mean[index_est]
 
# Observed response (escala original mm)
y_obs <- clim$rainfall
 
# Residuals
residuals <- y_obs - fitted_vals
 
# Attach residuals to dataset
clim$residuals <- residuals
 
###############################################################################
# 2. CREATE SPATIAL NEIGHBOUR STRUCTURE
###############################################################################
 
# Coordenadas únicas das 28 estações
# knearneigh requer coordenadas únicas — usar todas as linhas causa
# o aviso "identical points found" pois cada estação tem múltiplas
# observações temporais com as mesmas coordenadas
coords_unique <- clim %>%
  group_by(Longitude, Latitude) %>%
  summarise(.groups = "drop") %>%
  arrange(Longitude, Latitude)
 
coords_mat <- as.matrix(coords_unique[, c("Longitude", "Latitude")])
cat("Número de estações únicas:", nrow(coords_mat), "\n")
 
# 5-nearest neighbour structure sobre estações únicas
knn_nb <- knearneigh(coords_mat, k = 5)
nb     <- knn2nb(knn_nb)
lw     <- nb2listw(nb, style = "W")
 
###############################################################################
# 3. GLOBAL MORAN'S I TEST
###############################################################################
 
# Resíduo médio por estação — Moran's I testa dependência espacial
# entre locais, não entre observações temporais do mesmo local
resid_by_station <- clim %>%
  group_by(Longitude, Latitude) %>%
  summarise(mean_resid = mean(residuals, na.rm = TRUE), .groups = "drop") %>%
  arrange(Longitude, Latitude)
 
moran_test <- moran.test(resid_by_station$mean_resid, lw)
 
moran_test
 
###############################################################################
# 4. LOCAL MORAN'S I (LISA)
###############################################################################
 
local_moran <- localmoran(resid_by_station$mean_resid, lw)
 
resid_by_station$local_I <- local_moran[, 1]
resid_by_station$local_p <- local_moran[, 5]
 
###############################################################################
# 5. CONVERT DATA TO SPATIAL OBJECT
###############################################################################
 
resid_sf <- st_as_sf(
  resid_by_station,
  coords = c("Longitude", "Latitude"),
  crs    = 4326
)
 
###############################################################################
# 6. LOCAL MORAN'S I MAP
###############################################################################
 
moran.plot <- ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_sf(data = resid_sf,
          aes(color = local_I),
          size = 4) +
  scale_color_viridis(option = "plasma", direction = -1) +
  labs(
    color = "Local Moran's I"
  ) +
  theme_minimal()

ggsave(
  filename = "Moran_I_test.pdf", 
  plot = moran.plot,
  width = 8,           # Largura em polegadas
  height = 6,          # Altura em polegadas (ajuste conforme o formato do seu mapa)
  dpi = 300            # Resolução de nível de publicação (300 DPI)
)
###############################################################################
# 7. SIGNIFICANT SPATIAL CLUSTERS
###############################################################################
 
resid_sf$cluster <- ifelse(resid_sf$local_p < 0.05,
                           "Significant", "Not Significant")
 
ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_sf(data = resid_sf,
          aes(color = cluster),
          size = 4) +
  scale_color_manual(values = c("red", "grey40")) +
  labs(
    color = "Cluster"
  ) +
  theme_minimal()



fit_m3$summary.fixed


data.frame(
  Model = names(fits),
  DIC   = sapply(fits, function(x) x$dic$dic),
  WAIC  = sapply(fits, function(x) x$waic$waic),
  Marginal_LogLik = sapply(fits, function(x) x$mlik[1,1]),
  Mean_CPO = sapply(fits, function(x) mean(-log(x$cpo$cpo), na.rm=TRUE))
)



### Estatistica descritiva

# Rainfall
summary(clim$rainfall)
sd(clim$rainfall)
sum(clim$rainfall == 0) / nrow(clim) * 100  # % zeros

# Tmax
summary(clim$temp.max)
sd(clim$temp.max, na.rm=TRUE)

# Tmin
summary(clim$temp.min)
sd(clim$temp.min, na.rm=TRUE)

clim %>%
  group_by(Province) %>%
  summarise(
    Mean_rain = mean(rainfall, na.rm=TRUE),
    Mean_tmax = mean(temp.max, na.rm=TRUE),
    Mean_tmin = mean(temp.min, na.rm=TRUE)
  )


# Mann-Kendall

library("Kendall")
library("trend")

# Mann-Kendall por estação para Tmax e Tmin
mk_results <- clim %>%
  group_by(district, Province, Longitude, Latitude) %>%
  arrange(year, Months) %>%
  summarise(
    MK_tmax_tau   = MannKendall(temp.max)$tau,
    MK_tmax_pval  = MannKendall(temp.max)$sl,
    MK_tmin_tau   = MannKendall(temp.min)$tau,
    MK_tmin_pval  = MannKendall(temp.min)$sl,
    Sen_tmax      = sens.slope(temp.max)$estimates,
    Sen_tmin      = sens.slope(temp.min)$estimates,
    .groups = "drop"
  )

	mk_results


library(ggplot2)
library(sf)
mk_results <- mk_results %>%
  mutate(
    sig_tmax    = ifelse(MK_tmax_pval < 0.05, "p < 0.05", "p >= 0.05"),
    sig_tmin    = ifelse(MK_tmin_pval < 0.05, "p < 0.05", "p >= 0.05"),
    Sen_tmax_yr = Sen_tmax * 12,
    Sen_tmin_yr = Sen_tmin * 12
  )


mk_results


# Adicionar coordenadas dos não-significativos
nonsig_tmax <- mk_results %>% filter(MK_tmax_pval >= 0.05)
nonsig_tmin <- mk_results %>% filter(MK_tmin_pval >= 0.05)

p1 <- ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_point(data = mk_results,
             aes(x = Longitude, y = Latitude, fill = Sen_tmax_yr),
             shape = 21, size = 4, colour = "black") +
  geom_point(data = nonsig_tmax,
             aes(x = Longitude, y = Latitude),
             shape = 4, size = 3, colour = "black", stroke = 1.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Trend\n(°C/year)") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = "X = non-significant (p > 0.05)", size = 3) +
  labs(x = "Longitude", y = "Latitude",
       title = expression(T[max]~"Sen slope (°C/year)")) +
  theme_minimal()

p2 <- ggplot() +
  geom_sf(data = sa_shape, fill = "grey95", colour = "black") +
  geom_point(data = mk_results,
             aes(x = Longitude, y = Latitude, fill = Sen_tmin_yr),
             shape = 21, size = 4, colour = "black") +
  geom_point(data = nonsig_tmin,
             aes(x = Longitude, y = Latitude),
             shape = 4, size = 3, colour = "black", stroke = 1.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Trend\n(°C/year)") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = "X = non-significant (p > 0.05)", size = 3) +
  labs(x = "Longitude", y = "Latitude",
       title = expression(T[min]~"Sen slope (°C/year)")) +
  theme_minimal()

combined <- (p1 + p2) + plot_layout(guides = "collect")
combined & theme(legend.position = "right")


p1 <- p1 + theme(legend.position = "right")
p2 <- p2 + theme(legend.position = "right")

ggarrange(p1, p2, ncol=1)

#ggsave("MK_trend_map.pdf", width = 12, height = 6)

ggsave(
  filename = "MK_trend_map.pdf", 
#  plot = p_rain,
  width = 8,    
  height = 10,  
  dpi = 300   
)


# Mesh ggplot2


# ── 1. Extrair triângulos do mesh ──────────────────────────────
# Cada triângulo tem 3 vértices — extrair coordenadas
triangles <- mesh$graph$tv  # matriz N_tri x 3 com índices dos vértices
vertices  <- mesh$loc        # matriz N_vert x 2 com coordenadas (lon, lat)

# Construir data.frame com os 3 segmentos de cada triângulo
tri_df <- do.call(rbind, lapply(1:nrow(triangles), function(i) {
  v <- triangles[i, ]
  pts <- vertices[v, , drop = FALSE]
  # Fechar o triângulo repetindo o primeiro ponto
  pts <- rbind(pts, pts[1, ])
  data.frame(
    lon      = pts[, 1],
    lat      = pts[, 2],
    triangle = i
  )
}))

# ── 2. Coordenadas das estações ────────────────────────────────
stations <- data.frame(
  lon = coords[, 1],
  lat = coords[, 2]
)

# ── 3. Figura ─────────────────────────────────────────────────
plot.mesh <- ggplot() +
  # shapefile da Africa do Sul (fundo)
  geom_sf(data = sa_shape, fill = "grey95", colour = "grey60",
          linewidth = 0.3) +
  # triângulos do mesh
  geom_polygon(data = tri_df,
               aes(x = lon, y = lat, group = triangle),
               fill = NA, colour = "grey30", linewidth = 0.15) +
  # domínio interno (bounding box)
  geom_rect(aes(xmin = 13.01, xmax = 35.43,
                ymin = -37.73, ymax = -20.16),
            fill = NA, colour = "black", linewidth = 0.6) +
  # estações
  geom_point(data = stations,
             aes(x = lon, y = lat),
             colour = "red", size = 2.5, shape = 19) +
  # eixos e coordenadas
  scale_x_continuous(
    breaks = seq(14, 34, by = 4),
    labels = paste0(seq(14, 34, by = 4), "°E")
  ) +
  scale_y_continuous(
    breaks = seq(-36, -22, by = 4),
    labels = paste0(abs(seq(-36, -22, by = 4)), "°S")
  ) +
  coord_sf(xlim = c(11, 38), ylim = c(-40, -18)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )


plot.mesh
	ggsave("mesh_ggplot.pdf", width = 9, height = 9)
