###############################################################################

library(INLA)
library(sf)
library(sp)
library(spdep)
library(tmap)
library(dplyr)
library(tidyverse)
library(mapview)
library(ggplot2)
library(igraph)
library(splines)
library(survival)
library(viridis)
library(patchwork)
###############################################################################
# READ SHAPEFILE AND DATASET
###############################################################################

nig <- st_read("C:/Users/HP/Desktop/Research 2025/Flexible splines/shp/nga_admin1_em.shp")

data <- read.csv("C:/Users/HP/Desktop/Research 2025/Flexible splines/U5M24WMVs.csv")

###############################################################################
# FUNCTION TO SIMULATE WEIBULL-IMPUTED TIME VALUES
# (MODIFIED: now imputes BOTH missing AND zero values)
###############################################################################

simulate_weibull_times <- function(time_col, shape, scale, seed) {
  
  set.seed(seed)
  
  # Identify missing OR zero OR negative times
  bad_idx <- which(is.na(time_col) | time_col <= 0)
  
  imputed_times <- time_col
  
  # Generate Weibull values only for problematic entries
  simulated_vals <- rweibull(length(bad_idx), shape = shape, scale = scale)
  
  # Keep within DHS bounds (0–60 months)
  simulated_vals <- pmin(60, pmax(0.5, round(simulated_vals)))
  
  # Replace missing and zero values
  imputed_times[bad_idx] <- simulated_vals
  
  return(imputed_times)
}

###############################################################################
# IMPUTE MISSING AND ZERO TIMES
###############################################################################

DeleL <- data
DeleL$time <- simulate_weibull_times(data$time, shape = 0.5, scale = 10, seed = 202)

###############################################################################
# ROBUST SURVIVAL TIME CLEANING (UNCHANGED STRUCTURE)
###############################################################################

# Ensure numeric
DeleL$time <- as.numeric(DeleL$time)

# Ensure status is numeric and binary
DeleL$status <- as.numeric(DeleL$status)
DeleL$status[is.na(DeleL$status)] <- 0

#write.csv(DeleL,file = "DeleL.csv",row.names = FALSE)

###############################################################################
# DIAGNOSTIC CHECK
###############################################################################

cat("Minimum survival time:", min(DeleL$time, na.rm = TRUE), "\n")
cat("Any time <= 0?:", any(DeleL$time <= 0), "\n")
cat("Any missing time?:", any(is.na(DeleL$time)), "\n")

###############################################################################
# PIECEWISE EXPONENTIAL MODEL - DATA EXPANSION
###############################################################################

intervals <- c(6, 12, 18, 24, 30, 36, 42, 48, 54)
expanded_dataL<- survSplit(
  Surv(time, status) ~ .,
  data = DeleL,
  cut = intervals,
  episode = "interval"
)

# Calculate duration and add spline terms
expanded_dataL <- expanded_dataL %>%
  group_by(id = row_number()) %>%
  mutate(
    start = lag(time, default = 0),
    stop = time,
    duration = stop - start
  ) %>%
  ungroup() %>%
  filter(duration > 0)

expanded_dataL$time_mid <- (expanded_dataL$start + expanded_dataL$stop) / 2
spline_basis <- ns(expanded_dataL$time_mid, df = 3)
colnames(spline_basis) <- paste0("spline", 1:3)
expanded_dataL<- bind_cols(expanded_dataL, as.data.frame(spline_basis))

# Ensure factors are properly coded
expanded_dataL <- expanded_dataL %>%
  mutate(across(c(c_tf, c_pbi, m_edu, c_dop, c_nav,c_wid, region, interval), as.factor))

###############################################################################
# PREPARE SPATIAL STRUCTURE
#########################################################################################
# Ensure shapefile index
nig$state_id <- 1:nrow(nig)

# Create correct spatial index
expanded_dataL$state1 <- match(expanded_dataL$sstate, nig$state_id)

# Build adjacency
nb <- poly2nb(nig)
nb2INLA("nw.graph", nb)
nw.adj <- inla.read.graph("nw.graph")

# Verify match
stopifnot(max(expanded_dataL$state1, na.rm=TRUE) == length(nb))

#Model 1: Basic Piecewise Exponential
formula1 <- status ~ mab + sex + c_twin + tpr + c_sodw + c_dbf + c_cu + mq_net +
  factor(c_tf )+ factor(c_pbi) + factor(m_edu) + factor(c_dop) + factor(c_nav) + 
  factor(c_wid) + factor(region)

pem1_L <- inla(
  formula1,
  family = "poisson",
  data = expanded_dataL,
  E = expanded_dataL$duration,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 2: With Splines
formula2 <- update(formula1, ~ . + spline1 + spline2 + spline3)

pem2_L <- inla(
  formula2,
  family = "poisson",
  data = expanded_dataL,
  E = expanded_dataL$duration,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 3: With Splines and Interval RE
formula3 <- update(formula2, ~ . + f(interval, model = "iid"))

pem3_L <- inla(
  formula3,
  family = "poisson",
  data = expanded_dataL,
  E = expanded_dataL$duration,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 4: With Splines and Spatial Effects
formula4 <- update(formula2, ~ . + f(state1, model = "besag", graph = nw.adj))

pem4_L <- inla(
  formula4,
  family = "poisson",
  data = expanded_dataL,
  E = expanded_dataL$duration,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 5: With Splines, Interval RE and Spatial Effects
formula5 <- update(formula3, ~ . + f(state1, model = "besag", graph = nw.adj))

pem5_L <- inla(
  formula5,
  family = "poisson",
  data = expanded_dataL,
  E = expanded_dataL$duration,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Compare LOW hazard models
model_compare_L <- data.frame(
  Model = c("Basic", "+Splines", "+Interval RE", "+Spatial", "+Int+Space"),
  DIC = c(pem1_L$dic$dic, pem2_L$dic$dic, pem3_L$dic$dic, 
          pem4_L$dic$dic, pem5_L$dic$dic),
  WAIC = c(pem1_L$waic$waic, pem2_L$waic$waic, pem3_L$waic$waic,
           pem4_L$waic$waic, pem5_L$waic$waic),
  CPO = c(-mean(log(pem1_L$cpo$cpo)), -mean(log(pem2_L$cpo$cpo)),
          -mean(log(pem3_L$cpo$cpo)), -mean(log(pem4_L$cpo$cpo)),
          -mean(log(pem5_L$cpo$cpo)))
) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

# Print LOW hazard results
cat("\nLOW HAZARD MODEL COMPARISON (Weibull shape=0.5, scale=10)\n")
print(model_compare_L)
###############################################################################
# COMBINED DUAL-ROLE INTERVAL EFFECT PLOT
# Shows both baseline hazard and random effect deviation
###############################################################################


# Select model
model <- pem5_L   

###############################################################################
# EXTRACT AND PREPARE INTERVAL EFFECTS
###############################################################################

interval_re <- model$summary.random$interval %>%
  
  mutate(
    interval_id = as.numeric(as.character(ID)),
    
    mean_log  = mean,
    lower_log = `0.025quant`,
    upper_log = `0.975quant`
  ) %>%
  
  arrange(interval_id)

# Define cut points (must match survSplit)
cut_points <- c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60)

interval_re <- interval_re %>%
  
  mutate(
    start = cut_points[interval_id],
    stop  = cut_points[interval_id + 1],
    midpoint = (start + stop)/2,
    label = paste0(start, "–", stop)
  )

###############################################################################
# NORMALIZE HAZARD RATIOS RELATIVE TO REFERENCE
###############################################################################

reference_log <- interval_re$mean_log[1]

interval_re <- interval_re %>%
  
  mutate(
    mean_hr  = exp(mean_log  - reference_log),
    lower_hr = exp(lower_log - reference_log),
    upper_hr = exp(upper_log - reference_log)
  )

###############################################################################
# PANEL A: CONSTANT BASELINE HAZARD FUNCTION
###############################################################################

p1 <- ggplot(interval_re,
             aes(x = midpoint, y = mean_hr)) +
  
  geom_step(
    linewidth = 1.2,
    color = "black"
  ) +
  
  geom_point(size = 2.8) +
  
  geom_errorbar(
    aes(ymin = lower_hr, ymax = upper_hr),
    width = 1.5
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  
  labs(
    #title = "A. Piecewise Constant Baseline Hazard",
    x = "Time (Months)",
    y = "HR"
  ) +
  
  theme_minimal(base_size = 13)

###############################################################################
# PANEL B: RANDOM EFFECT DEVIATION (FRAILTY COMPONENT)
###############################################################################

p2 <- ggplot(interval_re,
             aes(x = midpoint, y = mean_log)) +
  
  geom_step(
    linewidth = 1.2,
    color = "blue"
  ) +
  
  geom_point(
    size = 2.8,
    color = "blue"
  ) +
  
  geom_errorbar(
    aes(ymin = lower_log, ymax = upper_log),
    width = 1.5,
    color = "blue"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  
  labs(
    #title = "B. Interval Random Effect Deviations",
    x = "Time (Months)",
    y = "Log Hazard Deviation"
  ) +
  
  theme_minimal(base_size = 13)

###############################################################################
# COMBINE PANELS 
###############################################################################

combined_plot <- p1 / p2 +
  
  plot_annotation(
    #title = "Dual-Role Parameterization of Interval Effects in Piecewise Exponential Model",
    #subtitle = "Top: constant baseline hazard estimator | Bottom: random effect deviations",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16)
    )
  )

###############################################################################
# DISPLAY
###############################################################################

combined_plot

###############################################################################
# SAVE HIGH-RESOLUTION VERSION
###############################################################################

ggsave(
  "dual_role_interval_effect_plot.png",
  combined_plot,
  width = 8,
  height = 10,
  dpi = 600
)

#SPATIAL ASSESSMENT

# spatial random effects
spatial_effects <- pem5_L $summary.random$state1

# create a mapping between shapefile states and random effect IDs
# assume shapefile is ordered the same as data (geometries in the same 1..N sequence)
# safest: assign ID sequence
nig$state1 <- 1:nrow(nig)

# match random effect means to shapefile
nig$spatial_mean <- spatial_effects$mean[match(nig$state1, spatial_effects$ID)]
# plot
ggplot(nig) +
  geom_sf(aes(fill=spatial_mean), color="white", size=0.3) +
  scale_fill_viridis(
    option = "plasma",
    direction = -1,   # darker = more deaths
    name = "Spatial Frailty effects",
    na.value = "grey90") +
  labs(
    #title="Posterior Spatial Frailty Effects",
    #subtitle="Piecewise Exponential + ICAR Model (Model 5)",
    caption="INLA results"
  ) +
  theme_minimal(base_size=14) +
  theme(
    legend.position="right",
    panel.grid=element_blank(),
    plot.title=element_text(face="bold", size=16),
    plot.subtitle=element_text(size=12)
  )
