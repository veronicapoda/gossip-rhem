################################################################################
############################# SIMULATION STUDY #################################
################################################################################
# 
# This script simulates gossip hyperevents among students in different classes 
# over multiple waves. 
#
# Two simulation experiments are included:
#   1. SIMU 1: Censored vs. uncensored models with smooth and linear effects.
#   2. SIMU 2: Evaluating time-varying covariates from partially observed data: 
#              compare the effects of past, current, and average values.
#
##


# Import some useful functions
base_folder <- getwd() 
source(file.path(base_folder, "utils", "covariates.R"))
source(file.path(base_folder, "gossipdata_simulation", "functions_sim.R"))
source(file.path(base_folder, "gossipdata_simulation", "plot_sim.R"))




# Input parameters
n_classes = 1
n_students_per_class = 8
n_years = 6
s_max = 3       # sender set size

# Simulation settings
num_simulations <- 100


# Construct the Risk Set: generate all possible sender-receiver-wave combinations
result <- construct_risk_set(
  n_classes = n_classes,
  n_students_per_class = n_students_per_class,
  n_years = n_years,
  s_max = s_max
)
Risk_Set <- result$risk_set
gender_map <- result$gender   # Map of students' genders








################################################################################
################## Simulation 1: Censored vs. Uncensored Data ##################
################################################################################

library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)


# Covariate selection
exogenous_covariates <- c("girl_alter", "age")
endogenous_covariates <- NULL 
nonlinear_covariates <- c("age")

# Compute exogenous covariates
girl.alter <- compute_covariates(RiskSet = Risk_Set, gender_map = gender_map, CovariateName = "girl.alter")
age <- compute_covariates(RiskSet = Risk_Set, CovariateName = "age")

Risk_Set$girl_alter <- girl.alter
Risk_Set$age <- age



# True coefficients
beta_0_true <- -8   # Intercept for the simulation

manual_betas <- c(girl_alter = 0.9)

# Match beta values to filtered covariates
exo_cov_linear <- setdiff(exogenous_covariates, nonlinear_covariates)
endo_cov_linear <- setdiff(endogenous_covariates, nonlinear_covariates)
all_covariates_linear <- c(exo_cov_linear, endo_cov_linear)
beta_true <- manual_betas[all_covariates_linear]



# Smooth effect function for age (decreasing function)
true_function <- function(age) {
  val <- 1 / (1 + exp(2 * (age - 16)))  
  val - mean(val)
}





# -------------------
# Run simulations
# -------------------
set.seed(143)
sim_small <- run_simulation_cens(lambda_0_val = 0.25, label = "\u03bb\u2080 small")
sim_large <- run_simulation_cens(lambda_0_val = 0.75, label = "\u03bb\u2080 large")


# Generate boxplots of estimated linear coefficients
p_box_small <- plot_box(sim_small$betas_cens_list, sim_small$betas_unc_list, expression("" * lambda[0] * " small"), beta_true)
p_box_large <- plot_box(sim_large$betas_cens_list, sim_large$betas_unc_list, expression("" * lambda[0] * " large"), beta_true)

all_values <- c(
  unlist(sim_small$betas_cens_list),
  unlist(sim_small$betas_unc_list),
  unlist(sim_large$betas_cens_list),
  unlist(sim_large$betas_unc_list)
)
y_min <- min(all_values)
y_max <- max(all_values)

p_box_small <- p_box_small + coord_cartesian(ylim = c(y_min, y_max))
p_box_large <- p_box_large + coord_cartesian(ylim = c(y_min, y_max))


# Plot smooth effect of age from GAM fits (with 95% CI across simulations)
p_smooth_small <- plot_smooth(sim_small$gam_fit_cens_list, sim_small$gam_fit_count_list, sim_small$covariates_df, expression("" * lambda[0] * " small"))
p_smooth_large <- plot_smooth(sim_large$gam_fit_cens_list, sim_large$gam_fit_count_list, sim_large$covariates_df, expression("" * lambda[0] * " large"))

data_small <- ggplot_build(p_smooth_small)$data[[1]]
data_large <- ggplot_build(p_smooth_large)$data[[1]]

y_min <- min(c(data_small$ymin, data_small$y, data_large$ymin, data_large$y), na.rm = TRUE)
y_max <- max(c(data_small$ymax, data_small$y, data_large$ymax, data_large$y), na.rm = TRUE)

y_limits <- c(y_min, y_max)

p_smooth_small <- p_smooth_small + coord_cartesian(ylim = y_limits)
p_smooth_large <- p_smooth_large + coord_cartesian(ylim = y_limits)

p_box_small    <- p_box_small + theme(plot.margin = margin(10, 20, 10, 10))
p_smooth_small <- p_smooth_small + theme(plot.margin = margin(10, 10, 10, 20))
p_box_large    <- p_box_large + theme(plot.margin = margin(20, 20, 10, 10))
p_smooth_large <- p_smooth_large + theme(plot.margin = margin(20, 10, 10, 20))



# Combine all four plots into a 2x2 panel
final_plot <- (p_box_small | p_smooth_small) /
  (p_box_large | p_smooth_large)

print(final_plot)









################################################################################
######### Simulation 2: Covariate Evaluation: Past, Current, Average ###########
################################################################################

library(mgcv)
library(dplyr)
library(ggplot2)


# Define global covariate and coefficients
global_coef <- c(gc = 0.8)            # True effect of global covariate
global_cov <- function(t) log(t+1)    # Global covariate as function of time


# Define baseline hazard
beta_0_true <- -8
lambda_0 <- function(t) {0.038}

# Initialize data frame to store estimated beta coefficients
beta_df <- data.frame(
  sim = integer(),
  model = character(),
  beta = numeric()
)




# -------------------
# Run simulations
# -------------------
set.seed(222)

for (s in 1:num_simulations) {
  
  cat("Simulation", s, "\n")
  
  # Simulate events
  events <- run_simulation(lambda_0_val = lambda_0,
                           n_years = n_years,
                           RiskSet = Risk_Set,
                           global_cov = global_cov,
                           global_coef = global_coef)
  
  # Extract events as a list with receiver, sender, and integer time
  event_list <- lapply(seq_along(events[[1]]), function(i) {
    list(
      time     = floor(events[[1]][[i]]$time + 1),
      receiver = events[[1]][[i]]$receiver,
      sender   = events[[1]][[i]]$sender
    )
  })
  
  # Create a data frame of events
  event_df <- data.frame(
    Time     = sapply(event_list, `[[`, "time"),
    Receiver = sapply(event_list, `[[`, "receiver"),
    Senders  = I(lapply(event_list, `[[`, "sender"))
  )
  event_df$Senders <- lapply(event_df$Senders, as.integer)
  
   # Aggregate RiskSet into wave-based summary
   agg_cov <- aggregate_cov(Risk_Set, n_years)
  
  # Evaluate global covariate for each wave
  agg_cov$current <- global_cov(agg_cov$wave)       # value at current wave
  agg_cov$past    <- global_cov(agg_cov$wave - 1)   # value at previous wave
  agg_cov$average     <- (global_cov(agg_cov$wave) + global_cov(agg_cov$wave - 1)) / 2     # mean of current and past
  
  # Define binary vector for event occurrence
   z_cens  <- define_vector_z(RiskSet = Risk_Set, EventList = unique(event_list))
  agg_cov$z_cens <- z_cens
  
  # Fit cloglog models for past, current, and average covariates
  fit_past_clog    <- gam(z_cens ~ past,    family = binomial(link = "cloglog"), data = agg_cov)
  fit_current_clog <- gam(z_cens ~ current, family = binomial(link = "cloglog"), data = agg_cov)
  fit_avg_clog     <- gam(z_cens ~ average,     family = binomial(link = "cloglog"), data = agg_cov)
  
  # Store estimated beta coefficients in beta_df
   beta_df <- rbind(beta_df,
                   data.frame(sim = s, model = "past",    beta = as.numeric(coef(fit_past_clog)["past"])),
                   data.frame(sim = s, model = "current", beta = as.numeric(coef(fit_current_clog)["current"])),
                   data.frame(sim = s, model = "average",     beta = as.numeric(coef(fit_avg_clog)["average"]))
  )
}


# Convert model column to factor with specific order
beta_df$model <- factor(beta_df$model, levels = c("past", "current", "average"))


# Plot boxplots of estimated beta coefficients
plot_beta_boxplot(beta_df, global_coef)

