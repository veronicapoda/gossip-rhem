################################################################################
################### Plot Functions for Simulation Results ######################
################################################################################


# ------------------------------------------------------------------------------
# Function: plot_box()
# Description: creates boxplots comparing estimated betas from censored and 
# uncensored models
# ------------------------------------------------------------------------------
plot_box <- function(betas_cens_list, betas_unc_list, title_text, beta_true = NULL) {
  
  # Helper: convert list of beta vectors to tidy data frame
  make_df <- function(betas_list, model_type) {
    do.call(rbind, lapply(betas_list, function(beta_vec) {
      tibble(term = names(beta_vec), estimate = as.numeric(beta_vec))
    })) %>%
      mutate(model = model_type)
  }
  
  #  Prepare data
  df_cens <- make_df(betas_cens_list, "censored")
  df_unc  <- make_df(betas_unc_list, "uncensored")
  df_all <- bind_rows(df_cens, df_unc)
  df_all <- df_all %>%
    mutate(term_model = interaction(term, model, sep = " - "))
  
  # Plot boxplots
  p <- ggplot(df_all, aes(x = term_model, y = estimate, fill = model)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    labs(title = title_text, x = "girl alter", y = "estimated \u03B2") +
    theme_minimal() +
    scale_fill_manual(values = c("censored" = "grey", "uncensored" = "grey")) +
    theme(plot.title = element_text(hjust = 0.5, size = 18), 
          legend.position = "none",
          axis.text.x = element_text(size = 18),  # larger text labels
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18, margin = margin(t = 10)), 
          axis.title.y = element_text(size = 18, margin = margin(r = 10)),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    scale_x_discrete(labels = function(x) {
      gsub(".*-", "", x) |> trimws() # mostra solo "Censored"/"Uncensored"
    })
  
  # Add true beta reference line if provided
  if (!is.null(beta_true)) {
    beta_df <- tibble(term = names(beta_true), true_val = beta_true)
    p <- p + geom_hline(data = beta_df, aes(yintercept = true_val),
                        linetype = "dashed", color = "black", size = 1)
  }
  
  return(p)
}











# ------------------------------------------------------------------------------
# Function: plot_smooth()
# Description: plots smoothed effects (s(age)) from multiple GAM models,
# comparing censored vs. uncensored fits. Includes 95% confidence intervals
# and true function reference
# ------------------------------------------------------------------------------
plot_smooth <- function(gam_fit_cens_list, gam_fit_unc_list, covariates_df, title_text) {
  library(ggplot2)
  library(dplyr)
  
  # Sequence of age values
  age_seq <- seq(min(covariates_df$age), max(covariates_df$age), length.out = 100)
  ref_row <- covariates_df[1, , drop = FALSE]
  
  # Prepare new data for prediction
  new_data <- ref_row[rep(1, length(age_seq)), ]
  new_data$age <- age_seq
  
  # Helper: compute mean + 95% CI for a GAM smooth term
  get_smooth_ci <- function(model) {
    pred <- predict(model, newdata = new_data, type = "terms", terms = "s(age)", se.fit = TRUE)
    fit <- pred$fit[, 1]     
    se  <- pred$se.fit[, 1]
    tibble(
      age = age_seq,
      mean = fit,
      lower = fit - 1.96 * se,  # 95% CI
      upper = fit + 1.96 * se
    )
  }
  
  # Apply to all models
  cens_df_list <- lapply(gam_fit_cens_list, get_smooth_ci)
  unc_df_list  <- lapply(gam_fit_unc_list, get_smooth_ci)
  
  # Average across replicates
  average_df <- function(df_list, label) {
    mat_mean  <- sapply(df_list, function(x) x$mean)
    mat_lower <- sapply(df_list, function(x) x$lower)
    mat_upper <- sapply(df_list, function(x) x$upper)
    tibble(
      age = age_seq,
      mean = rowMeans(mat_mean, na.rm = TRUE),
      lower = rowMeans(mat_lower, na.rm = TRUE),
      upper = rowMeans(mat_upper, na.rm = TRUE),
      type = label
    )
  }
  
  df_cens <- average_df(cens_df_list, "censored")
  df_unc  <- average_df(unc_df_list,  "uncensored")
  df_all  <- bind_rows(df_cens, df_unc)
  
  # True function (defined in simulation code)
  true_df <- data.frame(age = age_seq, value = true_function(age_seq))
  
  # Plot smooths
  ggplot_object <- ggplot() +
    geom_ribbon(data = df_all, aes(x = age, ymin = lower, ymax = upper, fill = type), alpha = 0.25) +
    geom_line(data = df_all, aes(x = age, y = mean, color = type), linewidth = 1) +
    geom_line(data = true_df, aes(x = age, y = value), linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = c("censored" = "salmon", "uncensored" = "deepskyblue")) +
    scale_fill_manual(values = c("censored" = "salmon", "uncensored" = "deepskyblue")) +
    labs(title = title_text, x = "age", y = "log-hazard contribution",
         color = NULL, fill = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18, margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, margin = margin(r = 10)),
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  return(ggplot_object)
  
}









# ------------------------------------------------------------------------------
# Function: plot_beta_boxplot()
# Description: simple boxplot for estimated betas with a global coefficient 
# reference line
# ------------------------------------------------------------------------------
plot_beta_boxplot <- function(beta_df, global_coef) {
  ggplot(beta_df, aes(x = model, y = beta)) +
    geom_boxplot(fill = 'grey', alpha = 0.6, outlier.shape = NA) +
    geom_hline(yintercept = global_coef, linetype = "dashed", color = "black", size = 1) +
    labs(
      y = "estimated \u03B2",
      x = ""
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "none",
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18, margin = margin(r = 15)),
      axis.title.y = element_text(size = 18, margin = margin(r = 15))
    )
}
