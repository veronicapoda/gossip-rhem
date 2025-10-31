################################################################################
############################# SIMULATION FUNCTIONS #############################
################################################################################


# This section contains helper functions for simulation.
# Functions included:
# 1. construct_risk_set       - Builds the full risk set with senders, receivers, and waves
# 2. lambda_cov               - Computes hazard rates including linear, nonlinear, and random effects
# 3. simulate_events          - Simulates hyperevents dynamically over multiple waves
# 4. run_simulation_cens      - Wrapper to run multiple simulations for censored vs count models
# 5. simulate_events_glob     - Simulates events using a global time-varying covariate
# 6. aggregate_cov            - Aggregates RiskSet by discrete waves for covariate evaluation
# 7. run_simulation           - Wrapper to run one full simulation with global covariate
# 8. define_vector_z          - Defines binary indicator vector for observed events
# 9. define_vector_z_count    - Defines count vector for number of times each event occurred





################################################################################
# 1. CONSTRUCT THE RISK SET
################################################################################
# Builds all possible sender-receiver combinations across classes and waves.
# Assigns random genders to students and returns both RiskSet and gender map.

construct_risk_set <- function(n_classes, n_students_per_class, n_years, s_max = 3, start_class_code = 1100) {
  full_risk_set <- list()
  full_gender_map <- list()
  
  for (i in 0:(n_classes - 1)) {
    class_code <- start_class_code + i * 100  # e.g., 1100, 1200, ...
    class_prefix <- class_code %/% 100        # e.g., 11, 12, ...
    student_suffixes <- 1:n_students_per_class
    students_in_class <- as.character(class_prefix * 100 + student_suffixes)
    cat("\nProcessing Class:", class_code, "\n")
    gender_map <- setNames(sample(c(0, 1), length(students_in_class), replace = TRUE), students_in_class)
    full_gender_map <- c(full_gender_map, gender_map)
    for (year in 1:n_years) {
      cat("  Processing Wave:", year, "\n")
      for (num_senders in 2:s_max) {
        sender_combinations <- combn(students_in_class, num_senders, simplify = FALSE)
        for (senders in sender_combinations) {
          possible_receivers <- setdiff(students_in_class, senders)
          for (receiver in possible_receivers) {
            full_risk_set <- append(full_risk_set, list(list(
              senders = senders,
              receiver = receiver,
              time = year,
              class = class_code
            )))
          }
        }
      }
    }
  }
  df <- do.call(rbind, lapply(full_risk_set, function(x) {
    senders <- paste(x$senders, collapse = ", ")
    data.frame(Senders = senders, Receiver = x$receiver, Time = x$time, Class = x$class, stringsAsFactors = FALSE)
  }))
  df$Time <- as.numeric(df$Time)
  df$Receiver <- as.numeric(df$Receiver)
  df$Class <- as.numeric(df$Class)
  df$Senders <- lapply(df$Senders, function(x) as.numeric(trimws(unlist(strsplit(x, ",")))))
  df <- unique(df)
  df <- df[order(df$Time, df$Class), ]
  cat("Total unique sender-receiver-time combinations:", nrow(df), "\n")
  return(list(risk_set = df, gender = full_gender_map))
}






################################################################################
# 2. HAZARD FUNCTION
################################################################################
# Computes event rates (hazard) incorporating:
#  - Linear covariates (exogenous and endogenous)
#  - Nonlinear transformations
#  - Random intercepts and slopes


lambda_cov <- function(lambda_0, beta, RiskSet, history,
                       endogenous_cov, exogenous_cov,
                       nonlinear_cov = NULL,
                       random_intercept_cov = NULL,
                       gamma = NULL,
                       gamma_slope = NULL,
                       random_slope_cov = NULL) {
  
  # Determine current wave/time
  last_time <- tail(history$time, 1)
  wave <- ifelse(length(last_time) == 0, 1, last_time + 1)
  riskset_wave <- RiskSet[RiskSet$Time == wave, ]
  
  ## Linear effects 
  covariates_linear <- NULL
  
  # Exogenous covariates
  if (length(exogenous_cov) > 0) {
    covariates_linear <- as.matrix(riskset_wave[, exogenous_cov, drop = FALSE])
    colnames(covariates_linear) <- exogenous_cov
  }
  
  # Endogenous covariates (computed from history)
  for (cov in endogenous_cov) {
    updated_covariate <- numeric(nrow(riskset_wave))
    for (i in seq_len(nrow(riskset_wave))) {
      event <- riskset_wave[i, ]
      e <- list(time = event$Time, receiver = event$Receiver, sender = event$Senders[[1]])
      updated_covariate[i] <- compute_covariates(CovariateName = cov, History = history, Event = e)
    }
    covariates_linear <- if (is.null(covariates_linear)) {
      matrix(updated_covariate, ncol = 1)
    } else {
      cbind(covariates_linear, updated_covariate)
    }
    colnames(covariates_linear)[ncol(covariates_linear)] <- cov
  }
  
  
  # Non-linear effects
  covariates_nonlinear_transformed <- NULL
  
  if (!is.null(nonlinear_cov) && length(nonlinear_cov) > 0) {
    covariates_nonlinear_transformed <- matrix(NA, nrow = nrow(riskset_wave), ncol = length(nonlinear_cov))
    colnames(covariates_nonlinear_transformed) <- nonlinear_cov
    for (j in seq_along(nonlinear_cov)) {
      cov <- nonlinear_cov[j]
      if (!is.null(covariates_linear) && cov %in% colnames(covariates_linear)) {
        covariates_nonlinear_transformed[, j] <- 1 / (1 + exp(2 * (covariates_linear[, cov] - 16))) 
        covariates_linear <- covariates_linear[, !(colnames(covariates_linear) %in% cov), drop = FALSE]
      } else {
        warning(paste("Covariate", cov, "not found for nonlinear transformation"))
        covariates_nonlinear_transformed[, j] <- NA
      }
    }
  }
  
  # Random effects (intercept and/or slope)
  random_effect_term <- rep(0, nrow(riskset_wave))
  random_effect_covariate_term <- 0
  if (!is.null(random_intercept_cov) && (!is.null(gamma) || !is.null(gamma_slope))) {
    group_ids <- switch(random_intercept_cov,
                        "ClassID" = as.character(floor(riskset_wave$Receiver / 100) * 100),
                        "Receiver" = as.character(riskset_wave$Receiver),
                        "Senders" = sapply(riskset_wave$Senders, function(s) paste(sort(as.character(s)), collapse = "_")),
                        as.character(riskset_wave[[random_intercept_cov]])
    )
    group_factor <- factor(group_ids)
    group_index <- as.integer(group_factor)
    group_levels <- levels(group_factor)
    
    # Random intercept
    if (!is.null(gamma)) {
      if (!all(group_levels %in% names(gamma))) stop("Some group IDs missing in gamma vector.")
      gamma_vector <- gamma[group_levels]
      random_effect_term <- gamma_vector[group_index]
    }
    
    # Random slope
    if (!is.null(gamma_slope)) {
      if (!all(group_levels %in% rownames(gamma_slope))) stop("Some group IDs missing in gamma_slope matrix.")
      if (!is.null(random_slope_cov)) {
        covariates_for_slope <- covariates_nonlinear_transformed[, random_slope_cov, drop = FALSE]
        gamma_slope_matrix <- gamma_slope[group_levels, random_slope_cov, drop = FALSE]
        random_slope <- gamma_slope_matrix[group_index, , drop = FALSE]
      } else {
        covariates_for_slope <- covariates_nonlinear_transformed
        gamma_slope_matrix <- gamma_slope[group_levels, , drop = FALSE]
        random_slope <- gamma_slope_matrix[group_index, , drop = FALSE]
      }
      random_effect_covariate_term <- rowSums(covariates_for_slope * random_slope)
    }
  }
  
  # Compute linear predictor
  relative_risk_function <- rep(0, nrow(riskset_wave))
  if (!is.null(covariates_linear)) {
    relative_risk_function <- relative_risk_function + covariates_linear %*% beta
  }
  if (!is.null(covariates_nonlinear_transformed)) {
    relative_risk_function <- relative_risk_function + covariates_nonlinear_transformed 
  }
  relative_risk_function <- relative_risk_function + random_effect_term + random_effect_covariate_term
  
  if (any(relative_risk_function < -4 | relative_risk_function > 4)) {
    warning("Warning: Linear predictor outside [-4, 4]")
  }
  
  # Compute event rates and return results
  rates <- lambda_0[wave] * exp(relative_risk_function)
  
  # Combine covariates for output
  if (is.null(covariates_linear)) covariates_linear <- matrix(nrow=nrow(riskset_wave), ncol=0)
  if (is.null(covariates_nonlinear_transformed)) covariates_nonlinear_transformed <- matrix(nrow=nrow(riskset_wave), ncol=0)
  
  covariates <- cbind(covariates_linear, covariates_nonlinear_transformed)
  
  return(list(rates = rates, covariates = covariates))
}










################################################################################
# 3. SIMULATE HYPEREVENTS DYNAMICALLY
################################################################################
# Simulates a sequence of hyperevents over multiple waves based on hazard rates.
# Stores event history, covariates, and true hazards for each wave.


simulate_events <- function(RiskSet, lambda_0, beta, n_years, 
                            endogenous_cov, exogenous_cov, 
                            nonlinear_cov = NULL,
                            random_intercept_cov = NULL,
                            gamma = NULL, 
                            gamma_slope = NULL,
                            random_slope_cov = NULL) {
  # Initialize storage
  history <- list(time = numeric(), receiver = numeric(), sender = list())
  covariates_list <- list()
  hazard_true_list <- list()
  t <- 0
  
  pb <- txtProgressBar(min = 0, max = n_years, style = 3)
  cat("Starting simulation...\n")
  
  # Main simulation loop over waves
  for (wave in seq_len(n_years)) {
    # Compute hazard rates and covariates for this wave
    result <- lambda_cov(lambda_0, beta, RiskSet, history,
                         endogenous_cov, exogenous_cov,
                         nonlinear_cov, random_intercept_cov,
                         gamma, gamma_slope, random_slope_cov)
    
    hazard_rates <- result$rates
    hazard_true_list[[wave]] <- hazard_rates
    
    # Compute total event intensity and sampling probabilities
    total_intensity <- sum(hazard_rates)
    event_probabilities <- hazard_rates / total_intensity
    n.atrisk <- length(hazard_rates)
    
    # Simulate inter-event times until the next wave
    while (t < wave) {
      delta_t <- rexp(1, total_intensity)    # waiting time
      
      # If the next event time exceeds the next wave boundary, move to next wave
      if (t + delta_t >= wave + 1) {
        t <- floor(t + delta_t)
      } else {
        t <- t + delta_t
        
        # Stop simulating within this wave if crossing boundary
        if (t >= wave) break 
        
        # Randomly choose one event based on hazard probabilities
        chosen_event <- sample(n.atrisk, size = 1, prob = event_probabilities)
        event_wave <- RiskSet[RiskSet$Time == wave, ][chosen_event, ]
        
        # Record event in history
        event <- list(
          time = t, 
          receiver = event_wave$Receiver, 
          sender = event_wave$Senders[[1]]
          )
        
        history$time <- c(history$time, floor(event$time) + 1)
        history$receiver <- c(history$receiver, as.numeric(event$receiver))
        history$sender <- append(history$sender, list(event$sender))
      }
    }
    
    # Save covariates from lambda_cov() for this wave
    for (cov_name in colnames(result$covariates)) {
      if (is.null(covariates_list[[cov_name]])) {
        covariates_list[[cov_name]] <- result$covariates[, cov_name]
      } else {
        covariates_list[[cov_name]] <- c(covariates_list[[cov_name]], result$covariates[, cov_name])
      }
    }
    setTxtProgressBar(pb, wave)
    if (wave %% 5 == 0 || wave == n_years) {
      cat(sprintf("\rProgress: %d%% (%d/%d waves) | Events simulated: %d      ", 
                  round(wave/n_years * 100), wave, n_years, length(history$time)))
    }
  }
  # Return results
  close(pb)
  cat("\nSimulation completed!\n")
  
  return(list(
    history = history, 
    covariates = covariates_list, 
    hazard_true = hazard_true_list
    ))
}







################################################################################
# 4. RUN SIMULATIONS FOR CENSORED VS COUNT MODELS
################################################################################
# Wrapper function to run multiple simulations, fit GAMs for:
#  - Binary outcome (censored model)
#  - Count outcome (Poisson model)
# Returns estimated coefficients and fitted models for each simulation.

run_simulation_cens <- function(lambda_0_val, label,
                                exogenous_covariates = c("girl_alter","age"),
                                endogenous_covariates = NULL,
                                nonlinear_covariates = c("age"),
                                beta = beta_true,
                                random_intercept_cov = NULL) {
  
  # Create baseline hazard vector replicated across years
  lambda_0 <- rep(lambda_0_val, n_years)
  
  # Initialize lists to store results for each simulation
  betas_cens_list <- list()       # coefficients from censored model
  betas_unc_list  <- list()       # coefficients from count (Poisson) model
  gam_fit_cens_list  <- list()    # full GAM models for censored
  gam_fit_count_list <- list()    # full GAM models for count
  
  # Loop over number of simulations
  for (i in 1:num_simulations) {
    cat("Simulation", i, "-", label, "\n")
    
    # Simulate events
    events <- simulate_events(Risk_Set, lambda_0, beta, n_years,
                              endogenous_cov = endogenous_covariates,
                              exogenous_cov = exogenous_covariates, 
                              nonlinear_cov = nonlinear_covariates,
                              random_intercept_cov = NULL,
                              gamma = NULL,
                              gamma_slope = NULL,
                              random_slope_cov = NULL)
    
    # Convert event history into list of events
    event_list <- lapply(1:length(events$history$time), function(i) {
      list(time = events$history$time[i], 
           receiver = events$history$receiver[i],  
           sender = events$history$sender[[i]])
    })
    
    # Binary indicator for censored model
    z_cens <- define_vector_z(RiskSet = Risk_Set, EventList = unique(event_list))
    # Count variable for Poisson model
    z_count <- define_vector_z_count(RiskSet = Risk_Set, EventList = event_list)
    
    # Prepare covariates data frame
    covariates_df <- data.frame(
      z_cens = z_cens,
      z_count = z_count,
      girl_alter = Risk_Set$girl_alter,
      age = Risk_Set$age
    )
    
    # Add endogenous covariates from simulated events
    for (cov in endogenous_covariates) {
      covariates_df[[cov]] <- events$covariates[[cov]]
    }
    
    # Random intercepts
    get_group_component <- function(covariate) {
      switch(covariate,
             "ClassID" = as.character(floor(Risk_Set$Receiver / 100) * 100),
             "Receiver" = as.character(Risk_Set$Receiver),
             "Senders" = sapply(Risk_Set$Senders, function(s) paste(sort(as.character(s)), collapse = "_")),
             as.character(Risk_Set[[covariate]]))
    }
    
    # Apply random intercept covariates
    if (length(random_intercept_cov) == 1) {
      covariates_df[[random_intercept_cov]] <- factor(get_group_component(random_intercept_cov))
    } else {
      for (covariate in random_intercept_cov) {
        covariates_df[[covariate]] <- factor(get_group_component(covariate))
      }
    }
    
    # Smooth terms for GAM: function to decide nb of knots based on unique values
    make_smooth_term <- function(cov_name, data) {
      n_unique <- length(unique(data[[cov_name]]))
      if (n_unique < 10) {
        paste0("s(", cov_name, ", k = ", n_unique, ")")
      } else {
        paste0("s(", cov_name, ")")
      }
    }
    smooth_terms <- sapply(nonlinear_covariates, make_smooth_term, data = covariates_df)
    linear_covariates <- setdiff(c(exogenous_covariates, endogenous_covariates), nonlinear_covariates)
    
    # Random effect term for GAM if applicable
    if (!is.null(random_intercept_cov) && nzchar(trimws(random_intercept_cov))) {
      random_effect_terms <- paste0("s(", random_intercept_cov, ", bs = 're')")
    } else {
      random_effect_terms <- NULL
    }
    rhs_parts <- c(linear_covariates, smooth_terms, random_effect_terms)
    rhs_parts <- rhs_parts[nzchar(rhs_parts)]
    if (length(rhs_parts) == 0) stop("Formula RHS is empty")
    
    rhs_combined <- paste(rhs_parts, collapse = " + ")
    
    # Fit GAM models
    
    # Censored model
    formula_obj_cens <- as.formula(paste("z_cens ~", rhs_combined))
    gam_fit_cens <- gam(formula = formula_obj_cens, family = binomial(link = "cloglog"),
                        data = covariates_df, method = "REML", select = TRUE)
    
    # Count (Poisson) model
    formula_obj_unc <- as.formula(paste("z_count ~", rhs_combined))
    gam_fit_count <- gam(formula = formula_obj_unc, family = poisson(link = "log"),
                         data = covariates_df, method = "REML", select = TRUE)
    
    # Store results for this simulation
    betas_cens_list[[i]] <- coef(gam_fit_cens)[linear_covariates]
    betas_unc_list[[i]]  <- coef(gam_fit_count)[linear_covariates]
    gam_fit_cens_list[[i]]  <- gam_fit_cens
    gam_fit_count_list[[i]] <- gam_fit_count
  }
  
  return(list(
    betas_cens_list = betas_cens_list,
    betas_unc_list = betas_unc_list,
    gam_fit_cens_list = gam_fit_cens_list,
    gam_fit_count_list = gam_fit_count_list,
    covariates_df = covariates_df
  ))
}








################################################################################
# 5. SIMULATE HYPEREVENTS USING GLOBAL COVARIATE
################################################################################
# Simulates continuous-time events where hazard depends on a global, time-varying covariate.

simulate_events_glob <- function(RiskSet, lambda_0, n_years, global_cov, global_coef, dt = 0.1) {
  # Initialize output
  history <- list(time = numeric(), receiver = numeric(), sender = list())
  covariates_list <- list()
  hazard_true_list <- list()
  
  tl <- 0             # Left endpoint of current time window
  event_count <- 0    # Counter for total number of simulated events
  
  cat("Starting simulation...\n")
  
  # Main loop: iterate over continuous time in small increments of dt
  while (tl < n_years) {
    tu <- tl + dt           # Right endpoint of the current interval
    mid_t <- tl + dt / 2    # Midpoint (for evaluating covariates)
    
    # Compute hazard for each sender–receiver pair
    #   λ_sr(t) = λ₀(t) * exp(β * X(t))
    l_sr <- rep(lambda_0(mid_t) * exp(global_coef * global_cov(mid_t)), nrow(RiskSet))
    total_intensity <- sum(l_sr)
    
    # Skip interval if no risk
    if (total_intensity <= 0) {
      tl <- tu
      next
    }
    
    # Draw inter-arrival time for the next possible event
    tm <- rexp(1, total_intensity)
    
    dt_r <- dt    # Remaining width in this sub-interval
    tl_r <- tl    # Current left edge
    
    # Inner loop: handle multiple potential events within one dt window
    while (tm < dt_r) {
      event_count <- event_count + 1
      
      # Draw one event proportional to hazard rates
      link <- rmultinom(1, 1, prob = l_sr / total_intensity)
      id <- which(link == 1)
      
      st <- tl_r + tm   # Actual event time
      
      # Stop if event happens beyond simulation horizon
      if (st >= n_years) break
      
      # Record event details
      chosen <- RiskSet[id, ]
      history$time <- c(history$time, st)
      history$receiver <- c(history$receiver, chosen$Receiver)
      history$sender <- append(history$sender, list(unlist(chosen$Senders)))
      
      # Update for next potential event
      tl_r <- st                   # Move start to new event time
      dt_r <- tu - tl_r            # Remaining time in the current interval 
      mid_t_r <- tl_r + dt_r / 2   # New midpoint
      
      # Recalculate hazards at updated time
      l_sr <- rep(lambda_0(mid_t_r) * exp(global_coef * global_cov(mid_t_r)), nrow(RiskSet))
      total_intensity <- sum(l_sr)
      
      if (total_intensity <= 0) break
      
      # Draw next inter-arrival time
      tm <- rexp(1, total_intensity)
    }
    
    # Move forward to next dt window
    tl <- tu
  }
  
  cat("Simulation completed!\n")
  
  return(list(
    history = history,
    covariates = covariates_list,
    hazard_true = hazard_true_list
  ))
}









################################################################################
# 6. AGGREGATE RISKSET BY WAVES
################################################################################
# Aggregates sender-receiver pairs by discrete wave for covariate evaluation.

aggregate_cov <- function(Risk_Set, n_years) {
  Risk_Set$Senders <- lapply(Risk_Set$Senders, function(x) as.integer(unlist(x)))
  agg_list <- list()
  
  for (w in seq_len(n_years)) {
    links_wave <- Risk_Set[Risk_Set$Time == w, ]
    agg_wave <- lapply(seq_len(nrow(links_wave)), function(i) {
      rec <- links_wave$Receiver[i]
      senders <- links_wave$Senders[[i]]
      data.frame(
        wave = w,
        Receiver = rec,
        Senders = paste(senders, collapse = ",")
      )
    })
    agg_list[[w]] <- do.call(rbind, agg_wave)
  }
  agg_df <- do.call(rbind, agg_list)
  rownames(agg_df) <- NULL
  return(agg_df)
}











################################################################################
# 7. RUN SINGLE SIMULATION WITH GLOBAL COVARIATE
################################################################################
# Wrapper function for a single simulation using a global covariate.

run_simulation <- function(lambda_0_val, n_years, RiskSet, global_cov, global_coef) {
  simulation_results <- list()
  lambda_0 <- lambda_0_val      # Constant baseline hazard
  
  # Run simulation
  events <- simulate_events_glob(Risk_Set, lambda_0, n_years, global_cov, global_coef)
  
  # Convert to event list format
  event_list <- lapply(1:length(events$history$time), function(j) {
    list(time = events$history$time[j],
         receiver = events$history$receiver[j],
         sender = events$history$sender[[j]])
  })
  simulation_results[[1]] <- event_list
  return(simulation_results)
}








################################################################################
# 8. DEFINE BINARY VECTOR Z
################################################################################
# Creates a 0/1 indicator vector marking which events in the RiskSet actually occurred.

define_vector_z <- function(RiskSet, EventList = NULL) {
  # Initialize 
  z_sim <- numeric(nrow(RiskSet))
  if (!is.null(EventList)) {
    # Loop over all possible risk events
    for (i in 1:nrow(RiskSet)) {
      # Extract information from the i-th row of the RiskSet
      risk_time <- RiskSet$Time[i]
      risk_receiver <- RiskSet$Receiver[i]
      risk_senders <- RiskSet$Senders[i][[1]]
      # Compare current risk event with all observed events in EventList
      z_sim[i] <- any(sapply(EventList, function(hist_event) {
        # Extract observed event components
        hist_time <- hist_event$time
        hist_receiver <- hist_event$receiver
        hist_senders <- hist_event$sender
        # Check equality of time, receiver, and senders 
        if (length(risk_senders) == length(hist_senders)) {
          return(risk_time == hist_time && risk_receiver == hist_receiver && all(risk_senders == hist_senders))
        } else {
          return(FALSE)  
        }
      }))
    }
  }
  
  # Return binary indicator vector
  if (!is.null(EventList)) {
    return(z_sim)
  } else {
    stop("EventList must be provided.")
  }
}









################################################################################
# 9. DEFINE COUNT VECTOR Z
################################################################################
# Similar to define_vector_z(), but counts number of times each event occurred.

define_vector_z_count <- function(RiskSet, EventList = NULL) {
  # Initialize
  z_sim <- numeric(nrow(RiskSet))
  if (!is.null(EventList)) {
    # Loop through each potential event in the RiskSet
    for (i in 1:nrow(RiskSet)) {
      # Extract candidate event
      risk_time <- RiskSet$Time[i]
      risk_receiver <- RiskSet$Receiver[i]
      risk_senders <- RiskSet$Senders[[i]]
      # Count how many times this event appears in EventList
      z_sim[i] <- sum(sapply(EventList, function(hist_event) {
        hist_time <- hist_event$time
        hist_receiver <- hist_event$receiver
        hist_senders <- hist_event$sender
        # Checks
        same_time <- (risk_time == hist_time)
        same_receiver <- (risk_receiver == hist_receiver)
        same_senders <- (length(risk_senders) == length(hist_senders) &&
                           all(risk_senders == hist_senders))
        return(same_time && same_receiver && same_senders)
      }))
    }
    return(z_sim)
    
  } else {
    stop("EventList must be provided.")
  }
}
