################################################################################
########################## COVARIATE DEFINITIONS ###############################
################################################################################
#
# Description:
#   Defines all covariates (exogenous and endogenous) used in the gossip network
#   hyperevent model. This file includes:
#     1. Utility functions
#     2. Exogenous covariates (gender, age)
#     3. Endogenous covariates (degree, repetition, reciprocity, balance, etc.)
#     4. Generic function to compute covariates
#
# Notes:
#   - Exogenous covariates depend on static attributes of nodes.
#   - Endogenous covariates depend on the current state of the network.
################################################################################



###############################################################################
# 1. UTILITY FUNCTIONS
###############################################################################

# Helper function to retrieve the appropriate gender dataset name
get_gender_data_name <- function(id) {
  id_str <- as.character(id)
  class_group <- paste0(substr(id_str, 1, 2), "00")
  return(paste0("gender_", class_group))
}













###############################################################################
# 2. EXOGENOUS COVARIATES
###############################################################################


### --------------------------------------------------------------------------
### GENDER
### --------------------------------------------------------------------------

# Gender proportion of sender(s): 0 to 1 depending on the fraction of female senders
gender_sender <- function(RiskSet, gender_map) {
  # Initialize a numeric vector to store the proportion for each row
  gender_sender_vec <- numeric(nrow(RiskSet)) 
  
  # Loop over each row in the RiskSet
  for (i in 1:nrow(RiskSet)) {
    # Get the senders for this row
    senders <- RiskSet$Senders[[i]]
    
    # Map sender IDs to their gender (1 = female, 0 = male)
    gender_senders <- gender_map[as.character(senders)]
    
    # Count how many female senders there are
    female_count <- sum(gender_senders == 1, na.rm = TRUE)
    
    # Total number of senders for this row
    total_senders <- length(senders)
    
    # Calculate the proportion of female senders
    # If there are no senders, return NA
    gender_sender_vec[i] <- ifelse(total_senders > 0, female_count / total_senders, NA)
  }
  
  # Return the vector of proportions
  return(gender_sender_vec)
}





# Gender of receiver: 1 = female, 0 = male
gender_receiver <- function(RiskSet, gender_map) {
  gender_receiver_vec <- numeric(nrow(RiskSet))  
  
  for (i in 1:nrow(RiskSet)) {
    receiver <- as.character(RiskSet$Receiver[i])
    gender_val <- gender_map[receiver]
    gender_receiver_vec[i] <- ifelse(gender_val == 1, 1, 0)  
  }
  
  return(gender_receiver_vec)
}







### --------------------------------------------------------------------------
### AGE
### --------------------------------------------------------------------------

# Computes a simple age covariate: base 13 years + wave number
age <- function(RiskSet) {
  age_covariate <- numeric(nrow(RiskSet))
  for (i in 1:nrow(RiskSet)) {
    wave <- RiskSet$Time[i]
    age_covariate[i] <- 13 + wave  # Wave 1 -> 14 years, Wave 2 -> 15 years , etc.
  }
  return(age_covariate)
}










###############################################################################
# 3. ENDOGENOUS COVARIATES
###############################################################################

### --------------------------------------------------------------------------
### DEGREE
### --------------------------------------------------------------------------

sender_degree <- function(history, event) {
  
  # Extract current event details
  time <- event$time
  S <- event$sender 
  if (is.null(history) || time <= 1) {
    return(0)
  }
  # Step 1: Identify all past events before the current event
  past_events <- which(history$time < time)
  if (length(past_events) == 0) {
    return(0)
  }
  # Initialize sender degree sum
  sender_degree_sum <- 0
  # Step 2: Loop through past events and count occurrences where senders match
  for (i in past_events) {
    S_i <- history$sender[[i]]  # Extract the sender list for past event i
    # Check if the sender set S is a subset of the past sender's set S_i
    if (all(S %in% S_i)) {
        sender_degree_sum <- sender_degree_sum + 1
    }
  }
  # Return sender degree (normalized by number of senders)
  return(sender_degree_sum / length(S))
}




receiver_degree <- function(history, event) {
  
  # Extract current event details
  time <- event$time
  r <- event$receiver 
  if (is.null(history) || time <= 1) {
    return(0)
  }
  # Step 1: Identify past events before the current event where the same receiver was involved
  past_events <- which((history$receiver == r) & (history$time < time))
  if (length(past_events) == 0) {
    return(0)
  }
  # Step 2: Count occurences
  receiver_degree_sum <- 0
  for (i in past_events) {
      receiver_degree_sum <- receiver_degree_sum + 1
  }
  return(receiver_degree_sum)
}






### --------------------------------------------------------------------------
### REPETITION
### --------------------------------------------------------------------------


repetition <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender
  r <- event$receiver
  if (is.null(history)) {return(0)}
  if (time <= 1) {return(0)}
  # Step 1: Identify past events where the same receiver interacted before time t
  past_events <- which(
    (history$receiver == r) &  # Same receiver
      (history$time < time)    # Before current event
  )
  # If there are no past events, return 0
  if (length(past_events) == 0) {return(0)}
  counted_nodes <- c()
  repetition_sum <- 0
  # Step 2: Count occurrences where senders match (order-independent)
  for (i in past_events) {
    S_i <- history$sender[[i]]  # Extract past sender list
    if (setequal(S, S_i)) {
      # Check whether the same event has already been recorded in the same wave
      A <- which(
        (history$sender %in% list(history$sender[[i]])) &
          (history$receiver == history$receiver[i]) &
          (history$time == history$time[i])
      )
      if (A[1] == i) {     # Ensure unique counting within the same wave
        repetition_sum <- repetition_sum + 1
      }
    }
  }
  # Normalize by the number of senders |S|
  return(repetition_sum / length(S))
}








library(combinat)  # for combn

subset_repetition <- function(history, event) {
  time <- event$time
  S <- event$sender
  r <- event$receiver
  
  if (is.null(history) || time <= 1) {return(0)}
  # Filter relevant past events: same receiver, before current time
  past_events <- which((history$receiver == r) & (history$time < time))
  if (length(past_events) == 0) {return(0)}
  S_size <- length(S)
  total_score <- 0
  # Loop through subset sizes: p = 1 to |S|
  for (p in 1:S_size) {
    subsets_p <- combn(S, p, simplify = FALSE)
    num_subsets <- length(subsets_p)
    match_count <- 0
    for (subset in subsets_p) {
      for (i in past_events) {
        S_i <- history$sender[[i]]
        if (all(subset %in% S_i)) {
          match_count <- match_count + 1
          break  # Count subset only once if matched in any past event
        }
      }
    }
    # Normalize by total number of subsets of size p
    if (num_subsets > 0) {
      total_score <- total_score + (match_count / num_subsets)
    }
  }
  return(total_score)
}








### --------------------------------------------------------------------------
### RECIPROCITY/RETALIATION
### --------------------------------------------------------------------------

reciprocity <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender
  r <- event$receiver
  # Initial check
  if (is.null(history)) {return(0)}
  if (time <= 1) {return(0)}
  # Step 1: Identify past events where the current receiver was a sender
  past_events <- which(
    sapply(history$sender, function(s) r %in% s) &  # Check if r is in sender list
      (history$time < time)                         # Before current event
  )
  # If there are no past events, return 0
  if (length(past_events) == 0) {return(0)}
  recip_sum <- 0
  # Step 2: Count occurrences where senders match (order-independent)
  for (i in past_events) {
    past_receiver <- history$receiver[i]  # Extract past receiver
    if (past_receiver %in% S) {           # Check if past receiver is in current senders
      recip_sum <- recip_sum + 1
    }
  }
  return(recip_sum / length(S))  # Normalize by the number of senders
}






### --------------------------------------------------------------------------
### TRIADIC CLOSURE
### --------------------------------------------------------------------------

transitive_closure <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender
  r <- event$receiver
  # If first event, return 0 (no past data)
  if (is.null(history) || time <= 1) {return(0)}
  # Initialize transitive closure sum
  tc_sum <- 0
  # Step 1: Find indexes in history of past interactions where S gossiped about some node a
  past_interactions_a <- which(
    !(history$receiver %in% c(S, r)) &   # Receiver should NOT be in S or r
      (history$sender %in% list(S)) &    # Senders should match exactly 
      (history$time < time)              # Make sure it happens before the current wave
  )
  counted_nodes <- c()
  # Step 2: Compute transitive closure sum over valid third-party nodes
  for (a_index in past_interactions_a) {
    # Extract third-party node a
    a <- history$receiver[a_index]
    # Ensure a ≠ S, r
    a <- a[!a %in% c(S, r)]
    # Check whether we already counted this sender[a_index] and receiver[a_index] in the current wave, if so, then skip
    A <- which(
      ((history$sender %in% history$sender[a_index]) &
         (history$receiver == history$receiver[a_index]) &
         (history$time == history$time[a_index])
      )
    )
    if (A[1] == a_index) {    # Ensure unique counting within the same wave
      # Compute past interactions between S and a
      past_interactions_S_a <- ((history$sender %in% list(S)) &     # Check if sender list matches
                                  (history$receiver == a) &         # Receiver should match with a
                                  (history$time < time))            # Ensure interaction happened before current time 
      # Compute hy_deg_t(S, a): Number of waves S gossiped about a
      hy_deg_S_a <- length(unique(history$time[past_interactions_S_a]))
      # Compute past interactions between a and r
      past_interactions_a_r <- (sapply(history$sender, function(x) a %in% x)) &  # Sender should contain a
        (history$receiver == r) &                                                # Receiver must be exactly r
        (history$time < time)                                                    # Interaction must happen before current time
      # Compute hy_deg_t(a, r): Number of waves `a` gossiped about r
      hy_deg_a_r <- length(unique(history$time[past_interactions_a_r]))
      # Ensure each node is counted only once when iterating over indexes, 
      # avoiding double counting for nodes appearing in multiple waves.
      for (i in which(past_interactions_S_a)) {  # Iterate only over TRUE indices
        node_id <- history$receiver[i]           # Get the actual node ID
        if (!(node_id %in% counted_nodes)) {     # Prevent duplicate counting
          tc <- min(hy_deg_S_a, hy_deg_a_r)
          tc_sum <- tc_sum + tc
          counted_nodes <- c(counted_nodes, node_id)  
        }
      }
    }
  }
  return(tc_sum / length(S))
}







cyclic_closure <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender
  r <- event$receiver
  # If first event, return 0 (no past data)
  if (is.null(history) || time <= 1) {return(0)}
  # Initialize cyclic closure sum
  cc_sum <- 0
  # Step 1: Find indexes in history of past interactions where some node a gossiped about a certain s in S
  past_interactions_a <- which(
    (history$receiver %in% S) &         # Ensure receiver stays in S
      !(history$sender %in% c(S, r)) &  # Exclude sender if in S or r
      (history$time < time)             # Ensure past interaction
  )
  counted_nodes <- c()
  # Step 2: Compute cyclic closure sum over valid third-party nodes
  for (a_index in past_interactions_a) {
    # Extract third-party node a
    vector_a <- history$sender[a_index]
    for (inner_list in vector_a) {
      for (a in inner_list) {
        # Ensure a ≠ S, r
        if (a %in% c(S, r)) next
        # Check whether we already counted this sender[a_index] and receiver[a_index] in the current wave, if so, then skip
        A <- which(
          ((history$sender %in% history$sender[a_index]) &
             (history$receiver == history$receiver[a_index]) &
             (history$time == history$time[a_index]))
        )
        # Ensure unique counting within the same wave
        if (A[1] == a_index) {
          # Compute hy_deg_t(a, S): Number of waves a gossiped about S
          past_interactions_a_S <- (sapply(history$sender, function(x) a %in% x)) &   # Sender should contain a
            (history$receiver %in% S) &                                               # Receiver should be in S
            (history$time < time)                                                     # Ensure past interaction
          hy_deg_a_S <- length(unique(history$time[past_interactions_a_S]))
          # Compute hy_deg_t(r, a): Number of waves r gossiped about a
          past_interactions_r_a <- (sapply(history$sender, function(x) r %in% x)) &   # Sender should contain r
            (history$receiver == a) &                                                 # Receiver must be exactly a
            (history$time < time)                                                     # Interaction must happen before current time
          hy_deg_r_a <- length(unique(history$time[past_interactions_r_a]))
          for (i in which(past_interactions_a_S)) {  # Iterate only over TRUE indices
            node_id <- a                             # Get the actual node ID
            if (!(node_id %in% counted_nodes)) {     # Prevent duplicate counting
              cc <- min(hy_deg_a_S, hy_deg_r_a)
              cc_sum <- cc_sum + cc
              counted_nodes <- c(counted_nodes, node_id) 
            }
          }
        }
      }
    }
  }
  return(cc_sum / length(S))
}










### --------------------------------------------------------------------------
### BALANCE
### --------------------------------------------------------------------------


sender_balance <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender 
  r <- event$receiver 
  # Initial check: If history is empty or it's the first event, return 0
  if (is.null(history) || time <= 1) {return(0)}
  # Initialize sender balance sum
  sb_sum <- 0
  # Step 1: Find past interactions where some node a sent gossip to some node in S
  past_interactions_a <- which(
    (history$receiver %in% S) &         # Receiver should stay in S
      !(history$sender %in% c(S, r)) &  # Exclude sender if in S or r
      (history$time < time)             # Ensure past interactions
  )
  counted_nodes <- c()
  # Step 2: Compute balance sum over valid third-party nodes
  for (a_index in past_interactions_a) {
    # Extract third-party node a
    vector_a <- history$sender[a_index]
    for (inner_list in vector_a) {
      for (a in inner_list) {
        # Ensure a ≠ S, r
        if (a %in% c(S, r)) next
        # Check whether we already counted this sender[a_index] and receiver[a_index] in the current wave, if so, then skip
        A <- which(
          ((history$sender %in% history$sender[a_index]) &
             (history$receiver == history$receiver[a_index]) &
             (history$time == history$time[a_index]))
        )
        if (A[1] == a_index){ 
          past_interactions_a_S <- (sapply(history$sender, function(x) a %in% x)) &  # Sender should contain a
            (history$receiver %in% S) &                                              # Receiver should be in S
            (history$time < time)                                                    # Ensure past interaction
          # Compute hy_deg_t(a, S): Number of waves a gossiped about S
          hy_deg_a_S <- length(unique(history$time[past_interactions_a_S]))
          past_interactions_a_r <- (sapply(history$sender, function(x) a %in% x)) &   # Sender should contain a
            (history$receiver == r) &                                                 # Receiver must be exactly r
            (history$time < time)                                                     # Interaction must happen before current time
          # Compute hy_deg_t(a, r): Number of waves a gossiped about r
          hy_deg_a_r <- length(unique(history$time[past_interactions_a_r]))
          for (i in which(past_interactions_a_S)) {  # Iterate only over TRUE indices
            node_id <- a                             # Get the actual node ID
            if (!(node_id %in% counted_nodes)) {     # Prevent duplicate counting
              sb <- min(hy_deg_a_S, hy_deg_a_r)
              sb_sum <- sb_sum + sb
              counted_nodes <- c(counted_nodes, node_id)  
            }
          }
        }
      }
    }
  }
  return(sb_sum / length(S))
}








receiver_balance <- function(history, event) {
  
  # Extract the current event details
  time <- event$time
  S <- event$sender
  r <- event$receiver
  # Initial check
  if (is.null(history)) {return(0)}
  if (time<=1){return(0)}
  # Initialize receiver balance sum
  balance_sum <- 0
  # Step 1: Find indexes in history of past interactions where S gossiped about some node a
  past_interactions_a <- which(
    !(history$receiver %in% c(S, r)) &   # Receiver should NOT be in S or r
      (history$sender %in% list(S)) &    # Senders should match exactly 
      (history$time < time))             # Make sure it happens before the current wave
  counted_nodes <- c()
  # Step 2: Compute balance sum over valid third-party nodes
  for (a_index in past_interactions_a) {
    # Extract all potential third-party nodes a
    a <- history$receiver[a_index]
    # Ensure a ≠ S, r
    a <- a[!a %in% c(S, r)]
    # Check whether we already counted this sender[a_index] and receiver[a_index] in the current wave, if so, then skip
    A <- which(((history$sender %in% history$sender[a_index]) &
                  (history$receiver == history$receiver[a_index]) &
                  #(floor(history$time) == floor(history$time[a_index])) 
                  (history$time == history$time[a_index])
    ))
    if (A[1] == a_index){  
      past_interactions_S_a <-((history$sender %in% list(S)) &     # Check if sender list matches
                               (history$receiver == a) &           # Receiver should match with a
                               (history$time < time))              # Ensure interaction happened before current time
      # Compute hy_deg_t(S, a): Number of waves S gossiped about a
      hy_deg_S_a <- length(unique(history$time[past_interactions_S_a])) 
      past_interactions_r_a <- (sapply(history$sender, function(x) r %in% x)) &  # Ensure r is inside the list of sender
        (history$receiver == a) &                                                # Ensure a is the receiver
        (history$time < time)                                                    # Ensure interaction happened before current wave
      # Compute hy_deg_t(r, a): Number of waves r gossiped about a 
      hy_deg_r_a <- length(unique(history$time[past_interactions_r_a])) 
      # Ensure each node is counted only once when iterating over indexes, 
      # avoiding double counting for nodes appearing in multiple waves.
      for (i in which(past_interactions_S_a)) {     # Iterate only over TRUE indices
        node_id <- history$receiver[i]              # Get the actual node ID
        if (!(node_id %in% counted_nodes)) {        # Prevent duplicate counting
          balance <- min(hy_deg_S_a, hy_deg_r_a)
          balance_sum <- balance_sum + balance
          counted_nodes <- c(counted_nodes, node_id)  
        }
      }
    }
  }
  return(balance_sum/length(S))
}











###############################################################################
# 4. GENERIC FUNCTION TO COMPUTE COVARIATES
###############################################################################

# -------------------------------------------------------------------------                                       
# Generic function to compute any covariate for the gossip hyperevent model.
# Handles both:
#   - Exogenous covariates (girl.ego, girl.alter, age, etc.)
#   - Endogenous covariates (degree, repetition, reciprocity, balance, etc.)
# -------------------------------------------------------------------------

compute_covariates <- function(RiskSet, CovariateName, History = NULL, Event = NULL, gender_map = NULL) {
  
  # Check if CovariateName is endogenous or exogenous
  if (CovariateName %in% c("sender.degree", "receiver.degree", "repetition", 
                           "subset.repetition", "reciprocity", "transitive.closure", 
                           "cyclic.closure", "sender.balance", "receiver.balance")) {
    # Endogenous covariates require History as input
    if (is.null(History)) {
      stop("For endogenous covariates, History must be provided.")
    }
    # If CovariateName is one of the endogenous covariates, process the event
    if (is.null(Event)) {
      stop("For endogenous covariates, Event must be provided.")
    }
    # Compute the endogenous covariate for a single event
    if (CovariateName == "sender.degree") {
      CovariateValues <- sender_degree(History, Event)
    } else if (CovariateName == "receiver.degree") {
      CovariateValues <- receiver_degree(History, Event)
    } else if (CovariateName == "repetition") {
      CovariateValues <- repetition(History, Event)
    } else if (CovariateName == "subset.repetition") {
      CovariateValues <- subset_repetition(History, Event)
    } else if (CovariateName == "reciprocity") {
      CovariateValues <- reciprocity(History, Event)
    } else if (CovariateName == "transitive.closure") {
      CovariateValues <- transitive_closure(History, Event)
    } else if (CovariateName == "cyclic.closure") {
      CovariateValues <- cyclic_closure(History, Event)
    } else if (CovariateName == "sender.balance") {
      CovariateValues <- sender_balance(History, Event)
    } else if (CovariateName == "receiver.balance") {
      CovariateValues <- receiver_balance(History, Event)
    } else {
      stop("Unknown endogenous covariate: ", CovariateName)
    }
  } else {
    # Handle exogenous covariates (whole RiskSet is used)
    if (CovariateName == "girl.ego") {
      CovariateValues <- gender_sender(RiskSet, gender_map)
    } else if (CovariateName == "girl.alter") {
      CovariateValues <- gender_receiver(RiskSet, gender_map)
    } else if (CovariateName == "age") {
      CovariateValues <- age(RiskSet)
    } else if (!(CovariateName %in% colnames(RiskSet))) {
      stop("Covariate not found in the dataset: ", CovariateName)
    } else {
      # If it's a simple column, just extract the values from RiskSet
      CovariateValues <- RiskSet[[CovariateName]]
    }
  }
  # Return the computed covariates
  return(CovariateValues)
}
