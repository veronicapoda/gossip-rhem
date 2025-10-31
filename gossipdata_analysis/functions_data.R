################################################################################
########################### DATA ANALYSIS FUNCTIONS ############################
################################################################################
# These functions are used to generate the risk set and define observed events
# for real gossip data.
################################################################################



################################################################################
# CONSTRUCT THE RISK SET
################################################################################
# Construct all possible sender–receiver combinations for a given wave
# considering friendship relations.

generate_risk_set <- function(max_senders, Wave) {
  clean_riskset <- list()
  
  # Helper: create a unique key to avoid duplicates
  make_key <- function(senders, receiver, wave) {
    senders <- as.character(senders)
    paste(sort(senders), as.character(receiver), as.character(wave), collapse = "_")
  }
  
  # Iterate over all filtered friendship pairs
  for (pair in filtered_friendships) {
    base_senders <- as.numeric(pair)
    
    # Skip invalid IDs
    if (!all(sapply(base_senders, is_valid_id))) next
    
    # Step 1: 2-person sender groups
    S2 <- base_senders
    possible_receivers <- setdiff(as.numeric(all_students_class), S2)
    possible_receivers <- possible_receivers[sapply(possible_receivers, is_valid_id)]
    
    # Create risk set entries for 2-person sender groups
    for (d in possible_receivers) {
      key <- make_key(S2, d, Wave)
      if (!exists(key, envir = existing_keys_global)) {
        clean_riskset <- append(clean_riskset, list(list(senders = S2, receiver = d)))
        assign(key, TRUE, envir = existing_keys_global)
      }
    }
    
    # Step 2: Expand sender groups up to max_senders
    current_senders <- list(S2)
    for (size in 3:max_senders) {
      next_senders_list <- list()
      for (S in current_senders) {
        potential_extra <- setdiff(as.numeric(students), S)
        potential_extra <- potential_extra[sapply(potential_extra, is_valid_id)]
        
        for (new_sender in potential_extra) {
          S_new <- as.numeric(c(S, new_sender))
          
          # Skip invalid combinations
          if (!all(sapply(S_new, is_valid_id))) next
          
          # Check if new sender is friend with all in current group
          is_friend_with_all <- all(sapply(S, function(sender) {
            any(apply(filtered_friendships_matrix, 1, function(x) all(sort(c(sender, new_sender)) == sort(as.numeric(x)))))
          }))
          
          if (is_friend_with_all) {
            possible_receivers <- setdiff(as.numeric(all_students_class), S_new)
            possible_receivers <- possible_receivers[sapply(possible_receivers, is_valid_id)]
            
            for (d in possible_receivers) {
              key <- make_key(S_new, d, Wave)
              if (!exists(key, envir = existing_keys_global)) {
                clean_riskset <- append(clean_riskset, list(list(senders = S_new, receiver = d)))
                assign(key, TRUE, envir = existing_keys_global)
              }
            }
            
            next_senders_list <- append(next_senders_list, list(S_new))
          }
        }
      }
      current_senders <- next_senders_list
      if (length(current_senders) == 0) break
    }
  }
  
  # Convert list to dataframe
  if (length(clean_riskset) == 0) {
    df <- data.frame(Senders = I(list()), Receiver = numeric(0), Time = numeric(0))
  } else {
    df <- do.call(rbind, lapply(clean_riskset, function(x) {
      data.frame(
        Senders = I(list(x$senders)),
        Receiver = x$receiver,
        Time = as.numeric(Wave),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  cat("Wave", Wave, "- Total risk set rows:", nrow(df), "\n")
  return(df)
}









################################################################################
# DEFINE VECTOR Z OF EVENTS
################################################################################
# Create a binary vector indicating which events in the risk set actually
# occurred (gossip observed).

define_vector_z <- function(RiskSet, GossipData = NULL, EventList = NULL) {
  
  z_dat <- numeric(nrow(RiskSet))   # observed gossip
  z_sim <- numeric(nrow(RiskSet))   # historical/simulated events
  
  # Compare with GossipData
  if (!is.null(GossipData)) {
    gossip_norm <- data.frame(
      senders = sapply(strsplit(GossipData$data$Gossipers, ",\\s*"), function(x) paste(sort(trimws(x)), collapse = ",")),
      target  = as.character(trimws(GossipData$data$Target)),
      time    = as.integer(gsub("Wave", "", GossipData$data$Time)),
      stringsAsFactors = FALSE
    )
    risk_norm <- data.frame(
      senders = sapply(RiskSet$Senders, function(x) paste(sort(trimws(as.character(x))), collapse = ",")),
      target  = as.character(trimws(RiskSet$Receiver)),
      time    = as.integer(RiskSet$Time),
      stringsAsFactors = FALSE
    )
    
    # Compare unique keys
    gossip_keys <- paste(gossip_norm$senders, gossip_norm$target, gossip_norm$time, sep = "_")
    risk_keys   <- paste(risk_norm$senders,  risk_norm$target,  risk_norm$time, sep = "_")
    
    z_dat <- as.numeric(risk_keys %in% gossip_keys)
  }
  
  # Compare with EventList
  if (!is.null(EventList)) {
    event_norm <- lapply(EventList, function(e) list(
      senders = paste(sort(as.character(e$sender)), collapse = ","),
      receiver = as.character(e$receiver),
      time = as.integer(e$time)
    ))
    for (i in seq_len(nrow(RiskSet))) {
      risk_key <- paste(
        paste(sort(as.character(RiskSet$Senders[[i]])), collapse = ","),
        as.character(RiskSet$Receiver[i]),
        as.integer(RiskSet$Time[i]),
        sep = "_"
      )
      
      z_sim[i] <- any(sapply(event_norm, function(ev) {
        paste(ev$senders, ev$receiver, ev$time, sep = "_") == risk_key
      }))
    }
  }
  
  # Return appropriate z vector
  if (!is.null(GossipData)) {
    return(z_dat)
  } else if (!is.null(EventList)) {
    return(z_sim)
  } else {
    stop("Either GossipData or EventList must be provided.")
  }
}









################################################################################
# IDENTIFY THE MAXIMUM FRIENDSHIP CLIQUE ACROSS WAVES AND CLASSES
################################################################################
# Scans all friendship matrices across waves and classes, identifies reciprocated
# friendships, and computes maximal cliques. Returns both the largest clique size
# and all clique sizes.
################################################################################

get_max_clique_size <- function(base_path = NULL) {
  library(igraph)
  
  if (is.null(base_path)) {
    base_path <- file.path(getwd(), 
                           "gossipdata_analysis",
                           "gossipdata",
                           "RECENS_Wired_into_each_other-Data_files",
                           "RECENS_networks_dataset_w1234")
  }
  
  all_clique_sizes <- c()
  waves <- c("Wave1", "Wave2", "Wave3", "Wave4")
  folders <- c("1000", "2000", "3000", "4000", "5000", "6000", "7000")
  
  for (wave in waves) {
    for (folder in folders) {
      folder_path <- paste0(base_path, "/", wave, "/", folder)
      classes_in_folder <- list.files(path = folder_path, full.names = FALSE)
      
      for (class in classes_in_folder) {
        class_path <- file.path(folder_path, class)
        if (length(list.files(path = class_path, full.names = TRUE)) == 0) next
        
        # Friendship matrix processing
        if (wave == "Wave1") {
          friend_scale <- "15"
          friend_scale_all <- "15_1"
        } else if (wave == "Wave2") {
          friend_scale <- "11"
          friend_scale_all <- "11_2h"
        } else if (wave == "Wave3") {
          friend_scale <- "19"
          friend_scale_all <- "19_3h"
        } else if (wave == "Wave4") {
          friend_scale <- "20"
          friend_scale_all <- "20_4h"
        }
        
        M <- read.csv(paste0(folder_path, "/", class, "/", friend_scale, "/", class, "_", friend_scale_all, ".csv"), row.names = 1)
        matrix_data <- as.matrix(M)
        rows <- rownames(matrix_data)
        cols <- colnames(matrix_data)
        
        # Identify all reported friendships (value = 2)
        indices <- which(matrix_data == 2, arr.ind = TRUE)
        pairs_of_friends <- t(apply(indices, 1, function(idx) { c(rows[idx[1]], cols[idx[2]]) }))
        pairs_of_friends[, 2] <- sub("^X", "", pairs_of_friends[, 2])
        
        # Sort each pair alphabetically and remove self-ties
        pairs_sorted <- t(apply(pairs_of_friends, 1, function(pair) { sort(pair) }))
        pairs_sorted <- pairs_sorted[pairs_sorted[, 1] != pairs_sorted[, 2], ]
        
        # Keep only reciprocated (mutual) friendships
        reciprocated_friends <- pairs_sorted[duplicated(pairs_sorted), ]
        if (is.vector(reciprocated_friends) && length(reciprocated_friends) == 2) {
          reciprocated_friends <- matrix(reciprocated_friends, ncol = 2, byrow = TRUE)
        } else if (!is.matrix(reciprocated_friends)) {
          reciprocated_friends <- as.matrix(reciprocated_friends)
        }
        
        # Build undirected graph and find all maximal cliques
        graph <- igraph::graph_from_edgelist(reciprocated_friends, directed = FALSE)
        cliques_list <- igraph::max_cliques(graph)
        
        # Filter cliques with at least 2 members
        cliques_filtered <- lapply(Filter(function(clique) length(clique) >= 2, cliques_list),
                                   function(clique) sort(igraph::V(graph)[clique]$name))
        
        # Store clique sizes for summary statistics
        clique_sizes <- sapply(cliques_filtered, length)
        all_clique_sizes <- c(all_clique_sizes, clique_sizes)
      }
    }
  }
  
  # Report the largest clique size in the entire dataset
  max_clique_size <- max(all_clique_sizes)
  cat("Maximum clique size in the dataset:", max_clique_size, "\n")
  
  return(list(
    max_clique_size = max_clique_size,
    all_clique_sizes = all_clique_sizes
  ))
}

