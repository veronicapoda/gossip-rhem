# ========================================================================
# DATA PRE-PROCESSING: GOSSIP HYPEREVENTS CONSTRUCTION
# ========================================================================


# Load necessary libraries
library(igraph)   
library(dplyr)   

# Define waves and folders
waves <- c("Wave1", "Wave2", "Wave3", "Wave4")
folders <- c("1000", "2000", "3000", "4000", "5000", "6000", "7000")


# Initialize empty data structures
Gossip_data <- data.frame(
  Gossiper1 = character(),
  Gossiper2 = character(),
  Target = character(),
  Time = character(),
  stringsAsFactors = FALSE
)
gender_data_list_all <- list()




# Main loop over waves and classes
for (wave in waves) {
  cat("\n")
  cat("\n")
  cat("\n")
  cat("Wave:", wave, "\n")  
  
  # Loop through each folder
  for (folder in folders) {
    # Print the current folder
    cat("\n")
    cat("\n")
    cat("\n")
    cat("Current folder:", folder, "\n")
    
    # Construct folder path for this wave and folder
    folder_path <- file.path("gossipdata_analysis",
                             "gossipdata",
                             "RECENS_Wired_into_each_other-Data_files",
                             "RECENS_networks_dataset_w1234",
                             wave,
                             folder)
    
    # Get subfolders (classes) inside the folder
    classes_in_folder <- list.files(path = folder_path, full.names = FALSE)
    
    # Keep only classes whose name starts with the same digit as the folder
    folder_prefix <- substr(folder, 1, 1)
    classes_in_folder <- classes_in_folder[substr(classes_in_folder, 1, 1) == folder_prefix]
    
    # Loop through each subfolder (class) inside the folder
    for (class in classes_in_folder) {
      class_path <- file.path(folder_path, class)
      
      # Skip empty folders
      if (length(list.files(path = class_path, full.names = TRUE)) == 0) {
        cat("\nSkipping empty class:", class, "\n")
        next
      }
      
      cat("\nInspecting class:", class, "\n")
      
      
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
      # Path
      csv_path <- file.path(folder_path, class, friend_scale, paste0(class, "_", friend_scale_all, ".csv"))
      M <- read.csv(csv_path, row.names = 1)
      matrix_data <- as.matrix(M)
      rows <- rownames(matrix_data)
      cols <- colnames(matrix_data)
      
      # Find reciprocated friendships (value == 2)
      indices <- which(matrix_data == 2, arr.ind = TRUE)
      pairs_of_friends <- t(apply(indices, 1, function(idx) {c(rows[idx[1]], cols[idx[2]])}))
      pairs_of_friends[, 2] <- sub("^X", "", pairs_of_friends[, 2])
      
      # Normalize pairs and remove self-loops
      pairs_sorted <- t(apply(pairs_of_friends, 1, function(pair) {sort(pair)}))
      pairs_sorted <- pairs_sorted[pairs_sorted[, 1] != pairs_sorted[, 2], ]
      
      # Keep only reciprocated pairs
      reciprocated_friends <- pairs_sorted[duplicated(pairs_sorted), ]
      
      # Ensure matrix format even for single pair
      if (is.vector(reciprocated_friends) && length(reciprocated_friends) == 2) {
        reciprocated_friends <- matrix(reciprocated_friends, ncol = 2, byrow = TRUE)
      } else if (!is.matrix(reciprocated_friends)) {
        reciprocated_friends <- as.matrix(reciprocated_friends)
      }
      
      # Build undirected friendship graph
      graph <- graph_from_edgelist(reciprocated_friends, directed = FALSE)
      
      # Find maximal cliques (fully connected subgroups)
      cliques_list <- max_cliques(graph)
      cliques_filtered <- lapply(Filter(function(clique) length(clique) >= 2, cliques_list), 
                                 function(clique) sort(V(graph)[clique]$name))
      # Summary
      cat("FRIENDSHIP \n")
      cat("Number of pairs of reciprocated friends:", nrow(reciprocated_friends), "\n")
      cat("Percentage of reciprocated friends:", round(nrow(reciprocated_friends)/nrow(pairs_sorted)*100, 0),"%", "\n")
      
      # Print cliques
      for (clique in cliques_filtered) {
        cat(paste(clique, collapse = ", "), "\n")
      }
      
      # Negative gossip matrix processing
      if (wave == "Wave1") {
        neg_gossip <- "33"
        neg_gossip_all <- "33_3"
      } else if (wave == "Wave2") {
        neg_gossip <- "32"
        neg_gossip_all <- "32_3_2h"
      } else if (wave == "Wave3") {
        neg_gossip <- "39"
        neg_gossip_all <- "39_3_3h"
      } else if (wave == "Wave4") {
        neg_gossip <- "38"
        neg_gossip_all <- "38_3_4h"
      }
      # Path
      ngossip_path <- file.path(folder_path, class, neg_gossip, paste0(class, "_", neg_gossip_all, ".csv"))
      NGossip <- read.csv(ngossip_path, row.names = 1)
      gossip_data <- as.matrix(NGossip)
      rows <- rownames(gossip_data)
      cols <- colnames(gossip_data)
      
      # Extract negative gossip pairs (value == 1)
      indices <- which(gossip_data == 1, arr.ind = TRUE)
      gossip_pairs <- t(apply(indices, 1, function(idx) {c(rows[idx[1]], cols[idx[2]])}))
      gossip_pairs[, 2] <- sub("^X", "", gossip_pairs[, 2])
      colnames(gossip_pairs) <- c("Gossiper", "Target")
      
      cat("NEGATIVE GOSSIP \n")
      cat("Total number of Gossiper-Target pairs:", nrow(gossip_pairs), "\n")
      
      # Group gossipers by target
      gossip_pairs <- as.data.frame(gossip_pairs)
      gossip_pairs_grouped <- gossip_pairs %>%
        group_by(Target) %>%
        summarise(Gossiper = paste(Gossiper, collapse = ", "), .groups = "drop")
      
      cat("Unique targets grouped:", nrow(gossip_pairs_grouped), "\n")
      
      # Gossip Hyperevents construction
      target_gossiper_data <- data.frame(gossip_pairs_grouped)
      cliques_sets <- lapply(cliques_filtered, function(clique) sort(unique(clique)))
      valid_cliques <- list()
      
      for (i in 1:nrow(target_gossiper_data)) {
        gossiper_vector <- sort(unique(unlist(strsplit(target_gossiper_data$Gossiper[i], ", "))))
        matching_cliques <- lapply(cliques_sets, function(clique) {if (all(clique %in% gossiper_vector)) return(clique) else return(NULL)})
        matching_cliques <- matching_cliques[!sapply(matching_cliques, is.null)]
        
        if (length(matching_cliques) > 0) {
          for (clique in matching_cliques) {
            valid_cliques <- append(valid_cliques, list(data.frame(
              Gossipers = paste(clique, collapse = ", "),
              Target = target_gossiper_data$Target[i],
              Time = wave,
              Class = class,      
              Folder = folder,    
              stringsAsFactors = FALSE
            )))
          }
        }
      }
      
      # Combine valid hyperevents
      if (length(valid_cliques) > 0) {
        Gossip_cliques_data <- do.call(rbind, valid_cliques)
        Gossip_data <- rbind(Gossip_data, Gossip_cliques_data)
        }
      
      # Process gender data
      gender_path <- file.path("gossipdata_analysis",
                               "gossipdata",
                               "RECENS_Wired_into_each_other-Data_files",
                               "RECENS_students_dataset_w1234.csv")
      gender <- read.csv(gender_path)
      process_gender_data <- function(class, gender_column) {
        class_prefix <- substr(as.character(class), 1, 2)  # Extract the first two digits from the class (e.g., 6100 -> "61")
        filtered_data <- gender %>%
          filter(substr(as.character(idcode), 1, 2) == class_prefix) %>% 
          select(idcode, all_of(gender_column))  
        colnames(filtered_data)[colnames(filtered_data) == gender_column] <- "gender"  
        return(filtered_data)
      }
      class_data <- process_gender_data(class, "gender_1")
      gender_data_list_all[[paste0("gender_", substr(as.character(class), 1, 2), "00")]] <- class_data
    }
  }
}






# Save final data
GossipData <- list(
  data = Gossip_data,
  info = gender_data_list_all
)

# Path
base_folder <- getwd()  
save_folder <- file.path(base_folder, "gossipdata_analysis")
if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

# Save gossip hyperevents as RData file
save(GossipData, file = file.path(save_folder, "GossipData.RData"))

cat("\n GossipData saved successfully!\n")
