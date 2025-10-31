################################################################################
################################ PLOT GOSSIPDATA ###############################
################################################################################


# Base folder where the script is located
base_folder <- getwd() 

# Load GossipData.RData (only if it was saved previously from pre_processing_data.R)
load(file.path(base_folder, "gossipdata_analysis", "GossipData.RData"))

# Folder to save the frequency plots
output_dir <- file.path(base_folder, "gossipdata_analysis", "frequency_plots")
# Create the folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}





################################################################################
# Gossiping frequency per month
################################################################################

library(dplyr)
library(ggplot2)

gossip_df <- GossipData$data
gossip_df$SenderSize <- sapply(strsplit(gossip_df$Gossipers, ","), length)

freq_wave <- gossip_df %>%
  mutate(WaveNum = as.numeric(gsub("Wave", "", Time))) %>%
  count(WaveNum)

# wave duration (in months)
wave_duration <- c(2, 6, 12, 12)

freq_wave <- freq_wave %>%
  mutate(
    months = wave_duration,            
    rate_per_month = n / months         
  )

p1 <- ggplot(freq_wave, aes(x = factor(WaveNum), y = rate_per_month)) +
  geom_col(fill = "grey", color = "black") +  
  labs(
    x = "wave",
    y = "gossiping frequency per month"
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, max(freq_wave$rate_per_month) * 1.1)) +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )

print(p1)

ggsave(file.path(output_dir, "waves.pdf"), p1, width = 6, height = 5)







################################################################################
# Gossiping frequency for each gossiper set size
################################################################################

freq_size <- gossip_df %>% count(SenderSize)

p2 <- ggplot(freq_size, aes(x = SenderSize, y = n)) +
  geom_col(fill = "grey", color = "black") +
  labs(x = "gossiper set size", y = "gossiping frequency") +
  scale_x_continuous(breaks = sort(unique(freq_size$SenderSize))) +
  scale_y_continuous(limits = c(0, 200)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),  
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
print(p2)


ggsave(file.path(output_dir, "gossipers_size.pdf"), p2, width = 6, height = 5)






################################################################################
# Gossiping frequency in each class
################################################################################
library(dplyr)
library(ggplot2)

# Load students dataset
base_folder <- getwd() 
students_df <- read.csv(
  file.path(base_folder, "gossipdata_analysis", 
            "gossipdata", 
            "RECENS_Wired_into_each_other-Data_files", 
            "RECENS_students_dataset_w1234.csv")
)

# Obtain the class for each gossip hyperevent
get_classgroup <- function(receiver_id, time) {
  student_row <- students_df[students_df$idcode == receiver_id, ]
  if(nrow(student_row) == 1) {
    return(switch(as.character(time),
                  "1" = student_row$class_1,
                  "2" = student_row$class_2,
                  "3" = student_row$class_3,
                  "4" = student_row$class_4,
                  NA))
  } else {
    return(NA)
  }
}

GossipData$data$Class <- mapply(get_classgroup, GossipData$data$Target, 
                                as.numeric(sub("Wave", "", GossipData$data$Time)))
GossipData$data$Class <- factor(GossipData$data$Class)


# Summarize number of gossip hyperevents per class
class_summary <- GossipData$data %>%
  group_by(Class) %>%
  summarise(n_events = n()) %>%
  arrange(desc(n_events))


# Plot histogram (bar chart)
p <- ggplot(class_summary, aes(x = reorder(Class, -n_events), y = n_events)) +
  geom_col(fill = "grey", color = "black") +
  labs(
    x = "class",
    y = "gossiping frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=10),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

print(p)
ggsave(file.path(output_dir, "classes.pdf"), p, width = 6, height = 5)








################################################################################
# Gossiping frequency for each mutual friend clique size
################################################################################

all_clique_sizes <- c()
waves <- c("Wave1", "Wave2", "Wave3", "Wave4")
folders <- c("1000", "2000", "3000", "4000", "5000", "6000", "7000")

for (wave in waves) {
  for (folder in folders) {
    base_folder <- getwd()  
    folder_path <- file.path(
      base_folder,
      "gossipdata_analysis",
      "gossipdata",
      "RECENS_Wired_into_each_other-Data_files",
      "RECENS_networks_dataset_w1234",
      wave,
      folder
    )
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
      indices <- which(matrix_data == 2, arr.ind = TRUE)
      pairs_of_friends <- t(apply(indices, 1, function(idx) {c(rows[idx[1]], cols[idx[2]])}))
      pairs_of_friends[, 2] <- sub("^X", "", pairs_of_friends[, 2])
      pairs_sorted <- t(apply(pairs_of_friends, 1, function(pair) {sort(pair)}))
      pairs_sorted <- pairs_sorted[pairs_sorted[, 1] != pairs_sorted[, 2], ]
      reciprocated_friends <- pairs_sorted[duplicated(pairs_sorted), ]
      if (is.vector(reciprocated_friends) && length(reciprocated_friends) == 2) {
        reciprocated_friends <- matrix(reciprocated_friends, ncol = 2, byrow = TRUE)
      } else if (!is.matrix(reciprocated_friends)) {
        reciprocated_friends <- as.matrix(reciprocated_friends)
      }
      
      # Graph e clique
      graph <- igraph::graph_from_edgelist(reciprocated_friends, directed = FALSE)
      cliques_list <- igraph::max_cliques(graph)
      cliques_filtered <- lapply(Filter(function(clique) length(clique) >= 2, cliques_list), 
                                 function(clique) sort(igraph::V(graph)[clique]$name))
      
      clique_sizes <- sapply(cliques_filtered, length)
      all_clique_sizes <- c(all_clique_sizes, clique_sizes)
    }
  }
}


# Count the number of cliques of dimension 2, 3, 4, 5 and 6
num_2 <- sum(all_clique_sizes == 2)
num_3 <- sum(all_clique_sizes == 3)
num_4 <- sum(all_clique_sizes == 4)
num_5 <- sum(all_clique_sizes == 5)
num_6 <- sum(all_clique_sizes == 6)

cat("Number of 2-friend cliques:", num_2, "\n")
cat("Number of 3-friend cliques:", num_3, "\n")
cat("Number of 4-friend cliques:", num_4, "\n")
cat("Number of 5-friend cliques:", num_5, "\n")
cat("Number of 6-friend cliques:", num_6, "\n")




# Plot histogram
library(ggplot2)

clique_df <- data.frame(
  Size = all_clique_sizes
)

p <- ggplot(clique_df, aes(x = factor(Size))) +
  geom_bar(fill = "grey", color = "black") +
  labs(
    x = "mutual friend clique size",
    y = "frequency"
  ) +
  scale_y_continuous(
    limits = c(0, 1800),
    breaks = seq(0, 1800, by = 300),      
    minor_breaks = seq(0, 1800, by = 20)  
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )
print(p)


ggsave(file.path(output_dir, "mutual_friends.pdf"), p, width = 6, height = 5)
