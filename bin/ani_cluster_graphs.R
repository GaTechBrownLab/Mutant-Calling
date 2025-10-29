#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)
library(svglite)
library(forcats)

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]
ani_cluster_input <- args[2]
ani_clusters <- args[3]
genome_annotations_classified <- args[4]
functions <- args[5]

# Read in environmental counts
ani_cluster_input_df <- read.csv(ani_cluster_input)

counts <- ani_cluster_input_df %>%
  count(Threshold, Cluster_Status, name = "Cluster_count")

ani_cluster_input_df <- ani_cluster_input_df %>%
  mutate(New_envs = case_when(
    Unique_Group_Count == 1 ~ "1",
    Unique_Group_Count == 2 ~ "2",
    Unique_Group_Count > 2  ~ "greater than 2"
  ))

env_counts <- ani_cluster_input_df %>%
  count(Threshold, Cluster_Status, New_envs, name = "Env_count")

merged_df <- left_join(env_counts, counts, by = c("Threshold", "Cluster_Status"))

merged_df$Proportion <- merged_df$Env_count / merged_df$Cluster_count

merged_df_no_mixed <- merged_df %>%
  filter(!grepl("Mixed", Cluster_Status))

prop_plot <- ggplot(merged_df_no_mixed,
                    aes(x = factor(Threshold), y = Proportion)) +
  geom_bar(aes(fill = New_envs),
           position = "stack",
           stat = "identity",
           color = "black") +
  geom_text(aes(label = Env_count,
                group = New_envs),
            position = position_stack(vjust = 0.5), # center in stack
            color = "black",
            size = 3) +
  facet_grid(~ Cluster_Status) +
  xlab("ANI clusters") +
  scale_fill_discrete(name = "N of unique envs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


prop_bar_plot <- paste0(gene, "_ANI_distributions_prop.svg")

ggsave(prop_bar_plot, plot = prop_plot, width = 10, height = 7.5, units = "in")


plot <- ggplot(ani_cluster_input_df, aes(x = factor(Threshold), y = Unique_Group_Count, group = Cluster_Status)) +
  geom_jitter(color = "black", pch = 21, aes(fill = Cluster_Status, alpha = Cluster_Status), width = 0.2, height = 0.2)  +
  scale_fill_manual(values = c("Functional" = "blue", "No function" = "red", "Mixed" = "#883268"), name = "Function Status") +
  scale_alpha_manual(
    values = c("No function" = 1, "Functional" = 0.3, "Mixed" = 0.3),
    guide = "none"
  ) +
  xlab("ANI clusters") +
  ylab("Number of unique environments") +
  theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

bar_plot <- paste0(gene, "_ANI_distributions.svg")

ggsave(bar_plot, plot = plot, width = 10, height = 7.5, units = "in")

check_no_function <- function(df) {
  tbl <- table(df$Cluster_Status, df$Unique_Group_Count)
  test <- chisq.test(tbl)
  
  resid <- test$residuals["No function", , drop = FALSE]

  tibble(
    Environment = colnames(tbl),
    residual = as.numeric(resid),
    enriched = residual > 2,
    depleted = residual < -2,
    statistic = test$statistic,
    df = test$parameter,
    p_value = test$p.value
  )
}

ani_cluster_input_df$Unique_Group_Count <- as.character(ani_cluster_input_df$Unique_Group_Count)

ani_cluster_results <- ani_cluster_input_df %>%
  group_by(Threshold) %>%
  group_modify(~ check_no_function(.x))

ani_cluster_results_name <- paste0(gene, "_ANI_clusters_chi.csv")

write.csv(ani_cluster_results, ani_cluster_results_name, row.names = FALSE)


# ----------------------------
# ANI clusters heatmap
# ----------------------------

ani_clusters_df <- read.csv(ani_clusters)

functions_df <- read.csv(functions)

functions_df_status <- functions_df %>%
  filter(Mut_Status == "No function")

functions_df_status_merge <- inner_join(functions_df_status, ani_clusters_df, by = c("Genome"))

ani_cluster_input_df_thresh <- ani_cluster_input_df %>%
  filter(Threshold == 99.955)

ani_cluster_input_df_status <- ani_cluster_input_df_thresh %>%
  filter(Cluster_Status == "No function")

merged_no_funct_df <- inner_join(functions_df_status_merge, ani_cluster_input_df_status, by = c("Cluster", "Threshold"))

test_name <- paste0(gene, "_cluster_status_test.csv")

write.csv(merged_no_funct_df, test_name, row.names = FALSE)

meta <- read.csv(genome_annotations_classified)

env_merged_no_funct <- left_join(merged_no_funct_df, meta, by = "Genome")

env_merged_no_funct_name <- paste0(gene, "_no_funct_clusters_env.csv")

write.csv(env_merged_no_funct, env_merged_no_funct_name, row.names = FALSE)

env_merged_no_funct_tbl <- table(env_merged_no_funct$Cluster, env_merged_no_funct$Group)

env_merged_no_funct_df <- as.data.frame(env_merged_no_funct_tbl)

cluster_totals_df <- env_merged_no_funct_df %>%
  group_by(Var1) %>%
  summarise(Total = sum(Freq))

env_merged_no_funct_df <- env_merged_no_funct_df %>%
  left_join(cluster_totals_df, by = "Var1")

env_merged_no_funct_df$Proportion <- env_merged_no_funct_df$Freq / env_merged_no_funct_df$Total

all_envs <- c("CF", "Clinical", "Clinical, Unknown", "Missing", "Environmental",
  "Human-associated environmental", "Lab", "Animal")

env_merged_no_funct_complete <- env_merged_no_funct_df %>%
  complete(
    Var1, Var2 = all_envs,
    fill = list(Freq = 0, Proportion = 0)
  )

env_heatmap <- ggplot(env_merged_no_funct_complete, aes(x = Var2, y = Var1, fill = Proportion)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    name = "Proportion"
  ) +
  geom_text(aes(label = Freq), color = "white", size = 3) +
  labs(x = "Environment", y = "Lineage") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

env_heatmap_file <- paste0(gene, "_environment_heatmap.svg")

ggsave(env_heatmap_file, plot = env_heatmap, width = 10, height = 10, units = "in")

# Functional

ani_cluster_input_df_status_funct <- ani_cluster_input_df_thresh %>%
  filter(Cluster_Status == "Functional")

functions_df_status_funct <- functions_df %>%
  filter(Mut_Status == "Functional")

functions_df_status_merge_funct <- inner_join(functions_df_status_funct, ani_clusters_df, by = c("Genome"))

merged_funct_df <- inner_join(functions_df_status_merge_funct, ani_cluster_input_df_status_funct, by = c("Cluster", "Threshold"))

env_merged_funct <- left_join(merged_funct_df, meta, by = "Genome")

env_merged_funct_tbl <- table(env_merged_funct$Cluster, env_merged_funct$Group)

env_merged_funct_df <- as.data.frame(env_merged_funct_tbl)

cluster_totals_df <- env_merged_funct_df %>%
  group_by(Var1) %>%
  summarise(Total = sum(Freq))

env_merged_funct_df <- env_merged_funct_df %>%
  left_join(cluster_totals_df, by = "Var1")

env_merged_funct_df$Proportion <- env_merged_funct_df$Freq / env_merged_funct_df$Total

env_merged_funct_complete <- env_merged_funct_df %>%
  complete(
    Var1, Var2 = all_envs,
    fill = list(Freq = 0, Proportion = 0)
  )

env_heatmap_funct <- ggplot(env_merged_funct_complete, aes(x = Var2, y = Var1, fill = Proportion)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    name = "Proportion"
  ) +
  geom_text(aes(label = Freq), color = "white", size = 3) +
  labs(x = "Environment", y = "Lineage") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

env_heatmap_funct_file <- paste0(gene, "_environment_heatmap_funct.svg")

ggsave(env_heatmap_funct_file, plot = env_heatmap_funct, width = 10, height = 10, units = "in")