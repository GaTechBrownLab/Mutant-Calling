#!/usr/bin/env Rscript

library(ape)
library(dplyr)
library(tibble)
library(geiger)
library(phylometrics)
library(extrafont)
library(randomcoloR)
library(stats)
library(phangorn)
library(tidyr)
library(ggplot2)
library(gt)

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
tree <- args[1]
tree_root <- args[2]
base <- args[3]
functions_incomplete <- args[4]
functions_pres_abs_incomplete <- args[5]
fastani <- args[6]
genome_annotations_classified <- args[7]
lineage_probability <- args[8]
no_funct_clusters_env <- args[9]

# Import fonts
font_import()

# Load fonts
loadfonts()

# ----------------------------
# Read in tree and perform basic tree analyses
# ----------------------------

if (base == "PA1430") {
  gene <- "LasR"
}

if (base == "PA3168") {
  gene <- "GyrA"
}

#Read in tree
tree <- ape::read.tree(tree)

#Root tree on PA7
rooted_tree <- root(tree, outgroup = tree_root, resolve.root = TRUE)

#Remove tips with branch length 0
tips_with_zero_length <- function(tree) {
  zero_length_tips <- character(0)
  
  for (i in tree$tip.label) {
    if (tree$edge.length[which(tree$tip.label == i)] == 0) {
      zero_length_tips <- c(zero_length_tips, i)
    }
  }
  
  return(zero_length_tips)
}

# Identify tips with branch length of 0
zero_length_tips <- tips_with_zero_length(rooted_tree)

# Subset trees
rooted_tree <- drop.tip(rooted_tree, zero_length_tips, trim.internal = TRUE)

tree <- drop.tip(tree, zero_length_tips, trim.internal = TRUE)

# Pull tip labels from tree
all_tips <- rooted_tree$tip.label
tips_df <- as.data.frame(all_tips)

# Rename columns
colnames(tips_df)[1] <- "Genome"

# Set seed
set.seed(123)

# Read in functions data
mol_data <- read.csv(functions_incomplete)

# Identify functions data not in tree
non_mol <- tips_df %>% filter(!Genome %in% mol_data$Genome)

tips_to_drop <- non_mol$Genome

# Drop tips that are not in functions data
pruned_tree_rooted <- drop.tip(rooted_tree, tips_to_drop, trim.internal = TRUE)

pruned_tree_unrooted <- drop.tip(tree, tips_to_drop, trim.internal = TRUE)

# ----------------------------
# Create an ultrametric tree using APE
# ----------------------------

# Calibrate the tree
mycalibration <- makeChronosCalib(pruned_tree_rooted, node="root")

relaxed_0 <- chronos(pruned_tree_rooted, lambda=0)

# Breakpoint - saving tree ##
relaxed_0_file <- paste0(base, "_tree.rds")
saveRDS(relaxed_0, relaxed_0_file)

# Capture output
chronogram_output <- capture.output(chronos(pruned_tree_rooted, lambda=0))

# Extract the log-likelihood value using regex
log_lik_line <- grep("log-Lik =", chronogram_output, value = TRUE)
log_lik_value <- sub(".*log-Lik = (-?[0-9.]+).*", "\\1", log_lik_line)

# Write to a file
log_lik_file <- paste0(base, "_log_likelihood.txt")

write(log_lik_value, file = log_lik_file)

# ----------------------------
# Create state table
# ----------------------------

# Create the traits table in the same order as the tips
tree_labels <- relaxed_0$tip.label
tree_labels <- as.data.frame(tree_labels)
colnames(tree_labels) <- c('Genome')
tree_labels$id  <- 1:nrow(tree_labels)

# Read in presence absence table
tree_data <- read.csv(functions_pres_abs_incomplete)

# Merge trait table with presence absence data
trait_table <- merge(tree_labels,tree_data, by = "Genome")
trait_table <- trait_table[order(trait_table$id), ]
row.names(trait_table) <- NULL

# Obtain all the positive states (non-zero)
pos_state <- subset(trait_table, Mut_Status != 0)

## Reformat this
traits <- trait_table[['Mut_Status']]

pos <- pos_state[['Genome']]

# Check the number of tips
length(relaxed_0$tip.label) == length(traits)

# ----------------------------
# Phylometrics values
# ----------------------------

# Values for PA tree
noto_cheat <- treestat(relaxed_0, stlist=pos, state=traits, func=noto)
tars_cheat <- treestat(relaxed_0, stlist=pos, state=traits, func=tars)
sscd_cheat <- treestat(relaxed_0, stlist=pos, state=traits, func=sscd)
fpd_cheat <- treestat(relaxed_0, stlist=pos, state=traits, func=fpd)

# Create a data frame with the results
phylometrics_table <- data.frame(
  Gene = base,
  noto = noto_cheat,
  tars = tars_cheat,
  sscd = sscd_cheat,
  fpd = fpd_cheat
)

phylometrics_file <- paste0(base, "_phylometrics_table.csv")
write.csv(phylometrics_table, phylometrics_file, row.names = FALSE)

# ----------------------------
# Identify environmental origins
# ----------------------------

# Remove 0 length tips
mol_data <- mol_data %>% filter(!Genome %in% zero_length_tips)

# Transform metadata
rownames(mol_data) <- mol_data[,1]
mol_data[,1] <- NULL
LasR_meta <- setNames(mol_data[,1],rownames(mol_data))

#Set metadata colors
color_map <- c("Functional" = "blue", "No function" = "red")

#Map colors to tips
tip_labels <- relaxed_0$tip.label

tip_states <- LasR_meta[tip_labels]

tip_colors <- color_map[as.character(tip_states)]

# ----------------------------
# Ancestral State Reconstruction
# ----------------------------

# Read in functions file
funct_no <- read.csv(functions_incomplete)

# Drop rows where "Genome" is in zero_length_tips
funct_no <- funct_no %>% filter(!Genome %in% zero_length_tips)

# Remove duplicate rows
funct_no <- funct_no %>% distinct()

# Set "Genome" as row names and drop the original "Genome" column
rownames(funct_no) <- funct_no$Genome
funct_no$Genome <- NULL

# Reformat
funct_no_ASR <- setNames(funct_no[,1],rownames(funct_no))

# Run ancestral state reconstruction under the ARD model

internal_nodes <- (length(relaxed_0$tip.label) + 1):(length(relaxed_0$tip.label) + relaxed_0$Nnode)
relaxed_0$node.label <- as.character(internal_nodes)

fit_ARD_ASR <- ace(funct_no_ASR, relaxed_0, model = 'ARD', type="discrete")

st_ASR <- fit_ARD_ASR$lik.anc

# Get the length of all internal nodes, and only keep probabilities for internal nodes
internal_nodes <- (length(relaxed_0$tip.label) + 1):(length(relaxed_0$tip.label) + relaxed_0$Nnode)
st_ASR <- st_ASR[as.character(internal_nodes), , drop = FALSE]

pa <- read.csv(functions_pres_abs_incomplete)

rownames(pa) <- pa$Genome
pa$Genome <- NULL

# Set colors
LasR_pa<-setNames(pa[,1],rownames(pa))

color_map <- c('0' = 'blue', '1' = 'red')

tip_labels <- relaxed_0$tip.label

tip_states <- LasR_pa[tip_labels]

tip_colors <- color_map[as.character(tip_states)]

# Plot
ASR_name_file <- paste0(base, "_ASR_tree.svg")

gene_label_wt <- paste0(gene, " WT")
gene_label_mut <- paste0(gene, " Mutant")

svg(filename = ASR_name_file, width = 10, height = 10)

plot.new()
plot(relaxed_0, type = "fan", lwd = 0.5, show.tip.label = FALSE)
nodelabels(pie = st_ASR, piecol = c("blue", "red"), cex = 0.15)
tiplabels(pch = 21, bg = tip_colors, col = "black", cex = 0.5)
legend("bottomright", legend=c(gene_label_wt, gene_label_mut),
        col=c("blue", "red"), pch=19, box.lty=0, text.font=1)
add.scale.bar(length = 0.1, x = -1, y = -1)

dev.off ()

# ----------------------------
# Identify lineages based on ancestral state reconstruction
# ----------------------------

# Covert st to dataframe
st_df <- as.data.frame(st_ASR)

# Identify internal nodes
num_tips <- length(relaxed_0$tip.label)

internal_nodes <- (num_tips + 1):max(relaxed_0$edge)

# Assign nodes to probabilities (should be in same order, as that is how the example tree is read in)
st_df$Node <- internal_nodes

# Rename columns
colnames(st_df) <- c("State0_Probability", "State1_Probability", "Node")

# Drop those with a state 1 (mutant) probability less than 0.9
lineage_probability <- as.numeric(lineage_probability)

identify_lineages <- function(state_prob, type) {
  st_df_lin <- st_df %>% filter(.data[[state_prob]] >= lineage_probability)

  # Save lineages
  node_tips_df_unique_file <- paste0(base, "_lineages_", type, ".csv")
  write.csv(st_df_lin, node_tips_df_unique_file, row.names = FALSE)

  # Get descedents (children) of all nodes in tree, merging with tip names to only get direct downstream tips
  node_tips <- list()

  for (node in st_df_lin$Node) {
  
  tips_for_node_list <- Descendants(relaxed_0, node, type = "children")

  tips_for_node <- unlist(tips_for_node_list)

  tips_for_node_labels <- relaxed_0$tip.label[tips_for_node][!is.na(relaxed_0$tip.label[tips_for_node])]

  node_tips[[paste0(node)]] <- tips_for_node_labels
  
  }

  # Convert the list to a dataframe
  node_tips_df <- stack(node_tips)

  # Subset dataframe
  colnames(node_tips_df) <- c("Genome", "Node")

  # Convert to character, and numeric
  node_tips_df$Node <- as.numeric(as.character(node_tips_df$Node))

  # Create clusters from original node list - if they are sequential, assign to same cluster
  st_df_lin$Node <- as.numeric(st_df_lin$Node)

  st_df_lin_clusters <- st_df_lin %>%
      arrange(Node) %>%
      mutate(Cluster = cumsum(c(1, diff(Node) > 2)))

  # Subset
  st_df_lin_clusters <- st_df_lin_clusters[c("Node","Cluster")]

  # Drop duplicate rows
  st_df_lin_clusters_distinct <- st_df_lin_clusters %>% distinct()

  # Left merge, assigning all genomes to a cluster based on the nodes
  node_tips_df_cluster <- left_join(node_tips_df, st_df_lin_clusters_distinct, by = "Node")

  # Remove duplicates
  node_tips_df_unique <- node_tips_df_cluster[!duplicated(node_tips_df_cluster), ]

  # Set lineage status to lineage
  node_tips_df_unique$Lineage_Status <- "Lineage"

  return(list(clusters = st_df_lin_clusters, tips = node_tips_df_unique))
}

cheat_lineages <- identify_lineages("State1_Probability", "Mutant")

WT_lineages <- identify_lineages("State0_Probability", "WT")

# ----------------------------
# Associate lineages with colors for cheat graphing
# ----------------------------

# Subset table
lineage_clusters_sub <- cheat_lineages$tips[c("Genome","Cluster")]

# Transform metadata
rownames(lineage_clusters_sub) <- lineage_clusters_sub[,1]
lineage_clusters_sub[,1] <- NULL

# Get number of clusters
num_clusters <- length(unique(lineage_clusters_sub$Cluster))

# Generate random, unrelated colors
cluster_colors <- setNames(distinctColorPalette(num_clusters), 
                            unique(lineage_clusters_sub$Cluster))

tip_labels <- relaxed_0$tip.label
tip_clusters <- lineage_clusters_sub[tip_labels, "Cluster"]  # Get clusters for the tip labels

# Map colors to the clusters
lineage_tip_colors <- cluster_colors[as.character(tip_clusters)]

# ----------------------------
# Node depth
# ----------------------------

# Node depths are proportional to the number of tips descending from each node
node_depth <- node.depth(relaxed_0, method = 1)

# Convert node_depth vector into a dataframe
node_depth_df <- data.frame(Node = 1:length(node_depth), Depth = node_depth)

# Ensure the Node column is numeric
node_depth_df$Node <- as.numeric(node_depth_df$Node)

# Keep the lowest (deepest) node per cluster
st_df_lin_clusters_lowest <- cheat_lineages$clusters %>%
  group_by(Cluster) %>%
  filter(Node == min(Node)) %>%
  ungroup()

# Get the depth of all of the deepest nodes
deep_node_depth <- left_join(st_df_lin_clusters_lowest, node_depth_df, by = "Node")

# Calculate average cluster depth
avg_depth <- mean(deep_node_depth$Depth)

# Save average node depth
node_depth_table <- data.frame(
  Gene = base,
  Average_Node_Depth = avg_depth
)

node_depth_file <- paste0(base, "_node_depth_table.csv")
write.csv(node_depth_table, node_depth_file, row.names = FALSE)

# ----------------------------
# ANI comparison
# ----------------------------

# Read in genomavars
ANI <- read.table(fastani, sep = ",")

colnames(ANI)[1] <- "Genome1"
colnames(ANI)[2] <- "Genome2"
colnames(ANI)[3] <- "ANI"

ANI$Genome1 <- gsub(".fna", "", ANI$Genome1)
ANI$Genome2 <- gsub(".fna", "", ANI$Genome2)

ANI <- ANI[c("Genome1", "Genome2", "ANI")]

merged_ANI1 <- left_join(ANI, cheat_lineages$tips, by = c("Genome1" = "Genome"))
merged_ANI2 <- left_join(merged_ANI1, cheat_lineages$tips, by = c("Genome2" = "Genome"))

merged_ANI2 <- merged_ANI2 %>% drop_na()

filtered_ANI <- merged_ANI2 %>%
  filter(Cluster.x == Cluster.y)

filtered_ANI_lowest <- filtered_ANI %>%
  group_by(Cluster.x) %>%
  slice_min(ANI, with_ties = FALSE) %>%
  ungroup()

filtered_ANI_lowest_sub <- filtered_ANI_lowest[c("Cluster.x", "Lineage_Status.x", "ANI")]

colnames(filtered_ANI_lowest_sub)[1] <- "Cluster"
colnames(filtered_ANI_lowest_sub)[2] <- "Lineage_Status"

ANI_low_file <- paste0(base, "_lowest_ANI_cluster.csv")
write.csv(filtered_ANI_lowest_sub, ANI_low_file, row.names = FALSE)

ani_plot <- ggplot(filtered_ANI_lowest_sub, aes(x = ANI)) +
  geom_histogram(fill = "red", color = "black", bins = 50) +
  geom_vline(xintercept = 99.99, color = "black", linetype = "dashed", size = 1) +
  xlab("ANI") +
  ylab("Count") +
  theme_bw()

ani_plot_name <- paste0(base, "_ani_distribution.svg")
ggsave(ani_plot_name, plot = ani_plot, width = 10, height = 7.5, units = "in")

avg_ANI <- mean(filtered_ANI_lowest$ANI)

# Save average ANI
avg_ANI_table <- data.frame(
  Gene = base,
  Average_ANI = avg_ANI
)

avg_ANI_file <- paste0(base, "_avg_ANI_table.csv")
write.csv(avg_ANI_table, avg_ANI_file, row.names = FALSE)

# ----------------------------
# Associate cheat and WT lineages with environmental origins
# ----------------------------

lin_env <- function(input, mut_status_int, type, clin) {
  # Join cheat lineages with tree data
  pres_abs_merge_lin <- left_join(input, tree_data, by = "Genome")

  # Drop those that are not correct group
  pres_abs_merge_lin_no_funct <- pres_abs_merge_lin[pres_abs_merge_lin$Mut_Status == mut_status_int, ]

  # Save mutant lineages
  node_tips_df_unique_file <- paste0(base, "_lineages_nodes_", type, ".csv")
  write.csv(pres_abs_merge_lin_no_funct, node_tips_df_unique_file, row.names = FALSE)

  # Read in environmental data
  meta <- read.csv(genome_annotations_classified)

  # Join lineages with environmental data
  lineage_env_merge <- left_join(pres_abs_merge_lin_no_funct, meta, by = "Genome")

  lineage_env_merge <- lineage_env_merge %>% 
    filter(!(Group %in% c("Animal", "Lab")))

  # If combining clinical, do
  if (clin == "clinical") {
    lineage_env_merge <- lineage_env_merge %>%
      mutate(
        Group = ifelse(Group %in% c("Clinical", "Clinical, Unknown"),
                      "Clinical all",
                      Group))

    print(lineage_env_merge)
  }
  
  # Make a table with counts per environment
  lineage_env_merge_tbl <- table(lineage_env_merge$Cluster, lineage_env_merge$Group)

  lineage_env_merge_df <- as.data.frame(lineage_env_merge_tbl)

  # Group
  cluster_totals_df <- lineage_env_merge_df %>%
  group_by(Var1) %>%
  summarise(Total = sum(Freq))

  lineage_env_merge_df <- lineage_env_merge_df %>%
  left_join(cluster_totals_df, by = "Var1")

  # Calculate proportion
  lineage_env_merge_df$Proportion <- lineage_env_merge_df$Freq / lineage_env_merge_df$Total

  if (clin == "clinical") {
    all_envs <- c("CF", "Clinical all", "Missing", "Environmental",
    "Human-associated environmental")
  } else {
    all_envs <- c("CF", "Clinical", "Clinical, Unknown", "Missing", "Environmental",
    "Human-associated environmental")
  }

  # Include all environments if missing set to 0
  lineage_env_merge_df_complete <- lineage_env_merge_df %>%
  complete(
      Var1, Var2 = all_envs,
      fill = list(Freq = 0, Proportion = 0)
  )

  # ----------------------------
  # Environmental heatmap for lineages
  # ----------------------------

  env_heatmap <- ggplot(lineage_env_merge_df_complete, aes(x = Var2, y = Var1, fill = Proportion)) +
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

      env_heatmap_file <- paste0(base, "_environment_heatmap_", type, "_", clin, ".svg")
      ggsave(env_heatmap_file, plot = env_heatmap, width = 10, height = 10, units = "in")

  lineage_env_merge_df_complete$Status <- type

  # ----------------------------
  # Calculating shannon diversity
  # ----------------------------

  lin_env_non_zero <- lineage_env_merge_df_complete %>% 
    filter(Proportion > 0)

  # Compute Shannon components
  lin_env_non_zero$H_component <- lin_env_non_zero$Proportion * log(lin_env_non_zero$Proportion)

  # Compute Shannon index for each lineage
  shannons_diversity <- lin_env_non_zero %>%
    group_by(Var1) %>%
    summarise(H = -sum(H_component))

  # Save mutant lineages
  shannons_diversity_file <- paste0(base, "_shannon_diversity_", type, "_", clin, ".csv")
  write.csv(shannons_diversity, shannons_diversity_file, row.names = FALSE)
  
  return(list(counts = lineage_env_merge_df_complete, shannon = shannons_diversity))

}

cheat <- lin_env(cheat_lineages$tips, 1, "Mutant", "all")
WT <- lin_env(WT_lineages$tips, 0, "WT", "all")

cheat_clin_combined <- lin_env(cheat_lineages$tips, 1, "Mutant", "clinical")
WT_clin_combined <- lin_env(WT_lineages$tips, 0, "WT", "clinical")

# Run Wilcoxon rank sum test
wilcox_results <- wilcox.test(cheat$shannon$H, WT$shannon$H)

wilcox_results_df <- data.frame(
  Test = wilcox_results$method,
  W = as.numeric(wilcox_results$statistic),
  p_value = wilcox_results$p.value,
  Alternative = wilcox_results$alternative
)

wilcox_file <- paste0(base, "_wilcox.csv")
write.csv(wilcox_results_df, wilcox_file, row.names = FALSE)

# Run Wilcoxon rank sum test for combined
wilcox_results_combined <- wilcox.test(cheat_clin_combined$shannon$H, WT_clin_combined$shannon$H)

wilcox_results_comined_df <- data.frame(
  Test = wilcox_results_combined$method,
  W = as.numeric(wilcox_results_combined$statistic),
  p_value = wilcox_results_combined$p.value,
  Alternative = wilcox_results_combined$alternative
)

wilcox_comined_file <- paste0(base, "_wilcox_clinical.csv")
write.csv(wilcox_results_comined_df, wilcox_comined_file, row.names = FALSE)


# ----------------------------
# Calculating odds ratio
# ----------------------------

combined_cheat_WT <- rbind(cheat$counts, WT$counts)

print(combined_cheat_WT)

# Sum frequency for environment across lineages of mutant or WT
env_summary <- combined_cheat_WT %>% 
  group_by(Var2, Status) %>% 
  summarise(Freq_sum = sum(Freq, na.rm = TRUE), .groups = "drop")

# Reshape to wide
env_wide <- env_summary %>%
  pivot_wider(names_from = Status,
              values_from = Freq_sum,
              values_fill = 0) %>%
  rename(Mutant = Mutant, WT = WT)

# Compute totals across all environments
total_mutant <- sum(env_wide$Mutant)
total_wt     <- sum(env_wide$WT)

# Compute odds ratio for each environment
odds_ratio_all <- env_wide %>%
  mutate(
    Em = Mutant,
    Ew = WT,
    NEm = total_mutant - Mutant,
    NEw = total_wt - WT,
    m_ratio = Em / NEm,
    w_ratio = Ew / NEw,
    OR = ( Em / NEm ) / ( Ew / NEw )
  ) %>%
  select(Var2, Em, Ew, NEw, NEm, m_ratio, w_ratio,OR)

colnames(odds_ratio_all)[1] <- "Environment"

odds_ratio_all_file <- paste0(base, "_odds_ratio.csv")
write.csv(odds_ratio_all, odds_ratio_all_file, row.names = FALSE)

or_gt <- odds_ratio_all %>%
  gt()

or_gt_file <- paste0(base, "_odds_ratio.html")
gtsave(or_gt, or_gt_file)

# ----------------------------
# Calculating odds ratio for combined clinical
# ----------------------------

combined_cheat_WT_clinical <- rbind(cheat_clin_combined$counts, WT_clin_combined$counts)
print(combined_cheat_WT_clinical)

# Sum frequency for environment across lineages of mutant or WT
env_summary_clinical <- combined_cheat_WT_clinical %>% 
  group_by(Var2, Status) %>% 
  summarise(Freq_sum = sum(Freq, na.rm = TRUE), .groups = "drop")

# Reshape to wide
env_wide_clinical <- env_summary_clinical %>%
  pivot_wider(names_from = Status,
              values_from = Freq_sum,
              values_fill = 0) %>%
  rename(Mutant = Mutant, WT = WT)

# Compute totals across all environments
total_mutant_clinical <- sum(env_wide_clinical$Mutant)
total_wt_clinical <- sum(env_wide_clinical$WT)

# Compute odds ratio for each environment
odds_ratio_all_clinical <- env_wide_clinical %>%
  mutate(
    Em = Mutant,
    Ew = WT,
    NEm = total_mutant_clinical - Mutant,
    NEw = total_wt_clinical - WT,
    m_ratio = Em / NEm,
    w_ratio = Ew / NEw,
    OR = ( Em / NEm ) / ( Ew / NEw )
  ) %>%
  select(Var2, Em, Ew, NEw, NEm, m_ratio, w_ratio,OR)

colnames(odds_ratio_all_clinical)[1] <- "Environment"

odds_ratio_all_clinical_file <- paste0(base, "_odds_ratio_clinical.csv")
write.csv(odds_ratio_all_clinical, odds_ratio_all_clinical_file, row.names = FALSE)

or_gt_clinical <- odds_ratio_all_clinical %>%
  gt()

or_gt_clinical_file <- paste0(base, "_odds_ratio_clinical.html")
gtsave(or_gt_clinical, or_gt_clinical_file)

# # ----------------------------
# # Tree figures for lineages
# # ----------------------------

# # Identify terminal edges
# tip_indices <- 1:Ntip(relaxed_0)  # Tips are numbered from 1 to Ntip
# terminal_edges <- which(relaxed_0$edge[, 2] %in% tip_indices)  # Find terminal edges

# # Create a dataframe that associates tip indices with terminal edges
# tip_to_edge_df <- data.frame(
#   tip_index = relaxed_0$edge[terminal_edges, 2],  # Tip index (descendant node)
#   edge_index = terminal_edges  # Corresponding edge index
# )

# # Associate tip labels with their terminal edges
# tip_to_edge_df$Genome <- relaxed_0$tip.label[tip_to_edge_df$tip_index]

# # Merge terminal edges with lineage status
# edge_merge <- merge(x = tip_to_edge_df, y = cheat_lineages$tips, by = "Genome", all = TRUE)

# # Fill missing lineage statuses as 'Non-lineage'
# edge_merge$Lineage_Status[is.na(edge_merge$Lineage_Status)] <- 'Non-lineage'

# # Subset table
# edge_merge <- edge_merge[c("edge_index", "Lineage_Status")]

# # Create a color map for the edges
# edge_map <- c('Lineage' = 'red', 'Non-lineage' = 'black')

# # Associate colors with lineage status, set default color to black
# edge_colors <- rep("black", nrow(relaxed_0$edge))

# # Match edge indices from the edge_merge dataframe
# edge_colors[edge_merge$edge_index] <- edge_map[as.character(edge_merge$Lineage_Status)]

# # Plot along with tip colors
# lin_tree_name_file <- paste0(base, "_new_lin_tree.svg")

# gene_label_wt <- paste0(gene, " WT")
# gene_label_mut <- paste0(gene, " Mutant")

# gene_label_non_lin <- paste0(gene, " Non-Lineage")
# gene_label_lin <- paste0(gene, " Lineage")

# svg(filename = lin_tree_name_file, width = 10, height = 10)

# plot.new()
# plot(relaxed_0,  edge.color = edge_colors, type = "fan", lwd = 0.5, show.tip.label = FALSE)
# tiplabels(pch = 19, col = tip_colors, cex = 0.5)
# par(family = "Times New Roman")
# legend("bottomright", inset = c(0.045, 0.07), legend = c(gene_label_wt, gene_label_mut),
#         col = c("blue", "red"), pch = 19, box.lty = 0, text.font = 1)
# legend("bottomright", legend = c(gene_label_non_lin, gene_label_lin),
#         col = c("black", "red"), lty = 1, box.lty = 0, text.font = 1)
# add.scale.bar(length = 0.1, x = -1, y = -1)

# dev.off ()

# # ----------------------------
# # Set lineages to ANI cutoff
# # ----------------------------

# ANI_clusters <- read.csv(no_funct_clusters_env)

# ANI_clusters <- ANI_clusters %>% filter(!Genome %in% zero_length_tips)

# ANI_clusters$Lineage_Status <- 'Lineage'

# # Merge terminal edges with lineage status
# edge_merge_ANI <- merge(x = tip_to_edge_df, y = ANI_clusters, by = "Genome", all = TRUE)

# print(edge_merge_ANI)

# # Fill missing lineage statuses as 'Non-lineage'
# edge_merge_ANI$Lineage_Status[is.na(edge_merge_ANI$Lineage_Status)] <- 'Non-lineage'

# # Subset table
# edge_merge_ANI <- edge_merge_ANI[c("edge_index", "Lineage_Status")]

# # Create a color map for the edges
# edge_map_ANI <- c('Lineage' = 'red', 'Non-lineage' = 'black')

# # Associate colors with lineage status, set default color to black
# edge_colors_ANI <- rep("black", nrow(relaxed_0$edge))

# # Match edge indices from the edge_merge dataframe
# edge_colors_ANI[edge_merge_ANI$edge_index] <- edge_map_ANI[as.character(edge_merge_ANI$Lineage_Status)]

# # Plot along with tip colors
# lin_tree_name_ANI_file <- paste0(base, "_new_lin_tree_ANI.svg")
# svg(filename = lin_tree_name_ANI_file, width = 10, height = 10)

# plot.new()
# plot(relaxed_0,  edge.color = edge_colors_ANI, type = "fan", lwd = 0.5, show.tip.label = FALSE)
# tiplabels(pch = 19, col = tip_colors, cex = 0.5)
# par(family = "Times New Roman")
# legend("bottomright", inset = c(0.045, 0.07), legend = c(gene_label_wt, gene_label_mut),
#         col = c("blue", "red"), pch = 19, box.lty = 0, text.font = 1)
# legend("bottomright", legend = c(gene_label_non_lin, gene_label_lin),
#         col = c("black", "red"), lty = 1, box.lty = 0, text.font = 1)
# add.scale.bar(length = 0.1, x = -1, y = -1)

# dev.off ()

# # Plot just lineages, without tips colors
# lin_tree_name_no_file <- paste0(base, "_new_lin_tree_no.svg")
# svg(filename = lin_tree_name_no_file, width = 10, height = 10)

# plot.new()
# plot(relaxed_0,  edge.color = edge_colors, type = "fan", lwd = 0.5, show.tip.label = FALSE)
# par(family = "Times New Roman")
# legend("bottomright", legend=c(gene_label_non_lin, gene_label_lin),
#         col=c("black", "red"), lty = 1, box.lty=0, text.font=1)
# add.scale.bar(length = 0.1, x = -1, y = -1)

# dev.off ()

# # Plot lineages with lineage tip colors
# lin_tree_name_tc_file <- paste0(base, "_new_lin_tree_tip_col.svg")

# gene_label_non_lin <- paste0(gene, " Non-Lineage")
# gene_label_lin <- paste0(gene, " Lineage")

# svg(filename = lin_tree_name_tc_file, width = 10, height = 10)

# plot.new()
# plot(relaxed_0,  edge.color = edge_colors, type = "fan", lwd = 0.5, show.tip.label = FALSE)
# tiplabels(pch = 19, col = lineage_tip_colors, cex = 0.5)
# par(family = "Times New Roman")
# legend("bottomright", legend=c(gene_label_non_lin, gene_label_lin),
#         col=c("black", "red"), lty = 1, box.lty=0, text.font=1)
# legend("right", legend = names(cluster_colors), 
#         col = cluster_colors, pch = 19, box.lty = 0, text.font = 1, cex = 0.7,
#         title = "Cluster")
# add.scale.bar(length = 0.1, x = -1, y = -1)

# dev.off ()
  
# # Plot lineages with environments
# # Read in metadata
# meta <- read.csv(genome_annotations_classified)

# funct <- read.csv(functions_incomplete)

# funct <- funct %>% filter(!Genome %in% zero_length_tips)

# # Merge metadata and functions
# meta <- merge(meta,funct, by = "Genome", all = TRUE)

# # Only keep mutants
# meta <- meta[grepl("No function", meta$Mut_Status), ]

# # Subset table
# meta <- meta[c("Genome","Group")]

# # Rename envs
# meta$Group <- gsub("Clinical, Unknown", "Clinical U", meta$Group)
# meta$Group <- gsub("Human-assocaited environmental", "Human Assoc", meta$Group)

# # Remove 0 length tips
# meta <- meta %>% filter(!Genome %in% zero_length_tips)

# # Transform metadata
# rownames(meta) <- meta[,1]
# meta[,1] <- NULL
# LasR_meta<-setNames(meta[,1],rownames(meta))

# #Set metadata colors
# color_map_env <- c('Clinical' = '#8B0000', 'Clinical U' = '#E74E00', 'CF' = '#DEAE21', 'Environmental' = '#3BB497', 'Human Assoc' = '#5DA2A7', 'Animal' = '#132157', 'Lab' = 'lightgray', 'Missing' = '#6C728C')

# #Map colors to tips
# tip_labels_env <- relaxed_0$tip.label

# tip_states_env <- LasR_meta[tip_labels]

# tip_colors_env <- color_map_env[as.character(tip_states_env)]

# env_tree_name_file <- paste0(base, "_env_tree.svg")
# svg(filename = env_tree_name_file, width = 10, height = 10)

# plot.new()
# plot(relaxed_0, edge.color = edge_colors, type = "fan", lwd = 0.5, show.tip.label = FALSE)
# tiplabels(pch = 19, col = tip_colors_env, cex = 0.5)
# par(family = "Times New Roman")
# legend("bottomright", inset = c(-0.065, 0.07), legend=c("Clinical", "Clinical U", "CF", "Environmental", "Human Assoc", "Animal", "Lab", "Missing"),
#         col=c("#8B0000", "#E74E00", "#DEAE21", "#3BB497", "#5DA2A7", "#132157", "lightgray", "#6C728C"), pch=19, box.lty=0, text.font=1, cex = 0.7)
# legend("bottomright", legend=c(gene_label_non_lin, gene_label_lin),
#         col=c("black", "red"), lty = 1, box.lty=0, text.font=1)
# add.scale.bar(length = 0.1, x = -1, y = -1)

# dev.off ()

# # ----------------------------
# # Chi-squared test for environmental association
# # ----------------------------

# # Read in environmental data
# env_data <- read.csv(genome_annotations_classified)

# # Merge with presence absence table (already read in)
# trait_table <- merge(env_data,tree_data, by = "Genome")

# # Subset table
# state_env <- trait_table[, c("Mut_Status", "Group")]

# # Get basic counts overview table, and save
# tbl <- table(state_env$Mut_Status, state_env$Group)

# tbl_file <- paste0(base, "_env.csv")
# write.csv(tbl, tbl_file)

# # Run a generalized linear model for the mutation status vs environment
# model <- glm(Mut_Status ~ Group, data = state_env, family = "binomial")

# # Obtain the output
# model_sum <- summary(model)

# coefficients <- model_sum$coefficients

# coefficients_df <- as.data.frame(coefficients)

# # If there are any significant associations, reformat the file
# if (any(coefficients_df$`Pr(>|z|)` < 0.05)) {
#   env_sig <- coefficients_df %>%
#     filter(`Pr(>|z|)` < 0.05) %>%
#     mutate(
#       Sig_env = paste(row.names(.), "in", ifelse(Estimate > 0, "No function", "Functional")),
#       Gene = paste(base)
#     ) %>%
#     dplyr::select(Gene, Sig_env, `Pr(>|z|)`, Estimate, `Std. Error`, `z value`)

#   env_sig$Sig_env <- as.character(env_sig$Sig_env)

#   env_sig$Sig_env <- gsub("Group", "", env_sig$Sig_env)

#   # Save
#   env_sig_file <- paste0(base, "_env_sig.csv")
#   write.csv(env_sig, env_sig_file, row.names=FALSE)

# } else {
#   env_sig <- NULL
# }

# # ----------------------------
# # Chi-squared test for lineage environmental association
# # ----------------------------

# # Only keep mutants on presence absence table
# tree_data_muts <- tree_data %>% filter(Mut_Status == 1)

# # Merge inner with lineage tips list to keep just mutant lineages
# muts_lineage_only <- merge(tree_data_muts,cheat_lineages$tips, by = "Genome", all = FALSE)

# muts_lineage_only_sub <- muts_lineage_only[c("Genome", "Node", "Lineage_Status")]

# # Merge outer with lineage tips list to keep just mutant lineages
# muts_lineage <- merge(tree_data_muts,muts_lineage_only_sub, by = "Genome", all = TRUE)

# # Subset
# muts_lineage <- muts_lineage[c("Genome", "Lineage_Status")]

# # If NA, replace with 0 (non lineage)
# muts_lineage <- muts_lineage %>%
#   mutate(Lineage_Status = ifelse(is.na(Lineage_Status), 0, Lineage_Status))

# # Merge environmental data with lineage table
# trait_table_lin <- merge(env_data,muts_lineage, by = "Genome")

# # Subset
# state_env_lin <- trait_table_lin[, c("Lineage_Status", "Group")]

# # Replace lineage with 1
# state_env_lin$Lineage_Status <- gsub("Lineage", 1, state_env_lin$Lineage_Status)

# # Make column numberic
# state_env_lin$Lineage_Status <- as.numeric(state_env_lin$Lineage_Status)

# # Get basic counts overview table, and save
# tbl_lin <- table(state_env_lin$Lineage_Status, state_env_lin$Group)

# tbl_lin_file <- paste0(base, "_lin_env.csv")
# write.csv(tbl_lin, tbl_lin_file)

# # Run a generalized linear model for the lineage status vs environment
# model_lin <- glm(Lineage_Status ~ Group, data = state_env_lin, family = "binomial")

# # Obtain the output
# model_sum_lin <- summary(model_lin)

# coefficients_lin <- model_sum_lin$coefficients

# coefficients_df_lin <- as.data.frame(coefficients_lin)

# # If there are any significant associations, reformat the file
# if (any(coefficients_df_lin$`Pr(>|z|)` < 0.05)) {
#   lin_env_sig <- coefficients_df_lin %>%
#     filter(`Pr(>|z|)` < 0.05) %>%  # Filter rows where A is less than 0.05
#     mutate(
#       Sig_env = paste(row.names(.), "in", ifelse(Estimate > 0, "lineage", "non-lineage")),
#       Gene = paste(base)
#     ) %>%
#     dplyr::select(Gene, Sig_env, `Pr(>|z|)`, Estimate, `Std. Error`, `z value`)

#   lin_env_sig$Sig_env <- as.character(lin_env_sig$Sig_env)

#   lin_env_sig$Sig_env <- gsub("Group", "", lin_env_sig$Sig_env)

#   lin_env_sig_file <- paste0(base, "_lin_env_sig.csv")
#   write.csv(lin_env_sig, lin_env_sig_file, row.names=FALSE)

# } else {
#   lin_env_sig <- NULL
# }
