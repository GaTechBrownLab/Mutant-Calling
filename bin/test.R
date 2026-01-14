
# ----------------------------
# Identify lineages based on ancestral state reconstruction
# ----------------------------

# Covert st to dataframe
st_df <- as.data.frame(st_ASR)

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
    st_df_lin <- st_df %>% filter(state_prob >= lineage_probability)

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

    return(node_tips_df_unique)
}

cheat_lineages <- identify_lineages(State1_Probability, "Mutant")
WT_lineages <- identify_lineages(State0_Probability, "WT")

# ----------------------------
# Associate lineages with colors for cheat graphing
# ----------------------------

# Subset table
lineage_clusters_sub <- cheat_lineages[c("Genome","Cluster")]

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
