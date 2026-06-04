library(ggplot2)
library(tidyr)
library(emmeans)
library(Ternary)
library(dplyr)
library(forcats)
library(tidytext)
library(ggrepel)
library(stringr)
library(nnet)
library(rcompanion)
library(vegan)
library(umap)
library(infotheo)
library(purrr)
library(ggh4x)
library(patchwork)
library(viridis)
library(tibble)
library(pheatmap)
library(grid)
library(factoextra)

# ----------------------------
# Initial data management
# ----------------------------

# Read in mapping file for genes
status <- read.csv("canonical_CF_gene_full_summary_new.csv")
status <- status[, c("gene_of_interest", "PA_locus_tag")]
colnames(status)[1] <- "Gene_funct"

colnames(status)[2] <- "Gene"

# Define the base directory
list_dir <- "E:/IFI/Brown_lab_research/Legacy_Data/Canonical_CF_genes/PAO1_new/outputs_for_graphing"

# Get a list of all directories in the base directory (excluding the base directory itself)
list_dirs <- list.dirs(path = list_dir, full.names = TRUE, recursive = FALSE)

run <- "new"

if (run == "first") {
  # IF RUNNING FOR THE FIRST TIME
  # Read in environmental data
  env_data <- read.csv("allGenomesAllNichesMetadata_BronchMapped.csv")

  # Remove clade C from data
  env_data <- env_data %>% 
    filter(!str_detect(Clade, "CladeC")) 

  # Reformat genome column
  colnames(env_data)[1] <- "Genome"

  env_data$Genome <- gsub("\\.", "_", env_data$Genome)

  # Map values
  map_vals <- c(
    'Bronchiectasis'= 'Bronchiectasis',
    'Pneumonia'= 'Pneumonia',
    'CF'= 'CF',
    'early.CF'= 'Early CF',
    'Urinary'= 'Human Infection',
    'Gastrointestinal'= 'Human Infection', # Not used
    'Wound'= 'Human Infection',
    'Burn'= 'Human Infection',
    'Blood'= 'Human Infection',
    'Eye'= 'Human Infection',
    'Terrestrial'= 'Environmental',
    'Aquatic'= 'Environmental',
    'Human_Environment'= 'Human Environment',
    'Animal'= 'Animal',
    'Ocean'= 'Environmental',
    'Built.Environment'= 'Human Environment',
    'Hospital.Environment'= 'Human Environment',
    'Rectal/Feces'= 'Human Infection',
    'Waste.Water'= 'Human Environment'
  )

  # replace values in niche
  env_data <- env_data %>%
    mutate(Niche_mapped = recode(Niche, !!!map_vals, .default = NA_character_)) %>%
    filter(!is.na(Niche_mapped))

  # Remove animal from environmental data
  env_data <- env_data %>% 
    filter(!str_detect(Niche, "Animal")) 

  # Drop clones
  env_data <- env_data %>%
    distinct(Clone_Cluster, Niche_mapped, .keep_all = TRUE)

  write.csv(env_data, "Environmental_data_dropped_animal_clones.csv", row.names = FALSE)
  
} else {
  
  env_data <- read.csv("Environmental_data_dropped_animal_clones.csv")
  
}

# ----------------------------
# Generate environmental counts tables
# ----------------------------

for (dir in list_dirs) {
  base <- basename(dir)
  
  print(base)
  
  # Read in presence absence table
  funct_data_file <- paste0(base, "_funct_no_funct_alt.csv")
  funct_data_path <- paste0(dir, "/", funct_data_file)
  funct_data <- read.csv(funct_data_path)
  
  # Merge with presence absence table (already read in)
  trait_table <- merge(env_data, funct_data, by = "Genome", all.x = TRUE)
  
  # Drop contig splits
  trait_table <- trait_table %>%
    filter(!grepl("Contig split", Mutation))
  
  state_env <- trait_table[, c("Mut_Status", "Niche_mapped")]
  
  # Get basic counts overview table, and save
  tbl <- table(state_env$Mut_Status, state_env$Niche_mapped)
  
  tbl_file <- paste0(base, "_env.csv")
  tbl_file_path <- paste0(dir, "/", tbl_file)
  write.csv(tbl, tbl_file_path)
  
}

# ----------------------------
# Looping through all gene output directories for overall counts
# ----------------------------

# Initializing lists
tbl_list <- list()

for (dir in list_dirs) {
  base <- basename(dir)
  
  print(base)
  
  # Read in environmental table counts
  tbl_file <- paste0(base, "_env.csv")
  tbl_file_path <- paste0(dir, "/", tbl_file)
  tbl_table <- read.csv(tbl_file_path)
  
  colnames(tbl_table)[1] <- "Mut_Status"
  
  # Make columns if missing
  required_cols <- c("Pneumonia", "Bronchiectasis", "CF", "Early.CF", 
                     "Human.Infection", "Environmental", "Human.Environment")
  
  # Add missing columns with value 0
  for (col in required_cols) {
    if (!col %in% colnames(tbl_table)) {
      tbl_table[[col]] <- 0
    }
  }
  
  required_status <- c("Functional", "No function", "Alternative function")
  
  tbl_table <- tbl_table %>%
    tidyr::complete(
      Mut_Status = required_status,
      fill = list(
        Bronchiectasis = 0,
        CF = 0,
        Early.CF = 0,
        Environmental = 0,
        Human.Environment = 0,
        Human.Infection = 0,
        Pneumonia = 0
      )
    )
  
  tbl_table$Gene <- base
  
  print(tbl_table)
  
  tbl_list[[length(tbl_list) + 1]] <- tbl_table
}

# Concatenate environmental association tables
combined_tbl <- do.call(rbind, tbl_list)

# Merge with social status file
tbl_status <- left_join(combined_tbl, status, by = "Gene")

# ----------------------------
# Figure 2: Graphing environments by clade distribution
# ----------------------------
env_data_counts_clade <- env_data %>%
  group_by(Clade, Niche_mapped) %>%
  summarise(env_n = n())

clade_counts <- env_data %>%
  group_by(Clade) %>%
  summarise(clade_n = n())

merge_env_clade_counts <- left_join(env_data_counts_clade, clade_counts, by = "Clade")

merge_env_clade_counts$Proportion <- merge_env_clade_counts$env_n / merge_env_clade_counts$clade_n

merge_env_clade_counts$Clade <- gsub("CladeA", "Clade A", merge_env_clade_counts$Clade)
merge_env_clade_counts$Clade <- gsub("CladeB", "Clade B", merge_env_clade_counts$Clade)


niche_colors <- c(
  "Human Environment"  = "#9CBED2",
  "Environmental"      = "#9CA780",
  "CF"                 = "#F3BA36",
  "Early CF"           = "#F7D074",
  "Bronchiectasis"     = "#E36B2E",
  "Pneumonia"          = "#F28D32",
  "Human Infection"    = "#C36557"
)

# Figure 2 panel A
clade_counts_raw <- ggplot(merge_env_clade_counts, aes(x = Clade, y = env_n, fill = Niche_mapped)) +
  geom_bar(stat = "identity", color = "black", aes(fill = Niche_mapped)) +
  xlab("Clade") +
  ylab("Count") +
  scale_fill_manual(values = niche_colors, name = "Environment") +
  theme_bw()

ggsave("Figures_4_21_2026/Clade_counts_no_lung_no_animal_no_clone.svg", plot = clade_counts_raw, width = 7.5, height = 7.5, units = "in")
ggsave("Figures_4_21_2026/Clade_counts_no_lung_no_animal_no_clone.png", plot = clade_counts_raw, width = 7.5, height = 7.5, units = "in")

# Figure 2 panel B
clade_proportion <- ggplot(merge_env_clade_counts, aes(x = Clade, y = Proportion, fill = Niche_mapped)) +
  geom_bar(stat = "identity", color = "black", aes(fill = Niche_mapped)) +
  labs(fill = "Niche") +
  xlab("Clade") +
  geom_text(aes(label = env_n),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  scale_fill_manual(values = niche_colors, name = "Environment") +
  theme_bw()

ggsave("Figures_4_21_2026/Clade_proportion_no_lung_no_animal_no_clone.svg", plot = clade_proportion, width = 7.5, height = 7.5, units = "in")
ggsave("Figures_4_21_2026/Clade_proportion_no_lung_no_animal_no_clone.png", plot = clade_proportion, width = 7.5, height = 7.5, units = "in")

# ----------------------------
# Figure 3: Heatmap of loss-of-function and alternate function mutants per gene
# ----------------------------

diff_from_mean <- function() {
  # Subset all environments
  gene_tbl_social <- tbl_status[c("Mut_Status", "Bronchiectasis", "CF", "Early.CF", "Pneumonia", 
                                  "Human.Infection", "Environmental", "Human.Environment", "Gene_funct")]
  
  # Transform table so that status becomes the new columns
  gene_tbl_social_t <- gene_tbl_social %>%
    pivot_longer(cols = !c(Mut_Status, Gene_funct), names_to = "Environments", values_to = "Value") %>%
    pivot_wider(names_from = Mut_Status, values_from = Value)
  
  gene_tbl_social_t <- as.data.frame(gene_tbl_social_t)
  
  # Calculate totals and proportions
  gene_props <- gene_tbl_social_t %>%
    group_by(Gene_funct) %>%
    summarise(
      total_alt  = sum(`Alternative function`, na.rm = TRUE),
      total_func = sum(Functional, na.rm = TRUE),
      total_no   = sum(`No function`, na.rm = TRUE),
      total_all  = total_alt + total_func + total_no,
      
      prop_alt  = (total_alt  / total_all) * 100,
      prop_func = (total_func / total_all) * 100,
      prop_no   = (total_no   / total_all) * 100,
      prop_both = ((total_no + total_alt)   / total_all) * 100,
      
      .groups = "drop"
    ) 
  
  # Calculate percentages and reformat
  gene_tbl_social_t <- gene_tbl_social_t %>%
    left_join(gene_props, by = "Gene_funct") %>%
    mutate(Sum = Functional + `No function` + `Alternative function`) %>%
    mutate(`Percent functional` = (Functional / Sum) * 100) %>%
    mutate(`Percent non-functional` = (`No function` / Sum) * 100) %>%
    mutate(`Percent potential alt function` = (`Alternative function` / Sum) * 100) %>%
    mutate(`Percent non and potential alt function` = ((`Alternative function` + `No function`) / Sum) * 100) %>%
    mutate(`Difference functional` = `Percent functional` - prop_func) %>%
    mutate(`Difference non-functional` = `Percent non-functional` - prop_no) %>%
    mutate(`Difference potential alt function` = `Percent potential alt function` - prop_alt) %>%
    mutate(`Difference non and potential alt function` = `Percent non and potential alt function` - prop_both) %>%
    mutate(
      Environments = str_replace_all(Environments, c(
        "Human\\.Environment" = "Human Environment",
        "Early\\.CF"          = "Early CF",
        "Human\\.Infection"   = "Human Infection"
      ))
    )
  
  # Re-order environements
  gene_tbl_social_t$Environments <- factor(gene_tbl_social_t$Environments, levels = c("CF", "Early CF", "Bronchiectasis",
                                                                                      "Pneumonia", "Human Infection",
                                                                                      "Human Environment", "Environmental"))
  
  # Re-order genes
  gene_tbl_social_t$Gene_funct <- factor(gene_tbl_social_t$Gene_funct, levels = rev(c("lasR", "mucA", "algU", "mutL", "mutS",
                                                                                      "rpoN", "gyrA", "gyrB", "mexY", "mexZ",
                                                                                      "aceE", "aceF", "algG", "ampC", "mexA",
                                                                                      "mexB", "mexX", "nfxB", "oprD", "pelA")))
  
  ## Loss-of-function results ##
  # Supplemental figure 2 panel A
  plot_no_funct <- ggplot(gene_tbl_social_t, aes(x = Environments, y = Gene_funct, fill = `Difference non-functional` )) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ylab("Gene") +
    xlab("Environment") +
    labs(
      fill = str_wrap("Difference non-functional", width = 20)
    ) +
    theme_minimal() +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_diff_no_clones.png", plot = plot_no_funct, width = 5, height = 7.5, units = "in")
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_diff_no_clones.svg", plot = plot_no_funct, width = 5, height = 7.5, units = "in")
  
  
  # Create df matrix precursor from df
  mat_df_no_funct <- gene_tbl_social_t %>%
    select(Environments, Gene_funct, `Difference non-functional`) %>%
    pivot_wider(names_from = Environments, values_from = `Difference non-functional`)

  # Create matrix
  mat_no_funct <- mat_df_no_funct %>%
    column_to_rownames("Gene_funct") %>%
    as.matrix()
  
  # Calculate maximum value
  max_val_no_funct <- max(abs(mat_no_funct), na.rm = TRUE)

  # Calculate breaks based on max value
  breaks_no_funct <- seq(-max_val_no_funct, max_val_no_funct, length.out = 101)

  set.seed(123)

  dev.off()

  svg("Figures_4_21_2026/Canonical_CF_genes_full_pheatmap_clustered_no_funct_diff_no_clones.svg", width = 5, height = 7.5)

  plot.new()

  # Figure 3
  pheatmap(
    mat_no_funct,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = breaks_no_funct
  )

  dev.off()

  svg(filename = "Figures_4_21_2026/Canonical_CF_genes_full_pheatmap_k2_no_funct_diff_no_clones.svg", width = 7.5, height = 10)
  
  transposed_mat_no_funct <- t(mat_no_funct)
  
  # Calculate kmeans gap start
  gap_start_no_funct <- fviz_nbclust(
    transposed_mat_no_funct,
    FUNcluster = kmeans, 
    k.max = 6,
    method = "gap_stat",
    nboot = 500,
    diss = get_dist(transposed_mat_no_funct, method = "euclidean")
  )
  
  ggsave("Figures_4_21_2026/gap_start_no_funct.svg", plot = gap_start_no_funct, width = 7.5, height = 7.5, units = "in")

  # Calculate kmeans silhouette
  silhouette_no_funct <- fviz_nbclust(
    transposed_mat_no_funct,
    FUNcluster = kmeans, 
    k.max = 6,
    method = "silhouette",
    nboot = 500,
    diss = get_dist(transposed_mat_no_funct, method = "euclidean")
  )
  
  ggsave("Figures_4_21_2026/silhouette_no_funct.svg", plot = silhouette_no_funct, width = 7.5, height = 7.5, units = "in")
  
  # Calculate kmeans WSS
  wss_no_funct <- fviz_nbclust(
    transposed_mat_no_funct,
    FUNcluster = kmeans, 
    k.max = 6,
    method = "wss",
    nboot = 500,
    diss = get_dist(transposed_mat_no_funct, method = "euclidean")
  )
  
  ggsave("Figures_4_21_2026/wss_no_funct.svg", plot = wss_no_funct, width = 7.5, height = 7.5, units = "in")
  
  # Get all results of all Ks for k-means clustering
  k_results_list_no_funct <- lapply(2:6, function(k) {

    ph <- pheatmap(transposed_mat_no_funct, kmeans_k = k)

    data.frame(
      K = k,
      sample = names(ph$kmeans$cluster),
      cluster = ph$kmeans$cluster
    )
  })

  k_results_list_no_funct_df <- bind_rows(k_results_list_no_funct)
  
  ## Combined loss-of-function and alternate function ##
  # Supplemental figure 2 panel B
  plot_both <- ggplot(gene_tbl_social_t, aes(x = Environments, y = Gene_funct, fill = `Difference non and potential alt function` )) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ylab("Gene") +
    xlab("Environment") +
    labs(
      fill = str_wrap("Difference non and potential alt function", width = 20)
    ) +
    theme_minimal() +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_alt_funct_diff_no_clones.png", plot = plot_both, width = 5, height = 7.5, units = "in")
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_alt_funct_diff_no_clones.svg", plot = plot_both, width = 5, height = 7.5, units = "in")
  
  # Create df matrix precursor from df
  mat_df_both <- gene_tbl_social_t %>%
    select(Environments, Gene_funct, `Difference non and potential alt function`) %>%
    pivot_wider(names_from = Environments, values_from = `Difference non and potential alt function`)
  
  # Create matrix
  mat_both <- mat_df_both %>%
    column_to_rownames("Gene_funct") %>%
    as.matrix()
  
  # Calculate maximum value
  max_val_both <- max(abs(mat_both), na.rm = TRUE)
  
  # Calculate breaks based on max value
  breaks_both <- seq(-max_val_both, max_val_both, length.out = 101)
  
  svg(filename = "Figures_4_21_2026/Canonical_CF_genes_full_pheatmap_clustered_both_diff_no_clones.svg", width = 5, height = 7.5)
  
  plot.new()
  
  # Supplemental figure x
  pheatmap(
    mat_both,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = breaks_no_funct
  )
  
  dev.off ()
  
  transposed_mat_both <- t(mat_both)
  
  # Get all results of all Ks for k-means clustering
  k_results_list_both <- lapply(2:6, function(k) {
    
    ph <- pheatmap(transposed_mat_both, kmeans_k = k)
    
    data.frame(
      K = k,
      sample = names(ph$kmeans$cluster),
      cluster = ph$kmeans$cluster
    )
  })
  
  k_results_list_both_df <- bind_rows(k_results_list_both)
  
  return(list(no_funct = k_results_list_no_funct_df, both = k_results_list_both_df, mat_no_funct = mat_no_funct))
}

# Barplot environments
diff_from_mean <- diff_from_mean()

# ----------------------------
# Supplemental figure X: Heatmaps of overall counts
# ----------------------------

# Graph all WT with mutant and non-deleterious mutants per gene
all_bar_heatmap <- function() {
  # Subset all environments
  gene_tbl_social <- tbl_status[c("Mut_Status", "Bronchiectasis", "CF", "Early.CF", "Pneumonia", 
                                  "Human.Infection", "Environmental", "Human.Environment", "Gene_funct")]
  
  # Transform table so that status becomes the new columns
  gene_tbl_social_t <- gene_tbl_social %>%
    pivot_longer(cols = !c(Mut_Status, Gene_funct), names_to = "Environments", values_to = "Value") %>%
    pivot_wider(names_from = Mut_Status, values_from = Value)
  
  gene_tbl_social_t <- as.data.frame(gene_tbl_social_t)
  
  # Calculate percentages and reformat
  gene_tbl_social_t <- gene_tbl_social_t %>%
    mutate(Sum = Functional + `No function` + `Alternative function`) %>%
    mutate(`Percent functional` = (Functional / Sum) * 100) %>%
    mutate(`Percent non-functional` = (`No function` / Sum) * 100) %>%
    mutate(`Percent potential alt function` = (`Alternative function` / Sum) * 100) %>%
    mutate(`Percent non and potential alt function` = ((`Alternative function` + `No function`) / Sum) * 100) %>%
    mutate(
      Environments = str_replace_all(Environments, c(
        "Human\\.Environment" = "Human Environment",
        "Early\\.CF"          = "Early CF",
        "Human\\.Infection"   = "Human Infection"
      ))
    )
  
  gene_tbl_social_t$Environments <- factor(gene_tbl_social_t$Environments, levels = c("CF", "Early CF", "Bronchiectasis",
                                                                                      "Pneumonia", "Human Infection",
                                                                                      "Human Environment", "Environmental"))
  
  gene_tbl_social_t$Gene_funct <- factor(gene_tbl_social_t$Gene_funct, levels = rev(c("lasR", "mucA", "algU", "mutL", "mutS",
                                                                                      "rpoN", "gyrA", "gyrB", "mexY", "mexZ",
                                                                                      "aceE", "aceF", "algG", "ampC", "mexA",
                                                                                      "mexB", "mexX", "nfxB", "oprD", "pelA")))
  
  
  plot_no_funct <- ggplot(gene_tbl_social_t, aes(x = Environments, y = Gene_funct, fill = `Percent non-functional` )) +
    geom_tile(color = "white") +
    # scale_fill_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = 50, limits = c(0, 100)) +
    scale_fill_viridis_c(limits = c(0, 100)) +
    ylab("Gene") +
    xlab("Environment") +
    labs(
      fill = str_wrap("Percent non-functional", width = 20)
    ) +
    theme_minimal() +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_no_clones.png", plot = plot_no_funct, width = 7.5, height = 10, units = "in")
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_no_clones.svg", plot = plot_no_funct, width = 7.5, height = 10, units = "in")
  
  
  plot_both <- ggplot(gene_tbl_social_t, aes(x = Environments, y = Gene_funct, fill = `Percent non and potential alt function` )) +
    geom_tile(color = "white") +
    # scale_fill_gradient(low = "blue", high = "red", limits = c(0, 100)) +
    scale_fill_viridis_c(limits = c(0, 100)) +
    ylab("Gene") +
    xlab("Environment") +
    labs(
      fill = str_wrap("Percent non and potential alt function", width = 20)
    ) +
    theme_minimal() +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_alt_funct_no_clones.png", plot = plot_both, width = 7.5, height = 10, units = "in")
  ggsave("Figures_4_21_2026/Canonical_CF_genes_full_heatmap_no_funct_alt_funct_no_clones.svg", plot = plot_both, width = 7.5, height = 10, units = "in")
  
  genes_to_remove <- c("mexY", "aceF", "ampC", "mexX", "oprD", "pelA")
  
  gene_tbl_social_t_filtered <- gene_tbl_social_t %>%
    filter(!Gene_funct %in% genes_to_remove)
  
  plot_both_filtered <- ggplot(gene_tbl_social_t_filtered, aes(x = Environments, y = Gene_funct, fill = `Percent non and potential alt function` )) +
    geom_tile(color = "white") +
    # scale_fill_gradient(low = "blue", high = "red", limits = c(0, 100)) +
    scale_fill_viridis_c() +
    ylab("Gene") +
    xlab("Environment") +
    labs(
      fill = str_wrap("Percent non and potential alt function", width = 20)
    ) +
    theme_minimal() +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  return(plot_both_filtered)
}

# Barplot environments
all_bar_heatmap <- all_bar_heatmap()

# ----------------------------
# Supplemental figure X: Counts across genes and environments
# ----------------------------
all_bar <- function() {
  # Subset all environments
  gene_tbl_social <- tbl_status[c("Mut_Status", "Bronchiectasis", "CF", "Early.CF", "Pneumonia", 
                                  "Human.Infection", "Environmental", "Human.Environment", "Gene_funct")]
  
  # Transform table so that status becomes the new columns
  gene_tbl_social_t <- gene_tbl_social %>%
    pivot_longer(cols = !c(Mut_Status, Gene_funct), names_to = "Environments", values_to = "Value") %>%
    pivot_wider(names_from = Mut_Status, values_from = Value)

  gene_tbl_social_t <- as.data.frame(gene_tbl_social_t)

  # Calcuate a sum of WT, mutant, and nd mutant
  gene_tbl_social_t$Sum <- gene_tbl_social_t$Functional + gene_tbl_social_t$`No function` +  gene_tbl_social_t$`Alternative function`

  # Get the percentages of WT, mutant, and nd mutant
  gene_tbl_social_t$Percent_Functional <- gene_tbl_social_t$Functional/gene_tbl_social_t$Sum
  gene_tbl_social_t$Percent_No_function <- gene_tbl_social_t$`No function`/gene_tbl_social_t$Sum
  gene_tbl_social_t$Percent_alt_mutant <- gene_tbl_social_t$`Alternative function`/gene_tbl_social_t$Sum
  
  # Subset just the counts
  gene_tbl_social_envs_count <- gene_tbl_social_t[c("Gene_funct", "Environments", "Functional", "No function", "Alternative function")]

  # Transform the dataframe by environmets
  gene_tbl_social_envs_count_t <- gene_tbl_social_envs_count %>%
    pivot_longer(cols = !c(Gene_funct, Environments), names_to = "Status", values_to = "Count")

  gene_tbl_social_envs_count_t <- as.data.frame(gene_tbl_social_envs_count_t)

  gene_tbl_social_envs_count_t$Status <- gsub('Alternative function', 'alt_mutant', gene_tbl_social_envs_count_t$Status)
  gene_tbl_social_envs_count_t$Status <- gsub('No function', 'No_function', gene_tbl_social_envs_count_t$Status)

  gene_tbl_social_envs_count_t <- gene_tbl_social_envs_count_t %>% mutate(Status = paste0("Percent_", Status))
  
  # Subset just the percentages
  gene_tbl_social_envs <- gene_tbl_social_t[c("Gene_funct", "Environments", "Percent_Functional", "Percent_No_function", "Percent_alt_mutant")]

  # Transform the dataframe by environmets
  gene_tbl_social_envs_t <- gene_tbl_social_envs %>%
    pivot_longer(cols = !c(Gene_funct, Environments), names_to = "Status", values_to = "Percent")

  gene_tbl_social_envs_t <- as.data.frame(gene_tbl_social_envs_t)
  
  gene_tbl_social_envs_t <-  gene_tbl_social_envs_t %>% inner_join(gene_tbl_social_envs_count_t, by = c("Gene_funct", "Environments", "Status"))

  merged_gene_tbl_social_envs_t <- gene_tbl_social_envs_t %>%
    mutate(Count = ifelse(Count == 0, "", Count))

  # Rename percentages
  merged_gene_tbl_social_envs_t$Status <- gsub("Percent_Functional", "Percent functional", merged_gene_tbl_social_envs_t$Status)
  merged_gene_tbl_social_envs_t$Status <- gsub("Percent_alt_mutant", "Percent potential alt function", merged_gene_tbl_social_envs_t$Status)
  merged_gene_tbl_social_envs_t$Status <- gsub("Percent_No_function", "Percent non-functional", merged_gene_tbl_social_envs_t$Status)
  
  # Rename environments
  merged_gene_tbl_social_envs_t$Environments <- gsub("Human.Environment", "Human Environment", merged_gene_tbl_social_envs_t$Environments)
  merged_gene_tbl_social_envs_t$Environments <- gsub("Early.CF", "Early CF", merged_gene_tbl_social_envs_t$Environments)
  merged_gene_tbl_social_envs_t$Environments <- gsub("Human.Infection", "Human Infection", merged_gene_tbl_social_envs_t$Environments)

  
  # Re-order
  merged_gene_tbl_social_envs_t$Status <- factor(
    merged_gene_tbl_social_envs_t$Status,
    levels = c("Percent functional",
               "Percent potential alt function",
               "Percent non-functional")
  )
  

  # Merge with stat names
  merged_gene_tbl_social_envs_t_sig <- merge(merged_gene_tbl_social_envs_t, merge_CF_chronic_status, by = c("Gene_funct"), all = TRUE)

  # Re-order bars and facets according to disease type and lit review order
  merged_gene_tbl_social_envs_t_sig$Environments <- factor(merged_gene_tbl_social_envs_t_sig$Environments, levels = c("CF", "Early CF", "Bronchiectasis",
                                                                                                              "Pneumonia", "Human Infection",
                                                                                                              "Human Environment", "Environmental"))
  # Re order levels
  merged_gene_tbl_social_envs_t_sig <- merged_gene_tbl_social_envs_t_sig %>%
    mutate(
      Significance = factor(
        Significance,
        levels = c(
          "Adult CF significant",
          "Adult CF and chronic lung significant",
          "Chronic lung significant",
          "Non-significant"
        )
      )
    )
  
  # Plot
  plot <- ggplot(merged_gene_tbl_social_envs_t_sig, aes(x = Environments, y = Percent, fill = Status)) +
    geom_bar(position="stack", stat="identity", color = "black")  +
    geom_text(aes(label = Count), color = "white", size = 2, position = position_stack(vjust = 0.5)) +
    # geom_text(aes(label = significance), color = "white", size = 3, position = position_stack(vjust = 0.7)) +
    facet_nested_wrap(vars(Significance, Gene_funct)) +
    scale_fill_manual(values = c("Percent functional" = "blue", "Percent potential alt function" = "lightpink", "Percent non-functional" = "red", "Percent contig split" = "darkgray"), name = "Mutant Status") +
    xlab("Environments") +
    ylab("Percent") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme_bw()
  
  return(plot)
}

# Barplot environments
all_bar <- all_bar()
all_bar