library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

# Read in gene list
gene_list <- read.csv("canonical_CF_gene_full_list.csv")

# Summarize, dropping those with a count less than 3
gene_list_summary <- gene_list %>%
  distinct(gene_of_interest, Paper) %>%
  group_by(gene_of_interest) %>%
  summarise(
    count = n(),
    papers = paste(unique(Paper), collapse = ",")
  ) %>%
  filter(count >= 3) %>%
  arrange(desc(count))

write.csv(gene_list_summary, "canonical_CF_gene_full_summary_new.csv", row.names = FALSE)

# Summarize, inclusing all genes
gene_list_full_counts <- gene_list %>%
  distinct(gene_of_interest, Paper) %>%
  group_by(gene_of_interest) %>%
  summarise(
    count = n(),  # total occurrences
    papers = paste(unique(Paper), collapse = ",")
  ) %>%
  arrange(desc(count))

write.csv(gene_list_full_counts, "canonical_CF_gene_full_summary_new_all_genes.csv", row.names = FALSE)


gene_list_full_counts$Status <- ifelse(gene_list_full_counts$count >= 3, "Kept", "Dropped")

status_colors <- c(
  "Kept"          = "#C36557",
  "Dropped"         = "#9CBED2"
)

gene_list_plot <- ggplot(gene_list_full_counts, aes(x = count, y = reorder(gene_of_interest, count))) +
  geom_bar(stat = "identity", aes(fill = Status), width = 0.5) +
  scale_fill_manual(values = status_colors, name = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  xlab("Number of publication mentions") +
  ylab("Gene")

gene_list_plot

ggsave("gene_list_plot.svg", plot = gene_list_plot, width = 7.5, height = 15, units = "in")


# Visualization for dropped list
gene_list_dropped_plot <- ggplot(gene_list_summary, aes(x = count, y = reorder(gene_of_interest, count))) +
  geom_bar(stat = "identity", fill = "#C36557", width = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  xlab("Number of publication mentions") +
  ylab("Gene")

ggsave("gene_list_dropped_plot.svg", plot = gene_list_dropped_plot, width = 5, height = 7.5, units = "in")

# Re-grouping genes
gene_list_summary_cat <- gene_list_summary %>%
  mutate(Function = case_when(
    gene_of_interest == "lasR" ~ "Quorum sensing",
    gene_of_interest %in% c("mucA", "algU", "algG", "pelA") ~ "Biofilm",
    gene_of_interest %in% c("mutS", "mutL") ~ "DNA repair",
    gene_of_interest %in% c("mexZ", "mexY", "mexX", "mexA", "mexB",
                            "gyrA", "gyrB", "oprD", "nfxB", "ampC") ~ "Antibiotic resistance",
    gene_of_interest == "rpoN" ~ "Motility",
    gene_of_interest %in% c("aceE", "aceF") ~ "Metabolism",
    TRUE ~ NA_character_
  ))

group_colors <- c(
  "Quorum sensing"         = "#9CBED2",
  "Biofilm"                = "#9CA780",
  "DNA repair"             = "#F3BA36",
  "Motility"               = "#F28D32",
  "Metabolism"             = "#E36B2E",
  "Antibiotic resistance"  = "#C36557"
)


gene_list_dropped_color_plot <- ggplot(gene_list_summary_cat, aes(x = count, y = reorder(gene_of_interest, count))) +
  geom_bar(stat = "identity", aes(fill = Function), width = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  scale_fill_manual(values = group_colors) +
  xlab("Number of publication mentions") +
  ylab("Gene")

ggsave("gene_list_dropped_color_plot.svg", plot = gene_list_dropped_color_plot, width = 6.5, height = 7.5, units = "in")

# Plotting the genes as groups
gene_list_summary_cat_plot <- gene_list_summary_cat %>%
  arrange(Function, count) %>%
  mutate(gene_of_interest = factor(gene_of_interest, levels = unique(gene_of_interest)))

totals <- gene_list_summary_cat_plot %>%
  group_by(Function) %>%
  summarize(total_count = sum(count), .groups = "drop")

gene_list_summary_cat_plot <- gene_list_summary_cat_plot %>%
  left_join(totals, by = "Function") %>%
  mutate(Function = reorder(Function, total_count))


# ----------------------------
# Figure 1: Plot of genes from publication review
# ----------------------------

gene_list_dropped_grouped_plot <- ggplot(gene_list_summary_cat_plot, aes(x = count, y = Function, fill = Function, group = gene_of_interest)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, linewidth = 0.5, color = "black") +
  scale_fill_manual(values = group_colors) +
  guides(fill = "none") +
  geom_text(aes(label = gene_of_interest),
            position = position_stack(vjust = 0.5),
            size = 4, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  xlab("Number of publication mentions") +
  ylab("Gene function")

ggsave("gene_list_dropped_grouped_plot.svg", plot = gene_list_dropped_grouped_plot, width = 10, height = 7.5, units = "in")
