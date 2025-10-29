#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)
library(svglite)
library(forcats)

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
base <- args[1]
env_counts <- args[2]
lineage_counts <- args[3]


# Read in environmental counts
environmental_table <- read.csv(env_counts)
colnames(environmental_table)[1] <- "Mut_Status"

# Read in lineage counts
lineage_table <- read.csv(lineage_counts)
colnames(lineage_table)[1] <- "Lineage_Status"

# Subset all environments
environmental_table_sub <- environmental_table[c("Mut_Status", "CF", "Clinical", "Clinical..Unknown", "Environmental",
                                      "Human.associated.environmental")]
# Subset just WT
WT_tbl <- environmental_table_sub %>%
  filter(grepl("0", Mut_Status))

WT_tbl$Status <- "WT"

# Drop Mut_Status column
WT_tbl <- WT_tbl[, !names(WT_tbl) %in% "Mut_Status"]
  
# Subset all environments
lineage_table_sub <- lineage_table[c("Lineage_Status", "CF", "Clinical", "Clinical..Unknown", "Environmental",
                                      "Human.associated.environmental")]
  
# Initialize status column
lineage_table_sub <- lineage_table_sub %>%
  mutate(Status = NA_character_)

# Label as lineage or non-lineage
lineage_table_sub <- lineage_table_sub %>%
  mutate(Status = ifelse(Lineage_Status == 0, "Non-lineage",
                                  ifelse(Lineage_Status == 1, "Lineage", Status)))

# Drop lineage_status column
lineage_table_sub <- lineage_table_sub[, !names(lineage_table_sub) %in% "Lineage_Status"]

# Concatenate WT table with lineage table
wt_lin_tbl <- rbind(WT_tbl, lineage_table_sub)
  
# Transform table so that status becomes the new columns
wt_lin_tbl_t <- wt_lin_tbl %>%
  pivot_longer(cols = -Status, names_to = "Environments", values_to = "Value") %>%
  pivot_wider(names_from =Status, values_from = Value)
  
wt_lin_tbl_t <- as.data.frame(wt_lin_tbl_t)

# Calcuate a sum of WT, lineage, and non-lineage
wt_lin_tbl_t$Sum <- wt_lin_tbl_t$WT + wt_lin_tbl_t$`Non-lineage` + wt_lin_tbl_t$Lineage
  
# Get the percentages of WT, lineage, and non-lineage
wt_lin_tbl_t$Percent_WT <- wt_lin_tbl_t$WT/wt_lin_tbl_t$Sum
wt_lin_tbl_t$`Percent_Non-lineage` <- wt_lin_tbl_t$`Non-lineage`/wt_lin_tbl_t$Sum
wt_lin_tbl_t$Percent_Lineage <- wt_lin_tbl_t$Lineage/wt_lin_tbl_t$Sum

# Subset just the counts
wt_lin_envs_count <- wt_lin_tbl_t[c("Environments", "WT", "Non-lineage", "Lineage")]

# Transform the dataframe by environmets
wt_lin_envs_count_t <- wt_lin_envs_count %>%
  pivot_longer(cols = -Environments, names_to = "Status", values_to = "Count")

wt_lin_envs_count_t <- as.data.frame(wt_lin_envs_count_t)

wt_lin_envs_count_t <- wt_lin_envs_count_t %>% mutate(Status = paste0("Percent_", Status))

# Subset just the percentages
wt_lin_envs <- wt_lin_tbl_t[c("Environments", "Percent_WT", "Percent_Non-lineage", "Percent_Lineage")]

# Transform the dataframe by environmets
wt_lin_envs_t <- wt_lin_envs %>%
  pivot_longer(cols = -Environments, names_to = "Status", values_to = "Percent")

wt_lin_envs_t <- as.data.frame(wt_lin_envs_t)

merged_wt_lin_envs_t <-  wt_lin_envs_t %>% inner_join(wt_lin_envs_count_t, by = c("Environments", "Status"))

merged_wt_lin_envs_t <- merged_wt_lin_envs_t %>%
  mutate(Count = ifelse(Count == 0, "", Count))
  
# Rename Clinical unknown and human associated environment
merged_wt_lin_envs_t$Environments <- gsub("Clinical..Unknown", "Clinical U", merged_wt_lin_envs_t$Environments)
merged_wt_lin_envs_t$Environments <- gsub("Human.associated.environmental", "Human Assoc", merged_wt_lin_envs_t$Environments)

merged_wt_lin_envs_t$Status <- gsub("_", " ", merged_wt_lin_envs_t$Status)

# Plot
plot <- ggplot(merged_wt_lin_envs_t, aes(x = Environments, y = Percent, fill = fct_rev(Status))) +
  geom_bar(position="stack", stat="identity", color = "black")  +
  geom_text(aes(label = Count), color = "white", size = 4, position = position_stack(vjust = 0.5), 
            family = "Times New Roman") +
  scale_fill_manual(values = c("Percent WT" = "blue", "Percent Lineage" = "darkred", "Percent Non-lineage" = "red"), name = "Lineage Status") +
  xlab("Environments") +
  ylab("Percent") +
  theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

bar_plot <- paste0(base, "_environment_bar.svg")

ggsave(bar_plot, plot = plot, width = 10, height = 7.5, units = "in")
