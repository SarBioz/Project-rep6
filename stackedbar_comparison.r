Class_0_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_0_Macaca.csv")
Class_1_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_1_Macaca.csv")
Class_2_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2_Macaca.csv")
Class_3_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3_Macaca.csv")
Class_4_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4_Macaca.csv")

chris_genes <- read.csv ("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/2/ClassCombinations_TissuesCells.macaque.csv")

Class_3_s$paralogs <- paste(Class_3_s$name, Class_3_s$paralogue_name, sep = ";")
Class_2_s$paralogs <- paste(Class_2_s$name, Class_2_s$paralogue_name, sep = ";")
Class_4_s$paralogs <- paste(Class_4_s$name, Class_4_s$paralogue_name, sep = ";")
Class_1_s$paralogs <- paste(Class_1_s$name, Class_1_s$paralogue_name, sep = ";")
Class_0_s$paralogs <- paste(Class_0_s$name, Class_0_s$paralogue_name, sep = ";")

library(dplyr)
Class_0_ch <-chris_genes %>%filter(cells == 'Class_0')
Class_1_ch <-chris_genes %>%filter(cells == 'Class_1')
Class_2_ch <-chris_genes %>%filter(cells == 'Class_2')
Class_3_ch <-chris_genes %>%filter(cells == 'Class_3')
Class_4_ch <-chris_genes %>%filter(cells == 'Class_4')

colnames(Class_0_ch)[1] <- "paralogs"
colnames(Class_1_ch)[1] <- "paralogs"
colnames(Class_2_ch)[1] <- "paralogs"
colnames(Class_3_ch)[1] <- "paralogs"
colnames(Class_4_ch)[1] <- "paralogs"


library(dplyr)

# Function to find common values and calculate percentages
find_common_values <- function(df1, df2, column_name) {
  common_values <- intersect(df1[[column_name]], df2[[column_name]])
  percentage_common <- length(common_values) / nrow(df2) * 100
  percentage_not_in_mine <- 100 - percentage_common
  return(list(common_values = common_values,
              percentage_common = percentage_common,
              percentage_not_in_mine = percentage_not_in_mine))
}

# Find common genes and percentages for each class
common_class_0 <- find_common_values(Class_0_s, Class_0_ch, "paralogs")
common_class_1 <- find_common_values(Class_1_s, Class_1_ch, "paralogs")
common_class_2 <- find_common_values(Class_2_s, Class_2_ch, "paralogs")
common_class_3 <- find_common_values(Class_3_s, Class_3_ch, "paralogs")
common_class_4 <- find_common_values(Class_4_s, Class_4_ch, "paralogs")

# Display results
print("Class 0:")
print(paste("Percentage of common genes:", common_class_0$percentage_common))
print(paste("Percentage of genes in Chris' data not in mine:", common_class_0$percentage_not_in_mine))

print("Class 1:")
print(paste("Percentage of common genes:", common_class_1$percentage_common))
print(paste("Percentage of genes in Chris' data not in mine:", common_class_1$percentage_not_in_mine))

print("Class 2:")
print(paste("Percentage of common genes:", common_class_2$percentage_common))
print(paste("Percentage of genes in Chris' data not in mine:", common_class_2$percentage_not_in_mine))

print("Class 3:")
print(paste("Percentage of common genes:", common_class_3$percentage_common))
print(paste("Percentage of genes in Chris' data not in mine:", common_class_3$percentage_not_in_mine))

print("Class 4:")
print(paste("Percentage of common genes:", common_class_4$percentage_common))
print(paste("Percentage of genes in Chris' data not in mine:", common_class_4$percentage_not_in_mine))


library(tidyverse)

# Create a data frame for percentages
data <- data.frame(
  sample = c("Class 0", "Class 1", "Class 2", "Class 3", "Class 4"),
  common_genes = c(common_class_0$percentage_common, 
                   common_class_1$percentage_common,
                   common_class_2$percentage_common,
                   common_class_3$percentage_common,
                   common_class_4$percentage_common),
  not_in_mine = c(common_class_0$percentage_not_in_mine, 
                   common_class_1$percentage_not_in_mine,
                   common_class_2$percentage_not_in_mine,
                   common_class_3$percentage_not_in_mine,
                   common_class_4$percentage_not_in_mine)
)

# Reshape data for plotting
data_long <- data %>% 
  pivot_longer(cols = c(common_genes, not_in_mine), 
               names_to = "status", 
               values_to = "percentage")

# Plot stacked bar chart
ggplot(data_long, aes(fill = status, x = sample, y = percentage)) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Percentage of Common Genes and Genes in Chris' Data Not in Mine",
       x = "Class",
       y = "Percentage") +
  scale_fill_manual(values = c("common_genes" = "blue", "not_in_mine" = "red")) +
  theme_minimal()

