Class_0_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_0.csv")
Class_1_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_1.csv")
Class_2_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2.csv")
Class_3_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3.csv")
Class_4_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4.csv")


Class_0_s_chimp <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_0_chimp.csv")
Class_1_s_chimp <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_1_chimp.csv")
Class_2_s_chimp <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2_chimp.csv")
Class_3_s_chimp <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3_chimp.csv")
Class_4_s_chimp <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4_chimp.csv")

Class_0_s_macaca <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_0_Macaca.csv")
Class_1_s_macaca <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_1_Macaca.csv")
Class_2_s_macaca <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2_Macaca.csv")
Class_3_s_macaca <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3_Macaca.csv")
Class_4_s_macaca <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4_Macaca.csv")
   

chris_genes_human <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/Chris_result/ClassCombinations_TissuesCells.human.csv")
chris_genes_macaca <- read.csv ("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/2/ClassCombinations_TissuesCells.macaque.csv")
chris_genes_chimp <- read.csv ("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/2/ClassCombinations_TissuesCells.chimp.csv")

# List of column names to remove
columns_to_remove <- c("tissues", "bulk","class","gtex","species","cac","cer","cn")

# Remove columns based on their names
chris_genes_human <- chris_genes_human[, !names(chris_genes_human) %in% columns_to_remove]
chris_genes_macaca <- chris_genes_macaca[, !names(chris_genes_macaca) %in% columns_to_remove]
chris_genes_chimp <- chris_genes_chimp[, !names(chris_genes_chimp) %in% columns_to_remove]
# Rename a column
names(chris_genes_human)[names(chris_genes_human) == "cells"] <- "chris"
names(chris_genes_macaca)[names(chris_genes_macaca) == "cells"] <- "chris"
names(chris_genes_chimp)[names(chris_genes_chimp) == "cells"] <- "chris"

add_paralogs_and_class <- function(data, name_col, paralogue_name_col, class_name) {
  paralogs <- paste(data[[name_col]], data[[paralogue_name_col]], sep = ";")
  class_col <- rep(class_name, nrow(data))
  result <- data.frame(paralogs = paralogs, class = class_col)
  return(result)
}

Class_3_s <- add_paralogs_and_class(Class_3_s, "name", "paralogue_name", "class 3")
Class_2_s <- add_paralogs_and_class(Class_2_s, "name", "paralogue_name", "class 2")
Class_1_s <- add_paralogs_and_class(Class_1_s, "name", "paralogue_name", "class 1")
Class_4_s <- add_paralogs_and_class(Class_4_s, "name", "paralogue_name", "class 4")
Class_0_s <- add_paralogs_and_class(Class_0_s, "name", "paralogue_name", "class 0")

Class_3_s_chimp <- add_paralogs_and_class(Class_3_s_chimp, "name", "paralogue_name", "class 3")
Class_2_s_chimp <- add_paralogs_and_class(Class_2_s_chimp, "name", "paralogue_name", "class 2")
Class_1_s_chimp <- add_paralogs_and_class(Class_1_s_chimp, "name", "paralogue_name", "class 1")
Class_4_s_chimp <- add_paralogs_and_class(Class_4_s_chimp, "name", "paralogue_name", "class 4")
Class_0_s_chimp <- add_paralogs_and_class(Class_0_s_chimp, "name", "paralogue_name", "class 0")

Class_3_s_macaca <- add_paralogs_and_class(Class_3_s_macaca, "name", "paralogue_name", "class 3")
Class_2_s_macaca <- add_paralogs_and_class(Class_2_s_macaca, "name", "paralogue_name", "class 2")
Class_4_s_macaca <- add_paralogs_and_class(Class_4_s_macaca, "name", "paralogue_name", "class 4")


S_data <- do.call(rbind, list(Class_3_s, Class_2_s, Class_4_s, Class_1_s, Class_0_s))
S_data_chimp <- do.call(rbind, list(Class_3_s_chimp, Class_2_s_chimp, Class_4_s_chimp, Class_1_s_chimp, Class_0_s_chimp))
S_data_macaca <- do.call(rbind, list(Class_3_s_macaca, Class_2_s_macaca, Class_4_s_macaca))
write.csv(S_data, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/human_class_merge_sar&chris.csv", row.names = FALSE)
write.csv(S_data_chimp, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/chimp_class_merge_sar&chris.csv", row.names = FALSE)
write.csv(S_data_macaca, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/macaca_class_merge_sar&chris.csv", row.names = FALSE)

colnames(S_data) <- c("genes", "sarv")
colnames(S_data_chimp) <- c("genes", "sarv")
colnames(S_data_macaca) <- c("genes", "sarv")

# Merge the two data frames based on the "genes" column
merged_data_human <- merge(chris_genes_human, S_data, by = "genes", all = TRUE)
merged_data_chimp <- merge(chris_genes_chimp, S_data_chimp, by = "genes", all = TRUE)
merged_data_macaca <- merge(chris_genes_macaca, S_data_macaca, by = "genes", all = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)

process_data <- function(data) {
  # Replace NA values with "NA"
  data[is.na(data)] <- "NA"
  
  # Create a contingency table
  contingency_table <- table(data$sarv, data$chris)
  
  # Convert the table to a dataframe
  matrix_df <- as.data.frame.matrix(contingency_table)
  
  # Rename the rows and columns
  rownames(matrix_df) <- c("class 0", "class 1", "class 2", "class 3", "class 4", "NA")
  colnames(matrix_df) <- c("class_0", "class_1", "class_2", "class_3", "class_4", "NA")
  
  # If there are NA values in the dataframe, replace them with 0
  matrix_df[is.na(matrix_df)] <- 0
  
  # Reshape the data frame and create the new columns
  new_df <- matrix_df %>%
    rownames_to_column(var = "rowname") %>%
    gather(key = "colname", value = "value", -rowname) %>%
    mutate(method = paste0(rowname, "_", colname),
           condition = gsub("_.*", "", method))
  
  return(new_df)
}


merged_data_human_processed <- process_data(merged_data_human)
merged_data_chimp_processed <- process_data(merged_data_chimp)

# Plotting
plot_human<- ggplot(merged_data_human_processed, aes(fill = method, y = value, x = condition)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Comparison of Gene Pair Classifications for human: Exploring Consistency Across Classifications ",
       x = "Sarvenaz Class",
       y = "number of genepiars",
       fill="Classification Comparison" )

# Plotting
plot_chimp<- ggplot(merged_data_chimp_processed, aes(fill = method, y = value, x = condition)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Comparison of Gene Pair Classifications for chimp: Exploring Consistency Across Classifications ",
       x = "Sarvenaz Class",
       y = "number of genepiars",
       fill="Classification Comparison" )

# Replace NA values with "NA"
merged_data_macaca[is.na(merged_data_macaca)] <- "NA"

# Create a contingency table
contingency_table <- table(merged_data_macaca$sarv, merged_data_macaca$chris)

# Convert the table to a dataframe
matrix_df <- as.data.frame.matrix(contingency_table)

# Rename the rows and columns
rownames(matrix_df) <- c("class 2", "class 3", "class 4", "NA")
colnames(matrix_df) <- c("class_0", "class_1", "class_2", "class_3", "class_4", "NA")

# If there are NA values in the dataframe, replace them with 0
matrix_df[is.na(matrix_df)] <- 0
# Reshape the data frame and create the new columns
new_df <- matrix_df %>%
  rownames_to_column(var = "rowname") %>%
  gather(key = "colname", value = "value", -rowname) %>%
  mutate(method = paste0(rowname, "_", colname),
         condition = gsub("_.*", "", method))

library(ggplot2)

plto_macaca <- ggplot(new_df, aes(fill = method, y = value, x = condition)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Comparison of Gene Pair Classifications for macaca: Exploring Consistency Across Classifications ",
       x = "Sarvenaz Class",
       y = "number of genepiars",
       fill="Classification Comparison" )


ggsave(filename = "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/plots/macaca_classcomparison_plot.png", plot = plto_macaca)
ggsave(filename = "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/plots/human_classcomparison_plot.png", plot = plot_human)
ggsave(filename = "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/plots/chimp_classcomparison_plot.png", plot = plot_chimp)

# Function to count occurrences of classes
count_classes <- function(data, class_column_name) {
  table(gsub("\\s+", "_", data[[class_column_name]]))  # Replace spaces with underscores for consistency
}

# Count occurrences of classes in each dataset
count_data_sarv <- count_classes(S_data, "sarv")
count_data_chris <- count_classes(chris_genes_human, "chris")

# Create a dataframe with counts
human_merged_data <- data.frame(nClass = names(count_data_sarv),
                          nchris = as.vector(count_data_chris),
                          nsarv = as.vector(count_data_sarv))


# Function to count occurrences of classes
count_classes <- function(data, class_column_name) {
  table(gsub("\\s+", "_", data[[class_column_name]]))  # Replace spaces with underscores for consistency
}

# Count occurrences of classes in each dataset
count_data_sarv <- count_classes(S_data_chimp, "sarv")
count_data_chris <- count_classes(chris_genes_chimp, "chris")

# Create a dataframe with counts
chimp_merged_data <- data.frame(nClass = names(count_data_sarv),
                                nchris = as.vector(count_data_chris),
                                nsarv = as.vector(count_data_sarv))



# Function to count occurrences of classes
count_classes <- function(data, class_column_name) {
  table(gsub("\\s+", "_", data[[class_column_name]]))  # Replace spaces with underscores for consistency
}

# Count occurrences of classes in each dataset
count_data_sarv <- count_classes(S_data_macaca, "sarv")
count_data_chris <- count_classes(chris_genes_macaca, "chris")

# Ensure consistent class names across datasets
all_classes <- union(names(count_data_sarv), names(count_data_chris))

# Create an empty dataframe with all classes
macaca_merged_data <- data.frame(nClass = all_classes, nchris = 0, nsarv = 0)

# Fill in counts for chris dataset
macaca_merged_data$nchris[match(names(count_data_chris), macaca_merged_data$nClass)] <- as.vector(count_data_chris)

# Fill in counts for sarv dataset
macaca_merged_data$nsarv[match(names(count_data_sarv), macaca_merged_data$nClass)] <- as.vector(count_data_sarv)

# Print the resulting dataframe
print(macaca_merged_data)


# Manually add 
macaca_merged_data <- data.frame(
  nClass = c("Class_0 ", "Class_1", "Class_2", "Class_3", "Class_4"),
  nchris = c(115, 396, 267, 1967, 398),
  nsarv = c(0, 0, 124, 2067, 1097)
)

