library(Seurat)
library(ggplot2)
library(stringr) 
library(patchwork)
library(tibble)
library(ggforce)
library(dplyr)
library(tidyverse)
library(scattermore)
library(patchwork)
geneanno <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/2/mart_export.txt")
umi <- readRDS(file.path("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/2", "human_SCT_UMI_expression_matrix.RDS"))
obj <- FindVariableFeatures(umi)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.75)
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
umap_coords <- as.data.frame(Embeddings(obj, reduction = "umap.unintegrated")[, 1:2])
clusters <- as.numeric(Idents(obj))

#genes expressed in more than 1 and more than 2 cell before clustring
expression_matrix <- umi$RNA@counts
binary_table <- rowSums(expression_matrix > 1)
# 1. How many genes have expression in at least one cell > 1
genes_at_least_one_cell_1 <-  rownames(expression_matrix)[binary_table >= 1]
# 2. How many genes have expression in at least two cells > 1
genes_at_least_two_cells_1 <- rownames(expression_matrix)[binary_table >= 2]

Genes_At_Least_One_Cell_name <- data.frame(
  gene_name = genes_at_least_one_cell_1)
Genes_At_Least_Two_Cells_name <- data.frame(
  gene_name = genes_at_least_two_cells_1)

# Access gene expression data for each cluster
cluster_gene_expression <- human_10x_data@assays$RNA@data
# Accessing cluster information
cluster_info <- human_10x_data$cluster
colnames(cluster_gene_expression) <- cluster_info
binary_table2 <- rowSums(cluster_gene_expression > 1)
genes_at_least1 <-  rownames(cluster_gene_expression)[binary_table2 >= 1]
# How many genes have expression in at least two cells > 1
genes_at_least_2 <- rownames(cluster_gene_expression)[binary_table2 >= 2]

Genes_At_Least_1Cell_name <- data.frame(
  gene_name = genes_at_least1)
Genes_At_Least_2Cells_name <- data.frame(
  gene_name = genes_at_least_2)

update_transcript_types <- function(gene_df, gene_data) {
  # Merge based on gene name
  merged_df <- merge(gene_df, gene_data, by.x = "gene_name", by.y = "Gene.name", all.x = FALSE)
  
  # If there are unmatched gene names, merge based on synonyms
  unmatched_rows <- is.na(merged_df$Transcript.type)
  if (any(unmatched_rows)) {
    merged_df_synonyms <- merge(gene_df[unmatched_rows, , drop = FALSE], gene_data, 
                                by.x = "gene_name", by.y = "Gene.Synonym", all.x = FALSE)
    
    # Update the Transcript.type for the unmatched rows
    merged_df$Transcript.type[unmatched_rows] <- merged_df_synonyms$Transcript.type
    
    # Identify still unmatched rows
    still_unmatched_rows <- is.na(merged_df$Transcript.type)
    
    # Save the unmatched rows in a separate data frame
    unmatched_data <- gene_df[unmatched_rows | still_unmatched_rows, , drop = FALSE]
  } else {
    unmatched_data <- NULL  # No unmatched data
  }
  
  return(list(merged_df, unmatched_data))
}

# Apply the function to each dataset
result1 <- update_transcript_types(Genes_At_Least_One_Cell_name, geneanno)
Genes_At_Least_One_Cell_names <- result1[[1]]
unmatched_data1 <- result1[[2]]

result2 <- update_transcript_types(Genes_At_Least_Two_Cells_name, geneanno)
Genes_At_Least_Two_Cells_names <- result2[[1]]
unmatched_data2 <- result2[[2]]

result3 <- update_transcript_types(Genes_At_Least_1Cell_name, geneanno)
Genes_At_Least_1Cell_names <- result3[[1]]
unmatched_data3 <- result3[[2]]

result4 <- update_transcript_types(Genes_At_Least_2Cells_name, geneanno)
Genes_At_Least_2Cells_names <- result4[[1]]
unmatched_data4 <- result4[[2]]
# Genes are classified as protein-coding if at least one of their transcript types falls under this category. For other transcript types, choose the one that appeared most frequently for each gene.
# save all of the gene synonyms for each gene and separate them by ","
process_data <- function(input_data) {
  result <- input_data %>%
    group_by(gene_name) %>%
    summarize(
      Transcript.type = ifelse(any(Transcript.type == "protein_coding"), "protein_coding", names(which.max(table(Transcript.type)))),
      Gene.Synonym = ifelse(all(is.na(Gene.Synonym)), NA, paste(na.omit(Gene.Synonym), collapse = ","))
    ) %>%
    ungroup()
  
  return(result)
}
# Apply the function to each dataset
At_Least_One_Cell_result <- process_data(Genes_At_Least_One_Cell_names)
Genes_At_Least_Two_Cells_result <- process_data(Genes_At_Least_Two_Cells_names)
Genes_At_Least_1Cell_result <- process_data(Genes_At_Least_1Cell_names)
Genes_At_Least_2Cells_result <- process_data(Genes_At_Least_2Cells_names)


# Function to create bar charts
create_bar_chart <- function(data, title) {
  ggplot(data, aes(x = Transcript.type, fill = Transcript.type)) +
    geom_bar() +
    geom_text(stat="count", aes(label = after_stat(count)), vjust=-0.5) +
    theme_minimal() +
    labs(title = title, 
         x = "Transcript Type", 
         y = "Number of Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create and display each bar chart separately
bar_chart_1 <- create_bar_chart(At_Least_One_Cell_result, "Genes At Least Once Before clustering")
print(bar_chart_1)

bar_chart_2 <- create_bar_chart(Genes_At_Least_Two_Cells_result, "Genes At Least Twice Before clustering")
print(bar_chart_2)

bar_chart_3 <- create_bar_chart(Genes_At_Least_1Cell_result, "Genes At Least 1 Cell after clustering")
print(bar_chart_3)

bar_chart_4 <- create_bar_chart(Genes_At_Least_2Cells_result, "Genes At Least 2 Cells after clustering")
print(bar_chart_4)

