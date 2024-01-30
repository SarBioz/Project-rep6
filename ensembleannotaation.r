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
# Run UMAP on the Seurat object
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
umap_coords <- as.data.frame(Embeddings(obj, reduction = "umap.unintegrated")[, 1:2])
clusters <- as.numeric(Idents(obj))
# Filter the Seurat object to include only cells where cluster_species is "human_10x"
human_10x_data <- subset(obj, subset = species_tech == "human_10x")
# Calculate average expression within each cluster for the filtered data
human_10x_cluster_means <- Seurat::AverageExpression(human_10x_data, by.group = "clusters")
human_10x_cluster_means_df <- as.data.frame(human_10x_cluster_means)
# Print the result
print(human_10x_cluster_means)
genes expressed in more than 1 and more than 2 cell before clustring

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
# 2. How many genes have expression in at least two cells > 1
genes_at_least_2 <- rownames(cluster_gene_expression)[binary_table2 >= 2]

Genes_At_Least_1Cell_name <- data.frame(
  gene_name = genes_at_least1)
Genes_At_Least_2Cells_name <- data.frame(
  gene_name = genes_at_least_2)

genes_at_least_one <- Genes_At_Least_One_Cell_name$Genes_At_Least_One_Cell
more_than_one <- Genes_At_Least_1Cell_name$Genes_At_Least_1Cell
genes_at_least_two <- Genes_At_Least_Two_Cells_name$Genes_At_Least_Two_Cells
more_than_two <- Genes_At_Least_2Cells_name$Genes_At_Least_2Cells

# Function to update transcript types------------------
update_transcript_types <- function(gene_df, gene_data) {
  transcript_types <- character(nrow(gene_df))

  for (i in seq_len(nrow(gene_df))) {
    gene_name <- gene_df$gene_name[i]
    match_index <- match(gene_name, gene_data$Gene.name)
    if (length(match_index) == 0){
        synonym_index <- match(gene_name, gene_data$Gene.Synonym)
        if (length(synonym_index) > 0) 
          transcript_types[i] <- gene_data$Transcript.type[synonym_index]
        else transcript_types[i] <- NA  
    }
    else 
        transcript_types[i] <- gene_data$Transcript.type[match_index]
  }
  
  gene_df$transcript_type <- transcript_types
  
  return(gene_df)
}

# Apply the function to each data frame
Genes_At_Least_One_Cell_names <- update_transcript_types(Genes_At_Least_One_Cell_name, geneanno)
Genes_At_Least_Two_Cells_names <- update_transcript_types(Genes_At_Least_Two_Cells_name, geneanno)
Genes_At_Least_1Cell_names <- update_transcript_types(Genes_At_Least_1Cell_name, geneanno)
Genes_At_Least_2Cells_names <- update_transcript_types(Genes_At_Least_2Cells_name, geneanno)
#--------------
# Function to create bar charts
create_bar_chart <- function(data, title) {
  ggplot(data, aes(x = transcript_type, fill = transcript_type)) +
    geom_bar() +
    geom_text(stat="count", aes(label = after_stat(count)), vjust=-0.5) +
    theme_minimal() +
    labs(title = title, 
         x = "Transcript Type", 
         y = "Number of Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create and display each bar chart separately
bar_chart_1 <- create_bar_chart(Genes_At_Least_One_Cell_names, "Genes At Least Once Before clustering")
print(bar_chart_1)

bar_chart_2 <- create_bar_chart(Genes_At_Least_Two_Cells_names, "Genes At Least Twice Before clustering")
print(bar_chart_2)

bar_chart_3 <- create_bar_chart(Genes_At_Least_1Cell_names, "Genes At Least 1 Cell after clustering")
print(bar_chart_3)

bar_chart_4 <- create_bar_chart(Genes_At_Least_2Cells_names, "Genes At Least 2 Cells after clustering")
print(bar_chart_4)

