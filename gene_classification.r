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
genepair <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/genepairs/Human_GeneParalogPairs_v104.csv")
genepair <- genepair[genepair$pairs=="pair",] 
human_10x_cluster_means_df$gene_name <- rownames(human_10x_cluster_means_df)
update_transcript_types <- function(gene_df, gene_data) {
  # Merge based on gene name
  merged_df <- merge(gene_df, gene_data, by.x = "gene_name", by.y = "Gene.name", all.x = FALSE)
  
  # If there are unmatched gene names, merge based on synonyms
  unmatched_rows <- setdiff(gene_df$gene_name, merged_df$gene_name) 
  
  if (length(unmatched_rows)) {
    unmatched_df <- as.data.frame(gene_df[gene_df$gene_name %in% unmatched_rows, ], stringsAsFactors = FALSE)
    colnames(unmatched_df) <- "gene_name"
    merged_df_synonyms <- merge(unmatched_df , gene_data, 
                                by.x = "gene_name", by.y = "Gene.Synonym", all.x = FALSE)
  }
  
  return (bind_rows(merged_df, merged_df_synonyms))
}




# Apply the function to each dataset
avgExpData <- update_transcript_types(human_10x_cluster_means_df, geneanno)
avgExpData<- subset(avgExpData, Transcript.type == "protein_coding")

avgExpData <- merge(human_10x_cluster_means_df, geneanno, by.x = "gene_name", by.y = "Gene.name", all.x = FALSE)
avgExpData<- subset(avgExpData, Transcript.type == "protein_coding")




df_no_duplicates <- avgExpData[!duplicated(avgExpData$gene_name), ]
df_no_duplicates <- df_no_duplicates[, !(names(df_no_duplicates) %in% c("Gene.Synonym"))]
df_no_duplicates$Gene.Synonym <- avgExpData_name$Gene.Synonym
rownames(df_no_duplicates) <- df_no_duplicates$gene_name
#Class 0 (neither paralog expressed)--------------------------
not_expressed <- data.frame(apply(df_no_duplicates[, 2:34], 1, sum),df_no_duplicates$Gene.Synonym,df_no_duplicates$Gene.stable.ID)
colnames(not_expressed) <- c("sum","Gene.Synonym","Gene.stable.ID")
Class_0_values <- data.frame(gene_name = rownames(not_expressed[not_expressed$sum == 0, , drop = FALSE]), Gene.stable.ID = not_expressed$Gene.stable.ID[not_expressed$sum == 0])



genepair_0 <-genepair
zero_sum_names <- Class_0_values$Gene.stable.ID
genepair_0$zero_sum <- genepair$ensembl_gene_id %in% zero_sum_names & genepair$paralogue_ensembl_gene_id %in% zero_sum_names
genepair_0 <- subset(genepair_0, zero_sum == TRUE)

#Class 1 (one paralog expressed)-----------------------------------------
genepair_1 <-genepair
names <- not_expressed[c("Gene.stable.ID","sum")]
non_zero_sum <- names[names$sum > 0, ]
zero_sum <- names[names$sum == 0, ]
genepair_1$class1 <- ((genepair$ensembl_gene_id %in% non_zero_sum$Gene.stable.ID) & (genepair$paralogue_ensembl_gene_id %in% zero_sum$Gene.stable.ID)) | ((genepair$ensembl_gene_id %in% zero_sum$Gene.stable.ID) & (genepair$paralogue_ensembl_gene_id %in% non_zero_sum$Gene.stable.ID))
genepair_1 <- subset(genepair_1, class1 == TRUE)

# Class 2
genepair_2 <-genepair
genepair_2$class2 <- ((genepair$ensembl_gene_id %in% non_zero_sum$Gene.stable.ID) & (genepair$paralogue_ensembl_gene_id %in% non_zero_sum$Gene.stable.ID))
genepair_2 <- subset(genepair_2, class2 == TRUE)
genepair_2 <- genepair_2[, c("name", "paralogue_name", "ensembl_gene_id", "paralogue_ensembl_gene_id")]
df_cluster <- df_no_duplicates[, 2:36]

 
pair_cluster <- df_cluster[match(genepair_2$ensembl_gene_id, df_cluster$Gene.stable.ID), ]
para_cluster <- df_cluster[match(genepair_2$paralogue_ensembl_gene_id, df_cluster$Gene.stable.ID), ]

ratio <- pair_cluster[, 1:34]/ para_cluster[, 1:34]
ratio <- ifelse(ratio >= 1/9 & ratio <= 9, 1, 0)

genepair_2$ratio <- rowSums(ratio)
genepair_2 <- subset(genepair_2, ratio == 34)


# Class 3-------------------------------------------------------------------------------------------------
genepair_3 <-genepair
genepair_3$class3 <- ((genepair$ensembl_gene_id %in% non_zero_sum$Gene.stable.ID) & (genepair$paralogue_ensembl_gene_id %in% non_zero_sum$Gene.stable.ID))
genepair_3 <- subset(genepair_3, class3 == TRUE)
genepair_3 <- genepair_3[, c("name", "paralogue_name", "ensembl_gene_id", "paralogue_ensembl_gene_id")]
df_cluster <- df_no_duplicates[, 2:36]


pair_cluster <- df_cluster[match(genepair_3$ensembl_gene_id, df_cluster$Gene.stable.ID), ]
para_cluster <- df_cluster[match(genepair_3$paralogue_ensembl_gene_id, df_cluster$Gene.stable.ID), ]

ratio <- pair_cluster[, 1:34]/ para_cluster[, 1:34]
ratio <- ifelse(ratio <= 1/9 | ratio >= 9, 1, 0)

genepair_3$ratio <- rowSums(ratio)
genepair_3 <- subset(genepair_3, ratio >= 1)

# Class 4------------------------------------------------------------------------------------------------
genepair_4 <-genepair
genepair_4$class4 <- ((genepair$ensembl_gene_id %in% non_zero_sum$Gene.stable.ID) & (genepair$paralogue_ensembl_gene_id %in% non_zero_sum$Gene.stable.ID))
genepair_4 <- subset(genepair_4, class4 == TRUE)
genepair_4 <- genepair_4[, c("name", "paralogue_name", "ensembl_gene_id", "paralogue_ensembl_gene_id")]
df_cluster <- df_no_duplicates[, 2:36]


pair_cluster <- df_cluster[match(genepair_4$ensembl_gene_id, df_cluster$Gene.stable.ID), ]
para_cluster <- df_cluster[match(genepair_4$paralogue_ensembl_gene_id, df_cluster$Gene.stable.ID), ]

ratio <- pair_cluster[, 1:34]/ para_cluster[, 1:34]
ratio <- ifelse(ratio <= 1/9, 1, 0)

genepair_4$ratio_l <- rowSums(ratio)

ratio <- pair_cluster[, 1:34]/ para_cluster[, 1:34]
ratio <- ifelse(ratio >= 9, 1, 0)

genepair_4$ratio_h <- rowSums(ratio)

genepair_4 <- subset(genepair_4, ratio_l > 0 & ratio_h > 0)

