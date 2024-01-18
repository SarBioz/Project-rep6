umi <- readRDS(file.path("human_SCT_UMI_expression_matrix.RDS"))
# 1) genes expressed in more than 1 and more than 2 cell before clustring

expression_matrix <- umi$RNA@counts
# 1. How many genes have expression in at least one cell > 1
genes_at_least_one_cell_1 <-  rownames(expression_matrix)[rowSums(expression_matrix > 1) > 0]
# 2. How many genes have expression in at least two cells > 1
genes_at_least_two_cells_1 <- rownames(expression_matrix)[rowSums(expression_matrix >= 1) > 1]

Genes_At_Least_One_Cell_name <- data.frame(
  Genes_At_Least_One_Cell = genes_at_least_one_cell_1)
Genes_At_Least_Two_Cells_name <- data.frame(
  Genes_At_Least_Two_Cells = genes_at_least_two_cells_1)

# 2) preprocessing step for seurat obj
obj <- FindVariableFeatures(umi)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.75)
# Run UMAP on the Seurat object
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
umap_coords <- as.data.frame(Embeddings(obj, reduction = "umap.unintegrated")[, 1:2])
clusters <- as.numeric(Idents(obj))

# Access gene expression data for each cluster
cluster_gene_expression <- human_10x_data@assays$RNA@data
# Accessing cluster information
cluster_info <- human_10x_data$cluster
# changing the colnames to the cluster 
colnames(cluster_gene_expression) <- cluster_info

binary_table2 <- rowSums(cluster_gene_expression > 1)
genes_at_least1 <-  rownames(cluster_gene_expression)[binary_table2 >= 1]
genes_at_least_2 <- rownames(cluster_gene_expression)[binary_table2 >= 2]

Genes_At_Least_1Cell_name <- data.frame(
  Genes_At_Least_1Cell = genes_at_least1)
Genes_At_Least_2Cells_name <- data.frame(
  Genes_At_Least_2Cells = genes_at_least_2)

genes_at_least_one <- Genes_At_Least_One_Cell_name$Genes_At_Least_One_Cell
more_than_one <- Genes_At_Least_1Cell_name$Genes_At_Least_1Cell
genes_at_least_two <- Genes_At_Least_Two_Cells_name$Genes_At_Least_Two_Cells
more_than_two <- Genes_At_Least_2Cells_name$Genes_At_Least_2Cells

library(VennDiagram)

y = list(genes_at_least_two, more_than_two , genes_at_least_one , more_than_one )
venn.plot <- venn.diagram(
  x = y ,
  category.names = c("at_least2_before" , "at_least2_after " , "at_least1_before" , "at_least1_after " ),
  filename = NULL,
  output=TRUE
)

grid.draw(venn.plot)


venn.plot <- venn.diagram(
  x = list(genes_at_least_two, more_than_two  ) ,
  category.names = c("at_least2_before" , "at_least2_after "),
  filename = NULL,
  output=TRUE
)

grid.draw(venn.plot)

venn.plot <- venn.diagram(
  x = list(genes_at_least_one , more_than_one ) ,
  category.names = c("at_least1_before","at_least1_after "),
  filename = NULL,
  output=TRUE
)

grid.draw(venn.plot)
