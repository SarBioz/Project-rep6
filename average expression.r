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

# Filter the Seurat object to include only cells where cluster_species is "human_10x"
# Subset based on species
human_10x_data <- subset(obj, species_tech == "human_10x")
# Calculate average expression within each cluster for the filtered data
human_10x_cluster_means <- Seurat::AverageExpression(human_10x_data, by.group = "clusters")
human_10x_cluster_means_df <- as.data.frame(human_10x_cluster_means)
# Print the result
print(human_10x_cluster_means)

# 3) genes expressed in more than 1 and more than 2 cell after clustring

genes_in_seurat <- rownames(human_10x_cluster_means_df)
genepir_new <- genepair %>% select(name, paralogue_name)
human_10x_cluster_means_df$gene_name_column_in_obj <- rownames(human_10x_cluster_means_df)
human_10x_cluster_means_df_name <- human_10x_cluster_means_df %>% select(gene_name_column_in_obj)
filtered_gene_pairs_data <- genepir_new[genepair$name %in% genes_in_seurat, ]

# Merge gene pairs data with Seurat object based on gene names
merged_data <- merge(x = human_10x_cluster_means_df_name, y = filtered_gene_pairs_data, by.x = "gene_name_column_in_obj", by.y = "name")
x_data <-human_10x_cluster_means_df[, -ncol(human_10x_cluster_means_df)]
value_data <- ifelse(x_data > 1, 1, 0)
# Calculate row sums
row_sums <- rowSums(value_data)
# Identify rows with sum greater than 1
more_than_1_rows <- names(row_sums[row_sums > 1])
# Create a new dataframe
more_than_1 <- data.frame(RowName = more_than_1_rows)
# Identify rows with sum greater than 2
more_than_2_rows <- names(row_sums[row_sums > 2])
# Create a new dataframe
more_than_2 <- data.frame(RowName = more_than_2_rows)


genes_at_least_one <- Genes_At_Least_One_Cell_name$Genes_At_Least_One_Cell
more_than_one <- more_than_1$RowName
genes_at_least_two <- Genes_At_Least_Two_Cells_name$Genes_At_Least_Two_Cells
more_than_two <- more_than_2$RowName



y = list(genes_at_least_two, more_than_two , genes_at_least_one , more_than_one )
venn.plot <- venn.diagram(
  x = y ,
  category.names = c("at_least2_before" , "at_least2_after " , "at_least1_before" , "at_least1_after " ),
  filename = NULL,
  output=TRUE
)

grid.draw(venn.plot)
