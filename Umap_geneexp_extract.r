library(tibble)
library(ggforce)
library(dplyr)
library(tidyverse)
library(scattermore)

umi <- readRDS(file.path("/project/scRNAseq/Sar/2", "human_SCT_UMI_expression_matrix.RDS"))

obj <- FindVariableFeatures(umi, selection.method = "vst")
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 2)

# Run UMAP on the Seurat object
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
umap_coords <- as.data.frame(Embeddings(obj, reduction = "umap.unintegrated")[, 1:2])

clusters <- as.numeric(Idents(obj))

metadata_info <- data.frame(
  species_tech = obj$species_tech,
  nCount_RNA = obj$nCount_RNA,
  sample_id = obj$sample_id,
  donor = obj$donor,
  cross_species_cluster = obj$cross_species_cluster,
  nFeature_RNA = obj$nFeature_RNA,
  neighborhood = obj$neighborhood,
  species = obj$species,
  cluster_color = obj$cluster_color,
  cross_species_cluster = obj$cross_species_cluster,
  cluster = obj$cluster, donor = obj$donor 
)


to_plot <- umap_coords %>%
  as_tibble(rownames = "sample_id") %>%
  filter(sample_id %in% metadata_info$sample_id) %>%
  left_join(metadata_info, by = "sample_id")

#plot 10x nuclei by cluster
tmp <- metadata_info %>% distinct(cluster, cluster_color)
cluster_colors <- tmp$cluster_color
names(cluster_colors) <- tmp$cluster

p1 <- to_plot %>%
  filter(species_tech %>% str_detect("human_10x")) %>%
  sample_n(n()) %>%
  ggplot() +
  geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = cluster), pointsize = 1) +
  scale_color_manual(values = cluster_colors) +
  
  labs(x = "",
       y = "") +
  
  theme(aspect.ratio = 1,
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p1



#---------------------------------------

# Filter the Seurat object to include only cells where cluster_species is "human_10x"
human_10x_data <- subset(obj, subset = species_tech == "human_10x")
# Calculate average expression within each cluster for the filtered data
human_10x_cluster_means <- Seurat::AverageExpression(human_10x_data, by.group = "clusters")

human_10x_cluster_means_df <- as.data.frame(human_10x_cluster_means)
# Print the result
print(human_10x_cluster_means)


