genes_in_seurat <- rownames(human_10x_cluster_means_df)
genepir_new <- genepair %>% select(name, paralogue_name)
human_10x_cluster_means_df$gene_name_column_in_obj <- rownames(human_10x_cluster_means_df)
human_10x_cluster_means_df_name <- human_10x_cluster_means_df %>% select(gene_name_column_in_obj)
filtered_gene_pairs_data <- genepir_new[genepair$name %in% genes_in_seurat, ]

# Merge gene pairs data with Seurat object based on gene names
merged_data <- merge(x = human_10x_cluster_means_df_name, y = filtered_gene_pairs_data, by.x = "gene_name_column_in_obj", by.y = "name")

x_data <-human_10x_cluster_means_df[, -ncol(human_10x_cluster_means_df)]


calculate_ratio <- function(pair) {
  gene1 <- pair[1]
  gene2 <- pair[2]
  ratio <- x_data[gene1, ] / x_data[gene2, ]
  return(ratio)
}
ratio_data <- apply(merged_data, 1, calculate_ratio)
combined_df <- do.call(rbind, ratio_data)
combined_df$gene_pairs <- merged_data$sorted_pair

last_column <- ncol(combined_df)

# Set the last column as row names
row.names(combined_df) <- combined_df[, last_column]

# Remove the last column (if you want)
combined_df <- combined_df[, -last_column]


binary_df <- ifelse(combined_df >= 1/9 & combined_df <= 9, 1, 0)


column_sums <- colSums(binary_df, na.rm = TRUE)
sum_data <- data.frame(Column_Name = names(column_sums), Sum_Value = column_sums)
library(ggplot2)

ggplot(sum_data, aes(x = Column_Name, y = Sum_Value)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Sum of Values in Each Column",
       x = "Column Name",
       y = "Sum Value") +
  theme_minimal()
