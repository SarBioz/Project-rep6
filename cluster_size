
subclass_data<- unique(obj$confirmed_subclass)
# Initialize an empty data frame to store results
cluster_subclass <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(cluster_subclass) <- c('cluster', 'num_of_cells', 'median_reads', 'mean_reads')

# Loop through each subclass
for (i in 1:length(subclass_data)) {
  tmp_df <- data.frame(matrix(nrow = 1, ncol = 4))
  colnames(tmp_df) <- c('cluster', 'num_of_cells', 'median_reads', 'mean_reads')
  
  # Subset the Seurat object by the current subclass
  cluster_exp <- subset(obj, confirmed_subclass == subclass_data[i])
  
  # Extract the expression matrix for the subset
  cluster_mtx <- as.matrix(cluster_exp@assays$RNA@data)
  
  # Calculate the required statistics
  tmp_df$cluster <- subclass_data[i]
  tmp_df$num_of_cells <- ncol(cluster_mtx)  # Number of cells is the number of columns
  tmp_df$median_reads <- median(cluster_mtx)
  tmp_df$mean_reads <- mean(cluster_mtx)
  
  # Append the results to the final data frame
  cluster_subclass <- rbind(cluster_subclass, tmp_df)
}

# Print the result
print(cluster_subclass)

write.csv(cluster_subclass, '/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/subclass.csv', row.names = FALSE)


plot <- ggplot(data = cluster_subclass, aes(x = subclass, y = num_of_cells)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
plot

plot <- ggplot(data = cluster_subclass, aes(x = subclass, y = mean_reads)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
plot
