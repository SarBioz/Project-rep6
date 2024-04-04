library(Seurat)
library(SeuratObject)
library(data.table)



# Init DataFrame to contain cell ID and clustering info
clustered_exp <- data.frame(sample_id = umi@meta.data$sample_id,
                            neighborhood = umi@meta.data$neighborhood,
                            subclass = umi@meta.data$confirmed_subclass,
                            cluster = umi@meta.data$cluster
)

# Extract expression matrix from SeuratObject and coerce sparse matrix to full
exp_matrix <- as.matrix(umi@assays[["RNA"]]@counts)
exp_matrix <- t(exp_matrix) # Transpose for formatting

# Create list of human genes from Ensembl data
gene_names <- geneanno$Gene.name
gene_names <- gene_names[nzchar(gene_names)] # Remove empty characters

# Create list of genes found in both Ensembl list and in data
genes_in_both <- intersect(gene_names, colnames(exp_matrix))

# Subset exp matrix to include only genes found in Ensembl
exp_matrix_df <- as.data.frame(exp_matrix)
exp_matrix_reduced <- exp_matrix_df[names(exp_matrix_df) %in% genes_in_both]

# Combine cluster data with exp matrix
clustered_exp <- cbind(clustered_exp, exp_matrix_reduced)
exp_data <- clustered_exp
# Expects: string containing target subclass name
# Returns: vector containing avg expression for each gene across target subclass
avg_exp_for_subclass <- function(tgt_subclass) {
  subclass_exp_data <- subset(exp_data, subclass == tgt_subclass)
  subclass_exp_data <- subclass_exp_data[,5:ncol(subclass_exp_data)] # Remove clustering data
  avg_exp <- apply(subclass_exp_data, 2, mean)
  return(avg_exp)
}

# Generate lists of gene and subclass names
gene_names <- colnames(exp_data)[5:length(colnames(exp_data))]
subclass_names <- unique(exp_data$subclass)

# Init dataframe for avg expression across subclasses
avg_exp_data <- data.frame(gene = gene_names)

# Iterate through subclasses, append avg expression data for genes
for(i in 1:length(subclass_names)) {
  avg_exp_data[subclass_names[i]] <- avg_exp_for_subclass(subclass_names[i])
}


# Read in paralog list
paralogs_full_list <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/genepairs/Human_GeneParalogPairs_v104.csv")
row.names(avg_exp_data) <- avg_exp_data[,1]
avg_exp_data <- avg_exp_data[, -1]
# Only use paired paralogs
paralogs <- paralogs_full_list[paralogs_full_list$pairs == 'pair', ]
paralogs <- paralogs[paralogs$name != "" & paralogs$paralogue_name != "", ] # Filter rows with empty names
paralogs <- subset(paralogs, name %in% rownames(avg_exp_data)
                   & paralogue_name %in%rownames(avg_exp_data)) # Remove paralog pairs not in dataset

# Init df for classification output
paralog_classifications <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(paralog_classifications) <- c('paralogs', 'class')

# Perform classification
for (i in 1:nrow(paralogs)) { # Iterate through list of paralogs
  para_A <- paralogs$name[i]
  para_B <- paralogs$paralogue_name[i]
  
  A_exp_by_subclass <- c() # Init empty vectors
  B_exp_by_subclass <- c()
  for (subclass in colnames(avg_exp_data)) { # For each paralog pair, iterate through subclasses
    A_exp <- avg_exp_data[para_A, subclass]
    B_exp <- avg_exp_data[para_B, subclass]
    
    A_exp_by_subclass <- c(A_exp_by_subclass, A_exp)
    B_exp_by_subclass <- c(B_exp_by_subclass, B_exp)
  }
  
  # Add arbitrary small value to prevent dividing by zero (*only use after checking for zero vals)
  exp_ratios <- (A_exp_by_subclass + 0.0000001) / (B_exp_by_subclass + 0.0000001)
  
  if(sum(A_exp_by_subclass) == 0 & sum(B_exp_by_subclass) == 0){
    class_tmp <- 0
  } else if(sum(A_exp_by_subclass) == 0 | sum(B_exp_by_subclass) == 0){
    class_tmp <- 1
  } else if(all(exp_ratios >= (1/9)) & all(exp_ratios <= 9)) {
    class_tmp <- 2
  } else if(xor(any(exp_ratios > 9), any(exp_ratios < (1/9)))){
    class_tmp <- 3
  } else if(any(exp_ratios > 9) & any(exp_ratios < (1/9))){
    class_tmp <- 4
  }
  row_tmp <- data.frame(paralogs = paste(para_A, para_B, sep=';'), class = class_tmp)
  paralog_classifications <- rbind(paralog_classifications, row_tmp)
}

# Find counts of classifications
class_0_count <- length(which(paralog_classifications$class == 0))
class_1_count <- length(which(paralog_classifications$class == 1))
class_2_count <- length(which(paralog_classifications$class == 2))
class_3_count <- length(which(paralog_classifications$class == 3))
class_4_count <- length(which(paralog_classifications$class == 4))
class_counts <- c(class_0_count, class_1_count, class_2_count, class_3_count, class_4_count)
names(class_counts) <- c('Class 0', 'Class 1', 'Class 2', 'Class 3', 'Class 4')

# Bar plot
x <- barplot(class_counts, ylim = c(0, max(class_counts) + 500))
y = as.matrix(class_counts)
text(x, y+2, labels = as.character(y), pos = 3)




Class_2_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2.csv")
Class_3_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3.csv")
Class_4_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4.csv")
library(dplyr)
Class_2_m <-paralog_classifications %>%
  filter(class == 2)
Class_3_m <-paralog_classifications %>%
  filter(class == 3)
Class_4_m <-paralog_classifications %>%
  filter(class == 4)

write.csv(Class_2_m, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2_Mat.csv", row.names = FALSE)
write.csv(Class_3_m, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3_Mat.csv", row.names = FALSE)
write.csv(Class_4_m, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4_Mat.csv", row.names = FALSE)




Class_3_s$paralogs <- paste(Class_3_s$name, Class_3_s$paralogue_name, sep = ";")
Class_2_s$paralogs <- paste(Class_2_s$name, Class_2_s$paralogue_name, sep = ";")
Class_4_s$paralogs <- paste(Class_4_s$name, Class_4_s$paralogue_name, sep = ";")

find_common_values <- function(df1, df2, column_name) {
  common_values <- intersect(df1[[column_name]], df2[[column_name]])
  return(data.frame(common_values))
}
common_class_2 <- find_common_values(Class_2_s, Class_2_m, "paralogs")
common_class_3 <- find_common_values(Class_3_s, Class_3_m, "paralogs")
common_class_4 <- find_common_values(Class_4_s, Class_4_m, "paralogs")

common_class_2_count <- nrow(common_class_2)
common_class_3_count <- nrow(common_class_3)
common_class_4_count <- nrow(common_class_4)
class_counts <- c(common_class_2_count, common_class_3_count, common_class_4_count)
names(class_counts) <- c('Class 2', 'Class 3', 'Class 4')

par(mar = c(5, 5, 2, 2))
# Bar plot
x <- barplot(class_counts, ylim = c(0, max(class_counts) + 500))
y = as.matrix(class_counts)
text(x, y+2, labels = as.character(y), pos = 3)
write.csv(common_class_2, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_2.csv", row.names = FALSE)
write.csv(common_class_3, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_3.csv", row.names = FALSE)
write.csv(common_class_4, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_4.csv", row.names = FALSE)
