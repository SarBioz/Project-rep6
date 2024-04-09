#chris and mine data comparison--------------------------------------------------------------------------------------------------------------
Class_0_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_0.csv")
Class_1_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_1.csv")
Class_2_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2.csv")
Class_3_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3.csv")
Class_4_s <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4.csv")
chris_genes <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/Chris_result/ClassCombinations_TissuesCells.human.csv")

Class_3_s$paralogs <- paste(Class_3_s$name, Class_3_s$paralogue_name, sep = ";")
Class_2_s$paralogs <- paste(Class_2_s$name, Class_2_s$paralogue_name, sep = ";")
Class_4_s$paralogs <- paste(Class_4_s$name, Class_4_s$paralogue_name, sep = ";")
Class_1_s$paralogs <- paste(Class_1_s$name, Class_1_s$paralogue_name, sep = ";")
Class_0_s$paralogs <- paste(Class_0_s$name, Class_0_s$paralogue_name, sep = ";")

library(dplyr)
Class_0_ch <-chris_genes %>%filter(cells == 'Class_0')
Class_1_ch <-chris_genes %>%filter(cells == 'Class_1')
Class_2_ch <-chris_genes %>%filter(cells == 'Class_2')
Class_3_ch <-chris_genes %>%filter(cells == 'Class_3')
Class_4_ch <-chris_genes %>%filter(cells == 'Class_4')

colnames(Class_0_ch)[1] <- "paralogs"
colnames(Class_1_ch)[1] <- "paralogs"
colnames(Class_2_ch)[1] <- "paralogs"
colnames(Class_3_ch)[1] <- "paralogs"
colnames(Class_4_ch)[1] <- "paralogs"

# Find counts of classifications
class_0_count_ch <- nrow(Class_0_ch)
class_1_count_ch <- nrow(Class_1_ch)
class_2_count_ch <- nrow(Class_2_ch)
class_3_count_ch <- nrow(Class_3_ch)
class_4_count_ch <- nrow(Class_4_ch)
class_counts <- c(class_0_count_ch, class_1_count_ch, class_2_count_ch, class_3_count_ch, class_4_count_ch)
names(class_counts) <- c('Class 0', 'Class 1', 'Class 2', 'Class 3', 'Class 4')

# Bar plot
x <- barplot(class_counts, ylim = c(0, max(class_counts) + 500))
y = as.matrix(class_counts)
text(x, y+2, labels = as.character(y), pos = 3)


find_common_values <- function(df1, df2, column_name) {
  common_values <- intersect(df1[[column_name]], df2[[column_name]])
  return(data.frame(common_values))
}
common_class_0 <- find_common_values(Class_0_s, Class_0_ch, "paralogs")
common_class_1 <- find_common_values(Class_1_s, Class_1_ch, "paralogs")
common_class_2 <- find_common_values(Class_2_s, Class_2_ch, "paralogs")
common_class_3 <- find_common_values(Class_3_s, Class_3_ch, "paralogs")
common_class_4 <- find_common_values(Class_4_s, Class_4_ch, "paralogs")

common_class_0_count <- nrow(common_class_0)
common_class_1_count <- nrow(common_class_1)
common_class_2_count <- nrow(common_class_2)
common_class_3_count <- nrow(common_class_3)
common_class_4_count <- nrow(common_class_4)
class_counts <- c(common_class_0_count,common_class_1_count,common_class_2_count, common_class_3_count, common_class_4_count)
names(class_counts) <- c('Class 0','Class 1','Class 2', 'Class 3', 'Class 4')

par(mar = c(5, 5, 2, 2))
# Bar plot
x <- barplot(class_counts, ylim = c(0, max(class_counts) + 500))
y = as.matrix(class_counts)
text(x, y+2, labels = as.character(y), pos = 3)

write.csv(common_class_0, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_0_sarv_chris.csv", row.names = FALSE)
write.csv(common_class_1, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_1_sarv_chris.csv", row.names = FALSE)
write.csv(common_class_2, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_2_sarv_chris.csv", row.names = FALSE)
write.csv(common_class_3, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_3_sarv_chris.csv", row.names = FALSE)
write.csv(common_class_4, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_4_sarv_chris.csv", row.names = FALSE)




#chris and matt data comparison--------------------------------------------------------------------------------------------------------------
Class_2_m <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_2_Mat.csv")
Class_3_m <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_3_Mat.csv")
Class_4_m <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/Class_4_Mat.csv")

find_common_values <- function(df1, df2, column_name) {
  common_values <- intersect(df1[[column_name]], df2[[column_name]])
  return(data.frame(common_values))
}

common_class_2_ <- find_common_values(Class_2_m, Class_2_ch, "paralogs")
common_class_3_ <- find_common_values(Class_3_m, Class_3_ch, "paralogs")
common_class_4_ <- find_common_values(Class_4_m, Class_4_ch, "paralogs")


common_class_2_count1 <- nrow(common_class_2_)
common_class_3_count1 <- nrow(common_class_3_)
common_class_4_count1 <- nrow(common_class_4_)

class_counts_ <- c(common_class_2_count1, common_class_3_count1, common_class_4_count1)
names(class_counts_) <- c('Class 2', 'Class 3', 'Class 4')

par(mar = c(5, 5, 2, 2))
# Bar plot
x <- barplot(class_counts_, ylim = c(0, max(class_counts_) + 500))
y = as.matrix(class_counts_)
text(x, y+2, labels = as.character(y), pos = 3)

write.csv(common_class_2_, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_2_matt_chris.csv", row.names = FALSE)
write.csv(common_class_3_, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_3_matt_chris.csv", row.names = FALSE)
write.csv(common_class_4_, "/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/gene_pair_classification/common_Genes_class_4_matt_chris.csv", row.names = FALSE)


   
