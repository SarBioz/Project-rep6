protein_1before<- subset(At_Least_One_Cell_result, Transcript.type == "protein_coding")
protein_1after <- subset(Genes_At_Least_1Cell_result, Transcript.type == "protein_coding")
protein_2before <- subset(Genes_At_Least_Two_Cells_result, Transcript.type == "protein_coding")
protein_2after <- subset(Genes_At_Least_2Cells_result, Transcript.type == "protein_coding")


genepair <- read.csv("/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/genepairs/Human_GeneParalogPairs_v104.csv")

filtered_df_1before_gene <- (protein_1before[protein_1before$gene_name%in%genepair$name, ])
uniq_1before_gene <- data.frame(unique(filtered_df_1before_gene$gene_name))
filtered_df_1before_para <- (protein_1before[protein_1before$gene_name%in%genepair$paralogue_name, ])
uniq_1before_para <- data.frame(unique(filtered_df_1before_para$gene_name))
filtered_df_1after_gene <- (protein_1after[protein_1after$gene_name%in%genepair$name, ])
uniq_1after_gene <- data.frame(unique(filtered_df_1after_gene$gene_name))
filtered_df_1after_para <- (protein_1after[protein_1after$gene_name%in%genepair$paralogue_name, ])
uniq_1after_para <- data.frame(unique(filtered_df_1after_para$gene_name))
filtered_df_2before_gene <- (protein_2before[protein_2before$gene_name%in%genepair$name, ])
uniq_2before_gene <- data.frame(unique(filtered_df_2before_gene$gene_name))
filtered_df_2before_para <- (protein_2before[protein_2before$gene_name%in%genepair$paralogue_name, ])
uniq_2before_para <- data.frame(unique(filtered_df_2before_para$gene_name))
filtered_df_2after_gene <- (protein_2after[protein_2after$gene_name%in%genepair$name, ])
uniq_2after_gene <- data.frame(unique(filtered_df_2after_gene$gene_name))
filtered_df_2after_para <- (protein_2after[protein_2after$gene_name%in%genepair$paralogue_name, ])
uniq_2after_para <- data.frame(unique(filtered_df_2after_para$gene_name))




update_paralogue_gene <- function(gene_df, pair_df, paralog_df) {
  merged_data <- merge(gene_df, pair_df[, c("name", "paralogue_name")], by.x = "gene_name", by.y = "name", all.x = FALSE)
  merged_data <- merge( merged_data, paralog_df, by.x = "paralogue_name", by.y = "gene_name", all.x = FALSE)
  return(merged_data)
}

pairgene1 <- update_paralogue_gene(filtered_df_1before_gene, genepair, filtered_df_1before_para)
pairgene2 <- update_paralogue_gene(filtered_df_1after_gene, genepair, filtered_df_1after_para )
pairgene3 <- update_paralogue_gene(filtered_df_2before_gene, genepair, filtered_df_2before_para )
pairgene4 <- update_paralogue_gene(filtered_df_2after_gene, genepair,filtered_df_2after_para )
