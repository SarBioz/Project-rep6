Class_0_s <- Class_0_s[, "paralogs", drop = FALSE]
Class_0_s$sarvenaz_class <- "Class_0"
Class_1_s <- Class_1_s[, "paralogs", drop = FALSE]
Class_1_s$sarvenaz_class <- "Class_1"
Class_2_s <- Class_2_s[, "paralogs", drop = FALSE]
Class_2_s$sarvenaz_class <- "Class_2"
Class_3_s <- Class_3_s[, "paralogs", drop = FALSE]
Class_3_s$sarvenaz_class <- "Class_3"
Class_4_s <- Class_4_s[, "paralogs", drop = FALSE]
Class_4_s$sarvenaz_class <- "Class_4"
combined_df <- rbind(Class_0_s, Class_1_s, Class_2_s, Class_3_s, Class_4_s)
colnames(combined_df)[colnames(combined_df) == "paralogs"] <- "genes"
merged_df <- merge(combined_df, chris_genes, by = "genes", all.y = TRUE)
write.csv(merged_df, '/project/pi_frederic_chain_uml_edu/scRNAseq/Sarvenaz/geneslist/ClassCombinations_merged_V.csv', row.names = FALSE)


if ("sarvenaz_class" %in% colnames(merged_df) & "cells" %in% colnames(merged_df)) {
  # Create a column to check if the classifications match
  merged_df <- merged_df %>% mutate(Concordant = sarvenaz_class == cells)
  
  # Calculate percentage of matching classifications
  concordance_percentage <- mean(merged_df$Concordant, na.rm = TRUE) * 100
  print(paste("Concordance Percentage: ", round(concordance_percentage, 2), "%", sep = ""))
} else {
  print("Column names do not match.")
}
install.packages("irr")
library(irr)
# Calculate Cohen's kappa
kappa_score <- kappa2(merged_df[,c("sarvenaz_class", "cells")])
print(paste("Cohen's Kappa Score: ", round(kappa_score$value, 2), sep = ""))


library(ggplot2)
library(reshape2)
# Create a contingency table
contingency_table <- table(merged_df$sarvenaz_class, merged_df$cells)

# Convert to data frame for ggplot
contingency_df <- as.data.frame(as.table(contingency_table))
colnames(contingency_df) <- c("sarvenaz_class", "cells", "Frequency")

# Plot heatmap
ggplot(data = contingency_df, aes(x = sarvenaz_class, y = cells, fill = Frequency)) +
  geom_tile() +
  geom_text(aes(label = Frequency), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap of Gene Paralog Classifications Concordance",
       x = "Study 2 Classes (sarvenaz_class)",
       y = "Study 1 Classes (cells)")
