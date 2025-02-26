# Load necessary libraries
library(sleuth)
library(dplyr)

# Read in the sample table
stab <- read.table("sample_table.txt", header = TRUE)

# Initialize Sleuth object using sleuth_prep function
so <- sleuth_prep(stab)

# Fit a model comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')

# Fit the reduced model (null model)
so <- sleuth_fit(so, ~1, 'reduced')

# Perform the likelihood ratio test for differential expression
so <- sleuth_lrt(so, 'reduced', 'full')

# Extract the results of the likelihood ratio test
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# Filter results for significant transcripts (FDR/qval < 0.05)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% 
  dplyr::arrange(pval)

# Print top 10 significant transcripts
head(sleuth_significant, n = 10)

# Write results for significant transcripts (FDR < 0.05) to a separate file
write.table(sleuth_significant, file = "sleuthfdr_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Need bootstrap results from kallisto for plotting
so = sleuth_prep(stab, extra_bootstrap_summary = TRUE,read_bootstrap_tpm = TRUE)
# Select the most significant transcript (smallest p-value)
top_target <- sleuth_significant$target_id[1]  # Take the first row (lowest pval)

# Check if a significant transcript was found
if (!is.na(top_target)) {
  
  # Generate the bootstrap expression plot
  topplot <- plot_bootstrap(so, top_target, units = "tpm", color_by = "condition")
  print(topplot)  # Display in RStudio
  
  # Save the plot dynamically with the transcript ID
  png_filename <- paste0(top_target, "_TPM.png")
  png(png_filename)
  print(topplot)
  dev.off()
  
  # Print confirmation message
  print(paste("Plot saved as:", png_filename))
  
} else {
  print("No significant transcripts found (qval < 0.05).")}

