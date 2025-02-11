# This R script accompanies the unpublished manuscript: 
# Van der Beek, J. G., Ibáñez-Justicia, A., Biesmeijer, J. C., Lizarazo-Forero, E., 
# Stroo, A., van de Vossenberg, B. T. L. H., Warbroek, T., & Schrama, M. J. J. 
# (n.d.). The differentiating power of mitochondrial genes: complete mitogenome 
# sequences of 27 mosquito species present in Europe. [Unpublished manuscript].  
#
# The script provides a comprehensive guide to the analysis and visualizations 
# conducted on pairwise distance matrices for all mitochondrial genes, rRNAs, 
# and tRNAs. These matrices were generated using Geneious Prime (v.2023.2.1) 
# based on MAFFT alignments (v.7.490) with default settings. 
#
# All distance matrices used in the analyses are available for download at:
# https://github.com/JordyvdB97/mosquito-genomes-pipeline/. This script also 
# automates the process of downloading the required files directly. 
#
# For questions or further inquiries, please contact: 
# Jordy van der Beek (jordy.vanderbeek@naturalis.nl).


# Install and load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
if (!requireNamespace("gridExtra", quietly = TRUE)) {install.packages("gridExtra")}
if (!requireNamespace("httr", quietly = TRUE)) {install.packages("httr")}
if (!requireNamespace("jsonlite", quietly = TRUE)) {install.packages("jsonlite")}
if (!requireNamespace("dplyr", quietly = TRUE)) {install.packages("dplyr")}
if (!requireNamespace("tidyr", quietly = TRUE)) {install.packages("tidyr")}
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(httr)
library(jsonlite)

# Function to find the smallest common classification level between two species
find_common_classification <- function(species1, species2) {
  metadata1 <- metadata[metadata$number == species1, ]
  metadata2 <- metadata[metadata$number == species2, ]
  
  # If metadata not found for either species, return NA
  if (nrow(metadata1) == 0 || nrow(metadata2) == 0) {
    return(NA)
  }
  
  # Compare at species level
  if (metadata1$species == metadata2$species) {
    return("species")
  }
  
  # Compare at subgenus level
  if (!is.na(metadata1$subgenus) && !is.na(metadata2$subgenus) && metadata1$subgenus == metadata2$subgenus) {
    return("subgenus")
  }
  
  # Compare at genus level
  if (metadata1$genus == metadata2$genus) {
    return("genus")
  }
  
  # Compare at tribe level
  if (!is.na(metadata1$tribe) && !is.na(metadata2$tribe) && metadata1$tribe == metadata2$tribe) {
    return("tribe")
  }
  
  # Compare at subfamily level
  if (metadata1$subfamily == metadata2$subfamily) {
    return("subfamily")
  }
  
  # Compare at family level
  if (metadata1$family == metadata2$family) {
    return("family")
  }
  
  # If no common classification level found, return NA
  return(NA)
}

# Define the URLs for the files and the GitHub API for the distance matrices folder
metadata_url <- "https://raw.githubusercontent.com/JordyvdB97/mosquito-genomes-pipeline/61505e3ac2b8f38d1ff85b5697476a99d5f54e2b/R_metadata_species.csv"
github_api_url <- "https://api.github.com/repos/JordyvdB97/mosquito-genomes-pipeline/contents/distance_matrices"
distance_matrix_base_url <- "https://raw.githubusercontent.com/JordyvdB97/mosquito-genomes-pipeline/main/distance_matrices/"


# Fetch the list of files in the distance_matrices folder using the GitHub API
response <- GET(github_api_url)
if (status_code(response) != 200) {
  stop("Failed to retrieve file list from GitHub. Check the repository or internet connection.")
}
file_info <- content(response, as = "parsed", type = "application/json")
file_names <- sapply(file_info, function(x) x$name)

# Download and read the metadata
metadata <- read.csv(metadata_url, na.strings = c("", "NA"))

# Create a dataframe to store plots and their names
plot_list <- c()

# Create an empty dataframe to store the combined similarity data
combined_data <- data.frame(spec_1 = character(3321),
                            spec_2 = character(3321),
                            stringsAsFactors = FALSE)

# Iterate over each file
for (file_name in file_names) {
  file_url <- paste0(distance_matrix_base_url, file_name)
  
  # Download and read the distance matrix
  similarity_matrix <- read.csv(file_url, row.names = 1)
  
  # Define the gene name
  gene_name <- gsub(".csv", "", basename(file_name))
  
  # Convert the distance matrix to a dataframe
  similarity_df <- as.data.frame(similarity_matrix)
  
  # Get row and column names
  row_names <- rownames(similarity_df)
  col_names <- colnames(similarity_df)
  
  # Create an empty dataframe to store results
  similarity_df <- data.frame(spec_1 = character(),
                              spec_2 = character(),
                              similarity = numeric(),
                              stringsAsFactors = FALSE)
  
  # Loop through the upper triangle of the distance matrix
  for (i in 1:(length(row_names) - 1)) {
    for (j in (i + 1):length(col_names)) {
      # Extract registration number and similarity score
      spec1 <- substr(row_names[i], 1, 16)
      spec2 <- substr(col_names[j], 1, 16)
      similarity <- similarity_matrix[i, j]
      
      # Add to the dataframe
      similarity_df <- rbind(similarity_df, data.frame(spec_1 = spec1, spec_2 = spec2, similarity = similarity))
    }
  }
  
  # Add similarity data to dataframe
  combined_data[[gene_name]] <- similarity_df$similarity
  combined_data$spec_1 <- similarity_df$spec_1
  combined_data$spec_2 <- similarity_df$spec_2
  
  # Add the smallest common classification level to the similarity_df
  similarity_df$common_level <- mapply(find_common_classification, similarity_df$spec_1, similarity_df$spec_2)
  
  # Specify the order of levels
  level_order <- c("species", "subgenus", "genus", "subfamily", "family")
  
  # Convert common_level to a factor with specified order
  similarity_df$common_level <- factor(similarity_df$common_level, levels = level_order)
  
  # Compute the range of data for 'species' and 'subgenus'
  species_range <- range(similarity_df$similarity[similarity_df$common_level == "species"])
  subgenus_range <- range(similarity_df$similarity[similarity_df$common_level == "subgenus"])
  
  # Check for overlap
  overlap <- !((species_range[2] < subgenus_range[1]) | (species_range[1] > subgenus_range[2]))
  
  # Create the main ggplot object
  main_plot <- ggplot(similarity_df, aes(x = common_level, y = similarity)) +
    geom_boxplot(size = 0.25, outlier.size = 1) +  # Set boxplot line size and reduce outlier point size
    labs(x = "Taxonomic level", y = "Similarity (%)", title = paste(gene_name)) +
    scale_y_continuous(limits = c(70, 100)) +
    scale_x_discrete(labels = c("spec.", "subg.", "gen.", "subf.", "fam.")) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 9),
      plot.title = element_text(size = 9, face = "bold.italic")  
    )
  
  # Add geom_hline if there's no overlap
  if (!overlap) {
    main_plot <- main_plot +
      geom_hline(
        yintercept = mean(c(species_range[1], subgenus_range[2])),
        linetype = "dashed",
        color = "red",
        size = 0.5
      )
  }
  
  # Save the arranged plot
  assign(paste("plot_", gene_name, sep = ""), main_plot)
  
  # Store plot and its name in the list
  plot_list <- append(plot_list, paste("plot_", gene_name, sep = ""))
}

# Print a list of plots that have been generated
print(paste(plot_list, collapse = ", "))

# Arrange plots of regular genes in a grid
arranged_PCGs_plot <- gridExtra::grid.arrange(plot_atp6, plot_atp8, plot_cox1, plot_cox2, plot_cox3, plot_cytb, plot_nad1, plot_nad2, plot_nad3, plot_nad4, plot_nad4L, plot_nad5, plot_nad6, plot_rrnL, plot_rrnS, ncol = 5)

# Save the plot
#ggplot2::ggsave("figure_2.jpg", arranged_PCGs_plot, width = 165, height = 110, units = "mm", dpi = 300)
#ggplot2::ggsave("figure_2.eps", arranged_PCGs_plot, width = 165, height = 110, units = "mm", dpi = 300, device = cairo_ps)

# Arrange plots of regular genes in a grid
arranged_tRNA_plot <- gridExtra::grid.arrange(plot_trnA, plot_trnC, plot_trnD, plot_trnE, plot_trnF, plot_trnG, plot_trnH, plot_trnI, plot_trnK, plot_trnL1, plot_trnL2, plot_trnM, plot_trnN, plot_trnP, plot_trnQ, plot_trnR, plot_trnS1, plot_trnS2, plot_trnT, plot_trnV, plot_trnW, plot_trnY, ncol = 5)

# Save the plot
#ggplot2::ggsave("supplementary_figure_S1.jpg", arranged_tRNA_plot, width = 165, height = 190, units = "mm", dpi = 300)
#ggplot2::ggsave("supplementary_figure_S1.eps", arranged_tRNA_plot, width = 165, height = 190, units = "mm", dpi = 300, device = cairo_ps)

###### constructing tables ######
# Add the taxonomic classification a pairwise comparison has in common
combined_data$common_level <- mapply(find_common_classification, combined_data$spec_1, combined_data$spec_2)

# Melt the data to long format
gene_data_long <- combined_data %>%
  pivot_longer(cols = -c(common_level, spec_1, spec_2), names_to = "Gene", values_to = "Similarity")

# Group by common_level and Gene, calculate average, min, and max
gene_summary <- gene_data_long %>%
  group_by(common_level, Gene) %>%
  summarize(
    Average = mean(Similarity),
    Min = min(Similarity),
    Max = max(Similarity)
  ) %>%
  mutate(Average_Range = paste(round(Average, 2), " (", round(Min, 2), " - ", round(Max, 2), ")", sep = "")) %>%
  select(-Min, -Max, -Average)

# Pivot the summary data to wide format for creating the table
gene_summary_wide <- tidyr::pivot_wider(gene_summary, names_from = Gene, values_from = Average_Range)

# Reorder the table
custom_order <- c("species", "subgenus", "genus", "subfamily", "family")
gene_summary_wide$common_level <- factor(gene_summary_wide$common_level, 
                                         levels = custom_order)
gene_summary_wide <- gene_summary_wide[order(gene_summary_wide$common_level), ]

# Save the data frame as a CSV file
#write.csv(gene_summary_wide, "supplementary_table_S2.csv", row.names = FALSE)

######

# Filter the combined data for only the species-level comparisons
species_level_similarity <- gene_data_long[gene_data_long$common_level == "species", ]

# Add species identification by joining species_level_similarity with metadata based on spec_1 and number
species_level_similarity$species <- paste(
  metadata$genus[match(species_level_similarity$spec_1, metadata$number)],
  metadata$species[match(species_level_similarity$spec_1, metadata$number)]
)

# Calculate the number of rows for each genus and species based on metadata
species_counts <- metadata %>%
  group_by(genus, species) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(genus_species = paste(genus, species))  # Create a new column combining genus and species

# Group by species and Gene, calculate average, min, and max
species_level_similarity_gene_summary <- species_level_similarity %>%
  group_by(species, Gene) %>%
  summarize(
    Average = mean(Similarity),
    Min = min(Similarity),
    Max = max(Similarity),
    .groups = "drop"
  ) %>%
  mutate(Average_Range = paste(round(Average, 2), " (", round(Min, 2), " - ", round(Max, 2), ")", sep = "")) %>%
  select(-Min, -Max, -Average)

# Add genus_species (genus + species) and species counts (n) from metadata to the species name in the summary table
species_level_similarity_gene_summary <- species_level_similarity_gene_summary %>%
  left_join(species_counts, by = c("species" = "genus_species")) %>%
  mutate(species = paste(species, "(n =", n, ")")) %>%
  select(-n)  # Optionally, remove the count column if it's no longer needed

# Pivot the summary data to wide format for creating the table
species_level_similarity_gene_summary_wide <- tidyr::pivot_wider(
  species_level_similarity_gene_summary,
  names_from = Gene,
  values_from = Average_Range
)

# Save the data frame as a CSV file
#write.csv(species_level_similarity_gene_summary_wide, "supplementary_table_S3.csv", row.names = FALSE)
