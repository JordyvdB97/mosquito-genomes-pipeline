# This R script accompanies the unpublished manuscript:  
# Van der Beek, J. G., Ibáñez-Justicia, A., Biesmeijer, J. C., Lizarazo-Forero, E.,  
# Stroo, A., van de Vossenberg, B. T. L. H., Warbroek, T., & Schrama, M. J. J.  
# (n.d.). The differentiating power of mitochondrial genes: complete mitogenome  
# sequences of 27 mosquito species present in Europe. [Unpublished manuscript].  
#
# This script calculates nucleotide diversity across the mitochondrial genomes  
# of the analyzed mosquito species. The complete mitochondrial genomes were 
# aligned using MAFFT (v.7.490) with default settings.  
#
# All data used in these analyses, including the aligned mitochondrial sequences,  
# can be accessed and downloaded from: https://github.com/JordyvdB97/mosquito-genomes-pipeline/.  
# The script also automates the process of downloading these files directly.  
#
# For questions or further inquiries, please contact:  
# Jordy van der Beek (jordy.vanderbeek@naturalis.nl).

# Install and load necessary packages
# Install and load necessary packages
if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
if (!requireNamespace("PopGenome", quietly = TRUE)) {devtools::install_github("pievos101/PopGenome")}
if (!requireNamespace("tidyverse", quietly = TRUE)) {install.packages("tidyverse")}
if (!requireNamespace("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
if (!requireNamespace("gggenes", quietly = TRUE)) {install.packages("gggenes")}
if (!requireNamespace("cowplot", quietly = TRUE)) {install.packages("cowplot")}
if (!requireNamespace("seqinr", quietly = TRUE)) {install.packages("seqinr")}
if (!requireNamespace("pegas", quietly = TRUE)) {install.packages("pegas")}
if (!requireNamespace("ape", quietly = TRUE)) {install.packages("ape")}
library(PopGenome)
library(tidyverse)
library(ggplot2)
library(gggenes)
library(cowplot)
library(seqinr)
library(pegas)
library(ape)

# Define URLs for GitHub-hosted files
alignments_base_url <- "https://raw.githubusercontent.com/JordyvdB97/mosquito-genomes-pipeline/main/mitogenome_alignments/"
annotations_url <- "https://raw.githubusercontent.com/JordyvdB97/mosquito-genomes-pipeline/main/R_gene_annotations.csv"
alignment_all_genomes_url <- "https://raw.githubusercontent.com/JordyvdB97/mosquito-genomes-pipeline/9c0252307b04fef1bcc9006f5f6c94a8a45837b6/mitogenome_alignments/R_whole_mitogenome_alignment.fasta"

# Files to download
alignment_files <- c(
  "R_whole_mitogenome_alignment_Aedes.fasta",
  "R_whole_mitogenome_alignment_Anopheles.fasta",
  "R_whole_mitogenome_alignment_Culex.fasta",
  "R_whole_mitogenome_alignment_Culiseta.fasta"
)

# Function to read FASTA directly from github
read_fasta_from_url <- function(url) {
  # Create a temporary directory and download the file
  temp_dir <- file.path(tempdir(), "fasta_files")
  dir.create(temp_dir, showWarnings = FALSE)
  temp_file <- file.path(temp_dir, "alignment.fasta")
  # Download the file to the temporary directory
  download.file(url, temp_file, mode = "wb")
  # Read the data from the temporary directory
  data <- PopGenome::readData(temp_dir, format = "fasta")
  # Clean up temporary directory if needed
  return(data)
}

# Set all window parameters
total_mito_size <- 18702 # set chromosome size
window_size <- 200 # set window size
window_jump <- 50 # set jump size
window_start <- seq(from = 1, to = total_mito_size, by = window_jump) # use seq to find the start points of each window
window_stop <- window_start + window_size # add the size of the window to each start point 
window_start <- window_start[which(window_stop < total_mito_size)] # remove windows from the start and stop vectors
window_stop <- window_stop[which(window_stop < total_mito_size)]
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2) # save as a data.frame

# Genera list
genera <- c("Aedes", "Anopheles", "Culex", "Culiseta")

# Initialize nucleotide diversity data frame
pi_mosquitoes <- data.frame(window_mid = windows$mid,
                            stringsAsFactors = FALSE)

for (genus in genera) {
  # Load the FASTA data from GitHub
  url <- paste0(alignments_base_url, "R_whole_mitogenome_alignment_", genus, ".fasta")
  mosquitoes <- read_fasta_from_url(url)
  
  # Generate sliding windows
  mosquitoes_sw <- sliding.window.transform(mosquitoes, width = window_size, jump = window_jump, type = 2)
  # Calculate pi in each window
  mosquitoes_sw <- diversity.stats(mosquitoes_sw, pi = TRUE)
  # Extract diversity stats and divide by window length
  mosquito_nuc_div_sw <- mosquitoes_sw@nuc.diversity.within
  mosquito_nuc_div_sw <- mosquito_nuc_div_sw/window_size
  # Combine into a data.frame and remove the row.names
  mosquito_nd_sw <- data.frame(windows$mid, mosquito_nuc_div_sw, row.names = NULL)
  
  # Add diversity data to the combined data frame
  pi_mosquitoes[[genus]] <- mosquito_nd_sw$pop.1
}

# Download and read the annotations
annotations <- read.csv(annotations_url, na.strings = c("", "NA"))

# Define colors for different gene types
color_map <- c(PCG = "black", rRNA = "black", tRNA = "gray")

# Filter out tRNA genes for labeling
annotations_no_TRNA <- annotations %>%
  filter(Type != "tRNA")

# Replace "reverse" with FALSE and "forward" with TRUE in Direction column
annotations$Direction <- ifelse(annotations$Direction == "reverse", FALSE, TRUE)

# Get the common x-axis limits
x_min <- min(c(min(annotations$Minimum), min(mosquito_nd_sw$windows.mid)))
x_max <- max(c(max(annotations$Maximum), max(mosquito_nd_sw$windows.mid)))

# Create plot 1: Gene Annotations
p1 <- ggplot(annotations, aes(xmin = Minimum, xmax = Maximum, y = 0, fill = Type, forward = Direction)) +
  geom_gene_arrow(aes(xmin = Minimum, xmax = Maximum, fill = Type, color = Type), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_label_repel(
    data = annotations_no_TRNA,
    aes(x = (Minimum + Maximum) / 2, label = Name),
    box.padding = 0.5,
    segment.color = '#cccccc',
    fill = "white",
    nudge_x = 0,
    nudge_y = 0.01,
    size = 7 / .pt,  # Convert to points for ggplot
    fontface = "italic"
  ) +
  theme_void() +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)) +
  xlim(x_min, x_max)

# Create plot 2: Line Plot
p2 <- ggplot(data = pi_mosquitoes, aes(x = window_mid)) +
  annotate("rect", xmin = 2524, xmax = 4029, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
  annotate("rect", xmin = 11158, xmax = 11680, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
  geom_line(aes(y = Aedes, color = "Aedes"), linewidth = 0.5) +
  geom_line(aes(y = Anopheles, color = "Anopheles"), linewidth = 0.5) +
  geom_line(aes(y = Culex, color = "Culex"), linewidth = 0.5) +
  geom_line(aes(y = Culiseta, color = "Culiseta"), linewidth = 0.5) +
  theme_light(base_size = 9) +
  ylab("Nucleotide Diversity (\u03c0)") +
  xlab("Nucleotide Position (bp)") +
  xlim(x_min, x_max) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9),
    legend.text = element_text(size = 9, face = "italic")
  ) +
  labs(color = "")

# Combine the plots
combined_plot <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v", axis = "l", rel_heights = c(2,8))

# Save the plot
#cowplot::ggsave2("figure_3.jpg", combined_plot, width = 165, height = 110, units = "mm", dpi = 300)
#cowplot::ggsave2("figure_3.eps", combined_plot, width = 165, height = 110, units = "mm", dpi = 300, device = cairo_ps)

###### Creating table format #####
# Read the alignment directly from the GitHub URL
temp_file <- tempfile(fileext = ".fasta")
download.file(alignment_all_genomes_url, temp_file, mode = "wb")

# Read genomic data from alignment folder in FASTA format
mosquitoes_all <- seqinr::read.alignment(temp_file, format = "fasta")

# Clean up the temporary file if no longer needed
unlink(temp_file)

# set parameters annotation of the barcoding region
cox1_barcoding_region_annotation <- data.frame(
  Name = "cox1_barcoding_region",
  Minimum = 2553,
  Maximum = 3210,
  Type = "region",
  Length = 658,  # Length = Maximum - Minimum + 1
  Direction = "forward",
  stringsAsFactors = FALSE
)

# Append COX1 annotation to the existing annotations table
annotations_table <- rbind(annotations, cox1_barcoding_region_annotation)

# Define the genera
genera2 <- c("Ae", "An", "Cx", "Cs")  # Use abbreviations for genus names

# Initialize results table
results <- data.frame(
  Markers = character(),
  Length = numeric(),
  Variable_sites_number = numeric(),
  Variable_sites_percentage = numeric(),
  Nucleotide_diversity = numeric(),
  Nucleotide_diversity_Ae = numeric(),
  Nucleotide_diversity_An = numeric(),
  Nucleotide_diversity_Cx = numeric(),
  Nucleotide_diversity_Cs = numeric(),
  stringsAsFactors = FALSE
)

# Process each annotation
for (i in 1:nrow(annotations_table)) {
  # Extract annotation details
  annotation_name <- annotations_table$Name[i]
  start <- annotations_table$Minimum[i]
  end <- annotations_table$Maximum[i]
  length <- end - start + 1
  
  cat("Processing marker:", annotation_name, "\n")
  cat("Start:", start, "End:", end, "Expected Length:", length, "\n")
  
  # Subset alignment for the region
  sub_sequences <- lapply(mosquitoes_all$seq, function(seq) {
    substr(seq, start, end)
  })
  
  # Convert to a matrix of characters
  sequence_matrix <- do.call(rbind, lapply(sub_sequences, strsplit, ""))
  
  # Convert to DNAbin format
  dna_bin <- ape::as.DNAbin(sequence_matrix)
  
  # Calculate variable sites
  variable_sites <- length(ape::seg.sites(dna_bin))
  variable_sites_percentage <- (variable_sites / length) * 100
  
  # Calculate nucleotide diversity for all sequences
  nucleotide_diversity <- pegas::nuc.div(dna_bin)
  
  # Initialize genus-specific nucleotide diversity
  nd_Ae <- 0
  nd_An <- 0
  nd_Cx <- 0
  nd_Cs <- 0
  
  # Filter sequences by genus based on sequence names and calculate ND
  genus_sequences <- mosquitoes_all$nam
  genus_Ae <- grep("Ae", genus_sequences)  # Aedes (Ae)
  genus_An <- grep("An", genus_sequences)  # Anopheles (An)
  genus_Cx <- grep("Cx", genus_sequences)  # Culex (Cx)
  genus_Cs <- grep("Cs", genus_sequences)  # Culiseta (Cs)
  
  # Calculate nucleotide diversity for each genus
  # Calculate nucleotide diversity for each genus
  if (length(genus_Ae) > 0) {
    dna_bin_Ae <- dna_bin[genus_Ae]  # Subset by genus (Aedes)
    nd_Ae <- pegas::nuc.div(dna_bin_Ae)
  }
  
  if (length(genus_An) > 0) {
    dna_bin_An <- dna_bin[genus_An]  # Subset by genus (Anopheles)
    nd_An <- pegas::nuc.div(dna_bin_An)
  }
  
  if (length(genus_Cx) > 0) {
    dna_bin_Cx <- dna_bin[genus_Cx]  # Subset by genus (Culex)
    nd_Cx <- pegas::nuc.div(dna_bin_Cx)
  }
  
  if (length(genus_Cs) > 0) {
    dna_bin_Cs <- dna_bin[genus_Cs]  # Subset by genus (Culiseta)
    nd_Cs <- pegas::nuc.div(dna_bin_Cs)
  }
  
  # Append to results table
  results <- rbind(
    results,
    data.frame(
      Markers = annotation_name,
      Length = length,
      Variable_sites_number = variable_sites,
      Variable_sites_percentage = variable_sites_percentage,
      Nucleotide_diversity = nucleotide_diversity,
      Nucleotide_diversity_Ae = nd_Ae,
      Nucleotide_diversity_An = nd_An,
      Nucleotide_diversity_Cx = nd_Cx,
      Nucleotide_diversity_Cs = nd_Cs
    )
  )
}

# Print the results table
print(results)

# Optionally, save the table to a CSV file
#write.csv(results, file = "table_2.csv", row.names = FALSE)
