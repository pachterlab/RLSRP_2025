#!/usr/bin/env Rscript

# Function to check and install missing packages
install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        options(repos = c(CRAN = "https://cloud.r-project.org"))  # Set CRAN mirror
        install.packages(pkg)
    }
}

# Check & install CRAN packages
install_if_missing("dplyr")
install_if_missing("jsonlite")

# Install missing Bioconductor package (depmap)
if (!requireNamespace("depmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")  # Ensure BiocManager installs correctly
    }
    BiocManager::install("depmap")
}

# Load the libraries
library(dplyr)
library(jsonlite)
library(depmap)


# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript update_ccle_metadata.R input_json output_json")
}

input_json <- args[1]
output_json <- args[2]

# Retrieve the DepMap metadata
metadata <- depmap_metadata()

# Load the JSON file as a list and convert to a data frame
json_data <- fromJSON(input_json, simplifyDataFrame = TRUE)
json_df <- as.data.frame(json_data)

# Loop through each row in the JSON data
for (i in 1:nrow(json_df)) {
  sample_title <- json_df$sample_title[i]
  
  # Filter the DepMap metadata for the given cell line
  sample_data <- metadata %>% filter(cell_line == sample_title)
  
  # Extract primary_disease and other metadata
  disease_info <- sample_data %>% select(primary_disease, subtype_disease, sex, age, lineage, lineage_subtype, Cellosaurus_NCIt_disease)
  
  # Update json_df with new metadata
  json_df$primary_disease[i] <- disease_info$primary_disease[1]
  json_df$subtype_disease[i] <- disease_info$subtype_disease[1]
  json_df$sex[i] <- disease_info$sex[1]
  json_df$age[i] <- disease_info$age[1]
  json_df$lineage[i] <- disease_info$lineage[1]
  json_df$lineage_subtype[i] <- disease_info$lineage_subtype[1]
  json_df$Cellosaurus_NCIt_disease[i] <- disease_info$Cellosaurus_NCIt_disease[1]
}

# Convert the updated data frame back to JSON
updated_json <- toJSON(json_df, pretty = TRUE, auto_unbox = TRUE)

# Save to the specified output file
write(updated_json, output_json)
