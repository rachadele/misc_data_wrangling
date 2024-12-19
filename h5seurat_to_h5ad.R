
library(reticulate)
library(Seurat)
library(SeuratDisk)

get_directory_path <- function() {
  message("Please enter the directory path containing h5seurat files:")
  path <- readline(prompt = "Directory path: ")
  return(path)
}

# Check for command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# If no directory path is provided, prompt the user interactively
if (length(args) < 1) {
  parent_dir <- get_directory_path()
} else {
  parent_dir <- args[1]
}

# Get file paths
filepaths <- list.files(path = parent_dir, full.names = TRUE, pattern = "h5seurat")

# Convert each file
for (filepath in filepaths) {
  Convert(filepath, dest = "h5ad")
}