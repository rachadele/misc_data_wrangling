projPath="/space/grp/rschwartz/rschwartz/single_cell_stuff"
source("/space/grp/rschwartz/rschwartz/single_cell_stuff/src/seurat_functions.r")
source("/space/grp/rschwartz/rschwartz/single_cell_stuff/src/malat1_function.R") 
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
load_required_packages()
library(gplots)
library(argparse)
library(cellxgene.census)
set.seed(123)

query_files = list.files("/space/grp/rschwartz/rschwartz/single_cell_stuff/rds", pattern="raw",full.names=TRUE)
query_files <- grep(".*subsample*", query_files, value = TRUE, invert=TRUE)

queries <- list()
for (filepath in query_files) {
  data <- readRDS(filepath)
  dataset_name <- strsplit(strsplit(filepath,"\\/")[[1]][8],"_")[[1]][1]
  queries[[dataset_name]] <- data
  #output_path=gsub("raw", "raw_", tissue, filepath)
}

velmeshev <- queries$velmeshev
#velmeshev@meta.data %>% colnames()
#velmeshev$region

lim <- queries$lim
#lim@meta.data %>% colnames()
#lim$Location

# Function to split a Seurat object by a specified column and write to files
split_seurat_and_save <- function(seurat_obj, split_column, output_dir, prefix) {
  # Split the Seurat object by the given column
  split_objs <- SplitObject(seurat_obj, split.by = split_column)
  
  # For each split object, save as an RDS file with the split value in the filename
  lapply(names(split_objs), function(split_val) {
    split_obj <- split_objs[[split_val]]
    file_name <- paste0(output_dir, "/", prefix, "_", split_val, "_raw.rds")
    saveRDS(split_obj, file_name)
    message("Saved: ", file_name)
  })
}

# Example usage:
output_dir <- paste0(projPath,"/rds")
split_seurat_and_save(velmeshev, "region", output_dir, prefix="velmeshev")
split_seurat_and_save(lim, "Location", output_dir, "lim")



query_files = list.files("/space/grp/rschwartz/rschwartz/single_cell_stuff/rds/ensembl_mapped", pattern="raw",full.names=TRUE)
query_files <- grep(".*subsample*", query_files, value = TRUE, invert=TRUE)

queries <- list()
for (filepath in query_files) {
  data <- readRDS(filepath)
  dataset_name <- strsplit(strsplit(filepath,"\\/")[[1]][8],"_")[[1]][1]
  queries[[dataset_name]] <- data
  #output_path=gsub("raw", "raw_", tissue, filepath)
}

# Example usage:
output_dir <- paste0(projPath,"/rds/ensembl_mapped")
split_seurat_and_save(velmeshev, "region", output_dir, prefix="velmeshev")
split_seurat_and_save(lim, "Location", output_dir, "lim")
