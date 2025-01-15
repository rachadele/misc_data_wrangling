projPath="/space/grp/rschwartz/rschwartz/single_cell_stuff"
source("/space/grp/rschwartz/rschwartz/single_cell_stuff/src/seurat_functions.r")
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
load_required_packages()
library(gplots)
library(argparse)
library(cellxgene.census)
set.seed(123)
topdir="/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/"


add_feature_id <- function(seurat_obj, mapping, gene_column, feature_column) {
  gene_names <- rownames(seurat_obj[["RNA"]])
  # Create a lookup table using the mapping file
  lookup_table <- setNames(mapping[[feature_column]], mapping[[gene_column]])
  # Match Seurat object genes to the lookup table
  feature_ids <- lookup_table[gene_names]
  # Add the feature_id to the Seurat object metadata
  seurat_obj[["RNA"]] <- AddMetaData(seurat_obj[["RNA"]], feature_ids, "feature_id")

  counts <- seurat_obj[["RNA"]]$counts
  meta <- seurat_obj[["RNA"]][[]]
  meta$feature_id <- sub("\\..*$", "", meta$feature_id)
  rownames(counts) <- meta$feature_id
  meta$hgnc <- rownames(meta)
  filtered_meta <- meta[!is.na(meta$feature_id), ]
  rownames(filtered_meta) <- filtered_meta$feature_id
  filtered_counts <- counts[!is.na(rownames(counts)), ]
  filtered_counts <- filtered_counts[rownames(filtered_meta), , drop = FALSE]

  new_counts <- CreateAssayObject(filtered_counts)
  seurat_obj[["RNA"]] <- new_counts
  seurat_obj[["RNA"]][[]] <- filtered_meta
  DefaultAssay(seurat_obj) <- "RNA"
  return(seurat_obj)
}

parse_gtf <- function(gtf_file) {
  # Read the GTF file
  gtf_data <- read.table(gtf_file, sep = "\t", comment.char = "#", header = FALSE, stringsAsFactors = FALSE)
  
  # Extract gene_name and ensembl_id from the 9th column (attributes) and filter for 'gene' entries
  gtf_data <- gtf_data %>%
   mutate(
    feature_id = str_replace(str_extract(V9, 'gene_id [^;]+'), 'gene_id ', ''),
    feature_name = str_replace(str_extract(V9, 'gene_name [^;]+'), 'gene_name ', '')
  ) %>% unique()

  return(gtf_data[, c("feature_id", "feature_name")])  # Return the relevant columns
}


query_files = list.files("/space/grp/rschwartz/rschwartz/single_cell_stuff/rds", pattern="raw",full.names=TRUE)
query_files <- grep(".*subsample*", query_files, value = TRUE, invert=TRUE)

queries <- list()
for (filepath in query_files) {
  data <- readRDS(filepath)
  dataset_name <- strsplit(strsplit(filepath,"\\/")[[1]][8],"_raw")[[1]][1]
  queries[[dataset_name]] <- data
  #output_path=gsub("raw", "raw_", tissue, filepath)
}

#lau_map <- read.table(file.path(topdir,"/lau/GSM4775574_NC7_features.tsv.gz"), sep="\t", header=FALSE)
lim_map <- parse_gtf("/space/grp/rschwartz/rschwartz/census-stuff/meta/genome_refs/gencode.v30.annotation.gtf")
nagy_map <- parse_gtf("/space/grp/rschwartz/rschwartz/census-stuff/meta/genome_refs/GRCh38-1.2.0_premrna/genes/genes.gtf")
rosmap_map <- parse_gtf("/space/grp/rschwartz/rschwartz/census-stuff/meta/genome_refs/gencode.v24.annotation.gtf")

lim <- add_feature_id(seurat_obj=queries$lim, mapping=lim_map, gene_column=2,feature_column=1)
saveRDS(lim, file.path(projPath,"/rds/ensembl_mapped/lim_raw.rds"))
nagy <- add_feature_id(seurat_obj=queries$nagy, mapping=nagy_map, gene_column=2,feature_column=1)
saveRDS(nagy, file.path(projPath,"/rds/ensembl_mapped/nagy_raw.rds")) 
rosmap <- add_feature_id(seurat_obj=queries$rosmap, mapping=rosmap_map, gene_column=2,feature_column=1)
saveRDS(rosmap, file.path(projPath,"/rds/ensembl_mapped/rosmap_raw.rds")) 
