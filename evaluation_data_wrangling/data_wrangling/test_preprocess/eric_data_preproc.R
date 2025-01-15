
projPath="/space/grp/rschwartz/rschwartz/single_cell_stuff"
source("/space/grp/rschwartz/rschwartz/single_cell_stuff/src/malat1_function.R")
source("/space/grp/rschwartz/rschwartz/single_cell_stuff/src/seurat_functions.r")
source("/space/grp/rschwartz/rschwartz/goldowitz/single_cell/scripts/sc_functions.R") 
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
load_required_packages()
library(knitr)
library(gplots)
set.seed(123)

add_feature_meta <- function(seurat_obj, mapping, gene_column, feature_column) {
  gene_names <- rownames(seurat_obj[["RNA"]])
  # Create a lookup table using the mapping file
  lookup_table <- setNames(mapping[[feature_column]], mapping[[gene_column]])
  # Match Seurat object genes to the lookup table
  feature_ids <- lookup_table[gene_names]
  # Add the feature_id to the Seurat object metadata
  seurat_obj[["RNA"]] <- AddMetaData(seurat_obj[["RNA"]], feature_ids, "feature_name")
  seurat_obj
}


topdir <- "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data"

#lim
counts <- read.csv(paste0(topdir,"/lim/GSE180928_filtered_cell_counts.csv.gz"), row.names=1)
meta <- read.table(paste0(topdir,"/lim/GSE180928_metadata.csv.gz"),sep=",",header=TRUE, row.names=1)
lim <- CreateSeuratObject(counts=counts, min.cells = 3, min.features = 200)
meta$cellid <- gsub("-", "\\.", rownames(meta))
lim$cellid <- colnames(lim)
newmeta <- left_join(lim@meta.data, meta, by="cellid")
rownames(newmeta) <- newmeta$cellid
lim@meta.data <- newmeta
saveRDS(lim, paste0(projPath,"/rds/lim_raw.rds"))

#nagy
mtx <- ReadMtx(paste0(topdir,"/nagy/GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz"),cells=paste0(topdir,"/nagy/GSE144136_CellNames.csv.gz"),features=paste0(topdir,"/nagy/GSE144136_GeneNames.csv.gz"), feature.column=2, feature.sep=",",cell.sep=",", cell.column=2, skip.cell=1, skip.feature=1)
dim(mtx)
#dim(mtx)
# [1] 30062 78886
nagy <- CreateSeuratObject(counts=mtx, min.cells = 3, min.features = 200)
cell_type <- as.data.frame(do.call(c, lapply(colnames(nagy), function(x) {
    label <- strsplit(x, "\\.")[[1]][1]
    return(label)
})))
rownames(cell_type) <- colnames(nagy)
nagy$cell_type <- cell_type
saveRDS(nagy, paste0(projPath,"/rds/nagy_raw.rds"))

#rosmap
mtx <- ReadMtx(paste0(topdir,"/rosmap/filtred_count_matrix.mtx"),cells=paste0(topdir,"/rosmap/filtered_column_metadata.txt"),features=paste0(topdir,"/rosmap/filtered_gene_row_names.txt"), feature.column=1, skip.cell=1)
meta <- read.table(paste0(topdir,"/rosmap/filtered_column_metadata.txt"),header=TRUE, row.names=1)
dim(mtx)
#dim(mtx)
# 65217 104559
rosmap <- CreateSeuratObject(counts=mtx, min.cells = 3, min.features = 200)
meta <- meta[colnames(rosmap), ]
meta$cellid <- rownames(meta)
rosmap$cellid <- colnames(rosmap)
identical(rownames(meta),colnames(rosmap))
newmeta <- left_join(rosmap@meta.data, meta, by="cellid")
rownames(newmeta) <- newmeta$cellid
newmeta$cellid <- NULL
rosmap@meta.data <- newmeta
saveRDS(rosmap, paste0(projPath,"/rds/rosmap_raw.rds"))



#pineda
#changed for ensembl mapping
pineda_dirs <- list.dirs(file.path(topdir, "/pineda/GSE174332_RAW"), full.names = TRUE, recursive = TRUE)
# Remove the topxlevel directory from the list
pineda_dirs <- pineda_dirs[pineda_dirs != file.path(topdir, "/pineda/GSE174332_RAW")]


load_seurat_data <- function(directory) {
  files_in_dir <- list.files(directory, full.names = TRUE)
  counts_file <- grep("counts_fil\\.mtx$", files_in_dir, value = TRUE)
  cell_metadata_file <- grep("col_metadata\\.tsv$", files_in_dir, value = TRUE)
  feature_metadata_file <- grep("row_metadata\\.tsv$", files_in_dir, value = TRUE)
  
  # Read the count matrix using the provided format
  counts <- ReadMtx(
    counts_file,
    cells = cell_metadata_file,
    features = feature_metadata_file,
    feature.column = 1, # Assumes the feature column is 1
    skip.cell = 1,
    skip.feature=1       # Assumes the first column in cell metadata should be skipped
  )
  # Create a Seurat object using the counts and cell metadata
  seurat_obj <- CreateSeuratObject(counts = counts, min.cells=3, min.features = 200)
  return(seurat_obj)
}

seurat_list <- lapply(pineda_dirs, function(dir) load_seurat_data(dir))
# Merge all Seurat objects into one
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
pineda <- JoinLayers(combined_seurat)
# View the combined object
cellmeta <- list.files(file.path(topdir, "/pineda/GSE174332_RAW"), pattern="col_metadata\\.tsv$", recursive=TRUE, full.names=TRUE)
# Read all files and combine them into a single data frame
combined_cellmeta <- do.call(rbind, lapply(cellmeta, function(x) {
    read.table(x, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    }))
rownames(combined_cellmeta) <- combined_cellmeta$Barcode
combined_cellmeta$Barcode <- NULL
pineda <- AddMetaData(pineda, metadata=combined_cellmeta)

pineda_map <- read.table(file.path(topdir,"pineda/GSE174332_RAW/191114_PN_324_snRNA-F7/row_metadata.tsv"), sep="\t", header=FALSE)
pineda <- add_feature_meta(seurat_obj=pineda,mapping=pineda_map, gene_column=1,feature_column=2)

saveRDS(pineda, file="/space/grp/rschwartz/rschwartz/single_cell_stuff/rds/ensembl_mapped/pineda_raw.rds")

mapping_df <- data.frame(
CellClass = pineda$CellClass,
  SubType = pineda$SubType,
  CellType = pineda$CellType
)
mapdf <- mapping_df %>% distinct()
write.table(mapdf, file=file.path(projPath,"/meta/pineda_mapdf.tsv"),sep="\t",row.names=FALSE)

#lau

library(stringr)
lau <- list()

filenames <- dir(paste0(topdir,"/lau"), pattern="^GSM")
# Extract filename prefixes using sub
samples <- sub("_.*", "", filenames) %>% unique()

for (sample in samples) {
  # Get the list of directories matching the sample prefix
  matching_files <- dir(paste0(topdir, "/lau"), pattern = sample, full.names=TRUE)
  orig.ident <- basename(matching_files)
  group <- strsplit(orig.ident,"_")[[1]][2]
  index <- gsub("^GSM[0-9]+_(?:AD|NC)([0-9]+).*", "\\1", orig.ident) %>% unique()
                         
  
  raw_counts <- ReadMtx(mtx=matching_files[3],features=matching_files[2],cells=matching_files[1]) # , feature.column=1)
  obj <- CreateSeuratObject(counts=raw_counts, min.cells = 3, min.features = 200)
  obj$active.ident=group
  obj$sample=sample
  obj$orig.ident <- gsub("-1",paste0("_",index),colnames(obj))
  lau[[group]] <- obj
}

lau <- merge(lau[[1]], lau[2:length(lau)], add.cell.ids =names(lau), project = "lau")
lau <- JoinLayers(lau)


meta <- read_xlsx(paste0(topdir,"/lau/Lau_et_al_PNAS_AD_snRNAseq_celltype_subcluster_ID.xlsx"))
meta.data <- left_join(lau@meta.data, meta, by=c("orig.ident"="ID"), keep=TRUE)
rownames(meta.data) <- rownames(lau@meta.data)
lau@meta.data <- meta.data

lau <- SetIdent(lau, value=lau$active.ident)
lau <- lau[ ,lau$Cell_type!="NA"]
#lau_map <- read.table(file.path(topdir,"/lau/GSM4775574_NC7_features.tsv.gz"), sep="\t", header=FALSE)
##lau <- add_feature_meta(seurat_obj=lau,mapping=lau_map, gene_column=1,feature_column=2)
saveRDS(lau,"./rds/lau_raw.rds")

#velmeshev
mtx <- ReadMtx(paste0(topdir,"/velmeshev/matrix.mtx"),cells=paste0(topdir,"/velmeshev/barcodes.tsv"),features=paste0(topdir,"/velmeshev/genes.tsv"))#, feature.column=1)
meta <- read.table(paste0(topdir,"/velmeshev/meta.txt"),sep="\t",header=TRUE)
dim(mtx)
#dim(mtx)
# 65217 104559
velmeshev <- CreateSeuratObject(counts=mtx, min.cells = 3, min.features = 200)
velmeshev@meta.data$cell <- rownames(velmeshev@meta.data)
newmeta <- left_join(velmeshev@meta.data, meta, by="cell") 
rownames(newmeta) <- colnames(velmeshev)
velmeshev@meta.data <- newmeta
Idents(velmeshev) <- velmeshev$sample
genemeta <- read.table(paste0(topdir,"/velmeshev/genes.tsv"), sep='\t',header=FALSE)
colnames(genemeta) <- c("feature_id","feature_name")
#feature_mapping <- genemeta[match(rownames(velmeshev), genemeta$feature_id),]
#velmeshev@assays$RNA[[]]$feature_name <- feature_mapping$feature_name
saveRDS(velmeshev, paste0(projPath,"/rds/velmeshev_raw.rds"))