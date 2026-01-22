# library(zellkonverter)
library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(SingleCellExperiment)
library(SeuratDisk)

library(BiocParallel)
register(MulticoreParam(64))

# Get the cell type argument from the command line.
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[1]
cat("Processing cell type:", cell_type, "\n")

####################################
# Load data
# Convert("/work/aliu10/AD_sc_Project/results/adata_GEM.h5ad", dest = "h5seurat", overwrite = TRUE, dense = TRUE)
seurat_obj <- LoadH5Seurat("/work/aliu10/AD_sc_Project/results/adata_GEM.h5seurat", assay = "RNA", meta.data = FALSE)
metadata <- read.csv("/work/aliu10/AD_sc_Project/results/adata_metadata.csv", row.names = 1) # Load metadata CSV
seurat_obj@meta.data = metadata
seurat_obj@assays$RNA$data = seurat_obj@assays$RNA$scale.data # Make sure "data" is normalized data

# Scale data; necessary for computing harmonized module eigengenes.
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

####################################
# Set up Seurat object for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

####################################
# Construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("CellType", "donor_id"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'CellType' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

####################################
# Set up expression data for WGCNA for this cell type.
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cell_type,  # use current cell type as group name
  group.by = 'CellType',
  assay = 'RNA',   # using the RNA assay
  layer = 'data'           # using normalized data
)

# Test different soft powers.
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'   # alternatively "unsigned" or "signed hybrid"
)

###################################################################
# Check soft power
output_dir <- "/work/aliu10/AD_sc_Project/results/hdWGCNA/soft_power/"

saveRDS(
  seurat_obj,
  file = paste0(output_dir, cell_type, "_soft_power.rds")
)

cat("Saved soft thresholding object for:", cell_type, "\n")
###################################################################

# library(patchwork)
# # plot the results:
# plot_list <- PlotSoftPowers(seurat_obj)
# # assemble with patchwork
# wrap_plots(plot_list, ncol=2)

# Define soft power per cell type; obtained from last step
soft_power_list <- list(
  "Astrocyte" = 8,
  "Excitatory" = 6,
  "Inhibitory" = 7,  
  "Microglia-PVM" = 7,  
  "Oligodendrocyte" = 6,    
  "OPC" = 6
)

seurat_obj <- readRDS(paste0("/work/aliu10/AD_sc_Project/results/hdWGCNA/soft_power/", cell_type, "_soft_power.rds"))

# Get the soft power for the current cell type
softpower <- soft_power_list[[cell_type]]

# Construct the co-expression network; TOM is written with the cell type name.
dir.create("/work/aliu10/AD_sc_Project/results/hdWGCNA/TOM/", recursive = TRUE, showWarnings = FALSE)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = softpower,
  tom_outdir = "/work/aliu10/AD_sc_Project/results/hdWGCNA/TOM", 
  tom_name = cell_type,
  overwrite_tom = TRUE
)

# Module Eigengenes and Connectivity
# Compute module eigengenes (MEs) for each donor in the entire snRNA-seq data.
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "donor_id"
)

# Compute eigengene-based connectivity (kME).
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'CellType',
  group_name = cell_type,
  TOM_use = paste0("/work/aliu10/AD_sc_Project/results/hdWGCNA/TOM/", cell_type, "_TOM.rda")
)

# Optionally, rename modules to reflect the cell type.
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(cell_type, ".M")
)

# Save the resulting object for this cell type.
dir.create("/work/aliu10/AD_sc_Project/results/hdWGCNA/hdWGCNA_object/", recursive = TRUE, showWarnings = FALSE)
output_file <- paste0('/work/aliu10/AD_sc_Project/results/hdWGCNA/hdWGCNA_object/hdWGCNA_object_', cell_type, '.rds')
saveRDS(seurat_obj, file = output_file)
cat("Saved:", output_file, "\n")
