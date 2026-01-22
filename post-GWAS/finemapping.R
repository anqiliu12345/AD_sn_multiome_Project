library(susieR)
library(dplyr)
library(tidyselect)
library(stringr)
library(GenomicRanges)
library(data.table)

set.seed(1)

window_size <- 1500000  # ±1.5Mb
name_group <- "female_AD_vs_CT"
compare_groups <- c("female_AD_vs_CT_GRN", "female_AD_vs_CT_GEM")

# #####################
## Genes of interest
cell_types <- c("Astrocyte", "Microglia-PVM", "Excitatory", "Inhibitory", "OPC", "Oligodendrocyte")  

for (compare_group in compare_groups) {
    # Initialize empty data frame
    gene_sources <- data.frame(Gene = character(), CellType = character(), Category = character(), stringsAsFactors = FALSE)

    genes <- c()

    for (cell_type in cell_types) {
        top_genes_path <- paste0("/work/aliu10/AD_sc_Project/results/FlowSig/", compare_group, "/", cell_type, "/Top_Genes_per_GEM.csv")
        in_out_path <- paste0("/work/aliu10/AD_sc_Project/results/FlowSig/", compare_group, "/", cell_type, "/inflow_gem_outflow_paths.csv")

        if (file.exists(top_genes_path) && length(readLines(top_genes_path, n = 2)) >= 2 &&
        file.exists(in_out_path) && length(readLines(in_out_path, n = 2)) >= 2) {
            top_genes <- read.csv(top_genes_path)
            in_out_genes <- read.csv(in_out_path) # Ligands and receptors involved in pathways

            top_gene_names <- unique(unlist(top_genes))
            gene_sources <- rbind(gene_sources, data.frame(Gene = top_gene_names, CellType = cell_type, Category = "GEM"))

            in_out_names <- unique(c(in_out_genes$Inflow, in_out_genes$Outflow))
            gene_sources <- rbind(gene_sources, data.frame(Gene = in_out_names, CellType = cell_type, Category = "InOut"))

            merged_genes <- unique(unlist(strsplit(c(top_gene_names, in_out_names), split = "\\+")))
            
            if (str_detect(compare_group, "GRN")) {
                # For TF-gene networks only; Extract TFs
                genes_extracted <- str_extract(colnames(top_genes), "^[A-Z0-9]+(?=_direct)")
                gene_sources <- rbind(gene_sources, data.frame(Gene = genes_extracted, CellType = cell_type, Category = "TF"))
                genes <- unique(c(genes, merged_genes, genes_extracted))
            } else {genes <- unique(c(genes, merged_genes))}
                
        } else {
            warning(paste("Missing file for", cell_type))
        }
    }

    write.csv(
      gene_sources,
      file = paste0("/work/aliu10/AD_sc_Project/results/FlowSig/", compare_group, "/genes_sources.csv"),
      row.names = FALSE
    )
} 

# Load gene source data from both comparison groups
compare_groups1 <- compare_groups[1]
compare_groups2 <- compare_groups[2]
gene_sources1 <- read.csv(paste0("/work/aliu10/AD_sc_Project/results/FlowSig/", compare_groups1, "/genes_sources.csv"))
gene_sources2 <- read.csv(paste0("/work/aliu10/AD_sc_Project/results/FlowSig/", compare_groups2, "/genes_sources.csv"))

# Combine both tables
combined_sources <- rbind(gene_sources1, gene_sources2)

# Identify genes that appear in both sources
shared_genes <- intersect(gene_sources1$Gene, gene_sources2$Gene)

# Filter the combined data to retain only shared genes
shared_sources <- combined_sources[combined_sources$Gene %in% shared_genes, ]

# Count how many times each Gene-CellType-Category combination appears
# If a combination appears twice, it means it exists in both sources
shared_summary <- shared_sources %>%
  group_by(Gene, CellType, Category) %>%
  tally(name = "Count") %>%
  filter(Count == 2)

# Keep only those shared in the "GEM" category
shared_summary <- shared_summary[shared_summary$Category == 'GEM', ]

# Write result to CSV
dir.create(paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group), showWarnings = FALSE, recursive = TRUE)
write.csv(shared_summary, paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/shared_genes_GEM_GRN.csv"), row.names = FALSE)

genes = shared_summary$Gene
# #####################
## eQTL summary statistics & extract gene regions
eqtl_genes <- read.table("/work/aliu10/AD_sc_Project/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eGenes.txt.gz", header = TRUE, sep = "\t")

eqtl_data <- eqtl_genes %>%
  mutate(chr = gsub("^chr", "", chr),
         SNP = paste0(chr, ":", variant_pos))

# eqtl_subset <- eqtl_subset[!duplicated(eqtl_subset$rs_id_dbSNP155_GRCh38p13), ]
eqtl_data$z_score <-eqtl_data$slope / eqtl_data$slope_se
eqtl_subset = eqtl_data[eqtl_data$gene_name %in% genes,]
eqtl_subset = eqtl_subset[eqtl_subset$gene_chr != "chrX",]


# #####################
## AD GWAS summary statistics
gwas5 <- fread("/work/aliu10/AD_sc_Project/susieR/AD_GWAS/GCST90027158.tsv.gz")
gwas5$SNP <- paste0(gwas5[[3]], ":", gwas5[[4]]) 
gwas5$z_score <-gwas5$beta / gwas5$standard_error

#####################
## Create ± windows around gene regions
gr <- GRanges(
  seqnames = eqtl_subset$gene_chr,
  ranges = IRanges(
    start = pmax(eqtl_subset$gene_start - window_size, 1),
    end = eqtl_subset$gene_end + window_size
  ),
  variant_id = eqtl_subset$rs_id_dbSNP155_GRCh38p13,
  gene = eqtl_subset$gene_name
)

#####################
# Convert eQTL/GWAS SNPs to GRanges
eqtl_snps_gr <- GRanges(
  seqnames = paste0("chr", eqtl_data$chr),
  ranges = IRanges(start = eqtl_data$variant_pos, end = eqtl_data$variant_pos),
  SNP = eqtl_data$SNP
)

gwas_snps_gr <- GRanges(
  seqnames = paste0("chr", gwas5[[3]]),
  ranges = IRanges(start = gwas5[[4]], end = gwas5[[4]]),
  SNP = gwas5$SNP
)

# Overlapping SNPs list
overlap_snp_list <- list()

for (i in seq_along(gr)) {
  gene_region <- gr[i]
  gene_name <- gene_region$gene
  
  # Get overlapping SNPs
  eqtl_hits <- subsetByOverlaps(eqtl_snps_gr, gene_region)
  gwas_hits <- subsetByOverlaps(gwas_snps_gr, gene_region)

  # Take intersection of SNPs
  shared_snps <- intersect(eqtl_hits$SNP, gwas_hits$SNP)
  
  # Store by gene name
  overlap_snp_list[[gene_name]] <- shared_snps
}

#####################
plink_path <- "/work/aliu10/conda/envs/plink/bin/plink"
plink2_path <- "/work/aliu10/conda/envs/plink2/bin/plink2"
plink_prefix_base <- "/work/aliu10/AD_sc_Project/susieR/1000Genome/Plink_files/"
ld_out_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/LD_matrix_1.5M")
dir.create(ld_out_dir, showWarnings = FALSE, recursive = TRUE)

# ==== Extract genotype data and compute LD for each gene region ====
for (i in seq_along(gr)) {
  chr <- gsub("chr", "", as.character(seqnames(gr[i])))
  region_start <- start(gr[i])
  region_end <- end(gr[i])
  gene <- gr$gene[i]
    
  snps <- overlap_snp_list[[gene]]
  if (length(snps) < 2) {
    message("Skipping ", gene, ": too few SNPs.")
    next
  }

  # Write SNP list to temporary file
  snp_list_file <- paste0(ld_out_dir, "/", gene, "_overlap_snps.txt")
  writeLines(snps, con = snp_list_file)
    
  bfile <- paste0(plink_prefix_base, "chr", chr, "_EUR")
  out_prefix <- file.path(ld_out_dir, paste0(gene, "_chr", chr, "_", region_start, "_", region_end))

  # Step 1: Use plink2 to extract genotype data for the overlaping SNPs in the gene region
  plink2_cmd <- paste(
    plink2_path,
    "--bfile", shQuote(bfile),
    "--chr", chr,
    "--extract", shQuote(snp_list_file),
    "--maf", 0.01, 
    "--set-missing-var-ids", "'@:#'",
    # "--new-id-max-allele-len", "23", "missing",
    "--make-bed",
    "--out", shQuote(out_prefix)
  )
  system(plink2_cmd)

  # Step 2: Use plink1.9 to compute pairwise LD matrix
  plink_ld_cmd <- paste(
    plink_path,
    "--bfile", shQuote(out_prefix),
    "--r square",
    "--out", shQuote(out_prefix)
  )
  system(plink_ld_cmd)

  # Step 3: Extract SNP ID list
  snp_list_cmd <- paste("cut -f2", shQuote(paste0(out_prefix, ".bim")), ">", shQuote(paste0(out_prefix, ".snp_list")))
  system(snp_list_cmd)

  message(gene, " done!")
}
##########################################################################################################################

## Fine-mapping with susieR using LD blocks
# Directory where LD matrix and SNP list files are stored
ld_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/LD_matrix_1.5M")
ld_files <- list.files(ld_dir, pattern = "\\.ld$", full.names = TRUE)
dir.create(file.path("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "susie", "susie_obj_eQTL"), recursive = TRUE)

# Iterate through each gene region and its corresponding LD file
for (ld_file in ld_files){
  
  # Derive the SNP list file name
  snp_file <- gsub(".ld$", ".snp_list", ld_file)
  
  # Read SNP list
  snps <- fread(snp_file, header = FALSE)[[1]]
    
  # Subset eQTL and LD matrix by SNP list (chr:position)
  eqtl_subset_sub <- eqtl_data[eqtl_data$SNP %in% snps, ]
  match_idx <- match(snps, eqtl_subset_sub$SNP)
  valid <- !is.na(match_idx)
  
  eqtl_subset_matched <- eqtl_subset_sub[match_idx[valid], ]

  R_full <- as.matrix(read.table(ld_file))
  R <- R_full[valid, valid]
    
  if (is.null(dim(R)) || any(dim(R) < 2)) {
      message("Skipping ", gene, ": LD matrix too small or invalid.")
      next
  }
    
  colnames(R) = eqtl_subset_matched$rs_id_dbSNP155_GRCh38p13
  rownames(R) = eqtl_subset_matched$rs_id_dbSNP155_GRCh38p13
    
  susie_fit = susie_rss(z = eqtl_subset_matched$z_score, R = R, n = 943, L = 10, niter = 100)
    
  # Save results
  gene <- gsub("\\.ld$", "", basename(ld_file))
  saveRDS(susie_fit, file = paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/susie/susie_obj_eQTL/", gene, ".rds"))
  
  message("Done: ", gene)
} 


## Fine-mapping with susieR using LD blocks
# Directory where LD matrix and SNP list files are stored
ld_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/LD_matrix_1.5M")
ld_files <- list.files(ld_dir, pattern = "\\.ld$", full.names = TRUE)
dir.create(file.path("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "susie", "susie_obj_GWAS"), recursive = TRUE)

# Iterate through each gene region and its corresponding LD file
for (ld_file in ld_files){
  
  # Derive the SNP list file name
  snp_file <- gsub(".ld$", ".snp_list", ld_file)
  
  # Read SNP list
  snps <- fread(snp_file, header = FALSE)[[1]]
    
  # Subset GWAS and LD matrix by SNP list (chr:position)
  gwas_sub <- gwas5[gwas5$SNP %in% snps, ]
  match_idx <- match(snps, gwas_sub$SNP)
  valid <- !is.na(match_idx)
  
  gwas_matched <- gwas_sub[match_idx[valid], ]

  R_full <- as.matrix(read.table(ld_file))
  R <- R_full[valid, valid]
    
  if (is.null(dim(R)) || any(dim(R) < 2)) {
      message("Skipping ", gene, ": LD matrix too small or invalid.")
      next
  }

  colnames(R) = gwas_matched$variant_id
  rownames(R) = gwas_matched$variant_id
    
  susie_fit = susie_rss(z = gwas_matched$z_score, R = R, n = 487511, L = 10, niter = 100)
    
  # Save results
  gene <- gsub("\\.ld$", "", basename(ld_file))
  saveRDS(susie_fit, file = paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", name_group, "/susie/susie_obj_GWAS/", gene, ".rds"))
  
  message("Done: ", gene)
} 
