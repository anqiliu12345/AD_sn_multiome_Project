library(coloc)

compare_groups <- "male_AD_vs_CT"

threshold <- 0.5

GWAS_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", compare_groups, "/susie/susie_obj_GWAS")
eQTL_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", compare_groups, "/susie/susie_obj_eQTL")
GWAS_rds_files <- list.files(GWAS_dir, pattern = "\\.rds$", full.names = TRUE)
eQTL_rds_files <- list.files(eQTL_dir, pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(GWAS_rds_files) == length(eQTL_rds_files))

output_dir <- paste0("/work/aliu10/AD_sc_Project/results/Post_GWAS/", compare_groups, "/coloc/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(GWAS_rds_files)) {
  gwas_data <- readRDS(GWAS_rds_files[i])
  eqtl_data <- readRDS(eQTL_rds_files[i])
  gene <- gsub("\\.rds$", "", basename(GWAS_rds_files[i]))
  
  susie.res <- coloc.susie(gwas_data, eqtl_data)
  
  if (!is.null(susie.res$summary) &&
      is.data.frame(susie.res$summary) &&
      "PP.H4.abf" %in% names(susie.res$summary)) {
    
    coloc_results <- susie.res$summary[susie.res$summary$PP.H4.abf > threshold, ]
    
    if (nrow(coloc_results) > 0) {
      write.csv(
        coloc_results,
        file = file.path(output_dir, paste0(gene, "_coloc_results.csv")),
        row.names = FALSE
      )
      message("Saved coloc results for: ", gene)
    } else {
      message("Skipping ", gene, ": no PP.H4.abf > ", threshold)
    }
  } else {
    message("Skipping ", gene, ": susie.res$summary is NULL or malformed")
  }
}
