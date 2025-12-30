library(tidyr)
library(tidyverse)
library(dplyr)
library(vcfR)
library(lsa)
library(ComplexHeatmap)

##########

scriptdir <- dirname(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))))

dir.create(path = paste0(scriptdir, '/Data/10x/Demultiplexing/02_VCF_Consensus'), showWarnings = F)

outdir <- paste0(scriptdir, '/Data/10x/Demultiplexing/02_VCF_Consensus')

########## Set up paths and variables ##########

souporcell_vcf <- list()
souporcell_chromR <- list()
souporcell_vcfR_filtered <- list()
souporcell_genotypes <- list()

souporcell_path <- '/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Demultiplexing/01_Souporcell'

list.files(souporcell_path)

chromium_lanes <- c("Run_S1081_L364A", "Run_S1146_L364B")

for (lane in chromium_lanes) {
  vcf_file <- file.path(souporcell_path, lane, 'cluster_genotypes.vcf')
  souporcell_vcf[[lane]] <- read.vcfR(vcf_file)
  souporcell_chromR[[lane]] <- vcfR::create.chromR(souporcell_vcf[[lane]])
  souporcell_vcfR_filtered[[lane]] <- chromR2vcfR(souporcell_chromR[[lane]], use.mask = TRUE)
  
  souporcell_genotypes[[lane]] <- as_tibble(
    extract.gt(element = "GT", souporcell_vcfR_filtered[[lane]], IDtoRowNames = FALSE)
  )
  
  # Identify genotype columns (excluding 'ID')
  gt_cols <- colnames(souporcell_genotypes[[lane]])
  
  if (!all(colSums(is.na(souporcell_genotypes[[lane]][, gt_cols])) == nrow(souporcell_genotypes[[lane]]))) {
    message("Found GT genotype format in query vcf. Will use that metric for query correlation.")
    
    # Detect separator.00000000000000000000000000000000000
    first_row <- souporcell_genotypes[[lane]][1, gt_cols]
    separator <- ifelse(any(grepl("\\|", first_row)), "|", 
                        ifelse(any(grepl("/", first_row)), "/", NA))
    
    if (is.na(separator)) {
      message("Can't identify a separator for the GT field in query vcf, moving on to using GP.")
    } else {
      message(paste("Detected", separator, "separator for GT genotype format in query vcf"))
      souporcell_genotypes[[lane]][, gt_cols] <- lapply(
        souporcell_genotypes[[lane]][, gt_cols],
        function(x) {
          x <- gsub(paste0("0", separator, "0"), "0", x)
          x <- gsub(paste0("0", separator, "1"), "1", x)
          x <- gsub(paste0("1", separator, "0"), "1", x)
          x <- gsub(paste0("1", separator, "1"), "2", x)
          as.numeric(x)
        }
      )
    }
  }
  
  # Add 'chr' prefix if missing
  if (!any(grepl('chr', souporcell_vcfR_filtered[[lane]]@fix[,'CHROM']))) {
    souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'] <- paste0('chr', souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'])
  }
  
  # Create ID column for SNPs
  souporcell_genotypes[[lane]]$ID <- if (all(is.na(souporcell_vcfR_filtered[[lane]]@fix[,'REF'])) & all(is.na(souporcell_vcfR_filtered[[lane]]@fix[,'ALT']))) {
    message("The REF and ALT categories are not provided for the reference and/or the query vcf. Will use just the chromosome and position to match SNPs.")
    paste0(souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'], ":", souporcell_vcfR_filtered[[lane]]@fix[,'POS'])
  } else if (all(is.na(souporcell_vcfR_filtered[[lane]]@fix[,'REF']))) {
    message("The REF categories are not provided for the reference and/or the query vcf. Will use the chromosome, position and ALT to match SNPs.")
    paste0(souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'], ":", souporcell_vcfR_filtered[[lane]]@fix[,'POS'], "_", souporcell_vcfR_filtered[[lane]]@fix[,'REF'])
  } else if (all(is.na(souporcell_vcfR_filtered[[lane]]@fix[,'ALT']))) {
    message("The ALT categories are not provided for the reference and/or the query vcf. Will use the chromosome, position and REF to match SNPs.")
    paste0(souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'], ":", souporcell_vcfR_filtered[[lane]]@fix[,'POS'], "_", souporcell_vcfR_filtered[[lane]]@fix[,'ALT'])
  } else {
    message("Found REF and ALT in both query and reference genotype vcfs. Will use chromosome, position, REF and ALT to match SNPs.")
    paste0(souporcell_vcfR_filtered[[lane]]@fix[,'CHROM'], ":", souporcell_vcfR_filtered[[lane]]@fix[,'POS'], "_", souporcell_vcfR_filtered[[lane]]@fix[,'REF'], "_", souporcell_vcfR_filtered[[lane]]@fix[,'ALT'])
  }
  
  # Rename columns (except ID)
  ncols <- ncol(souporcell_genotypes[[lane]])
  colnames(souporcell_genotypes[[lane]])[1:(ncols-1)] <- paste('Lane', lane, colnames(souporcell_genotypes[[lane]])[1:(ncols-1)], sep = '_')
  
  # Remove rows with NA IDs
  souporcell_genotypes[[lane]] <- souporcell_genotypes[[lane]][!str_detect(souporcell_genotypes[[lane]]$ID, "NA:NA_NA_NA"), ]
}

merged_data_souporcell <- Reduce(function(x, y) right_join(x, y, by="ID"), souporcell_genotypes) %>%
  select(ID, everything())

################################################################################

pearson_correlation <- function(df, ref_df, query_df) {
  for (col in 1:ncol(df)) {
    for (row in 1:nrow(df)) {
      query_row <- as.numeric(query_df[, row][[1]])
      ref_col <- as.numeric(ref_df[, col][[1]])
      query_sd <- sd(query_row, na.rm = TRUE)
      ref_sd <- sd(ref_col, na.rm = TRUE)
      if (all(is.na(query_row)) || all(is.na(ref_col)) || 
          is.na(query_sd) || query_sd == 0 || 
          is.na(ref_sd) || ref_sd == 0 || 
          sum(complete.cases(query_row, ref_col)) == 0) {
        df[row, col] <- 0
      } else {
        df[row, col] <- cor(query_row, ref_col, method = "pearson", use = "complete.obs")
      }
    }
  }
  return(df)
}

########## Correlate all the query genotypes with the individuals genotyped ##########

merged_data_souporcell_noid <- merged_data_souporcell[, colnames(merged_data_souporcell)]
merged_data_souporcell_noid <- merged_data_souporcell_noid[, setdiff(colnames(merged_data_souporcell), 'ID')]

pearson_correlations <- as.data.frame(matrix(nrow = ncol(merged_data_souporcell_noid), ncol = ncol(merged_data_souporcell_noid)))
colnames(pearson_correlations) <- colnames(merged_data_souporcell_noid)
rownames(pearson_correlations) <- colnames(merged_data_souporcell_noid)
pearson_correlations <- pearson_correlation(pearson_correlations, merged_data_souporcell_noid, merged_data_souporcell_noid)

S1081_L364A_samples <- grep("S1081_L364A", rownames(pearson_correlations), value = TRUE)
S1146_L364B_samples <- grep("S1146_L364B", rownames(pearson_correlations), value = TRUE)

# Create subset: S1081 as rows, S1146 as columns
pearson_correlations <- pearson_correlations[S1081_L364A_samples, S1146_L364B_samples]

pearson_with_names <- cbind(Samples = rownames(pearson_correlations), 
                            pearson_correlations)

write.table(pearson_with_names, 
            file = paste0(outdir, "/Souporcell_vs_Souporcell_clusters_SetA__Pearson_correlations.tsv"), 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")

########## Create correlation figures ##########
col_fun <- colorRampPalette(c("white", "red"))(10)
pPearsonCorrelations <- Heatmap(as.matrix(pearson_correlations), cluster_rows = T, col = col_fun)

########## Save the correlation figures ##########
Cairo::CairoPNG(filename = paste0(outdir, "/Souporcell_vs_Souporcell_clusters_SetA__Pearson_correlations.png"), width = 1500, height = 1500)
ComplexHeatmap::draw(pPearsonCorrelations)
dev.off()
