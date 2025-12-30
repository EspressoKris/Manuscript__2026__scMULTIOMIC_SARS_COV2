library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(scCustomize)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(hdf5r)

# Get the script directory
projectdir <- dirname(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))))

dir.create(file.path(projectdir, "Data/10x/Per_Sample_Demultiplexed_Mapped_Data"), showWarnings = F)

destination_dir <- file.path(projectdir, "Data/10x/Per_Sample_Demultiplexed_Mapped_Data")

# Define the per_sample_outs directory
Run_S1081_L364A_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1081_L364A/outs/per_sample_outs")
Run_S1146_L364B_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364B/outs/per_sample_outs")
Run_S1146_L364C_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364C/outs/per_sample_outs")
Run_S1146_L364D_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364D/outs/per_sample_outs")

# Get all sample directories
Run_S1081_L364A_samples <- basename(list.dirs(Run_S1081_L364A_dir, full.names = TRUE, recursive = FALSE))
Run_S1146_L364B_samples <- basename(list.dirs(Run_S1146_L364B_dir, full.names = TRUE, recursive = FALSE))
Run_S1146_L364C_samples <- basename(list.dirs(Run_S1146_L364C_dir, full.names = TRUE, recursive = FALSE))
Run_S1146_L364D_samples <- basename(list.dirs(Run_S1146_L364D_dir, full.names = TRUE, recursive = FALSE))

# Function to extract patient ID (remove _ASB, _GB, _Ambiguous)
get_patient_id <- function(sample_name) {
  gsub("_(ASB|GB|Ambiguous)$", "", sample_name)
}

# Remove suffixes from sample names
Run_S1081_L364A_samples <- unique(sapply(Run_S1081_L364A_samples, get_patient_id))
Run_S1146_L364B_samples <- unique(sapply(Run_S1146_L364B_samples, get_patient_id))
Run_S1146_L364C_samples <- unique(sapply(Run_S1146_L364C_samples, get_patient_id))
Run_S1146_L364D_samples <- unique(sapply(Run_S1146_L364D_samples, get_patient_id))

# Print cleaned sample names
cat("Run_S1081_L364A samples:", paste(Run_S1081_L364A_samples, collapse = ", "), "\n")
cat("Run_S1146_L364B samples:", paste(Run_S1146_L364B_samples, collapse = ", "), "\n")
cat("Run_S1146_L364C samples:", paste(Run_S1146_L364C_samples, collapse = ", "), "\n")
cat("Run_S1146_L364D samples:", paste(Run_S1146_L364D_samples, collapse = ", "), "\n")

# Load demultiplexed barcodes
Demultiplexed_barcodes <- readRDS(file.path(projectdir, 'Data/10x/Demultiplexing/03_Barcodes_Assignment/corrected_sample_barcodes.rds'))

# Create list of run sample vectors with their names
run_samples_list <- list(
  Run_S1081_L364A = Run_S1081_L364A_samples,
  Run_S1146_L364B = Run_S1146_L364B_samples,
  Run_S1146_L364C = Run_S1146_L364C_samples,
  Run_S1146_L364D = Run_S1146_L364D_samples
)

# Create filtered barcode lists for each run
Demultiplexed_barcodes_by_run <- list()

for (run_name in names(run_samples_list)) {
  samples <- run_samples_list[[run_name]]
  
  # Create pattern to match any of the samples in this run
  pattern <- paste0("^(", paste(samples, collapse = "|"), ")_")
  
  # Filter Demultiplexed_barcodes
  filtered_barcodes <- Demultiplexed_barcodes[str_detect(names(Demultiplexed_barcodes), pattern)]
  
  # Store in list
  Demultiplexed_barcodes_by_run[[run_name]] <- filtered_barcodes
  
  cat("Created filtered barcodes for", run_name, "with", length(filtered_barcodes), "samples\n")
  cat("  Samples:", paste(names(filtered_barcodes), collapse = ", "), "\n")
}

# Define the per_sample_outs directory of non HTO demultiplexed mapped data
run_dirs <- list(
  Run_S1081_L364A = file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1081_L364A_No_HTO/outs/per_sample_outs/Run_S1081_L364A_No_HTO/count"),
  Run_S1146_L364B = file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364B_No_HTO/outs/per_sample_outs/Run_S1146_L364B_No_HTO/count"),
  Run_S1146_L364C = file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364C_No_HTO/outs/per_sample_outs/Run_S1146_L364C_No_HTO/count"),
  Run_S1146_L364D = file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364D_No_HTO/outs/per_sample_outs/Run_S1146_L364D_No_HTO/count")
)

# Process each run
for (run_name in names(run_dirs)) {
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Processing", run_name, "\n")
  cat(rep("=", 80), "\n", sep = "")
  
  run_dir <- run_dirs[[run_name]]
  
  # Check if directory exists
  if (!dir.exists(run_dir)) {
    warning("Directory not found: ", run_dir)
    next
  }
  
  # Load count matrix
  h5_file <- file.path(run_dir, 'sample_filtered_feature_bc_matrix.h5')
  if (!file.exists(h5_file)) {
    warning("H5 file not found: ", h5_file)
    next
  }
  
  cat("Loading count matrix from:", h5_file, "\n")
  count_mat <- scCustomize::Read_CellBender_h5_Mat(h5_file)
  
  # Identify GEX and ADT features
  all_features <- rownames(count_mat)
  gex_features <- all_features[1:38606]
  adt_features <- all_features[38607:length(all_features)]
  all_barcodes <- colnames(count_mat)
  
  cat("Total features:", length(all_features), "\n")
  cat("GEX features:", length(gex_features), "\n")
  cat("ADT features:", length(adt_features), "\n")
  cat("Total barcodes:", length(all_barcodes), "\n")
  
  # Get filtered barcodes for this run
  run_filtered_barcodes <- Demultiplexed_barcodes_by_run[[run_name]]
  
  # Process each sample in this run
  for (sample_name in names(run_filtered_barcodes)) {
    
    cat("\n  Processing sample:", sample_name, "\n")
    
    # Get barcodes for this sample
    sample_barcodes <- run_filtered_barcodes[[sample_name]]
    
    # Find common barcodes between sample and count matrix
    common_barcodes <- intersect(sample_barcodes, all_barcodes)
    
    if (length(common_barcodes) == 0) {
      warning("    No common barcodes found for ", sample_name)
      next
    }
    
    cat("    Found", length(common_barcodes), "barcodes\n")
    
    # Create sample directory
    sample_dir <- file.path(destination_dir, sample_name)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Remove old files if they exist
    for (file in c('features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz', 'filtered_feature_bc_matrix.h5')) {
      fpath <- file.path(sample_dir, file)
      if (file.exists(fpath)) file.remove(fpath)
    }
    
    # Subset count matrix for this sample (includes both GEX and ADT)
    sample_count_mat <- count_mat[all_features, common_barcodes]
    
    cat("    Writing combined GEX+ADT data...\n")
    
    # Write features (both GEX and ADT)
    features_df <- data.frame(
      feature_id = all_features,
      feature_name = all_features,
      feature_type = c(rep("Gene Expression", length(gex_features)), 
                       rep("Antibody Capture", length(adt_features))),
      stringsAsFactors = FALSE
    )
    write.table(features_df, file.path(sample_dir, 'features.tsv'), 
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    system2("gzip", file.path(sample_dir, 'features.tsv'))
    
    # Write barcodes
    gz_barcodes <- gzfile(file.path(sample_dir, 'barcodes.tsv.gz'), "w")
    writeLines(common_barcodes, gz_barcodes)
    close(gz_barcodes)
    
    # Write matrix
    mtx_path <- file.path(sample_dir, 'matrix.mtx')
    Matrix::writeMM(sample_count_mat, mtx_path)
    system2("gzip", mtx_path)
    
    cat("    Creating 10X H5 file...\n")
    
    # Create 10X H5 file using DropletUtils
    tryCatch({
      h5_output <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
      
      # Use DropletUtils to write H5
      DropletUtils::write10xCounts(
        path = h5_output,
        x = sample_count_mat,
        gene.id = features_df$feature_id,
        gene.symbol = features_df$feature_name,
        gene.type = features_df$feature_type,
        barcodes = common_barcodes,
        overwrite = TRUE,
        type = "HDF5",
        genome = "GRCh38",
        version = "3"
      )
      
      cat("    Successfully created H5 file for", sample_name, "\n")
    }, error = function(e) {
      cat("    Failed to create H5 file for", sample_name, ": ", e$message, "\n")
    })
    
    cat("    Completed", sample_name, "\n")
  }
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("Processing complete!\n")
cat(rep("=", 80), "\n", sep = "")