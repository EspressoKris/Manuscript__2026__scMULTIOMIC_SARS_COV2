library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Get the script directory
projectdir <- dirname(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))))

# Define the per_sample_outs directory
Run_S1081_L364A_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1081_L364A/outs/per_sample_outs")
Run_S1146_L364B_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364B/outs/per_sample_outs")
Run_S1146_L364C_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364C/outs/per_sample_outs")
Run_S1146_L364D_dir <- file.path(projectdir, "Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364D/outs/per_sample_outs")

# Get all sample directories
Run_S1081_L364A_samples <- list.dirs(Run_S1081_L364A_dir, full.names = TRUE, recursive = FALSE)
Run_S1146_L364B_samples <- list.dirs(Run_S1146_L364B_dir, full.names = TRUE, recursive = FALSE)
Run_S1146_L364C_samples <- list.dirs(Run_S1146_L364C_dir, full.names = TRUE, recursive = FALSE)
Run_S1146_L364D_samples <- list.dirs(Run_S1146_L364D_dir, full.names = TRUE, recursive = FALSE)

# Initialize list to store barcodes
Samples_barcodes <- list()

# Loop through each sample directory
for (Run_dir in c(Run_S1081_L364A_samples, Run_S1146_L364B_samples,
                  Run_S1146_L364C_samples, Run_S1146_L364D_samples)) {
  # Get sample name
  Sample_name <- basename(Run_dir)
  
  # Path to filtered barcodes
  barcode_file <- file.path(Run_dir, "count", "sample_filtered_feature_bc_matrix", "barcodes.tsv.gz")
  
  # Check if file exists
  if (file.exists(barcode_file)) {
    # Read barcodes from gzipped file
    barcodes <- read.table(gzfile(barcode_file), header = FALSE, stringsAsFactors = FALSE)$V1
    
    # Store in list
    Samples_barcodes[[Sample_name]] <- barcodes
    
    cat("Loaded", length(barcodes), "barcodes for", Sample_name, "\n")
  } else {
    warning("Barcode file not found for ", Sample_name)
  }
}

# View summary
cat("\nTotal samples:", length(Samples_barcodes), "\n")

# Define the demultiplexed directories
Run_S1081_L364A_snp_demultiplex_dir <- file.path(projectdir, "Data/10x/Demultiplexing/01_Souporcell/Run_S1081_L364A")
Run_S1146_L364B_snp_demultiplex_dir <- file.path(projectdir, "Data/10x/Demultiplexing/01_Souporcell/Run_S1146_L364B")
Run_S1146_L364C_snp_demultiplex_dir <- file.path(projectdir, "Data/10x/Demultiplexing/01_Souporcell/Run_S1146_L364C")
Run_S1146_L364D_snp_demultiplex_dir <- file.path(projectdir, "Data/10x/Demultiplexing/01_Souporcell/Run_S1146_L364D")

# Create mapping lists
run_to_snp_dir <- list(
  "Run_S1081_L364A" = Run_S1081_L364A_snp_demultiplex_dir,
  "Run_S1146_L364B" = Run_S1146_L364B_snp_demultiplex_dir,
  "Run_S1146_L364C" = Run_S1146_L364C_snp_demultiplex_dir,
  "Run_S1146_L364D" = Run_S1146_L364D_snp_demultiplex_dir
)

run_to_samples <- list(
  "Run_S1081_L364A" = Run_S1081_L364A_samples,
  "Run_S1146_L364B" = Run_S1146_L364B_samples,
  "Run_S1146_L364C" = Run_S1146_L364C_samples,
  "Run_S1146_L364D" = Run_S1146_L364D_samples
)

# Function to extract patient ID (remove _ASB, _GB, _Ambiguous)
get_patient_id <- function(sample_name) {
  gsub("_(ASB|GB|Ambiguous)$", "", sample_name)
}

# Initialize results storage
all_results_detailed <- list()
all_results_corrected <- list()
cluster_assignments <- list()

# SYSTEMATIC ANALYSIS BY RUN
for (run_name in names(run_to_snp_dir)) {
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("DEMULTIPLEX RUN:", run_name, "\n")
  cat(rep("=", 80), "\n", sep = "")
  
  snp_dir <- run_to_snp_dir[[run_name]]
  samples <- run_to_samples[[run_name]]
  
  # Read Souporcell clusters (ALL cells)
  clusters_file <- file.path(snp_dir, "clusters.tsv")
  if (!file.exists(clusters_file)) {
    warning("clusters.tsv not found for ", run_name)
    next
  }
  
  SNPs_demultiplex <- read.table(clusters_file, row.names = 'barcode', header = TRUE, stringsAsFactors = FALSE)
  
  # Rename assignments with "/" to "doublet"
  SNPs_demultiplex$assignment <- ifelse(grepl("/", SNPs_demultiplex$assignment), 
                                        "doublet", 
                                        SNPs_demultiplex$assignment)
  
  # Get unique clusters (excluding doublet)
  clusters <- sort(unique(SNPs_demultiplex$assignment))
  clusters <- clusters[clusters != "doublet"]
  
  cat("\nNumber of clusters identified:", length(clusters), "\n")
  cat("Clusters:", paste(clusters, collapse = ", "), "\n")
  
  # Get counts of doublets (ignore unassigned status, use assignment column)
  n_doublets <- sum(SNPs_demultiplex$assignment == "doublet")
  
  cat("Doublets:", n_doublets, "\n")
  
  # Initialize data frame for this run
  run_results_detailed <- data.frame()
  run_results_corrected <- data.frame()
  
  # Track doublets to be moved to Doublet cluster
  total_doublets_from_clusters <- 0
  
  # Analyze each cluster
  cat("\n", rep("-", 80), "\n", sep = "")
  for (cluster in clusters) {
    cat("\nCLUSTER", cluster, ":\n")
    
    # Get barcodes in this cluster
    cluster_barcodes <- rownames(SNPs_demultiplex[SNPs_demultiplex$assignment == cluster, ])
    total_cluster_cells <- length(cluster_barcodes)
    
    cat("  Total cells in cluster:", total_cluster_cells, "\n")
    cat("  Overlapping barcodes by sample:\n")
    
    # Check overlap with each sample - aggregate by patient ID for assignment
    patient_counts <- list()
    assigned_cells <- 0  # Track cells assigned to samples
    sample_details <- list()  # Store details for each sample
    
    for (sample_dir in samples) {
      sample_name <- basename(sample_dir)
      patient_id <- get_patient_id(sample_name)
      sample_barcodes <- Samples_barcodes[[sample_name]]
      
      if (!is.null(sample_barcodes)) {
        # Find common barcodes
        common_barcodes <- intersect(cluster_barcodes, sample_barcodes)
        count <- length(common_barcodes)
        
        # Store detailed counts (with _ASB/_GB)
        if (count > 0) {
          assigned_cells <- assigned_cells + count
          percentage <- (count / total_cluster_cells) * 100
          
          # Store for detailed plotting
          run_results_detailed <- rbind(run_results_detailed, data.frame(
            Run = run_name,
            Cluster = cluster,
            Sample = sample_name,
            Count = count,
            Percentage = percentage,
            stringsAsFactors = FALSE
          ))
          
          # Store sample details
          sample_details[[sample_name]] <- list(patient_id = patient_id, count = count)
        }
        
        # Aggregate counts by patient ID for assignment
        if (count > 0) {
          if (is.null(patient_counts[[patient_id]])) {
            patient_counts[[patient_id]] <- 0
          }
          patient_counts[[patient_id]] <- patient_counts[[patient_id]] + count
        }
      }
    }
    
    # Print aggregated counts by patient
    for (patient_id in names(patient_counts)) {
      count <- patient_counts[[patient_id]]
      percentage <- (count / total_cluster_cells) * 100
      cat(sprintf("    %s: %d cells (%.2f%%)\n", patient_id, count, percentage))
    }
    
    # Calculate unassigned cells (cells in cluster not matched to any sample)
    n_cluster_unassigned <- total_cluster_cells - assigned_cells
    
    if (n_cluster_unassigned > 0) {
      pct_unassigned <- (n_cluster_unassigned / total_cluster_cells) * 100
      cat(sprintf("    Unassigned: %d cells (%.2f%%)\n", n_cluster_unassigned, pct_unassigned))
      
      run_results_detailed <- rbind(run_results_detailed, data.frame(
        Run = run_name,
        Cluster = cluster,
        Sample = "Unassigned",
        Count = n_cluster_unassigned,
        Percentage = pct_unassigned,
        stringsAsFactors = FALSE
      ))
    }
    
    # Determine assignment (patient with most cells, excluding unassigned)
    if (length(patient_counts) > 0) {
      max_patient <- names(which.max(unlist(patient_counts)))
      max_count <- max(unlist(patient_counts))
      max_percentage <- (max_count / total_cluster_cells) * 100
      
      cat(sprintf("  → Assignment: Souporcell Cluster %s is %s (%.2f%%, %d cells)\n", 
                  cluster, max_patient, max_percentage, max_count))
      
      # Store assignment
      if (is.null(cluster_assignments[[run_name]])) {
        cluster_assignments[[run_name]] <- list()
      }
      cluster_assignments[[run_name]][[as.character(cluster)]] <- list(
        sample = max_patient,
        percentage = max_percentage,
        count = max_count
      )
      
      # CREATE CORRECTED DATASET
      # Keep cells from assigned patient's samples
      # Move cells from other patients to Doublet cluster
      # Mark unassigned cells as Ambiguous
      
      corrected_assigned_count <- 0
      corrected_doublet_count <- 0
      corrected_cluster_samples <- list()
      
      for (sample_name in names(sample_details)) {
        detail <- sample_details[[sample_name]]
        
        if (detail$patient_id == max_patient) {
          # Keep this sample - store for later percentage calculation
          corrected_assigned_count <- corrected_assigned_count + detail$count
          corrected_cluster_samples[[sample_name]] <- detail$count
        } else {
          # Move to Doublet cluster
          corrected_doublet_count <- corrected_doublet_count + detail$count
        }
      }
      
      # Track total doublets to add to Doublet cluster
      total_doublets_from_clusters <- total_doublets_from_clusters + corrected_doublet_count
      
      # Calculate new total for this cluster (after removing doublets)
      corrected_cluster_total <- corrected_assigned_count + n_cluster_unassigned
      
      # Add corrected samples with recalculated percentages
      for (sample_name in names(corrected_cluster_samples)) {
        count <- corrected_cluster_samples[[sample_name]]
        run_results_corrected <- rbind(run_results_corrected, data.frame(
          Run = run_name,
          Cluster = cluster,
          Sample = sample_name,
          Count = count,
          Percentage = (count / corrected_cluster_total) * 100,
          stringsAsFactors = FALSE
        ))
      }
      
      # Rename unassigned to Ambiguous with assigned patient ID (recalculated percentage)
      if (n_cluster_unassigned > 0) {
        ambiguous_sample <- paste0(max_patient, "_Ambiguous")
        run_results_corrected <- rbind(run_results_corrected, data.frame(
          Run = run_name,
          Cluster = cluster,
          Sample = ambiguous_sample,
          Count = n_cluster_unassigned,
          Percentage = (n_cluster_unassigned / corrected_cluster_total) * 100,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # NOW HANDLE DOUBLETS AS A SEPARATE "CLUSTER"
  doublet_barcodes <- rownames(SNPs_demultiplex[SNPs_demultiplex$assignment == "doublet", ])
  if (length(doublet_barcodes) > 0) {
    cat("\nDOUBLET CELLS:\n")
    cat("  Total doublet cells:", length(doublet_barcodes), "\n")
    
    # Check which samples these doublets overlap with
    doublet_assigned <- 0
    for (sample_dir in samples) {
      sample_name <- basename(sample_dir)
      sample_barcodes <- Samples_barcodes[[sample_name]]
      
      if (!is.null(sample_barcodes)) {
        common_barcodes <- intersect(doublet_barcodes, sample_barcodes)
        count <- length(common_barcodes)
        
        if (count > 0) {
          doublet_assigned <- doublet_assigned + count
          percentage <- (count / length(doublet_barcodes)) * 100
          
          run_results_detailed <- rbind(run_results_detailed, data.frame(
            Run = run_name,
            Cluster = "Doublet",
            Sample = sample_name,
            Count = count,
            Percentage = percentage,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Unassigned doublets
    doublet_unassigned <- length(doublet_barcodes) - doublet_assigned
    if (doublet_unassigned > 0) {
      pct_unassigned <- (doublet_unassigned / length(doublet_barcodes)) * 100
      run_results_detailed <- rbind(run_results_detailed, data.frame(
        Run = run_name,
        Cluster = "Doublet",
        Sample = "Unassigned",
        Count = doublet_unassigned,
        Percentage = pct_unassigned,
        stringsAsFactors = FALSE
      ))
    }
    
    # For corrected dataset: Doublet cluster = souporcell doublets + doublets from other clusters
    total_doublet_cells <- length(doublet_barcodes) + total_doublets_from_clusters
    
    run_results_corrected <- rbind(run_results_corrected, data.frame(
      Run = run_name,
      Cluster = "Doublet",
      Sample = "Doublet",
      Count = total_doublet_cells,
      Percentage = 100,
      stringsAsFactors = FALSE
    ))
  }
  
  # Store results for this run
  all_results_detailed[[run_name]] <- run_results_detailed
  all_results_corrected[[run_name]] <- run_results_corrected
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("\nFINAL ASSIGNMENTS FOR", run_name, ":\n")
  cat(rep("-", 80), "\n", sep = "")
  for (cluster in names(cluster_assignments[[run_name]])) {
    assignment <- cluster_assignments[[run_name]][[cluster]]
    cat(sprintf("Souporcell Cluster %s is %s based on %.1f%% overlap (%d cells)\n",
                cluster, assignment$sample, assignment$percentage, assignment$count))
  }
}

# CREATE STACKED BARPLOTS
cat("\n", rep("=", 80), "\n", sep = "")
cat("GENERATING PLOTS\n")
cat(rep("=", 80), "\n", sep = "")

# Create output directory
outdir <- file.path(projectdir, "Data/10x/Demultiplexing/03_Barcodes_Assignment")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Plot for each run (ORIGINAL - detailed version with _ASB/_GB separate)
for (run_name in names(all_results_detailed)) {
  
  run_data <- all_results_detailed[[run_name]]
  
  if (nrow(run_data) == 0) next
  
  # Set cluster order (numeric clusters first, then "Doublet")
  cluster_levels <- c(sort(unique(run_data$Cluster[run_data$Cluster != "Doublet"])), "Doublet")
  run_data$Cluster <- factor(run_data$Cluster, levels = cluster_levels)
  
  # Set factor levels to control stacking order (patients first, then Unassigned)
  sample_levels <- unique(run_data$Sample)
  sample_levels <- c(setdiff(sample_levels, "Unassigned"), "Unassigned")
  run_data$Sample <- factor(run_data$Sample, levels = sample_levels)
  
  # Create color palette
  n_samples <- length(setdiff(sample_levels, "Unassigned"))
  sample_colors <- setNames(
    colorRampPalette(brewer.pal(min(9, max(3, n_samples)), "Set1"))(n_samples),
    setdiff(sample_levels, "Unassigned")
  )
  
  all_colors <- c(sample_colors, "Unassigned" = "gray90")
  
  # Create stacked barplot
  p <- ggplot(run_data, aes(x = Cluster, y = Percentage, fill = Sample)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_fill_manual(values = all_colors) +
    labs(
      title = paste("Souporcell Cluster Composition -", run_name),
      x = "Souporcell Cluster",
      y = "Proportion (%)",
      fill = "Sample"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # Add assignment labels on top (only for numeric clusters, not doublet)
  assignments <- cluster_assignments[[run_name]]
  if (length(assignments) > 0) {
    label_data <- data.frame(
      Cluster = factor(names(assignments), levels = cluster_levels),
      Label = sapply(assignments, function(x) x$sample),
      y = 102
    )
    
    p <- p + geom_text(data = label_data, 
                       aes(x = Cluster, y = y, label = Label, fill = NULL),
                       inherit.aes = FALSE, fontface = "bold", size = 3.5)
  }
  
  # Save plot
  ggsave(file.path(outdir, paste0(run_name, "_cluster_composition.png")),
         plot = p, width = 10, height = 6, dpi = 300)
  ggsave(file.path(outdir, paste0(run_name, "_cluster_composition.pdf")),
         plot = p, width = 10, height = 6)
  
  cat("Saved plot for", run_name, "\n")
}

# Plot for each run (CORRECTED VERSION)
for (run_name in names(all_results_corrected)) {
  
  run_data <- all_results_corrected[[run_name]]
  
  if (nrow(run_data) == 0) next
  
  # Set cluster order (numeric clusters first, then "Doublet")
  cluster_levels <- c(sort(unique(run_data$Cluster[run_data$Cluster != "Doublet"])), "Doublet")
  run_data$Cluster <- factor(run_data$Cluster, levels = cluster_levels)
  
  # Set factor levels to control stacking order (patients first, then Ambiguous, then Doublet)
  sample_levels <- unique(run_data$Sample)
  ambiguous_samples <- grep("_Ambiguous$", sample_levels, value = TRUE)
  regular_samples <- setdiff(sample_levels, c(ambiguous_samples, "Doublet"))
  sample_levels <- c(regular_samples, ambiguous_samples, "Doublet")
  run_data$Sample <- factor(run_data$Sample, levels = sample_levels)
  
  # Create color palette
  n_regular <- length(regular_samples)
  if (n_regular > 0) {
    sample_colors <- setNames(
      colorRampPalette(brewer.pal(min(9, max(3, n_regular)), "Set1"))(n_regular),
      regular_samples
    )
  } else {
    sample_colors <- c()
  }
  
  # Add colors for ambiguous (lighter shades) and doublet
  ambiguous_colors <- setNames(rep("gray80", length(ambiguous_samples)), ambiguous_samples)
  all_colors <- c(sample_colors, ambiguous_colors, "Doublet" = "gray40")
  
  # Create stacked barplot
  p <- ggplot(run_data, aes(x = Cluster, y = Percentage, fill = Sample)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_fill_manual(values = all_colors) +
    labs(
      title = paste("Souporcell Cluster Composition (Corrected) -", run_name),
      x = "Souporcell Cluster",
      y = "Proportion (%)",
      fill = "Sample"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # Add assignment labels on top (only for numeric clusters, not doublet)
  assignments <- cluster_assignments[[run_name]]
  if (length(assignments) > 0) {
    label_data <- data.frame(
      Cluster = factor(names(assignments), levels = cluster_levels),
      Label = sapply(assignments, function(x) x$sample),
      y = 102
    )
    
    p <- p + geom_text(data = label_data, 
                       aes(x = Cluster, y = y, label = Label, fill = NULL),
                       inherit.aes = FALSE, fontface = "bold", size = 3.5)
  }
  
  # Save plot
  ggsave(file.path(outdir, paste0(run_name, "_cluster_composition_corrected.png")),
         plot = p, width = 10, height = 6, dpi = 300)
  ggsave(file.path(outdir, paste0(run_name, "_cluster_composition_corrected.pdf")),
         plot = p, width = 10, height = 6)
  
  cat("Saved corrected plot for", run_name, "\n")
}

# Save summary table (original)
summary_df <- data.frame()
for (run_name in names(cluster_assignments)) {
  for (cluster in names(cluster_assignments[[run_name]])) {
    assignment <- cluster_assignments[[run_name]][[cluster]]
    summary_df <- rbind(summary_df, data.frame(
      Run = run_name,
      Cluster = cluster,
      Assigned_Sample = assignment$sample,
      Confidence_Percent = round(assignment$percentage, 2),
      Total_Cells = assignment$count,
      stringsAsFactors = FALSE
    ))
  }
}

write.table(summary_df, 
            file = file.path(outdir, "souporcell_assignments_summary.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nSummary table saved to:", file.path(outdir, "souporcell_assignments_summary.tsv"), "\n")

# Save corrected summary table
corrected_summary_df <- data.frame()
for (run_name in names(cluster_assignments)) {
  for (cluster in names(cluster_assignments[[run_name]])) {
    assignment <- cluster_assignments[[run_name]][[cluster]]
    assigned_sample <- assignment$sample
    
    # Get corrected data for this cluster
    cluster_data <- all_results_corrected[[run_name]]
    cluster_data <- cluster_data[cluster_data$Cluster == cluster, ]
    
    # Calculate counts
    # Assigned cells = cells from assigned patient's samples (ASB + GB)
    assigned_count <- sum(cluster_data$Count[grepl(paste0("^", assigned_sample, "_"), cluster_data$Sample)])
    
    # Ambiguous cells = cells assigned to {patient}_Ambiguous
    ambiguous_sample_name <- paste0(assigned_sample, "_Ambiguous")
    ambiguous_count <- sum(cluster_data$Count[cluster_data$Sample == ambiguous_sample_name])
    
    # Total cells in cluster (no doublets here - they moved to Doublet cluster)
    total_cells <- assigned_count + ambiguous_count
    
    # Confidence = 100% since all cells in this cluster are accounted for
    corrected_summary_df <- rbind(corrected_summary_df, data.frame(
      Run = run_name,
      Cluster = cluster,
      Assigned_Sample = assigned_sample,
      Assigned_Cells = assigned_count,
      Ambiguous_Cells = ambiguous_count,
      Total_Cells = total_cells,
      Confidence_Percent = 100,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add Doublet cluster summary
  doublet_data <- all_results_corrected[[run_name]]
  doublet_data <- doublet_data[doublet_data$Cluster == "Doublet", ]
  if (nrow(doublet_data) > 0) {
    total_doublet_cells <- sum(doublet_data$Count)
    
    corrected_summary_df <- rbind(corrected_summary_df, data.frame(
      Run = run_name,
      Cluster = "Doublet",
      Assigned_Sample = "Doublet",
      Assigned_Cells = total_doublet_cells,
      Ambiguous_Cells = 0,
      Total_Cells = total_doublet_cells,
      Confidence_Percent = 100,
      stringsAsFactors = FALSE
    ))
  }
}

write.table(corrected_summary_df, 
            file = file.path(outdir, "souporcell_reassignments_summary.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("Corrected summary table saved to:", file.path(outdir, "souporcell_reassignments_summary.tsv"), "\n")

# Save corrected run_data for each run
for (run_name in names(all_results_corrected)) {
  run_data_corrected <- all_results_corrected[[run_name]]
  
  if (nrow(run_data_corrected) > 0) {
    # Save to TSV file
    output_file <- file.path(outdir, paste0(run_name, "_corrected_run_data.tsv"))
    write.table(run_data_corrected, 
                file = output_file,
                row.names = FALSE, quote = FALSE, sep = "\t")
    
    cat("Saved corrected run_data for", run_name, "to:", output_file, "\n")
  }
}

# CREATE CORRECTED BARCODE ASSIGNMENTS (excluding doublets)
cat("\n", rep("=", 80), "\n", sep = "")
cat("CREATING CORRECTED BARCODE ASSIGNMENTS\n")
cat(rep("=", 80), "\n", sep = "")

corrected_sample_barcodes <- list()

for (run_name in names(run_to_snp_dir)) {
  
  cat("\nProcessing", run_name, "...\n")
  
  snp_dir <- run_to_snp_dir[[run_name]]
  
  # Read Souporcell clusters
  clusters_file <- file.path(snp_dir, "clusters.tsv")
  if (!file.exists(clusters_file)) {
    warning("clusters.tsv not found for ", run_name)
    next
  }
  
  SNPs_demultiplex <- read.table(clusters_file, row.names = 'barcode', header = TRUE, stringsAsFactors = FALSE)
  
  # Rename assignments with "/" to "doublet"
  SNPs_demultiplex$assignment <- ifelse(grepl("/", SNPs_demultiplex$assignment), 
                                        "doublet", 
                                        SNPs_demultiplex$assignment)
  
  # Get cluster assignments for this run
  assignments <- cluster_assignments[[run_name]]
  
  # Process each cluster (excluding doublet cluster)
  for (cluster in names(assignments)) {
    assignment <- assignments[[cluster]]
    assigned_patient <- assignment$sample
    
    # Get all barcodes in this cluster
    cluster_barcodes <- rownames(SNPs_demultiplex[SNPs_demultiplex$assignment == cluster, ])
    
    # Get samples for this run
    samples <- run_to_samples[[run_name]]
    
    # Check each barcode and assign to sample
    for (sample_dir in samples) {
      sample_name <- basename(sample_dir)
      patient_id <- get_patient_id(sample_name)
      
      # Only process samples belonging to the assigned patient
      if (patient_id == assigned_patient) {
        sample_barcodes <- Samples_barcodes[[sample_name]]
        
        if (!is.null(sample_barcodes)) {
          # Find barcodes in both cluster and sample
          common_barcodes <- intersect(cluster_barcodes, sample_barcodes)
          
          if (length(common_barcodes) > 0) {
            # Add to corrected list
            if (is.null(corrected_sample_barcodes[[sample_name]])) {
              corrected_sample_barcodes[[sample_name]] <- character(0)
            }
            corrected_sample_barcodes[[sample_name]] <- c(corrected_sample_barcodes[[sample_name]], common_barcodes)
            
            cat("  Added", length(common_barcodes), "barcodes to", sample_name, "\n")
          }
        }
      }
    }
    
    # Handle unassigned/ambiguous barcodes
    # Find barcodes in cluster that don't match any sample
    assigned_barcodes <- character(0)
    for (sample_dir in samples) {
      sample_name <- basename(sample_dir)
      sample_barcodes <- Samples_barcodes[[sample_name]]
      if (!is.null(sample_barcodes)) {
        assigned_barcodes <- c(assigned_barcodes, intersect(cluster_barcodes, sample_barcodes))
      }
    }
    
    unassigned_barcodes <- setdiff(cluster_barcodes, assigned_barcodes)
    
    if (length(unassigned_barcodes) > 0) {
      # Create ambiguous sample name
      ambiguous_sample_name <- paste0(assigned_patient, "_Ambiguous")
      
      if (is.null(corrected_sample_barcodes[[ambiguous_sample_name]])) {
        corrected_sample_barcodes[[ambiguous_sample_name]] <- character(0)
      }
      corrected_sample_barcodes[[ambiguous_sample_name]] <- c(corrected_sample_barcodes[[ambiguous_sample_name]], 
                                                              unassigned_barcodes)
      
      cat("  Added", length(unassigned_barcodes), "barcodes to", ambiguous_sample_name, "\n")
    }
  }
}

# Print summary
cat("\n", rep("-", 80), "\n", sep = "")
cat("CORRECTED BARCODE ASSIGNMENT SUMMARY:\n")
cat(rep("-", 80), "\n", sep = "")
for (sample_name in sort(names(corrected_sample_barcodes))) {
  n_barcodes <- length(corrected_sample_barcodes[[sample_name]])
  cat(sprintf("  %s: %d barcodes\n", sample_name, n_barcodes))
}

# Save corrected barcode assignments to RDS file
output_rds <- file.path(outdir, "corrected_sample_barcodes.rds")
saveRDS(corrected_sample_barcodes, file = output_rds)
cat("\nSaved corrected barcode assignments to:", output_rds, "\n")

# Also save as individual files per sample
barcode_dir <- file.path(outdir, "corrected_barcodes")
dir.create(barcode_dir, recursive = TRUE, showWarnings = FALSE)

for (sample_name in names(corrected_sample_barcodes)) {
  sample_barcode_file <- file.path(barcode_dir, paste0(sample_name, "_barcodes.txt"))
  writeLines(corrected_sample_barcodes[[sample_name]], sample_barcode_file)
  cat("Saved", sample_name, "barcodes to:", sample_barcode_file, "\n")
}

cat("\nAnalysis complete!\n")