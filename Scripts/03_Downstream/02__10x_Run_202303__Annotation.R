
### Defining work dir
setwd('/project/meadlab/kgurashi') # <- change

options(Seurat.object.assay.version = 'v3')
options(future.globals.maxSize= 891289600*10)

# Set the bitmap type to "cairo" to avoid X11 dependency
options(bitmapType = "cairo")

### Loading libraries

library(Seurat)
library(dplyr)
library(SingleR)
library(BiocParallel)
library(stringr)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(ggplot2)
library(harmony)

## - Options

BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-2)

### Defining seed

set.seed(seed = 123, kind = "Mersenne-Twister", normal.kind = "Inversion")

### Loading query data

data_folder = './CG__250425_10x/Data/10x'

dir.create('./CG__250425_10x/Results/10x/Annotation', showWarnings = F)

results_folder = './CG__250425_10x/Results/10x/Annotation'

Query.SeuratObj <- readRDS(file = paste(data_folder, 'SC_Objects', 'Integrated.rds', sep = '/'))

### Loading reference data

REF_DB_Triana.SeuratObj <- readRDS(paste(data_folder, "/Resources/Reference_10x_Datasets/TrianaS_VeltenL_Healthy_Reference_ALL_dataset.rds", sep = ''))
REF_DB_Stewart.SeuratObj <- readRDS(paste(data_folder, "/Resources/Reference_10x_Datasets/Stewart2021_BCells.rds", sep = ''))
REF_DB_Stewart.SeuratObj <- UpdateSeuratObject(REF_DB_Stewart.SeuratObj)

### Processing reference data

### Triana S.

REF_DB_Triana.SeuratObj_label <- paste('REF:TrianaS  ', str_replace(REF_DB_Triana.SeuratObj$ct, "\\d+_", ""),sep = '')
names(REF_DB_Triana.SeuratObj_label) <- names(REF_DB_Triana.SeuratObj$ct)
REF_DB_Triana.SeuratObj <- SCTransform(REF_DB_Triana.SeuratObj, 
                                       assay = 'RNA', 
                                       residual.features = REF_DB_Triana.SeuratObj@assays$RNA@var.features)
DefaultAssay(REF_DB_Triana.SeuratObj) <- 'SCT'

### Stewart A.

REF_DB_Stewart.SeuratObj_label <- paste('REF:StewartA  ', str_replace(REF_DB_Stewart.SeuratObj$Population, "\\d+_", ""),sep = '')
names(REF_DB_Stewart.SeuratObj_label) <- names(REF_DB_Stewart.SeuratObj$Population)
REF_DB_Stewart.SeuratObj <- SCTransform(REF_DB_Stewart.SeuratObj, 
                                        assay = 'RNA', 
                                        residual.features = REF_DB_Stewart.SeuratObj@assays$integrated@var.features)
DefaultAssay(REF_DB_Stewart.SeuratObj) <- 'SCT'

### Running SingleR annotation

singler_results <- SingleR(
  test = as.matrix(Query.SeuratObj@assays$SCT@data), 
  ref = list(TrianaS=REF_DB_Triana.SeuratObj@assays$SCT@data,
             StewartA=REF_DB_Stewart.SeuratObj@assays$SCT@data),
  labels = list(REF_DB_Triana.SeuratObj_label,
                REF_DB_Stewart.SeuratObj_label),
  #BPPARAM = BPPARAM, #Activate BPPARAM for parallelisation
  de.method = 'wilcox'
)

### Saving SingleR annotation output

saveRDS(singler_results, 
        file = paste(results_folder, '/SingleR_annotation_results.RDS', sep=''))

png(filename = paste(results_folder, "/Annotation_Combined_Score_Heatmap.png", sep = ''),
    width = 14,height = 10,units = 'in',res = 300, type='cairo')
plotScoreHeatmap(singler_results,scores.use = 0, show.labels = F, max.labels = 70)
dev.off()

png(filename = paste(results_folder, "/Annotation_REF_TrianaS_Score_Heatmap.png", sep = ''),
    width = 12,height = 8,units = 'in',res = 300, type='cairo')
plotScoreHeatmap(singler_results,scores.use = 1, show.labels = F)
dev.off()

png(filename = paste(results_folder, "/Annotation_REF_StewartA_Score_Heatmap.png", sep = ''),
    width = 12,height = 8,units = 'in',res = 300, type='cairo')
plotScoreHeatmap(singler_results,scores.use = 2, show.labels = F)
dev.off()

### Adding annotation metadata to query object

Query.SeuratObj$Combined_celltype_annotation <- singler_results$labels
Query.SeuratObj$Triana_celltype_annotation <- singler_results$orig.results$TrianaS$labels
Query.SeuratObj$Stewart_celltype_annotation <- singler_results$orig.results$StewartA$labels

table(Query.SeuratObj$Combined_celltype_annotation)

Query.SeuratObj$Celltype_draft <- singler_results$labels

### Refining combined annotation by merging highly correlating annotations

Expression_mat <- AverageExpression(Query.SeuratObj, assays = 'SCT', group.by = 'Combined_celltype_annotation')$SCT
Corr <- cor(as.matrix(Expression_mat), method = 'spearman')
ComplexHeatmap::Heatmap(Corr)

### Removing cells poorly represented or outliers

labels_with_enough_cells = names(table(Query.SeuratObj$Combined_celltype_annotation))[grep('FALSE', table(Query.SeuratObj$Combined_celltype_annotation) < 20)]

Query.SeuratObj = subset(Query.SeuratObj,
                         cells = WhichCells(Query.SeuratObj, expression = Combined_celltype_annotation %in% labels_with_enough_cells))

Query.SeuratObj$Combined_celltype_annotation <- sub("REF:[^ ]*  ", "", Query.SeuratObj$Combined_celltype_annotation)

### Plotting final results

png(filename = paste(results_folder, "/Annotated_sample_umap.png", sep = ''),
    width = 8.5,height = 6,units = 'in',res = 300, type='cairo')

DimPlot(Query.SeuratObj, 
        reduction = 'h_umap', 
        group.by = 'Combined_celltype_annotation', 
        pt.size = 1) + scale_color_manual(values = c('Naive' = "#B15928",
                                                     "Classical_Memory" = "#6A3D9A",
                                                     "DN" = "#FDBF6F",
                                                     "IgM_Memory" = "#33A02C",
                                                     "CD56dimCD16+ NK cells" = 'pink',
                                                     "Mature naive B cells" = 'blue',
                                                     "Nonswitched memory B cells" = "darkgreen",
                                                     "Class switched memory B cells" = "#E31A1C",
                                                     "CD4+ naive T cells" = "yellow",
                                                     "CD4+ memory T cells" = "yellow",
                                                     "CD8+ naive T cells" = "cyan",
                                                     "CD8+ effector memory T cells" = "cyan",
                                                     "Non-classical monocytes" = "orange",
                                                     "Classical Monocytes" = "orange",
                                                     "Transitional" = "#A6CEE3")) + ylab('UMAP-2') + xlab('UMAP-1') + ggtitle('Cell types')

dev.off()

######################
### Checks and tidy up
######################

# Finding and removing contamination cells

# Define the file path
#file_path <-  paste(results_folder, 'CB_Noise_Contamination.csv', sep = '/')

# Check if the file exists
#if (file.exists(file_path)) {
  # Import unwanted_cell_names from the file
#  unwanted_cell_names <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)[, 1]
#} else {
#  # Extract the UMAP embeddings and metadata from the Seurat object
#  umap_data <- as.data.frame(Query.SeuratObj@reductions$h_umap@cell.embeddings)
#  umap_data$celltype <- Query.SeuratObj$Combined_celltype_annotation
  
#  # Rename columns for consistency (if needed)
#  colnames(umap_data) <- c("UMAP_1", "UMAP_2", "celltype")
  
#  # Filter cells for group1
#  group1 <- subset(umap_data, UMAP_2 >= -3 & UMAP_2 <= 3 & UMAP_1 >= -11 & UMAP_1 <= -7)
  
#  # Filter cells for group2
#  group2 <- subset(umap_data, UMAP_2 >= -3 & UMAP_2 <= 5 & UMAP_1 >= 11 & UMAP_1 <= 15)
  
#  # Combine both groups into a single data frame
#  unwanted_cells <- rbind(group1, group2)
  
#  # Get the cell names of unwanted cells
#  unwanted_cell_names <- rownames(unwanted_cells)
  
#  # Write the unwanted cell names to the file
#  write.table(x = unwanted_cell_names, 
#              file = file_path, 
#              quote = FALSE, row.names = FALSE, sep = ',', col.names = FALSE)
#}

## Subset the Seurat object to exclude unwanted cells
#Query.SeuratObj_filtered <- subset(Query.SeuratObj, cells = setdiff(Cells(Query.SeuratObj), unwanted_cell_names))

## Verify the subsetted object
#print(Query.SeuratObj_filtered)

## Plot the filtered Seurat object
#DimPlot(Query.SeuratObj_filtered, 
#        reduction = 'h_umap', 
#        group.by = 'Combined_celltype_annotation', 
#        pt.size = 1) + scale_color_manual(values = c('Naive' = "#B15928",
#                                                     "Classical_Memory" = "#6A3D9A",
#                                                     "DN" = "#FDBF6F",
#                                                     "IgM_Memory" = "#33A02C",
#                                                     "CD56dimCD16+ NK cells" = 'pink',
#                                                     "Mature naive B cells" = 'blue',
#                                                     "Nonswitched memory B cells" = "darkgreen",
#                                                     "Class switched memory B cells" = "#E31A1C",
#                                                     "Transitional" = "#A6CEE3")) + 
#  ylab('UMAP-2') + xlab('UMAP-1') + ggtitle('Cell types')

## RNA normalisation

#DefaultAssay(Query.SeuratObj_filtered) <- 'RNA'
#Query.SeuratObj_filtered@assays$SCT <- NULL

#Query.SeuratObj_filtered <- SCTransform(Query.SeuratObj_filtered, verbose = FALSE)
#Query.SeuratObj_filtered <- FindVariableFeatures(Query.SeuratObj_filtered, assay = 'SCT')
#Query.SeuratObj_filtered <- ScaleData(Query.SeuratObj_filtered, assay = 'SCT')

## UMAP non-integrated

#Query.SeuratObj_filtered <- RunUMAP(Query.SeuratObj_filtered, dims=1:12, reduction = 'pca', 
#                                    reduction.name = 'umap', reduction.key = 'umap', 
#                                    spread=1.5, min.dist=0.3, n.neighbors=50, 
#                                    metric="cosine", n.epochs = 200, 
#                                    negative.sample.rate = 5, return.model = T, seed.use = 123)
#
## Harmony integration

#Query.SeuratObj_filtered <- RunPCA(Query.SeuratObj_filtered, 
#                                   verbose = FALSE, assay = 'SCT')

## These are now standard steps in the Seurat workflow for visualization and clustering
#elbow_plot <- ElbowPlot(object = Query.SeuratObj_filtered, ndims = 40)

#print(elbow_plot)

#Query.SeuratObj_filtered <- RunHarmony(Query.SeuratObj_filtered, 
#                                       group.by.vars = c('Identity'), 
#                                       dims.use = 1:12, max_iter = 50)

## UMAP integrated

#Query.SeuratObj_filtered <- RunUMAP(Query.SeuratObj_filtered, dims=1:12, 
#                                    reduction = 'harmony', reduction.name = 'h_umap', 
#                                    reduction.key = 'h_umap', spread=1.5, 
#                                   min.dist=0.3, n.neighbors=50, metric="cosine", 
#                                    n.epochs = 200, negative.sample.rate = 5, 
#                                    return.model = T, seed.use = 123)
#
## Check umap results

#cells <- colnames(Query.SeuratObj_filtered)
#cell_types <- Query.SeuratObj_filtered@meta.data[cells, 'Combined_celltype_annotation']
#cells_highlight <- split(cells, cell_types)

#DimPlot(Query.SeuratObj_filtered, reduction = 'h_umap', cells.highlight = cells_highlight,
#        group.by='Combined_celltype_annotation', 
#        label = F, pt.size = 1.4, sizes.highlight = 0.9, alpha=0.8, order = T,
#        repel = T) + theme_void() + scale_color_manual(values = c('Naive' = "#c97c5d",
#                                                                  "Classical_Memory" = "#bf0603",
#                                                                  "DN" = "grey40",
#                                                                  "IgM_Memory" = "pink",
#                                                                  "CD56dimCD16+ NK cells" = 'orange',
#                                                                  "Mature naive B cells" = '#593d3b',
#                                                                  "Nonswitched memory B cells" = "#96c7c6",
#                                                                  "Class switched memory B cells" = "#004e89",
#                                                                  "Transitional" = "darkgreen"),
#                                                       na.value = '#DCDADB') + 
#  theme(plot.title = element_text(hjust=0.5, vjust = 0.5),
#        plot.subtitle = element_text(hjust=0.5, vjust = 0.5),
#        legend.position = 'bottom') +
#  ggtitle(label = 'ALL', subtitle = paste('N = ', length(cells), sep = ''))

## Result with a chosen patient

#sample <- 'VA013_ASB'

#cells <- WhichCells(Query.SeuratObj_filtered, expression = Identity %in% sample)
#cell_types <- Query.SeuratObj_filtered@meta.data[cells, 'Combined_celltype_annotation']
#cells_highlight <- split(cells, cell_types)

#DimPlot(Query.SeuratObj_filtered, reduction = 'h_umap', cells.highlight = cells_highlight,
#        group.by='Combined_celltype_annotation', 
#        label = F, pt.size = 1.4, sizes.highlight = 0.9, alpha=0.8, order = T,
#        repel = T) + theme_void() + scale_color_manual(values = c('Naive' = "#c97c5d",
#                                                                  "Classical_Memory" = "#bf0603",
#                                                                  "DN" = "grey40",
#                                                                  "IgM_Memory" = "pink",
#                                                                  "CD56dimCD16+ NK cells" = 'orange',
#                                                                  "Mature naive B cells" = '#593d3b',
#                                                                  "Nonswitched memory B cells" = "#96c7c6",
#                                                                  "Class switched memory B cells" = "#004e89",
#                                                                  "Transitional" = "darkgreen"),
#                                                       na.value = '#DCDADB') + 
#  theme(plot.title = element_text(hjust=0.5, vjust = 0.5),
#        plot.subtitle = element_text(hjust=0.5, vjust = 0.5),
#        legend.position = 'none') +
#  ggtitle(label = sprintf('%s', sample), subtitle = paste('N = ', length(cells), sep = ''))

#### Saving final results

saveRDS(object = Query.SeuratObj, 
        file = paste(data_folder, 'SC_Objects', 'Annotated.rds', sep = '/'))

## - Session Info

writeLines(capture.output(sessionInfo()), file.path(dirname(rstudioapi::getSourceEditorContext()$path), '02_Session_Info.txt'))

