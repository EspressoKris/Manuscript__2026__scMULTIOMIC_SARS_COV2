
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
library(ggthemes)
library(SCpubr)

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
REF_DB_Terekhova.SeuratObj <- readRDS(paste(data_folder, "/Resources/Reference_10x_Datasets/Terekhova_PBMCs2023/b_cells_rna.rds", sep = ''))
REF_DB_Terekhova.UMAP <- read.csv(paste(data_folder, "/Resources/Reference_10x_Datasets/Terekhova_PBMCs2023/b_cells_umap.csv", sep = ''), row.names = 'X')
REF_DB_Terekhova.UMAP <- CreateDimReducObject(embeddings = as.matrix(REF_DB_Terekhova.UMAP), assay = 'RNA', global = T, key = 'UMAP')
REF_DB_Terekhova.Labels <- read.csv(paste(data_folder, "/Resources/Reference_10x_Datasets/Terekhova_PBMCs2023/b_cells_metadata.csv", sep = ''), row.names = 'X')
REF_DB_Terekhova.SeuratObj <- Seurat::CreateSeuratObject(counts = REF_DB_Terekhova.SeuratObj, assay = 'RNA', meta.data = REF_DB_Terekhova.Labels)
REF_DB_Terekhova.SeuratObj@reductions$UMAP <- REF_DB_Terekhova.UMAP

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

### Terekhova M.

REF_DB_Terekhova.SeuratObj_label <- paste('REF:TerekhovaM  ', str_replace(REF_DB_Terekhova.SeuratObj$Cluster_names, "\\d+_", ""),sep = '')
names(REF_DB_Terekhova.SeuratObj_label) <- names(REF_DB_Terekhova.SeuratObj$Cluster_names)
REF_DB_Terekhova.SeuratObj <- SCTransform(REF_DB_Terekhova.SeuratObj, 
                                        assay = 'RNA')
DefaultAssay(REF_DB_Terekhova.SeuratObj) <- 'SCT'

### Running SingleR annotation

singler_results <- SingleR(
  test = as.matrix(Query.SeuratObj@assays$SCT@data), 
  ref = list(TrianaS=REF_DB_Triana.SeuratObj@assays$SCT@data,
             StewartA=REF_DB_Stewart.SeuratObj@assays$SCT@data,
             TerekhovaM=REF_DB_Terekhova.SeuratObj@assays$SCT@data),
  labels = list(REF_DB_Triana.SeuratObj_label,
                REF_DB_Stewart.SeuratObj_label,
                REF_DB_Terekhova.SeuratObj_label),
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

png(filename = paste(results_folder, "/Annotation_REF_TerekhovaM_Score_Heatmap.png", sep = ''),
    width = 12,height = 8,units = 'in',res = 300, type='cairo')
plotScoreHeatmap(singler_results,scores.use = 3, show.labels = F)
dev.off()

### Adding annotation metadata to query object

Query.SeuratObj$Combined_celltype_annotation <- singler_results$labels
Query.SeuratObj$Triana_celltype_annotation <- singler_results$orig.results$TrianaS$labels
Query.SeuratObj$Stewart_celltype_annotation <- singler_results$orig.results$StewartA$labels
Query.SeuratObj$Terekhova_celltype_annotation <- singler_results$orig.results$TerekhovaM$labels

table(Query.SeuratObj$Combined_celltype_annotation)

Query.SeuratObj$Celltype_draft <- singler_results$orig.results$TrianaS$labels
Query.SeuratObj@meta.data[grep(TRUE,
                               str_detect(pattern = "REF:TerekhovaM  Non-switched memory", 
                                          string = singler_results$orig.results$TerekhovaM$labels)),
                          'Celltype_draft'] <- "Non-switched memory"
Query.SeuratObj@meta.data[grep(TRUE,
                               str_detect(pattern = "REF:TerekhovaM  Switched memory", 
                                          string = singler_results$orig.results$TerekhovaM$labels)),
                          'Celltype_draft'] <- "Switched memory"
Query.SeuratObj@meta.data[grep(TRUE,
                               str_detect(pattern = "REF:TerekhovaM  Atypical memory", 
                                          string = singler_results$orig.results$TerekhovaM$labels)),
                          'Celltype_draft'] <- "Atypical memory"
Query.SeuratObj@meta.data[grep(TRUE,
                               str_detect(pattern = "REF:StewartA  Naive", 
                                          string = singler_results$orig.results$TerekhovaM$labels)),
                          'Celltype_draft'] <- "Naive"
Query.SeuratObj@meta.data[grep(TRUE,
                               str_detect(pattern = "REF:StewartA  Transitional", 
                                          string = singler_results$orig.results$StewartA$labels)),
                          'Celltype_draft'] <- "Transitional"

### Refining combined annotation by merging highly correlating annotations

Query.SeuratObj$Celltype_draft <- sub("REF:[^ ]*  ", "", Query.SeuratObj$Celltype_draft)

digest::digest(Query.SeuratObj$Celltype_draft)
#4e4036d5b598fe02c54a799b5629918c

###

set.seed(123)

Query.SeuratObj <- FindNeighbors(Query.SeuratObj, assay = 'harmony')

set.seed(123)

Query.SeuratObj <- FindClusters(Query.SeuratObj, algorithm=4, resolution = 1, random.seed = 123)

DimPlot(Query.SeuratObj, group.by = 'seurat_clusters', label = T, reduction = 'h_umap')

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = seurat_clusters %in% c('1', '3', '4', '8', '12', '10', '15')),
                          'Celltype_draft'] <- 'B naive'

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = seurat_clusters %in% c('9')),
                          'Celltype_draft'] <- 'B transitional'

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = seurat_clusters %in% c('5', '6')),
                          'Celltype_draft'] <- 'B switched memory'

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = seurat_clusters %in% c('11')),
                          'Celltype_draft'] <- 'B atypical memory'

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = seurat_clusters %in% c('13', '2', '7')),
                          'Celltype_draft'] <- 'B non switched memory'

Idents(Query.SeuratObj) <- Query.SeuratObj$seurat_clusters

set.seed(123)

DEG <- FindAllMarkers(Query.SeuratObj, assay = 'SCT')
DEG$pct2_vs_pct1 <- DEG$pct.2 - DEG$pct.1
DEG_top10 <- DEG %>%
  group_by(cluster) %>%
  arrange(desc(pct2_vs_pct1), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

DimPlot(Query.SeuratObj, group.by = 'Celltype_draft', label = T, reduction = 'h_umap')

###

T_NK_subset.SeuratObj <- subset(Query.SeuratObj, 
                             cells = WhichCells(Query.SeuratObj, 
                                                expression = seurat_clusters %in% c('14')))

set.seed(123)

T_NK_subset.SeuratObj <- FindNeighbors(T_NK_subset.SeuratObj, assay = 'SCT', 'harmony')
T_NK_subset.SeuratObj <- FindClusters(T_NK_subset.SeuratObj, algorithm = 4, 
                                      resolution = 0.5, graph.name = 'SCT_snn')

T_NK_subset.SeuratObj$temp_cluster <- 'Other'
T_NK_subset.SeuratObj@meta.data[names(T_NK_subset.SeuratObj$seurat_clusters),
                                'temp_cluster'] <- T_NK_subset.SeuratObj$seurat_clusters

DimPlot(T_NK_subset.SeuratObj, group.by = 'temp_cluster', label = T)
DimPlot(T_NK_subset.SeuratObj, group.by = 'Celltype_draft', label = T)
FeaturePlot(T_NK_subset.SeuratObj, features = 'CD8A')

Query.SeuratObj@meta.data[WhichCells(T_NK_subset.SeuratObj, 
                                      expression = temp_cluster %in% c('4')),
                           'Celltype_draft'] <- 'CD8+ T'

Query.SeuratObj@meta.data[WhichCells(T_NK_subset.SeuratObj, 
                                     expression = temp_cluster %in% c('3', '1')),
                          'Celltype_draft'] <- 'CD4+ T'

Query.SeuratObj@meta.data[WhichCells(T_NK_subset.SeuratObj, 
                                     expression = temp_cluster %in% c('2')),
                          'Celltype_draft'] <- 'NK'

DimPlot(Query.SeuratObj, group.by = 'Celltype_draft', label = T, reduction = 'h_umap')

###

Myeloid_subset.SeuratObj <- subset(Query.SeuratObj, 
                                cells = WhichCells(Query.SeuratObj, 
                                                   expression = seurat_clusters %in% c('16')))

set.seed(123)

Myeloid_subset.SeuratObj <- FindNeighbors(Myeloid_subset.SeuratObj, assay = 'SCT', 'harmony')
Myeloid_subset.SeuratObj <- FindClusters(Myeloid_subset.SeuratObj, algorithm = 4, 
                                      resolution = 1, graph.name = 'SCT_snn')

Myeloid_subset.SeuratObj$temp_cluster <- 'Other'
Myeloid_subset.SeuratObj@meta.data[names(Myeloid_subset.SeuratObj$seurat_clusters),
                                'temp_cluster'] <- Myeloid_subset.SeuratObj$seurat_clusters

DimPlot(Myeloid_subset.SeuratObj, group.by = 'temp_cluster', label = T)
DimPlot(Myeloid_subset.SeuratObj, group.by = 'Celltype_draft', label = T)
FeaturePlot(Myeloid_subset.SeuratObj, features = 'CD14')

Query.SeuratObj@meta.data[WhichCells(Myeloid_subset.SeuratObj, 
                                     expression = temp_cluster %in% c('2', '1')),
                          'Celltype_draft'] <- 'CD14+ monocyte'

Query.SeuratObj@meta.data[WhichCells(Myeloid_subset.SeuratObj, 
                                     expression = temp_cluster %in% c('3')),
                          'Celltype_draft'] <- 'CD16+ monocyte'

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Triana_celltype_annotation %in% c('REF:TrianaS  Conventional dendritic cell 2',
                                                                                    'REF:TrianaS  Conventional dendritic cell 1',
                                                                                    'REF:TrianaS  Plasmacytoid dendritic cells')),
                          'Celltype_draft'] <- 'Dendritic cells'

###

B_subcluster.SeuratObj <- subset(Query.SeuratObj, 
                                   cells = WhichCells(Query.SeuratObj, 
                                                      expression = seurat_clusters %in% c('15')))

set.seed(123)

B_subcluster.SeuratObj <- FindNeighbors(B_subcluster.SeuratObj, assay = 'SCT', 'harmony')
B_subcluster.SeuratObj <- FindClusters(B_subcluster.SeuratObj, algorithm = 4, 
                                         resolution = 1, graph.name = 'SCT_snn')

B_subcluster.SeuratObj$temp_cluster <- 'Other'
B_subcluster.SeuratObj@meta.data[names(B_subcluster.SeuratObj$seurat_clusters),
                                   'temp_cluster'] <- B_subcluster.SeuratObj$seurat_clusters

DimPlot(B_subcluster.SeuratObj, group.by = 'temp_cluster', label = T)
DimPlot(B_subcluster.SeuratObj, group.by = 'Celltype_draft', label = T)
FeaturePlot(B_subcluster.SeuratObj, features = 'TCL1A')
VlnPlot(B_subcluster.SeuratObj, features = 'IGHD', group.by = 'temp_cluster')

Query.SeuratObj@meta.data[WhichCells(B_subcluster.SeuratObj, 
                                     expression = temp_cluster %in% c('1', '2', '3', '4')),
                          'Celltype_draft'] <- 'B naive'

Query.SeuratObj@meta.data[WhichCells(B_subcluster.SeuratObj, 
                                     expression = temp_cluster %in% c('5')),
                          'Celltype_draft'] <- 'B non switched memory'

### Plotting final results

Query.SeuratObj$Celltype <- Query.SeuratObj$Celltype_draft

Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('B naive')),
                                     'Celltype'] <- "Naive B-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('B non switched memory')),
                          'Celltype'] <- "IgM memory B-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('B switched memory')),
                          'Celltype'] <- "Class switched memory B-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('B transitional')),
                          'Celltype'] <- "Transitional B-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('B atypical memory')),
                          'Celltype'] <- "Atypical memory B-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('Dendritic cells')),
                          'Celltype'] <- "Dendritic cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('CD4+ T')),
                          'Celltype'] <- "CD4+ T-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('CD14+ monocyte')),
                          'Celltype'] <- "CD14+ monocytes"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('CD8+ T')),
                          'Celltype'] <- "CD8+ T-cells"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('CD16+ monocyte')),
                          'Celltype'] <- "CD16+ monocytes"
Query.SeuratObj@meta.data[WhichCells(Query.SeuratObj, 
                                     expression = Celltype %in% c('NK')),
                          'Celltype'] <- "NK cells"

### 

png(filename = paste(results_folder, "/Annotated_sample_umap.png", sep = ''),
    width = 7.5,height = 8.5,units = 'in',res = 300, type='cairo')

cols <- c(
  "Naive B-cells"                = "#E15759",
  "Transitional B-cells"         = "#EDC948",
  "Atypical memory B-cells"      = "#F28E2B",
  "IgM memory B-cells"           = "#76B7B2",
  "Class switched memory B-cells"= "#4E79A7",
  "CD4+ T-cells"                 = "#9C755F",
  "CD8+ T-cells"                 = "#BAB0AC",
  "NK cells"                     = "#59A14F",
  "CD14+ monocytes"              = "#FF9DA7",
  "CD16+ monocytes"              = "#B07AA1",
  "Dendritic cells"              = "#B6992D"
)

SCpubr::do_DimPlot(
  Query.SeuratObj,
  reduction = "h_umap",
  group.by = "Celltype",
  border.size = 2,
  pt.size = 1
) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(
      ncol = 3,
      override.aes = list(
        shape  = 21,              # fillable legend key
        fill   = unname(cols),    # make fill match your annotation colours
        colour = "black",         # black stroke
        stroke = 0.8,
        size   = 3
      )
    )
  )

dev.off()

Query.SeuratObj$SCT_snn_res.1 <- NULL
Query.SeuratObj$Celltype_draft <- NULL

###

Query.B.SeuratObj <- subset(Query.SeuratObj, cells = WhichCells(Query.SeuratObj,
                                                                expression = Celltype %in% c('Naive B-cells',
                                                                                             'Transitional B-cells',
                                                                                             'Atypical memory B-cells',
                                                                                             'IgM memory B-cells',
                                                                                             'Class switched memory B-cells')))

Query.B.SeuratObj <- PrepSCTFindMarkers(Query.B.SeuratObj)

set.seed(123)

# Run PCA on SCT assay
Query.B.SeuratObj <- RunPCA(
  Query.B.SeuratObj, 
  assay = 'SCT', 
  features = VarGenes.No.IG,
  verbose = TRUE
)

digest::digest(Query.SeuratObj@reductions$pca)
#3246ac49a82241982b6c261a9f658532

seed.use = 123

Query.B.SeuratObj <- RunHarmony(Query.B.SeuratObj, 
                               group.by.vars = c('Timepoint', 'Donor', 'Identity'), 
                               lambda=c(0.5,1,1),
                               dims.use = 1:11, max_iter = 50)

seed.use = 123

Query.B.SeuratObj <- RunUMAP(Query.B.SeuratObj, dims=1:11, reduction = 'harmony', 
                            reduction.name = 'h_umap', reduction.key = 'h_umap', 
                            spread=0.35, min.dist=0.75, n.neighbors=50, 
                            n.epochs = 200, return.model = T, seed.use = 123)

DimPlot(Query.B.SeuratObj, group.by = 'Celltype', label = T, reduction = 'h_umap')

#### Saving final results

saveRDS(object = Query.SeuratObj, 
        file = paste(data_folder, 'SC_Objects', 'Annotated.rds', sep = '/'))

## - Session Info

library(tools)

script_path <- rstudioapi::getSourceEditorContext()$path
script_path_no_ext <- file_path_sans_ext(script_path)

writeLines(capture.output(sessionInfo()), paste0(script_path_no_ext, '__Session_Info.txt'))
