# Enter commands in R (or R studio, if installed)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
  }
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages('ggplot2')
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages('dplyr')
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages('tidyr')
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages('reshape2')
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")
}
if (!requireNamespace("bluster", quietly = TRUE)) {
  BiocManager::install("bluster")
}
if (!requireNamespace("multtest", quietly = TRUE)) {
  BiocManager::install("multtest")
}
if (!requireNamespace("metap", quietly = TRUE)) {
  install.packages('metap')
}
if (!requireNamespace("mutoss", quietly = TRUE)) {
  install.packages('mutoss')
}
if (!requireNamespace("qqconf", quietly = TRUE)) {
  install.packages('qqconf')
  }
if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
  BiocManager::install("ExperimentHub")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}
if (!requireNamespace("harmony", quietly = TRUE)) {
  install.packages("harmony")
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}
if (!"MuSiC" %in% rownames(installed.packages())) {
  devtools::install_github('xuranw/MuSiC')
}
if (!"xbioc" %in% rownames(installed.packages())) {
  devtools::install_github('crhisto/xbioc')
}
if (!"MuSiC2" %in% rownames(installed.packages())) {
  devtools::install_github('Jiaxin-Fan/MuSiC2')
}
if (!requireNamespace("TOAST", quietly = TRUE)) {
  BiocManager::install("TOAST")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
if (!requireNamespace("celldex", quietly = TRUE)) {
  BiocManager::install("celldex")
}
if (!requireNamespace("SingleR", quietly = TRUE)) {
  BiocManager::install("SingleR")
}
if (!"presto" %in% rownames(installed.packages())) {
  devtools::install_github('immunogenomics/presto')
}
library(presto)
library(SingleR)
library(celldex)
library(Biobase)
library(MuSiC2)
library(MuSiC)
library(TOAST)
library(biomaRt)
library(SingleCellExperiment)
library(harmony)
library(ExperimentHub)
library(DESeq2)
library(mutoss)
library(qqconf)
library(multtest)
library(metap)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scDblFinder)
library(dplyr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)




TD0 = ReadMtx( mtx = "GSE310060_RAW/T/matrix.mtx",
                features = "GSE310060_RAW/T/features.tsv",
                cells = "GSE310060_RAW/T/barcodes.tsv")

SauretTD0 = CreateSeuratObject(counts = TD0,
                             project = "Day0",
                             min.cells = 3,
                             min.features = 200)

TD7 = ReadMtx( mtx = "GSE310060_RAW/TiMKD7/matrix.mtx",
                features = "GSE310060_RAW/TiMKD7/features.tsv",
                cells = "GSE310060_RAW/TiMKD7/barcodes.tsv")

SauretTD7 = CreateSeuratObject(counts = TD7,
                               project = "Day7",
                               min.cells = 3,
                               min.features = 200)

TD14 = ReadMtx( mtx = "GSE310060_RAW/TiMKD14/matrix.mtx",
               features = "GSE310060_RAW/TiMKD14/features.tsv",
               cells = "GSE310060_RAW/TiMKD14/barcodes.tsv")

SauretTD14 = CreateSeuratObject(counts = TD14,
                                project = "Day14",
                                min.cells = 3,
                                min.features = 200)

merged = merge(SauretTD0, y = c(SauretTD7, SauretTD14), add.cell.ids = ls()[1:3], project = "Chemical_Reporogramming")


merged$sample = rownames(merged@meta.data)

merged@meta.data = separate(merged@meta.data, col = 'sample', into = c("sample", "barcode"), sep = '_')

merged$percent.mt = PercentageFeatureSet(merged, pattern = '^MT-')

VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave("Violin_Raw.png", width = 12, height = 6, dpi = 300)

FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") +
  geom_smooth(method = 'lm') + ggtitle("Feature Scatter Plot - RAW")
ggsave("Scatter_Raw.png", width = 7, height = 6, dpi = 300)

merged_filtered = subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                      percent.mt < 20)

merged_filtered = NormalizeData(object = merged_filtered)
merged_filtered = FindVariableFeatures(object = merged_filtered)
merged_filtered = ScaleData(object = merged_filtered)
merged_filtered = RunPCA(object = merged_filtered)

VlnPlot(merged_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave("Violin_Norm.png", width = 12, height = 6, dpi = 300)

FeatureScatter(merged_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") +
  geom_smooth(method = 'lm') + ggtitle("Feature Scatter Plot - Normalized")
ggsave("Scatter_Norm.png", width = 7, height = 6, dpi = 300)

ElbowPlot(merged_filtered) + ggtitle("Elbow Plot")
ggsave("Elbow_Plot.png", width = 7, height = 6, dpi = 300)

merged_filtered = FindNeighbors(object = merged_filtered, dims = 1:15)
merged_filtered = FindClusters(object = merged_filtered)
merged_filtered = RunUMAP(object = merged_filtered, dims = 1:15)

DimPlot(merged_filtered, reduction = 'umap', group.by = 'orig.ident', cols = c('red', 'green', 'blue')) + ggtitle("Dimensional Reduction Plot")
ggsave("DimPlot.png", width = 10, height = 8, dpi = 300)



prop = prop.table(table(merged_filtered$seurat_clusters, merged_filtered$orig.ident), margin = 2)
df = reshape2::melt(prop)
df$Var2 = factor(df$Var2, levels = c("Day0", "Day7", "Day14"))
ggplot(df, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Timepoint", y = "Proportion", fill = "Cluster") +
  ggtitle("Cluster Proportions Across Timepoints")
ggsave("Barplot-Clusters.png", width = 12, height = 9, dpi = 300)


# run Harmony -----------
merged.harmony <- merged_filtered %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = FALSE)

merged.harmony@reductions

merged.harmony.embed <- Embeddings(merged.harmony, "harmony")
merged.harmony.embed[1:10,1:10]
merged.harmony <- merged.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

after <- DimPlot(merged.harmony, reduction = 'umap', group.by = 'orig.ident')

saveRDS(merged.harmony, file = "merged.harmony.rds")


hpca <- HumanPrimaryCellAtlasData()
dice <- DatabaseImmuneCellExpressionData()
hpca <- HumanPrimaryCellAtlasData()
merged <- readRDS("merged.harmony.rds")
names(merged[["RNA"]]@layers)
merged <- JoinLayers(object = merged, layers = "data")
pbmc_counts <- GetAssayData(merged, layer = "data")
# run SingleR
pred <- SingleR(test = pbmc_counts,
                ref = hpca,
                labels = hpca$label.main)
merged$singleR.labels <- pred$labels
DimPlot(merged, reduction = 'umap', group.by = 'singleR.labels', label = TRUE)



url <- "https://drive.google.com/uc?export=download&id=1ucKin0dHFmc3sdXcz3NR49dpR4DzlVCJ"
df <- read.delim(url, header = TRUE, row.names = 1, sep = "\t")
colnames(df)
# 1) If there’s a very specific long name you want to force to a given sample ID:
colnames(df)[colnames(df) == 
               "featureCounts.on.dataset.224.and.251..Counts_RNA.STAR.on.dataset.188..223..and.224..mapped.bam"] <- "T_iMK_D7_3"
colnames(df)

# 2) General clean-up:
#    - Remove everything up to and including 'Counts' plus separator ('.' or '_')
#    - Then remove any trailing RNA/STAR/mapped suffixes
names_clean <- colnames(df)
names_clean <- sub("^.*Counts[._]", "", names_clean)       # keep content AFTER Counts[._]
names_clean <- sub("_RNA.*$", "", names_clean)             # drop _RNA... and the rest
names_clean <- sub("\\.RNA.*$", "", names_clean)           # also handle '.RNA...' variants
names_clean <- sub("\\.STAR.*$", "", names_clean)          # drop .STAR... if present
names_clean <- sub("\\.mapped\\.bam$", "", names_clean)    # drop .mapped.bam if present

# 3) Optional: remove 'featureCounts...' prefix if somehow still present
names_clean <- sub("^featureCounts\\..*", "", names_clean)

# 4) Apply any known casing fixes (your example: T_IMK -> T_iMK)
names_clean <- sub("^T_IMK_", "T_iMK_", names_clean)

# Assign back
colnames(df) <- names_clean

# Remove version suffix (everything after the dot)
rownames(df) = sub("\\..*", "", rownames(df))


# Use the latest Ensembl release
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

IDs <- rownames(df)

# Query Ensembl
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = IDs,
  mart = mart
)

# Merge back to your data
mapping <- annot$hgnc_symbol[match(IDs, annot$ensembl_gene_id)]

# Replace NAs with original IDs
mapping[is.na(mapping) | mapping == ""] <- IDs[is.na(mapping) | mapping == ""]

# Make unique
mapping_unique <- make.unique(mapping)

# Set rownames
rownames(df) <- mapping_unique

rownames(df)

length(rownames(df))
raw_mapping <- annot$hgnc_symbol[match(IDs, annot$ensembl_gene_id)]


countData = as.matrix(sapply(df, as.numeric))
rownames(countData) = rownames(df)
# Keep genes with count > 3 in at least one sample
countData <- countData[rowSums(countData > 3) >= 3, , drop = FALSE]

# Check dimensions after filtering
dim(countData)

# Optional: confirm no all-zero rows remain
stopifnot(all(rowSums(countData) > 0))

mean(apply(countData, 1, function(x) median(x) == 0))

#check  if the numeric matrix is successfully created and ready for DEG
if(is.numeric(countData)){
  cat("The matrix is numeric and ready to use\n")
  cat("------------------------------------------------\n")
  cat("The summary of raw count data:\n\n")
  summary(countData)
} else{
  stop("The matrix must be numeric")
}



DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)  # good step before conversion

# Convert Seurat -> SCE
sce <- as.SingleCellExperiment(merged, assay = "RNA")

# Map your Seurat metadata into colData(sce).
# Replace these with the actual column names in merged@meta.data:
# e.g., merged$cell_type or merged$seurat_clusters for cell types,
#       merged$subject / donor_id for per-cell subject IDs.

colnames(merged@meta.data)


colData(sce)$cell_type <- as.factor(merged$singleR.labels)  # <-- change if name differs
colData(sce)$subject   <- as.factor(merged$orig.ident)    # <-- change if name differs

# Gene intersection between sc and bulk
common_genes <- intersect(rownames(sce), rownames(countData))
sce       <- sce[common_genes, ]
countData <- countData[common_genes, ]



# IMPORTANT: Do not pass sc.mtx or sc.cluster when using sc.sce

Idents(merged) <- merged$singleR.labels

markers_tbl <- FindAllMarkers(
  merged,
  only.pos = TRUE,        # keep positive markers
  min.pct  = 0.10,        # present in ≥10% of cells in the type
  logfc.threshold = 0.25  # adjust to your sensitivity
)

saveRDS(markers_tbl, file = "FindAllMarkers_output.rds")



# Optional: select top N markers per cell type by adjusted p-value or avg_log2FC
markers_tbl <- markers_tbl %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj) %>%  # prioritise strong and significant
  ungroup()

# Convert to a named list: names = cell types, values = gene vectors
markers_list <- split(markers_tbl$gene, markers_tbl$cluster)

# 2) Convert to SCE
sce <- as.SingleCellExperiment(merged, assay = "RNA")
colData(sce)$cell_type <- as.factor(merged$singleR.labels)
colData(sce)$subject   <- as.factor(merged$orig.ident)

# 3) Align gene sets with bulk (critical for marker matching)
common_genes <- intersect(rownames(sce), rownames(countData))
sce       <- sce[common_genes, ]
countData <- as.matrix(countData[common_genes, ])

# 4) Clean markers to the common gene universe (avoid mismatched IDs)
markers_list <- lapply(markers_list, function(g) intersect(g, common_genes))
# Optionally remove empty marker sets (cell types with 0 common markers)
markers_list <- markers_list[vapply(markers_list, length, 1L) > 0]

#pseudobulk analysis
#pseudobulk should be run on raw data
SeuratTD0$sample <- "Day0"
SeuratTD7$sample <- "Day7"
SeuratTD14$sample <- "Day14"
merged_pb <- merge(SeuratTD0, y = list(SeuratTD7, SeuratTD14))
# Normalize and find variable features
merged_pb <- NormalizeData(merged_pb)
merged_pb <- FindVariableFeatures(merged_pb)

# Scale and run PCA
merged_pb <- ScaleData(merged_pb)
merged_pb <- RunPCA(merged_pb)

# Build nearest-neighbor graph
merged_pb <- FindNeighbors(merged_pb, dims = 1:20)

# clustering
merged_pb <- FindClusters(merged_pb, resolution = 0.5)
# 4. Aggregate counts per sample × cluster
pb_counts <- AggregateExpression(         #retrieves raw data for DE and seurat clusters from normalized data
  merged_pb,
  group.by = c("sample","seurat_clusters"),
  assays = "RNA",
  slot = "counts"
)$RNA

# View results for Day7 vs Day0
res_D7 <- results(dds, name="condition_D7_vs_D0")

# View results for Day14 vs Day0
res_D14 <- results(dds, name="condition_D14_vs_D0")

# Inspect top DE genes
head(res_D7)
head(res_D14)

# Order by adjusted p-value
res_D7 <- res_D7[order(res_D7$padj), ]
res_D14 <- res_D14[order(res_D14$padj), ]
# MA plot: Day7 vs Day0
plotMA(res_D7, ylim = c(-5, 5), main = "Day7 vs Day0")

# MA ploy: Day14 vs Day0
plotMA(res_D14, ylim = c(-5, 5), main = "Day14 vs Day0")
# Convert results to data frame
df_D7 <- as.data.frame(res_D7)
df_D14 <- as.data.frame(res_D14)

# Volcano plot for Day14 vs Day0
ggplot(df_D7, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Day7 vs Day0",
       x = "Log2 Fold Change",
       y = "-log10(p-value)")
ggplot(df_D14, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Day14 vs Day0",
       x = "Log2 Fold Change",
       y = "-log10(p-value)")
# Variance stabilizing transformation
vsd <- vst(dds)

# PCA plot
plotPCA(vsd, intgroup = "condition")
# Install pheatmap from CRAN
install.packages("pheatmap")
# Select top 30 genes by fold change (Day14 vs Day0)
topgenes <- head(order(res_D7$padj, na.last = NA), 30)

# Extract normalized counts
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)

# Heatmap
pheatmap(mat, annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
         show_rownames = TRUE,
         main = "Top DE genes (Day7 vs Day0)")
# Select top 30 genes by fold change (Day14 vs Day0)
topgenes <- head(order(res_D14$padj, na.last = NA), 30)

# Extract normalized counts
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)

# Heatmap
pheatmap(mat, annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
         show_rownames = TRUE,
         main = "Top DE genes (Day14 vs Day0)")
res_D7  <- res_D7[order(res_D7$log2FoldChange, decreasing = TRUE), ]
res_D14 <- res_D14[order(res_D14$log2FoldChange, decreasing = TRUE), ]
head(res_D14, 20)   # top 20 genes at Day14
# Use BiocManager to install clusterProfiler
BiocManager::install("clusterProfiler")

# Load the package
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#downregulation enrichment at Day7
down_genes <- rownames(subset(res_D7, log2FoldChange < -1))
ego <- enrichGO(gene = down_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")

dotplot(ego, showCategory = 20)
#downregulation enrichment at Day14
down_genes <- rownames(subset(res_D14, log2FoldChange < -1))
ego <- enrichGO(gene = down_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")

dotplot(ego, showCategory = 20)

#upregulation enrichment at Day7
up_genes <- rownames(subset(res_D7, log2FoldChange > 1))
ego <- enrichGO(gene = up_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")

dotplot(ego, showCategory = 20)
#upregulation enrichment at Day14
up_genes <- rownames(subset(res_D14, log2FoldChange > 1))
ego <- enrichGO(gene = up_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")

dotplot(ego, showCategory = 20)

# 5. Build colData
meta <- merged_pb[[]]
meta$condition <- case_when(
  meta$sample == "Day0" ~ "D0",
  meta$sample == "Day7" ~ "D7",
  meta$sample == "Day14" ~ "D14"
)
meta$condition <- factor(meta$condition, levels=c("D0","D7","D14"))
meta$sample_cluster_id <- paste(meta$sample, meta$seurat_clusters, sep="_")
colData <- unique(meta[,c("sample","seurat_clusters","condition","sample_cluster_id")])
rownames(colData) <- colData$sample_cluster_id
colData <- colData[colnames(pb_counts),,drop=FALSE]

# 6. DESeq2

dds <- DESeqDataSetFromMatrix(countData=round(pb_counts),
                              colData=colData,
                              design=~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds)
res_D7  <- results(dds, name="condition_D7_vs_D0")
res_D14 <- results(dds, name="condition_D14_vs_D0")

# 5) Run MuSiC with markers (marker-based)
res <- music_prop(
  bulk.mtx = countData,
  sc.sce   = sce,
  clusters = "cell_type",
  samples  = "subject",
  markers  = markers_list   # <-- key change
)

write.csv(res$Est.prop.weighted, "MuSiC_celltype_proportions.csv")



prop_mat <- res$Est.prop.weighted
prop_mat <- as.matrix(prop_mat)


res_long <- reshape2::melt(
  prop_mat,
  varnames = c("Sample", "CellType"),
  value.name = "Proportion"
)

colnames(res_long) <- c("Sample", "CellType", "Proportion")
ggplot(res_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "Estimated Cell Proportion", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
res_long$TimePoint <- c("D0","D7","D14")[match(res_long$Sample, c("sample1","sample2","sample3"))]
ggplot(res_long, aes(x=TimePoint, y=Proportion, fill=CellType)) +
  geom_boxplot() +
  facet_wrap(~CellType) +
  theme_classic()



colnames(res_long)[colnames(res_long) == ""] <- "Unknown"
head(res_long)




# Derive TimePoint from Sample names that look like T_iMK_D7_3 etc., using base R
res_tp <- within(res_long, {
  TimePoint <- ifelse(grepl("^T_iMK_", Sample),
                      sub(".*(D[0-9]+).*", "\\1", Sample),
                      NA_character_)
})

# Keep only T_iMK samples and desired time points
res_tp <- subset(res_tp, grepl("^T_iMK_", Sample) & TimePoint %in% c("D0", "D7", "D14"))

# ---- ANOVA per CellType across D0/D7/D14 ----
anova_base <- function(df) {
  # Fit one-way ANOVA
  fit <- tryCatch(aov(Proportion ~ TimePoint, data = df), error = function(e) NULL)
  
  if (!is.null(fit)) {
    # Extract ANOVA p-value
    anova_p <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    # Tukey post-hoc comparisons
    tukey_res <- tryCatch(TukeyHSD(fit)$TimePoint, error = function(e) NULL)
    
    # Collect means per group
    means <- tapply(df$Proportion, df$TimePoint, mean, na.rm = TRUE)
    
    data.frame(
      CellType   = unique(df$CellType),
      mean_D0    = means["D0"],
      mean_D7    = means["D7"],
      mean_D14   = means["D14"],
      anova_p    = anova_p,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(CellType = unique(df$CellType),
               mean_D0 = NA, mean_D7 = NA, mean_D14 = NA,
               anova_p = NA)
  }
}

# Apply to each CellType
cell_types <- unique(res_tp$CellType)
anova_results_list <- lapply(cell_types, function(ct) {
  df_ct <- res_tp[res_tp$CellType == ct, ]
  anova_base(df_ct)
})

anova_results <- do.call(rbind, anova_results_list)

# Adjust p-values (Benjamini-Hochberg)
anova_results$p_adj <- p.adjust(anova_results$anova_p, method = "BH")

# View
head(anova_results, 20)

# Optionally save
write.csv(anova_results, "MuSiC_TiMK_timepoint_ANOVA.csv", row.names = FALSE)






















pheatmap(
  prop_mat,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("#313695", "#4575B4", "#74ADD1",
                             "#ABD9E9", "#E0F3F8", "#FFFFBF",
                             "#FEE090", "#FDAE61", "#F46D43",
                             "#D73027", "#A50026"))(100),
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Cell-type proportion heatmap (MuSiC, row-scaled)"
)
ggsave("MuSiC_heatmap_row_scaledr.png", width = 12, height = 6, dpi = 300)


# ---- 2) Stacked bar chart per sample ----
df_long <- reshape2::melt(prop_mat)
colnames(df_long) <- c("Sample", "CellType", "Proportion")

# (Optional) reorder samples (e.g., alphabetical or custom)
df_long$Sample <- factor(df_long$Sample, levels = sort(unique(df_long$Sample)))

# (Optional) order cell types by overall mean proportion
ct_order <- df_long %>% group_by(CellType) %>%
  summarize(mean_prop = mean(Proportion, na.rm = TRUE)) %>%
  arrange(desc(mean_prop)) %>%
  pull(CellType)
df_long$CellType <- factor(df_long$CellType, levels = ct_order)

# (Optional) a color palette for cell types
n_ct <- length(levels(df_long$CellType))
ct_colors <- scales::hue_pal()(n_ct)

# If sample names encode groups (e.g., Con, AZD, SM, T_iMK), derive a Condition column:
df_long <- df_long %>%
  mutate(
    Condition = case_when(
      grepl("^Con_", Sample) ~ "Control",
      grepl("^AZD_", Sample) ~ "AZD",
      grepl("^SM_",  Sample) ~ "SM",
      grepl("^T_iMK_", Sample) ~ "T_iMK",
      TRUE ~ "Other"
    )
  )

ggplot(df_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = ct_colors) +
  labs(
    title = "Cell-type composition by sample and condition",
    x = "Bulk sample",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  ) +
  facet_wrap(~ Condition, scales = "free_x", nrow = 1)

# ---- 4) Save outputs (optional) ----
ggsave("MuSiC_stacked_bar.png", width = 12, height = 6, dpi = 300)
ggsave("MuSiC_stacked_bar.pdf", width = 12, height = 6)

