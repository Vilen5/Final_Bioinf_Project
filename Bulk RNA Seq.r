if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
if (!requireNamespace("EDASeq", quietly = TRUE)) {
  BiocManager::install("EDASeq")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("DESeq2")
}
if (!requireNamespace("apeglm", quietly = TRUE)) {
  install.packages("apeglm")
}
if (!requireNamespace("ggplotify", quietly = TRUE)) {
  install.packages("ggplotify")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  devtools::install_github("kevinblighe/EnhancedVolcano")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi", force = TRUE)
}
if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  BiocManager::install("ReactomePA", force = TRUE)
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot", force = TRUE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(biomaRt)
library(EDASeq)
library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggplotify)
library(patchwork)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)
library(ReactomePA)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

#The dataset is read directly from a Google Drive link, with the first column used as row names.
#Renames long, complicated featureCounts column names to simpler, interpretable names (, T_iMK_D7_3).
#Uses sub with a regex to systematically simplify other column names by extracting the sample identifier.
#Fixes specific column names for consistency and to avoid typos (T_IMK_D14_2 → T_iMK_D14_2).
url = "https://drive.google.com/uc?export=download&id=1ucKin0dHFmc3sdXcz3NR49dpR4DzlVCJ"
df = read.delim(url, header=TRUE, row.names=1)
colnames(df)[colnames(df) == "featureCounts.on.dataset.224.and.251..Counts_RNA.STAR.on.dataset.188..223..and.224..mapped.bam"] = "T_iMK_D7_3"
colnames(df) = sub('.*Counts[._](.*?)(?:_RNA.*|$)', '\\1', colnames(df))
colnames(df)[colnames(df) == "T_IMK_D14_2"] = "T_iMK_D14_2"


#Removes version suffixes from Ensembl gene IDs to ensure proper matching with annotation databases.
rownames(df) = sub("\\..*", "", rownames(df))


# Use the latest Ensembl release
#Connects to the latest Ensembl release for human genes.
#Retrieves hgnc_symbol (gene names) corresponding to each Ensembl gene ID.
mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

IDs = rownames(df)

# Query Ensembl
#Maps HGNC gene symbols back to the data. #Replaces missing symbols with original Ensembl IDs.
#Ensures all row names are unique, preventing issues in downstream analyses.
annot = getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = IDs,
  mart = mart
)

# Merge back to your data
mapping = annot$hgnc_symbol[match(IDs, annot$ensembl_gene_id)]

# Replace NAs with original IDs
mapping[is.na(mapping) | mapping == ""] = IDs[is.na(mapping) | mapping == ""]

# Make unique
mapping_unique = make.unique(mapping)

# Set rownames
#dfAll: contains samples excluding D0, D7, D14 (other treatments).
#df5M: contains only 5M treatment samples at D0, D7, D14.
rownames(df) = mapping_unique


raw_mapping = annot$hgnc_symbol[match(IDs, annot$ensembl_gene_id)]
dfAll = df[, !grepl("D0|D7|D14", colnames(df))]
df5M = df[, grepl("D0|D7|D14", colnames(df))]
#Converts data to numeric for differential expression analysis (DESeq2).
countDataAll = as.matrix(sapply(dfAll, as.numeric))
rownames(countDataAll) = rownames(dfAll)

# Keep genes with count > 3 in at least three samples
countDataAll = countDataAll[rowSums(countDataAll > 3) >= 3, , drop = FALSE]

#Checks to ensure no all-zero rows remain
stopifnot(all(rowSums(countDataAll) > 0))

#check  if the numeric matrix is successfully created and ready for DEG
if(is.numeric(countDataAll)){
  cat("The Different treatments matrix is numeric and ready to use\n")
  cat("------------------------------------------------\n")
  cat("The summary of Different treatments raw count data:\n\n")
  summary(countDataAll)
  cat("------------------------------------------------\n")
  cat("------------------------------------------------\n")
} else{
  stop("The Different treatments matrix must be numeric")
}
#Same process is repeated for countData5M (5M treatment across time points).
#Verifies that matrices are numeric and ready for downstream analysis.
countData5M = as.matrix(sapply(df5M, as.numeric))
rownames(countData5M) = rownames(df5M)
# Keep genes with count > 3 in at least three samples
countData5M = countData5M[rowSums(countData5M > 3) >= 3, , drop = FALSE]

# confirm no all-zero rows remain
stopifnot(all(rowSums(countData5M) > 0))

#check  if the numeric matrix is successfully created and ready for DEG
if(is.numeric(countData5M)){
  cat("\n\nThe different timepoints of 5M matrix is numeric and ready to use\n")
  cat("------------------------------------------------\n")
  cat("The summary of different timepoints of 5M raw count data:\n\n")
  summary(countData5M)
} else{
  stop("The different timepoints of 5M matrix must be numeric")
}

#ColData Preparation

#Lets prepare the sample names separately
samplesAll = colnames(countDataAll)
samples5M = colnames(countData5M)

#Check whether the sample column is created or not
#These checks confirm that the number of sample names matches the number of columns in the corresponding count matrices.
#This is a critical validation step to ensure that metadata and count data are aligned, preventing downstream errors in differential expression analysis.
if(length(samplesAll) == length(colnames(countDataAll))){
  cat("The samples for different treatments metadata have been generated successfully\n")
} else{
  stop("Failed generating all samples for different treatments from the countData")
}
if(length(samples5M) == length(colnames(countData5M))){
  cat("The samples for diffferent timepoints of 5M metadata have been generated successfully\n")
} else{
  stop("Failed generating all samples for different timepoints of 5M from the countData")
}

#Lets prepare the source column
#Assigns a biological source to each sample based on naming patterns
#For all treatments: samples are classified as Cultured T Cells.
#For 5M timepoints: samples are classified as CD3 or iMK. Converting these vectors to factors with predefined levels ensures consistent reference levels in statistical models.
sourcesAll = ifelse(grepl("^Con", samplesAll), "Cultured_T_Cells", "Cultured_T_Cells")
sourcesAll = factor(sourcesAll, levels = c("Cultured_T_Cells"))
sources5M = ifelse(grepl("^T_iMK_D0", samples5M), "CD3", "iMK")
sources5M = factor(sources5M, levels = c("CD3", "iMK"))


#Confirms that the source variables are correctly created as factors and match the sample length.
#This step ensures the experimental design is properly encoded.
if(is.factor(sourcesAll) && length(sourcesAll) == length(samplesAll)){
  cat("The sources for different treatments metadata have been generated successfully\n")
} else{
  stop("Failed generating all sources for different treatments from the countData")
}
if(is.factor(sources5M) && length(sources5M) == length(samples5M)){
  cat("The sources for different timepoints of 5M metadata have been generated successfully\n")
} else{
  stop("Failed generating all sources for different timepoints of 5M from the countData")
}

#lets prepare the groups
#For all treatments, samples are categorized into experimental groups such as Control, AZD4205, 4M_AZD4205, and 4M.
#For the 5M dataset, samples are grouped by time point (D0, D7, D14).

# Converting these to factors with ordered levels ensures correct baseline and contrasts in DEG analysis.
groupsAll = ifelse(grepl("^Con", samplesAll), "Control", ifelse(grepl("^AZD", samplesAll), "AZD4205", ifelse(grepl("^SM_AZD", samplesAll), "4M_AZD4205", "4M")))
groupsAll = factor(groupsAll, levels = c("Control", "AZD4205", "4M_AZD4205", "4M"))
groups5M = ifelse(grepl("D0", samples5M), "D0", ifelse(grepl("D7", samples5M), "D7", ifelse(grepl("D14", samples5M), "D14", "D7")))
groups5M = factor(groups5M, levels = c("D0", "D7", "D14"))

#Check if its prepared
#Confirms that all samples have been successfully assigned to groups.
#Protects against mislabeling or missing metadata.
if(is.factor(groupsAll) && length(groupsAll) == length(samplesAll)){
  cat("The groupes for different treatments metadata have been generated successfully\n")
} else{
  stop("Failed generating all groupes for different treatments from the countData")
}
if(is.factor(groups5M) && length(groups5M) == length(samples5M)){
  cat("The groupes for different timepoints of 5M metadata have been generated successfully\n")
} else{
  stop("Failed generating all groupes for different timepoints of 5M from the countData")
}


#Now, lets prepare the colData
#Constructs colData data frames required by DESeq2 and similar tools.
#Each row corresponds to a sample, and each column describes a biological or experimental variable.
#Row names are set to sample IDs to ensure exact matching with count matrix columns.
colDataAll = data.frame(sample = samplesAll, source = sourcesAll, group = groupsAll, row.names = samplesAll, stringsAsFactors = FALSE)
colData5M = data.frame(sample = samples5M, source = sources5M, group = groups5M, row.names = samples5M, stringsAsFactors = FALSE)

#Lets order the samples based on their groups
#reorders both the count matrices and metadata so samples are grouped logically.
#This improves clarity in downstream visualizations (PCA, heat maps, MA plots).
ordAll = order(colDataAll$group, colDataAll$sample)
countDataAll = countDataAll[, ordAll]
colDataAll   = colDataAll[ordAll, ]

ord5M = order(colData5M$group, colData5M$sample)
countData5M = countData5M[, ord5M]
colData5M   = colData5M[ord5M, ]

#Ensures perfect alignment between colData row names and count matrix column names.
#Confirms the expected dimensions of metadata
#These checks guarantee that the dataset is ready for differential expression analysis.
if(all(rownames(colDataAll) %in% colnames(countDataAll)) && all(rownames(colDataAll) == colnames(countDataAll)) && all(dim(colDataAll) == c(12, 3))){
  cat("The colData for different treatments is generated successfully\n")
} else{
  stop("Failed to generate the colData for different treatments")
}
if(all(rownames(colData5M) %in% colnames(countData5M)) && all(rownames(colData5M) == colnames(countData5M)) && all(dim(colData5M) == c(9, 3))){
  cat("The colData for different timepoints of 5M is generated successfully\n")
} else{
  stop("Failed to generate the colData for different timepoints of 5M")
}



#plotRLE() is a quality control plot that shows how consistent the expression distributions are between samples.
#countDataAll contains raw RNA-seq counts for different treatments.
#outline=FALSE , it removes outlier whiskers for cleaner visualization.
#Samples are colored according to their experimental group, allowing visual comparison between treatments. .
# main='Raw Counts - Different treatments' → title for the plot.
plotRLE(countDataAll, outline=FALSE,col= as.numeric(colDataAll$group),
        main = 'Raw Counts - Different treatments')

#Save the plot to the local files "Or colab runtime files"
#Opens a PNG device to save the plot as a high-resolution image.
#The dimensions (12 × 8 inches) and resolution (300 dpi) make the figure suitable for publications or reports.
png("RLE_Raw_Different_Treatments.png", width = 12, height = 8, units = "in", res = 300)
plotRLE(countDataAll, outline=FALSE,col= as.numeric(colDataAll$group),
        main = 'Raw Counts - Different treatments')
dev.off() #Finalizes and saves the PNG file to disk and Ensures the plot is written correctly and prevents graphics device conflicts.

plotRLE(countData5M, outline=FALSE,col= as.numeric(colData5M$group),
        main = 'Raw Counts - Different timepoints of 5M')

#Save the plot to the local files "Or colab runtime files"
png("RLE_Raw_Different_Times_5M.png", width = 12, height = 8, units = "in", res = 300)
plotRLE(countData5M, outline=FALSE,col= as.numeric(colData5M$group),
        main = 'Raw Counts - Different timepoints of 5M')
dev.off()


#Performs Principal Component Analysis (PCA) on the transposed count matrix.
#The transpose (t(countDataAll)) is necessary because PCA expects samples as rows and genes as columns — opposite of typical count matrices where genes are rows and samples are columns.
#PCA reduces the dimensionality of the dataset, summarizing the variation across thousands of genes into a few main components
pcaAll = prcomp(t(countDataAll))
pca5M = prcomp(t(countData5M))


#Calculates the percentage of variance explained by each principal component.
#pca$sdev gives the standard deviation for each component. Squaring and dividing by the total gives the proportion of total variance each component explains.
#This helps identify how much of the data’s variation is captured by PC1 and PC2.
percentVarAll = pcaAll$sdev^2 / sum(pcaAll$sdev^2)
percentVar5M = pca5M$sdev^2 / sum(pca5M$sdev^2)


#Creates a data frame containing the first two principal components (PC1 and PC2) for each sample, along with combined which that corresponding to group and conditions information
# each color represent unique treatment and each shape represent unique time of treatment

pcaDataAll = data.frame(
  PC1 = pcaAll$x[, 1],
  PC2 = pcaAll$x[, 2],
  group = colDataAll$group
  
)
pcaData5M = data.frame(
  PC1 = pca5M$x[, 1],
  PC2 = pca5M$x[, 2],
  group = colData5M$group
  
)

#This is used for visualization with ggplot2.
#Each point in the plot represents one sample positioned according to its gene expression profile.
#The xlab and ylab include the percentage of variance explained by each axis, giving context to the separation seen
ggplot(pcaDataAll, aes(x = PC1, y = PC2, color = group)) +
  geom_point( size = 3) +
  xlab(paste0("PC1: ", round(percentVarAll[1]*100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVarAll[2]*100), "% variance")) +
  theme_bw() + ggtitle("PCA Plot for different treatments - Raw") # title of PCA plot
ggsave("PCA_Raw_Different_Treatments.png", width = 12, height = 6, dpi = 300) # the dimensions of plot

ggplot(pcaData5M, aes(x = PC1, y = PC2, color = group)) +
  geom_point( size = 3) +
  xlab(paste0("PC1: ", round(percentVar5M[1]*100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar5M[2]*100), "% variance")) +
  theme_bw() + ggtitle("PCA Plot for different timepoints of 5M - Raw") # title of PCA plot
ggsave("PCA_Raw_Different_Times_5M.png", width = 12, height = 6, dpi = 300) # the dimensions of plot



#Sets the reference (baseline) groups for differential expression analysis:
#Control for all treatments.
#D0 for 5M time-course analysis.
#Identifies all comparison groups relative to these references, which will be tested automatically in the subsequent loops.

ref_groupAll  = "Control"
ref_group5M  = "D0"
treatmentsAll = setdiff(levels(colDataAll$group), ref_groupAll)
time5M = setdiff(levels(colData5M$group), ref_group5M)

#constructs DESeqDataSet objects that integrate which are Raw count matrices ,Sample metadata (colData), andExperimental design formula (~ group).
#These objects are the core input for DESeq2’s statistical modeling.
# this made for two analysis
ddsAll = DESeqDataSetFromMatrix(
  countData = countDataAll,
  colData   = colDataAll,
  design    = ~ group
)


dds5M = DESeqDataSetFromMatrix(
  countData = countData5M,
  colData   = colData5M,
  design    = ~ group
)

#Removes genes with very low total counts across all samples.
#This filtering step improves statistical power and reduces noise from uninformative features.
ddsAll = ddsAll[rowSums(counts(ddsAll)) > 3, ]
dds5M = dds5M[rowSums(counts(dds5M)) > 3, ]

#Performs the full DESeq2 workflow Library size normalization,Dispersion estimation,Fitting of negative binomial models,and Hypothesis testing for differential expression.
ddAll = DESeq(ddsAll)
dd5M = DESeq(dds5M)

#Retrieves normalized expression values, which are suitable for:
#Visualization (heat maps, PCA).
#Exploratory analysis.

#These counts are not used for DEG testing but for interpretation and plotting.
norm_countsAll = counts(ddAll, normalized = TRUE)
norm_counts5M = counts(dd5M, normalized = TRUE)

#Explicitly sets the reference levels to ensure that log2 fold changes are calculated relative to:
#Control for all treatments and D0 for the 5M time course.
#This guarantees biologically meaningful and consistent contrasts.

colDataAll$group = relevel(colDataAll$group, ref = ref_groupAll)
colData5M$group  = relevel(colData5M$group, ref = as.character(ref_group5M))

#Lists all coefficients generated by the DESeq2 model.
#Useful for validating contrast names and ensuring correct comparisons.
cat("Coefficient names for all treatments:\n")
print(resultsNames(ddAll))
cat("Coefficient names for all timepoints of 5M:\n")
print(resultsNames(dd5M))
#Iterates through each treatment group and compares it against Control.
res_listAll = list()

for (tr in treatmentsAll) {
  
  message("\nRunning: ", tr, " vs ", ref_groupAll)
  
  # Extract DESeq2 results
  resAll = results(ddAll,
                   contrast = c("group", tr, ref_groupAll),
                   alpha = 0.05)
  
  # Clean the results (remove rows with NA p-values)
  resAll = resAll[complete.cases(resAll$pvalue), ]
  
  # Order by P-values
  resAll = resAll[order(resAll$pvalue), ]
  
  # Save into list using a clean name
  name = paste0(tr, "_vs_", ref_groupAll)
  res_listAll[[name]] = resAll
  
  # Write CSV
  write.csv(as.data.frame(resAll),
            paste0("results_", name, ".csv"))
  
  message("Saved: results_", name, ".csv")
}
#Performs pairwise comparisons of D7 vs D0 and D14 vs D0.
res_list5M = list()
for (tr in time5M) {
  
  message("\nRunning: ", tr, " vs ", ref_group5M)
  
  # Extract DESeq2 results
  res5M = results(dd5M,
                  contrast = c("group", tr, ref_group5M),
                  alpha = 0.05)
  
  # Clean the results (remove rows with NA p-values)
  res5M = res5M[complete.cases(res5M$pvalue), ]
  
  # Order by P-values
  res5M = res5M[order(res5M$pvalue), ]
  
  # Save into list using a clean name
  name = paste0(tr, "_vs_", ref_group5M)
  res_list5M[[name]] = res5M
  
  # Write CSV
  write.csv(as.data.frame(res5M),
            paste0("results_", name, ".csv"))
  
  message("Saved: results_", name, ".csv")
}

vstAll = vst(ddAll, blind = TRUE)
vst5M = vst(dd5M, blind = TRUE)


#plotRLE() is a quality control plot that shows how consistent the expression distributions are between samples.
#as.matrix(countDataAll) → raw gene count matrix (each column = a sample).
#outline=FALSE → removes outlier whiskers for cleaner visualization.
#col=as.numeric(colDataAll$group) → each treatment have unique color.
# main='RLE Plot for different treatments - Raw Data' → title for the plot.
options(repr.plot.width = 14, repr.plot.height =8) # adjust the dimensions of the plot
par(mfrow = c(1, 2))
EDASeq::plotRLE(countDataAll, outline=FALSE, col= as.numeric(colDataAll$group), main = 'RLE Plot for different treatments - Raw Data')
#And now, for normalized data
#plotRLE() is a quality control plot that shows how consistent the expression distributions are between samples.
#as.matrix(norm_countsAll) → normalized gene count matrix (each column = a sample).
#outline=FALSE → removes outlier whiskers for cleaner visualization.
#col=as.numeric(colDataAll$group) → each treatment have unique color.
# main='RLE Plot for different treatments - Normalized Data' → title for the plot.
EDASeq::plotRLE(norm_countsAll, outline=FALSE, col= as.numeric(colDataAll$group), main = 'RLE Plot for different treatments - Normalized Data')

#Save the plot to the local files "Or colab runtime files"
png("RLE_Raw-Norm_Different_Treatments.png", width = 14, height = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
EDASeq::plotRLE(countDataAll, outline=FALSE, col= as.numeric(colDataAll$group), main = 'RLE Plot for different treatments - Raw Data')
EDASeq::plotRLE(norm_countsAll, outline=FALSE, col= as.numeric(colDataAll$group), main = 'RLE Plot for different treatments - Normalized Data')
dev.off()


#as.matrix(countData5M) → raw gene count matrix (each column = a sample).
#outline=FALSE → removes outlier whiskers for cleaner visualization.
#col=as.numeric(colData5M$group) → each timepoint have unique color.
# main='RLE Plot for different timepoints of 5M - Raw Data' → title for the plot.

options(repr.plot.width = 14, repr.plot.height =8) # adjust the dimensions of the plot
par(mfrow = c(1, 2))
EDASeq::plotRLE(countData5M, outline=FALSE, col= as.numeric(colData5M$group), main = 'RLE Plot for different timepoints of 5M - Raw Data')
#And now, for normalized data
#plotRLE() is a quality control plot that shows how consistent the expression distributions are between samples.
#as.matrix(norm_counts5M) → normalized gene count matrix (each column = a sample).
#outline=FALSE → removes outlier whiskers for cleaner visualization.
#col=as.numeric(colData5M$group) → each timepoint have unique color.
# main='RLE Plot for different treatments - Normalized Data' → title for the plot.
EDASeq::plotRLE(norm_counts5M, outline=FALSE, col= as.numeric(colData5M$group), main = 'RLE Plot for different timepoints of 5M - Normalized Data')

#Save the plot to the local files "Or colab runtime files"
png("RLE_Raw-Norm_Different_Times_5M.png", width = 14, height = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
EDASeq::plotRLE(countData5M, outline=FALSE, col= as.numeric(colData5M$group), main = 'RLE Plot for different timepoints of 5M - Raw Data')
EDASeq::plotRLE(norm_counts5M, outline=FALSE, col= as.numeric(colData5M$group), main = 'RLE Plot for different timepoints of 5M - Normalized Data')
dev.off()

# Performs Principal Component Analysis (PCA) using variance-stabilized transformed (VST) count data to explore global expression patterns across samples.

# Samples are grouped according to their experimental condition (diferent treatment) "group" to assess clustering and separation between treatments.

# Coloring samples by group allows visualization of treatment-

# A clean black-and-white theme is applied to enhance readability,
# descriptive title is added to clearly indicate to the plot
# PCA is based on normalized data from different treatments.

DESeq2::plotPCA(vstAll, intgroup = c("group")) +
  aes(color = group) +
  theme_bw(base_size = 12) +
  ggtitle("PCA for different treatments - Normalized")

#Save the plot to the local files "Or colab runtime files"
ggsave("PCA_Norm_Different_Treatments.png", width = 12, height = 6, dpi = 300)

DESeq2::plotPCA(vst5M, intgroup = c("group")) +
  aes(color = group) +
  theme_bw(base_size = 12) +
  ggtitle("PCA for different timepoints of 5M - Normalized")

#Save the plot to the local files "Or colab runtime files"
ggsave("PCA_Norm_Different_Times_5M.png", width = 12, height = 6, dpi = 300)


# plotDispEsts() generates a dispersion plot for the DESeq2 dataset (ddAll), which shows how gene-wise dispersion estimates relate to mean expression.
# The plot displays:
# - Black points: dispersion estimates for individual genes
# - A red line: the fitted dispersion trend across all genes
# the title of plot is Dispersion Plot for Different Treatments

plotDispEsts(ddAll, main="Dispersion Plot for Different Treatments")

#Save the plot to the local files "Or colab runtime files"
# Save the dispersion plot as a high-resolution PNG file
# Opens a PNG graphics device with defined width (12 inches), height (6 inches), and resolution (300 dpi) for high-quality export.
# The same dispersion plot is generated and written to "Dispersion_Treatments.png".
png("Dispersion_Treatments.png", width = 12, height = 6, units = "in", res = 300)
plotDispEsts(ddAll, main="Dispersion Plot for Different Treatments")
dev.off()

plotDispEsts(dd5M, main="Dispersion Plot for Different Timepoints of 5M")

#Save the plot to the local files "Or colab runtime files"
png("Dispersion_Times.png", width = 12, height = 6, units = "in", res = 300)
plotDispEsts(dd5M, main="Dispersion Plot for Different Timepoints of 5M")
dev.off()

#plotMA generate MA plot for diferent treatments
# ylim=c(-5,5), it adjust the dimensions of plot
#main="Pre-Shrinkage MA Plot for Different Treatments") is the title of plot
plotMA(ddAll, ylim=c(-5,5), main="Pre-Shrinkage MA Plot for Different Treatments")

#Save the plot to the local files "Or colab runtime files"
png("MA_Treatments.png", width = 12, height = 6, units = "in", res = 300)
plotMA(ddAll, ylim=c(-5,5), main="Pre-Shrinkage MA Plot for Different Treatments")
dev.off()

plotMA(dd5M, ylim=c(-5,5), main="Pre-Shrinkage MA Plot for different Timepoints of 5M")

#Save the plot to the local files "Or colab runtime files"
png("MA_Times.png", width = 12, height = 6, units = "in", res = 300)
plotMA(dd5M, ylim=c(-5,5), main="Pre-Shrinkage MA Plot for different Timepoints of 5M")
dev.off()



# Initialize lists to store results and MA plots
# res_shrink_listAll: will hold log2 fold change results after shrinkage for each treatment comparison.
# ma_plotsAll: will store MA plots corresponding to each treatment contrast for visualization.


# Identify coefficients for DESeq2 comparisons

# resultsNames(ddAll) returns all coefficients in the DESeq2 model.
# The first coefficient is typically the intercept, which is not used for treatment contrasts, so it is excluded with [-1].
# coef_namesAll now contains only the treatment vs Control contrasts.


# Loop through each treatment contrast

# For each coefficient (treatment comparison):
# Print a message indicating which contrast is currently being processed. This provides runtime feedback during long analyses.
# The loop will later be used to perform log2 fold change shrinkage (apeglm) and generate MA plots for each comparison.


res_shrink_listAll = list()
ma_plotsAll = list()
coef_namesAll = resultsNames(ddAll)[-1]
for (coef in coef_namesAll) {
  
  message("Processing: ", coef)
  
  # LFC shrinkage
  shr = lfcShrink(ddAll, coef = coef, type = "normal")
  res_shrink_listAll[[coef]] = shr
  
  # Convert plotMA to ggplot
  p = ggplotify::as.ggplot(
    ~ DESeq2::plotMA(
      shr,
      ylim = c(-5, 5),
      main = paste("Shrunken MA:", coef)
    )
  )
  
  ma_plotsAll[[coef]] = p
}

res_shrink_list5M = list()
ma_plots5M = list()
coef_names5M = resultsNames(dd5M)[-1]
for (coef in coef_names5M) {
  
  message("Processing: ", coef)
  
  # LFC shrinkage
  shr = lfcShrink(dd5M, coef = coef, type = "normal")
  res_shrink_list5M[[coef]] = shr
  
  # Convert plotMA to ggplot
  p = ggplotify::as.ggplot(
    ~ DESeq2::plotMA(
      shr,
      ylim = c(-5, 5),
      main = paste("Shrunken MA:", coef)
    )
  )
  
  ma_plots5M[[coef]] = p
}

# Combine individual MA plots into a single figure

# wrap_plots() from the patchwork package is used to combine all individual MA plots stored in 'ma_plotsAll' into a single composite figure.
# ncol = 2 specifies that the plots should be arranged in two columns, creating a clear and organized layout.

# Adjust plotting dimensions
# options(repr.plot.width, repr.plot.height) sets the dimensions of the combined figure when displayed.

# Display the combined MA plot figure

combined_maAll = wrap_plots(ma_plotsAll, ncol = 2)
combined_maAll
combined_ma5M = wrap_plots(ma_plots5M, ncol = 2)
combined_ma5M

#Save the plot to the local files "Or colab runtime files"
ggsave(
  filename = "Shrunken_MA_Different_Treatments.png",
  plot = combined_maAll,
  width = 12,
  height = 10,
  dpi = 300
)
ggsave(
  filename = "Shrunken_MA_Different_Times.png",
  plot = combined_ma5M,
  width = 12,
  height = 6,
  dpi = 300
)

# Clean and filter
# Removes rows with missing p-values or adjusted p-values (NA) to ensure only valid statistical results are considered.
res_cleanAll = resAll[!is.na(resAll$pvalue), ]
res_cleanAll = resAll[!is.na(resAll$padj), ]
res_clean5M = res5M[!is.na(res5M$pvalue), ]
res_clean5M = res5M[!is.na(res5M$padj), ]

# Filters for statistically significant DEGs (adjusted p-value < 0.05) and with meaningful expression changes (|log2FoldChange| > 1).
filtered_DEGsAll = res_cleanAll[res_cleanAll$padj < 0.05, ]
filtered_DEGs5M = res_clean5M[res_clean5M$padj < 0.05, ]
filtered_DEGsAll = filtered_DEGsAll[abs(filtered_DEGsAll$log2FoldChange) > 1, ]
filtered_DEGs5M = filtered_DEGs5M[abs(filtered_DEGs5M$log2FoldChange) > 1, ]

# Select top DEGs for visualization
# Limits the number of top DEGs to a maximum of 50 based on lowest p-values to focus on the most significant genes.
n_topAll = min(50, nrow(filtered_DEGsAll))
n_top5M = min(50, nrow(filtered_DEGs5M))

# Extract row names (gene identifiers) of the top DEGs ordered by statistical significance.
top_genesAll = rownames(filtered_DEGsAll)[order(filtered_DEGsAll$pvalue)][1:n_topAll]

top_genes5M = rownames(filtered_DEGs5M)[order(filtered_DEGs5M$pvalue)][1:n_top5M]

# Normalized counts for heatmap
top_countsAll = norm_countsAll[top_genesAll, ]
top_counts5M = norm_counts5M[top_genes5M, ]

# Subset DESeq2 results
res_topAll = resAll[top_genesAll, ]
res_top5M = res5M[top_genes5M, ]

#Group annotation for unique colors in heatmap
annotation_colAll = data.frame(Group = ddAll$group)
rownames(annotation_colAll) = colnames(ddAll)

annotation_col5M = data.frame(Group = dd5M$group)
rownames(annotation_col5M) = colnames(dd5M)

heatmapAll = ggplotify::as.ggplot(
  ~ pheatmap(
  top_countsAll,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_colAll,
  fontsize_row = 6,
  main = paste("Top", n_topAll, "DEGs with Sample Annotation")
))
ggsave(
  filename = "Heatmap_Different_Treatments.png",
  plot = heatmapAll,
  width = 12,
  height = 8,
  dpi = 300
)

# Initialize list to store volcano plots
# volcano_plotsAll will hold all volcano plots generated for each treatment vs Control comparison. Storing plots in a list allows easy manipulation, combination, and export of multiple figures.

# Loop through all DESeq2 result objects
# For each treatment comparison in 'res_listAll':
# Print a message indicating which comparison is being processed, providing runtime feedback during figure generation,Extract the DESeq2 results for the current comparison.
#   This includes log2 fold changes, p-values, and adjusted p-values,which are required to generate volcano plots.
volcano_plotsAll <- list()

for (name in names(res_listAll)) {
  
  message("Plotting volcano: ", name)
  
  res <- res_listAll[[name]]
  
  # remove NA padj / pvalue safely
  res <- res[!is.na(res$pvalue), ]
  
  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    title = name,
    subtitle = NULL,
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 0.2,
    labSize = 3.0,
    colAlpha = 1,
    axisLabSize = 8.0,
    legendLabSize = 5.0,
    legendIconSize = 3.0,
    legendPosition = "right"
  )
  
  volcano_plotsAll[[name]] <- p
}

# reapt pervioues codes but for 5M treatments in diferent time intervels
volcano_plots5M <- list()

for (name in names(res_list5M)) {
  
  message("Plotting volcano: ", name)
  
  res <- res_list5M[[name]]
  
  # remove NA padj / pvalue safely
  res <- res[!is.na(res$pvalue), ]
  
  #   - lab = labels each point with its gene name
  #   - x = uses log2FoldChange for the X-axis (magnitude of expression change)
  #   - y = uses pvalue for the Y-axis (statistical significance)
  #   - title = gives the plot a descriptive title
  #   - pCutoff = sets the p-value threshold for statistical significance
  #   - FCcutoff = sets the minimum fold change threshold to mark DEGs
  #   - pointSize = controls the size of each gene point
  #   - labSize = adjusts the text size of gene labels
  #   - colAlpha = sets the transparency of the points
  #   - ylim = limits the Y-axis range
  #   - axisLabSize = sets font size for axis labels
  #   - legendLabSize = adjusts font size for legend text
  #   - legendIconSize = adjusts the size of legend symbols
  #   - legendPosition = places the legend to the right of the plot
  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    title = name,
    subtitle = NULL,
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 0.5,
    labSize = 3.0,
    colAlpha = 1,
    axisLabSize = 8.0,
    legendLabSize = 5.0,
    legendIconSize = 3.0,
    legendPosition = "right"
  )
  
  volcano_plots5M[[name]] <- p
}

# Combine individual volcano plots into a single figure

# wrap_plots() from the patchwork package is used to combine all individual volcano plots stored in 'volcano_plotsAll' into a single composite figure.
# ncol = 2 specifies that the plots should be arranged in two columns, creating a clear and organized layout.

# Adjust plotting dimensions
# options(repr.plot.width, repr.plot.height) sets the dimensions of the combined figure when displayed.
combined_volcanoAll <- wrap_plots(volcano_plotsAll, ncol = 2)
combined_volcanoAll
combined_volcano5M <- wrap_plots(volcano_plots5M, ncol = 2)
combined_volcano5M

ggsave(
  filename = "Volcano_Different_Treatments.png",
  plot = combined_volcanoAll,
  width = 12,
  height = 10,
  dpi = 300
)
ggsave(
  filename = "Volcano_Different_Times.png",
  plot = combined_volcano5M,
  width = 12,
  height = 5,
  dpi = 300
)

# Generate a heatmap for the top DEGs across all treatments

# pheatmap() is used to visualize the expression patterns of
# the top DEGs (top_countsAll) across all samples.
# ggplotify::as.ggplot() converts the pheatmap object into a ggplot object, which allows further manipulation, integration with other ggplot figures, or saving with ggplot-compatible functions.
# Parameters:
# - scale = "row": standardizes each gene (row) to have mean 0 and standard deviation 1, highlighting relative expression differences across samples.
# - cluster_rows = TRUE: hierarchical clustering is applied to genes to group similar expression patterns together.
# - cluster_cols = TRUE: hierarchical clustering is applied to samples to identify patterns related to treatment groups.
# - annotation_col = annotation_colAll: adds color-coded annotations for sample groups, enhancing interpretability.
# - fontsize_row = 6: adjusts row label font size for readability.
# - main: provides a descriptive title indicating the number of top DEGs visualized.
heatmapAll = ggplotify::as.ggplot(
  ~ pheatmap(
    top_countsAll,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_colAll,
    fontsize_row = 6,
    main = paste("Top", n_topAll, "DEGs for different treatments")
  ))
ggsave(
  filename = "Heatmap_Different_Treatments.png",
  plot = heatmapAll,
  width = 12,
  height = 8,
  dpi = 300
)

heatmap5M = ggplotify::as.ggplot(
  ~ pheatmap(
    top_counts5M,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col5M,
    fontsize_row = 6,
    main = paste("Top", n_top5M, "DEGs for different timepoints of 5M")
  ))
ggsave(
  filename = "Heatmap_Different_Times.png",
  plot = heatmap5M,
  width = 12,
  height = 8,
  dpi = 300
)


## Enrichment Analysis
# In DESeq results, the genes were already mapped as gene symbols, and some of them was left as ENSEMBL IDs, so we had to reprocess it


all_results = c(res_listAll, res_list5M)

processed_results = list()


for (name in names(all_results)) {
  
  message("Processing: ", name)
  
  # Get all results from the list at ewach loop run
  res = all_results[[name]]
  ids = rownames(res)
  clean_ids = sub("\\.\\d+$", "", ids)  # Remove the suffix used to make genes mapped to the same id unique
  
  # Get the gene IDs
  res_df = as.data.frame(res)
  res_df$gene_ID = clean_ids
  
  
  # Take the original annotation vector used fot initial mapping as a reference
  res_annot = res_df %>%
    left_join(annot, by = c("gene_ID" = "hgnc_symbol")) %>%
    filter(!is.na(ensembl_gene_id))
  
  # Map the original ENSEMBL IDs to Entrez
  mapping = getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = res_annot$ensembl_gene_id,
    mart = mart
  )
  
  # Remove any genes with no ENTREZ
  mapping = mapping[!is.na(mapping$entrezgene_id), ]
  
  
  # Merge mapping results to the annotations
  res_annot = merge(res_annot, mapping, by = "ensembl_gene_id")
  
  # Remoce duplicate EntrezID
  res_annot = res_annot[!duplicated(res_annot$entrezgene_id), ]
  
  # We're only interested in the gene ids, log FC and Padj columns
  res2 = res_annot[c("entrezgene_id", "log2FoldChange", "padj")]
  names(res2)[1] = "entrezID"
  
  # Store the results of each tun to a list
  processed_results[[name]] = res2
}

#Make lists to store the results
ORA_up_list = list()
ORA_down_list = list()
dotplots_up = list()
dotplots_down = list()

# Loop to produce dot plots
for (name in names(processed_results)) {
  
  message("Processing: ", name)
  
  # Get all results from the list at ewach loop run
  res = processed_results[[name]]
  
  # Get downregulated and up regulated genes
  down = res %>% filter(!is.na(padj) & padj < 0.05 & log2FoldChange < 0)
  up   = res %>% filter(!is.na(padj) & padj < 0.05 & log2FoldChange > 0)
  
  #Run ORA with reactome for upregulated genes
  # Reactome was chosen as its the most suitable for blood samples and pathways
  if(nrow(up) > 0){
    ORA_Reactome_Up = enrichPathway(gene = up$entrezID,
                                    organism = "human",
                                    readable = TRUE,
                                    qvalueCutoff = 0.05)
    # Store ORA results in its list
    ORA_up_list[[name]] = ORA_Reactome_Up
    
    #Plot the dot plot
    updot = dotplot(ORA_Reactome_Up,showCategory=20,font.size=10,label_format=70)+
      scale_size_continuous(range=c(1, 7))+
      theme_minimal() +
      ggtitle(paste0("Reactome Up-regulated: ", name))
    #Save the plot
    ggsave(paste0("Dot_Up_", name, ".png"), plot = updot, width = 10, height = 8, dpi = 300)
    dotplots_up[[name]] = updot
  }
  #The same for the downregulated
  if(nrow(down) > 0){
    ORA_Reactome_Down = enrichPathway(gene = down$entrezID,
                                      organism = "human",
                                      readable = TRUE,
                                      qvalueCutoff = 0.05)
    ORA_down_list[[name]] = ORA_Reactome_Down
    downdot = dotplot(ORA_Reactome_Down,showCategory=20,font.size=10,label_format=40)+
      scale_size_continuous(range=c(1, 7))+
      theme_minimal() +
      ggtitle(paste0("Reactome Down-regulated: ", name))
    ggsave(paste0("Dot_Down_", name, ".png"), plot = downdot, width = 10, height = 8, dpi = 300)
    dotplots_down[[name]] = downdot
  }
}

#Make lists to store the results
GSEA_list <- list()
gseaplots <- list()

#Loop to produce GSEA plots
for (name in names(processed_results)) {
  
  message("Processing: ", name)
  
  # Get all results from the list at ewach loop run
  res <- processed_results[[name]]
  
  # make a vector of the log FC with the entrez ID as names
  entire_gene_list= res %>% pull(log2FoldChange,name="entrezID")
  
  # Sort the vector
  entire_gene_list <- sort(entire_gene_list, decreasing = TRUE)
  
  #Set seeds for reproducibility as gsePathway can go through randomization
  set.seed(123)
  GSEA_Reactome  = gsePathway(entire_gene_list,
                              organism="human",
                              seed=TRUE,
                              pvalueCutoff = 0.05)
  
  #Search for the most relevant pathway
  rank_1= which(GSEA_Reactome$Description == "Factors involved in megakaryocyte development and platelet production")
  rank_2= which(GSEA_Reactome$Description == "Kinesins")
  
  #If both pathways are present, plot both of them in one plot
  if(length(rank_1) > 0 && length(rank_2) > 0){
    GSEA_list[[name]] = GSEA_Reactome
    gseaplot = enrichplot::gseaplot2(GSEA_Reactome,
                                     geneSetID = c(rank_1, rank_2),
                                     title = paste0("GSEA Reactome: ", name))
    ggsave(paste0("Enrichment_", name, ".png"), plot = gseaplot, dpi = 300)
    gseaplots[[name]] = gseaplot
    
    # If only one of them is present, plot it on its own
  } else if(length(rank_1) > 0 && length(rank_2) == 0){
    GSEA_list[[name]] = GSEA_Reactome
    gseaplot = enrichplot::gseaplot2(GSEA_Reactome,
                                     geneSetID = rank_1,
                                     title = paste0("GSEA Reactome: ", name, " ", GSEA_Reactome$Description[[rank_1]]))
    ggsave(paste0("Enrichment_", name, ".png"), plot = gseaplot, dpi = 300)
    gseaplots[[name]] = gseaplot
  } else if(length(rank_1) == 0 && length(rank_2) > 0){
    GSEA_list[[name]] = GSEA_Reactome
    gseaplot = enrichplot::gseaplot2(GSEA_Reactome,
                                     geneSetID = rank_2,
                                     title = paste0("GSEA Reactome: ", name, " ", GSEA_Reactome$Description[[rank_2]]))
    ggsave(paste0("Enrichment_", name, ".png"), plot = gseaplot, dpi = 300)
    gseaplots[[name]] = gseaplot
    
    #If both of them are absent, get the top ranked pathway
  } else if(length(rank_1) == 0 && length(rank_2) == 0){
    GSEA_list[[name]] = GSEA_Reactome
    gseaplot = enrichplot::gseaplot2(GSEA_Reactome,
                                     geneSetID = 1,
                                     title = paste0("GSEA Reactome: ", name, " ", GSEA_Reactome$Description[[1]]))
    ggsave(paste0("Enrichment_", name, ".png"), plot = gseaplot, dpi = 300)
    gseaplots[[name]] = gseaplot
  }
}
