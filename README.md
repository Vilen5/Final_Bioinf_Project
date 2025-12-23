OAS1 Emerges as a Critical Node in the Chemical Reprogramming of T-Cells to Megakaryocytes
Authors: Rawan Salah, Alhussein Ahmed, Nour A. Mohamed, Shahd Alaa, Maryam M. Anwar, Hana H. Khalaf, Youssef K. Attia

Affiliation: Biotechnology Biomolecular Chemistry Department, Faculty of Science, Cairo University

📌 Project Overview
This project investigates the computational evidence for chemically reprogramming human CD3+ T-cells into Induced Megakaryocytes (iMKs) and platelets using a five-molecule cocktail (5M). By utilizing integrative multi-omics analysis (Bulk and Single-cell RNA-Seq) and molecular dynamics simulations, we identified OAS1 as a potential critical node mediating this lineage transition.

📂 Repository Structure
Final_Report.pdf: The complete scientific report detailing the theoretical basis, methodology, and discussion of results.

Bulk RNA Seq.r: R script containing the pipeline for bulk RNA-Seq analysis, including preprocessing, Differential Expression Analysis (DESeq2), and Pathway Enrichment (GSEA/ORA).

scRNA-Seq Analysis.R: R script for single-cell RNA-Seq analysis, covering Quality Control, Dimensionality Reduction (PCA/UMAP), Clustering, Harmony Integration, and Cell Type Annotation.

🛠️ Methodology & Workflow
1. Bulk RNA-Seq Analysis
Data: 21 samples (Control, AZD4205, 4M, 4M+AZD4205, and 5M time-course at Day 0, 7, and 14).

Pipeline:

Alignment: RNA-STAR via Galaxy.

Quantification: FeatureCounts.

Differential Expression: DESeq2 (v1.50.2) to identify DEGs (|log2FC| > 1.0, p-adj < 0.05).

Visualization: PCA, Volcano plots (EnhancedVolcano), and Heatmaps (pheatmap).

Enrichment: ClusterProfiler and ReactomePA for ORA and GSEA.

2. Single-Cell RNA-Seq Analysis
Data: 3 samples (Day 0, Day 7, Day 14).

Pipeline:

QC & Filtering: Seurat package; filtering cells with high mitochondrial content or low feature counts.

Integration: Harmony to correct for batch effects.

Annotation: SingleR with Human Primary Cell Atlas (HPCA) reference.

Deconvolution: MuSiC2 to estimate cell proportions in bulk samples using scRNA-seq signatures.

3. Structural Bioinformatics
Target: OAS1 protein.

Tools: AutoDock Vina for molecular docking and GROMACS (v2024.4) for Molecular Dynamics simulations (50ns).

📦 Dependencies
To reproduce the analysis, ensure the following R packages are installed:

```
# Core Analysis
install.packages(c("ggplot2", "dplyr", "tidyr", "Seurat", "patchwork", "devtools"))

# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "EDASeq", "biomaRt", "clusterProfiler", 
                       "ReactomePA", "enrichplot", "org.Hs.eg.db", "SingleCellExperiment", 
                       "scDblFinder", "glmGamPoi"))

# Visualization
install.packages("pheatmap")
devtools::install_github("kevinblighe/EnhancedVolcano")
```
📊 Key Findings
Treatment Efficacy: The 5M cocktail induced the most significant expression shifts compared to 4M or AZD4205 alone.

Lineage Transition: Day 7 represents a critical "transition state," while Day 14 samples show clear upregulation of platelet-specific pathways (hemostasis, wound healing) and downregulation of T-cell immune markers.

OAS1 Role: Identified as a key protein in the transition; MD simulations confirm stable binding potential with cocktail ligands.

<img width="3000" height="3000" alt="Graphs full" src="https://github.com/user-attachments/assets/ce631806-aa48-47f6-819b-935dd8f4aa61" />
