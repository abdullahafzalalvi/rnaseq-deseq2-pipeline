# =============================================================================
# Script: 00_install_packages.R
# Purpose: Install all required packages for the RNA-seq DE pipeline
# Author: Abdullah Afzal Alvi
# =============================================================================

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define required packages
cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "ggVennDiagram",
  "scales",
  "gridExtra",
  "corrplot",
  "here"
)

bioc_packages <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "clusterProfiler",
  "EnhancedVolcano",
  "org.Hs.eg.db",       # Human annotation (swap for org.At.tair.db for Arabidopsis)
  "AnnotationDbi",
  "airway",             # Example dataset
  "SummarizedExperiment",
  "BiocParallel"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
installed_cran <- rownames(installed.packages())
for (pkg in cran_packages) {
  if (!pkg %in% installed_cran) {
    install.packages(pkg, dependencies = TRUE)
    cat(sprintf("  Installed: %s\n", pkg))
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
installed_bioc <- rownames(installed.packages())
for (pkg in bioc_packages) {
  if (!pkg %in% installed_bioc) {
    BiocManager::install(pkg, ask = FALSE)
    cat(sprintf("  Installed: %s\n", pkg))
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

cat("\nAll packages installed. Session info:\n")
sessionInfo()
