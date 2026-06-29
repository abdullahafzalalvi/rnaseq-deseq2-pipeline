pdf("results/plots/qc_sample_correlation_heatmap.pdf", width = 8, height = 7)
pheatmap(
sample_cor,
annotation_col    = annotation_col,
annotation_colors = ann_colors,
color             = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
main              = "Sample-to-Sample Pearson Correlation",
fontsize          = 10,
display_numbers   = TRUE,
number_format     = "%.3f"
)
dev.off()
cat("Sample correlation heatmap saved.\n")
# =============================================================================
# SECTION 6: PCA on Raw Log-Counts
# =============================================================================
pca_res <- prcomp(t(log_counts), scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample    <- rownames(pca_df)
pca_df$condition <- col_data$condition[match(pca_df$sample, rownames(col_data))]
var_explained <- unname(round(summary(pca_res)$importance[2, 1:2] * 100, 1))
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
geom_point(size = 5, alpha = 0.85) +
ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20) +
scale_color_manual(values = c(untrt = "#2196F3", trt = "#F44336")) +
labs(
title = "PCA of Raw Log2-Transformed Counts",
x     = paste0("PC1 (", var_explained[1], "% variance)"),
y     = paste0("PC2 (", var_explained[2], "% variance)"),
color = "Condition"
) +
theme_classic(base_size = 13) +
theme(legend.position = "right")
ggsave("results/plots/qc_pca_raw.pdf", p_pca, width = 7, height = 6)
ggsave("results/plots/qc_pca_raw.png", p_pca, width = 7, height = 6, dpi = 300)
# =============================================================================
# SECTION 7: Save Filtered Count Matrix and QC Summary
# =============================================================================
write.csv(count_filtered, "data/raw/counts_filtered.csv")
qc_summary <- data.frame(
Metric = c("Total genes (raw)", "Genes after filtering", "Genes removed",
"Number of samples", "Min library size", "Max library size"),
Value  = c(nrow(count_matrix), nrow(count_filtered),
nrow(count_matrix) - nrow(count_filtered),
ncol(count_filtered), min(lib_sizes), max(lib_sizes))
)
write.csv(qc_summary, "results/tables/qc_summary.csv", row.names = FALSE)
cat("\n=== QC Complete ===\n")
cat("Filtered count matrix saved to: data/raw/counts_filtered.csv\n")
cat("QC summary saved to: results/tables/qc_summary.csv\n")
cat("All QC plots saved to: results/plots/\n")
suppressPackageStartupMessages({
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(airway)
library(SummarizedExperiment)
library(tidyverse)
})
# =============================================================================
# SECTION 1: Load Data and Build DESeqDataSet
# =============================================================================
data("airway")
se <- airway
count_matrix <- assay(se, "counts")
col_data <- as.data.frame(colData(se))
col_data$condition <- col_data$dex  # "trt" vs "untrt"
# Pre-filter low-count genes
keep <- rowSums(count_matrix) >= 10
count_filtered <- count_matrix[keep, ]
# Build DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
countData = count_filtered,
colData   = col_data,
design    = ~ condition
)
# Set reference level (control group must be reference)
dds$condition <- relevel(dds$condition, ref = "untrt")
cat("DESeqDataSet created.\n")
cat(sprintf("Genes: %d | Samples: %d\n", nrow(dds), ncol(dds)))
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)
# =============================================================================
# SECTION 3: Extract Results
# =============================================================================
# Default contrast: trt vs untrt
res <- results(dds,
contrast      = c("condition", "trt", "untrt"),
alpha         = 0.05,          # FDR threshold for independent filtering
pAdjustMethod = "BH")          # Benjamini-Hochberg correction
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)
# Convert to data frame and sort by adjusted p-value
res_df <- as.data.frame(res) %>%
rownames_to_column("gene_id") %>%
arrange(padj) %>%
filter(!is.na(padj))
# Apply significance thresholds
res_df <- res_df %>%
mutate(
significance = case_when(
padj < 0.05 & log2FoldChange >  1.5 ~ "UP",
padj < 0.05 & log2FoldChange < -1.5 ~ "DOWN",
TRUE                                 ~ "NS"
)
)
cat(sprintf("\nTotal tested genes: %d\n", nrow(res_df)))
cat(sprintf("Significant DEGs (FDR < 0.05, |LFC| > 1.5): %d\n",
sum(res_df$significance != "NS")))
cat(sprintf("  Upregulated: %d\n",   sum(res_df$significance == "UP")))
cat(sprintf("  Downregulated: %d\n", sum(res_df$significance == "DOWN")))
# Save results
write.csv(res_df, "results/tables/DESeq2_results.csv", row.names = FALSE)
# Save significant DEGs only
sig_df <- res_df %>% filter(significance != "NS")
write.csv(sig_df, "results/tables/DESeq2_significant_DEGs.csv", row.names = FALSE)
cat("\nResults saved.\n")
# =============================================================================
# SECTION 4: Normalized Count Matrix (VST)
# =============================================================================
# Variance Stabilizing Transformation for visualization purposes
# VST is preferred over rlog for large datasets (>30 samples)
vst_data <- vst(dds, blind = FALSE)
# Save normalized counts
vst_matrix <- assay(vst_data)
write.csv(as.data.frame(vst_matrix), "results/tables/DESeq2_VST_normalized_counts.csv")
# =============================================================================
# SECTION 5: PCA on VST Data (Post-normalization)
# =============================================================================
pca_res <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
var_exp <- round(100 * attr(pca_res, "percentVar"), 1)
p_pca_vst <- ggplot(pca_res, aes(x = PC1, y = PC2, color = condition, label = name)) +
geom_point(size = 5, alpha = 0.85) +
ggrepel::geom_text_repel(size = 3.5) +
scale_color_manual(values = c(untrt = "#2196F3", trt = "#F44336")) +
labs(
title = "PCA of VST-Normalized Counts (DESeq2)",
x     = paste0("PC1: ", var_exp[1], "% variance"),
y     = paste0("PC2: ", var_exp[2], "% variance"),
color = "Condition"
) +
theme_classic(base_size = 13)
ggsave("results/plots/deseq2_pca_vst.pdf", p_pca_vst, width = 7, height = 6)
ggsave("results/plots/deseq2_pca_vst.png", p_pca_vst, width = 7, height = 6, dpi = 300)
# =============================================================================
# SECTION 6: MA Plot
# =============================================================================
pdf("results/plots/deseq2_ma_plot.pdf", width = 7, height = 5)
plotMA(res,
alpha = 0.05,
main  = "MA Plot — DESeq2 (trt vs untrt)",
colSig = "#F44336",
colNonSig = "grey60",
ylim  = c(-5, 5))
abline(h = c(-1.5, 1.5), lty = 2, col = "darkblue", lwd = 1.2)
dev.off()
# =============================================================================
# SECTION 7: Volcano Plot
# =============================================================================
p_volcano <- EnhancedVolcano(
res_df,
lab            = res_df$gene_id,
x              = "log2FoldChange",
y              = "padj",
pCutoff        = 0.05,
FCcutoff       = 1.5,
title          = "DESeq2 Differential Expression",
subtitle       = "trt vs untrt | FDR < 0.05 | |LFC| > 1.5",
xlab           = expression(paste("Log"[2], " Fold Change")),
ylab           = expression(paste("-Log"[10], " Adjusted p-value")),
col            = c("grey60", "grey40", "#2196F3", "#F44336"),
colAlpha       = 0.75,
legendLabels   = c("NS", "LFC only", "FDR only", "FDR + LFC"),
legendPosition = "right",
pointSize      = 2.5,
labSize        = 3.0,
drawConnectors = TRUE,
max.overlaps   = 20
)
ggsave("results/plots/deseq2_volcano.pdf", p_volcano, width = 10, height = 8)
ggsave("results/plots/deseq2_volcano.png", p_volcano, width = 10, height = 8, dpi = 300)
# =============================================================================
# SECTION 8: Heatmap of Top 50 DEGs
# =============================================================================
top50_genes <- sig_df %>%
arrange(padj) %>%
head(50) %>%
pull(gene_id)
top50_vst <- vst_matrix[top50_genes, , drop = FALSE]
# Z-score scaling per gene (row)
top50_scaled <- t(scale(t(top50_vst)))
annotation_col <- data.frame(
Condition = col_data$condition[match(colnames(top50_scaled), rownames(col_data))]
)
rownames(annotation_col) <- colnames(top50_scaled)
ann_colors <- list(Condition = c(untrt = "#2196F3", trt = "#F44336"))
pdf("results/plots/deseq2_heatmap_top50.pdf", width = 9, height = 14)
pheatmap(
top50_scaled,
annotation_col    = annotation_col,
annotation_colors = ann_colors,
color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
cluster_rows      = TRUE,
cluster_cols      = TRUE,
show_rownames     = TRUE,
show_colnames     = TRUE,
fontsize_row      = 7,
fontsize_col      = 9,
main              = "Top 50 DEGs — Z-Scored VST Counts (DESeq2)"
)
dev.off()
# =============================================================================
# SECTION 9: Dispersion Plot
# =============================================================================
pdf("results/plots/deseq2_dispersion.pdf", width = 7, height = 5)
plotDispEsts(dds,
main = "Gene-Wise Dispersion Estimates (DESeq2)",
cex  = 0.5)
dev.off()
# =============================================================================
# SECTION 10: Save Session Info for Reproducibility
# =============================================================================
sink("results/tables/DESeq2_session_info.txt")
cat("=== DESeq2 Analysis Session Info ===\n\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
sessionInfo()
sink()
cat("\n=== DESeq2 Analysis Complete ===\n")
cat("All outputs saved to results/\n")
suppressPackageStartupMessages({
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(airway)
library(SummarizedExperiment)
library(tidyverse)
})
# =============================================================================
# SECTION 1: Load Data
# =============================================================================
data("airway")
se <- airway
count_matrix <- assay(se, "counts")
col_data <- as.data.frame(colData(se))
col_data$condition <- col_data$dex
# Same pre-filter as DESeq2 for fair comparison
keep <- rowSums(count_matrix) >= 10
count_filtered <- count_matrix[keep, ]
# =============================================================================
# SECTION 2: Build DGEList and Normalize
# =============================================================================
group <- factor(col_data$condition[match(colnames(count_filtered), rownames(col_data))])
group <- relevel(group, ref = "untrt")
dge <- DGEList(counts = count_filtered, group = group)
# TMM normalization (Trimmed Mean of M-values)
# This corrects for differences in RNA composition between samples
dge <- calcNormFactors(dge, method = "TMM")
cat("=== edgeR DGEList Created ===\n")
cat(sprintf("Genes: %d | Samples: %d\n", nrow(dge), ncol(dge)))
cat("\nNormalization factors:\n")
print(dge$samples)
# =============================================================================
# SECTION 3: Design Matrix and Dispersion Estimation
# =============================================================================
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "trt_vs_untrt")
# Estimate dispersions (three steps):
#   1. Common dispersion
#   2. Trended dispersion
#   3. Tagwise (gene-wise) dispersion
dge <- estimateDisp(dge, design, robust = TRUE)
cat(sprintf("\nCommon BCV (Biological Coefficient of Variation): %.4f\n",
sqrt(dge$common.dispersion)))
# Save BCV plot
pdf("results/plots/edger_bcv_plot.pdf", width = 7, height = 5)
plotBCV(dge, main = "Biological Coefficient of Variation (edgeR)")
dev.off()
fit <- glmQLFit(dge, design, robust = TRUE)
# Save QLDisp plot
pdf("results/plots/edger_qldisp_plot.pdf", width = 7, height = 5)
plotQLDisp(fit, main = "Quasi-Likelihood Dispersions (edgeR)")
dev.off()
# Test the contrast: trt vs untrt
qlf <- glmQLFTest(fit, coef = "trt_vs_untrt")
cat("\n=== edgeR QLF Test Summary ===\n")
summary(decideTests(qlf, p.value = 0.05, lfc = 1.5))
# =============================================================================
# SECTION 5: Extract and Format Results
# =============================================================================
res_edger <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")
res_df_edger <- as.data.frame(res_edger$table) %>%
rownames_to_column("gene_id") %>%
rename(
log2FoldChange = logFC,
pvalue         = PValue,
padj           = FDR
) %>%
filter(!is.na(padj))
# Significance classification
res_df_edger <- res_df_edger %>%
mutate(
significance = case_when(
padj < 0.05 & log2FoldChange >  1.5 ~ "UP",
padj < 0.05 & log2FoldChange < -1.5 ~ "DOWN",
TRUE                                 ~ "NS"
)
)
cat(sprintf("\nTotal tested genes: %d\n", nrow(res_df_edger)))
cat(sprintf("Significant DEGs (FDR < 0.05, |LFC| > 1.5): %d\n",
sum(res_df_edger$significance != "NS")))
cat(sprintf("  Upregulated: %d\n",   sum(res_df_edger$significance == "UP")))
cat(sprintf("  Downregulated: %d\n", sum(res_df_edger$significance == "DOWN")))
write.csv(res_df_edger, "results/tables/edgeR_results.csv", row.names = FALSE)
sig_edger <- res_df_edger %>% filter(significance != "NS")
write.csv(sig_edger, "results/tables/edgeR_significant_DEGs.csv", row.names = FALSE)
# =============================================================================
# SECTION 6: Volcano Plot
# =============================================================================
p_volcano_edger <- EnhancedVolcano(
res_df_edger,
lab            = res_df_edger$gene_id,
x              = "log2FoldChange",
y              = "padj",
pCutoff        = 0.05,
FCcutoff       = 1.5,
title          = "edgeR Differential Expression",
subtitle       = "trt vs untrt | FDR < 0.05 | |LFC| > 1.5",
xlab           = expression(paste("Log"[2], " Fold Change")),
ylab           = expression(paste("-Log"[10], " Adjusted p-value")),
col            = c("grey60", "grey40", "#4CAF50", "#FF9800"),
colAlpha       = 0.75,
legendPosition = "right",
pointSize      = 2.5,
labSize        = 3.0,
drawConnectors = TRUE,
max.overlaps   = 20
)
ggsave("results/plots/edger_volcano.pdf", p_volcano_edger, width = 10, height = 8)
ggsave("results/plots/edger_volcano.png", p_volcano_edger, width = 10, height = 8, dpi = 300)
# =============================================================================
# SECTION 7: MD Plot (equivalent to MA plot in edgeR)
# =============================================================================
pdf("results/plots/edger_md_plot.pdf", width = 7, height = 5)
plotMD(qlf,
main   = "MD Plot — edgeR (trt vs untrt)",
status = decideTests(qlf, p.value = 0.05, lfc = 1.5),
values = c(-1, 0, 1),
col    = c("#2196F3", "grey70", "#F44336"),
cex    = 0.5)
abline(h = c(-1.5, 1.5), lty = 2, col = "darkgreen", lwd = 1.2)
dev.off()
# =============================================================================
# SECTION 8: Heatmap of Top 50 edgeR DEGs
# =============================================================================
# Get TMM-normalized log-CPM values
logcpm <- cpm(dge, log = TRUE, normalized.lib.sizes = TRUE)
top50_edger <- sig_edger %>%
arrange(padj) %>%
head(50) %>%
pull(gene_id)
top50_logcpm <- logcpm[top50_edger, , drop = FALSE]
top50_scaled_edger <- t(scale(t(top50_logcpm)))
annotation_col <- data.frame(
Condition = as.character(group[match(colnames(top50_scaled_edger),
rownames(col_data))])
)
rownames(annotation_col) <- colnames(top50_scaled_edger)
ann_colors <- list(Condition = c(untrt = "#2196F3", trt = "#F44336"))
pdf("results/plots/edger_heatmap_top50.pdf", width = 9, height = 14)
pheatmap(
top50_scaled_edger,
annotation_col    = annotation_col,
annotation_colors = ann_colors,
color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
cluster_rows      = TRUE,
cluster_cols      = TRUE,
show_rownames     = TRUE,
fontsize_row      = 7,
fontsize_col      = 9,
main              = "Top 50 DEGs — Z-Scored logCPM (edgeR)"
)
dev.off()
# =============================================================================
# SECTION 9: Session Info
# =============================================================================
sink("results/tables/edgeR_session_info.txt")
cat("=== edgeR Analysis Session Info ===\n\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
sessionInfo()
sink()
cat("\n=== edgeR Analysis Complete ===\n")
cat("All outputs saved to results/\n")
suppressPackageStartupMessages({
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)
library(RColorBrewer)
})
# =============================================================================
# SECTION 1: Load Results from Both Tools
# =============================================================================
deseq2_res <- read.csv("results/tables/DESeq2_results.csv")
edger_res  <- read.csv("results/tables/edgeR_results.csv")
# Significant DEG sets
deseq2_sig <- deseq2_res %>% filter(significance != "NS") %>% pull(gene_id)
edger_sig  <- edger_res  %>% filter(significance != "NS") %>% pull(gene_id)
cat("=== Cross-Tool DEG Count Summary ===\n")
cat(sprintf("DESeq2 significant DEGs: %d\n",  length(deseq2_sig)))
cat(sprintf("edgeR significant DEGs:  %d\n",  length(edger_sig)))
cat(sprintf("Consensus (overlap):     %d\n",  length(intersect(deseq2_sig, edger_sig))))
cat(sprintf("DESeq2 unique:           %d\n",  length(setdiff(deseq2_sig, edger_sig))))
cat(sprintf("edgeR unique:            %d\n",  length(setdiff(edger_sig, deseq2_sig))))
# =============================================================================
# SECTION 2: Venn Diagram of DEG Overlap
# =============================================================================
deg_list <- list(DESeq2 = deseq2_sig, edgeR = edger_sig)
p_venn <- ggVennDiagram(
deg_list,
label_alpha = 0,
set_color   = c("#1565C0", "#C62828")
) +
scale_fill_gradient(low = "#FFFFFF", high = "#EF9A9A") +
labs(title = "DEG Overlap: DESeq2 vs edgeR",
subtitle = "FDR < 0.05 | |LFC| > 1.5") +
theme(legend.title = element_blank())
ggsave("results/plots/venn_deseq2_edger_overlap.pdf", p_venn, width = 7, height = 6)
ggsave("results/plots/venn_deseq2_edger_overlap.png", p_venn, width = 7, height = 6, dpi = 300)
# =============================================================================
# SECTION 3: LFC Correlation Scatter Plot
# =============================================================================
# Merge on gene_id — only genes present in both result tables
merged <- inner_join(
deseq2_res %>% select(gene_id, log2FoldChange, padj) %>% rename(LFC_DESeq2 = log2FoldChange, padj_DESeq2 = padj),
edger_res  %>% select(gene_id, log2FoldChange, padj) %>% rename(LFC_edgeR  = log2FoldChange, padj_edgeR  = padj),
by = "gene_id"
)
merged <- merged %>%
mutate(
category = case_when(
gene_id %in% intersect(deseq2_sig, edger_sig) ~ "Consensus DEG",
gene_id %in% setdiff(deseq2_sig, edger_sig)   ~ "DESeq2 only",
gene_id %in% setdiff(edger_sig, deseq2_sig)   ~ "edgeR only",
TRUE                                            ~ "Non-significant"
)
)
# Calculate Pearson correlation
lfc_cor <- cor(merged$LFC_DESeq2, merged$LFC_edgeR, use = "complete.obs", method = "pearson")
p_scatter <- ggplot(merged, aes(x = LFC_DESeq2, y = LFC_edgeR, color = category)) +
geom_point(alpha = 0.5, size = 1.5) +
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
scale_color_manual(values = c(
"Consensus DEG"   = "#F44336",
"DESeq2 only"     = "#2196F3",
"edgeR only"      = "#4CAF50",
"Non-significant" = "grey80"
)) +
annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,
label = sprintf("Pearson r = %.4f", lfc_cor), size = 4.5) +
labs(
title    = "LFC Correlation: DESeq2 vs edgeR",
subtitle = "Each point represents one gene",
x        = expression(paste("Log"[2], " Fold Change (DESeq2)")),
y        = expression(paste("Log"[2], " Fold Change (edgeR)")),
color    = "Category"
) +
theme_classic(base_size = 13)
ggsave("results/plots/lfc_correlation_scatter.pdf", p_scatter, width = 8, height = 7)
ggsave("results/plots/lfc_correlation_scatter.png", p_scatter, width = 8, height = 7, dpi = 300)
cat(sprintf("\nPearson r (LFC correlation, DESeq2 vs edgeR): %.4f\n", lfc_cor))
# =============================================================================
# SECTION 4: Save Consensus DEG Table
# =============================================================================
consensus_genes <- intersect(deseq2_sig, edger_sig)
consensus_df <- merged %>%
filter(gene_id %in% consensus_genes) %>%
select(gene_id, LFC_DESeq2, padj_DESeq2, LFC_edgeR, padj_edgeR) %>%
mutate(
mean_LFC       = (LFC_DESeq2 + LFC_edgeR) / 2,
direction      = ifelse(mean_LFC > 0, "UP", "DOWN")
) %>%
arrange(padj_DESeq2)
write.csv(consensus_df, "results/tables/consensus_DEGs.csv", row.names = FALSE)
cat(sprintf("\nConsensus DEG table saved: %d genes\n", nrow(consensus_df)))
cat(sprintf("  Upregulated in consensus:   %d\n", sum(consensus_df$direction == "UP")))
cat(sprintf("  Downregulated in consensus: %d\n", sum(consensus_df$direction == "DOWN")))
# =============================================================================
# SECTION 5: Comparison Table for Manuscript
# =============================================================================
comparison_table <- data.frame(
Metric               = c("Total genes tested", "Significant DEGs",
"Upregulated", "Downregulated",
"Unique to tool", "Shared with other tool"),
DESeq2               = c(nrow(deseq2_res), length(deseq2_sig),
sum(deseq2_res$significance == "UP", na.rm = TRUE),
sum(deseq2_res$significance == "DOWN", na.rm = TRUE),
length(setdiff(deseq2_sig, edger_sig)),
length(intersect(deseq2_sig, edger_sig))),
edgeR                = c(nrow(edger_res), length(edger_sig),
sum(edger_res$significance == "UP", na.rm = TRUE),
sum(edger_res$significance == "DOWN", na.rm = TRUE),
length(setdiff(edger_sig, deseq2_sig)),
length(intersect(deseq2_sig, edger_sig)))
)
write.csv(comparison_table, "results/tables/tool_comparison_summary.csv", row.names = FALSE)
cat("\n=== Tool Comparison Table ===\n")
print(comparison_table)
cat("\nComparison table saved to results/tables/tool_comparison_summary.csv\n")
cat("\n=== Cross-Tool Comparison Complete ===\n")
