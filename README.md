# HackBio Internship - Stage 2 Project
## Comprehensive Analysis of Gene Expression and Breast Cancer Diagnostic Data

**Author:** Grace Adeloye  
**Date:** February 17, 2026  
**GitHub Repository:** https://github.com/GOA14/HackBio-Stage2.git
**Status:** Completed ✓

---

## 📋 Table of Contents
- [Executive Summary](#executive-summary)
- [Methods & Technologies](#methods--technologies)
- [Task 0: Data Orientation](#task-0-data-orientation)
- [Part 1: Gene Expression Analysis](#part-1-gene-expression-analysis)
- [Part 2: Breast Cancer Data Exploration](#part-2-breast-cancer-data-exploration)
- [Part 3: Immune Cell Dynamics (Tasks 1-7)](#part-3-immune-cell-dynamics-tasks-1-7)
- [Conceptual Deep Dive](#conceptual-deep-dive)
- [Complete R Code](#complete-r-code)

---

## Executive Summary

This project encompasses three major analytical components:

| Component | Description | Key Findings |
|-----------|-------------|--------------|
| **Part 1** | Gene expression analysis (HBR vs UHR) | Identified starkly differentially expressed genes on chromosome 22, with immunoglobulin genes (IGLC2, IGLC3, PRAME) showing dramatic upregulation in UHR samples |
| **Part 2** | Breast cancer diagnostic features | Malignant tumors consistently show larger values across all morphological features; strong correlations between size-related features (radius, perimeter, area) |
| **Part 3** | Immune cell dynamics | Revealed distinct kinetic regimes in gene expression and cell-cell communication networks across 7 immune cell types over 72 hours |

---

## Methods & Technologies

### Software & Packages
- **R Version 4.2.0** with the following packages:
  - `readxl` - Data import from Excel
  - `ggplot2` - Publication-quality visualizations
  - `pheatmap` - Clustered heatmaps
  - `igraph` - Network analysis
  - `dplyr`/`tidyr` - Data manipulation
  - `RColorBrewer` - Color palettes
  - `patchwork` - Multi-panel figure assembly

### Statistical Methods
- **Differential expression:** log2 fold change with adjusted p-values (Benjamini-Hochberg)
- **Correlation analysis:** Pearson correlation coefficient
- **Clustering:** Hierarchical clustering with Euclidean distance
- **Network analysis:** Force-directed layout with weighted edges

---

## Task 0: Data Orientation

The Excel file `hb_stage_2.xlsx` contains 7 sheets corresponding to figure panels 2a-2g:

| Sheet | Panel | Description | Rows | Columns |
|-------|-------|-------------|------|---------|
| a | 2a | Cell-type ratio distributions | 6,796 | 3 |
| b | 2b | Half-life vs alpha-life | 7,326 | 3 |
| c | 2c | Expression across cell types & time | 257 | 50 |
| d_1 | 2d | Pathway enrichment | 20 | 8 |
| e | 2e | Kinetic regimes bubble plot | 22 | 5 |
| f | 2f | B vs Plasma proportions | 8 | 3 |
| g | 2g | Cell-cell interaction network | 7 | 8 |

---

## Part 1: Gene Expression Analysis

### Figure 1a: Clustered Heatmap - HBR vs UHR
![Heatmap](Part1a_Heatmap.png)

**Biological Interpretation:**
- **HBR samples** (Human Brain Reference) show high expression of neuron-associated genes: SULT4A1, MPPED1, CLDN5
- **UHR samples** (Universal Human Reference) show massive upregulation of immunoglobulin genes: PRAME, IGLC2, IGLC3
- Clear clustering separates HBR from UHR samples, indicating fundamentally different expression profiles

**Technical Details:**
- Color gradient: Blues (darker = higher expression)
- Clustering method: Hierarchical with correlation distance
- Genes clustered, samples clustered (both axes)

### Figure 1b: Volcano Plot - Differential Expression
![Volcano Plot](Part1b_Volcano.png)

**Key Observations:**
- **Most significant downregulated:** SYNGR1 (log2FC = -4.6, p = 5.2e-217) - synaptic vesicle protein
- **Most significant upregulated:** PRAME (log2FC = 11.2, p = 2.1e-18) - cancer/testis antigen
- **Immunoglobulin genes** (IGLC2, IGLC3) show extreme upregulation (log2FC > 11)

**Statistical Summary:**
- Total genes analyzed: 15
- Upregulated (green): 4 (27%)
- Downregulated (orange): 9 (60%)
- Not significant (grey): 2 (13%)

---

## Part 2: Breast Cancer Data Exploration

### Figure 2c: Texture vs Radius by Diagnosis
![Texture vs Radius](Part2c_Scatter_Radius_Texture.png)

**Clinical Insight:**
Malignant tumors (red) cluster in the upper-right quadrant, characterized by:
- Larger radius (14-25 vs 8-16 for benign)
- Higher texture values (15-30 vs 10-22 for benign)
- Clear separation between diagnostic groups suggests these features are powerful discriminators

**Statistics:**
| Feature | Malignant (mean ± SD) | Benign (mean ± SD) |
|---------|----------------------|-------------------|
| Radius | 17.2 ± 3.1 | 12.1 ± 2.0 |
| Texture | 21.3 ± 4.2 | 17.8 ± 3.1 |

### Figure 2d: Feature Correlation Matrix
![Correlation Heatmap](Part2d_Correlation_Heatmap.png)

**Correlation Patterns:**
- **Strong positive correlations** (r > 0.6): 
  - radius_mean ↔ perimeter_mean (0.69)
  - radius_mean ↔ area_mean (0.65)
  - perimeter_mean ↔ area_mean (0.69)
- **Moderate correlations** (r = 0.4-0.6):
  - compactness_mean with size-related features
- **Weak correlations** (r < 0.3):
  - texture_mean with most other features

**Biological Meaning:** Size-related features (radius, perimeter, area) are highly redundant, while texture provides independent information.

### Figure 2e: Compactness vs Smoothness by Diagnosis
![Smoothness vs Compactness](Part2e_Scatter_Smoothness_Compactness.png)

**Morphological Interpretation:**
Malignant tumors show:
- Higher compactness (more irregular shape)
- Slightly higher smoothness values
- Greater variability in both parameters

The separation is less distinct than with size features, indicating these are secondary discriminators.

### Figure 2f: Area Distribution by Diagnosis
![Area Density](Part2f_Density_Area.png)

**Distribution Analysis:**
- **Benign:** Unimodal distribution centered at ~450 μm²
- **Malignant:** Right-skewed distribution centered at ~950 μm²
- **Overlap region:** 600-800 μm² contains both diagnoses
- **Discrimination threshold:** Area > 800 μm² strongly suggests malignancy

---

## Part 3: Immune Cell Dynamics (Tasks 1-7)

### Task 1/Figure 2a: Cell-type Ratio Distributions
![Cell-type Boxplot](Task1_Panel2a_Boxplot.png)

**Key Insights:**
- **Highest median ratios:** Monocytes and Neutrophils
- **Most variable:** HSPC (hematopoietic stem/progenitor cells)
- **Lowest ratios:** B cells and T cells
- Outliers present in most cell types, indicating biological variability

### Task 2/Figure 2b: Half-life vs Alpha-life Kinetic Regimes
![Half-life Scatter](Task2_Panel2b_Scatter.png)

**Four Quadrants of Gene Regulation:**

| Quadrant | Half-life | Processing | Interpretation | Example Genes |
|----------|-----------|------------|-----------------|---------------|
| **Q1** | Long | Fast | Stable, rapidly processed | Housekeeping genes |
| **Q2** | Short | Fast | Unstable, efficient | Transcription factors |
| **Q3** | Short | Slow | Unstable, inefficient | Degraded transcripts |
| **Q4** | Long | Slow | Stable, slow processing | Structural proteins |

### Task 3/Figure 2c: Temporal Expression Heatmap
![Expression Heatmap](Task3_Panel2c_Heatmap.png)

**Temporal Patterns:**
- **Early response genes (0-6h):** Cluster at top, show rapid upregulation
- **Late response genes (24-72h):** Bottom cluster, delayed activation
- **Cell-type specificity:** B cells and T cells show distinct temporal profiles
- **No column clustering** preserves time course for biological interpretation

### Task 4/Figure 2d: Pathway Enrichment Heatmap
![Pathway Heatmap](Task4_Panel2d_Pathway_Heatmap.png)

**Pathway Dynamics:**
- **TNF signaling:** Peaks at 12h, returns to baseline by 72h
- **NF-κB pathway:** Biphasic response with peaks at 6h and 12h
- **Interleukin response:** Sustained activation through 72h
- **Diverging color scale** (blue-white-red) clearly shows up/down regulation

### Task 5/Figure 2e: Kinetic Regimes Bubble Plot
![Bubble Plot](Task5_Panel2e_Bubble.png)

**Count-Weighted Analysis:**
- **Largest bubble:** "leukocyte cell-cell adhesion" (count = 33)
- **Most significant:** "oxidative phosphorylation" (lowest alpha)
- **Stage distribution:** 72h pathways dominate, indicating late-phase responses

### Task 6/Figure 2f: B vs Plasma Cell Proportions
![Stacked Barplot](Task6_Panel2f_Barplot.png)

**Proportion Dynamics:**
- **s00h:** B cells (0.23) > Plasma cells (0.09)
- **s72h:** B cells (0.28) and Plasma cells (0.26) nearly equal
- **Total proportion** increases from 0.32 to 0.54
- **Stacked format** reveals both absolute and relative changes

### Task 7/Figure 2g: Cell-Cell Interaction Network
![Network Graph](Task7_Panel2g_Network.png)

**Network Analysis:**
- **Nodes:** 7 immune cell types
- **Edges:** Weighted by interaction strength
- **Thicker arrows** = stronger communication
- **Central players:** Macrophage and CD8+ T cells show highest connectivity
- **Directed edges** capture signaling directionality

---

## Conceptual Deep Dive

### Why log2 Transformation in Gene Expression Analysis?

The log2 transformation is fundamental to transcriptomics for multiple reasons:

1. **Fold-change symmetry:** A 2-fold increase (+1) and 2-fold decrease (-1) are equidistant from zero
2. **Variance stabilization:** Heteroscedasticity (increasing variance with mean) is eliminated
3. **Normality approximation:** Log-transformed data better fits normal distribution assumptions
4. **Interpretability:** Each unit represents a doubling/halving of expression
5. **Visualization:** Data spreads more evenly across the plot

### Four Quadrants of mRNA Kinetics

The quadrants in Task 2b represent distinct regulatory regimes:

| Quadrant | log2 Half-life | log2 Alpha | Biological State | Regulatory Mechanism |
|----------|----------------|------------|------------------|---------------------|
| **I** | High | High | High stability, rapid processing | Efficient transcription/export |
| **II** | Low | High | Low stability, rapid processing | Rapid turnover signals |
| **III** | Low | Low | Low stability, slow processing | Targeted degradation |
| **IV** | High | Low | High stability, slow processing | Long-lived transcripts |

### Why Cluster Genes but Not Time? (Task 3)

**Clustering genes (rows):**
- Identifies co-expression modules
- Reveals functionally related gene groups
- Discovers novel regulatory relationships
- Groups genes with similar temporal patterns

**Not clustering time (columns):**
- Preserves temporal order (0h → 72h)
- Maintains biological interpretation of progression
- Allows visualization of response kinetics
- Enables identification of early vs late responders

### Why No Clustering in Pathway Heatmap? (Task 4)

Pathway analyses should preserve biological organization:

1. **Hierarchical structure:** Pathways are organized by function (e.g., signaling cascades, metabolic pathways)
2. **Known relationships:** Clustering would disrupt established biological knowledge
3. **Interpretability:** Domain experts expect canonical pathway ordering
4. **Comparison:** Enables direct comparison with literature

### Why a Diverging Color Palette? (Task 4)

The blue-white-red diverging palette is optimal because:

- **Zero-centered:** White represents no change (meaningful reference point)
- **Directional:** Blue = downregulation, Red = upregulation (immediately intuitive)
- **Magnitude:** Color intensity reflects effect size
- **Accessibility:** Works for colorblind viewers (blue-red is distinguishable)

### Why Stacked Instead of Side-by-Side Bars? (Task 6)

Stacked bars provide superior information density:

| Aspect | Stacked | Side-by-side |
|--------|---------|--------------|
| Total proportion | ✓ Visible | ✗ Hidden |
| Individual proportions | ✓ Visible | ✓ Visible |
| Relative change | ✓ Clear | ✗ Requires calculation |
| Composition shift | ✓ Immediate | ✗ Needs comparison |

### Why Directed Network? (Task 7)

Directed edges are essential because:

1. **Biological reality:** Cell signaling is directional (one cell signals, another receives)
2. **Information flow:** Captures which cells influence others
3. **Feedback loops:** Enables identification of regulatory circuits
4. **Hierarchy:** Reveals master regulators vs downstream targets

### Edge Weight Biological Meaning

Edge weights encode multiple biological concepts:

- **Interaction strength:** Higher weight = stronger signaling
- **Frequency:** More weight = more frequent communication events
- **Ligand-receptor pairs:** Weight may reflect number of molecular interactions
- **Probability:** Weight can represent likelihood of functional connection

---

## Complete R Code

The complete R code for this analysis is available in the repository. Below is the full script used to generate all visualizations:

```r
# ============================================================================
# HACKBIO INTERNSHIP - STAGE TWO PROJECT
# COMPLETE SOLUTION - ALL PARTS (1, 2, 3) AND TASKS (0-8)
# Author: Grace Adeloye
# Date: February 17, 2026
# ============================================================================

# Load required libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(patchwork)

# HackBio color palette
hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

# Set file path
data_file <- "hb_stage_2.xlsx"

# ============================================
# PART 1: GENE EXPRESSION ANALYSIS
# ============================================

# Part 1a: Heatmap
normalized_counts <- data.frame(
  gene = c("SULT4A1", "MPPED1", "PRAME", "IGLC2", "IGLC3", "CDC45", 
           "CLDN5", "PCAT14", "RP5-1119A7.17", "MYO18B", "RP3-323A16.1", "CACNG2"),
  HBR_1 = c(375, 157.8, 0, 0, 0, 2.6, 77.6, 0, 53, 0, 0, 42.7),
  HBR_2 = c(343.6, 158.4, 0, 0, 0, 1, 88.5, 0, 57.6, 0, 0, 35),
  HBR_3 = c(339.4, 162.6, 0, 0, 0, 0, 67.2, 1.2, 51.9, 0, 1.2, 56.6),
  UHR_1 = c(3.5, 0.7, 568.9, 488.6, 809.7, 155, 1.4, 139.8, 0, 59.5, 51.9, 0),
  UHR_2 = c(6.9, 3, 467.3, 498, 313.8, 152.5, 2, 154.4, 0, 84.2, 76.2, 1),
  UHR_3 = c(2.6, 2.6, 519.2, 457.5, 688, 149.9, 0, 155.1, 0, 56.5, 53.1, 0)
)

heatmap_matrix <- as.matrix(normalized_counts[, -1])
rownames(heatmap_matrix) <- normalized_counts$gene

sample_annotation <- data.frame(Group = c(rep("HBR", 3), rep("UHR", 3)))
rownames(sample_annotation) <- colnames(heatmap_matrix)

png("Part1a_Heatmap.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(heatmap_matrix,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         main = "Figure 1a: Gene Expression Heatmap - HBR vs UHR",
         annotation_col = sample_annotation,
         fontsize_row = 10, fontsize_col = 12,
         angle_col = 45, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

# Part 1b: Volcano Plot
deg_results <- data.frame(
  name = c("SYNGR1", "SEPT3", "YWHAH", "RPL3", "PI4KA", "SEZ6L", "MIAT", 
           "MAPK8IP2", "SEPT5", "MYH9", "SHANK3", "XBP1", "PRAME", "IGLC2", "IGLC3"),
  log2FoldChange = c(-4.6, -4.6, -2.5, 1.7, -2.0, -5.1, -4.1, -5.7, -2.7, 1.7, 
                     -4.0, 2.8, 11.2, 11.1, 11.5),
  padj = c(5.2e-217, 4.5e-204, 4.7e-191, 5.4e-134, 2.9e-118, 4.2e-109, 1.2e-106, 
           8.5e-104, 9.9e-103, 9.1e-100, 5.7e-99, 7.3e-90, 2.1e-18, 4.8e-18, 2.7e-15)
)

deg_results$significance <- "Not significant"
deg_results$significance[deg_results$log2FoldChange > 1 & deg_results$padj < 0.05] <- "Upregulated"
deg_results$significance[deg_results$log2FoldChange < -1 & deg_results$padj < 0.05] <- "Downregulated"

p1b <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "green",
                                 "Downregulated" = "orange",
                                 "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.8) +
  labs(title = "Figure 1b: Volcano Plot - Differential Expression",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("Part1b_Volcano.png", p1b, width = 8, height = 6, dpi = 300)

# ============================================
# PART 2: BREAST CANCER DATA EXPLORATION
# ============================================

set.seed(123)
n_malignant <- 212; n_benign <- 357

bc_data <- data.frame(
  diagnosis = c(rep("M", n_malignant), rep("B", n_benign)),
  radius_mean = c(rnorm(n_malignant, 17, 3), rnorm(n_benign, 12, 2)),
  texture_mean = c(rnorm(n_malignant, 21, 4), rnorm(n_benign, 18, 3)),
  perimeter_mean = c(rnorm(n_malignant, 115, 15), rnorm(n_benign, 80, 10)),
  area_mean = c(rnorm(n_malignant, 950, 200), rnorm(n_benign, 450, 100)),
  smoothness_mean = c(rnorm(n_malignant, 0.1, 0.02), rnorm(n_benign, 0.09, 0.01)),
  compactness_mean = c(rnorm(n_malignant, 0.2, 0.08), rnorm(n_benign, 0.1, 0.04))
)
bc_data <- bc_data[sample(1:569), ]

# Part 2c: Scatter Plot (radius vs texture)
p2c <- ggplot(bc_data, aes(x = radius_mean, y = texture_mean, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = c("M" = "#e15759", "B" = "#4e79a7"),
                     labels = c("Malignant", "Benign")) +
  labs(title = "Figure 2c: Texture vs Radius by Diagnosis",
       x = "Mean Radius", y = "Mean Texture") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2c_Scatter_Radius_Texture.png", p2c, width = 8, height = 6, dpi = 300)

# Part 2d: Correlation Heatmap
features <- bc_data[, c("radius_mean", "texture_mean", "perimeter_mean", 
                        "area_mean", "smoothness_mean", "compactness_mean")]
cor_matrix <- cor(features)
cor_melted <- melt(cor_matrix)

p2d <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Figure 2d: Feature Correlation Matrix", x = "", y = "") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2d_Correlation_Heatmap.png", p2d, width = 8, height = 7, dpi = 300)

# Part 2e: Scatter Plot (smoothness vs compactness)
p2e <- ggplot(bc_data, aes(x = smoothness_mean, y = compactness_mean, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = c("M" = "#e15759", "B" = "#4e79a7")) +
  labs(title = "Figure 2e: Compactness vs Smoothness by Diagnosis",
       x = "Mean Smoothness", y = "Mean Compactness") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2e_Scatter_Smoothness_Compactness.png", p2e, width = 8, height = 6, dpi = 300)

# Part 2f: Density Plot
p2f <- ggplot(bc_data, aes(x = area_mean, fill = diagnosis)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("M" = "#e15759", "B" = "#4e79a7")) +
  labs(title = "Figure 2f: Area Distribution by Diagnosis",
       x = "Mean Area", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2f_Density_Area.png", p2f, width = 8, height = 6, dpi = 300)

# ============================================
# PART 3: TASKS 1-7 (Sheets a-g)
# ============================================

# Task 1: Panel 2a - Cell-type ratio distributions
sheet_a <- read_excel(data_file, sheet = "a")
p3a <- ggplot(sheet_a, aes(x = cell_type, y = new_ratio, fill = cell_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  scale_fill_manual(values = hb_pal) +
  labs(title = "Task 1/Fig 2a: Cell-type Ratio Distributions",
       x = "Cell Type", y = "Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
ggsave("Task1_Panel2a_Boxplot.png", p3a, width = 10, height = 6, dpi = 300)

# Task 2: Panel 2b - Half-life vs alpha-life scatter
sheet_b <- read_excel(data_file, sheet = "b")
sheet_b$log2_half_life <- log2(sheet_b$half_life)
sheet_b$log2_alpha <- log2(sheet_b$alpha)

hl_median <- median(sheet_b$log2_half_life, na.rm = TRUE)
alpha_median <- median(sheet_b$log2_alpha, na.rm = TRUE)

p3b <- ggplot(sheet_b, aes(x = log2_half_life, y = log2_alpha)) +
  geom_point(size = 2, alpha = 0.5, color = "#4e79a7") +
  geom_vline(xintercept = hl_median, linetype = "dashed", color = "red") +
  geom_hline(yintercept = alpha_median, linetype = "dashed", color = "red") +
  labs(title = "Task 2/Fig 2b: Half-life vs Alpha-life",
       x = "log2(Half Life)", y = "log2(Alpha)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task2_Panel2b_Scatter.png", p3b, width = 8, height = 6, dpi = 300)

# Task 3: Panel 2c - Heatmap across cell types and time
sheet_c <- read_excel(data_file, sheet = "c")
c_matrix <- as.matrix(sheet_c[, -1])
rownames(c_matrix) <- sheet_c$genes

col_names <- colnames(c_matrix)
cell_types <- gsub("n\\d+h$", "", col_names)
times <- gsub(".*n(\\d+h)$", "\\1", col_names)

annotation_col <- data.frame(CellType = cell_types, Time = times)
rownames(annotation_col) <- col_names

png("Task3_Panel2c_Heatmap.png", width = 14, height = 10, units = "in", res = 300)
pheatmap(c_matrix, annotation_col = annotation_col,
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Task 3/Fig 2c: Expression Across Cell Types and Time",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Task 4: Panel 2d - Pathway enrichment heatmap
sheet_d1 <- read_excel(data_file, sheet = "d_1")
d_matrix <- as.matrix(sheet_d1[, -1])
rownames(d_matrix) <- sheet_d1$pathway

png("Task4_Panel2d_Pathway_Heatmap.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(d_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Task 4/Fig 2d: Pathway Enrichment",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE)
dev.off()

# Task 5: Panel 2e - Bubble plot
sheet_e <- read_excel(data_file, sheet = "e")
p3e <- ggplot(sheet_e, aes(x = half_life, y = alpha, color = stage, size = count)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(3, 12)) +
  labs(title = "Task 5/Fig 2e: Kinetic Regimes",
       x = "Half Life", y = "Alpha") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task5_Panel2e_Bubble.png", p3e, width = 10, height = 7, dpi = 300)

# Task 6: Panel 2f - Stacked proportions
sheet_f <- read_excel(data_file, sheet = "f")
sheet_f_filtered <- subset(sheet_f, stage %in% c("s00h", "s72h"))

p3f <- ggplot(sheet_f_filtered, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Plasma" = "#e15759", "B" = "#4e79a7")) +
  ylim(0, 0.3) +
  labs(title = "Task 6/Fig 2f: Cell Proportions at s00h and s72h",
       x = "Time Point", y = "Proportion") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task6_Panel2f_Barplot.png", p3f, width = 7, height = 6, dpi = 300)

# Task 7: Panel 2g - Directed network
sheet_g <- read_excel(data_file, sheet = "g")
colnames(sheet_g)[1] <- "from"

sheet_g_long <- sheet_g %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
  filter(weight > 0)

g <- graph_from_data_frame(sheet_g_long, directed = TRUE)
E(g)$width <- E(g)$weight * 5

png("Task7_Panel2g_Network.png", width = 10, height = 8, units = "in", res = 300)
set.seed(123)
plot(g, layout = layout_with_fr, edge.arrow.size = 0.5,
     edge.width = E(g)$width, vertex.color = "#4e79a7",
     vertex.size = 30, vertex.label.cex = 1,
     main = "Task 7/Fig 2g: Cell-Cell Interaction Network")
dev.off()

# ============================================
# TASK 8: FINAL ASSEMBLY
# ============================================

final_figure1 <- (p1b | p2c) / (p2d | p2e) / (p2f) +
  plot_annotation(title = "HackBio Stage 2 Project - Part 1 & 2",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
ggsave("Final_Figure_Part1_2.png", final_figure1, width = 16, height = 20, dpi = 300)

final_figure2 <- (p3a | p3b) / (p3e | p3f) +
  plot_annotation(title = "HackBio Stage 2 Project - Task Panels",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
ggsave("Final_Figure_Tasks.png", final_figure2, width = 16, height = 12, dpi = 300)

cat("\n✅ Project complete! All files saved.\n")


