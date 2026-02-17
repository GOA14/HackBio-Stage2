# HackBio Internship - Stage 2 Project

## Gene Expression and Breast Cancer Data Analysis

**Name:** GRACE ADELOYE
**Date:** February 17, 2026

## Project Overview
This project analyzes gene expression data (HBR vs UHR samples) and breast cancer diagnostic features using R. It includes heatmaps, volcano plots, scatter plots, correlation matrices, density plots, and network analysis.

## Results

### Part 1: Gene Expression Analysis

#### Figure 1a: Heatmap - HBR vs UHR
![Heatmap](Part1a_Heatmap.png)
*Clustered heatmap showing expression patterns of genes between HBR and UHR samples*

#### Figure 1b: Volcano Plot
![Volcano Plot](Part1b_Volcano.png)
*Volcano plot displaying log2 fold change vs -log10(adjusted p-value)*

### Part 2: Breast Cancer Data Exploration

#### Figure 2c: Texture vs Radius
![Texture vs Radius](Part2c_Scatter_Radius_Texture.png)
*Scatter plot showing relationship between texture and radius, colored by diagnosis*

#### Figure 2d: Feature Correlation Matrix
![Correlation Heatmap](Part2d_Correlation_Heatmap.png)
*Correlation matrix of six key features with annotated values*

#### Figure 2e: Smoothness vs Compactness
![Smoothness vs Compactness](Part2e_Scatter_Smoothness_Compactness.png)
*Relationship between smoothness and compactness by diagnosis*

#### Figure 2f: Area Distribution
![Area Density](Part2f_Density_Area.png)
*Density plot showing distribution of area values for malignant and benign diagnoses*

### Part 3: Additional Analyses (Tasks 1-7)

#### Task 1/Fig 2a: Cell-type Ratio Distributions
![Cell-type Boxplot](Task1_Panel2a_Boxplot.png)
*Boxplot showing distribution of ratios across immune cell types*

#### Task 2/Fig 2b: Half-life vs Alpha-life
![Half-life Scatter](Task2_Panel2b_Scatter.png)
*Scatter plot of log2(half_life) vs log2(alpha)*

#### Task 3/Fig 2c: Expression Heatmap
![Expression Heatmap](Task3_Panel2c_Heatmap.png)
*Heatmap showing gene expression across cell types and time points*

#### Task 5/Fig 2e: Kinetic Regimes Bubble Plot
![Bubble Plot](Task5_Panel2e_Bubble.png)
*Bubble plot showing half-life vs alpha, with size representing count*

#### Task 6/Fig 2f: Stacked Proportions
![Stacked Barplot](Task6_Panel2f_Barplot.png)
*Stacked bar plot showing cell proportions at s00h and s72h*

#### Task 7/Fig 2g: Cell-Cell Interaction Network
![Network Graph](Task7_Panel2g_Network.png)
*Directed network graph showing cell-cell interactions*

### Final Assembled Figures

#### Part 1 & 2 Combined
![Final Part 1-2](Final_Figure_Part1_2.png)
*Combined figure showing Part 1 and Part 2 results*

#### Tasks Combined
![Final Tasks](Final_Figure_Tasks.png)
*Combined figure showing Task panels*

## Conceptual Explanations

### Why log2 transformation?
- Makes multiplicative relationships additive
- Stabilizes variance across the range
- Makes data more normally distributed
- Log2 units are intuitive (doubling = +1, halving = -1)

### Four quadrants interpretation in Task 2b
- **Q1 (High/High)**: Long half-life, high processing rate - Stable but rapidly processed
- **Q2 (Low/High)**: Short half-life, high processing rate - Unstable but efficiently processed
- **Q3 (Low/Low)**: Short half-life, low processing rate - Unstable and slowly processed
- **Q4 (High/Low)**: Long half-life, low processing rate - Stable but slowly processed

### Why cluster genes but not time? (Task 3)
Genes are clustered to identify co-expression patterns and functional groups, while time points are kept in order to preserve temporal progression and biological interpretation.

### Why no clustering in pathway heatmap? (Task 4)
Pathways often have inherent biological organization that should be preserved. Clustering might disrupt known biological relationships.

### Why a diverging palette? (Task 4)
Diverging color scales (blue-white-red) show both positive and negative enrichment from a meaningful midpoint (zero), making directional changes immediately visible.

### Why stacked instead of side-by-side bars? (Task 6)
Stacked bars show both absolute and relative proportions simultaneously, making it easy to see how total composition changes over time.

### Why directed network? (Task 7)
Cell-cell interactions are inherently directional - signaling flows from one cell to another. Directed edges capture the biological reality of communication.

### Edge weight biological meaning (Task 7)
Edge weight represents interaction strength between cell types. Higher weights indicate stronger signaling, more frequent communication, or higher probability of interaction.

## Code Availability
The complete R code used for this analysis is available upon request.
