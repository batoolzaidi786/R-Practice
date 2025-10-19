# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================

# ---------------------------------------------------------------------
#              AI and Omics Research Internship (2025)
# ---------------------------------------------------------------------
#             Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------------------
#                     Microarray Data Analysis
# =====================================================================

# Topics:
# 1. Probe IDs to gene mapping
# 2. Differential Gene Expression Analysis
# 3. Data Visualization


gc()  # Clear memory to free up space before analysis

#### Install and Load Required Packages ####
# Check if BiocManager is installed; install if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages required for microarray analysis
BiocManager::install(c("limma", "AnnotationDbi", "hgu133plus2.db"))

# Each microarray platform has its own annotation package. 
# For example, for the Affymetrix Human Genome U133 Plus 2.0 Array, 
# the package name is 'hgu133plus2.db'.
# Find the appropriate annotation package for your platform 
# (e.g., hgu133a.db, hugene10sttranscriptcluster.db, etc.)
# and install/load it accordingly.

# Install CRAN packages for data manipulation and visualization
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

# Load Bioconductor packages
library(AnnotationDbi)   # Handles annotation and probe–gene mapping
library(hgu133plus2.db)  # Annotation database for Affymetrix HG-U133 Plus 2.0 array
library(limma)           # Performs linear modeling and differential expression
library(dplyr)           # Simplifies data manipulation tasks
library(tibble)          
library(ggplot2)         # Used for plotting and visualization
library(pheatmap)        # Generates heatmaps for gene expression data

# Note: dplyr and ggplot2 belong to the tidyverse collection

# -------------------------------------------------------------
#### Probe-to-Gene Mapping using AnnotationDbi ####
# -------------------------------------------------------------

# Load preprocessed expression and phenotype data
load("Class3_Assignment.RData")

# Alternative annotation resources include:
# Ensembl Genome Browser and DAVID Functional Annotation Tool

# check annotation slot of your dataset
annotation(raw_data1)
#My dataset (GSE43322) uses the NuGO_Hs1a520180 platform
#not the standard hgu133plus2 chip.

library(nugohs1a520180.db)
BiocManager::install("nugohs1a520180.db")

raw_data1

# Display objects available in the annotation package
ls("package:nugohs1a520180.db")

columns(nugohs1a520180.db)
keytypes(nugohs1a520180.db)

# Example of organism-specific annotation packages:
# org.Hs.eg.db for human
# org.Mm.eg.db for mouse

# -------------------------------------------------------------
# Extract probe IDs from processed microarray data
# -------------------------------------------------------------
probe_IDs <- rownames(processed_Data1)

# Map probe IDs to gene symbols using the platform annotation database
gene_symbols <- mapIds(
  nugohs1a520180.db,          # Database used for mapping
  keys = probe_IDs,        # Input probe IDs
  keytype = "PROBEID",     # Probe ID key type
  column = "SYMBOL",       # Desired annotation column (gene symbols)
  multiVals = "first"      # Return first match if multiple exist
)

# Convert mapping to a data frame and rename columns
Gene_Map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)

# -------------------------------------------------------------
# Handle multiple probes mapping to a single gene
# -------------------------------------------------------------
# Several strategies exist:
# 1. Retain probe with highest expression or variance
# 2. Average or summarize probe signals
# 3. Remove duplicate probes to maintain one row per gene

# Summarize number of probes per gene symbol
Duplicate_summary <- Gene_Map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# Identify genes associated with multiple probes
duplicate_genes <- Duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)
3511-1943
nrow(duplicate_genes)

# -------------------------------------------------------------
# Merge annotation mapping with expression data
# -------------------------------------------------------------

# Verify if probe IDs in mapping correspond to expression data
all(Gene_Map_df$PROBEID == row.names(processed_Data1))

# Merge annotation (SYMBOL) with expression matrix
processed_data_df <- processed_Data1 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

# Remove probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

# Select only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

# -------------------------------------------------------------
# Collapse multiple probes per gene using average expression
# -------------------------------------------------------------
# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

# Example to demonstrate how avereps works
x <- matrix(rnorm(8*3), 8, 3)
colnames(x) <- c("S1", "S2", "S3")
rownames(x) <- c("b", "a", "a", "c", "c", "b", "b", "b")
head(x)
avereps(x)  # Collapses duplicated row names by averaging

dim(averaged_data)

# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

# -------------------------------------------------------------
#### Differential Gene Expression Analysis ####
# -------------------------------------------------------------

# Define sample groups based on phenotype data
# Adjust group labels according to dataset annotation
groups <- factor(phenotype_data1$source_name_ch1,
                 levels = c("Control", "PCOS"),
                 label = c("Normal", "PCOS"))

class(groups)
levels(groups)

# -------------------------------------------------------------
# Create design matrix for linear modeling
# -------------------------------------------------------------
# Using no intercept (~0 + groups) allows each group to have its own coefficient
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit linear model to expression data
fit_1 <- lmFit(data, design)

# Define groups properly
groups <- factor(phenotype_data1$source_name_ch1,
                 levels = c("Subcutaneous adipose tissue; Control",
                            "Subcutaneous adipose tissue; PCOS",
                            "Subcutaneous adipose tissue; PCOS; placebo",
                            "Subcutaneous adipose tissue; PCOS; LC n-3 PUFA"),
                 labels = c("Control", "PCOS", "PCOS_placebo", "PCOS_treated"))

# Simplify to just Control vs PCOS if desired
groups_simple <- ifelse(grepl("Control", phenotype_data1$source_name_ch1), 
                        "Control", "PCOS")

groups_simple <- factor(groups_simple)

# Ensure sample order matches expression data
all(names(groups_simple) == colnames(data))  # should return TRUE

# Create design matrix
design <- model.matrix(~0 + groups_simple)
colnames(design) <- levels(groups_simple)

# Fit model
fit_1 <- lmFit(data, design)


# Define contrast to compare PCOS vs normal samples
contrast_matrix <- makeContrasts(PCOS_vs_Normal = PCOS - Control,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# -------------------------------------------------------------
# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
deg_results <- topTable(fit_2,
                        coef = "PCOS_vs_Normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

# -------------------------------------------------------------
# Classify DEGs into Upregulated, Downregulated, or Not Significant
# -------------------------------------------------------------
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))

# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

# Combine both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "Results/DEGs_Results.csv")
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Results/Updown_DEGs.csv")


# -------------------------------------------------------------
#### Data Visualization ####
# -------------------------------------------------------------

# -------------------------------------------------------------
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values
# -------------------------------------------------------------
# Note: x-axis = log2 fold change, y-axis = -log10 adjusted p-value

# Save volcano plot as PNG
png("Results/volcano_plot.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

dev.off()

# -------------------------------------------------------------
# Heatmap of Top Differentially Expressed Genes
# -------------------------------------------------------------

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

# Save heatmap as PNG
png("Heatmap_1st_25_DEGs.png", width = 2000, height = 1500, res = 300)

# Generate heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              # Disable row clustering
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 25 Differentially Expressed Genes"
)

dev.off()
save.image(file = "Class_3c_Assignment.RData")
