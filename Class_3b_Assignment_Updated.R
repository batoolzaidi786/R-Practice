# Assignment
# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

#### 0. Install and Load Packages ####

#Bioconductor provides R Packages for Omics Data Analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install Bioconductor Packages
BiocManager::install ("GEOquery" , force = TRUE) 
BiocManager::install ("affy" , force = TRUE) 
BiocManager::install ("arrayQualityMetrics" , force = TRUE) 



#Load required libraries
library(BiocGenerics)
?BiocGenerics
BiocManager:: install ("BiocGenerics" , force = TRUE)

library(GEOquery)

library(affy)

library(arrayQualityMetrics)

library(dplyr)

# GEOQuery is for Downloading GEO datasets (series matrix or raw CEL files)
# affy for preprocessing of the Affymatrix microarray data sets (RMA normalization)
# arrayQualityMetrics for QC reports of the microarray data (Raw or processed)
#dplyr is a CRAN package for data manipulation... graphs/figures etc



### Install Series Matrix files ####
#Preprocessed text files having
#expression values, sample annotations and probe info
#useful for quick analysis exploration when raw CEL Files are not needed
library(getGEO)

getGEO("GSE43322" , GSEMatrix = TRUE)
gse_data1 <- getGEO("GSE43322" , GSEMatrix = TRUE)

####Extract Expression Data Matrix (Gene/probes * samples) 
# Rows correspondes to probes and columns to samples

expression_data1 <- exprs(gse_data1$GSE43322_series_matrix.txt.gz)
?exprs

#### Extract Feature Data (probe annotation)
# Rows correspondes to probes and columns to samples
feature_data1 <- fData(gse_data1$GSE43322_series_matrix.txt.gz)

###Extract phenotype data (sample metadata)
# Rows correspondes to probes and columns to samples
phenotype_data1 <- pData(gse_data1$GSE43322_series_matrix.txt.gz)


#checking any missing values in sample annotations
sum(is.na(phenotype_data1$source_name_ch1))

####  Download Raw Files "(CEL Files) ####
# CEL files contain raw prob level intensity values for affymetrix platforms
#Raw data requires full preprocessing (RMA normalization, QC etc)

#CEL files are large and download may even fail with good connections so recommended to download directly

#Fetching the raw/CEL GEO supplementary files
getGEOSuppFiles("GSE43322" , baseDir = "Raw_Data1" , makeDirectory = TRUE)

#Important note: The pipline for preprocessing of Affymetrix microarray data is same whehter the raw CEL files are downloaded from NCBI GEO or ArrayExpress


#unzip files if as .zip

setwd("D:/BZ/AI for Omics Internship/Class3_Assignment")
file.exists("E-GEOD-43322.zip")

unzip("E-GEOD-43322.zip", exdir = "E_GEOD43322")


#Read CEL files  as Affybatch object
raw_data1 <- ReadAffy(celfile.path = "E-GEOD-43322")

raw_data1 #gives the basic info about the dataset. Notwdown the annotation (hgu133plus2) 
#from the output. you will need this for selection and installation correctly in next step
# the annotation package (hgu133plus2.db ) to map probe IDs to genes
dim(raw_data1)

#### QC before Pre-processing ####

#QC identifies outliers, hybridization problems and/or technical biases
#arrayQualityMetrics generates automated QC reports for microarray data
#It applies multiple complementary methods for detecting technical issues:
#Density pplots and Box plots: check distribution of intensities
# MA plots: visualize systemic biasis between arrays
#Heat maps and distance matrices: identify clustering/outliers
#PCA: detect unusual variations/ detecting outliers or batch effects

#Output is interactive html report index.html summarizing QC results

arrayQualityMetrics(expressionset = raw_data1 ,
                    outdir = "Results/QC_raw_data" ,
                    force = TRUE , 
                    do.logtransform = TRUE)

#RMA (Robust Multi Array Average ) Normalization
#RMA is a popular method for normalizing Affymetrix microarray data by: 
#1. background noise correcting
#2. normalizing probe intensities using quantile normalization
#3. summarizing them into gene level expression values using the robust median polish algorithm

#this method reduces experimetnal variation across multiple arrays
#Producing more symmetrical and normalized expression data
#compared to others

normalized_data1 <- rma(raw_data1)

#QC of the normalized data
arrayQualityMetrics(expressionset = normalized_data1 ,
                    outdir = "Results/QC_normalized_data" ,
                    force = TRUE)

#extract normalized expression values into a data frame
processed_Data1 <- as.data.frame(exprs(normalized_data1))

class(normalized_data1)
class(processed_Data1)

dim(processed_Data1) #Dimensions: no. of probes * samples

####Filter Low Variance Transcripts ("Soft" intensity based filtering) ####

#filtering removes low or uninformative expression signals probes
#Reasons: reduce noise and improves statistical power in differential expression
# or ML

#Calculate median intensity per probe across samples
Row_median1 <- rowMedians(as.matrix(processed_Data1))

#Visualize distribution of probe median intensity 
dev.new()      # Opens a new plotting device
hist(Row_median1,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
summary(Row_median1)

#Set a threshold value to remove low variance probes (dataset-specific, set accordingly)
threshold <- 5.5

abline(v = threshold , col = "Black" , lwd = 2) #abline a line created at the threshold value with color specified and width

# Select probes value above threshold

indx1 <- Row_median1 > threshold
filtered_data <- processed_Data1 [indx1, ]


#rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

#overwrite processed data with filtered data
processed_Data1 <- filtered_data

#Checking if subgroups withing PCOS groups (Treated and untreated) cluster together
#so that we can simplify groups as only PCOS vs. Control
library(limma)
exprs_mat <- exprs(normalized_data1)    # expression matrix: probes x samples
pca <- prcomp(t(exprs_mat), scale.=TRUE)
groups_detailed <- ifelse(grepl("Control", phenotype_data1$title, ignore.case = TRUE), "Control",
                          ifelse(grepl("LC n-3 PUFA", phenotype_data1$title, ignore.case = TRUE), "PCOS_LC_n3",
                                 ifelse(grepl("placebo", phenotype_data1$title, ignore.case = TRUE), "PCOS_placebo",
                                        ifelse(grepl("PCOS", phenotype_data1$title, ignore.case = TRUE), "PCOS_untreated", NA))))
groups_detailed <- factor(groups_detailed)

plot(pca$x[,1], pca$x[,2], col=as.factor(groups_detailed), pch=19)
legend("bottom", legend=levels(as.factor(groups_detailed)), col=1:length(levels(as.factor(groups_detailed))), pch=19)

# Phenotype data contains sample-level metadata such as condition, 
# tissue type, or disease status.
# Required to define experimental groups for statistical analysis.

class(phenotype_data$source_name_ch1) 

# Define experimental groups (normal vs PCOS)

simplified_group <- ifelse(grepl("Control", phenotype_data1$title, ignore.case = TRUE),
                       "Control",
                       ifelse(grepl("PCOS", phenotype_data1$title, ignore.case = TRUE),
                              "PCOS", NA))
simplified_group <- factor(simplified_group, levels = c("Control", "PCOS"))

groups <- factor(phenotype_data1$source_name_ch1,
                 levels = c("Control", "PCOS"),
                 labels = c("Control" , "PCOS"))
str(simplified_group)

class(groups)
levels(groups)

#Also for simplified groups 
groups <- ifelse(grepl("Control", phenotype_data1$source_name_ch1, ignore.case = TRUE),
                 "Control", "PCOS")

groups <- factor(groups, levels = c("Control", "PCOS"))

# Check what you got
table(groups)
levels(groups)

unique(phenotype_data1$source_name_ch1)
unique(phenotype_data1$title)

save.image(file = "Class3_Assignment.RData")
