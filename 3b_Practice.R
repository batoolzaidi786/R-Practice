#Module 2: Intro to Genomics Data Analysis
#Microarray Data Analysis
###Quality Control
###RMA Normalization
###Pre-Prcoessing and Filtering

#### 0. Install and Load Packages ####

#Bioconductor provides R Packages for Omics Data Analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install Bioconductor Packages
BiocManager::install (c("GEOquery" , "affy" , "arrayQualityMetrics"))

 #Load required libraries
# GEOQuery is for Downloading GEO datasets (series matrix or raw CEL files)
library(GEOquery)

# affy for preprocessing of the Affymatrix microarray data sets (RMA normalization)
library(affy)
# arrayQualityMetrics for QC reports of the microarray data (Raw or processed)
library(arrayQualityMetrics)

#dplyr is a CRAN package for data manipulation... graphs/figures etc
library(dplyr)

### Install Series Matrix files ####
#Preprocessed text files having
#expression values, sample annotations and probe info
#useful for quick analysis exploration when raw CEL Files are not needed

getGEO("GSE79973" , GSEMatrix = TRUE)
gse_data <- getGEO("GSE79973" , GSEMatrix = TRUE)

####Extract Expression Data MAtrix (Gene/probes * samples) 
# Rows correspondes to probes and columns to samples

expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)
?exprs

#### Extract Feature Data (probe annotation)
# Rows correspondes to probes and columns to samples
feature_data <- fData(gse_data$GSE79973_series_matrix.txt.gz)

###Extract phenotype data (sample metadata)
# Rows correspondes to probes and columns to samples
phenotype_data <- pData(gse_data$GSE79973_series_matrix.txt.gz)


#checking any missing values in sample annotations
sum(is.na(phenotype_data$source_name_ch1))

####  Download Raw Files "(CEL Files) ####
# CEL files contain raw prob level intensity values for affymetrix platforms
#Raw data requires full preprocessing (RMA normalization, QC etc)

#CEL files are large and download may even fail with good connections so recommended to download directly

#Fetching the raw/CEL GEO supplementary files
getGEOSuppFiles("GSE79973" , baseDir = "Raw_Data" , makeDirectory = TRUE)

#Important note: The pipline for preprocessing of Affymetrix microarray data is same whehter the raw CEL files are downloaded from NCBI GEO or ArrayExpress

#untar files if compressed as .tar 
untar("Raw_Data/GSE79973_RAW.tar" , exdir = "Raw_Data/CEL_Files")
#unzip files if as .zip
unzip("Raw_Data/E-GEOD-79973.zip", exdir = "Raw_Data/E_GEOD79973")

file.exists("Raw_Data/E-GEOD-79973.zip")
#So we dont have E-GEOD.zip file. it comes from ArrayExpress not from NCBI GEO

#Read CEL files  as Affybatch object
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")
getwd()
raw_data #gives the basic info about the dataset. Notwdown the annotation (hgu133plus2) 
#from the output. you will need this for selection and installation correctly in next step
# the annotation package (hgu133plus2.db ) to map probe IDs to genes


#### QC before Pre-processing ####

#QC identifies outliers, hybridization problems and/or technical biases
#arrayQualityMetrics generates automated QC reports for microarray data
#It applies multiple complementary methods for detecting technical issues:
  #Density pplots and Box plots: check distribution of intensities
  # MA plots: visualize systemic biasis between arrays
  #Heat maps and distance matrices: identify clustering/outliers
  #PCA: detect unusual variations/ detecting outliers or batch effects

#Output is interactive html report index.html summarizing QC results

arrayQualityMetrics(expressionset = raw_data ,
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

normalized_data <- rma(raw_data)

#QC of the normalized data
arrayQualityMetrics(expressionset = normalized_data ,
                    outdir = "Raw_Data/QC_normalized_data" ,
                    force = TRUE)

#extract normalized expression values into a data frame
processed_Data <- as.data.frame(exprs(normalized_data))

class(normalized_data)

dim(processed_Data) #Dimensions: no. of probes * samples

####Filter Low Variance Transcripts ("Soft" intensity based filtering) ####

#filtering removes low or uninformative expression signals probes
#Reasons: reduce noise and improves statistical power in differential expression
# or ML

#Calculate median intensity per probe across samples
Row_median <- rowMedians(as.matrix(processed_Data))

#Visualize distribution of probe median intensity 
hist(Row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

#Set a threshold value to remove low variance probes (dataset-specific, set accordingly)
threshold <- 3.5

abline(v = threshold , col = "Black" , lwd = 2) #abline a line created at the threshold value with color specified and width

# Select probes value above threshold

indx <- Row_median > threshold
filtered_data <- processed_Data [indx, ]

#rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

#overwrite processed data with filtered data
processed_Data <- filtered_data

#### Phenotype Data Preparation ####
#phenotype data contains sample-level metadata like conditions,
# tissue type or disease status
# required for defining exp group for stats analysis

class(phenotype_data$source_name_ch1)

#Change exp groups as normal vs disease

groups <- factor(phenotype_data$source_name_ch1 , 
                 levels = c("gastric mucosa" , "gastric adenocarcinoma") , 
                 label = c("normal" , "cancer"))
class(groups)
levels(groups)

save.image( file = "3b.RData")
