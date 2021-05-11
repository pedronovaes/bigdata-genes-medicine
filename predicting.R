# In this code, we are going to to practice predictions based on gene expressions. This implies some data preprocessing then building
# a classification model to predict a disease. Here we are continuing to work on a breast cancer dataset downloaded from Firehose
# (https://gdac.broadinstitute.org/). These are the next generation sequencing data that are provided already normalized.

# Preparing the environment
memory.limit(size = 3500)
library(randomForest)
library(class)

# Loading the data
# We create a dataframe 'mrnaNorm' with the gene expression values and the first column being the gene names. The second dataframe
# 'mrnaIDs' contains the IDs of the patients.
mrnaNorm <- read.table(
    "data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
    header = FALSE,
    fill = TRUE,
    skip = 2
)
mrnaIDs <- read.table(
    "data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
    header = FALSE,
    fill = TRUE,
    nrows = 1
)
mrnaIDs <- mrnaIDs[, -1][, -1]

# Data preprocessing
# 'mrnaClass' and 'mrnaClassNum' are created and contain the diagnostic class - 0 for normal and 1 for tumor
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
sampleType <- as.data.frame(samp)
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)

dim(mrnaNorm)
dim(mrnaIDs)
dim(mrnaClass)
table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum)

# We also created a dataframe with only the gene names, called 'geneNames', which are located in the first column of 'mrnaNorm'
geneNames <- mrnaNorm[1]
dim(geneNames)

mrnaData <- t(mrnaNorm[, -1])
rm(samp)
rm(sampClass)
rm(mrnaNorm)
gc()