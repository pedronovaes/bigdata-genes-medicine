# Loading the data

# This dataset has genes as rows and patients as columns. Patient ID's refer to TCGA IDs such as TCGA.3C.AAAU.01A.11R.A41B.07. THe top
# two lines of the file are header information, and we direct 'read.table' to skip these two lines since we only want the gene
# expressions to be placed in a dataframe we call 'mrnaNorm'. The first column of 'mrnaNorm' contains the list of genes measured. Since
# this file is a large dataset, I didn't upload it to Github
mrnaNorm <- read.table(
    "data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
    header = FALSE,
    fill = TRUE,
    skip = 2
)

# We also want the patient IDs in the first row, which is why we read only this first line and place it in a second dataframe called
# 'mrnaIDs'.
mrnaIDs <- read.table(
    "data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
    header = FALSE,
    fill = TRUE,
    nrows = 1
)

mrnaIDs <- mrnaIDs[, -1][, -1]


# Data preprocessing

# Extract first 5 patients IDs
mrnaIDs5 <- mrnaIDs[, c(1, 2, 3, 4, 5)]

# First 5 rows and columns
mrnaNorm5x5 <- mrnaNorm[1:5, 1:5]

head(mrnaIDs5)
head(mrnaNorm5x5)
summary(mrnaNorm5x5)

# The following instruction analyzes the patient IDs because the IDs can tell us whether a particular column is from a normal sample
# or a tumor sample. TCGA has data both from tumor samples and from normal samples taken from the same patient. This information can
# be found in the label of the column. For example in id: TCGA.3C.AAAU.01A.11R.A41B.07 the type of sample is indicated in the 4th
# group: 01A. Tumor types range from 01-09, normal types from 10-19, and control samples from 20-29
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
sampleType <- as.data.frame(samp)
tab <- table(unlist(sampleType))
tab

# We are going to associate a class of '1' for the tumor samples and of '0' for the normal samples
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum)


# Feature Selection with BSS/WSS

# Returns a list with two elements: a list of bss/wss ratios in decreasing order bss[1] and a list of features bss[2]
bss <- bssWssFast(X = t(mrnaNorm[, -1]), givenClassArr = t(mrnaClassNum), numClass = 2)

# Show list of 50 ranked gene indexes
bss$ix[1:50]

# We can then list the genes by their name as provided in the 'mrnaNorm' dataset. The names have two codes separeted by
# character '|': HUGO code Entrez code
# For example, FHL1|2273 is gene with HUGO code FHL1 and Entrez code 2273
genes <- mrnaNorm[bss$ix[1:50], 1]
genes
