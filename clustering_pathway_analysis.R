# In this code, we are going to practice clustering and pathway analysis. Here we are continuing to work on a breast cancer dataset
# downloaded from Firehose. We are going to compare three clustering methods (kmeans, pam, and dbscan). We will also perform gene
# pathway enrichment on the set of genes previously identified for this dataset


# Preparing the environment
memory.limit(size = 3500)
library(cluster)
library(ReactomePA)
source("bss_wss.R")

# Loading the data
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
# mrnaClass and mrnaClassNum are created and contain the diagnostic class - 0 for normal and 1 for tumor
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
sampleType <- as.data.frame(samp)
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)

dim(mrnaNorm)
dim(mrnaClass)
table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum)

# We also create a dataframe with only the genes names, called genesNames, which are located in the first column of mrnaNorm
genesNames <- mrnaNorm[1]

# We transpose the mrnaNorm because we want to select genes
mrnaData <- t(mrnaNorm[, -1])

# Because we are working with large datasets, we free space from memory by removing objects we will not be using anymore
rm(samp)
rm(sampClass)
rm(mrnaNorm)
gc()


# Feature selection with BSS/WSS
# We select the top 100 genes and create a subset of mrnaData with these genes, called mrnaDataReduced. This dataset has 1212
# patients and 100 genes
bss <- bssWssFast(X = mrnaData, givenClassArr = t(mrnaClassNum), numClass = 2)
mrnaDataReduced <- mrnaData[, bss$ix[1:100]]
dim(mrnaDataReduced)

trainClasses <- unlist(mrnaClassNum[1, ], use.names = FALSE)


# Clustering with KMeans
set.seed(20)
kmeans.clusters <- kmeans(mrnaDataReduced, 3, nstart = 20)
table(kmeans.clusters$cluster)
table(kmeans.clusters$cluster, trainClasses)

# We can plot the results and see the intersection between the two clusters. We also see that the clusters are ellipsoid or
# circular, which is probably why the results are not perfect
clusplot(mrnaDataReduced, kmeans.clusters$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)


# Clustering with Pam
# Pam is a improvement over KMeans
pam.clusters <- pam(mrnaDataReduced, 3)
table(pam.clusters$clustering)
table(pam.clusters$clustering, trainClasses)
clusplot(mrnaDataReduced, pam.clusters$clustering, color = TRUE, shade = TRUE, labels = 2, lines = 0)

# In conclusion, we can say that Pam performs better at grouping 0's and 1's in different clusters. However, the ellipsoid shape of
# the clusters limits the the ability to separate the two classes. Maybe other methods would be good to try - such as density based
# clustering or agglomerative clustering


# Pathway analysis
# The genes are provided in TCGA by the hugo_code|entrez_code. Therefore we want to separate the two and split the name into hugo
# name and entrez name. We then pass the list of genes Entrez codes into enrichPathway from ReactomePA package. ReactomePA searches
# the reactome pathway database for pathways having these genes somewhere along their paths
genes <- genesNames[bss$ix[1:100], 1]
genes <- as.character(genes)
hugoNames <- lapply(genes, function(t) substr(t, 1, regexpr("\\|", t) - 1))
entrezNames <- lapply(genes, function(t) substr(t, regexpr("\\|", t) + 1, nchar(t)))
paths <- enrichPathway(unlist(entrezNames), pvalueCutoff = 1)
head(summary(paths))