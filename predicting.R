# In this code, we are going to to practice predictions based on gene expressions. This implies some data preprocessing then building
# a classification model to predict a disease. Here we are continuing to work on a breast cancer dataset downloaded from Firehose
# (https://gdac.broadinstitute.org/). These are the next generation sequencing data that are provided already normalized.


# Preparing the environment
memory.limit(size = 3500)
library(randomForest)
library(class)
source("bss_wss.R")


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


# Feature selection with BSS/WSS
# We will select top genes with bss/wss method
dim(mrnaData)
dim(mrnaClass)
dim(geneNames)
bss <- bssWssFast(X = mrnaData, givenClassArr = t(mrnaClassNum), numClass = 2)


# Classification on selected features on training set
# We extract the top 100 gene expressions and place them in 'mrnaDataReduced'
mrnaDataReduced <- mrnaData[, bss$ix[1:100]]
dim(mrnaDataReduced)

trainSet <- mrnaDataReduced
testSet <- mrnaDataReduced

trainClasses <- unlist(mrnaClassNum[1,], use.names = FALSE)
testClasses <- unlist(mrnaClassNum[1,], use.names = FALSE)

# KNN on selected features on training set
knn.predict <- knn(trainSet, testSet, trainClasses, testClasses, k = 1)
knn.predict <- as.vector(knn.predict)

# Build the confusion matrix
table(knn.predict, testClasses)
tab <- table(knn.predict, t(testClasses))

# Calculate accuracy
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100 - (error * 100 / length(testClasses)))
print(paste("accuracy =", as.character(accuracy), "%"), quote = FALSE)

# Random Forest on selected features on training set
trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,])))
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,])))
colnames(trainSetClass)[101] <- "class"

trainSetClass$class <- as.factor(trainSetClass$class)
class(trainSetClass$class)

rf <- randomForest(class ~ ., trainSetClass, ntree = 100, importance = TRUE)

colnames(testSetClass)[101] <- "class"
testSetClass$class <- as.factor(testSetClass$class)

rf.predict <- predict(rf, testSetClass)
rf.predict <- as.factor(rf.predict)

table(rf.predict, testClasses)
tab <- table(rf.predic, t(testClasses))

error <- sum(tab) - sum(diag(tab))
accuracy <- round(100 - (error * 100 / length(testClasses)))
print(paste("accuracy =", as.character(accuracy), "%"), quote = FALSE)


# Classification on selected features on independent test set
# In this experiment, we are first going to separate the dataset into two sets: (i) a training set composed of 70% of the data
# samples; and (ii) a test set composed of 30% of the data samples
nbRows <- nrow(mrnaDataReduced)
set.seed(33)

trainRows <- sample(1:nbRows, 0.7 * nbRows)
trainSet <- mrnaDataReduced[trainRows,]
testSet <- mrnaDataReduced[-trainRows,]

dim(trainSet)
dim(testSet)

# KNN on selected features on independent test set
trainClasses <- unlist(mrnaClassNum[1,trainRows], use.names = FALSE)
testClasses <- unlist(mrnaClassNum[1,-trainRows], use.names = FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses, k = 1)
knn.predic = as.vector(knn.predic)

table(knn.predic, testClasses)
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100 - (error * 100 / length(testClasses)))
print(paste("accuracy = ", as.character(accuracy), "%"), quote = FALSE)

# Random Forest on selected features on independent test set
trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1, trainRows])))
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1, -trainRows])))
colnames(trainSetClass)[101] <- "class"
trainSetClass$class <- as.factor(trainSetClass$class)
class(trainSetClass$class)

rf <- randomForest(class ~., trainSetClass, ntree = 100, importance = TRUE)
colnames(testSetClass)[101] <- "class"
testSetClass$class <- as.factor(testSetClass$class)
rf.predic <- predict(rf, testSetClass)
rf.predic <- as.vector(rf.predic)

table(rf.predic, testClasses)
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100 - (error * 100 / length(testClasses)))
print(paste("accuracy =", as.character(accuracy), "%"), quote = FALSE)
