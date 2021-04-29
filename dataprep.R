# Libs
library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(ggplot2)


# Loading the data
dataset <- read.table("data/heart-h.txt", header = TRUE, sep = ",", quote = "")

nrow(dataset)
ncol(dataset)
dim(dataset)
head(dataset)
summary(dataset)


# Dealing with missing values

# md.pattern function displays a table of the number of missing values in each column, where '0' represents a missing value, and '1'
# represents a present value.
md.pattern(dataset)

# aggr function ranks the variables by decreasing number of missing values. The variables at the top of the list are those having
# the most missing values. Here are the top, we have variable 'ca' with almost 100% of missing values
mice_plot <- aggr(
    dataset,
    col = c("green", "red"),
    numbers = TRUE,
    sortVars = TRUE,
    labels = names(dataset),
    cex.axis = .7,
    gap = 3,
    ylab = c("Missing data", "Pattern")
)

# To deal with missing values, one main method is to impute these values, for example with the mean of the column where the missing
# value is located.
dataset$chol_imputed <- with(dataset, impute(chol, mean))
summary(dataset)


# Normalization

# A function 'Normalize' exists in RWeka package to normalize using ZScore normalization
dataset_n <- Normalize(~., data = dataset)
summary(dataset_n)


# Discretization
# The main methods for discretization are equal-width and equal-depth

# Discretization can be accomplished by the 'cut2' function in several manners:
# - with 'g=3', we create 3 bins with approximately the same number of elements in each (equal-depth)
# - with 'm=100', we create bins of depth 100 (equal-depth)
dataset$chol_bin <- as.numeric(cut2(dataset$chol_imputed, g = 3))
summary(dataset)
