# Here I provide a function that calculates a ranking of features (genes in this case) with the BSS/WSS method
bssWssFast <- function(X, givenClassArr, numClass=2) {
    classVec <- matrix(0, numClass, length(givenClassArr))
    
    for (k in 1:numClass) {
        temp <- rep(0, length(givenClassArr))
        temp[givenClassArr == (k - 1)] <- 1
        classVec[k, ] <- temp
    }
    
    classMeanArr <- rep(0, numClass)
    ratio <- rep(0, ncol(X))
    
    for (j in 1:ncol(X)) {
        overallMean <- sum(X[, j]) / length(X[, j])
        
        for (k in 1:numClass) {
            classMeanArr[k] <- sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
        }
        
        classMeanVec <- classMeanArr[givenClassArr + 1]
        bss <- sum((classMeanVec - overallMean) ^ 2)
        wss <- sum((X[, j] - classMeanVec) ^ 2)
        ratio[j] <- bss / wss
    }
    
    sort(ratio, decreasing = TRUE, index = TRUE)
}