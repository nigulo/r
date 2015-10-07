library(MASS)
library(foreach)
library(doMC)
registerDoMC(7)

storms <- read.table("aa_gt50.csv", header=FALSE, sep="\t")
#storms <- read.table("kakioka_gt40.csv", header=FALSE, sep="\t")
#storms <- read.table("aaindexmean.csv", header=FALSE, sep="\t")
epochs <- scan("epochs.csv")
#epochs <- epochs[epochs > 1924] # Kakioka

#thresholds = c(25, 50, 100, 200, 300, 400) # aa-index mean
thresholds = c(0, 100, 200, 300, 400, 500, 600) # aa-index
#thresholds = c(50, 100, 150, 200, 250, 300, 350) # Kakioka

ascStartInds <- seq(1, length(epochs), 4)
maxStartInds <- seq(2, length(epochs), 4)
decStartInds <- seq(3, length(epochs), 4)
minStartInds <- seq(4, length(epochs), 4)

ascStarts <- epochs[ascStartInds]
maxStarts <- epochs[maxStartInds]
decStarts <- epochs[decStartInds]
minStarts <- epochs[minStartInds]
j = 1
#ascs <- matrix(nrow = 0, ncol = 2)
#maxs <- matrix(nrow = 0, ncol = 2)
#decs <- matrix(nrow = 0, ncol = 2)
#mins <- matrix(nrow = 0, ncol = 2)
ascLen = 0
maxLen = 0
decLen = 0
minLen = 0
numCycles <- length(ascStartInds)
for (i in 1:numCycles) {
    ascLen = ascLen + maxStarts[i] - ascStarts[i]
    maxLen = maxLen + decStarts[i] - maxStarts[i]
    decLen = decLen + minStarts[i] - decStarts[i]
    if (i < length(ascStartInds)) {
        minLen = minLen + ascStarts[i+1] - minStarts[i]
    }
}

print(sprintf("num cycles: %d", numCycles))
ascLen/numCycles
maxLen/numCycles
decLen/numCycles
minLen/numCycles

setClass("stormsInPhase", slots = c(ascs="data.frame", maxs="data.frame", decs="data.frame", mins="data.frame"))

stormsInPhase <- function(storms) {
    ascs <- matrix(nrow = 0, ncol = 2)
    maxs <- matrix(nrow = 0, ncol = 2)
    decs <- matrix(nrow = 0, ncol = 2)
    mins <- matrix(nrow = 0, ncol = 2)
    for (i in 1:numCycles) {
        ascs <- rbind(ascs, storms[storms[,1] >= ascStarts[i] & storms[,1] < maxStarts[i],])
        maxs <- rbind(maxs, storms[storms[,1] >= maxStarts[i] & storms[,1] < decStarts[i],])
        decs <- rbind(decs, storms[storms[,1] >= decStarts[i] & storms[,1] < minStarts[i],])
        if (i < length(ascStartInds)) {
            mins <- rbind(mins, storms[storms[,1] >= minStarts[i] & storms[,1] < ascStarts[i+1],])
        }
    }
    retVal <- new("stormsInPhase", ascs = ascs, maxs = maxs, decs = decs, mins = mins)
    return (retVal)
}

bsStorms <- function (storms) {
    bsIndices <- sample(1:nrow(storms),nrow(storms),replace=TRUE)
    #print(bsIndices[1:10])
    return (storms[bsIndices,])
}

filterStorms <- function(storms, threshold) {
    return (storms[storms[,2] >= threshold,])
}

countStorms <- function(storms) {
    return (nrow(storms))
}

calcFreqs <- function(storms, len, print=TRUE) {
    stormCounts = sapply(thresholds, function(threshold) countStorms(filterStorms(storms, threshold)))
    #if (print) {print(stormCounts)}
    stormFreqs = stormCounts / len
    return (stormFreqs)
}

printFreqs <- function(ascs, maxs, decs, mins) {
    ascFreqs = calcFreqs(ascs, ascLen)
    maxFreqs = calcFreqs(maxs, maxLen)
    decFreqs = calcFreqs(decs, decLen)
    minFreqs = calcFreqs(mins, minLen)
    print (rbind(ascFreqs, maxFreqs, decFreqs, minFreqs))
}

stormsPerPhase <- stormsInPhase(storms)
printFreqs(stormsPerPhase@ascs, stormsPerPhase@maxs, stormsPerPhase@decs, stormsPerPhase@mins)

setClass("stormsFreqsInPhase", slots = c(ascs="numeric", maxs="numeric", decs="numeric", mins="numeric"))

bsCount = 1000
allAscFreqs <- c()
allMaxFreqs <- c()
allDecFreqs <- c()
allMinFreqs <- c()
allStormFreqsInPhase <- foreach (i=1:bsCount) %dopar% {
#for (i in 1:bsCount) {
    #if (i %% 10 == 0) {print(i)}
    bs <- bsStorms(storms)
    stormsPerPhase <- stormsInPhase(bs)
    ascFreqs = calcFreqs(stormsPerPhase@ascs, ascLen)
    #allAscFreqs <- rbind(allAscFreqs, ascFreqs)
    maxFreqs = calcFreqs(stormsPerPhase@maxs, maxLen)
    #allMaxFreqs <- rbind(allMaxFreqs, maxFreqs)
    decFreqs = calcFreqs(stormsPerPhase@decs, decLen)
    #allDecFreqs <- rbind(allDecFreqs, decFreqs)
    minFreqs = calcFreqs(stormsPerPhase@mins, minLen)
    #allMinFreqs <- rbind(allMinFreqs, minFreqs)
    retVal <- new("stormsFreqsInPhase", ascs = ascFreqs, maxs = maxFreqs, decs = decFreqs, mins = minFreqs)
    return (retVal)
}

for (i in 1:bsCount) {
    allAscFreqs <- rbind(allAscFreqs, allStormFreqsInPhase[[i]]@ascs)
    allMaxFreqs <- rbind(allMaxFreqs, allStormFreqsInPhase[[i]]@maxs)
    allDecFreqs <- rbind(allDecFreqs, allStormFreqsInPhase[[i]]@decs)
    allMinFreqs <- rbind(allMinFreqs, allStormFreqsInPhase[[i]]@mins)
}

rbind(apply(allAscFreqs, 2, sd), apply(allMaxFreqs, 2, sd), apply(allDecFreqs, 2, sd), apply(allMinFreqs, 2, sd))

