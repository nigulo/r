library(class)
#library(stringi)

attributes <- c("age", "workclass", "fnlwgt", "education", "education-num", "marital-status", "occupation", "relationship", "race", "sex", "capital-gain", "capital-loss", "hours-per-week", "native-country")
attrTypes <- c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1) # 0 - continuous, 1 - discrete
data <- read.csv('data/adult.data', header=FALSE, strip.white = TRUE);

#data <- data[data[,] != '?',]
attrVals <- list()
for (col in 1:ncol(data)) {
    attrVals <- append(attrVals, list(as.list(unique(data[,col]))))
}

minLeafSize <- 1000
minEntropy <- 0.001
classAttr <- 15
classLabel1 <- unlist(attrVals[[classAttr]][1])
classLabel2 <- unlist(attrVals[[classAttr]][2])
sensitiveAttr <- 10
sensitiveAttrVal1 <- unlist(attrVals[[sensitiveAttr]][1])

nodeEntropy <- function(x, level) {
    entropy <- 0
    len <- nrow(x)
    for (classAttrVal in attrVals[[classAttr]]) {
        p <- nrow(x[x[,classAttr]==classAttrVal,]) / len
        if (p > 0) { # assuming 0 log 0 = 0
            entropy <- entropy - p * log2(p)
        }
    }
    #print(paste(level, sprintf("Entropy %f", entropy)))
    return (entropy)
}

splitEntropy <- function(subsets, level) {
    entropy <- 0
    len <- 0
    for (subset in subsets) {
        subsetLen <- nrow(subset)
        if (subsetLen == 0) {
            next;
        }
        subEntropy <- 0
        for (classAttrVal in attrVals[[classAttr]]) {
            p <- nrow(subset[subset[,classAttr]==classAttrVal,]) / subsetLen
            if (p > 0) { # assuming 0 log 0 = 0
                subEntropy <- subEntropy + p * log2(p)
            }
        }
        entropy <- entropy - subEntropy * subsetLen
        len <- len + subsetLen
    }
    entropy <- entropy / len
    #print(paste(level, sprintf("Entropy %f", entropy)))
    return (entropy)
}

splitFunc <- function(x, node) {
    attrIndex <- node$attr
    if (attrTypes[attrIndex] == 1) { # attribute is discrete
        #print("Discrete");
        return (match(x[attrIndex], attrVals[[attrIndex]]))
    } else { # attribute is continuous
        splitPoint <- node$splitParams;
        #print(paste("Continuous with split point", format(splitPoint)));
        if (x[attrIndex] <= splitPoint) {
            return (1)
        } else {
            return (2)
        }
    }
}

split <- function(x, level) {
    minEnt <- NULL
    bestSplit <- NULL
    bestSplitAttr <- NULL
    for (attrIndex in 1:(ncol(x) - 1)) {
        if (attrIndex == classAttr || attrIndex == sensitiveAttr) {
            next # skip class attribute and sensitive attribute
        }
        if (length(attrVals[[attrIndex]]) == 1) {
            next # skip attributes possibly with only one unique value
        }
        subsets <- list()
        splitParams <- NULL
        if (attrTypes[attrIndex] == 1) { # attribute is discrete
            #print(paste(level, sprintf("Discrete %d", attrIndex)))
            for (attrVal in attrVals[[attrIndex]]) {
                subset <- x[x[,attrIndex] == attrVal,]
                subsets <- append(subsets, list(subset))
            }
        } else { # attribute is continuous
            #print(paste(level, sprintf("Continuous %d", attrIndex)))
            splitPoint <- mean(x[,attrIndex])
            splitParams <- splitPoint
            subset1 <- x[x[,attrIndex] <= splitPoint,]
            subset2 <- x[x[,attrIndex] > splitPoint,]
            subsets <- append(subsets, list(subset1))
            subsets <- append(subsets, list(subset2))
        }
        e <- splitEntropy(subsets, level)
        if (is.null(minEnt) || e < minEnt) {
            minEnt <- e
            bestSplit <- subsets
            bestSplitAttr <- attrIndex
            bestSplitParams <- splitParams
        }
    }
    #print(paste(level, sprintf("Splitting by attribute %d into %d nodes", bestSplitAttr, length(bestSplit))))
    return (list(attr=bestSplitAttr, splitParams=bestSplitParams, subnodes=bestSplit))
}

createLeaf <- function(x, totalNumSensitive1, totalNumSensitive2, level) {
    maxCount <- NULL
    classLabel <- NULL
    for (classAttrVal in attrVals[[classAttr]]) {
        count <- nrow(x[x[,classAttr] == classAttrVal,])
        if (is.null(maxCount) || count > maxCount) {
            maxCount <- count
            classLabel <- classAttrVal
        }
    }
    numPositive <- nrow(x[x[,classAttr] == classLabel1,])
    numNegative <- nrow(x) - numPositive
    numSensitive1 <- nrow(x[x[,sensitiveAttr] == sensitiveAttrVal1,])
    numSensitive2 <- nrow(x) - numSensitive1
    deltaAcc <- numNegative - numPositive
    deltaDisc <- numSensitive1 / totalNumSensitive1 - numSensitive2 / totalNumSensitive2
    if (deltaAcc > 0) {
        deltaAcc <- -deltaAcc
        deltaDisc <- -deltaDisc
    }
    #print(paste(level, "Creating leaf with class label", format(classLabel), ", dAcc=", deltaAcc, ", dDisc=", deltaDisc))
    return (list(class=classLabel, dAcc=deltaAcc, dDisc=deltaDisc))
}

generateTree <- function(x, totalNumSensitive1, totalNumSensitive2, level) {
    if (nrow(x) <= minLeafSize || nodeEntropy(x, level) < minEntropy) {
        return (createLeaf(x, totalNumSensitive1, totalNumSensitive2, level))
    } else {
        node <- split(x, level)
        nodes <- list()
        i <- 1
        for (subnode in node$subnodes) {
            nodes <- append(nodes, list(generateTree(subnode, totalNumSensitive1, totalNumSensitive2, paste0(level, ".", i))))
            i <- i + 1
        }
        return (list(attr=node$attr, splitParams=node$splitParams, subnodes=nodes))
    }
}

# Relabel one random leaf in tree which has dDisc negative
relabelLeaf <- function(tree) {
    if (is.null(tree$attr)) { # leaf node
        if (tree$dDisc < 0) {
            classLabel <- classLabel1
            if (tree$class == classLabel) {
                classLabel <- classLabel2
            }
            return (list(class=classLabel, dAcc=tree$dAcc, dDisc=-tree$dDisc))
        } else {
            # do nothing as the discrimination would increase
            return (list(class=tree$class, dAcc=0, dDisc=0))
        }
    } else {
        path <- sample(1:length(tree$subnodes), 1)
        #print(sprintf("Using attribute %d and path %d", tree$attr, path))
        newSubNode <- relabelLeaf(tree$subnodes[[path]])
        if (path == 1) {
            return (list(attr=tree$attr, splitParams=tree$splitParams, dAcc=newSubNode$dAcc, dDisc=newSubNode$dDisc, subnodes=c(list(newSubNode), tree$subnodes[(path + 1):length(tree$subnodes)])))
        } else if (path == length(tree$subnodes)) {
            return (list(attr=tree$attr, splitParams=tree$splitParams, dAcc=newSubNode$dAcc, dDisc=newSubNode$dDisc, subnodes=c(tree$subnodes[1:(path-1)], list(newSubNode))))
        } else {
            return (list(attr=tree$attr, splitParams=tree$splitParams, dAcc=newSubNode$dAcc, dDisc=newSubNode$dDisc, subnodes=c(tree$subnodes[1:(path-1)], list(newSubNode), tree$subnodes[(path + 1):length(tree$subnodes)])))
        }
    }
}

# We use an approach similar to Monte-Carlo method
# and just generate random new trees with relabeled leaves
# and discrimination drop more than specified value
relabel <- function(tree, discDrop) {
    dAcc <- 0
    dDisc <- 0
    numLeaves <- 0
    newTree <- tree
    while (dDisc < discDrop) {
        newTree <- relabelLeaf(newTree)
        dAcc <- dAcc + newTree$dAcc
        dDisc <- dDisc + newTree$dDisc
        numLeaves <- numLeaves + 1
    }
    print(paste("Number of leaves relabeled:", numLeaves, ", dAcc:", dAcc))
    return (newTree)
}

predict <- function(x, tree) {
    if (is.null(tree$attr)) { # leaf node
        return (tree$class)
    } else {
        path <- splitFunc(x, tree)
        #print(sprintf("Using attribute %d and path %d", tree$attr, path))
        return (predict(x, tree$subnodes[[path]]))
    }
}

# Creating baseline datasets with the most correlated attributes removed
corrs <- c()
for (i in 1:ncol(data)) {
	if (i == sensitiveAttr || i == classAttr) {
		next;
	}
	corrs <- c(corrs, list(attr=i, corr=abs(correlation(data[,i], data[,sensitiveAttr]))))
}
baselineData <- list()
reducedData <- data
while (length(corrs) > 0) {
	maxCorr <- 0
	maxAttr <- NULL
	for (j in 1:length(corrs)) {
		if (corrs[j]$corr > maxCorr) {
			maxCorr <- corrs[j]$corr
			maxAttr <- corrs[j]$attr
		}
	}
	lastData <- baselineData[length(baselineData)]
	reducedData <- lastData[,1:maxAttr]
	baselineData <- append(baselineData, newData)
	lastData <- reducedData
}
# End creating baseline datasets

chunkSize <- ceiling(nrow(data) / 10)
meanAcc <- 0
meanDisc1 <- 0
meanDisc2 <- 0
relabMeanAccs <- NULL
relabMeanDiscs1Relab <- NULL
relabMeanDiscs2Relab <- NULL
for (i in 0:9) { # 10-fold cross validation
    print(paste("Building tree", i + 1))
    validation <- data[(i*chunkSize + 1):(min((i+1)*chunkSize,nrow(data))),] # 10% of data for validation
    ranges <- c()
    if (i > 0) {
        ranges <- c(ranges, c(1:(i*chunkSize)))
    }
    if (i < 9) {
        ranges <- c(ranges, c(((i+1)*chunkSize + 1):nrow(data)))
    }
    train <- data[ranges,] # 90% of data for training
    numSensitive1 <- nrow(train[train[,sensitiveAttr] == unlist(sensitiveAttrVal1),])
    numSensitive2 <- nrow(train) -  numSensitive1
    tree <- generateTree(train, numSensitive1, numSensitive2, "1")
    acc <- 0
    disc1 <- 0
    disc2 <- 0
    print(paste("Validating tree", i + 1))
    for (j in 1:nrow(validation)) {
        x <- validation[j,]
        class <- predict(x, tree)
        #print(sprintf("Prediction: %s, true label: %s", format(label), format(x[classAttr])))
        if (class == unlist(x[classAttr])) {
            acc <- acc + 1
        }
        if (class == classLabel1) {
            if (unlist(x[sensitiveAttr]) == sensitiveAttrVal1) {
                disc1 <- disc1 + 1
            } else {
                disc2 <- disc2 + 1
            }
        }
    }
	numSensitive1 <- nrow(validation[validation[,sensitiveAttr] == sensitiveAttrVal1,])
	numSensitive2 <- nrow(validation) -  numSensitive1
	disc <- disc1 / numSensitive1 - disc2 / numSensitive2
	print(paste("Accuracy:", acc / nrow(validation)))
	print(paste("Discrimination:", disc))
	dDiscs <- seq(abs(disc) / 10, abs(disc) - abs(disc) / 10, by = abs(disc) / 10)
	#if (is.null(meanDDiscs)) {
	#	meanDDiscs <- dDiscs
	#} else {
	#	meanDDiscs <- meanDDiscs + dDiscs
	#}
	relabAccs <- c()
	relabDiscs1 <- c()
	relabDiscs2 <- c()
	for (dDisc in dDiscs) {
		print(paste("Relabeling tree with dDisc", dDisc))
		treeRelab <- relabel(tree, dDisc)
		accRelab <- 0
		disc1Relab <- 0
		disc2Relab <- 0
		print(paste("Validating relabeled tree", i + 1))
		for (j in 1:nrow(validation)) {
			x <- validation[j,]
			classRelab <- predict(x, treeRelab)
			if (classRelab == unlist(x[classAttr])) {
				accRelab <- accRelab + 1
			}
			if (classRelab == classLabel1) {
				if (unlist(x[sensitiveAttr]) == sensitiveAttrVal1) {
					disc1Relab <- disc1Relab + 1
				} else {
					disc2Relab <- disc2Relab + 1
				}
			}
		}
		relabAccs <- c(relabAccs, accRelab)
		relabDiscs1 <- c(relabDiscs1, disc1Relab)
		relabDiscs2 <- c(relabDiscs2, disc2Relab)
		
	}
		
	#print(paste("Accuracy after relabeling:", accRelab / nrow(validation)))
	#print(paste("Discrimination after relabeling:", disc1Relab / numSensitive1 - disc2Relab / numSensitive2))
    meanAcc <- meanAcc + acc
    meanDisc1 <- meanDisc1 + disc1
    meanDisc2 <- meanDisc2 + disc2
	if (is.null(relabMeanAccs)) {
		relabMeanAccs <- relabAccs
		relabMeanDiscs1 <- relabDiscs1
		relabMeanDiscs2 <- relabDiscs2
	} else {
		relabMeanAccs <- relabMeanAccs + relabAccs
		relabMeanDiscs1 <- relabMeanDiscs1 + relabDiscs1
		relabMeanDiscs2 <- relabMeanDiscs2 + relabDiscs2
	}
}
numSensitive1 <- nrow(data[data[,sensitiveAttr] == unlist(sensitiveAttrVal1),])
numSensitive2 <- nrow(data) -  numSensitive1
print(paste("10-Fold cross validation accuracy:", meanAcc / nrow(data)))
print(paste("10-Fold cross validation discrimination:", meanDisc1 / numSensitive1 - meanDisc2 / numSensitive2))
relabMeanAccs <- relabMeanAccs / nrow(data)
relabMeanDiscs <- relabMeanDiscs1 / numSensitive1 - relabMeanDiscs2 / numSensitive2
print(paste("10-Fold cross validation accuracy after relabeling:", relabMeanAccs))
print(paste("10-Fold cross validation discrimination after relabeling:", relabMeanDiscs))
plot(relabMeanDiscs, relabMeanAccs)