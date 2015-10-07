library(MASS)
library(foreach)
library(doMC)
registerDoMC(7)

#spots <- read.table("yearssn.dat", header=FALSE)
spots <- read.table("dayssnv0.dat", header=FALSE, sep="\t")
spots <- spots[spots[,2] != 999,]
storms <- read.table("aa_gt50.csv", header=FALSE, sep="\t")
#storms <- read.table("kakioka_gt40.csv", header=FALSE, sep="\t")
epochs <- scan("minima_epochs.csv")
epochs <- epochs[epochs > 1878] # for aa-index
#epochs <- epochs[epochs > 1923] # for Kakioka
thresholds <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600) # for aa-index
#thresholds <- c(50, 100, 150, 200, 250, 300, 350) # for Kakioka

filterStorms <- function(storms, threshold) {
    return (storms[storms[,2] >= threshold,])
}

findNearest <- function(x, dat) {
    nearest <- unlist(dat[1,])
    nearestDiff <- abs(nearest[1] - x)
    for (i in 1:nrow(dat)) {
	val <- unlist(dat[i,])
        diff <- val[1] - x
        if (abs(diff) < nearestDiff) {
            nearest <- val
            nearestDiff <- abs(diff)
        }
        if (diff > 0) {
            break
        }
    }
    return (nearest)
    #unlist(dat[which.min(lapply(abs(unlist(dat[,1]) - x),min)),])
}

normalize <- function(dat) {
    minY <- min(unlist(dat[,2]))
    maxY <- max(unlist(dat[,2]))
    scaleY <- 1 / (maxY - minY)
    norm <- function(x) {
        (x - minY) * scaleY
    } 
    return (cbind(dat[,1], lapply(dat[,2], norm)))
}

int <- function(count, dat, valueOrCount = TRUE) {
    xMin <- unlist(dat[1,1])
    xMax <- unlist(dat[nrow(dat),1])
    step <- (xMax - xMin) / count
    res <- matrix(nrow=0, ncol = 2)
    for (x in seq(xMin, xMax - step, by = step)) {
        bin <- dat[dat[,1] >= x & dat[,1] < x + step, 2]
        if (valueOrCount) {
            res <- rbind(res, c(x + step / 2, sum(unlist(bin))))
        } else {
            res <- rbind(res, c(x + step / 2, length(bin)))
        }
    }
    return (normalize(res))
}


calcPhaseShift <- function(spots, storms1, threshold = 0, bootstrap = FALSE) {

    storms <- filterStorms(storms1, threshold)
    spotDist <- matrix(nrow=0, ncol = 2)
    stormDist <- matrix(nrow=0, ncol = 2)

    cycleIndices <- 1:(length(epochs)-1)
    if (bootstrap) {
        cycleIndices <- sample(cycleIndices, length(cycleIndices),replace=TRUE)
    }

    for (i in cycleIndices) {
        cycleStart <- epochs[i]
        cycleEnd <- epochs[i+1]
        norm <- function(x) (x - cycleStart) / (cycleEnd - cycleStart)

        spotsInCycle <- spots[spots[,1] >= cycleStart & spots[,1] < cycleEnd,]
        spotsInCycle <- cbind(lapply(spotsInCycle[,1], norm), spotsInCycle[,2])
        spotDist <- rbind(spotDist, normalize(spotsInCycle))

        stormsInCycle <- storms[storms[,1] >= cycleStart & storms[,1] < cycleEnd,]
        stormsInCycle <- cbind(lapply(stormsInCycle[,1], norm), stormsInCycle[,2])
        stormDist <- rbind(stormDist, stormsInCycle)
    }

    spotDist <- spotDist[order(unlist(spotDist[,1])),]
    stormDist <- stormDist[order(unlist(stormDist[,1])),]

    count1 = 50
    count2 <- min(count1, nrow(stormDist)/2)
    spotInt <- int(count1, spotDist)
    stormInt <- int(count2, stormDist, FALSE)

    if (!bootstrap) {
        if (threshold < 100) {
            meanSpotInt <- smooth.spline(spotInt)
            meanStormInt <- smooth.spline(stormInt)
            postscript("meandist.ps")
            #plot(spotInt, xlab="Phase", ylab="Mean distribution", type ="lines", lty=1)
            spotMin = min(meanSpotInt$y)
            spotMax = max(meanSpotInt$y)
            plot(meanSpotInt$x, (meanSpotInt$y - spotMin) / (spotMax - spotMin), xlab="Phase", ylab="Normalized value", type ="lines", lty=1)
            #lines(stormInt, lty=2)
            stormMin = min(meanStormInt$y)
            stormMax = max(meanStormInt$y)
            lines(meanStormInt$x, (meanStormInt$y - stormMin) / (stormMax - stormMin), lty=2)
            legend("topright", c("Sunspot number","AA-index"), lty=c(1,2))
        }
        #plot(stormDist)
        #plot(stormInt)
    }
    phaseshift <- matrix(nrow=0, ncol = 2)
    for (phase in seq(-0.5, 0.5 - 1/count1, by = 1/count1)) {
        dist <- 0
        for (i in 1:nrow(stormInt)) {
            y <- unlist(stormInt[i,2])
            if (y > 0) {
                x <- unlist(stormInt[i,1]) + phase
                if (x >= 1) x <- x - 1
                if (x < 0) x <- x + 1
                nearest <- findNearest(x, spotInt)
                dist <- dist + (nearest[2] - y) ^2
            }
        }
        phaseshift <- rbind(phaseshift, c(phase, dist))
        #print(sprintf("%f %f", phase, dist))
    }
    #disp <- approx(unlist(phaseshift[,1]), unlist(phaseshift[,2]), n=2)
    #print(disp)
    s <- smooth.spline(phaseshift)
    #plot(phaseshift)
    p <- predict(s, seq(-0.5, 0.5, by = 1/100))
    #lines(p)
    #return (unlist(phaseshift[which.min(lapply(phaseshift[,2],min)),])[1])
    shift <- p$x[which.min(lapply(p$y,min))]
    return (shift)
}

#test <- seq(-0.5, 0.5, by = 1/20)
#test <- lapply(test + 0.5, function (x) if (x >= 0.5) {x - 1} else {x})
#print(test)

shifts <- c()
shiftErrors  <- c()
bsShiftMeans <- c()
for (threshold in thresholds) {
    print(threshold)
    shift <- calcPhaseShift(spots, storms, threshold)
    bsCount <- 0#100
    #bsShifts <- c()
    bsShifts <- foreach (i=1:bsCount, .combine=rbind) %dopar% {
        print(i)
        bsShift <- calcPhaseShift(spots, storms, threshold, TRUE)
        #bsShifts <- rbind(bsShifts, bsShift)
        return(bsShift)
    }
    minKsRes <- -1
    minDist <- c()
    minDelta <- 0
    for (delta in seq(0, 0.5, by = 0.01)) {
        dist <- bsShifts + delta
        dist <- sapply(dist, function (x) if (x >= 0.5) {x - 1} else {x})
        #print(dist)
        fit=fitdistr(dist, "normal")
        ksRes=ks.test(dist,rnorm(length(dist), fit$estimate[1], fit$estimate[2]))$statistic
        #var <- var(dist)
        if (minKsRes < 0 || ksRes < minKsRes) {
            minKsRes <- ksRes
            minDist <- dist
            minDelta <- delta
        }
    }
    #hist(minDist)
    #print(sprintf("Threshold shift bsShift sd: %f %f %f %f", threshold, shift, mean(bsShifts), sd(bsShifts)))
    #bsShifts1 <- apply(cbind(bsShifts + 0.5, 0.5 - bsShifts), 1, min) #c(bsShifts[bsShifts >= 0],bsShifts[bsShifts < 0] + 1)
    shifts <- rbind(shifts, shift)
    shiftErrors <- rbind(shiftErrors, 2*sd(minDist))
    #meanShift <- mean(minDist) - minDelta
    #if (meanShift < -0.5) {
    #    meanShift <- meanShift + 1;
    #}
    #bsShiftMeans <- rbind(bsShiftMeans, meanShift)
}
postscript("mygraph.ps")
plot(thresholds, shifts,
    ylim=range(c(shifts-shiftErrors, shifts+shiftErrors)), yaxt="n",
    pch=19, xlab="Storm threshold", ylab="Phase shift"
)
phaseTics <- c(-0.5,-0.25,0,0.25,0.5)
axis(2, at=phaseTics,labels=phaseTics, las=2)
#points(thresholds, bsShiftMeans)
arrows(thresholds, shifts-shiftErrors, thresholds, shifts+shiftErrors, length=0.05, angle=90, code=3)
abline(h=phaseTics, lty=2)

