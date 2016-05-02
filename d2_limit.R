#install.packages("signal")
library(signal)
library(pracma)
testDataNo = 2
if (testDataNo == 1) {
  freq <- 1
  count <- 10000
  timeSpan <- 200
  x <- seq(-timeSpan/2,timeSpan/2, length.out= count)
  s <- sin(2*pi*freq*x)
  dat <- data.frame(x, s)
  write.table(dat, 'testdata1.csv', row.names = FALSE, col.names = FALSE)
} else if (testDataNo == 2) {
  dat <- read.table('testdata2.csv', sep='', header = FALSE)
  x <- dat[,1]
  s <- dat[,2]
  count <- length(x)
  timeSpan <- x[length(x)] - x[1]
}
sampling <- timeSpan/count
cohlen <- 5
ln2 <- log(2)
#plot(x, s, type='l')

######################
#p <- Mod(fft(s))
#freqSampling <- 1/sampling
#f <- freqSampling/2*seq(-1, 1, length.out= count)
#print(length(p))
#p1 <- p[1:(length(p)/2)]
#p2 <- p[(length(p)/2+1):(length(p))]
#p <- c(p2, p1)
#print(length(p))
#print(length(f))
##p <- p[1:length(f)]
##plot(f,p, 'l')
#g <- cohlen * sqrt(pi/ln2) * exp(-pi^2*cohlen^2/ln2*f^2)
##g <- dnorm(f, 0, 0.1)
#c <- conv(g,p)
#print(length(c))
#f1 <- seq(f[1]*2, f[length(f)]*2, length.out= length(c))
##plot(f1, c, 'l')
######################



cff <- acf(s, lag.max=count, plot=TRUE)$acf
#plot(cff$lag*sampling, cff1acf)

term1 <- function(tau, freq) {
  exp(-ln2*tau^2/cohlen^2)*(cos(2*pi*tau*freq))*cff[round(tau/sampling)+1]
}

term2 <- function(tau) {
  exp(-ln2*tau^2/cohlen^2)*cff[round(tau/sampling)+1]
}

norm <- function(tau, freq) {
  exp(-ln2*tau^2/cohlen^2)*(cos(2*pi*tau*freq) + 1)
}

freqs <- seq(0.5, 1.5, 0.01)
d2 <- c()
taus <- seq(0, (length(cff) - 1)*sampling, sampling)
sampledTerm2 <- c()
for (tau in taus) {
  sampledTerm2 <- c(sampledTerm2, term2(tau))
}
intTerm2 <- 2 * trapz(taus, sampledTerm2) #integrate(term2, 0, 3*cohlen, subdivisions=1000)$value
print(intTerm2)
for (freq in freqs) {
  
  sampledTerm1 <- c()
  for (tau in taus) {
    sampledTerm1 <- c(sampledTerm1, term1(tau, freq))
  }
  #plot(taus, sampledTerm1, 'l')
  intTerm1 <- 2*trapz(taus, sampledTerm1) #integrate(term1, 0, 3*cohlen, subdivisions=1000)$value
  intNorm <- integrate(function(tau) norm(tau, freq), lower=-Inf, upper=Inf, subdivisions=1000)$value
  d2 <- c(d2, 1 - (intTerm1 + intTerm2) / intNorm)
  #d2 <- c(d2, intTerm1)
}


d2_expected <- read.table(paste('d2_expected_',toString(testDataNo),'_',toString(cohlen),'.csv', sep=''), header = FALSE)
plot(d2_expected, type='l', col = 'red')
lines(freqs, d2, col = 'green')
