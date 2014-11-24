##########################################
# Collective Risk classes (S3) and methods
##########################################

# distribution methods --------------------------------------------

# mean
mean <- function(d) UseMethod("mean", d)
# limited average severity
LAS <- function(d, ...) UseMethod("LAS", d)
# limited second moment of severity
SECM <- function(d, ...) UseMethod("SECM", d)
# limited variance
limVar <- function(d, ...) UseMethod("limVar", d)
# cummulative density function
cdf <- function(d, ...) UseMethod("cdf", d)
# inverse cummulative density function
icdf <- function(d, ...) UseMethod("icdf", d)
# discretize a continuous distribution
discretize <- function(d, ...) UseMethod("discretize", d)


# default methods -------------------------------------------------

limVar.default <- function(d,x) {
  SECM(d,x) - (LAS(d,x))^2
}

discretize.default <- function(d,h,N){
  X<-rep(0,times=n)
  X[1] <- 1 - LAS(d,h)/h
  for (k in 2:n)
  {
    X[k] <- (2*LAS(d, (h*(k-1)))-LAS(d, (h*k))-LAS(d, (h*(k-2))))/h
  }
  X  
}

# mixed exponential distribution ----------------------------------

# constructor
MixedExponential <- function(weights=numeric(), means=numeric()) {
  x <- list(weights=weights, means=means)
  class(x) <- "MixedExponential"
  x 
}

# methods
mean.MixedExponential <- function (d) { 
  sum(d$weights * d$means) 
}

LAS.MixedExponential <- function (d, x) { 
  sum(d$weights*d$means*(1-exp(-x/d$means))) 
}

cdf.MixedExponential <- function(d, x) {
  1-sum(d$weights*exp(-x/d$means))
}

icdf.MixedExponential <- function(d,p) {
  f <- function(x){p-cdf(d,x)}
  xmax <- 100000
  while (f(xmax)>0) {xmax <- xmax * 10}
  uniroot(f,c(0,xmax))$root
}

SECM.MixedExponential <- function(d,x){
  sum(2 * d$weights * (d$means)^2 * (1 - (1 + x/d$means) * exp(-x/d$means)))
}

# layer distribution ---------------------------------------------

# constructor
LayerDistribution <- function(distribution=list(), lowerLimit=numeric(), upperLimit=numeric()) {
  if (upperLimit < lowerLimit | lowerLimit < 0.) {stop('not a valid distribution')}
  x <- list(distribution=distribution, lowerLimit=lowerLimit, upperLimit=upperLimit)
  class(x) <- "LayerDistribution"
  x 
}

mean.LayerDistribution <- function(d) {
  LAS(d$distribution, d$upperLimit) - LAS(d$distribution, d$lowerLimit)
}

LAS.LayerDistribution <- function(d, x){
  if (x>(d$upperLimit - d$lowerLimit)) {
    mean(d)
  } else {
    LAS(d$distribution, x+d$lowerLimit) - LAS(d$distribution, d$lowerLimit)
  }
}

SECM.LayerDistribution <- function(d, x){
  if (x > (d$upperLimit - d$lowerLimit)) {
    (SECM(d$distribution, d$upperLimit) - SECM(d$distribution, d$lowerLimit)) -
      2 * d$lowerLimit * (LAS(d$distribution, d$upperLimit) - LAS(d$distribution, d$lowerLimit))
  } else {
    (SECM(d$distribution, (x+d$lowerLimit)) - SECM(d$distribution, d$lowerLimit)) -
      2 * d$lowerLimit * (LAS(d$distribution, (x + d$lowerLimit)) - LAS(d$distribution, d$lowerLimit))
  }
}

cdf.LayerDistribution <- function(d, x){
  if (x>(d$upperLimit-d$lowerLimit)){
    1
  } else {
    cdf(d$distribution, (x+d$lowerLimit))
  }
}



# test ----------------------------------------------

# d <- 
#   MixedExponential(
#     weights = c(0.345741,0.332724,0.177363,0.106868,0.029365,0.005924,
#                 0.001421,0.000399,0.000155,0.000040),
#     means = c(2597,9185,34521,129180,591936,1679763,4272482,9716829,
#               23319137,100000000))
# 
# # discretization parameters
# N <- 2^12               #length
# h <- 750                #step-size
# limit <- 25000
# 
# dLimited <- LayerDistribution(d, 0, limit)
# 
# probs <- discretize(d, h, N)


# aggregate methods -------------------------------------

#Heckman-Meyer's probability generating function
HMpgf <- function(t, count, contagion)
{
  if (contagion == 0)
    exp(count*(t-1))
  else
  {
    alpha <- 1 / contagion
    beta <- count / alpha
    (1 - beta * (t-1))^(-alpha)
  }
}

#scaled & real-value inverse fast fourier transform
ifft <- function(x)
{
  temp <- Re(fft(x,inverse=TRUE))
  temp / sum(temp)
}