#########################################
#Mixed Exponential class (S3) and methods
#########################################

#class constructor
MixedExponential <- function(weights=numeric(), means=numeric()) {
  x <- list(weights=weights, means=means)
  class(x) <- "MixedExponential"
  x 
}

#mean method
mean <- function(d) UseMethod("mean", d)
mean.MixedExponential <- function (d) {
  sum(d$weights * d$means)
}

#limited average severity method
LAS <- function(d, ...) UseMethod("LAS", d)
LAS.MixedExponential <- function (d, x) {
  sum(d$weights*d$means*(1-exp(-x/d$means)))
}

#cummulative density function
cdf <- function(d, ...) UseMethod("cdf", d)
cdf.MixedExponential <- function(d, x) {
  1-sum(d$weights*exp(-x/d$means))
}

#inverse cummulative density function
icdf <- function(d, ...) UseMethod("icdf", d)
icdf.MixedExponential <- function(d,p) {
  f <- function(x){p-cdf(d,x)}
  xmax <- 100000
  while (f(xmax)>0) {xmax <- xmax * 10}
  uniroot(f,c(0,xmax))$root
}

#mixed exponential discretized probability density function
pdf <- function(d, ...) UseMethod("pdf", d)
pdf.MixedExponential <- function(d, n, h, maxLoss){
  LAS <- function (x) 
    {
     if ((x>maxLoss)&&(maxLoss != 0))
       {x <- maxLoss}
     sum(d$weights*d$means*(1-exp(-x/d$means)))
    }
  X<-rep(0,times=n)
  X[1] <- 1 - LAS(h)/h
  for (k in 2:n)
  {
    X[k] <- (2*LAS(h*(k-1))-LAS(h*k)-LAS(h*(k-2)))/h
  }
  X
}

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