#####################
#function definitions
#####################

#setClass("People",
#         representation=representation(
#           firstNames="character",
#          ages="numeric"),
#         validity=function(object) {
#           if (length(object@firstNames) != length(object@ages))
#             "'firstNames' and 'ages' must have same length"
#           else TRUE
#         })
#
#People = function(firstNames=character(), ages=numeric(), ...)
#  new("People", firstNames=firstNames, ages=ages, ...)

setClass("MixedExponential",
         representation=representation(
           weights="numeric",
           means="numeric"),
         validity=function(object) {
           if (length(object@weights) != length(object@means))
             "'weights' and 'means' must have same length"
           else TRUE
         })

weights <- c(0.8, 0.2)
means <- c(10, 1000)

MixedExponential <- function(weights=numeric(), means=numeric()) {
  x <- list(weights=weights, means=means)
  class(x) <- "MixedExponential"
  x 
}

x <- MixedExponential(weights, means)

#mean method
mean <- function(x) UseMethod("mean", x)
mean.MixedExponential <- function (x) {
  sum(x$weights * x$means)
}

#limited average severity method

MixedExponential <- 
  function(weights=numeric(), means=numeric(), ...)
  new("MixedExponential", weights=weights, means=means, ...)                

setMethod("mean", signature("MixedExponential"),
          function(dist = "MixedExponential")
          {
            sum(dist@weights * dist@means)
          })

x <- MixedExponential(weights=c(0.8,0.2), means=c(10,1000))




LAS <- function (weight, mean, x) 
{sum(weight*mean*(1-exp(-x/mean)))}

#mixed exponential cummulative density function
mixexp.cdf <- function(weight, mean, x) {
  1-sum(weight*exp(-x/mean))
}

#secant method
secant <- function(fun, x0, x1, tol=1e-07, niter=500){
  for ( i in 1:niter ) {
    x2 <- x1-fun(x1)*(x1-x0)/(fun(x1)-fun(x0))
    if (abs(fun(x2)) < tol)
      return(x2)
    x0 <- x1
    x1 <- x2
  }
  stop("exceeded allowed number of iteractions")
}

#mixed exponential inverse cummulative density function
mixexp.icdf <- 
  
  
#   function(weight, mean, prob) {
# 
#   f <- function(x){prob - mixexp.cdf(weight, mean, x)}
# 
#   tol <- 1e-4
#   n0 <- 100
#   
#   p1 <- 10000
#   
#   #need to find upper value
#   while (f(p1) > 0) {
#     p1 <- p1 * 2
#   }
# 
#   p0 <- p1
#   
#   #need to find lower value
#   while (f(p0) < 0) {
#     p0 <- p0 / 2
#   }
#   
#   n <- 1
#   
#   good.enough <- FALSE
# 
#   while ((good.enough==FALSE) & (n <= n0)) {
#   
#   q0 <- f(p0)
#   q1 <- f(p1)
#   
#   p <- p1 - q1 * (p1 - p0) / (q1 - q0)
#   q <- f(p)
#   
#   if (abs(p-p1) < tol){
#     good.enough <- TRUE
#   }  
#   
#   if (q*q1 < 0) {p0 <- p1}
# 
#   p1 <- p
#   n <- n+1
#   }
#   
#   if (n <= n0) {
#     p
#   }
# }


# /////generalized method of false position (Numerical Analysis, Burden & Faires)
# /////solves for f(x)=0 with initial guesses x= p0 & p1 
# /////tolerance = tol and maximum iterations = n0 
# //let SolveFP f p0 p1 tol n0 =
#   //  let rec pNext p0 p1 n =
#   //    if n > n0 then failwith "could not converge to a solution"
# //    let q0 = f p0
# //    let q1 = f p1
# //    let p = p1 - q1 * (p1 - p0) / (q1 - q0)
# //    let q = f p
# //    match (abs(p-p1) < tol), ((q*q1) < 0.) with
# //    |true, _      -> p
# //    |false, true  -> pNext p1 p (n+1) 
# //    |false, false -> pNext p0 p (n+1)
# //  pNext p0 p1 1  


#mixed exponential discretized probability density function
pdf <- function(weight, mean, n, h, maxLoss)
{
  LAS <- function (x) 
    {
     if ((x>maxLoss)&&(maxLoss != 0))
       {x <- maxLoss}
     sum(weight*mean*(1-exp(-x/mean)))
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