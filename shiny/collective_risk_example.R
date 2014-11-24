source("CollectiveRiskS3.R")
library(distrEx)

# Collective Risk Model Example
# -----------------------------
# example using a mixed exponential distribution, but 
# approach could be applied to other distribution types

distribution <- 
  MixedExponential(
    weights = c(0.345741,0.332724,0.177363,0.106868,0.029365,0.005924,
                0.001421,0.000399,0.000155,0.000040),
    means = c(2597,9185,34521,129180,591936,1679763,4272482,9716829,
              23319137,100000000))

# discretization parameters
n <- 2^12               #length
h <- 750                #step-size

# Frequency Assumptions
count <- 10             #expected claim count
c <- 0.15               #frequency contagion factor

#Retention Limit
maxloss <- 1000000      #maximum per occurrence loss


#calculate aggregate probabilities (via FFT)
X <- pdf(distribution,n,h,maxloss)   #discrete severity probabilities
F_X <- fft(X)                        #transformed severity
F_S <- HMpgf(F_X,count,c)            #transformed aggregate loss
S <- ifft(F_S)                       #discrete aggregate probabilities


#
# aggregate loss cummulative probability and insurance charges
#

loss <- 0:(n-1) * h

#calculate insurance charges (entry ratios 0-10 by 0.01 increments)
distS <- DiscreteDistribution(prob = S, supp = 0:(n-1)*h)

E_S <- E(distS)
ER <- (0:1000)/100
charge <- sapply(ER, function(x){1-E(distS,upp=E_S*x)/E_S})

aggloss_graph <- function()
{  
  xpos1 <- seq(0,h*n,by=h*n/4)
  plot(loss,S,type = "h",main = "Aggregate Loss Distribution",ylab="Probability", xlab="Loss (in $m)",col = "dark blue", xaxt="n")
  axis(1, at=xpos1, labels=sprintf("%.2f", xpos1/1000000)) 
}

inscharge_graph <- function()
{
  xpos2 <- seq(0,10,by=1)
  plot(ER,charge, type="l",xaxt = "n", col = "red", lwd = 2, main = "Insurance Charges", xlab = "Entry Ratios")
  axis(1, at=xpos2, labels=sprintf("%.1f",xpos2))
}

#plot aggregate loss distribution & insurance charges
attach(mtcars)
par(mfrow=c(2,1)) 
aggloss_graph()
inscharge_graph()