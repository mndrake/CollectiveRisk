NoOverdraft <- setRefClass("NoOverdraft",
                           contains = "Account",
                           methods = list(
                             withdraw = function(x) {
                               if (balance < x) stop("Not enough money")
                               balance <<- balance - x
                             }
                           )
)


# /// base class for all distributions
# [<AbstractClass>]
# type DistributionBase() as this = 
#   abstract Mean : float with get
# abstract LAS : float -> float
# abstract SECM : float -> float
# abstract LimVar : float -> float
# abstract CDF : float -> float
# override this.LimVar limit = this.SECM limit - this.LAS limit ** 2.
# 
# // discrete approximation using mean-matching method
# // notation based on S. Wang paper  
# // h=step size, N=total steps, X=severity    
# member this.Discretize h N = 
#   let X = 
#   [| for k = 0 to N do
#    yield this.LAS(h * float k) |]
# [| yield 1.0 - X.[1] / h
#  for k = 1 to N - 1 do
#  yield (2. * X.[k] - X.[k + 1] - X.[k - 1]) / h |]
# 
# member this.toIDistribution = this :> IDistribution
# interface IDistribution with
# member I.Mean = this.Mean
# member I.LAS limit = this.LAS limit
# member I.SECM limit = this.SECM limit
# member I.LimVar limit = this.LimVar limit
# member I.CDF limit = this.CDF limit
# member I.Discretize h N = this.Discretize h N
# member I.ToString() = this.ToString()