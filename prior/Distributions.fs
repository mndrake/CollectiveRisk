namespace CATALINA.UDF.Library

open System
open System.Linq
open System.Text
open MathNet.Numerics

/// probability distributions --
/// all distributions inherit from distribution base class
/// to add new distributions implement DistributionBase and 
/// override methods Mean, LAS, SECM, and CDF
module Distributions = 
    [<Interface>]
    type IDistribution = 
        /// expected value
        abstract Mean : float with get
        /// limited average severity
        abstract LAS : float -> float
        /// limited second moment
        abstract SECM : float -> float
        /// limited variance
        abstract LimVar : float -> float
        /// cummulative distribution function
        abstract CDF : float -> float
        /// create a discrete distribution
        abstract Discretize : float -> int -> float []
        abstract ToString : unit -> string
    
    /// base class for all distributions
    [<AbstractClass>]
    type DistributionBase() as this = 
        abstract Mean : float with get
        abstract LAS : float -> float
        abstract SECM : float -> float
        abstract LimVar : float -> float
        abstract CDF : float -> float
        override this.LimVar limit = this.SECM limit - this.LAS limit ** 2.
        
        // discrete approximation using mean-matching method
        // notation based on S. Wang paper  
        // h=step size, N=total steps, X=severity    
        member this.Discretize h N = 
            let X = 
                [| for k = 0 to N do
                       yield this.LAS(h * float k) |]
            [| yield 1.0 - X.[1] / h
               for k = 1 to N - 1 do
                   yield (2. * X.[k] - X.[k + 1] - X.[k - 1]) / h |]
        
        member this.toIDistribution = this :> IDistribution
        interface IDistribution with
            member I.Mean = this.Mean
            member I.LAS limit = this.LAS limit
            member I.SECM limit = this.SECM limit
            member I.LimVar limit = this.LimVar limit
            member I.CDF limit = this.CDF limit
            member I.Discretize h N = this.Discretize h N
            member I.ToString() = this.ToString()
    
    type ExponentialDistribution(mean) = 
        inherit DistributionBase()
        member this.ExpectedSecondMoment = 2. * mean ** 2.
        override this.Mean = mean
        override this.LAS limit = mean * (1. - exp(-limit / mean))
        override this.SECM limit = 2. * mean ** 2. * (1. - (1. + limit / mean) * exp(-limit / mean))
        override this.CDF limit = 1. - exp(-limit / mean)
        override this.ToString() = sprintf "ExponentialDistribution(%g)" mean
    
    type MixedExponentialDistribution(weights, means) = 
        inherit DistributionBase()
        
        do 
            if (abs(weights |> Array.sum) - 1.0) > 0.0001 then invalidArg "weights" "must sum to 1"
            if means |> Array.exists(fun m -> m < 0.0) then 
                invalidArg "means" "means must be greater than 0."
        
        let sum f = 
            Array.zip weights means
            |> Array.filter(fun (_, m) -> m > 0.)
            |> Array.sumBy f
        
        member this.ExpectedSecondMoment = sum(fun (w, m) -> 2. * w * m ** 2.)
        override this.Mean = sum(fun (w, m) -> w * m)
        override this.LAS limit = sum(fun (w, m) -> w * m * (1. - exp(-limit / m)))
        override this.SECM limit = 
            sum(fun (w, m) -> 2. * w * m ** 2. * (1. - (1. + limit / m) * exp(-limit / m)))
        override this.CDF limit = 1. - sum(fun (w, m) -> w * exp(-limit / m))
        override this.ToString() = 
            sprintf "MixedExponentialDistribution(%s,%s)" (weights |> Array.toString) 
                (means |> Array.toString)
    
    type LogNormalDistribution(mu, sigma) = 
        inherit DistributionBase()
        let stdNormal = Distributions.Normal(0., 1.).CumulativeDistribution
        let invStdNormal = Distributions.Normal(0., 1.).InverseCumulativeDistribution
        member this.ICDF prob = exp(mu + sigma * invStdNormal(prob))
        override this.Mean = exp(mu + 0.5 * sigma ** 2.)
        override this.LAS limit = 
            exp(mu + 0.5 * sigma ** 2.) * stdNormal((log limit - mu - sigma ** 2.) / sigma) 
            + limit * (1. - stdNormal((log limit - mu) / sigma))
        override this.SECM limit = raise <| NotImplementedException()
        override this.CDF limit = stdNormal((log limit - mu) / sigma)
        override this.ToString() = sprintf "LogNormalDistribution(%g,%g)" mu sigma
    
    type ParetoDistribution(alpha, theta) = 
        inherit DistributionBase()
        override this.Mean = theta / (alpha - 1.)
        override this.LAS limit = 
            theta / (alpha - 1.) * (1. - (theta / (limit + theta)) ** (alpha - 1.))
        override this.SECM limit = raise <| NotImplementedException()
        override this.CDF limit = 1. - (theta / (limit + theta)) ** alpha
        override this.ToString() = sprintf "ParetoDistribution(%g,%g)" alpha theta
    
    /// distribution based on the mixture of two other distributions  
    type MixedDistribution(weights : float [], distributions : IDistribution []) = 
        inherit DistributionBase()
        
        do 
            if (abs(weights |> Array.sum) - 1.) > 1e-4 then invalidArg "weights" "must sum to 1"
        
        let sum f = 
            Array.zip weights distributions
            |> Array.filter(fun (w, _) -> w > 0.)
            |> Array.sumBy f
        
        override this.Mean = sum(fun (w, d) -> w * d.Mean)
        override this.LAS limit = sum(fun (w, d) -> w * d.LAS(limit))
        override this.CDF limit = sum(fun (w, d) -> w * d.CDF(limit))
        override this.SECM limit = sum(fun (w, d) -> w * d.SECM(limit))
        override this.ToString() = 
            sprintf "MixedDistribution(%s,%s)" (weights |> Array.toString) 
                (String.Join(";", distributions |> Array.map(fun d -> d.ToString())) 
                 |> sprintf "[%s]")
    
    /// distribution adjustment for excess of loss attachment
    type ExcessDistribution(dist : IDistribution, attachLimit : float) = 
        inherit DistributionBase()
        override this.LAS limit = dist.LAS(limit + attachLimit) - dist.LAS(attachLimit)
        override this.Mean = dist.Mean - dist.LAS(attachLimit)
        override this.SECM limit = raise <| NotImplementedException()
        override this.CDF limit = dist.CDF(limit + attachLimit)
        override this.ToString() = sprintf "ExcessDistribution(%s,%g)" (dist.ToString()) attachLimit
    
    /// distribution adjustment for layer of loss from lowerLimit to upperLimit
    type LayerDistribution(dist : IDistribution, lowerLimit, upperLimit) = 
        inherit DistributionBase()
        do 
            if (upperLimit < lowerLimit) || (lowerLimit < 0.) then 
                invalidArg "dist" "not a valid distribution"
        
        member this.ExpectedSecondMoment = 
            (dist.SECM upperLimit - dist.SECM lowerLimit) 
            - 2. * lowerLimit * (dist.LAS upperLimit - dist.LAS lowerLimit)
        override this.Mean = dist.LAS upperLimit - dist.LAS lowerLimit
        
        override this.LAS limit = 
            if limit > (upperLimit - lowerLimit) then this.Mean
            else dist.LAS(limit + lowerLimit) - dist.LAS lowerLimit
        
        override this.SECM limit = 
            if limit > (upperLimit - lowerLimit) then this.ExpectedSecondMoment
            else 
                (dist.SECM(limit + lowerLimit) - dist.SECM lowerLimit) 
                - 2. * lowerLimit * (dist.LAS(limit + lowerLimit) - dist.LAS lowerLimit)
        
        override this.CDF limit = 
            match limit > (upperLimit - lowerLimit) with
            | true -> 1.
            | false -> dist.CDF(limit + lowerLimit)
        
        override this.ToString() = 
            sprintf "LayerDistribution(%s,%g,%g)" (dist.ToString()) lowerLimit upperLimit
    
    /// linear model for allocated loss adjustment expense
    type AlaeBand = 
        { low : float
          high : float
          ratio : float }
        static member Create(lBound : float [], alaeRatio : float [], alaeAdj) = 
            let n = lBound.Length
            [| for i = 0 to n - 1 do
                   yield { low = lBound.[i]
                           high = 
                               if i = n - 1 then infinity
                               else lBound.[i + 1]
                           ratio = alaeRatio.[i] * alaeAdj } |]
    
    /// allocated loss adjustment treatment for limit
    type Alae = 
        | Included
        | Prorata
        | Unlimited
    
    /// distribution adjusted for allocated loss adjument treatment
    type AlaeDistribution(d : IDistribution, a : AlaeBand [], alae : Alae) = 
        inherit DistributionBase()
        let n = a.Length
        let cumLosses = a |> Array.map(fun x -> x.low)
        
        let cumAlaes = 
            a
            |> Array.scan (fun acc x -> acc + (x.high - x.low) * x.ratio) 0.
            |> fun array -> array.[0..a.Length - 1]
        
        let cumTotals = Array.zip cumLosses cumAlaes |> Array.map(fun (loss, alae) -> loss + alae)
        let lasLosses = a |> Array.map(fun x -> d.LAS x.low)
        
        let lasAlaes = 
            a
            |> Array.scan (fun acc x -> acc + (d.LAS x.high - d.LAS x.low) * x.ratio) 0.
            |> fun array -> array.[0..a.Length - 1]
        
        let meanTotal = a |> Array.sumBy(fun x -> (d.LAS x.high - d.LAS x.low) * (1. + x.ratio))
        let meanAlae = meanTotal - d.Mean
        
        let getLossLimit limit = 
            match alae with
            | Included -> 
                let idx = Array.FindLastIndex(cumTotals, fun x -> x < limit) |> max 0
                a.[idx].low + (limit - cumTotals.[idx]) / (1. + a.[idx].ratio)
            | Prorata | Unlimited -> limit
        
        override this.Mean = meanTotal
        override this.SECM limit = raise <| NotImplementedException()
        override this.CDF limit = d.CDF(getLossLimit limit)
        override this.LAS limit = 
            let lossLAS = d.LAS(getLossLimit limit)
            
            let alaeLAS = 
                match alae with
                | Included -> 
                    let idx = Array.FindLastIndex(cumTotals, fun x -> x < limit) |> max 0
                    lasAlaes.[idx] + (lossLAS - lasLosses.[idx]) * a.[idx].ratio
                | Prorata -> 
                    let idx = Array.FindLastIndex(cumLosses, fun x -> x < limit) |> max 0
                    lasAlaes.[idx] + (lossLAS - lasLosses.[idx]) * a.[idx].ratio
                | Unlimited -> meanAlae
            lossLAS + alaeLAS
    
    /// parser for IDistribution classes
    module Language = 
        open Regex
        /// active pattern to parse distributions from a string input
        let rec (|Distribution|_|) str : IDistribution option = 
            match str with
            // mixed exponential distribution
            | Match "\AMixedExponentialDistribution\((.*)\,(.*)\)" 
              [ FloatArray weights; FloatArray means ] -> 
                MixedExponentialDistribution(weights, means) :> IDistribution |> Some
            | _ -> None
    
    /// functional methods for the distribution class
    [<RequireQualifiedAccess>]
    module Distribution = 
        /// expected value
        let mean(d : IDistribution) = d.Mean
        
        /// limited average severity
        let las limit (d : IDistribution) = d.LAS limit
        
        /// cumulative distribution function        
        let cdf limit (d : IDistribution) = d.CDF limit

        /// excess ratio
        let excessRatio limit (d : IDistribution) = 1.0 - (d.LAS limit / d.Mean)
        
        /// discrete approximation using mean-matching method        
        let discretize h N (d : IDistribution) = d.Discretize h N
        
        /// parses a string to an IDistribution option
        let tryParse str = 
            match str with
            | Language.Distribution d -> d |> Some
            | _ -> None
        
        /// parses a string to an IDistribution otherwise return error
        let parse str = 
            tryParse str |> function 
            | Some d -> d
            | None -> invalidArg "str" "could not parse string to IDistribution"
        
        let toString(d : IDistribution) = d.ToString()