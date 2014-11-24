namespace CATALINA.UDF.Library

open MathNet.Numerics

module CollectiveRisk = 
    open CATALINA.UDF.Library.Distributions

    /// discritization steps used by collective risk model
    [<Literal>]
    let DiscreteSteps =  1024

    /// tolerance used to determine if there is tail wrap in
    /// the collective risk model
    let DiscreteTolerance = 1e-04
    
    // methods based on 
    // "Aggregation of Correlated Risk Portfolios", Shaun S. Wang
    type Line = 
        { SevWeights : float[];
          SevDists : IDistribution[];
          FreqCount : float;
          FreqContagion : float;
          ExpectedLoss : float }
    
    type Aggregate = 
        { StepSize : float;
          CumProbs : float[];
          Charges : float[];
          EntryRatios : float[] }
          
    /// operations on Line                             
    [<CompilationRepresentation(CompilationRepresentationFlags.ModuleSuffix)>]
    [<RequireQualifiedAccess>]
    module Line = 
        let init() = { SevWeights = [||]; SevDists = [||]; FreqCount = 0.0; FreqContagion = 0.0; ExpectedLoss = 0.0 }
        let create(sevDists, counts, contagion) = 
            // filter out distributions with no counts
            let (_sevDists, _counts) = 
                Array.zip sevDists counts
                |> Array.filter(fun a -> snd a <> 0.0)
                |> Array.unzip
            let _freqCount = Array.sum _counts
            // create Line
            { SevWeights = _counts |> Array.map(fun c -> c / _freqCount);
              SevDists = _sevDists;
              FreqCount = _freqCount;
              FreqContagion = contagion;
              ExpectedLoss = 
                  _sevDists
                  |> Array.map Distribution.mean
                  |> Array.zip _counts
                  |> Array.sumBy(fun a -> fst a * snd a) }
        
        /// Discrete Fourier transform of the aggregate loss
        /// with step size h and total steps N
        let aggTransform h N (line : Line) = 
            // marginal probabilities
            let pdfs = line.SevDists |> Array.map(Distribution.discretize h N)
            // weighted probability
            let prob k = 
                pdfs
                |> Array.map(fun d -> d.[k])
                |> Array.zip line.SevWeights
                |> Array.sumBy(fun (p, w) -> p * w)
                |> fun r -> complex r 0.
            // FFT(aggregate loss)
            [|for k = 0 to N - 1 do yield prob k|]
            |> FFT.forward
            |> FFT.pgf line.FreqCount line.FreqContagion

    /// operations on Aggregate                                                    
    [<CompilationRepresentation(CompilationRepresentationFlags.ModuleSuffix)>]
    [<RequireQualifiedAccess>]
    module Aggregate = 
        let create N (lines : Line[]) = 
            // filter lines that do not have count=0
            let lines = lines |> Array.filter(fun l -> l.FreqCount <> 0.0)
            // expected loss for combined lines
            let expLoss = lines |> Array.sumBy(fun l -> l.ExpectedLoss)
            // tolerance used for testing discretization
            let tol = DiscreteTolerance
            // calculates discrete aggregate probability
            let aggProb N = 
                // initial step size
                let initH = (floor(expLoss * 12.0 / float N / 2000.0) + 1.0) * 2000.0
                let rec loop h = 
                    let dist = 
                        lines
                        |> Array.map(Line.aggTransform h N)
                        |> Array.fold FFT.convolution (FFT.zeroTransform N)
                        |> FFT.inverse
                        |> Array.map(fun x -> x.Real)
                    let actLoss = 
                        dist
                        |> Array.mapi(fun i p -> p * float i)
                        |> Array.sum
                        |> (*) h
                    match abs(actLoss / expLoss - 1.0) < tol || actLoss = expLoss with
                    | true -> (h, dist |> Array.toList)
                    | false -> loop(h * 2.0)
                loop initH
            // discretization results
            let stepSize, probs = aggProb N
            let cumProbs = 
                probs
                |> (List.scan (+) 0.0 >> List.tail)
                |> List.toArray
            let cumExpected = 
                probs
                |> (List.mapi(fun i p -> float i * p)
                    >> List.scan (+) 0.0
                    >> List.tail)
                |> List.toArray
            let maxExpected = 
                probs
                |> List.mapi(fun i p -> float i * p)
                |> List.sum
            // create Aggregate
            { StepSize = stepSize;
              CumProbs = cumProbs;
              Charges = 
                  Array.zip cumProbs cumExpected 
                  |> Array.mapi(fun i (p, e) -> (1.0 - (e + (1.0 - p) * float i) / maxExpected));
              EntryRatios = 
                  [|for i = 0 to N - 1 do
                        yield float i / maxExpected|] }

        // initialize with default total steps
        let init = create DiscreteSteps
        
        /// insurance charge
        let charge entryRatio (agg : Aggregate) = 
            Interpolate.Linear(agg.EntryRatios, agg.Charges).GetYofX(entryRatio)
        
        /// insurance savings
        let savings entryRatio (agg : Aggregate) = entryRatio + (agg |> charge entryRatio) - 1.0
        
        /// cummulative distribution function
        let cdf entryRatio (agg : Aggregate) = 
            Interpolate.Linear(agg.EntryRatios, agg.CumProbs).GetYofX(entryRatio)

        let maxEntryRatio (agg : Aggregate) = agg.EntryRatios |> fun a -> a.[a.Length-1]

        let findRoot agg f = 
            FindRoots.OfFunction((fun x -> f x), 0., maxEntryRatio agg, 1e-7, 100)

        let tryFindRoot agg f =
            try Some(findRoot agg f)
            with |_ -> None

        /// function to determine Entry Ratio for a given Savings
        let savingsEntryRatio sav (agg : Aggregate) =
            let f x = savings x agg - sav
            findRoot agg f

        /// function to determine Entry Ratio for a given Charge
        let chargeEntryRatio chg (agg : Aggregate) =
            let f x = charge x agg - chg
            findRoot agg f

        /// function to balance retro based on selected ELG, Charge & ER differences
        let balanceRetro chgDiff erDiff (agg : Aggregate) =
            let f x = (charge x agg) - (charge (x + erDiff) agg) - chgDiff  
            findRoot agg f
            |> fun r -> [| r; r + erDiff; charge r agg; charge (r + erDiff) agg |]

        /// function to determine entry ratio at min for min/no-max retro
        let minRetro chgMin (agg : Aggregate) =
            let f x = (charge x agg) - chgMin  
            findRoot agg f
            |> fun r -> [| r; (r + (charge r agg) - 1.) |]

        /// function to determine entry ratio & charge at max for max/no-min retro
        let maxRetro savMax (agg : Aggregate) =
            let f x = savings x agg - savMax
            findRoot agg f
            |> fun r -> [| r; charge r agg |]


        let tryMaxRetro savMax (agg : Aggregate) =
            match savMax < 0. with
            |true  -> None
            |false -> 
                let f x = savings x agg - savMax
                findRoot agg f
                |> fun r -> Some [| r; charge r agg |]

        /// function to balance retro based on selected ELG, Charge & ER differences
        let tryBalanceRetro chgDiff erDiff (agg : Aggregate) : obj array option =
            match (chgDiff >= 0.) && (chgDiff <= 1.) && (erDiff > 0.) with
            |true  ->
                // verify that desired charge difference is between max/min charge difference
                if (1. - charge erDiff agg) < chgDiff then
                    // cannot balance suggest max/no-min retro
                    Some [| "MAX" ; 0.; chargeEntryRatio(1. - chgDiff) agg; 1.; 1. - chgDiff |]
                elif (charge ((maxEntryRatio agg) - erDiff) agg - charge (maxEntryRatio agg) agg) > chgDiff then
                    // cannot balance suggest no-max/min retro
                    Some [| "MIN"; chargeEntryRatio chgDiff agg; 0.; chgDiff; 0 |]
                else
                    let f x = (charge x agg) - (charge (x + erDiff) agg) - chgDiff  
                    match (tryFindRoot agg f) with 
                    |Some r -> Some [| "MINMAX"; r; r + erDiff; charge r agg; charge (r + erDiff) agg |]
                    |None -> None
            |false ->
                None