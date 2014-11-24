namespace CATALINA.UDF.Library

open System

module Interpolate = 
    open MathNet.Numerics.Interpolation
    
    type Linear(xVals : float [], yVals : float []) = 
        member val xVals = xVals |> Array.toIList
        member val yVals = yVals |> Array.toIList
        member this.GetYofX(x : float) = Interpolate.Common(this.xVals, this.yVals).Interpolate(x)
    
    /// linear interpolation between points
    let LinearInterpolation2D(tableVals : float [,], rowVals, colVals, rowVal, colVal) = 
        let findIndex value array = 
            match (array |> Array.tryFindIndex(fun i -> i > value)) with
            | Some 0 -> 1
            | Some x -> x
            | None -> array.Length - 1
        let rIndex = rowVals |> findIndex rowVal
        let cIndex = colVals |> findIndex colVal
        let r1 = rowVals.[rIndex - 1]
        let r2 = rowVals.[rIndex]
        let c1 = colVals.[cIndex - 1]
        let c2 = colVals.[cIndex]
        let f11 = tableVals.[rIndex - 1, cIndex - 1]
        let f21 = tableVals.[rIndex, cIndex - 1]
        let f12 = tableVals.[rIndex - 1, cIndex]
        let f22 = tableVals.[rIndex, cIndex]
        f11 * ((r2 - rowVal) * (c2 - colVal) / ((r2 - r1) * (c2 - c1))) 
        + f21 * ((rowVal - r1) * (c2 - colVal) / ((r2 - r1) * (c2 - c1))) 
        + f12 * ((r2 - rowVal) * (colVal - c1) / ((r2 - r1) * (c2 - c1))) 
        + f22 * ((rowVal - r1) * (colVal - c1) / ((r2 - r1) * (c2 - c1)))
    
    type FrischCarlson(xVals, yVals) = 
        let x : float [] = xVals
        let y : float [] = yVals
        
        let d = 
            [|for i = 0 to x.Length - 2 do
                  yield (y.[i + 1] - y.[i]) / (x.[i + 1] - x.[i])|]
        
        let m = 
            [|yield d.[0]
              for k = 1 to x.Length - 2 do
                  yield (d.[k - 1] + d.[k]) / 2.0
              yield d.[x.Length - 2]|]
        
        let findIndex (value : float) (A : float []) = 
            let test x = 
                // check to see if data ascending or descending         
                if A.[0] > A.[1] then x < value
                else x > value
            A
            |> Array.tryFindIndex test
            |> function 
               | None -> A.Length - 2
               | Some x -> max (x - 1) 0
        
        do 
            for k = 0 to x.Length - 2 do
                if d.[k] = 0.0 then 
                    m.[k] <- 0.0
                    m.[k + 1] <- 0.0
                else 
                    let a = m.[k] / d.[k]
                    let b = m.[k + 1] / d.[k]
                    if a * a + b * b > 9.0 then 
                        let t = 3.0 / sqrt(a * a + b * b)
                        m.[k] <- t * a * d.[k]
                        m.[k + 1] <- t * b * d.[k]

        new (table : float[,]) =
                let vals = table |> Array2D.getColumns
                FrischCarlson(vals.[0], vals.[1])

        
        member this.GetYofX(xVal : float) = 
            let i = findIndex xVal x
            let h = x.[i + 1] - x.[i]
            let t = (xVal - x.[i]) / h
            let h00 = 2.0 * t ** 3.0 - 3.0 * t ** 2.0 + 1.0
            let h10 = t ** 3.0 - 2.0 * t ** 2.0 + t
            let h01 = -2.0 * t ** 3.0 + 3.0 * t ** 2.0
            let h11 = t ** 3.0 - t ** 2.0
            y.[i] * h00 + h * m.[i] * h10 + y.[i + 1] * h01 + h * m.[i + 1] * h11
        
        member this.Interpolate(xs : float []) = xs |> Array.map(this.GetYofX)
    
    type FrischCarlson2D(rowVals, colVals, tableVals) = 
        let rowVals : float [] = rowVals
        let colVals : float [] = colVals
        let table : float [,] = tableVals
        
        member this.Interpolate(r : float [], c : float []) = 
            table
            |> Array2D.getColumns
            |> Array.map(fun tableCol -> 
                   let fc = new FrischCarlson(rowVals, tableCol)
                   r |> Array.map(fc.GetYofX))
            |> array2D
            |> Array2D.getColumns
            |> Array.map(fun tableRow -> 
                   let fc = new FrischCarlson(colVals, tableRow)
                   c |> Array.map(fc.GetYofX))
            |> array2D
