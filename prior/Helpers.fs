namespace CATALINA.UDF.Library

open System
open System.Numerics
open System.Collections.Generic
open System.Text

[<AutoOpen>]
module GlobalHelpers = 
    // Save computed results by using an internal dictionary.
    // Note that memoize is inferred to have type
    // memoize: ('a -> 'b) -> ('a -> 'b)
    let memoize f = 
        let cache = Dictionary<_, _>()
        fun x -> 
            if cache.ContainsKey(x) then cache.[x]
            else 
                let res = f x
                cache.[x] <- res
                res
                
    let inline isNull value = obj.ReferenceEquals(value, null)

    let isTuple value = 
        match box value with
        | null -> false
        | _ -> Microsoft.FSharp.Reflection.FSharpType.IsTuple(value.GetType())

[<RequireQualifiedAccess>]
module Array = 
    let toIList(x : _ []) = x :> System.Collections.Generic.IList<_>
    let toString(A : float[]) = String.Join(";", A |> Array.map(sprintf "%g")) |> sprintf "[%s]"

[<RequireQualifiedAccess>]
module Array2D = 
    let cast<'a>(A : obj [,]) : 'a [,] = A |> Array2D.map unbox
    let flatten<'a>(A : 'a [,]) = A |> Seq.cast<'a>
    let getColumn c (A : _ [,]) = flatten A.[*, c..c] |> Seq.toArray
    let getColumns(A : _ [,]) = [| for i = 0 to A.GetLength(1) - 1 do yield A |> getColumn i |]  
    let getRow<'T> r (A : 'T [,]) = flatten<'T> A.[r..r, *] |> Seq.toArray
    let getRows<'T>(A : 'T [,]) = [| for i = 0 to A.GetLength(0) - 1 do yield A |> getRow<'T> i |]
    let toArray<'a>(A : obj [,]) = A |> Seq.cast<'a> |> Seq.toArray
    let ofColumnArray(A : 'a []) = Array2D.init (A.Length) 1 (fun i j -> A.[i])

[<RequireQualifiedAccess>]
module Float = 
    let toComplex(x : float) = Complex(x, 0.)

module Regex =
    open System.Linq
    open System.Text.RegularExpressions

    let (|Match|_|) pattern input = 
        let m = Regex.Match(input, pattern)
        if m.Success then 
            Some(List.tail [ for g in m.Groups -> g.Value ])
        else None

    let (|Split|_|) regex str = 
        let rec loop (acc : string list) (m : Match) = 
            if m.Success then loop ((m.Value) :: acc) (m.NextMatch())
            elif acc.Length = 0 then None
            else Some(List.rev acc)
        Regex(regex).Match(str) |> loop []

    let (|Float|_|) s = 
        match s with
        |Match "([0-9.eE+-]+)" [res] -> Some(Double.Parse(res))
        |_ -> None

    let (|FloatArray|_|) s =
        match s with
        |Match "\A\[(.*)\]" [Split "[0-9.eE+-]+" res] -> Some(res.Select(Double.Parse).ToArray())
        |_ -> None