#the majority of the following code is inspired by/pulled from 
#the DataFrames.jl package (https://github.com/JuliaStats/DataFrames.jl)
#I didn't want to force DataFrames as a dependency for just a few functions.
#DataFrames.jl is licensed under the MIT License, a copy of which
#can be found in this repo in LICENSE.md.

#This whole thing contains functions that are useful for crossjoining 
#arbitrary arrays, and are thus useful for building a tensor product
#structure. 

#At some point, this code should be refactored to utilize Base.Cartesian
#since, AFAIK, julia does not yet feature tail-call optimization. 

function rep{T <: Integer}(x::AbstractVector, lengths::AbstractVector{T})
    if length(x) != length(lengths)
        error("vector lengths must match")
    end
    res = similar(x, sum(lengths))
    i = 1
    for idx in 1:length(x)
        tmp = x[idx]
        for kdx in 1:lengths[idx]
            res[i] = tmp
            i += 1
        end
    end
    return res
end

function rep(x::AbstractVector, times::Integer = 1, each::Integer = 1)
    res = similar(x, each * times * length(x))
    i = 1
    for jdx in 1:times
        for idx in 1:length(x)
            tmp = x[idx]
            for kdx in 1:each
                res[i] = tmp
                i += 1
            end
        end
    end
    return res
end

function rep(x::AbstractVector; times::Integer = 1, each::Integer = 1)
    rep(x, times, each)
end

rep(x::Any, times::Integer) = fill(x, times)

function crossjoin(A::Array, B::Vector)
    r1, r2 = size(A, 1), size(B, 1)
    columns = [[rep(A[:,c], 1, r2) for c=1:size(A,2)],
               [rep(B[:,c], r1, 1) for c=1:size(B,2)]]
    hcat(columns...)
end

function crossjoin(A::Array, B::Array)
    result = A;
    for i=1:length(B[1,:])
        result = crossjoin(result,B[:,i])
    end
    return result
end

function crossjoin(arr::Array...)
    if length(arr) == 2
        return crossjoin(arr[1], arr[2])
    else
        crossjoin(crossjoin(arr[1], arr[2]), arr[3:end]...)
    end
end