#the majority of the following code is inspired by
#the DataFrames.jl package (https://github.com/JuliaStats/DataFrames.jl)
#DataFrames.jl is licensed under the MIT License, a copy of which
#can be found in this package's repo in LICENSE.md.

#This file contains functions that are useful for crossjoining 
#arbitrary arrays, and are thus useful for building a tensor product
#structure. 

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

function rep(x::AbstractArray, times::Integer = 1, each::Integer = 1)
    res = similar(x, each * times * size(x,1), size(x,2))
    i = 1
    for jdx in 1:times
        for idx in 1:size(x,1)
            tmp = x[idx,:]
            for kdx in 1:each
                res[i,:] = tmp
                i += 1
            end
        end
    end
    return res
end

crossjoin(A::Array, B::Array) = hcat(rep(A, 1, size(B,1)), rep(B, size(A,1), 1))

crossjoin(A::Array...) = reduce(crossjoin, A)
