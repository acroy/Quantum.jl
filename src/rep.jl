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