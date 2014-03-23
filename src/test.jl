#######################################################################################
include("basistries.jl")
include("rep.jl")
using BasisTries;

function tensor(arr::Array...)
	if length(arr) == 2
		return crossjoin(arr[1], arr[2])
	else
		tensor(crossjoin(arr[1], arr[2]), arr[3:end]...)
	end
end

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

function toTrie(A::Array)
	res = BasisTrie()
	for i=1:size(A,1)
		res[A[i,:]...] = i
	end
	return res
end	

function toDict(A::Array)
	res = Dict{Tuple,Int}()
	sizehint(res, size(A,1))
	for i=1:size(A,1)
		res[A[i,:]...] = i
	end
	return res
end	

tensor([1:3],[1:3])
bob = tensor(["A", "B", "C"],[1:100],[[1:4],[1:20]])

@time t = toTrie(bob);
@time d = toDict(bob);
println("d")
