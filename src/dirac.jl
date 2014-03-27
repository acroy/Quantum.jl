module m
include("rep.jl")

const lang = "\u27E8"
const rang = "\u27E9"

import Base.show,
	   Base.getindex,
	   Base.ndims,
	   Base.size,
	   Base.length,
	   Base.slice,
	   Base.+,
	   Base.-,
	   Base.*,
	   Base.in,
	   Base.setdiff
#####Utility##########################################################
function todict(A::Array; d=Dict{eltype(A),Int}())
	sizehint(d, size(A,1))
	for i=1:size(A,1)
		d[A[i,:]...] = i
	end
	return d
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

tensor()=error("tensor() needs an argument of the form A::Array...")

function tensor(arr::Array...)
	if length(arr) == 2
		return crossjoin(arr[1], arr[2])
	else
		tensor(crossjoin(arr[1], arr[2]), arr[3:end]...)
	end
end

#####Basis##########################################################
abstract AbstractBasis

immutable Basis{T} <: AbstractBasis
	name::String
	bra_sym::String
	ket_sym::String
	label_map::Dict{T, Int}
	label_arr::Vector{T}
end

	#label_arr is of the form [label1,label2...]
	function Basis(name::String, label_arr::Vector; bra_sym=lang, ket_sym=rang)
		label_arr = unique(label_arr)
		label_map = todict(label_arr)
		Basis(name, bra_sym, ket_sym, label_map, label_arr)
	end

immutable TensorBasis{B<:Basis} <: AbstractBasis
	name::String
	label_map::Dict{Tuple, Int}
	label_arr::Array #memoization for future tensor products
	basis_arr::Vector{B}
	bra_sym::String
	ket_sym::String
end

	function TensorBasis{B<:Basis}(name::String, label_arr::Array, basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		label_map = todict(label_arr, d=Dict{Tuple, Int}())
		TensorBasis(name, label_map, label_arr, basis_arr, bra_sym, ket_sym)
	end

	function TensorBasis{B<:Basis}(name::String, basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		label_arr = tensor([i.label_arr for i in basis_arr]...)
		TensorBasis(name, label_arr, basis_arr, bra_sym=bra_sym, ket_sym=ket_sym)
	end

	function TensorBasis{B<:Basis}(basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		name = string(["_$(i.name)" for i in basis_arr]...)
		TensorBasis(name[2:end], basis_arr, bra_sym=bra_sym, ket_sym=ket_sym)
	end

	function tensor(bases::AbstractBasis...)
		TensorBasis(vcat([components(i) for i in bases]...))
	end

#####BasisFunctions##########################################################
function components(B::TensorBasis)
	return B.basis_arr
end

function components(B::Basis)
	return B
end

function extract(f::Function, A::Array)
    res = Array(eltype(A), 0)
    for i=1:size(A,1)
	    if f(A[i, :])
            push!(res, A[i,:])
        end
    end
    return vcat(res...)
end

function extract(f::Function, B::Basis; name=B.name)
	sub_arr = filter(f, B.label_arr)
	return Basis(name, sub_arr)
end

function extract(f::Function, B::TensorBasis; name=B.name)
	sub_arr = extract(f, B.label_arr)
	return TensorBasis(name, sub_arr, B.basis_arr, bra_sym=B.bra_sym, ket_sym=B.ket_sym)
end

# #####State##########################################################
type State
	label::String
	content::Dict{String, SparseMatrixCSC}
	basis::Dict{String, Basis}
end

# 	function State{S<:String, N<:Number}(coeff_dict::Dict{S, Array{N,1}})
# 		warn("No basis provided! Assuming coefficient arrays match input bases")
# 		bases_dict = Dict{String, Dict{Any, Int64}}()
# 		new_coeffs = Dict{String, SparseMatrixCSC}()
# 		for i in keys(coeff_dict)
# 			bases_dict[i] = Dict{Any, Int64}()
# 			new_coeffs[i] = sparse(coeff_dict[i])
# 		end
# 		State(new_coeffs, bases_dict)
# 	end
	###########coeff_dict, bases_dict input	
	# function State{N<:Number, K<:Any, V<:Integer}(coeff_dict::Dict{String, Array{N, Int64}}, bases_dict::Dict{String, Dict{K, V}})
	# 	warn("Assuming coefficient arrays match input bases")
	# 	new_coeffs = Dict{String, SparseMatrixCSC}()
	# 	for i in keys(coeff_dict)
	# 		new_coeffs[i] = sparse(coeff_dict[i])
	# 	end
	# 	State(new_coeffs, bases_dict)
	# end
	# #########Single Basis Input
	# function State{N<:Number}(basis::Basis, coeffs::Array{N})
	# 	warn("Assuming coefficient array matches input basis")
	# 	State([basis.name=>sparse(coeffs)], [basis.name=>basis.labels])
	# end
	# function State{N<:Number}(basis::Basis, coeffs::SparseMatrixCSC{N})
	# 	warn("Assuming sparse coefficient array matches input basis")
	# 	State([basis.name=>coeffs], [basis.name=>basis.labels])
	# end
	# function State{K<:Any, V<:Number}(basis::Basis, coeffs::Dict{K, V})
	# 	coeff_array = zeros(Number, basis.max_index)
	# 	for i in coeffs
	# 		coeff_array[basis.labels[i[1]]] = i[2]
	# 	end
	# 	State([basis.name=>coeff_array], [basis.name=>basis.labels])
	# end

	#Multiple Bases Input
	# function State{K<:Any, V<:Number}(basis::Array{Basis}, coeffs::Dict{K, V})
	# 	basis_dict = 
	# end
	# function State{K<:Any, V<:Number}(basis::(dirac.Basis...,), coeffs::Dict{K, V})
	# 	basis_dict = 
	# ends
#####StateFunctions##########################################################

function magnitude(A::Number...)
	if length(A)==2
		return hypot(A[1], A[2])
	end
	return magnitude(hypot(A[1], A[2]), A[3:end]...)
end

magnitude{N<:Number}(A::Vector{N}) = magnitude(A...)

function normalize{N<:Number}(A::Vector{N})
	return (1/magnitude(A))*A
end

#####Converter##########################################################

# type Converter
# 	neighbors::Dict{String, (String...,)}
# 	edges::Dict{String, SparseMatrixCSC}
# end

# 	function Converter()
# 		Converter(Dict{String, (String...,)}(), Dict{String, SparseMatrixCSC}())
# 	end

#   	function Converter{S<:String}(b1::S, b2::S, coeffs::SparseMatrixCSC)
#   		neighbors = Dict{String, (String...,)}()
#   		neighbors[b1]=(b2,)
#   		neighbors[b2]=(b1,)
  		
#   		edges = Dict{String, SparseMatrixCSC}()
#   		edges["$b1_$b2"]=coeffs
  		
#   		Converter(neighbors, edges)
#   	end
  	
#   	function Converter(b1::Basis, b2::Basis, coeffs::SparseMatrixCSC)
#   		Converter(b1.name, b2.name, coeffs)
#   	end

#   	function Converter(b1::Basis, b2::Basis, coeffs::Matrix)
#   		Converter(b1.name, b2.name, sparse(coeffs))
#   	end
#   	###Adding conversions
#   	function add_conversion!{S<:String}(b1::S, b2::S, coeffs::SparseMatrixCSC, conv::Converter)
#   		if !in(b1, collect(keys(conv.neighbors)))
#   			conv.neighbors[b1] = (b2)
#   		elseif !in(b1, conv.neighbors[b2])
#   			push!(conv.neighbors[b2], b1)
#   		end

#   		if !in(b2, collect(keys(conv.neighbors)))
#   			conv.neighbors[b2] = (b1)
#   		elseif !in(b2, conv.neighbors[b1])
#   			push!(conv.neighbors[b1], b2)
#   		end

#   		conv.edges["$b1_$b2"]=coeffs;
#   	end

#   	function add_conversion!(b1::Basis, b2::Basis, coeffs::SparseMatrixCSC, conv::Converter)
#   		add_conversion!(b1.name, b2.name, coeffs, conv)
#   	end
#   	function add_conversion!(b1::Basis, b2::Basis, coeffs::Matrix, conv::Converter)
#   		add_conversion!(b1.name, b2.name, sparse(coeffs), conv)
#   	end
#   	###Performing conversions

#   	function generate_conv{S<:String}(start::S, dest::S, conv::Converter)
#   		if start == dest
#   			error("Source and target bases are the same")
#   		end
#   		if !in(start, collect(keys(conv.neighbors))) || !in(dest, collect(keys(conv.neighbors)))
#   			error("Either the source or target basis is missing from the given Converter")
#   		end
#   		init = conv.neighbors[start][1]
#   		coeffs = spones(conv.edges["$start_$init"])
#   		dest_reached = false
#   		visited = {}
#   		current = start
#   		while !dest_reached
#   			shared = intersect(conv.neighbors[current], conv.neighbors[dest]) 
#   			if !isempty(shared)
# 				coeffs = coeffs * conv.edges["$current_$(shared[1])"] * conv.edges["$(shared[1])_$dest"] 
#   				add_conversion!(start, dest, coeffs, conv)
#   				dest_reached=true
#   			end
#   			push!(visited, current)
#   			current = setdiff(conv.neighbors[current], visited)[1]
#   		end
#   	end

#   	function convert{S<:String}(state::State, new_basis::S, conv::Converter)
#   		# old_bases = collect(keys(state.contents))
#   		for i in state.contents
# 	  		try #if conversion already directly exists
# 	  			state.contents[new_basis] = conv.edges["$(i[1])_$(new_basis)"]*i[2]
# 	  			delete!(state.contents, i[1])
# 	  		catch

# 	  		end
# 	  	end


#   	end



#####Base.function Overloading############################################
function show(io::IO, b::AbstractBasis)
	println("$(typeof(b))")
	println("Name: $(b.name)")
	println("$(length(b.label_arr[:,1])) Basis States:")
	for i=1:length(b.label_arr[:,1])
		println("| $(repr([b.label_arr[i, :]...])[2:end-1]) $(b.ket_sym)")
	end
end

size(b::AbstractBasis, x::Int...) = size(b.label_arr, x...)
length(b::AbstractBasis) = size(b, 1)
ndims(b::AbstractBasis) = ndims(b.label_arr)
in(a, b::Basis)=in(a, collect(keys(b.label_map)))
in(a, b::TensorBasis)=in(tuple(a...), collect(keys(b.label_map)))
getindex(b::TensorBasis, x::Int) = TensorBasis("$(b.name)_$x", b.label_arr[x,:], b.basis_arr, bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::TensorBasis, x::Range1{Int}) = TensorBasis("$(b.name)_$(x[1]) to $(b.name)_$(last(x))", b.label_arr[x,:], b.basis_arr, bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::Basis, x::Int) = Basis("$(b.name)_$x", b.label_arr[x,:], bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::Basis, x::Range1{Int}) = Basis("$(b.name)_$(x[1]) to $(b.name)_$(last(x))", b.label_arr[x,:], bra_sym=b.bra_sym, ket_sym=b.ket_sym)
+(a::Basis,b::Basis) = Basis("$(a.name)+$(b.name)", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym)
+(a::TensorBasis,b::TensorBasis) = size(a,2)==size(b,2) ? TensorBasis("$(a.name)+$(b.name)", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym) : error("dimension mismatch; label lengths differ")
-(a::Basis,b::Basis) = extract(x->!in(x,b), a, name="$(a.name)-$(b.name)")
-(a::TensorBasis, b::TensorBasis) = size(a,2)==size(b,2) ? extract(x->!in(x,b),a, name="$(a.name)-$(b.name)") : error("dimension mismatch; label lengths differ")
*(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
setdiff(a::AbstractBasis,b::AbstractBasis) = a-b

end #module
#####Tests##########################################################
using m

function main()
	A = m.Basis("A", [1:3])
	B = m.Basis("B", ["$i" for i=1:3])
	C = m.Basis("C", ["%", "\$", "#"])
	D = m.tensor(A, B, C)
	return (A,B,C,D)
end

res = main();
a = res[1]
b = res[2]
c = res[3]
d = res[4]
println("")

