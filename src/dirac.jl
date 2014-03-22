module m

const lang = "\u27E8"
const rang = "\u27E9"

import Base.show

#####Basis##########################################################
abstract AbstractBasis

immutable Basis{T} <: AbstractBasis
	name::String
	bra_sym::String
	ket_sym::String
	label_dict::Dict{T, Int}
	label_arr::Array{T}
end

	#labels is of the form [label1=>index1,label2=>index2...]
	function Basis{T}(name::String, label_dict::Dict{T, Int}; bra_sym=lang, ket_sym=rang, errors=true) 
		if errors #errors==true will cause the function to validate indices of label_dict
			sorted_vals = sort(collect(values(label_dict)))
			test_vals = [1:last(sorted_vals)]
			missing = setdiff(test_vals, sorted_vals)
			if length(missing)!=0
				error("Can't form label array - the following indices are unaccounted for: \n", repr(missing))
			end
		end	
		label_arr = Array(Labels, last(sorted_vals))
		for (k,v) in label_dict
			label_arr[v] = k
		end
		Basis(name, bra_sym, ket_sym, label_dict, label_arr)
	end

	#label_arr is of the form {label1,label2...}
	function Basis(name::String, label_arr::Array; bra_sym=lang, ket_sym=rang)
		label_dict = Dict{eltype(label_arr),Int}()
		len=length(label_arr)
		sizehint(label_dict, len)
		for i=1:len
			label_dict[label_arr[i]] = i
		end
		Basis(name, bra_sym, ket_sym, label_dict, label_arr)
	end

immutable TensorBasis{B<:Basis, V<:Vector} <: AbstractBasis
	name::String
	label_dict::Dict{V, Int}
	label_arr::Array #memoization for future tensor products
	basis_arr::Vector{B}
	bra_sym::String
	ket_sym::String
end

	function TensorBasis{B<:Basis, V<:Vector}(name::String, label_dict::Dict{V, Int}, basis_arr::Array{B}; bra_sym=lang, ket_sym=rang, errors=true)
		if errors #errors==true will cause the function to validate properties label_dict
			if B==Basis  
				vec_type = Any
			else
				vec_type = eltype(basis_arr[1].label_arr)
			end
			if eltype(label_dict)[1]!=Vector{vec_type}	
	  			error("parameters of keys in label_dict must match parameters of Bases in basis_arr")
	  		end
			sorted_vals = sort(collect(values(label_dict)))
			test_vals = [1:last(sorted_vals)]
			missing = setdiff(test_vals, sorted_vals)
			if length(missing)!=0
				error("Can't form label array - the following indices are unaccounted for: \n", repr(missing))
			end
		end	
		width = length(first(label_dict)[1])
		label_arr = cell(length(label_dict), width)
		for (k,v) in label_dict
			for i=1:width
				label_arr[v, i] = k[i]
			end
		end
		TensorBasis(name, label_dict, label_arr, basis_arr, bra_sym, ket_sym)
	end

	function TensorBasis{B<:Basis}(name::String, label_arr::Array, basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		if B==Basis  
			vec_type = Any
		else
			vec_type = eltype(basis_arr[1].label_dict)[1]
		end
		label_dict = Dict{Vector{vec_type}, Int}()
		len=length(label_arr[:, 1])
		sizehint(label_dict, len)
		nest_arr = nest_rows(label_arr)
		for i=1:len
			label_dict[nest_arr[i]] = i
		end
		TensorBasis(name, label_dict, label_arr, basis_arr, bra_sym, ket_sym)
	end

#####BasisFunctions##########################################################
function components(B::TensorBasis)
	return B.basis_arr
end

function components(B::Basis)
	return B
end

function nest_rows{T}(arr::Array{T,2})
	len::Int = length(arr[:,1])
	flipped::Array{T,2} =arr.'
	result::Vector{Vector} = Array(Vector, len)
	for i=1:len
		result[i] = flipped[:,i]
	end
	return result
end

function vectensor(A::Array, B::Vector)
	width_A::Int = length(A[1, :])
	result::Array = cell(length(A[:,1])*length(B), width_A+1)
	count::Int = 1
	for i=1:length(A[:,1])
		for j=1:length(B)
			for k=1:width_A
				result[count, k]= A[i, k]
			end
			result[count, end]=B[j]
			count+=1
		end
	end
	return result
end

function vectensor(A::Array, B::Array)
	result = A;
	for i=1:length(B[1,:])
		result = vectensor(result,B[:,i])
	end
	return result
end

function tensor() #preventing ambiguities
	error("tensor() needs arguments of type A::Vector... or AbstractBasis...")
end

function tensor(A::Array...)
	if length(A)==2
		return vectensor(A[1], A[2])
	else
		tensor(vectensor(A[1],A[2]), A[3:end]...)
	end
end

function tensor(bases::AbstractBasis...; bra_sym=lang, ket_sym=rang)
	name = string(["_$(i.name)" for i in bases]...)
	return tensor(name[2:end], bases..., bra_sym=bra_sym, ket_sym=ket_sym)
end

function tensor(name::String, bases::AbstractBasis...; bra_sym=lang, ket_sym=rang)
	label_arr = tensor([i.label_arr for i in bases]...)
	return TensorBasis(name, label_arr, vcat([components(i) for i in bases]...);  bra_sym=bra_sym, ket_sym=ket_sym)
end

function tensor(name::String, bases::AbstractBasis...; bra_sym=lang, ket_sym=rang)
	label_arr = tensor([i.label_arr for i in bases]...)
	return TensorBasis(name, label_arr, vcat([components(i) for i in bases]...);  bra_sym=bra_sym, ket_sym=ket_sym)
end

function extract(f::Function, B::Basis, name=B.name)
	sub_arr = filter(f, B.label_arr)
	return Basis(name, sub_arr)
end

function extract(f::Function, B::TensorBasis, name=B.name)
	sub_arr = filter(f, collect(keys(B.label_dict)))
	sub_dict = Dict{Vector{Any}, Int}() 
	for i=1:length(sub_arr)#insertion sort FTW
		j=i
		while j>1 && B.label_dict[sub_arr[j-1]]>B.label_dict[sub_arr[j]]
			el = sub_arr[j]
			sub_arr[j] = sub_arr[j-1]
			sub_arr[j-1] = el		
		end
	end
	for i=1:length(sub_arr)
		sub_dict[sub_arr[i]] = i
	end
	return TensorBasis(name, sub_dict, B.basis_arr, errors=false)
end

# #####State##########################################################
type State
	label::String
	content::Dict{String, SparseMatrixCSC}
	bases::Dict{String, Basis}
end

# 	#########content only
	function State{S1<:String, S2<:String, M<:SparseMatrixCSC}(label::S1, content::Dict{S2,M})
		warn("Basis is implicit; assuming sparse coefficient arrays match given basis names")
		bases = Dict{String, Basis}()
		for i in getkeys(content)
			bases[i] = Basis(i, Dict{Any, Int}(), cell(length(content)))
		end
		State(label, content, bases_dict)
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

#####REPL Representation############################################
function show(io::IO, b::AbstractBasis)
	println("$(typeof(b))")
	println("Name: $(b.name)")
	println("$(length(b.label_arr[:,1])) Basis States:")
	for i=1:length(b.label_arr[:,1])
		println("| $(repr([b.label_arr[i, :]...])[2:end-1]) $(b.ket_sym)")
	end
end


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

