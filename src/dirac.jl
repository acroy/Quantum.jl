module d

#####Basis##########################################################
abstract AbstractBasis

type Basis{S<:String, Labels<:Any} <: AbstractBasis
	name::S
	label_dict::Dict{Labels, Int}
	label_array::Array
end

	#labels is of the form [label1=>index1,label2=>index2...]
	function Basis{S<:String, Labels<:Any}(name::S, label_dict::Dict{Labels, Int64})
		sorted_vals = sort(collect(values(label_dict)))
		max_index = last(sorted_vals)
		test_vals = Set([1:max_index]...)
		setdiff!(test_vals, Set(sorted_vals...))
		if length(test_vals)==0
			label_array = Array(Labels, max_index)
			for i=1:max_index
				label_array[i] = label_dict[i]
			end
			Basis(name, label_dict, label_array)
		else
			error("The following indices are unaccounted for: \n", string(sort([test_vals...])))
		end
	end

	#label_array is of the form {label1,label2...}
	function Basis{S<:String, A<:Array}(name::S, label_array::A)
		label_dict = Dict{Any,Int}()
		len=length(label_array)
		sizehint(label_dict, len)
		for i=1:len
			label_dict[label_array[i]] = i
		end
		Basis(name, label_dict, label_array)
	end

type TensorBasis{S<:String, B<:Basis} <: AbstractBasis
	name::S
	label_dict::Dict{Array{Any}, Int}
	label_array::Array
	basis_array::Array{B}
end

	function TensorBasis{S<:String, B<:AbstractBasis}(name::S, label_dict::Dict{Array, Int}, basis_array::Array{B})
		label_array = Array(Any, length(label_dict))
		for i in label_dict
			label_array[i[2]] = i[1]
		end
		TensorBasis(name, label_dict, label_array, basis_array)
	end

	function TensorBasis{S<:String, B<:AbstractBasis}(name::S, label_array::Array, basis_array::Array{B})
		label_dict = Dict{Array{Any}, Int}()
		len=length(label_array[:, 1])
		sizehint(label_dict, len)
		for i=1:len
			label_dict[label_array[i]] = i
		end
		TensorBasis(name, label_dict, label_array, basis_array)
	end

#####TensorFunctions##########################################################

#  function vectensor(A::Vector, B::Vector)
# 	result = cell(length(A)*length(B))
# 	count=1
# 	for i=1:length(A)
# 		for j=1:length(B)
# 			result[count]= {A[i], B[j]}
# 			count+=1
# 		end
# 	end
# 	return result
# end

function vectensor(A::Array{Int}, B::Vector{Int})
	width_A = length(A[1, :])
	result = Array(Int, length(A[:,1])*length(B), width_A+1)
	count=1
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

function vectensor(A::Array{Int}...)
	if length(A)==2
		return vectensor(A[1], A[2])
	else
		vectensor(vectensor(A[1],A[2]), A[3:end]...)
	end
end

function map_inds(inds::Array{Int}, arr::Array)
	result = Array(Any, length(arr))
	for i=1:length(inds)
		result[i] = arr[i][inds[i]]
	end
	return result
end

function tensor(A::Array{Any})
	lens = [length(i) for i in A]
	result = cell(*(lens...))
	inds_array = vectensor([[1:i] for i in lens]...)
	for i=1:length(inds_array[:,1])
		result[i] = map_inds(inds_array[i, :], A)
	end
	return result
end

# function tensor()
# 	error("tensor() needs arguments of type A::Vector... or B<:AbstractBasis...")
# end

# function tensor(A::Vector...)
# 	if length(A)==2
# 		return vectensor(A[1], A[2])
# 	else
# 		tensor({vectensor(A[1],A[2]), A[3:end]...}...)
# 	end
# end

# function tensor{S<:String, B<:AbstractBasis}(name::S, bases::B...)
# 	tensor_array = tensor({i.label_array for i in bases}...)
# 	return TensorBasis(name, tensor_array, [bases...]) 
# end

# function tensor{B<:AbstractBasis}(bases::B...)
# 	name = ["_$(i.name)" for i in bases]
# 	return tensor(string(name...)[2:end], bases...)
# end

function separate(labels::Vector)
	basis_num = length(labels[1])
	bases = {}
	for i=1:basis_num
		push!(bases, unique({labels[j][i] for j=1:length(labels)}))
	end
	return bases
end

function components(B::TensorBasis)
	return B.basis_array
end
# #####State##########################################################
# type State
# 	label::String
# 	content::Dict{String, SparseMatrixCSC}
# 	bases::Dict{String, Basis}
# end

# 	#########content only
# 	function State{S1<:String, S2<:String, M<:SparseMatrixCSC}(label::S1, content::Dict{S2,M})
# 		warn("No basis provided! Assuming sparse coefficient arrays match given basis names")
# 		bases = Dict{String, Basis}()
# 		for i in keys(content)

# 			bases[i] = Basis(i, )
# 		end
# 		State(coeff_dict, bases_dict)
# 	end

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

#end module
end
#####Tests##########################################################
using d
inds = [1:4]
bob = {[1:5], [6:11], [12:18], [19:24]}

# println("Creating bases")
# tic();
# A = d.Basis("A", [1:100]);
# B = d.Basis("B", ["$i" for i=1:100]);
# C = d.Basis("C", ["%", "\$", "#"]);
# toc()
# println("Taking tensor product of bases")
# tic();
# D = d.tensor(A, B, C);
# toc()

# coeffs1 = [1:100 1:100];
# coeffs2 = [100:200 100:200];
# conv = d.Converter(A, B, coeffs1);
# d.add_conversion!(A, B, sparse(coeffs2), conv);
# dave = {{"A","B"},{1,2,3},{"#",";"},{4,5,6}};
# bob = {"A"=>("B", "E"),
# 	   "B"=>("A", "C"),
# 	   "C"=>("B", "E", "D"),
# 	   "D"=>("C",),
# 	   "E"=>("A","F","C"),
# 	   "F"=>("E",)};
# ;
