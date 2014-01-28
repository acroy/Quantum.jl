module d

#export Basis, Converter
#export Basis, State, Converter
#####Basis##########################################################
type Basis{Labels<:Any}
	name::String
	labels::Dict{Labels, Int}
	length::Int
end
	#labels is of the form [label1=>index1,label2=>index2...]
	function Basis{S<:String, K<:Any}(name::S, labels::Dict{K, Int64})
		sorted_vals = sort(collect(values(labels)))
		max_index = last(sorted_vals)
		test_vals = Set([1:max_index]...)
		setdiff!(test_vals, Set(sorted_vals...))
		if length(test_vals)==0
			Basis(name, labels, max_index)
		else
			error("The following indices are unaccounted for: \n", string(sort([test_vals...])))
		end
	end

	#label_array is of the form {label1,label2...}
	function Basis{S<:String, A<:Array}(name::S, label_array::A)
		labels = Dict{Any,Int64}()
		labels["A"] = 1
		len=length(label_array)
		sizehint(labels, len)
		for i=1:len
			labels[label_array[i]] = i
		end
		Basis(name, labels, len)
	end
#####Tensor##########################################################
function vectensor(A::Vector, B::Vector)
	result_length = length(A)*length(B)
	result = {{} for i=1:result_length}
	count = 1
	for i in A
		for j in B
			push!(result[count], i...)
			push!(result[count], j...)
			count+=1
		end
	end
	return result
end

function tensor(A::Array)
	if length(A)==2
		return vectensor(A[1], A[2])
	else
		a=A[1]
		b=A[2]
		tensor({vectensor(a,b), A[3:end]...})
	end
end

function separate(A::Array)
	basis_num = length(A[1])
	bases = {}
	for i=1:basis_num
		push!(bases, unique({A[j][i] for j=1:length(A)}))
	end
	return bases
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

A = d.Basis("A", [1:100]);
B = d.Basis("B", [1:100]);
# coeffs1 = [1:100 1:100];
# coeffs2 = [100:200 100:200];
# conv = d.Converter(A, B, coeffs1);
# d.add_conversion!(A, B, sparse(coeffs2), conv);
dave = {{"A","B"},{1,2,3},{"#",";"},{4,5,6}};
bob = {"A"=>("B", "E"),
	   "B"=>("A", "C"),
	   "C"=>("B", "E", "D"),
	   "D"=>("C",),
	   "E"=>("A","F","C"),
	   "F"=>("E",)};
;
