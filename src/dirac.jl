module d
include("rep.jl")
#includes,imports,consts#############################
const lang = "\u27E8"
const rang = "\u27E9"
const otimes = "\u2297"
import Base.show,
	   Base.repr,
	   Base.norm,	
	   Base.getindex,
	   Base.setindex!,
	   Base.ndims,
	   Base.size,
	   Base.length,
	   Base.slice,
	   Base.(.+),
	   Base.(.^),
	   Base.(.-),
	   Base.^,
	   Base.*,
	   Base.in,
	   Base.setdiff,
	   Base.get,
	   Base.!,
	   Base.exp,
	   Base.map,
	   Base.map!,
	   Base.filter,
	   Base.isequal,
	   Base.copy,
	   Base.hash,
	   Base.isequal,
	   Base.endof,
	   Base.start,
	   Base.find,
	   Base.trace
#Utility######################################
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

#AbstractTypes###################################
abstract DiracObj

abstract BraKet <: DiracObj
abstract Bra <: BraKet
abstract Ket <: BraKet
!(K::Type{Ket}) = Bra
!(B::Type{Bra}) = Ket
abstract AbstractStateRep{K<:BraKet} <: DiracObj
abstract AbstractBasis{K<:BraKet} <: DiracObj

#Basisless States###################################

immutable State{K<:BraKet} <: DiracObj
  label::Vector
  kind::Type{K}
end

State(label::Vector) = State(label, Ket)
State{K<:BraKet}(label, kind::Type{K}=Ket) = State([label], kind)
State{K<:BraKet}(label...; kind::Type{K}=Ket) = State([label...], kind)

function statevec{K<:BraKet}(v::Vector, kind::Type{K}=Ket)
	svec = Array(State{kind}, length(v))
	for i=1:length(v)
		svec[i] = State(v[i], kind)
	end
	return svec
end

function statevec{K<:BraKet}(arr::Array, kind::Type{K}=Ket)
	svec = Array(State{kind}, size(arr,1))
	for i=1:size(arr, 1)
		svec[i] = State(vec(arr[i,:]), kind)
	end
	return svec
end

tensor() = nothing
tensor{K<:BraKet}(s::State{K}...) = State(vcat([i.label for i in s]...), K)
*{K<:BraKet}(s1::State{K}, s2::State{K}) = tensor(s1, s2)
*(s1::State{Bra}, s2::State{Ket}) = s1.label==s2.label ? 1 : 0 

statejoin{S<:State}(v::Vector{S}) = tensor(v...)
statejoin{S<:State}(v::Vector{S}...) = broadcast(tensor, v...)
function statejoin{S<:State}(state_arr::Array{S}) 
	result = statejoin(state_arr[:,1], state_arr[:,2])
	for i=3:size(state_arr, 2)
		result = statejoin(result, state_arr[:,i])
	end
	return result
end

tensor{S<:State}(state_arrs::Array{S}...) = statejoin(crossjoin(state_arrs...))
separate(s::State) = statevec(s.label)
separate{S<:State}(v::Vector{S}) = hcat(map(separate, v)...).'

size(s::State) = size(s.label)
ndims(s::State) = 1
isequal(a::State,b::State) = a.label==b.label && a.kind==b.kind
hash(a::State) = hash(a.label)+hash(a.kind)

kind(s::State) = s.kind
ctranspose(s::State) = State(s.label, !s.kind)
getindex(s::State, x) = s.label[x]
setindex!(s::State, y, x) = setindex!(s.label, y, x)
endof(s::State) = endof(s.label)
repr(s::State{Bra}, extra="") = isempty(s.label) ? "$lang #undef$extra |" : "$lang $(repr(s.label)[2:end-1])$extra |"
repr(s::State{Ket}, extra="") = isempty(s.label) ? "| #undef$extra $rang" : "| $(repr(s.label)[2:end-1])$extra $rang"
label(s::State) = s.label

function show(io::IO, s::State)
	print(io, repr(s))
end

#Bases#############################################

function maplabels{K<:BraKet}(svec::Vector{State{K}})
	dict = Dict{Vector, Int}()
	sizehint(dict, length(svec))
	for i=1:length(svec)
		dict[svec[i].label] = i
	end
	return dict
end	

immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label
	states::Vector{State{K}}
	label_map::Dict{Vector, Int}
end

Basis{K<:BraKet}(label, states::Vector{State{K}}) = Basis(label, unique(states), maplabels(states))
Basis{K<:BraKet}(label, states::State{K}...) = Basis(label, vcat(states...))
Basis(label, label_vec::Vector) = Basis(label, statevec(label_vec))														

immutable TensorBasis{K<:BraKet} <: AbstractBasis{K}
	bases::Vector{Basis{K}}
	states::Vector{State{K}}
	label_map::Dict{Vector, Int}
end

function TensorBasis{K<:BraKet}(bases::Vector{Basis{K}}, states::Vector{State{K}})
	states = unique(states)
	label_map = maplabels(states)
	TensorBasis(bases, states, label_map)
end

function tensor{K<:BraKet}(bases::AbstractBasis{K}...)
	states = tensor([i.states for i in bases]...)
	TensorBasis(vcat([separate(i) for i in bases]...), states)
end

separate(b::TensorBasis) = b.bases
separate(b::Basis)=b
isequal(a::AbstractBasis,b::AbstractBasis) = a.states==b.states

label(b::Basis) = b.label
function label(b::TensorBasis)
	labels = [label(i) for i in b.bases]
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str $otimes $(labels[i])"
	end
	return str
end

function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) $(label(b)):")
	for i in b.states
		println(repr(i))
	end
end

ctranspose(b::Basis) = Basis(b.label, map(ctranspose, b.states), b.label_map)
ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states), b.label_map)

states(b::AbstractBasis) = b.states
filter(f::Function, b::Basis) = Basis(b.label, filter(f, b.states))
filter(f::Function, b::TensorBasis) = TensorBasis(b.bases, filter(f, b.states))
getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)
endof(b::AbstractBasis) = endof(b.states)
length(b::AbstractBasis) = length(b.states)
size(b::TensorBasis) = (length(b.states), length(b.bases))
in(s::State, b::AbstractBasis)=in(s, b.states)
kind(b::AbstractBasis) = kind(b[1])

*{K<:BraKet}(a::AbstractBasis{K}, b::AbstractBasis{K}) = tensor(a,b)
+{K<:BraKet}(a::Basis{K},b::Basis{K}) = Basis([a.label, b.label], vcat(a.states,b.states))
-{B<:AbstractBasis}(a::B,b::B) = filter(x->!in(x,b), a)
setdiff{B<:AbstractBasis}(a::B,b::B) = a-b
get(b::AbstractBasis, label...) = b.label_map[[label...]]
get(b::AbstractBasis, label::Vector) = b.label_map[label]
get(b::AbstractBasis, s::State) = b.label_map[s.label]

#StateRepresentation###############################
type StateRep{K<:BraKet} <: AbstractStateRep{K}
	state::State{K}
	coeffs::Array{Complex{Float64}}
	basis::AbstractBasis{K}
	function StateRep(s::State{K}, coeffs::Array{Complex{Float64}}, basis::AbstractBasis{K})
		if length(basis)==length(coeffs)
			new(s, coeffs, basis)
		elseif length(basis)>length(coeffs)
			error("coefficients unspecified for $(length(basis)-length(coeffs)) basis states")
		else
			error("basis states unspecified for $(length(coeffs)-length(basis)) coefficients")
		end	
	end	
end

StateRep{N<:Number}(s::State, coeffs::Vector{N}, basis::AbstractBasis) = StateRep(s, convert(Array{Complex{Float64}},coeffs), basis)
StateRep{N<:Number}(s::State{Ket}, coeffs::Vector{N}, basis::AbstractBasis) = StateRep{Ket}(s, convert(Array{Complex{Float64}},coeffs), basis)
StateRep{N<:Number}(s::State{Bra}, coeffs::Vector{N}, basis::AbstractBasis) = error("Dimensions of coefficient array does not match type $K")
function StateRep{N<:Number, K<:BraKet}(s::State{K}, coeffs::Array{N}, basis::AbstractBasis)
	if size(coeffs)[2]==1 && K==Ket
		StateRep{Ket}(s, convert(Array{Complex{Float64}},vec(coeffs)), basis)
	elseif K==Bra
		StateRep{Bra}(s, convert(Array{Complex{Float64}},coeffs), basis)
	else
		error("Dimensions of coefficient array does not match type $K")
	end
end

StateRep{N<:Number}(label, coeffs::Vector{N}, basis::AbstractBasis) = StateRep{Ket}(State(label, Ket), convert(Vector{Number},coeffs), basis)
function StateRep{N<:Number}(label, coeffs::Array{N}, basis::AbstractBasis)
	if size(coeffs)[2]==1
		StateRep{Ket}(State(label, Ket), convert(Array{Complex{Float64}},vec(coeffs)), basis)
	else
		StateRep{Bra}(State(label, Bra), label, convert(Array{Complex{Float64}},coeffs), basis)
	end
end

kind(s::StateRep) = s.state.kind
label(s::StateRep) = s.state.label
repr(s::StateRep) = repr(s.state, " ; $(label(s.basis))")

copy(s::StateRep, coeffs=copy(s.coeffs)) = StateRep(s.state, coeffs, s.basis)
find(s::StateRep) = find(s.coeffs)
length(s::StateRep) = length(s.coeffs)
endof(s::StateRep) = length(s.coeffs)
getindex(s::StateRep, x) = s.coeffs[x]
setindex!(s::StateRep, y, x) = setindex!(s.coeffs, y, x)
get(s::StateRep, label) = s[get(s.basis, label)]

norm(s::StateRep) = norm(s.coeffs)

function normalize!(s::StateRep) 
	s.coeffs=(1/norm(s))*s.coeffs
	return s
end

normalize(s::StateRep) = normalize!(copy(s))

ctranspose(s::StateRep) = StateRep(s.state', s.coeffs', s.basis')

function map!(f::Function, s::StateRep)
	s.coeffs = map!(f, s.coeffs)
	return s
end 

map(f::Function, s::StateRep) = map!(f, copy(s))

function mapmatch!(f_coeffs::Function, f_states::Function, s::StateRep)
	matched_states = filter(f_states, s.basis)	
	for i in matched_states
		s[get(s.basis, i)] = apply(f_coeffs, get(s, i))
	end
	return s
end

mapmatch(f_coeffs::Function, f_states::Function, s::StateRep) = mapmatch!(f_coeffs, f_labels, copy(s))

filter(f::Function, s::StateRep) = mapmatch((x)->0, f, s)
filter!(f::Function, s::StateRep) = mapmatch!((x)->0, f, s)

function show(io::IO, s::StateRep)
	println("$(typeof(s)) $(repr(s)):")
	if length(s)!=0
		filled = find(s.coeffs)
		table = cell(length(filled), 2)	
		if length(filled)>=52
			for i=1:25
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= s.basis[filled[i]]
			end
			table[26:(length(filled)-25),:] = 0 # prevents access to undefined reference
			for i=(length(filled)-25):length(filled)
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= s.basis[filled[i]]
			end
		else
			for i=1:length(filled)
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= s.basis[filled[i]]
			end
		end
		temp_io = IOBuffer()
		if kind(s)==Ket
			show(temp_io, table)
		else
			show(temp_io, [table[:,2]', table[:,1]'])
		end
		io_str = takebuf_string(temp_io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		print(io_str)
	else
		println("(all coefficients are zero)")
	end
end

*(n::Number, s::StateRep) = copy(s, n*s.coeffs) 
*(s::StateRep, n::Number) = copy(s, s.coeffs*n) 
.*(n::Number, s::StateRep) = n*s
.*(s::StateRep, n::Number) = s*n 
.+(s::StateRep, n::Number) = copy(s, s.coeffs.+n)
.+(n::Number, s::StateRep) = copy(s, n.+s.coeffs)
.-(s::StateRep, n::Number) = copy(s, s.coeffs.-n)
.-(n::Number, s::StateRep) = copy(s, n.-s.coeffs)
/(s::StateRep, n::Number) = copy(s, s.coeffs/n)
./(s::StateRep, n::Number) = s/n
./(n::Number, s::StateRep) = copy(s, n./s.coeffs)
.^(n::Number, s::StateRep) = copy(s, n.^s.coeffs)
.^(s::StateRep, n::Number) = copy(s, s.coeffs.^n)

*(arr::Array, s::StateRep) = copy(s, arr*s.coeffs)
*(s::StateRep, arr::Array) = copy(s, s.coeffs*arr)

*(sr::StateRep{Bra}, s::State{Ket}) = get(sr, s')
*(sr::State{Bra}, s::StateRep{Ket}) = get(sr, s')

*(a::StateRep{Bra}, b::StateRep{Ket}) = (a.coeffs*b.coeffs)[1]
*(a::StateRep{Ket}, b::StateRep{Ket}) = StateRep(a.state*b.state, kron(a.coeffs, b.coeffs), a.basis*b.basis)
*(a::StateRep{Bra}, b::StateRep{Bra}) = StateRep(a.state*b.state, kron(a.coeffs, b.coeffs), a.basis*b.basis)
*(a::StateRep{Ket}, b::StateRep{Bra}) = Operator(a.coeffs*b.coeffs, a.basis, b.basis)
#Operator#####################################################
type Operator{} 
	coeffs::Matrix{Complex{Float64}}
	row_basis::AbstractBasis{Ket}
	col_basis::AbstractBasis{Bra}
	function Operator(coeffs::Matrix{Complex{Float64}}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra})
		if size(coeffs)==(length(row_basis), length(col_basis))
			new(coeffs, row_basis, col_basis)
		else
			throw(DimensionMismatch)
		end
	end
end

Operator{N<:Number}(coeffs::Matrix{N}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), row_basis, col_basis) 
Operator{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Ket}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), b, b') 
Operator{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Bra}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), b', b) 

function Operator(label_func::Function, coeff_func::Function, b::AbstractBasis)
	coeffs = Array(Number, length(b), length(b))
	for i=1:length(b)
		for j=1:length(b)
			coeffs[i,j] = coeff_func(b[i]) * (b[j]'*label_func(b[i]))
		end
	end
	return Operator(coeffs, b)
end

samebasis(a::Operator, b::Operator) = isequal(a.row_basis,b.row_basis) && isequal(a.col_basis, b.col_basis)
isequal(a::Operator, b::Operator) = a.coeffs==b.coeffs && samebasis(a,b)

copy(op::Operator, coeffs=copy(op.coeffs)) = Operator(coeffs, op.row_basis, op.col_basis)
expm(op::Operator) = copy(op, expm(op.coeffs))

ndims(op::Operator) = 2
size(op::Operator) = size(op.coeffs)
length(op::Operator) = length(op.coeffs)
endof(op::Operator) = length(op)
find(op::Operator) = find(op.coeffs)
find(f::Function, op::Operator) = find(f::Function, op.coeffs)

function findstates(f::Function, op::Operator) #doesn't work perfectly just saying
	lens = size(op)
	inds = find(f, op)
	result = Array(State, length(inds), 2)
	for i=1:length(inds)
		mod_row = inds[i]<lens[1] ? 0 : 1
		mod_col = inds[i]<lens[2] ? 0 : 1
		result[i, 1] = op.row_basis[inds[i]%lens[1]+mod_row]
		result[i, 2] = op.col_basis[inds[i]%lens[2]+mod_col]
	end
	return result
end

size(op::Operator, i::Int) = size(op.coeffs, i)
ctranspose(op::Operator) = Operator(op.coeffs', op.col_basis',op.row_basis')
getindex(op::Operator, x...) = op.coeffs[x...]
setindex!(op::Operator, y, x) = setindex!(op.coeffs,y,x)
get(op::Operator, ket_label, bra_label) = op[get(op.row_basis, ket_label), get(op.col_basis, bra_label)]
get(op::Operator, s::State{Ket}) = op[get(op.row_basis, s), :]
get(op::Operator, s::State{Bra}) = op[:, get(op.col_basis, s)]

*(op::Operator, n::Number) = copy(op, op.coeffs*n)
*(n::Number, op::Operator) = copy(op, n*op.coeffs)
/(op::Operator, n::Number) = copy(op, op.coeffs/n)
.+(op::Operator, n::Number) = copy(op, op.coeffs.+n)
.+(n::Number, op::Operator) = copy(op, n.+op.coeffs)
.-(op::Operator, n::Number) = copy(op, op.coeffs.-n)
.-(n::Number, op::Operator) = copy(op, n.-op.coeffs)
.^(op::Operator, n::Number) = copy(op, op.coeffs.^n)
.^(n::Number, op::Operator) = copy(op, n.^op.coeffs)

*(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs*b.coeffs) : error("BasesMismatch")
+(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs+b.coeffs) : error("BasesMismatch")
-(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs-b.coeffs) : error("BasesMismatch")

commutator(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, (a.coeffs*b.coeffs)-(b.coeffs*a.coeffs)) : error("BasesMismatch")

*(op::Operator, s::State{Ket}) = StateRep(State([]), get(op, s'), op.row_basis)
*(s::State{Bra}, op::Operator) = StateRep(State([], Bra), get(op, s'), op.col_basis)
*(op::Operator, s::StateRep{Ket}) = op.row_basis == s.basis ? StateRep(s.state, op.coeffs*s.coeffs, op.row_basis) : error("BasesMismatch")
*(s::StateRep{Bra}, op::Operator) = op.col_basis == s.basis ? StateRep(s.state, s.coeffs*op.coeffs, op.col_basis) : error("BasesMismatch")
*(arr::Array, op::Operator) = copy(op, arr*op.coeffs)
*(op::Operator, arr::Array) = copy(op, op.coeffs*arr)

trace(op::Operator) = trace(op.coeffs)

function ptrace(op::Operator, ind::Int) #doesn't produce a normalized matrix
	if isequal(op.col_basis, op.row_basis')
		tr_row = tensor(vcat(separate(op.row_basis)[1:ind-1], separate(op.row_basis)[ind+1:end])...)
		tr_col = tr_row'
		coeffs = Complex{Float64}[sum([op[i,j]*(op[i,k])' for i=1:length(separate(op.row_basis)[ind])]) for j=1:length(tr_row), k=1:length(tr_col)] 
	else
		error("BasesMismatch")
	end
	return Operator(coeffs, tr_row, tr_col) 
end

function show(io::IO, op::Operator)
	println("$(typeof(op)):")
	table = cell(length(op.row_basis)+1, length(op.col_basis)+1)	
	for i = 1:length(op.row_basis)
		table[i+1,1] = op.row_basis[i]
	end
	for j = 1:length(op.col_basis)
		table[1,j+1] = op.col_basis[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = op.coeffs	
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	io_str = io_str[searchindex(io_str, "\n")+3:end]
	print(io_str)
end

#module ends#######################################
end

using d
s = d.State("sr")
a = d.Basis("a", [1:10])
b = d.Basis("b", ["$i" for i=1:4]);
c = d.Basis("c", [1:20]);
r = d.tensor(a,b,c)
sr = d.normalize!(d.StateRep(s, [1 for i=1:length(r)], r))
op = (sr*sr')
print("")
