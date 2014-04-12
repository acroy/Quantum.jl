module d
include("rep.jl")
#includes,imports,consts#############################
const lang = "\u27E8"
const rang = "\u27E9"
const vert_ell = "\u22EE"
const horiz_ell = "\u2026"
const otimes = "\u2297"
import Base.show,
	   Base.repr,	
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
	   Base.start  
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

#Basisless States###################################
abstract AbstractState
abstract BraKet
abstract Bra <: BraKet
abstract Ket <: BraKet
!(K::Type{Ket}) = Bra
!(B::Type{Bra}) = Ket

immutable State{K<:BraKet}
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

ctranspose(s::State) = State(s.label, !s.kind)
getindex(s::State, x) = s.label[x]
setindex!(s::State, y, x) = setindex!(s.label, y, x)
endof(s::State) = endof(s.label)
repr(s::State{Bra}, bra_sym=lang) = "$bra_sym $(repr(s.label)[2:end-1]) |"
repr(s::State{Ket}, ket_sym=rang) = "| $(repr(s.label)[2:end-1]) $ket_sym"
label(s::State) = s.label

function show(io::IO, s::State)
	print(io, repr(s))
end

#Bases#############################################
abstract AbstractBasis

function mapstates{K<:BraKet}(svec::Vector{State{K}})
	dict = Dict{State{K}, Int}()
	sizehint(dict, length(svec))
	for i=1:length(svec)
		dict[svec[i]] = i
	end
	return dict
end	

immutable Basis{K<:BraKet} <: AbstractBasis
	label
	states::Vector{State{K}}
	state_map::Dict{State{K}, Int}
	bra_sym::String
	ket_sym::String
end

Basis{K<:BraKet}(label, states::Vector{State{K}}; bra_sym=lang, ket_sym=rang) = Basis(label, unique(states), mapstates(states), bra_sym, ket_sym)
Basis{K<:BraKet}(label, states::State{K}...; bra_sym=lang, ket_sym=rang) = Basis(label, vcat(states...), bra_sym=bra_sym, ket_sym=ket_sym)
Basis(label, label_vec::Vector; bra_sym=lang, ket_sym=rang) = Basis(label, statevec(label_vec), bra_sym=bra_sym, ket_sym=ket_sym)														

immutable TensorBasis{K<:BraKet} <: AbstractBasis
	bases::Vector{Basis{K}}
	states::Vector{State{K}}
	state_map::Dict{State{K}, Int}
	bra_sym::String
	ket_sym::String
end

function TensorBasis{K<:BraKet}(bases::Vector{Basis{K}}, states::Vector{State{K}}; bra_sym=lang, ket_sym=rang)
	states = unique(states)
	state_map = mapstates(states)
	TensorBasis(bases, states, state_map, bra_sym, ket_sym)
end

function tensor(bases::AbstractBasis...; bra_sym=lang, ket_sym=rang)
	states = tensor([i.states for i in bases]...)
	TensorBasis(vcat([separate(i) for i in bases]...), states, bra_sym=bra_sym, ket_sym=ket_sym)
end

separate(b::TensorBasis) = b.bases
separate(b::Basis)=b
kind(b::AbstractBasis)=b.states[1].kind
isequal(a::AbstractBasis,b::AbstractBasis) = a.states==b.states

label(b::Basis) = repr(b.label)

function label(b::TensorBasis)
	labels = [label(i) for i in b.bases]
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str $otimes $(labels[i])"
	end
	str
end

function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) $(label(b)):")
	for i in b.states
		println(repr(i))
	end
end

ctranspose(b::Basis) = Basis(b.label, map(ctranspose, b.states))
ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states))

states(b::AbstractBasis) = b.states
filter(f::Function, b::Basis) = Basis(b.label, filter(f, b.states), bra_sym=b.bra_sym, ket_sym=b.ket_sym)
filter(f::Function, b::TensorBasis) = TensorBasis(b.bases, filter(f, b.states), bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)
endof(b::AbstractBasis) = endof(b.states)
length(b::AbstractBasis) = length(b.states)
size(b::TensorBasis) = (length(b.states), length(b.bases))
in(s::State, b::AbstractBasis)=in(s, b.states)

*(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
+(a::Basis,b::Basis) = Basis([a.label, b.label], vcat(a.states,b.states), bra_sym=a.bra_sym, ket_sym=a.ket_sym)
-{B<:AbstractBasis}(a::B,b::B) = filter(x->!in(x,b), a)
setdiff{B<:AbstractBasis}(a::B,b::B) = a-b
get(b::AbstractBasis, label...) = b.state_map[State(label...)]
get(b::AbstractBasis, label::Vector) = b.state_map[State(label)]
get(b::AbstractBasis, state_key::State) = b.state_map[state_key]

#StateRepresentation###############################
type StateRep{K<:BraKet} <: AbstractState
	state::State{K}
	coeffs::Array{Complex{Float64}}
	basis::AbstractBasis
	function StateRep(s::State{K}, coeffs::Array{Complex{Float64}}, basis::AbstractBasis)
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

repr(s::StateRep) = repr(s.state)
copy(s::StateRep) = StateRep(s.state, copy(s.coeffs), s.basis)
copy{N<:Number}(s::StateRep, coeffs::Array{N}) = StateRep(s.state, coeffs, s.basis)

length(s::StateRep) = length(find(s.coeffs))
getindex(s::StateRep, x) = s.coeffs[x]
setindex!(s::StateRep, y, x) = setindex!(s.coeffs, y, x)
get(s::StateRep, label...) = s[get(s.basis, label...)]
get(s::StateRep, label::Vector) = s[get(s.basis, label)]
get(s::StateRep, state_key::State) = s[get(s.basis, state_key)]

magnitude(A::Number...)=sqrt(sum(A.^2))
magnitude{N<:Number}(A::Array{N}) = sqrt(sum(A.^2))
magnitude(s::StateRep) = magnitude(s.coeffs)

function normalize!(s::StateRep) 
	s.coeffs=(1/magnitude(s))*s.coeffs
	return s
end

normalize(s::StateRep) = normalize!(copy(s))
isnorm(s::StateRep) = magnitude(s)==1

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
		kind(s)==Ket ? show(temp_io, table) : show(temp_io, [table[:,2].', table[:,1].'])
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

*(a::StateRep{Bra}, b::StateRep{Ket}) = (a.coeffs*b.coeffs)[1]
*(a::StateRep{Ket}, b::StateRep{Ket}) = StateRep(a.state*b.state, mapslices(prod, crossjoin(a.coeffs, b.coeffs), 2), a.basis*b.basis)
*(a::StateRep{Bra}, b::StateRep{Bra}) = StateRep(a.state*b.state, mapslices(prod, crossjoin(a.coeffs', b.coeffs'), 2)', a.basis*b.basis)
*(a::StateRep{Ket}, b::StateRep{Bra}) = Operator("$(repr(a))$(repr(b))", a.coeffs*b.coeffs, a.basis, b.basis)
#Operator#####################################################
type Operator 
	label
	coeffs::Matrix{Complex{Float64}}
	row_basis::AbstractBasis
	col_basis::AbstractBasis
	function Operator(label, coeffs::Matrix{Complex{Float64}}, row_basis::AbstractBasis, col_basis::AbstractBasis)
		if kind(row_basis)==Ket && kind(col_basis)==Bra
			new(label, coeffs, row_basis, col_basis)
		else
			error("input bases have wrong orientation")
		end
	end
end

Operator{N<:Number}(label, coeffs::Matrix{N}, row_basis::AbstractBasis, col_basis::AbstractBasis) = Operator(label, convert(Matrix{Complex{Float64}}, coeffs), row_basis, col_basis) 
Operator{N<:Number}(label, coeffs::Matrix{N}, b::AbstractBasis) = Operator(label, convert(Matrix{Complex{Float64}}, coeffs), b, b') 

function Operator(label, label_func::Function, coeff_func::Function, b::AbstractBasis)
	coeffs = Array(Number, length(b), length(b))
	if kind(b)==Bra
		b = b'
	end
	for i=1:length(b)
		for j=1:length(b)
			coeffs[i,j] = coeff_func(b[i]) * (b[j]'*label_func(b[i]))
		end
	end
	return Operator(label, coeffs, b)
end

isequal(a::Operator, b::Operator) = isequal(a.label, b.label) && coeffs==coeffs && isequal(a.row_basis,b.row_basis) && isequal(a.col_basis, b.col_basis)

copy(op::Operator) = Operator(copy(op.label), copy(op.coeffs), op.row_basis, op.col_basis)
copy{N<:Number}(op::Operator, coeffs::Matrix{N}) = Operator(copy(op.label), coeffs, op.row_basis, op.col_basis)
expm(op::Operator) = copy(op, exp(op.coeffs))

*(op::Operator, n::Number) = copy(op, op.coeffs*n)
*(n::Number, op::Operator) = copy(op, n*op.coeffs)
/(op::Operator, n::Number) = copy(op, op.coeffs/n)
.+(op::Operator, n::Number) = copy(op, op.coeffs.+n)
.+(n::Number, op::Operator) = copy(op, n.+op.coeffs)
.-(op::Operator, n::Number) = copy(op, op.coeffs.-n)
.-(n::Number, op::Operator) = copy(op, n.-op.coeffs)
.^(op::Operator, n::Number) = copy(op, op.coeffs.^n)
.^(n::Number, op::Operator) = copy(op, n.^op.coeffs)

*(op::Operator, s::StateRep{Ket}) = op.row_basis == s.basis ? StateRep(s.state, op.coeffs*s.coeffs, op.row_basis) : error("Bases don't match")
*(s::StateRep{Bra}, op::Operator) = op.col_basis == s.basis ? StateRep(s.state, s.coeffs*op.coeffs, op.col_basis) : error("Bases don't match")

function show(io::IO, op::Operator)
	println("$(typeof(op)) $(op.label):")
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

s = d.State([1:3])
a = d.Basis("a", [1:10])
b = d.Basis("b", ["$i" for i=1:4]);
c = d.tensor(a,b)
sr = d.StateRep(s, [1:10], a)

raiselabel(s::d.State)=d.State(s.label[1]+1, s.kind)
raisecoeff(s::d.State)=sqrt(s[1]+1)
op = d.Operator("dag", raiselabel, raisecoeff, a)
print("")
