module d
include("rep.jl")
#includes,imports,consts#############################
const lang = "\u27E8"
const rang = "\u27E9"
const vert_ell = "\u22EE"
const horiz_ell = "\u2026"
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
	   Base.endof  
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

immutable State{K<:BraKet,T}
  label::Vector{T}
  kind::Type{K}
end

State(label::Vector, kind=Ket) = State(label, kind)
State{T}(label::T, kind=Ket) = State(T[label], kind)
State(label...; kind=Ket) = State([label...], kind)

function statevec(v::Vector)
	svec = Array(State, length(v))
	for i=1:length(v)
		svec[i] = State(v[i])
	end
	return svec
end

function statevec(arr::Array)
	svec = Array(State, size(arr,1))
	for i=1:size(arr, 1)
		svec[i] = State(vec(arr[i,:]))
	end
	return svec
end

tensor() = nothing
tensor{K<:BraKet}(s::State{K}...) = State(vcat([i.label for i in s]...), K)
statejoin(v::Vector{State}...) = broadcast(tensor, v...)
function statejoin(state_arr::Array{State}) 
	result = statejoin(state_arr[:,1], state_arr[:,2])
	for i=3:size(state_arr, 2)
		result = statejoin(result, state_arr[:,i])
	end
	return result
end
tensor(state_arrs::Array{State}...) = statejoin(crossjoin(state_arrs...))
separate(s::State) = statevec(s.label)
separate(v::Vector{State}) = hcat(map(separate, v)...).'

ctranspose(s::State) = State(s.label, !s.kind)
getindex(s::State, x) = s.label[x]
setindex!(s::State, y, x) = setindex!(s.label, y, x)
endof(s::State) = endof(s.label)
repr(s::State{Bra}, bra_sym=lang) = "$bra_sym $(repr(s.label)[2:end-1]) |"
repr(s::State{Ket}, ket_sym=rang) = "| $(repr(s.label)[2:end-1]) $ket_sym"

function show(io::IO, s::State)
	print(io, repr(s))
end

#Bases#############################################
abstract AbstractBasis

function mapstates{S<:State}(svec::Vector{S})
	dict = Dict{S, Int}()
	sizehint(dict, length(svec))
	for i=1:length(svec)
		dict[svec[i]] = i
	end
	return dict
end	

immutable Basis{S<:State} <: AbstractBasis
	label
	states::Vector{S}
	state_map::Dict{S, Int}
	bra_sym::String
	ket_sym::String
end

function Basis{S<:State}(label, states::Vector{S}; bra_sym=lang, ket_sym=rang)
	state_map = mapstates(states)
	Basis(label, states, state_map, bra_sym, ket_sym)
end

function Basis(label, label_vec::Vector; bra_sym=lang, ket_sym=rang)
	states = statevec(unique(label_vec))
	Basis(label, states, bra_sym=bra_sym, ket_sym=ket_sym)
end

Basis{S<:State}(label, states::S...; bra_sym=lang, ket_sym=rang) = Basis(label, vcat(states...), bra_sym=bra_sym, ket_sym=ket_sym)

immutable TensorBasis{B<:Basis} <: AbstractBasis
	label
	bases::Vector{B}
	states::Vector{State}
	state_map::Dict{State, Int}
	bra_sym::String
	ket_sym::String
end

function TensorBasis{B<:Basis, S<:State}(label, bases::Vector{B}, states::Vector{S}; bra_sym=lang, ket_sym=rang)
	state_map = mapstates(states)
	TensorBasis(label, bases, states, state_map, bra_sym, ket_sym)
end

function TensorBasis{B<:Basis}(label, bases::Vector{B}; bra_sym=lang, ket_sym=rang)
	states = tensor([i.states for i in bases]...)
	TensorBasis(label, bases, states, bra_sym=bra_sym, ket_sym=ket_sym)
end

function TensorBasis{B<:Basis}(bases::Vector{B}; bra_sym=lang, ket_sym=rang)
	label = [i.label for i in bases]
	TensorBasis(label, bases, bra_sym=bra_sym, ket_sym=ket_sym)
end

function tensor(bases::AbstractBasis...)
	TensorBasis(vcat([components(i) for i in bases]...))
end

function components(b::TensorBasis)
	return b.bases
end

function components(b::Basis)
	return b
end

function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) $(repr(b.label)):")
	for i in b.states
		println(repr(i))
	end
end

filter(f::Function, b::Basis) = Basis(b.label, filter(f, b.states), bra_sym=b.bra_sym, ket_sym=b.ket_sym)
filter(f::Function, b::TensorBasis) = TensorBasis(b.label, b.bases, filter(f, b.states), bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)
endof(b::AbstractBasis) = endof(b.states)
length(b::AbstractBasis) = length(b.states)
size(b::TensorBasis) = (length(b.states), length(b.bases))

in(s::State, b::AbstractBasis)=in(s, b.states)

*(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
+(a::Basis,b::Basis) = Basis("$(repr(a.label))+$(repr(b.label))", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym)
+(a::TensorBasis,b::TensorBasis) = size(a,2)==size(b,2) ? TensorBasis("$(repr(a.label))+$(repr(b.label))", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym) : error("dimension mismatch; label lengths differ")
-{B<:AbstractBasis}(a::B,b::B) = filter(x->!in(x,b), a, name="$(repr(a.label))-$(repr(b.label))")
setdiff(a::AbstractBasis,b::AbstractBasis) = a-b
get(b::AbstractBasis, label...) = b.state_map[State(label...)]
get(b::AbstractBasis, label::Vector) = b.state_map[State(label)]
get(b::AbstractBasis, state_key::State) = b.state_map[state_key]



#StateRepresentation###############################
type StateRep{B<:AbstractBasis, K<:BraKet, N<:Number} <: AbstractState
	state::State{K}
	coeffs::Array{N}
	basis::B
	function StateRep(s::State{K}, coeffs::Array{N}, basis::B)
		if length(basis)==length(coeffs)
			new(s, coeffs, basis)
		elseif length(basis)>length(coeffs)
			error("coefficients unspecified for $(length(basis)-length(coeffs)) basis states")
		else
			error("basis states unspecified for $(length(coeffs)-length(basis)) coefficients")
		end	
	end	
end


StateRep{B<:AbstractBasis,N<:Number}(s::State{Ket}, coeffs::Vector{N}, basis::B) = StateRep{B, Ket, N}(s, coeffs, basis)
function StateRep{B<:AbstractBasis,N<:Number, K<:BraKet}(s::State{K}, coeffs::Array{N}, basis::B)
	if size(coeffs)[2]==1 && K==Ket
		StateRep{B, Ket, N}(s, vec(coeffs), basis)
	elseif K==Bra
		StateRep{B, Bra, N}(s, coeffs, basis)
	else
		error("Dimensions of coefficient array does not match type $K")
	end
end

function StateRep{B<:AbstractBasis,N<:Number}(label, coeffs::Array{N}, basis::B)
	if size(coeffs)[2]==1
		StateRep{B, Ket, N}(State(label, Ket), vec(coeffs), basis)
	else
		StateRep{B, Bra, N}(State(label, Bra), label, coeffs, basis)
	end
end

kind(s::StateRep) = s.state.kind
label(s::StateRep) = s.state.label

repr(s::StateRep) = repr(s.state)
copy(s::StateRep) = StateRep(s.state, copy(s.coeffs), s.basis)
copy{N<:Number}(S::StateRep, coeffs::Array{N}) = StateRep(s.state, coeffs, s.basis)

length(s::StateRep) = length(find(s.coeffs))
getindex(s::StateRep, x) = s.coeffs[x]
setindex!(s::StateRep, y, x) = setindex!(s.coeffs, y, x)
get(s::StateRep, label...) = s[get(s.basis, label...)]
get(s::StateRep, label::Vector) = s[get(s.basis, label)]
get(s::StateRep, state_key::State) = s[get(s.basis, state_key)]

function magnitude(A::Number...)
	if length(A)==2
		return hypot(A[1], A[2])
	end
	return magnitude(hypot(A[1], A[2]), A[3:end]...)
end

magnitude{N<:Number}(A::Array{N}) = magnitude(A...)
magnitude(s::StateRep) = magnitude(s.coeffs)

function normalize!(s::StateRep) 
	s.coeffs=(1/magnitude(s))*s.coeffs
	return s
end

normalize(s::StateRep) = normalize!(copy(s))
isnorm(s::StateRep) = magnitude(s)==1

ctranspose(s::StateRep) = State(s.state', s.coeffs', s.basis)

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
		show(temp_io, table)
		io_str = takebuf_string(temp_io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		print(io_str)
	else
		println("(all coefficients are zero)")
	end
end


#module ends#######################################
end

using d

s = d.State([1:3])
a = d.Basis("a", [1:10])
b = d.Basis("b", ["$i" for i=1:4]);
print("")



