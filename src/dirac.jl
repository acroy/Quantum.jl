module d
include("rep.jl")
#includes,imports,consts#############################
const lang = "\u27E8"
const rang = "\u27E9"
const vert_ell = "\u22EE"
const horiz_ell = "\u2026"
import Base.show,
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

State(label::Vector) = State(label, Ket)
State{T}(label::T) = State(T[label], Ket)
State(label...) = State([label...], Ket)

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

ctranspose(s::State) = State(s.label, !s.kind)
getindex(s::State, x) = s.label[x]
endof(s::State) = endof(s.label)
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

function components(B::TensorBasis)
	return B.bases
end

function components(B::Basis)
	return B
end



#module ends#######################################
end

using d

a = d.Basis("a", [1:10])
b = d.Basis("b", ["$i" for i=1:4])



