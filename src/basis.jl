#####################################
#Basis###############################
#####################################

immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label::String
	states::Vector{State{K}}
	statemap::Dict{State{K}, Int}
end

makebasis{K<:BraKet}(label::String, states::Vector{State{K}}) = Basis(label, unique(states), (State{K}=>Int)[states[i]=>i for i=1:length(states)])

Basis{K<:BraKet}(label::String, labelvec::Vector, kind::Type{K}=Ket) = makebasis(label, statearr(labelvec, label, kind))
Basis{K<:BraKet}(labelvec::Vector, kind::Type{K}=Ket) = Basis("?", labelvec, kind)														
Basis{K<:BraKet}(s::State{K}...) = Basis(vcat(s...))
function Basis{K<:BraKet}(s::Vector{State{K}}) 
	bases = unique(map(basislabel, s))
	if length(bases)>1
		makebasis("?", s)
	else
		makebasis(bases[1], s)
	end
end

#####################################
#TensorBasis#########################
#####################################
immutable TensorBasis{K<:BraKet} <: AbstractBasis{K}
	bases::Vector{Basis{K}}
	states::Vector{TensorState{K}}
	statemap::Dict{TensorState{K}, Int}
end

TensorBasis{K<:BraKet}(bases::Vector{Basis{K}}, states::Vector{TensorState{K}}) = TensorBasis(bases, unique(states), (TensorState{K}=>Int)[states[i]=>i for i=1:length(states)])

function tensor{K<:BraKet}(bases::AbstractBasis{K}...)
	TensorBasis(vcat([separate(i) for i in bases]...), convert(Vector{TensorState{K}}, statejoin(tensorarr([i.states for i in bases]...))))
end

function tensor{K<:BraKet}(basis::AbstractBasis{K})
	return basis
end

#####################################
#Functions###########################
#####################################

#utility#############################

label(b::Basis) = "$(b.label)"
function label(b::TensorBasis)
	labels = [label(i) for i in b.bases]
	#terrible way to grow a string
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str$otimes$(labels[i])"
	end
	return str
end

#imported############################
isequal(a::AbstractBasis,b::AbstractBasis) = a.states==b.states
in(s::State, b::AbstractBasis)=in(s, b.states)

filter(f::Function, b::Basis) = makebasis(b.label, filter(f, b.states))
function filter(f::Function, b::TensorBasis) 
	states = filter(f, b.states)
	b_arr = Array(Basis{kind(b)}, size(b,2))
	for i=1:size(b)[2]
		b_arr[i] = Basis("?", unique(map(x->getindex(x,i), states)))
    end
	TensorBasis(b_arr, states)
end

ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states))
ctranspose(b::Basis) = Basis(map(ctranspose, b.states))

size(b::TensorBasis) = (length(b.states), length(b.bases))
length(b::AbstractBasis) = length(b.states)
endof(b::AbstractBasis) = endof(b.states)

*{K<:BraKet}(a::AbstractBasis{K}, b::AbstractBasis{K}) = tensor(a,b)
+{K<:BraKet}(a::Basis{K}, b::Basis{K}) = Basis(vcat(a.states,b.states))
setdiff{B<:AbstractBasis}(a::B,b::B) = setdiff(a.states, b.states)

getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

get(b::AbstractBasis, s::AbstractState, notfound) = get(b.statemap, s, notfound)
get(b::Basis, s::State) = b.statemap[s]
get(b::TensorBasis, s::TensorState) = b.statemap[s]
get(b::Basis, label) = get(b, State(label))
get(b::TensorBasis, label) = get(b, TensorState(label))

showcompact(io::IO, b::AbstractBasis) = print(io, "$(typeof(b)) $(label(b))")

function show(io::IO, b::AbstractBasis)
	showcompact(io, b)
	println(", $(length(b)) states:")
	if length(b)>20
		for i=1:10
			println(io, b.states[i])
		end
		println(vdots)
		for i=length(b)-10:length(b)
			println(io, b.states[i])
		end
	else
		for i in b.states
			println(io, i)
		end	
	end
end

# #exported############################
statetobasis(s::State) = Basis(s)
statetobasis(s::TensorState) = TensorBasis(map(Basis, separate(s)), [s])

separate(b::Basis)=b
separate(b::TensorBasis) = b.bases
kind(b::AbstractBasis) = kind(b[1])
