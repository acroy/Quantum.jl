#####################################
#Basis###############################
#####################################

immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label
	states::Vector{State{K}}
	label_map::Dict{Vector, Int}
end

Basis{K<:BraKet}(label, states::Vector{State{K}}) = Basis(label, unique(states), maplabels(states))
Basis{K<:BraKet}(label, states::State{K}...) = Basis(label, vcat(states...))
Basis(label, label_vec::Vector) = Basis(label, statevec(label_vec))														

#####################################
#TensorBasis#########################
#####################################
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

function tensor{K<:BraKet}(basis::AbstractBasis{K})
	return basis
end

#####################################
#Functions###########################
#####################################

#utility#############################

function maplabels{K<:BraKet}(svec::Vector{State{K}})
	dict = Dict{Vector, Int}()
	sizehint(dict, length(svec))
	for i=1:length(svec)
		dict[svec[i].label] = i
	end
	return dict
end	


label(b::Basis) = b.label
function label(b::TensorBasis)
	labels = [label(i) for i in b.bases]
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str $otimes $(labels[i])"
	end
	return str
end

#imported############################
isequal(a::AbstractBasis,b::AbstractBasis) = a.states==b.states
in(s::State, b::AbstractBasis)=in(s, b.states)

filter(f::Function, b::Basis) = Basis(b.label, filter(f, b.states))
function filter(f::Function, b::TensorBasis) 
	states = filter(f, b.states)
	b_arr = Array(Basis{kind(b)}, size(b)[2])
	for i=1:size(b)[2]
		b_arr[i] = Basis("sub_($(label(b)))_$i", unique(map(x->getindex(x,i), states)))
    end
	TensorBasis(b_arr, states)
end

ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states), b.label_map)
ctranspose(b::Basis) = Basis(b.label, map(ctranspose, b.states), b.label_map)

size(b::TensorBasis) = (length(b.states), length(b.bases))
length(b::AbstractBasis) = length(b.states)
endof(b::AbstractBasis) = endof(b.states)

*{K<:BraKet}(a::AbstractBasis{K}, b::AbstractBasis{K}) = tensor(a,b)
+{K<:BraKet}(a::Basis{K},b::Basis{K}) = Basis([a.label, b.label], vcat(a.states,b.states))
-{B<:AbstractBasis}(a::B,b::B) = filter(x->!in(x,b), a)
setdiff{B<:AbstractBasis}(a::B,b::B) = a-b

getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

get(b::AbstractBasis, label::Vector, notfound) = get(b.label_map, label, notfound)
get(b::AbstractBasis, s::State, notfound) = get(b.label_map, s.label, notfound)

function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) $(label(b)):")
	for i in b.states
		println(repr(i))
	end
end

#exported############################
separate(b::Basis)=b
separate(b::TensorBasis) = b.bases

states(b::AbstractBasis) = b.states
kind(b::AbstractBasis) = kind(b[1])
