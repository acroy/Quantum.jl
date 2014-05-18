#####################################
#Basis###############################
#####################################

immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label::Symbol
	states::Vector{State{K}}
	label_map::Dict{Any, Int}
end

makebasis{K<:BraKet}(label::Symbol, states::Vector{State{K}}) = Basis(label, unique(states), maplabels(states)) #hidden function; not exported

Basis{K<:BraKet}(label::Symbol, labelvec::Vector, kind::Type{K}=Ket) = makebasis(label, statearr(labelvec, label, kind))
Basis{K<:BraKet}(labelvec::Vector, kind::Type{K}=Ket) = Basis(:?, label_vec, kind)														
Basis{K<:BraKet}(svec::Vector{State{K}}) = Basis(svec...) #this seems inefficient
function Basis{K<:BraKet}(s::State{K}...) 
	bases = unique(map(basis, s))
	if length(bases)>1
		makebasis(:?, s)
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
	label_map::Dict{Vector, Int}
end

TensorBasis{K<:BraKet}(bases::Vector{Basis{K}}, states::Vector{TensorState{K}}) = TensorBasis(bases, unique(states), maplabels(states))

function tensor{K<:BraKet}(bases::AbstractBasis{K}...)
	TensorBasis(vcat([separate(i) for i in bases]...), convert(Array{TensorState{K}}, statejoin(tensorarr([i.states for i in bases]...))))
end

function tensor{K<:BraKet}(basis::AbstractBasis{K})
	return basis
end

#####################################
#Functions###########################
#####################################

#utility#############################

function maplabels{K<:BraKet}(svec::Vector{TensorState{K}})
	dict = Dict{Vector, Int}()
	for i=1:length(svec)
		dict[label(svec[i])]=i
	end
	return dict
end

maplabels{K<:BraKet}(svec::Vector{State{K}}) = [label(svec[i])=>i for i=1:length(svec)]

label(b::Basis) = "$b.label"
function label(b::TensorBasis)
	labels = [label(i) for i in b.bases]
	#terrible way to grow a string
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
get(b::AbstractBasis, key) = get(b, key, nothing)

showcompact(io::IO, b::AbstractBasis) = print(io, "$(typeof(b)) $(label(b))")

function show(io::IO, b::AbstractBasis)
	showcompact(io, b)
	println(":")
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
statetobasis(s::State) = length(s)>1 ? TensorBasis(map(Basis, separate(s)), [s]) : Basis(s)
separate(b::Basis)=b
separate(b::TensorBasis) = b.bases
samelabels(b1::AbstractBasis, b2::AbstractBasis) = collect(keys(b1.label_map))==collect(keys(b2.label_map))
kind(b::AbstractBasis) = kind(b[1])
