#####################################
#Basis###############################
#####################################

immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label::String
	states::Vector{State{K}}
	statemap::Dict{(Any,String), Int}
end

function makebasis{K}(label::String, states::Vector{State{K}})
	states = unique(states)
	return Basis(label, states, ((Any,String)=>Int)[(states[i].label,states[i].basislabel)=>i for i=1:length(states)])
end

Basis{K<:BraKet}(labelvec::Vector, label::String, kind::Type{K}=Ket) = makebasis(label, statearr(labelvec, label, kind))
Basis{K}(s::State{K}...) = Basis(vcat(s...))
function Basis{K<:BraKet}(s::Vector{State{K}}) 
	bases = unique(map(basislabel, s))
	if length(bases)>1
		makebasis("?", s)
	else
		makebasis(bases[1], s)
	end
end

copy(b::Basis) = Basis(copy(b.label), copy(b.states), copy(b.statemap))

#####################################
#TensorBasis#########################
#####################################
immutable TensorBasis{K<:BraKet} <: AbstractBasis{K}
	bases::Vector{Basis{K}}
	states::Vector{TensorState{K}}
	statemap::Dict{(Vector, Vector{String}), Int}
end

function TensorBasis{K}(bases::Vector{Basis{K}}, states::Vector{TensorState{K}})
	return TensorBasis(bases, unique(states), ((Vector, Vector{String})=>Int)[(label(states[i]), basislabel(states[i]))=>i for i=1:length(states)])
end

function TensorBasis{K}(states::Vector{TensorState{K}})
	sepstates = hcat(map(separate, states)...).'
	bases = Array(Basis{K}, size(sepstates, 2))
	for i=1:size(sepstates, 2)
		bases[i] = Basis(sepstates[:, i])
	end
	return TensorBasis(bases, states)
end

copy(b::TensorBasis) = TensorBasis(copy(b.bases), copy(b.states), copy(b.statemap))

function tensor{K}(bases::AbstractBasis{K}...)
	TensorBasis(vcat([separate(i) for i in bases]...), convert(Vector{TensorState{K}}, statejoin(tensorarr([i.states for i in bases]...))))
end

function tensor{K}(basis::AbstractBasis{K})
	return basis
end

#####################################
#Functions###########################
#####################################

#utility#############################

label(b::Basis) = b.label
label(b::TensorBasis) = map(label, b.bases)
basislabel(b::AbstractBasis) = label(b)
reprlabel(b::Basis) = label(b)
function reprlabel(b::TensorBasis)
	labels = label(b)
	#terrible way to grow a string
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str$otimes$(labels[i])"
	end
	return str
end

#imported############################
isequal(a::AbstractBasis, b::AbstractBasis) = isequal(a.states, b.states) && label(a)==label(b)
==(a::AbstractBasis, b::AbstractBasis) = a.states==b.states && label(a)==label(b)
in(s::AbstractState, b::AbstractBasis)=in(s, b.states)
find(f::Function, b::Basis) = find(f, b.states)
filter(f::Function, b::Basis) = makebasis(b.label, filter(f, b.states))
function filter(f::Function, b::TensorBasis) 
	states = filter(f, b.states)
	b_arr = Array(Basis{kind(b)}, size(b,2))
	for i=1:size(b)[2]
		b_arr[i] = Basis("?", unique(map(x->getindex(x,i), states)))
    end
	TensorBasis(b_arr, states)
end

function map(f::Function, b::AbstractBasis) 
	newstates = map(f, b.states)
	if eltype(newstates) <: TensorState
		return TensorBasis(newstates)
	else
		return Basis(newstates)
	end
end

ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states), b.statemap)
ctranspose(b::Basis) = Basis(b.label, map(ctranspose, b.states), b.statemap)

size(b::TensorBasis) = (length(b.states), length(b.bases))
length(b::AbstractBasis) = length(b.states)
endof(b::AbstractBasis) = endof(b.states)

function basisjoin{K}(b::Basis{K}, s::State{K})
	resmap = copy(b.statemap)
	resmap[(label(s), basislabel(s))] = length(b)+1
	if samebasis(b, s)
		return Basis(b.label, vcat(b.states, s), resmap)
	else
		return Basis("?", vcat(b.states, s), resmap)
	end
end

basisjoin{K}(s::State{K}, b::Basis{K}) = Basis(vcat(s, b.states))

function basisjoin{K}(a::Basis{K}, b::Basis{K})
	resmap = copy(a.statemap)
	for i=1:length(b)
		resmap[(label(b[i]), basislabel(b[i]))] = length(b)+i
	end
	if samebasis(a,b)
		return Basis(a.label, vcat(a.states, b.states), resmap)
	else
		return Basis("?", vcat(a.states, b.states), resmap)
	end
end

function basisjoin{K}(b::TensorBasis{K}, s::TensorState{K})
	if samebasis(b, s)
		resmap = copy(b.statemap)
		resmap[(label(s), basislabel(s))] = length(b)+1
		return TensorBasis(b.bases, vcat(b.states, s), resmap)
	else
		return TensorBasis(vcat(b.states, s))
	end
end

basisjoin{K}(s::TensorState{K}, b::TensorBasis{K}) = TensorBasis(vcat(s, b.states))

function basisjoin{K}(a::TensorBasis{K}, b::TensorBasis{K})
	if samebasis(a, b)
		resmap = copy(b.statemap)
		for i=1:length(b)
			resmap[(label(b[i]), basislabel(b[i]))] = length(b)+i
		end
		return TensorBasis(a.bases, vcat(a.states, b.states), resmap)
	else
		return TensorBasis(vcat(a.states, b.states))
	end
end


*{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = tensor(a,b)
+{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = basisjoin(a,b)

setdiff{B<:AbstractBasis}(a::B,b::B) = setdiff(a.states, b.states)

getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

get(b::AbstractBasis, s::AbstractState, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
get(b::Basis, s::AbstractState, notfound::String) = get(b.statemap, (label(s), basislabel(s)), notfound) #this definition resolves ambiguities
get(b::AbstractBasis, s::AbstractState) = b.statemap[(label(s), basislabel(s))]
get(b::Basis, label, basislabel::String) = b.statemap[(label, basislabel)]
get(b::TensorBasis, label::Vector, basislabel::Vector{String}) = b.statemap[(label, basislabel)]

showcompact(io::IO, b::AbstractBasis) = print(io, "$(typeof(b)) $(reprlabel(b))")

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
tobasis(s::State) = Basis(s)
tobasis(s::TensorState) = TensorBasis(map(Basis, separate(s)), [s])
tobasis{S<:State}(v::Vector{S}) = Basis(v)
tobasis{S<:TensorState}(v::Vector{S}) = TensorBasis(v)

isdual(a::Basis{Ket}, b::Basis{Bra}) = label(a)==label(b) && a.statemap==b.statemap
isdual(a::Basis{Bra}, b::Basis{Ket}) = isdual(b,a)
isdual(a::TensorBasis{Ket}, b::TensorBasis{Bra}) = label(a)==label(b) && a.statemap==b.statemap
isdual(a::TensorBasis{Bra}, b::TensorBasis{Ket}) = isdual(b,a)
isdual(a::AbstractBasis,b::AbstractBasis)=false #default to false



separate(b::Basis)=b
separate(b::TensorBasis) = b.bases
kind(b::AbstractBasis) = kind(b[1])
