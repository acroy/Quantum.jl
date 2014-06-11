#####################################
#Basis###############################
#####################################
immutable Basis{K<:BraKet} <: AbstractBasis{K}
	label::String
	states::Vector{State{K}}
	statemap::Dict{(Any,String), Int}
	function Basis(label, states, statemap, errcheck=true)
		if errcheck
			@assert length(unique(states))==length(states) "Basis states must be uniquely labeled"
				new(label, states, statemap)
		else
			new(label, states, statemap)
		end
	end
end

function makebasis{K}(label::String, states::Array{State{K}})
	states = unique(states)
	return Basis{K}(label, states, ((Any,String)=>Int)[(states[i].label,states[i].basislabel)=>i for i=1:length(states)], false)
end

Basis{K<:BraKet}(labelvec::Array, label::String, kind::Type{K}=Ket) = makebasis(label, statearr(labelvec, label, kind))
Basis{K}(s::State{K}...) = Basis(vcat(s...))
function Basis{K<:BraKet}(s::Array{State{K}}) 
	bases = unique(map(basislabel, s)) 
	@assert length(bases)>1 bmm
		makebasis(bases[1], s)
	end
end

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

#####################################
#Misc Functions######################
#####################################

copy{K}(b::Basis{K}) = Basis{K}(copy(b.label), copy(b.states), copy(b.statemap))
copy(b::TensorBasis) = TensorBasis(copy(b.bases), copy(b.states), copy(b.statemap))

kind(b::AbstractBasis) = kind(b[1])

label(b::Basis) = b.label
label(b::TensorBasis) = map(label, b.bases)
basislabel(b::AbstractBasis) = label(b)

isequal(a::AbstractBasis, b::AbstractBasis) = isequal(a.states, b.states) && label(a)==label(b)
==(a::AbstractBasis, b::AbstractBasis) = a.states==b.states && label(a)==label(b)

ctranspose(b::TensorBasis) = TensorBasis(map(ctranspose, b.bases), map(ctranspose, b.states), b.statemap)
ctranspose{K}(b::Basis{K}) = Basis{!K}(b.label, map(ctranspose, b.states), b.statemap, false)

isdual(a::Basis{Ket}, b::Basis{Bra}) = label(a)==label(b) && a.statemap==b.statemap
isdual(a::Basis{Bra}, b::Basis{Ket}) = isdual(b,a)
isdual(a::TensorBasis{Ket}, b::TensorBasis{Bra}) = label(a)==label(b) && a.statemap==b.statemap
isdual(a::TensorBasis{Bra}, b::TensorBasis{Ket}) = isdual(b,a)
isdual(a::AbstractBasis,b::AbstractBasis)=false #default to false

size(b::TensorBasis) = (length(b.states), length(b.bases))
length(b::AbstractBasis) = length(b.states)
endof(b::AbstractBasis) = endof(b.states)

getindex(b::AbstractBasis, x) = b.states[x]
setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

get(b::AbstractBasis, s::AbstractState, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
get(b::Basis, s::AbstractState, notfound::String) = get(b.statemap, (label(s), basislabel(s)), notfound) #this definition resolves ambiguities
get(b::AbstractBasis, s::AbstractState) = b.statemap[(label(s), basislabel(s))]
get(b::Basis, label, basislabel::String) = b.statemap[(label, basislabel)]
get(b::TensorBasis, label::Vector, basislabel::Vector{String}) = b.statemap[(label, basislabel)]

in(s::AbstractState, b::AbstractBasis)= get(b,s,"FALSE")=="FALSE" ? false : true


#####################################
#Show Functions######################
#####################################

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
#####################################
#Function-Passing Functions##########
#####################################

find(f::Function, b::Basis) = find(f, b.states)
filter(f::Function, b::Basis) = Basis(filter(f, b.states))
filter(f::Function, b::TensorBasis) = TensorBasis(filter(f, b.states))

function map(f::Function, b::AbstractBasis) 
	newstates = map(f, b.states)
	if eltype(newstates) <: TensorState
		return TensorBasis(newstates)
	else
		return Basis(newstates)
	end
end

#####################################
#Joining/Separating Functions########
#####################################

function tensor{K}(bases::AbstractBasis{K}...)
	TensorBasis(vcat([separate(i) for i in bases]...), convert(Vector{TensorState{K}}, statejoin(tensorarr([i.states for i in bases]...))))
end

function tensor{K}(basis::AbstractBasis{K})
	return basis
end

function basisjoin{K}(b::Basis{K}, s::State{K})
	if in(s, b)
		return b
	else
		@assert samebasis(b, s) bmm
		resmap = copy(b.statemap)
		resmap[(label(s), basislabel(s))] = length(b)+1
		return Basis{K}(b.label, vcat(b.states, s), resmap, false)
	end
end

basisjoin{K}(s::State{K}, b::Basis{K}) = Basis(vcat(s, b.states))

function basisjoin{K}(a::Basis{K}, b::Basis{K})
	@assert samebasis(a,b) bmm
	resmap = copy(a.statemap)
	for i=1:length(b)
		if !in(b[i], a)
			resmap[(label(b[i]), basislabel(b[i]))] = length(b)+i
		end
	end
	return Basis{K}(a.label, unique(vcat(a.states, b.states)), resmap, false)
end

function basisjoin{K}(b::TensorBasis{K}, s::TensorState{K})
	@assert samebasis(b, s) bmm
	resmap = copy(b.statemap)
	resmap[(label(s), basislabel(s))] = length(b)+1
	return TensorBasis(b.bases, vcat(b.states, s), resmap)
end

basisjoin{K}(s::TensorState{K}, b::TensorBasis{K}) = TensorBasis(vcat(s, b.states))

function basisjoin{K}(a::TensorBasis{K}, b::TensorBasis{K})
	@assert samebasis(a, b) bmm
	resmap = copy(b.statemap)
	for i=1:length(b)
		resmap[(label(b[i]), basislabel(b[i]))] = length(b)+i
	end
	return TensorBasis(a.bases, vcat(a.states, b.states), resmap)
end

tensor{K}(a::AbstractBasis{K}, b::AbstractState{K}) = map(s->s*b, d.basis)
tensor{K}(a::AbstractState{K}, b::AbstractBasis{K}) = map(s->b*s, d.basis)

*{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = tensor(a,b)
*{K}(a::AbstractBasis{K}, b::AbstractState{K}) = tensor(a,b)
*{K}(a::AbstractState{K}, b::AbstractBasis{K}) = tensor(a,b)
+{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = basisjoin(a,b)
+{K}(a::AbstractBasis{K}, b::AbstractState{K}) = basisjoin(a,b)
+{K}(a::AbstractState{K}, b::AbstractBasis{K}) = basisjoin(a,b)

setdiff{B<:AbstractBasis}(a::B,b::B) = setdiff(a.states, b.states)

tobasis(s::State) = Basis(s)
tobasis(s::TensorState) = TensorBasis(map(Basis, separate(s)), [s])
tobasis(s::AbstractState...) = tobasis(vcat(s...))
tobasis{S<:State}(v::Array{S}) = Basis(vec(v))
tobasis{S<:TensorState}(v::Array{S}) = TensorBasis(vec(v))

separate(b::Basis)=[b]
separate(b::TensorBasis) = b.bases
