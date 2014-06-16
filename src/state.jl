#####################################
#State/TensorState###################
#####################################

immutable State{K<:BraKet, T} <: AbstractState{K}
  label::T
  basislabel::Symbol
end

State{K<:BraKet,T}(label::T, basislabel::Symbol, kind::Type{K}=Ket) = State{kind, T}(label, basislabel)

immutable TensorState{K<:BraKet} <: AbstractState{K}
  states::Vector
  TensorState{S<:State{K}}(v::Vector{S}) = new(v)
  TensorState(v::Vector) = new(convert(Vector{State{K}},v))  
end

TensorState{K<:BraKet}(labels::Vector, basislabel::Symbol, kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, basislabel, kind))

#####################################
#Misc Functions######################
#####################################

copy(s::State) = State(copy(s.label), copy(s.basislabel))
copy(s::TensorState) = TensorState(copy(s.states))

#hashing isn't working right now :(
# hash(s::State) = ?
# hash(s::TensorState) = ?

isequal{K,T}(a::State{K,T},b::State{K,T}) = isequal(a.label, b.label) && a.basislabel==b.basislabel
=={K,T}(a::State{K,T},b::State{K,T}) = a.label==b.label && a.basislabel==b.basislabel

isequal{K}(a::TensorState{K}, b::TensorState{K}) = isequal(a.states, b.states)
=={K}(a::TensorState{K}, b::TensorState{K}) = a.states==b.states

isequal(a::AbstractState, b::AbstractState) = false #default to false
==(a::AbstractState, b::AbstractState) = false #default to false

ctranspose{K,T}(s::State{K,T}) = State{!K,T}(s.label, s.basislabel)
ctranspose{K}(s::TensorState{K}) = TensorState{!K}([ctranspose(i) for i in s.states])

getindex(s::TensorState, x) = s.states[x]
separate(s::TensorState) = s.states

kind{K}(s::AbstractState{K}) = K

basislabel(s::State) = s.basislabel
label(s::State) = s.label
basislabel(s::TensorState) = map(basislabel, s.states)
label(s::TensorState) = [label(i) for i in s.states]

labeldelta(a::State, b::State) = label(a)==label(b) ? 1 : 0
labeldelta(a::TensorState, b::TensorState) = label(a)==label(b) ? 1 : 0
labeldelta(a::AbstractState, b::AbstractState) = 0 #default to 0

isdual{T}(a::State{Bra, T}, b::State{Ket, T}) = label(a)==label(b) && samebasis(a,b)
isdual{T}(a::State{Ket, T}, b::State{Bra, T}) = isdual(b,a)
isdual(a::TensorState{Ket}, b::TensorState{Bra}) = label(a)==label(b) && samebasis(a,b)
isdual(a::TensorState{Bra}, b::TensorState{Ket}) = isdual(b,a)
isdual(a::AbstractState, b::AbstractState) = false #default to false

for op=(:length, :endof)
	@eval ($op)(s::TensorState) = $(op)(s.states)
end

#####################################
#Show Functions######################
#####################################

reprlabel(s::State) = "$(repr(s.label)):$(s.basislabel)"
function reprlabel(s::TensorState)
	str = "$(reprlabel(s.states[1]))"
	for i=2:length(s.states)
		str = "$(str), $(reprlabel(s.states[i]))"
	end
	return str
end

show(io::IO, s::AbstractState{Ket}) = print(io, "| $(reprlabel(s)) $rang")
show(io::IO, s::AbstractState{Bra}) = print(io, "$lang $(reprlabel(s)) |")

#####################################
#Arithmetic Operations###############
#####################################
tensor() = error("tensor takes arguments of states")
tensor{K}(a::State{K}, b::State{K}) = TensorState{K}(vcat(a,b))
tensor{K}(a::TensorState{K}, b::State{K}) = TensorState{K}(vcat(a.states, b))
tensor{K}(a::State{K}, b::TensorState{K}) = TensorState{K}(vcat(a, b.states))
tensor{K}(a::TensorState{K}, b::TensorState{K}) = TensorState{K}(vcat(a.states, b.states))
tensor{K}(s::AbstractState{K}...) = reduce(tensor,s) 

*{K}(a::AbstractState{K}, b::AbstractState{K}) = tensor(a,b)

*(c::DiracCoeff, s::AbstractState) = c==0 ? 0 : (c==1 ? s : DiracVector([c], tobasis(s)))
*(s::AbstractState, c::DiracCoeff) = *(c,s)

function inner(a::State{Bra}, b::State{Ket})
	if samebasis(a,b)
		label(a)==label(b) ? 1 : 0
	else
		return InnerProduct(a, b)
	end
end
inner(a::State{Bra}, b::TensorState{Ket}, i::Int) = inner(a, b[i])*reduce(tensor,vcat(b[1:i-1], b[i+1:end]))
inner(a::TensorState{Bra}, b::State{Ket}, i::Int) = reduce(tensor,vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
inner(a::TensorState{Bra}, b::State{Ket}) = inner(a[1],b)*reduce(tensor,a[2:end])
inner(a::State{Bra}, b::TensorState{Ket}) = inner(a,b[1])*reduce(tensor,b[2:end])

function inner{K<:BraKet}(a::TensorState{Bra}, b::TensorState{Ket}, i::Int, target::Type{K}=Ket)
	if target==Ket
		return inner(a, b[i])*reduce(tensor,vcat(b[1:i-1], b[i+1:end]))
	else 
		return reduce(tensor,vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
	end
end

function inner(a::TensorState{Bra}, b::TensorState{Ket})
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		for s in a.states
			b = *(s,b)
		end
		return b
	end
end

*(a::AbstractState{Bra}, b::AbstractState{Ket}) = inner(a,b)
*(a::AbstractState{Ket}, b::AbstractState{Bra}) = OuterProduct(a,b)

#####################################
#Functions###########################
#####################################

function statearr{K<:BraKet, T}(arr::Array{T}, basislabel::Symbol, kind::Type{K}=Ket) 
	if is(T,Any)
		return convert(Array{State{kind}}, map(i->State(i, basislabel, kind), arr))
	else
		return map(i->State(i, basislabel, kind), arr)
	end
end

statejoin{S<:AbstractState}(state_arr::Array{S,2}) = TensorState{kind(S)}[reduce(tensor,state_arr[i, :]) for i=1:size(state_arr, 1)]

