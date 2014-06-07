#####################################
#State/TensorState###################
#####################################

immutable State{K<:BraKet} <: AbstractState{K}
  label
  basislabel::String
  kind::Type{K}
end

State{K<:BraKet}(label, basislabel::String = "?", kind::Type{K}=Ket) = State{kind}(label, basislabel, kind)
State{K<:BraKet}(label, kind::Type{K}) = State{K}(label, "?", kind)

immutable TensorState{K<:BraKet} <: AbstractState{K}
  states::Vector{State{K}}
  kind::Type{K}
end

TensorState{K<:BraKet}(states::Vector{State{K}}, kind::Type{K}=Ket) = TensorState{kind}(states, kind)
TensorState{K<:BraKet}(labels::Vector, basislabel::String="?", kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, basislabel, kind), kind)
TensorState{K<:BraKet}(labels::Vector, kind::Type{K}=Ket) = TensorState(labels, "?", kind)

#####################################
#Misc Functions######################
#####################################

copy(s::State) = State(copy(s.label), copy(s.basislabel), copy(s.kind))
copy(s::TensorState) = TensorState(copy(s.states), copy(s.kind))

hash(s::State) = hash((s.label, s.basislabel, s.kind))
hash(s::TensorState) = hash((s.states, s.kind))

isequal(a::State,b::State) = isequal(a.label, b.label) && a.basislabel==b.basislabel && a.kind==b.kind
==(a::State,b::State) = a.label==b.label && a.basislabel==b.basislabel && a.kind==b.kind

isequal(a::TensorState, b::TensorState) = isequal(a.states, b.states) && a.kind==b.kind
==(a::TensorState, b::TensorState) = a.states==b.states && a.kind==b.kind

ctranspose(s::State) = State(s.label, s.basislabel, !s.kind)
ctranspose(s::TensorState) = TensorState(map(ctranspose, s.states), !s.kind)

getindex(s::TensorState, x) = s.states[x]

kind(s::AbstractState) = s.kind
basislabel(s::State) = s.basislabel
label(s::State) = s.label
basislabel(s::TensorState) = map(basislabel, s.states)
label(s::TensorState) = map(label, s.states)

isdual(a::State{Ket}, b::State{Bra}) = label(a)==label(b) && samebasis(a,b)
isdual(a::State{Bra}, b::State{Ket}) = isdual(b,a)
isdual(a::TensorState{Ket}, b::TensorState{Bra}) = label(a)==label(b) && samebasis(a,b)
isdual(a::TensorState{Bra}, b::TensorState{Ket}) = isdual(b,a)
isdual(a::AbstractState, b::AbstractState) = false #default to false

for op=(:length, :endof)
	@eval ($op)(s::TensorState) = $(op)(s.states)
end

#####################################
#Show Functions######################
#####################################

reprlabel(s::State) = "$(repr(s.label))_$(s.basislabel)"
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

*{K}(a::State{K}, b::State{K}) = TensorState([a,b], K)
*{K}(a::TensorState{K}, b::State{K}) = TensorState(vcat(a.states, b), K)
*{K}(a::State{K}, b::TensorState{K}) = TensorState(vcat(a, b.states), K)
*{K}(a::TensorState{K}, b::TensorState{K}) = TensorState(vcat(a.states, b.states), K)
*(c::DiracCoeff, s::AbstractState) = c==0 ? 0 : (c==1 ? s : DiracVector([c], tobasis(s)))
*(s::AbstractState, c::DiracCoeff) = *(c,s)

function *(a::State{Bra}, b::State{Ket})
	if samebasis(a, b) && !samebasis(a, "?") && !samebasis(b, "?")
		if a.label == b.label
			return 1
		else
			return 0
		end
	else
		return InnerProduct(a, b)
	end
end

function *(a::TensorState{Bra}, b::State{Ket}) 
	ind = findfirst(s->samebasis(s,b), reverse(a.states))
	ind==0 ? prod(a[1:end-1])*(last(a)*b) : inner(a,b,ind)
end

function *(a::State{Bra}, b::TensorState{Ket}) 
	ind = findfirst(s->samebasis(s,a), b.states)
	ind==0 ? (a*b[1])*prod(b[2:end]) : inner(a,b,ind)
end

function *(a::TensorState{Bra}, b::TensorState{Ket})
	for s in a.states
		b = s*b
	end
	return b
end

*(a::AbstractState{Ket}, b::AbstractState{Bra}) = OuterProduct(a,b)

#####################################
#Product Functions###################
#####################################

function inner{K<:BraKet}(a::TensorState{Bra}, b::TensorState{Ket}, i::Int, target::Type{K}=Ket)
	if target==Ket
		return prod(vcat(a*b[i], b[1:i-1], b[i+1:end]))
	else 
		return prod(vcat(a[1:i-1], a[i+1:end], a[i]*b))
	end
end

inner(a::State{Bra}, b::TensorState{Ket}, i::Int) = prod(vcat(a*b[i], b[1:i-1], b[i+1:end]))
inner(a::TensorState{Bra}, b::State{Ket}, i::Int) = prod(vcat(a[1:i-1], a[i+1:end], a[i]*b))

tensor() = nothing
tensor(s::AbstractState) = s
tensor{K}(s::AbstractState{K}...) = prod(s) 
tensor{S<:AbstractState}(arr::Array{S}) = tensor(arr...)

tensorarr(arrs::Array...) = crossjoin(arrs...)

statearr{K<:BraKet}(arr::Array, basislabel::String="?", kind::Type{K}=Ket) = State{K}[State(i, basislabel, kind) for i in arr]
statearr{K<:BraKet}(arr::Array, kind::Type{K}) = statearr(arr, "?", kind)
statejoin{S<:AbstractState}(state_arr::Array{S,2}) = [prod(state_arr[i, :]) for i=1:size(state_arr, 1)]

separate(s::TensorState) = s.states