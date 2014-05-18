immutable State{K<:BraKet} <: AbstractState{K}
  label
  basislabel::Symbol
  kind::Type{K}
end

State{K<:BraKet}(label, basislabel::Symbol=:?, kind::Type{K}=Ket) = State{kind}(label, basislabel, kind)
State{K<:BraKet}(label, kind::Type{K}) = State{K}(label, :?, kind)

immutable TensorState{K<:BraKet} <: AbstractState{K}
  states::Vector{State{K}}
  kind::Type{K}
end

TensorState{K<:BraKet}(states::Vector{State{K}}, kind::Type{K}=Ket) = TensorState{kind}(states, kind)
TensorState{K<:BraKet}(labels::Vector, basislabel::Symbol=:?, kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, basislabel, kind), kind)
TensorState{K<:BraKet}(labels::Vector, kind::Type{K}=Ket) = TensorState(labels, :?, kind)

immutable InnerProduct
	bra::AbstractState{Bra}
	ket::AbstractState{Ket}
end

isequal(a::State,b::State) = isequal(a.label, b.label) && a.basislabel==b.basislabel && a.kind==b.kind
==(a::State,b::State) = a.label==b.label && a.basislabel==b.basislabel && a.kind==b.kind

hash(s::State) = hash(s.label)+hash(s.basislabel)+hash(s.kind)
ctranspose(s::State) = State(s.label, s.basislabel, !s.kind)
ctranspose(s::TensorState) = TensorState(map(ctranspose, s.states), !s.kind)
getindex(s::TensorState, x) = s.states[x]

reprlabel(s::State) = "$(repr(s.label))_$(repr(s.basislabel))"
function reprlabel(s::TensorState)
	str = "$(reprlabel(s.states[1]))"
	for i=2:length(s.states)
		str = "$(str) , $(reprlabel(s.states[i]))"
	end
	return str
end

show(io::IO, s::AbstractState{Ket}) = print(io, "| $(reprlabel(s)) $rang")
show(io::IO, s::AbstractState{Bra}) = print(io, "$lang $(reprlabel(s)) |")
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");


kind(s::AbstractState) = s.kind
basislabel(s::State) = s.basislabel
label(s::State) = s.label
basislabel(s::TensorState) = map(basislabel, s.states)
label(s::TensorState) = map(label, s.states)

for op=(:length, :endof)
	@eval ($op)(s::TensorState) = $(op)(s.states)
end

*{K<:BraKet}(a::State{K}, b::State{K}) = TensorState([a,b], K)
*{K<:BraKet}(a::TensorState{K}, b::State{K}) = TensorState(vcat(a.states, b), K)
*{K<:BraKet}(a::State{K}, b::TensorState{K}) = TensorState(vcat(a, b.states), K)
*{K<:BraKet}(a::TensorState{K}, b::TensorState{K}) = TensorState(vcat(a.states, b.states), K)

function *(s::AbstractState, n::Number)
	if n==1 
		return s
	elseif n==0 
		return 0
	else 
		error("implement StateRep")
	end
end

*(n::Number, s::AbstractState) = *(s, n)

function *(a::State{Bra}, b::State{Ket})
	if a.basislabel==b.basislabel && a.basislabel!=:?
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
	ind = findfirst(s->s.basislabel==b.basislabel, reverse(a.states))
	ind==0 ? InnerProduct(a,b) : sprod(a,b,ind)
end

function *(a::State{Bra}, b::TensorState{Ket}) 
	ind = findfirst(s->s.basislabel==a.basislabel, b.states)
	ind==0 ? InnerProduct(a,b) : sprod(a,b,ind)
end

*(a::TensorState{Bra}, b::TensorState{Ket}) = InnerProduct(a, b)

function sprod{K<:BraKet}(a::TensorState{Bra}, b::TensorState{Ket}, i::Int, target::Type{K}=Ket)
	if target==Ket
		return reduce(*, vcat(a*b[i], b[1:i-1], b[i+1:end]))
	else 
		return reduce(*, vcat(a[1:i-1], a[i+1:end], a[i]*b))
	end
end
sprod(a::State{Bra}, b::TensorState{Ket}, i::Int) = reduce(*, vcat(a*b[i], b[1:i-1], b[i+1:end]))
sprod(a::TensorState{Bra}, b::State{Ket}, i::Int) = reduce(*, vcat(a[1:i-1], a[i+1:end], a[i]*b))

tensor() = nothing
tensor(s::AbstractState) = s
tensor{K<:BraKet}(s::AbstractState{K}...) = reduce(*, s) 
tensor{S<:AbstractState}(arr::Array{S}) = tensor(arr...)

tensorarr(arrs::Array...) = crossjoin(arrs...)

statearr{K<:BraKet}(arr::Array, basislabel::Symbol=:?, kind::Type{K}=Ket) = State{K}[State(arr[i], basislabel, kind) for i in arr]
statearr{K<:BraKet}(arr::Array, kind::Type{K}) = statearr(arr, :?, kind)
statejoin{S<:AbstractState}(state_arr::Array{S,2}) = [reduce(*, state_arr[i, :]) for i=1:size(state_arr, 1)]

separate{K<:BraKet}(s::TensorState{K}) = s.states
