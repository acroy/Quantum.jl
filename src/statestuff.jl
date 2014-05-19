immutable State{K<:BraKet} <: AbstractState{K}
  label
  basislabel::String
  kind::Type{K}
end

State{K<:BraKet}(label, basislabel::String="?", kind::Type{K}=Ket) = State{kind}(label, basislabel, kind)
State{K<:BraKet}(label, kind::Type{K}) = State{K}(label, "?", kind)

immutable TensorState{K<:BraKet} <: AbstractState{K}
  states::Vector{State{K}}
  kind::Type{K}
end

TensorState{K<:BraKet}(states::Vector{State{K}}, kind::Type{K}=Ket) = TensorState{kind}(states, kind)
TensorState{K<:BraKet}(labels::Vector, basislabel::String="?", kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, basislabel, kind), kind)
TensorState{K<:BraKet}(labels::Vector, kind::Type{K}=Ket) = TensorState(labels, "?", kind)

immutable InnerProduct <: Dirac
	bra::AbstractState{Bra}
	ket::AbstractState{Ket}
end

immutable OuterProduct <: Dirac
	ket::AbstractState{Ket}
	bra::AbstractState{Bra}
end

isequal(a::State,b::State) = isequal(a.label, b.label) && a.basislabel==b.basislabel && a.kind==b.kind
==(a::State,b::State) = a.label==b.label && a.basislabel==b.basislabel && a.kind==b.kind

ctranspose(s::State) = State(s.label, s.basislabel, !s.kind)
ctranspose(s::TensorState) = TensorState(map(ctranspose, s.states), !s.kind)
ctranspose(i::InnerProduct) = InnerProduct(i.ket', i.bra')
ctranspose(o::OuterProduct) = OuterProduct(o.bra', o.ket')

getindex(s::TensorState, x) = s.states[x]

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
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");
show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");


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
	if a.basislabel==b.basislabel && a.basislabel!="?" && b.basislabel!="?"
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
	ind==0 ? InnerProduct(a,b) : inner(a,b,ind)
end

function *(a::State{Bra}, b::TensorState{Ket}) 
	ind = findfirst(s->s.basislabel==a.basislabel, b.states)
	ind==0 ? InnerProduct(a,b) : inner(a,b,ind)
end

function *(a::TensorState{Bra}, b::TensorState{Ket})
	for s in reverse(a.states)
		b = s*b
	end
	return b
end
*(a::AbstractState{Bra}, b::InnerProduct) = InnerProduct(a*b.bra, b.ket)
*(a::InnerProduct, b::AbstractState{Ket}) = InnerProduct(a.bra, a.ket*b)

function inner{K<:BraKet}(a::TensorState{Bra}, b::TensorState{Ket}, i::Int, target::Type{K}=Ket)
	if target==Ket
		return reduce(*, vcat(a*b[i], b[1:i-1], b[i+1:end]))
	else 
		return reduce(*, vcat(a[1:i-1], a[i+1:end], a[i]*b))
	end
end
inner(a::State{Bra}, b::TensorState{Ket}, i::Int) = reduce(*, vcat(a*b[i], b[1:i-1], b[i+1:end]))
inner(a::TensorState{Bra}, b::State{Ket}, i::Int) = reduce(*, vcat(a[1:i-1], a[i+1:end], a[i]*b))

*(a::AbstractState{Ket}, b::AbstractState{Bra}) = OuterProduct(a,b)

tensor() = nothing
tensor(s::AbstractState) = s
tensor{K<:BraKet}(s::AbstractState{K}...) = reduce(*, s) 
tensor{S<:AbstractState}(arr::Array{S}) = tensor(arr...)

tensorarr(arrs::Array...) = crossjoin(arrs...)

statearr{K<:BraKet}(arr::Array, basislabel::String="?", kind::Type{K}=Ket) = State{K}[State(i, basislabel, kind) for i in arr]
statearr{K<:BraKet}(arr::Array, kind::Type{K}) = statearr(arr, "?", kind)
statejoin{S<:AbstractState}(state_arr::Array{S,2}) = [reduce(*, state_arr[i, :]) for i=1:size(state_arr, 1)]

separate{K<:BraKet}(s::TensorState{K}) = s.states
