immutable State{K<:BraKet} <: AbstractState{K}
  label
  eigop::Symbol
  kind::Type{K}
end

State{K<:BraKet}(label, eigop::Symbol=:?, kind::Type{K}=Ket) = State{kind}(label, eigop, kind)
State{K<:BraKet}(label, kind::Type{K}) = State{K}(label, :?, kind)

#if a state's eigop is :?, then it will be converted to a number 
#corresponding to its position in the TensorState
immutable TensorState{K<:BraKet} <: AbstractState{K}
  states::Vector{State{K}}
  kind::Type{K}
end

TensorState{K<:BraKet}(states::Vector{State{K}}, kind::Type{K}=Ket) = TensorState{kind}(states, kind)
TensorState{K<:BraKet}(labels::Vector, eigop::Symbol=:?, kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, eigop, kind), kind)
TensorState{K<:BraKet}(labels::Vector, kind::Type{K}=Ket) = TensorState(labels, :?, kind)

isequal(a::State,b::State) = isequal(a.label, b.label) && a.eigop==b.eigop && a.kind==b.kind
==(a::State,b::State) = a.label==b.label && a.eigop==b.eigop && a.kind==b.kind

hash(s::State) = hash(s.label)+hash(s.eigop)+hash(s.kind)
ctranspose(s::State) = State(s.label, s.eigop, !s.kind)
ctranspose(s::TensorState) = TensorState(map(ctranspose, s.states), !s.kind)


getindex(s::TensorState, x) = s.states[x]

reprlabel(s::State) = "$(repr(s.label))_$(repr(s.eigop))"
function reprlabel(s::TensorState)
	str = "$(reprlabel(s.states[1]))"
	for i=2:length(s.states)
		str = "$(str) , $(reprlabel(s.states[i]))"
	end
	return str
end
show(io::IO, s::AbstractState{Ket}) = print(io, "| $(reprlabel(s)) $rang")
show(io::IO, s::AbstractState{Bra}) = print(io, "$lang $(reprlabel(s)) |")

kind(s::AbstractState) = s.kind
eigop(s::State) = s.eigop
label(s::State) = s.label
eigop(s::TensorState) = map(eigop, s.states)
label(s::TensorState) = map(label, s.states)



*(s::State{Bra}, op::Symbol) = op==eigop(s) ? error("implement StateRep") : error("$s is not an eigenstate of operator $(repr(op))")
*(op::Symbol, s::State{Ket}) = op==eigop(s) ? error("implement StateRep") : error("$s is not an eigenstate of operator $(repr(op))")
*(s::AbstractState{Ket}, op::Symbol) = error("Cannot apply Ket to operator. Perhaps you meant $(repr(op)) * $s")
*(op::Symbol, s::AbstractState{Bra}) = error("Cannot apply operator to bra. Perhaps you meant $s * $(repr(op))")

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
	if a.eigop==b.eigop && a.eigop!=:?
		if a.label == b.label
			return 1
		else
			return 0
		end
	else
		error("implement InnerProduct")
	end
end

function *(a::TensorState{Bra}, b::State{Ket}) 
	ind = findfirst(s->s*b==1, reverse(a.states))-1
	tensor(vcat(a.states[1:(length(a.states)-ind)-1], a.states[(length(a.states)-ind)+1:end]))
end

function *(a::State{Bra}, b::TensorState{Ket}) 
	ind = findfirst(s->s*b==1, a.states)
	tensor(vcat(a.states[1:ind-1], a.states[ind+1:end]))
end

*(a::TensorState{Bra}, b::TensorState{Ket}) = reduce(*, vcat(a.states, b.states))


tensor() = 1 #in case inner products happen between singlet TensorStates
tensor(s::AbstractState) = s
tensor{K<:BraKet}(s::AbstractState{K}...) = reduce(*, s) 
tensor{S<:AbstractState}(arr::Array{S}) = tensor(arr...)

tensorarr(arrs::Array...) = crossjoin(arrs...)

statearr{K<:BraKet}(arr::Array, eigop::Symbol=:?, kind::Type{K}=Ket) = State{K}[State(arr[i], eigop, kind) for i in arr]
statearr{K<:BraKet}(arr::Array, kind::Type{K}) = statearr(arr, :?, kind)
statejoin{S<:AbstractState}(state_arr::Array{S,2}) = [reduce(*, state_arr[i, :]) for i=1:size(state_arr, 1)]

separate{K<:BraKet}(s::TensorState{K}) = s.states
