#####################################
#State###############################
#####################################
abstract State{S} #S<:Single

immutable Ket <: State{Ket}
	label
	bsym::Symbol
end

immutable Bra <: State{Bra}
	label
	bsym::Symbol
end

typealias Single Union(Ket,Bra)

immutable Tensor{S<:Single} <: State{S}
	states::Vector{S}
end

#####################################
#Misc Functions######################
#####################################

copy{S<:Single}(s::S) = S(copy(s.label), copy(s.bsym))
copy{S<:Tensor}(s::S) = S(copy(s.states))

eltype{K}(t::Type{Tensor{K}}) = K
dual(t::Type{Ket}) = Bra
dual(t::Type{Bra}) = Ket
dual{K}(t::Type{Tensor{K}}) = Tensor{dual(K)}

isequal{S<:Single}(a::S,b::S) = isequal(a.label, b.label) && a.bsym==b.bsym
=={S<:Single}(a::S,b::S) = a.label==b.label && a.bsym==b.bsym

isequal{S<:Tensor}(a::S, b::S) = isequal(a.states, b.states)
=={S<:Tensor}(a::S, b::S) = a.states==b.states

hash(s::Tensor) = sum(map(hash, s.states))

ctranspose{S<:Single}(s::S) = dual(S)(s.label, s.bsym)
ctranspose{S<:Tensor}(s::S) = dual(S)([ctranspose(i) for i in s.states])

getindex(s::Tensor, x) = s.states[x]
separate(s::Tensor) = s.states

bsym(s::Single) = s.bsym
label(s::Single) = s.label
bsym(s::Tensor) = map(bsym, s.states)
label(s::Tensor) = [label(i) for i in s.states]

labeldelta(a::Single, b::Single) = label(a)==label(b) ? 1 : 0
labeldelta(a::Tensor, b::Tensor) = label(a)==label(b) ? 1 : 0
labeldelta(a::State, b::State) = 0 #default to 0

isdual(a::Bra, b::Ket) = label(a)==label(b) && samebasis(a,b)
isdual(a::Ket, b::Bra) = isdual(b,a)
isdual(a::Tensor{Ket}, b::Tensor{Bra}) = label(a)==label(b) && samebasis(a,b)
isdual(a::Tensor{Bra}, b::Tensor{Ket}) = isdual(b,a)
isdual(a::State, b::State) = false #default to false

for op=(:length, :endof, :eltype)
	@eval ($op)(s::Tensor) = $(op)(s.states)
end

statearr(arr::Array, bsym::Symbol, K::Type{Single}=Ket) = map(i->K(i, bsym), arr)

#####################################
#Show Functions######################
#####################################

reprlabel(s::Single) = "$(repr(s.label)):$(s.bsym)"
function reprlabel(s::Tensor)
	str = "$(reprlabel(s.states[1]))"
	for i=2:length(s.states)
		str = "$(str), $(reprlabel(s.states[i]))"
	end
	return str
end

show(io::IO, s::State{Ket}) = print(io, "| $(reprlabel(s)) $rang")
show(io::IO, s::State{Bra}) = print(io, "$lang $(reprlabel(s)) |")

#####################################
#Arithmetic Operations###############
#####################################

# *(c::DiracCoeff, s::AbstractState) = c==0 ? 0 : (c==1 ? s : DiracVector([c], tobasis(s)))
# *(s::AbstractState, c::DiracCoeff) = *(c,s)
tensor(s::State) = error("cannot perform tensor operation on one state")
tensor{S<:Single}(s::Vector{S}) = Tensor(s)
tensor(s::Array) = tensor(vec(s))
tensor{S<:Single}(a::S, b::S) = Tensor([a,b]) 
tensor{S<:Single}(a::S, b::Tensor{S}) = Tensor([a,b.states]) 
tensor{S<:Single}(a::Tensor{S}, b::S) = Tensor([a.states,b]) 
tensor{S<:Single}(a::Tensor{S}, b::Tensor{S}) = Tensor([a.states,b.states]) 
tensor{S<:Single}(s::State{S}...) = reduce(tensor,s) 
*{S<:Single}(a::State{S}, b::State{S}) = tensor(a,b)

function inner(a::Bra, b::Ket)
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		return InnerProduct(a, b)
	end
end
inner(a::Bra, b::Tensor{Ket}, i::Int) = inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
inner(a::Tensor{Bra}, b::Ket, i::Int) = tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
inner(a::Tensor{Bra}, b::Ket) = inner(a[1],b)*tensor(a[2:end])
inner(a::Bra, b::Tensor{Ket}) = inner(a,b[1])*tensor(b[2:end])

function inner(a::Tensor{Bra}, b::Tensor{Ket}, i::Int, target::Type{Single}=Ket)
	if target<:Ket
		return inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
	else 
		return tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
	end
end

function inner(a::Tensor{Bra}, b::Tensor{Ket})
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		for s in a.states
			b = *(s,b)
		end
		return b
	end
end

*(a::State{Bra}, b::State{Ket}) = inner(a,b)
*(a::State{Ket}, b::State{Bra}) = OuterProduct(a,b)
