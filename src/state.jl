#####################################
#State###############################
#####################################
abstract State{S} #S<:Single

immutable Ket{T} <: State{Ket{T}}
	label::T
	bsym::Symbol
end

immutable Bra{T} <: State{Bra{T}}
	label::T
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

labeltype{T}(t::Type{Ket{T}}) = T
labeltype{T}(t::Type{Bra{T}}) = T
labeltype(t::Type{Ket}) = Any
labeltype(t::Type{Bra}) = Any
labeltype{T}(t::Type{Tensor{Ket{T}}}) = T
labeltype{T}(t::Type{Tensor{Bra{T}}}) = T
labeltype(t::Type{Tensor{Ket}}) = Vector{Any}
labeltype(t::Type{Tensor{Bra}}) = Vectr{Any}
labeltype(s::State) = labeltype(typeof(s))

eltype{S<:Single}(t::Type{S}) = S
eltype{K}(t::Type{Tensor{K}}) = K
eltype(s::State) = eltype(typeof(s))

dual(t::Type{Ket}) = Bra
dual(t::Type{Bra}) = Ket
dual{T}(t::Type{Ket{T}}) = Bra{T}
dual{T}(t::Type{Bra{T}}) = Ket{T}
dual{K}(t::Type{Tensor{K}}) = Tensor{dual(K)}

isequal{S<:Single}(a::S, b::S) = isequal(a.label, b.label) && a.bsym==b.bsym
=={S<:Single}(a::S,b::S) = a.label==b.label && a.bsym==b.bsym

isequal{S<:Tensor}(a::S, b::S) = isequal(a.states, b.states)
=={S<:Tensor}(a::S, b::S) = a.states==b.states

hash(s::Tensor) = sum(map(hash, s.states))

ctranspose{S<:Single}(s::S) = dual(S)(s.label, s.bsym)
ctranspose{S<:Tensor}(s::S) = dual(S)([ctranspose(i) for i in s.states])

getindex(s::Tensor, x) = s.states[x]
separate(s::Single) = s
separate(s::Tensor) = s.states

bsym(s::Single) = s.bsym
label(s::Single) = s.label
bsym(s::Tensor) = map(bsym, s.states)
label(s::Tensor) = [label(i) for i in s.states]

labeldelta(a::Single, b::Single) = label(a)==label(b) ? 1 : 0
labeldelta(a::Tensor, b::Tensor) = label(a)==label(b) ? 1 : 0
labeldelta(a::State, b::State) = 0 #default to 0

isdual{B<:Bra, K<:Ket}(a::State{B},b::State{K}) = label(a)==label(b) && samebasis(a,b)
isdual{K<:Ket, B<:Bra}(a::State{K},b::State{B}) = isdual(b,a)
isdual(a::State, b::State) = false #default to false

for op=(:length, :endof, :eltype)
	@eval ($op)(s::Tensor) = $(op)(s.states)
end

statevec{S<:Single}(arr::Array, bsym::Symbol, K::Type{S}=Ket) = [K(i, bsym) for i in arr]

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

show{K<:Ket}(io::IO, s::State{K}) = print(io, "| $(reprlabel(s)) $rang")
show{B<:Bra}(io::IO, s::State{B}) = print(io, "$lang $(reprlabel(s)) |")

#####################################
#Arithmetic Operations###############
#####################################

# *(c::DiracCoeff, s::AbstractState) = c==0 ? 0 : (c==1 ? s : DiracVector([c], tobasis(s)))
# *(s::AbstractState, c::DiracCoeff) = *(c,s)
tensor(s::State) = error("cannot perform tensor operation on one state")
for t=(:Ket, :Bra)
	@eval begin
	tensor(a::($t), b::($t)) = Tensor([a,b]) 
	tensor{S<:($t)}(a::($t), b::Tensor{S}) = Tensor([a,b.states]) 
	tensor{S<:($t)}(a::Tensor{S}, b::($t)) = Tensor([a.states,b]) 
	tensor{S1<:($t), S2<:($t)}(a::Tensor{S1}, b::Tensor{S2}) = Tensor([a.states,b.states]) 
	*{S1<:($t), S2<:($t)}(a::State{S1}, b::State{S2}) = tensor(a,b)
	end
end
tensor{S<:Single}(s::Vector{S}) = Tensor(s)
tensor{S<:State}(s::Vector{S}) = Tensor([[separate(i) for i in s]...])
tensor(s::State...) = tensor(collect(s)) 
tensor{S<:State}(s::Array{S}) = tensor(vec(s))

function inner(a::Bra, b::Ket)
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		return InnerProduct(a, b)
	end
end
inner{K<:Ket}(a::Bra, b::Tensor{K}, i::Int) = inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
inner{B<:Bra}(a::Tensor{B}, b::Ket, i::Int) = tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
inner{K<:Ket}(a::Bra, b::Tensor{K}) = inner(a,b[1])*tensor(b[2:end])
inner{B<:Bra}(a::Tensor{B}, b::Ket) = inner(a[1],b)*tensor(a[2:end])

function inner{B<:Bra,K<:Ket,S<:Single}(a::Tensor{B}, b::Tensor{K}, i::Int, target::Type{S}=Ket)
	if target<:Ket
		return inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
	else 
		return tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
	end
end

function inner{B<:Bra,K<:Ket}(a::Tensor{B}, b::Tensor{K})
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		for s in a.states
			b = *(s,b)
		end
		return b
	end
end

*{B<:Bra,K<:Ket}(a::State{B}, b::State{K}) = inner(a,b)
*{K<:Ket,B<:Bra}(a::State{K}, b::State{B}) = OuterProduct(a,b)
