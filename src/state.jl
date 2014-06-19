#####################################
#State###############################
#####################################
abstract Single
abstract Tensor
abstract State{S<:Union(Single, Tensor)} <: Dirac
abstract AbstractKet{S<:Union(Single, Tensor)} <: State{S}
abstract AbstractBra{S<:Union(Single, Tensor)} <: State{S}

immutable Ket{T} <: AbstractKet{Single}
	label::T
	bsym::Symbol
end

immutable Bra{T} <: AbstractBra{Single}
	label::T
	bsym::Symbol
end

immutable TensorKet{K<:Ket} <: AbstractKet{Tensor}
	states::Vector{K}
end

immutable TensorBra{B<:Bra} <: AbstractBra{Tensor}
	states::Vector{B}
end

#####################################
#Misc Functions######################
#####################################

copy{S<:State{Single}}(s::S) = S(copy(s.label), copy(s.bsym))
copy{S<:State{Tensor}}(s::S) = S(copy(s.states))

eltype{K}(t::Type{TensorKet{K}}) = K
eltype{B}(t::Type{TensorBra{B}}) = B
dual(t::Type{Ket}) = Bra
dual(t::Type{Bra}) = Ket
dual(t::Type{TensorKet}) = TensorBra
dual(t::Type{TensorBra}) = TensorKet
dual{T}(t::Type{Ket{T}}) = Bra{T}
dual{T}(t::Type{Bra{T}}) = Ket{T}
dual{T}(t::Type{TensorKet{T}}) = TensorBra{dual(T)}
dual{T}(t::Type{TensorBra{T}}) = TensorKet{dual(T)}

isequal{S<:State{Single}}(a::S,b::S) = isequal(a.label, b.label) && a.bsym==b.bsym
=={S<:State{Single}}(a::S,b::S) = a.label==b.label && a.bsym==b.bsym

isequal{S<:State{Tensor}}(a::S, b::S) = isequal(a.states, b.states)
=={S<:State{Tensor}}(a::S, b::S) = a.states==b.states

ctranspose{S<:State{Single}}(s::S) = dual(S)(s.label, s.bsym)
ctranspose{S<:State{Tensor}}(s::S) = dual(S)([ctranspose(i) for i in s.states])

getindex(s::State{Tensor}, x) = s.states[x]
separate(s::State{Tensor}) = s.states

bsym(s::State{Single}) = s.bsym
label(s::State{Single}) = s.label
bsym(s::State{Tensor}) = map(bsym, s.states)
label(s::State{Tensor}) = [label(i) for i in s.states]

labeldelta{S<:State{Single}}(a::S, b::S) = label(a)==label(b) ? 1 : 0
labeldelta{S<:State{Tensor}}(a::S, b::S) = label(a)==label(b) ? 1 : 0
labeldelta{T}(a::Ket{T}, b::Bra{T}) = label(a)==label(b) ? 1 : 0
labeldelta{T}(a::Bra{T}, b::Ket{T}) = label(a)==label(b) ? 1 : 0
labeldelta{T}(a::TensorKet{Ket{T}}, b::TensorBra{Bra{T}}) = label(a)==label(b) ? 1 : 0
labeldelta{T}(a::TensorBra{Bra{T}}, b::TensorKet{Ket{T}}) = label(a)==label(b) ? 1 : 0
labeldelta(a::State, b::State) = 0 #default to 0

isdual{T}(a::Bra{T}, b::Ket{T}) = label(a)==label(b) && samebasis(a,b)
isdual{T}(a::Ket{T}, b::Bra{T}) = label(a)==label(b) && samebasis(a,b)
isdual{T}(a::TensorKet{Ket{T}}, b::TensorBra{Bra{T}}) = label(a)==label(b) && samebasis(a,b)
isdual{T}(a::TensorBra{Bra{T}}, b::TensorKet{Ket{T}}) = label(a)==label(b) && samebasis(a,b)
isdual(a::TensorKet{Ket}, b::TensorBra{Bra}) = label(a)==label(b) && samebasis(a,b)
isdual(a::TensorBra{Bra}, b::TensorKet{Ket}) = label(a)==label(b) && samebasis(a,b)
isdual(a::State, b::State) = false #default to false

for op=(:length, :endof, :eltype)
	@eval ($op)(s::State{Tensor}) = $(op)(s.states)
end

#####################################
#Show Functions######################
#####################################

reprlabel(s::State{Single}) = "$(repr(s.label)):$(s.bsym)"
function reprlabel(s::State{Tensor})
	str = "$(reprlabel(s.states[1]))"
	for i=2:length(s.states)
		str = "$(str), $(reprlabel(s.states[i]))"
	end
	return str
end

show(io::IO, s::AbstractKet) = print(io, "| $(reprlabel(s)) $rang")
show(io::IO, s::AbstractBra) = print(io, "$lang $(reprlabel(s)) |")

#####################################
#Arithmetic Operations###############
#####################################

# *(c::DiracCoeff, s::AbstractState) = c==0 ? 0 : (c==1 ? s : DiracVector([c], tobasis(s)))
# *(s::AbstractState, c::DiracCoeff) = *(c,s)
tensor(s::State) = error("cannot perform tensor operation on one state")
tensor(s::Array) = tensor(vec(s))
tensor{K<:Ket}(s::Vector{K}) = TensorKet(s)
tensor{B<:Bra}(s::Vector{B}) = TensorBra(s)
tensor(a::Ket, b::Ket) = TensorKet([a,b]) 
tensor(a::Ket, b::TensorKet) = TensorKet([a,b.states]) 
tensor(a::TensorKet, b::Ket) = TensorKet([a.states,b]) 
tensor(a::TensorKet, b::TensorKet) = TensorKet([a.states,b.states]) 
tensor(a::Bra, b::Bra) = TensorBra([a,b]) 
tensor(a::Bra, b::TensorBra) = TensorBra([a,b.states]) 
tensor(a::TensorBra, b::Bra) = TensorBra([a.states,b]) 
tensor(a::TensorBra, b::TensorBra) = TensorBra([a.states,b.states]) 

tensor(s::State...) = reduce(tensor,s) 

*(a::AbstractKet, b::AbstractKet) = tensor(a,b)
*(a::AbstractBra, b::AbstractBra) = tensor(a,b)

function inner(a::Bra, b::Ket)
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		return InnerProduct(a, b)
	end
end
inner(a::Bra, b::TensorKet, i::Int) = inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
inner(a::TensorBra, b::Ket, i::Int) = tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
inner(a::TensorBra, b::Ket) = inner(a[1],b)*tensor(a[2:end])
inner(a::Bra, b::TensorKet) = inner(a,b[1])*tensor(b[2:end])

function inner{S<:State{Single}}(a::TensorBra, b::TensorKet, i::Int; target::Type{S}=Ket)
	if target<:Ket
		return inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end]))
	else 
		return tensor(vcat(a[1:i-1], a[i+1:end]))*inner(a[i], b)
	end
end

function inner(a::TensorBra, b::TensorKet)
	if samebasis(a,b)
		return labeldelta(a,b)
	else
		for s in a.states
			b = *(s,b)
		end
		return b
	end
end

*(a::AbstractBra, b::AbstractKet) = inner(a,b)
*(a::AbstractKet, b::AbstractBra) = OuterProduct(a,b)

#####################################
#Utility Functions###################
#####################################

function statearr{T,S<:State{Single}}(arr::Array{T}, bsym::Symbol, K::Type{S}=Ket) 
	if is(T,Any)
	 	return convert(Array{K}, map(i->K(i, bsym), arr))
	else
		return map(i->K(i, bsym), arr)
	end
end

statejoin{K<:Ket}(sarr::Array{K,2}) = TensorKet{K}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]
statejoin{B<:Bra}(sarr::Array{B,2}) = TensorBra{B}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]
statejoin{K<:AbstractKet}(sarr::Array{K,2}) = TensorKet{Ket}[reduce(tensor,sarr[i, :]) for i=1:size(sarr, 1)]
statejoin{B<:AbstractBra}(sarr::Array{B,2}) = TensorBra{Bra}[reduce(tensor,sarr[i, :]) for i=1:size(sarr, 1)]
