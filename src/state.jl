#####################################
#State###############################
#####################################
abstract State{S} <: Dirac

immutable Ket{b,T} <: State{Ket{b,T}}
	label::T
end

immutable Bra{b,T} <: State{Bra{b,T}}
	label::T
end

ket{T}(b,label::T) = Ket{b,T}(label)
bra{T}(b,label::T) = Bra{b,T}(label)

typealias Single Union(Bra, Ket)

immutable Tensor{S<:Single} <: State{S}
	states::Vector{S}
end 

tensor{S<:Single}(s::Vector{S}) = Tensor(s)
tensor{S<:State}(s::Vector{S}) = tensor([[separate(i) for i in s]...])
tensor{S<:State}(s::Array{S}) = tensor(vec(s))
tensor(s::Vector) = tensor([[separate(i) for i in s]...])
tensor(s::State...) = tensor(collect(s)) 
tensor(s::State) = s 

#####################################
#Misc Functions######################
#####################################

in(s::Single, t::Tensor) = false
for kb=(:Ket,:Bra)
	@eval begin
	in{S1<:($kb),S2<:($kb)}(s::S1, t::Tensor{S2}) = in(s, t.states)

	bsym{b,T}(s::Type{($kb){b,T}}) = b
	bsym{b,T}(s::($kb){b,T}) = b
	bsym{b}(s::Type{($kb){b}}) = b
	bsym{b}(s::($kb){b}) = b
	bsym(s::Type{($kb)}) = Any
	bsym(s::($kb)) = Any

	labeltype{b,T}(t::Type{($kb){b,T}}) = T
	labeltype{b}(t::Type{($kb){b}}) = Any
	labeltype(t::Type{($kb)}) = Any
	
	kind{b,T}(s::Type{State{($kb){b,T}}}) = ($kb){b}
	kind{b}(s::Type{State{($kb){b}}}) = ($kb){b}
	kind(s::Type{State{($kb)}}) = ($kb)
	kind{b,T}(s::State{($kb){b,T}}) = ($kb){b}
	kind{b}(s::State{($kb){b}}) = ($kb){b}
	kind(s::State{($kb)}) = ($kb)
	end
end

labeltype{S}(s::Type{Tensor{S}}) = Vector{labeltype(S)}
labeltype(s::State) = labeltype(typeof(s))

label(s::Single) = s.label
bsym(s::Tensor) = map(bsym, s.states)
label(s::Tensor) = [label(i) for i in s.states]

copy{S<:Single}(s::S) = S(copy(s.label))
copy{S<:Tensor}(s::S) = S(copy(s.states))

labeltype(s::State) = labeltype(typeof(s))
in(s::State, t::State) = in(s, t.states)

eltype{S<:Single}(t::Type{S}) = S
eltype{S}(t::Type{Tensor{S}}) = S
eltype(s::State) = eltype(typeof(s))

dual(t::Type{Ket}) = Bra
dual(t::Type{Bra}) = Ket
dual{b,T}(t::Type{Ket{b,T}}) = Bra{b,T}
dual{b,T}(t::Type{Bra{b,T}}) = Ket{b,T}
dual{S}(t::Type{Tensor{S}}) = Tensor{dual(S)}

isequal{S<:Single}(a::S, b::S) = isequal(label(a), label(b))
=={S<:Single}(a::S,b::S) = label(a)==label(b)

isequal{S<:Tensor}(a::S, b::S) = isequal(a.states, b.states)
=={S<:Tensor}(a::S, b::S) = a.states==b.states

hash(s::Tensor) = sum(map(hash, s.states))

ctranspose{S<:Single}(s::S) = dual(S)(label(s))
ctranspose{S<:Tensor}(s::S) = dual(S)([ctranspose(i) for i in s.states])

getindex(s::Tensor, x) = s.states[x]
separate(s::Single) = s
separate(s::Tensor) = s.states

labeldelta(a::State, b::State) = 0 #default to 0
labeldelta(a::Single, b::Single) = label(a)==label(b) ? 1 : 0
labeldelta(a::Tensor, b::Tensor) = label(a)==label(b) ? 1 : 0

isdual(a::State, b::State) = false #default to false
isdual{b,T}(x::Bra{b,T},y::Ket{b,T}) = label(x)==label(y)
isdual{b,T}(x::Ket{b,T},y::Bra{b,T}) = isdual(y,x)
isdual{b,T}(x::Tensor{Bra{b,T}},y::Tensor{Ket{b,T}}) = label(x)==label(y)
isdual{b,T}(x::Tensor{Ket{b,T}},y::Tensor{Bra{b,T}}) = isdual(y,x)
isdual{b}(x::Tensor{Bra{b}},y::Tensor{Ket{b}}) = label(x)==label(y)
isdual{b}(x::Tensor{Ket{b}},y::Tensor{Bra{b}}) = isdual(y,x)
isdual(x::Tensor{Bra},y::Tensor{Ket}) = label(x)==label(y) && samebasis(x,y)
isdual(x::Tensor{Ket},y::Tensor{Bra}) = isdual(y,x)

for op=(:length, :endof)
	@eval ($op)(s::Tensor) = $(op)(s.states)
end

svec{T}(bsym::Symbol, arr::Array{T}) = map(Ket{bsym,T}, arr)

#####################################
#Show Functions######################
#####################################

reprlabel(s::Single) = "$(repr(label(s))):$(bsym(s))"
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

inner(x::Bra, y::Ket) = InnerProduct(x, y)
inner{b}(x::Bra{b}, y::Ket{b}) = labeldelta(x,y)

inner{K<:Ket}(a::Bra, b::Tensor{K}, i::Int) = inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end])...)
inner{B<:Bra}(a::Tensor{B}, b::Ket, i::Int) = tensor(vcat(a[1:i-1], a[i+1:end])...)*inner(a[i], b)
inner{K<:Ket}(a::Bra, b::Tensor{K}) = inner(a,b[1])*tensor(b[2:end]...)
inner{B<:Bra}(a::Tensor{B}, b::Ket) = inner(a[1],b)*tensor(a[2:end]...)

function inner{B<:Bra,K<:Ket,S<:Single}(a::Tensor{B}, b::Tensor{K}, i::Int, target::Type{S}=Ket)
	if target<:Ket
		return inner(a, b[i])*tensor(vcat(b[1:i-1], b[i+1:end])...)
	else 
		return tensor(vcat(a[1:i-1], a[i+1:end])...)*inner(a[i], b)
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

for t=(:Ket, :Bra)
	@eval begin
	kron{A<:($t),B<:($t)}(a::State{A}, b::State{B}) = tensor(a,b)
	end
end
kron{K<:Ket,B<:Bra}(a::State{K}, b::State{B}) = OuterProduct(a,b)
kron{B<:Bra,K<:Ket}(a::State{B}, b::State{K}) = OuterProduct(b,a)

*{B<:Bra,K<:Ket}(a::State{B}, b::State{K}) = inner(a,b)
*{K<:Ket,B<:Bra}(a::State{K}, b::State{B}) = kron(a,b)
