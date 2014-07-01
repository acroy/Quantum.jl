#####################################
#State###############################
#####################################
abstract State{S} <: Dirac

immutable Ket{T,b} <: State{Ket{T,b}}
	label::T
end

immutable Bra{T,b} <: State{Bra{T,b}}
	label::T
end

Ket{T}(label::T, b) = Ket{T,b}(label)
Bra{T}(label::T, b) = Bra{T,b}(label)

typealias Single Union(Bra, Ket)

immutable Tensor{S<:Single} <: State{S}
	states::Vector{S}
end 

#####################################
#Misc Functions######################
#####################################

in(s::Single, t::Tensor) = false
for kb=(:Ket,:Bra)
	@eval begin
	in{S1<:($kb),S2<:($kb)}(s::S1, t::Tensor{S2}) = in(s, t.states)
	bsym{T}(s::Type{($kb){T}}) = Any
	bsym{T,b}(s::Type{($kb){T,b}}) = b
	bsym{T,b}(s::($kb){T,b}) = b
	labeltype{T,b}(t::Type{($kb){T,b}}) = T
	labeltype(t::Type{($kb)}) = Any
	labeltype{T,b}(t::Type{Tensor{($kb){T,b}}}) = Vector{T}
	labeltype(t::Type{Tensor{($kb)}}) = Vector{Any}
	end
end

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

kind{K<:Ket}(b::State{K}) = Ket
kind{B<:Bra}(b::State{B}) = Bra

dual(t::Type{Ket}) = Bra
dual(t::Type{Bra}) = Ket
dual{T,b}(t::Type{Ket{T,b}}) = Bra{T,b}
dual{T,b}(t::Type{Bra{T,b}}) = Ket{T,b}
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
isdual{T1,T2,b}(x::State{Bra{T1,b}},y::State{Ket{T2,b}}) = label(x)==label(y)
isdual{T1,T2,b}(x::State{Ket{T1,b}},y::State{Bra{T2,b}}) = isdual(y,x)

for op=(:length, :endof)
	@eval ($op)(s::Tensor) = $(op)(s.states)
end

svec{T}(arr::Array{T}, bsym::Symbol) = map(Ket{T,bsym}, arr)

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

tensor{S<:Single}(s::Vector{S}) = Tensor(s)
tensor{S<:State}(s::Vector{S}) = tensor([[separate(i) for i in s]...])
tensor{S<:State}(s::Array{S}) = tensor(vec(s))
tensor(s::Vector) = tensor([[separate(i) for i in s]...])
tensor(s::State...) = tensor(collect(s)) 
tensor(s::State) = s 


inner(x::Bra, y::Ket) = InnerProduct(x, y)
inner{T1,T2,b}(x::Bra{T1,b}, y::Ket{T2,b}) = labeldelta(x,y)

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

*{K1<:Ket,K2<:Ket}(a::State{K1}, b::State{K2}) = tensor(a,b)
*{B1<:Bra,B2<:Bra}(a::State{B1}, b::State{B2}) = tensor(a,b)
*{B<:Bra,K<:Ket}(a::State{B}, b::State{K}) = inner(a,b)
*{K<:Ket,B<:Bra}(a::State{K}, b::State{B}) = OuterProduct(a,b)