#####################################
#Basis###############################
#####################################

abstract AbstractBasis{S<:Single} <: Dirac

immutable Basis{S<:Single} <: AbstractBasis{S}
	states::Vector{S}
	statemap::Dict
	adjhash::Uint #hash statemap for quick isdual() functionality 
end

statemapper{S<:State}(sv::Vector{S}) = (labeltype(S)=>Int)[label(sv[i])=>i for i=1:length(sv)]
hashbasis{S<:Single}(sv::Vector{S}, statemap::Dict) = Basis(sv, statemap, hash(statemap))
consbasis{S<:Single}(sv::Vector{S}) = hashbasis(sv, statemapper(sv))
consbasis(bsym::Symbol,labelvec::Array) = consbasis(svec(bsym,labelvec))

function makebasis{S<:Single}(sv::Vector{S})
	@assert bsym(eltype(sv))!=Any "all states must have same bsym to construct Basis"
	@assert labeltype(eltype(sv))!=Any "all states must have same labeltype to construct Basis"	
	consbasis(sv)
end

basis{S<:Single}(states::Vector{S}) = makebasis(unique(states))
basis{S<:Single}(states::Array{S}) = basis(vec(states))

basis(labelvec::Array, bsym::Symbol) = basis(svec(bsym, labelvec))

basis{S<:Single}(states::Array{S}) = basis(vec(states))
basis(bsym::Symbol, labelvec::Array) = basis(svec(bsym, labelvec))

#####################################
#TensorBasis#########################
#####################################

immutable TensorBasis{S<:Single,B<:Basis} <: AbstractBasis{S}
	bases::Vector{B}
	states::Vector{Tensor{S}}
	statemap::Dict
end

consbasis{B<:Basis,T<:Tensor}(bv::Vector{B}, sv::Vector{T}) = TensorBasis(bv, sv, statemapper(sv))

separate(b::Basis) = [b]
separate(b::TensorBasis) = b.bases

kind{K<:Ket}(b::AbstractBasis{K}) = Ket
kind{B<:Bra}(b::AbstractBasis{B}) = Bra

sepstates(sv) = hcat(map(separate, sv)...)

function consbvec!(bases, seps)
	for i=1:size(seps, 1)
	 	bases[i] = basis([seps[i, :]...]) #forcing correct typing with vcat
	end
	bases
end

function consbvec!{B<:Basis}(bases, seps::Vector{B})
	for i=1:length(seps)
	 	bases[i] = seps[i]
	end
	bases
end

prepbvec{S<:Single}(sv::Vector{Tensor{S}}) = Array(Basis{S}, length(sv[1]))
prepbvec{S<:Tensor}(sv::Vector{S}) = Array(Basis, length(sv[1]))

maketensorbasis{S<:Tensor}(sv::Vector{S}) = consbasis(consbvec!(prepbvec(sv),sepstates(sv)), sv)

basis{S<:Tensor}(states::Vector{S}) = maketensorbasis(unique(states))
basis{S<:Tensor}(states::Array{S}) = basis(vec(states))

#improved type inferencing for statejoin; 
#based on the idea that each column of sarr 
#has its own type due to the behavior of 
#crossjoin
sjointype(sarr::Array)=typejoin(map(eltype, sarr[1,:])...)

statejoin{S<:Single}(sarr::Array{S,2}) = Tensor{S}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]
statejoin{S<:State}(sarr::Array{S,2}) = Tensor{sjointype(sarr)}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]

statecross(v) = statejoin(reduce(crossjoin,v))

tensor() = error("tensor needs arguments of type T<:Basis or T<:State")
tensor(b::AbstractBasis) = b

tensor(b::AbstractBasis...) = consbasis([[separate(i) for i in b]...], statecross([i.states for i in b]))

for t=(:Ket,:Bra)
	@eval begin
	tensor{S1<:($t), S2<:($t)}(s::State{S1}, b::AbstractBasis{S2}) = tensor(basis(s), b)
	tensor{S1<:($t), S2<:($t)}(b::AbstractBasis{S1}, s::State{S2}) = tensor(b, basis(s))
	end
end

#####################################
#Misc Functions######################
#####################################
basis{S<:State}(s::S...) = basis(collect(s))
basis(s::State...) = basis(kind(s[1])[s...])
basis(arr::Array) = basis([arr...]) #vcat forces retyping


copy{B<:Basis}(b::B) = B(copy(b.states), copy(b.statemap), copy(b.adjhash))
copy{B<:TensorBasis}(b::B)= B(copy(b.bases), copy(b.states), copy(b.statemap))
for t=(:Ket,:Bra)
	@eval bsym{b,T}(bs::Basis{($t){b,T}}) = b
end
bsym(b::TensorBasis) = map(bsym, b.bases)
adjhash(b::Basis) = b.adjhash
adjhash(b::TensorBasis) = map(adjhash, b.bases)

for op=(:isequal, :(==))
	@eval ($op){S<:Single}(a::Basis{S}, b::Basis{S}) = adjhash(a)==adjhash(b)
	@eval ($op){b,T}(x::TensorBasis{Bra{b,T}}, y::TensorBasis{Bra{b,T}}) = adjhash(x)==adjhash(y)
	@eval ($op){b,T}(x::TensorBasis{Ket{b,T}}, y::TensorBasis{Ket{b,T}}) = adjhash(x)==adjhash(y)
	@eval ($op){S1<:Ket,S2<:Ket}(a::TensorBasis{S1}, b::TensorBasis{S2}) = adjhash(a)==adjhash(b) && samebasis(a,b)
	@eval ($op){S1<:Bra,S2<:Bra}(a::TensorBasis{S1}, b::TensorBasis{S2}) = adjhash(a)==adjhash(b) && samebasis(a,b)
end

for op=(:length, :endof)
	@eval ($op)(b::AbstractBasis) = ($op)(b.states)
end

dual{S<:Single}(t::Type{Basis{S}}) = Basis{dual(S)}
dual{S<:Single}(t::Type{TensorBasis{S}}) = TensorBasis{dual(S)}

ctranspose(b::Basis) = Basis([ctranspose(i) for i in b.states], b.statemap, b.adjhash)
ctranspose(b::TensorBasis) = TensorBasis([ctranspose(i) for i in b.bases], [ctranspose(i) for i in b.states], b.statemap)

isdual{b,T}(x::Basis{Ket{b,T}}, y::Basis{Bra{b,T}}) = adjhash(x)==adjhash(y)
isdual{b,T}(x::Basis{Bra{b,T}}, y::Basis{Ket{b,T}}) = isdual(y,x)
isdual{b,T}(x::TensorBasis{Ket{b,T}}, y::TensorBasis{Bra{b,T}}) = adjhash(x)==adjhash(y)
isdual{b,T}(x::TensorBasis{Bra{b,T}}, y::TensorBasis{Ket{b,T}}) = isdual(y,x)

isdual{S1<:Ket,S2<:Bra}(a::TensorBasis{S1}, b::TensorBasis{S2}) = adjhash(a)==adjhash(b) && samebasis(a,b)
isdual{S1<:Bra,S2<:Ket}(a::TensorBasis{S1}, b::TensorBasis{S2}) = isdual(b,a)


isdual(a::AbstractBasis,b::AbstractBasis)=false 

getindex(b::AbstractBasis, x) = b.states[x]

for t=(:Ket,:Bra)
	@eval begin
	get{b,T}(bs::Basis{($t){b,T}}, s::($t){b,T}, notfound) = get(bs.statemap, label(s), notfound)
	get{b,T}(bs::Basis{($t){b,T}}, s::($t){b,T}) = bs.statemap[label(s)]
	get{b,T}(bs::TensorBasis{($t){b,T}}, s::Tensor{($t){b,T}}, notfound) = get(bs.statemap, label(s), notfound)
	get{b,T}(bs::TensorBasis{($t){b,T}}, s::Tensor{($t){b,T}}) = bs.statemap[label(s)]

	get{S1<:($t),S2<:($t)}(b::TensorBasis{S1}, s::Tensor{S2}, notfound) = samebasis(b,s) ? get(b.statemap, label(s), notfound) : notfound
	get{S1<:($t),S2<:($t)}(b::TensorBasis{S1}, s::Tensor{S2}) = samebasis(b,s) ? b.statemap[label(s)] : throw(KeyError(s))

	in{b,T}(s::($t){b,T}, bs::Basis{($t){b,T}})= get(bs,s,"FALSE")=="FALSE" ? false : true
	in{S1<:($t),S2<:($t)}(s::Tensor{S1}, b::TensorBasis{S2})= get(b,s,"FALSE")=="FALSE" ? false : true
	end
end

get(b::AbstractBasis, s::State, notfound) = notfound
get(b::AbstractBasis, s::State) = throw(KeyError(s))	

getpos(b::AbstractBasis, args...) = get(b, ars...)

in(s::State, b::AbstractBasis)= false	
in(s, b::AbstractBasis)= error("in(s,b::AbstractBasis) is defined for s of type State, not $(typeof(s))")

######################################
##Show Functions######################
######################################

reprlabel(b::Basis) = bsym(b)
function reprlabel(b::TensorBasis)
	labels = bsym(b)
	#terrible way to grow a string
	str = "$(labels[1])"
	for i=2:length(labels)
		str = "$str$otimes$(labels[i])"
	end
	return str
end

summary(io::IO, b::AbstractBasis)= print(io, "$(typeof(b))")

showcompact(io::IO, b::AbstractBasis) = (summary(io, b); show(io, b.states))

function show(io::IO, b::AbstractBasis)
	summary(io, b)
	println(", $(length(b)) states:")
	if length(b)>20
		for i=1:10
			println(io, b.states[i])
		end
		println(vdots)
		show(io, b.states[length(b)-10])
		for i=length(b)-9:length(b)
			print('\n')
			show(io, b.states[i])
		end
	else
		show(io, b.states[1])
		for i=2:length(b.states)
			print('\n')
			show(io, b.states[i])
		end	
	end
end

######################################
##Function-Passing Functions##########
######################################

find(f::Function, b::Basis) = find(f, b.states)
filter(f::Function, b::Basis) = consbasis(filter(f, b.states))
filter(f::Function, b::TensorBasis) = maketensorbasis(filter(f, b.states))
map(f::Function, b::AbstractBasis) = basis(map(f, b.states))

######################################
##Joining/Separating Functions########
######################################

function stateappend(d::Dict,s::State)
	resmap = copy(d)
	resmap[label(s)] = length(d)+1
	return resmap
end

function bcat{S}(a::Basis{S}, b::Basis{S})
	consbasis(unique(vcat(a.states, b.states)))
end

function bcat{S}(b::Basis{S}, s::S)
	if in(s, b)
		return b
	else
		return hashbasis(vcat(b.states, s), stateappend(b.statemap,s))
	end
end

function bcat{S}(s::S, b::Basis{S})
	if in(s, b)
		return b
	else
		return consbasis(vcat(s,b.states))
	end
end

function bcat{S}(b::TensorBasis{S}, s::Tensor{S})
	if in(s, b)
		return b
	else
		@assert samebasis(b, s) "BasisMismatch"
		return TensorBasis(b.bases, vcat(b.states, s), stateappend(b.statemap,s))
	end
end

function bcat{S}(s::Tensor{S},b::TensorBasis{S})
	if in(s, b)
		return b
	else
		@assert samebasis(b, s) "BasisMismatch"
		return consbasis(b.bases, vcat(s,b.states))
	end
end

bcat{S}(a::TensorBasis{S}, b::TensorBasis{S}) = basis(vcat(a.states, b.states))

kron(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
kron(a::AbstractBasis, b::State) = tensor(a,b)
kron(a::State, b::AbstractBasis) = tensor(a,b)

setdiff(a::AbstractBasis,b::AbstractBasis) = setdiff(a.states, b.states)
