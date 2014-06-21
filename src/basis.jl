#####################################
#Basis###############################
#####################################

abstract AbstractBasis{S<:Single} <: Dirac

immutable Basis{S<:Single} <: AbstractBasis{S}
	bsym::Symbol
	states::Vector{S}
	statemap::Dict
end

statemapper{S<:State}(sv::Vector{S}) = (labeltype(S)=>Int)[label(sv[i])=>i for i=1:length(sv)]

consbasis{S<:Single}(sv::Vector{S}) = Basis(bsym(sv[1]), sv, statemapper(sv))

function makebasis{S<:Single}(sv::Vector{S})
	@assert length(unique(map(bsym, sv)))==1 "BasisMismatch"
	consbasis(sv)
end

basis{S<:Single}(states::Vector{S}) = makebasis(unique(states))
basis{S<:Single}(states::Array{S}) = basis(vec(states))
basis{S<:Single}(s::S...) = basis(collect(s))

#####################################
#TensorBasis#########################
#####################################

immutable TensorBasis{S<:Single} <: AbstractBasis{S}
	bases::Vector
	states::Vector{Tensor{S}}
	statemap::Dict
end

consbasis{S<:Single,B<:Basis}(bv::Vector{B}, sv::Vector{Tensor{S}}) = TensorBasis(bv, sv, statemapper(sv))

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

prepbvec{T}(sv::Vector{Tensor{Ket{T}}}) = Array(Basis{Ket{T}}, length(sv[1]))
prepbvec{T}(sv::Vector{Tensor{Bra{T}}}) = Array(Basis{Bra{T}}, length(sv[1]))
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

statecross(v::Vector) = statejoin(reduce(crossjoin, v))

tensor() = error("tensor needs arguments of type T<:Basis or T<:State")
tensor(b::AbstractBasis) = error("cannot perform tensor operation on one basis")

function tensor(b::AbstractBasis...)
	@assert length(unique(map(kind, b))) == 1 "arguments to tensor must be of uniform kind (B<:Bra or K<:Ket)"
	sv = statecross([i.states for i in b])
	TensorBasis([[separate(i) for i in b]...], sv, statemapper(sv))
end

#####################################
#Misc Functions######################
#####################################

copy{B<:Basis}(b::B) = B(copy(b.bsym), copy(b.states), copy(b.statemap))
copy{B<:TensorBasis}(b::B)= B(copy(b.bases), copy(b.states), copy(b.statemap))

bsym(b::Basis) = b.bsym
bsym(b::TensorBasis) = map(bsym, b.bases)

isequal{B<:AbstractBasis}(a::B, b::B) = isequal(a.states, b.states)
=={B<:AbstractBasis}(a::B, b::B) = a.states==b.states

for op=(:length, :endof)
	@eval ($op)(b::AbstractBasis) = ($op)(b.states)
end

dual{S<:Single}(t::Type{Basis{S}}) = Basis{dual(S)}
dual{S<:Single}(t::Type{TensorBasis{S}}) = TensorBasis{dual(S)}

ctranspose(b::Basis) = Basis(b.bsym, [ctranspose(i) for i in b.states], b.statemap)
ctranspose(b::TensorBasis) = TensorBasis([ctranspose(i) for i in b.bases], [ctranspose(i) for i in b.states], b.statemap)

for t=(:TensorBasis,:Basis)
	@eval begin
	isdual{T}(a::($t){Ket{T}}, b::($t){Bra{T}})= a.statemap==b.statemap && bsym(a)==bsym(b)
	isdual(a::($t){Ket}, b::($t){Bra})= a.statemap==b.statemap && bsym(a)==bsym(b)
	isdual{B<:Bra, K<:Ket}(a::($t){B}, b::($t){K}) = isdual(b,a)
	end
end
isdual(a::AbstractBasis,b::AbstractBasis)=false 

getindex(b::AbstractBasis, x) = b.states[x]

for t=(:Ket,:Bra)
	@eval begin
	get{S1<:($t),S2<:($t)}(b::Basis{S1}, s::S2, notfound) = samebasis(b,s) ? get(b.statemap, label(s), notfound) : notfound
	get{S1<:($t),S2<:($t)}(b::Basis{S1}, s::S2) = samebasis(b,s) ? b.statemap[label(s)] : throw(KeyError(s))
	get{S1<:($t),S2<:($t)}(b::TensorBasis{S1}, s::Tensor{S2}, notfound) = samebasis(b,s) ? get(b.statemap, label(s), notfound) : notfound
	get{S1<:($t),S2<:($t)}(b::TensorBasis{S1}, s::Tensor{S2}) = samebasis(b,s) ? b.statemap[label(s)] : throw(KeyError(s))
	in{S1<:($t),S2<:($t)}(s::S1, b::Basis{S2})= get(b,s,"FALSE")=="FALSE" ? false : true
	in{S1<:($t),S2<:($t)}(s::Tensor{S1}, b::TensorBasis{S2})= get(b,s,"FALSE")=="FALSE" ? false : true
	end
end

get(b::AbstractBasis, s::State, notfound) = notfound
get(b::AbstractBasis, s::State) = throw(KeyError(s))	
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

showcompact(io::IO, b::AbstractBasis) = print(io, "$(typeof(b)) $(reprlabel(b))")

function show(io::IO, b::AbstractBasis)
	showcompact(io, b)
	println(", $(length(b)) states:")
	if length(b)>20
		for i=1:10
			println(io, b.states[i])
		end
		println(vdots)
		for i=length(b)-10:length(b)
			println(io, b.states[i])
		end
	else
		for i in b.states
			println(io, i)
		end	
	end
end

######################################
##Function-Passing Functions##########
######################################

find(f::Function, b::Basis) = find(f, b.states)
filter(f::Function, b::Basis) = consbasis(filter(f, b.states))
filter(f::Function, b::TensorBasis) = consbasis(b.bases, filter(f, b.states))
map(f::Function, b::AbstractBasis) = basis(map(f, b.states))

######################################
##Joining/Separating Functions########
######################################

function stateappend(d::Dict,s::State)
	resmap = copy(d)
	resmap[label(s)] = length(d)+1
	return resmap
end

for t=(:Ket,:Bra)
	@eval begin
	function bjoin{S1<:($t), S2<:($t)}(b::Basis{S1}, s::S2)
		if in(s, b)
			return b
		else
			@assert samebasis(b, s) "BasisMismatch"
			return Basis(b.bsym, vcat(b.states, s), stateappend(b.statemap,s))
		end
	end

	function bjoin{S1<:($t), S2<:($t)}(s::S2,b::Basis{S1})
		if in(s, b)
			return b
		else
			@assert samebasis(b, s) "BasisMismatch"
			return consbasis(vcat(s,b.states))
		end
	end

	function bjoin{S1<:($t), S2<:($t)}(a::Basis{S1}, b::Basis{S2}) 
		@assert samebasis(a,b) "BasisMismatch"
		return consbasis(vcat(a.states, b.states))
	end
	
	function bjoin{S1<:($t), S2<:($t)}(b::TensorBasis{S1}, s::Tensor{S2})
		if in(s, b)
			return b
		else
			@assert samebasis(b, s) "BasisMismatch"
			return TensorBasis(b.bases, vcat(b.states, s), stateappend(b.statemap,s))
		end
	end

	function bjoin{S1<:($t), S2<:($t)}(s::Tensor{S2},b::TensorBasis{S1})
		if in(s, b)
			return b
		else
			@assert samebasis(b, s) "BasisMismatch"
			return consbasis(b.bases, vcat(s,b.states))
		end
	end

	function bjoin{S1<:($t), S2<:($t)}(a::TensorBasis{S1}, b::TensorBasis{S2}) 
		@assert samebasis(a,b) "BasisMismatch"
		return consbasis(b.bases, vcat(a.states, b.states))
	end

	tensor{S1<:($t), S2<:($t)}(s::State{S1},b::AbstractBasis{S2}) = map(s->tensor(s,b), b)
	tensor{S1<:($t), S2<:($t)}(b::AbstractBasis{S1}, s::State{S2}) = map(s->tensor(b,s), b)
	end
end

*(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
*(a::AbstractBasis, b::State) = tensor(a,b)
*(a::State, b::AbstractBasis) = tensor(a,b)
+(a::AbstractBasis, b::AbstractBasis) = bjoin(a,b)
+(a::AbstractBasis, b::State) = bjoin(a,b)
+(a::State, b::AbstractBasis) = bjoin(a,b)

setdiff(a::AbstractBasis,b::AbstractBasis) = setdiff(a.states, b.states)
