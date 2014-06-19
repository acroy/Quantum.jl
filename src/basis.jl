#####################################
#Basis###############################
#####################################

abstract Basis{S<:Union(Single, Tensor)} <: Dirac
abstract AbstractKetBasis{S<:Union(Single, Tensor)} <: Basis{S}
abstract AbstractBraBasis{S<:Union(Single, Tensor)} <: Basis{S}

immutable KetBasis{K<:Ket} <: AbstractKetBasis{Single}
	label::Symbol
	states::Vector{K}
	statemap::Dict{K,Int}
end

immutable BraBasis{B<:Bra} <: AbstractBraBasis{Single}
	label::Symbol
	states::Vector{B}
	statemap::Dict{B,Int}
end

#makebasis is for internal use only; 
#useful for when we KNOW that the input vector of
#states contains only unique elements

statemapper{S<:State}(sv::Vector{S}) = (S=>Int)[sv[i]=>i for i=1:length(sv)]

function makebasis{K<:Ket}(sv::Vector{K})
	@assert length(unique(map(bsym, sv)))==1 "BasisMismatch"
	KetBasis(bsym(sv[1]), sv, statemapper(sv))
end

function makebasis{B<:Bra}(sv::Vector{B})
	@assert length(unique(map(bsym, sv)))==1 "BasisMismatch"
	BraBasis(bsym(sv[1]), sv, statemapper(sv))
end

basis{S<:State{Single}}(states::Vector{S}) = makebasis(unique(states))
basis{S<:State{Single}}(states::Array{S}) = basis(vec(states))
basis(s::State...) = basis(collect(s))

#####################################
#TensorBasis#########################
#####################################

immutable TensorKetBasis{K<:Ket, KB<:KetBasis} <: AbstractKetBasis{Tensor}
	bases::Vector{KB}
	states::Vector{TensorKet{K}}
	statemap::Dict{TensorKet{K}, Int}
end

immutable TensorBraBasis{B<:Bra, BB<:BraBasis} <: AbstractBraBasis{Tensor}
	bases::Vector{BB}
	states::Vector{TensorBra{B}}
	statemap::Dict{TensorBra{B}, Int}
end

separate(b::Basis{Single}) = [b]
separate(b::Basis{Tensor}) = b.bases

sepstates(sv) = hcat(map(separate, sv)...)

function consbvec!(bases, seps::Array)
	for i=1:size(seps, 1)
	 	bases[i] = basis(seps[i, :])
	end
	bases
end

function consbvec!(bases, seps::Vector)
	for i=1:length(seps)
	 	bases[i] = seps[i]
	end
	bases
end

prepbvec{K<:TensorKet}(sv::Vector{K}) = Array(KetBasis{eltype(eltype(sv))}, length(sv[1]))
prepbvec{B<:TensorBra}(sv::Vector{B}) = Array(BraBasis{eltype(eltype(sv))}, length(sv[1]))

maketensorbasis{K<:TensorKet}(sv::Vector{K}) = TensorKetBasis(consbvec!(prepbvec(sv),sepstates(sv)), sv, statemapper(sv))
maketensorbasis{B<:TensorBra}(sv::Vector{B}) = TensorBraBasis(consbvec!(prepbvec(sv),sepstates(sv)), sv, statemapper(sv))

basis{S<:State{Tensor}}(states::Vector{S}) = maketensorbasis(unique(states))
basis{S<:State{Tensor}}(states::Array{S}) = basis(vec(states))

statecross(v::Vector) = statejoin(reduce(crossjoin, v))

tensor() = error("tensor needs arguments of type T<:Basis or T<:State")
tensor(b::Basis) = error("cannot perform tensor operation on one basis")

function tensor(b::AbstractKetBasis...)
	sv = statecross([i.states for i in b])
	TensorKetBasis([[separate(i) for i in b]...], sv, statemapper(sv))
end

function tensor(b::AbstractBraBasis...)
	sv = statecross([i.states for i in b])
	TensorBraBasis([[separate(i) for i in b]...], sv, statemapper(sv))
end

# statecross{K<:BraKet,T}(bases::Vector{Basis{K,T}}) = statejoin(crossjoin([i.states for i in bases]...))

# statecross{K<:BraKet}(bases::Vector{Basis{K}}) = statejoin(crossjoin([i.states for i in bases]...))
# statecross(bases::Vector{Basis}) = statecross(Basis{kind(bases[1])}[bases])

# statemapper{K<:BraKet,T}(sv::Vector{TensorState{K,State{K,T}}}) = ((Vector{T}, Vector{Symbol})=>Int)[(label(sv[i]), basislabel(sv[i]))=>i for i=1:length(sv)]

# function btensor{K}(bases::AbstractBasis{K}...) #utilizes memoization of already tensored states
# 	sv = statecross(bases)
# 	TensorBasis{K}(vcat([separate(i) for i in bases]...), sv, statemapper(sv))
# end

# btensor{K}(basis::AbstractBasis{K}) = basis
# TensorBasis{K}(bases::AbstractBasis{K}...) = btensor(bases...)

# #makebasis but for TensorBasis
# function maketensorbasis{K}(sv::Vector{TensorState{K}})
# 	sepstates = hcat(map(separate, sv)...)
# 	bases = Array(Basis{K}, size(sepstates,1))
# 	for i=1:size(sepstates, 1)
# 	 	bases[i] = Basis([sepstates[i, :]...]) #vcat trick to force correct typing of subarray
# 	end
# 	return TensorBasis{K}(bases, sv, statemapper(sv))
# end

# #in case Bases are already explicitly known;
# #once again, no checking is performed 
# #since this function is for internal use only
# maketensorbasis{B<:Basis,K}(b::Vector{B}, sv::Vector{TensorState{K}})=TensorBasis{K}(b, sv, statemapper(sv))
# TensorBasis{K}(sv::Vector{TensorState{K}}) = maketensorbasis(unique(sv))

# #####################################
# #Misc Functions######################
# #####################################

# separate(b::Basis)=[b]
# separate(b::TensorBasis) = b.bases

# copy{K,T}(b::Basis{K,T}) = Basis{K,T}(copy(b.label), copy(b.states), copy(b.statemap))
# copy{K}(b::TensorBasis{K}) = TensorBasis{K}(copy(b.bases), copy(b.states), copy(b.statemap))

# kind{K}(b::AbstractBasis{K}) = K

# label(b::Basis) = b.label
# label(b::TensorBasis) = map(label, b.bases)
# basislabel(b::AbstractBasis) = label(b)

# isequal{K,T}(a::Basis{K,T}, b::Basis{K,T}) = isequal(a.states, b.states) && label(a)==label(b)
# =={K,T}(a::Basis{K,T}, b::Basis{K,T}) = a.states==b.states && label(a)==label(b)
# isequal{K}(a::TensorBasis{K}, b::TensorBasis{K}) = isequal(a.states, b.states) && label(a)==label(b)
# =={K}(a::TensorBasis{K}, b::TensorBasis{K}) = a.states==b.states && label(a)==label(b)
# isequal(a::AbstractBasis, b::AbstractBasis) = false
# ==(a::AbstractBasis, b::AbstractBasis) = false 

# ctranspose{K}(b::TensorBasis{K}) = TensorBasis{!K}(map(ctranspose, b.bases), map(ctranspose, b.states), b.statemap)
# ctranspose{K,T}(b::Basis{K,T}) = Basis{!K,T}(b.label, map(ctranspose, b.states), b.statemap)

# isdual{T}(a::Basis{Ket,T}, b::Basis{Bra,T}) = label(a)==label(b) && a.statemap==b.statemap
# isdual{T}(a::Basis{Bra,T}, b::Basis{Ket,T}) = isdual(b,a)
# isdual(a::TensorBasis{Ket}, b::TensorBasis{Bra}) = label(a)==label(b) && a.statemap==b.statemap
# isdual(a::TensorBasis{Bra}, b::TensorBasis{Ket}) = isdual(b,a)
# isdual(a::AbstractBasis,b::AbstractBasis)=false 

# size(b::TensorBasis) = (length(b.states), length(b.bases))
# length(b::AbstractBasis) = length(b.states)
# endof(b::AbstractBasis) = endof(b.states)

# getindex(b::AbstractBasis, x) = b.states[x]
# setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

# get{K,T}(b::Basis{K,T}, s::State{K,T}, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
# get{K}(b::TensorBasis{K}, s::TensorState{K}, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
# get(b::AbstractBasis, s::AbstractState, notfound) = notfound
# get{K,T}(b::Basis{K,T}, s::State{K,T}) = b.statemap[(label(s), basislabel(s))]
# get{K}(b::TensorBasis{K}, s::TensorState{K}) = b.statemap[(label(s), basislabel(s))]
# get(b::AbstractBasis, s::AbstractState) = throw(KeyError(s))

# in(s::AbstractState, b::AbstractBasis)= get(b,s,"FALSE")=="FALSE" ? false : true

# ######################################
# ##Show Functions######################
# ######################################

# reprlabel(b::Basis) = label(b)
# function reprlabel(b::TensorBasis)
# 	labels = label(b)
# 	#terrible way to grow a string
# 	str = "$(labels[1])"
# 	for i=2:length(labels)
# 		str = "$str$otimes$(labels[i])"
# 	end
# 	return str
# end

# showcompact(io::IO, b::AbstractBasis) = print(io, "$(typeof(b)) $(reprlabel(b))")

# function show(io::IO, b::AbstractBasis)
# 	showcompact(io, b)
# 	println(", $(length(b)) states:")
# 	if length(b)>20
# 		for i=1:10
# 			println(io, b.states[i])
# 		end
# 		println(vdots)
# 		for i=length(b)-10:length(b)
# 			println(io, b.states[i])
# 		end
# 	else
# 		for i in b.states
# 			println(io, i)
# 		end	
# 	end
# end
# ######################################
# ##Function-Passing Functions##########
# ######################################

# find(f::Function, b::Basis) = find(f, b.states)
# filter(f::Function, b::Basis) = makebasis(filter(f, b.states))
# filter(f::Function, b::TensorBasis) = maketensorbasis(b.bases, filter(f, b.states))

# function map(f::Function, b::AbstractBasis) 
# 	newstates = map(f, b.states)
# 	if eltype(newstates) <: TensorState
# 		return TensorBasis(newstates)
# 	else
# 		return Basis(newstates)
# 	end
# end

# ######################################
# ##Joining/Separating Functions########
# ######################################

# function basisjoin{K,T}(b::Basis{K,T}, s::State{K,T})
# 	if in(s, b)
# 		return b
# 	else
# 		@assert samebasis(b, s) "BasisMismatch"
# 		resmap = copy(b.statemap)
# 		resmap[(label(s), basislabel(s))] = length(b)+1
# 		return Basis{K,T}(b.label, vcat(b.states, s), resmap)
# 	end
# end

# function basisjoin{K,T}(s::State{K,T}, b::Basis{K,T})
# 	if in(s, b)
# 		return b
# 	else
# 		return makebasis(vcat(s,b.states))
# 	end
# end

# basisjoin{K,T}(a::Basis{K,T}, b::Basis{K,T}) = Basis(vcat(a.states, b.states))

# function basisjoin{K}(b::TensorBasis{K}, s::TensorState{K})
# 	if in(s, b)
# 		return b
# 	else
# 		@assert samebasis(b, s) "BasisMismatch"
# 		resmap = copy(b.statemap)
# 		resmap[(label(s), basislabel(s))] = length(b)+1
# 		return TensorBasis{K}(b.bases, vcat(b.states, s), resmap)
# 	end
# end

# function basisjoin{K}(s::TensorState{K}, b::TensorBasis{K})
# 	if in(s, b)
# 		return b
# 	else
# 		@assert samebasis(b, s) "BasisMismatch"
# 		return tensormakebasis(s.basis, vcat(s,b.states))
# 	end
# end

# basisjoin{K}(a::TensorBasis{K}, b::TensorBasis{K}) = TensorBasis(vcat(a.states, b.states))

# btensor{K}(a::AbstractBasis{K}, b::AbstractState{K}) = map(s->s*b, d.basis)
# btensor{K}(a::AbstractState{K}, b::AbstractBasis{K}) = map(s->b*s, d.basis)

# *{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = btensor(a,b)
# *{K}(a::AbstractBasis{K}, b::AbstractState{K}) = btensor(a,b)
# *{K}(a::AbstractState{K}, b::AbstractBasis{K}) = btensor(a,b)
# +{K}(a::AbstractBasis{K}, b::AbstractBasis{K}) = basisjoin(a,b)
# +{K}(a::AbstractBasis{K}, b::AbstractState{K}) = basisjoin(a,b)
# +{K}(a::AbstractState{K}, b::AbstractBasis{K}) = basisjoin(a,b)

# setdiff{B<:AbstractBasis}(a::B,b::B) = setdiff(a.states, b.states)

# tobasis(s::State) = Basis(s)
# tobasis(s::TensorState) = TensorBasis([s])
# tobasis(s::AbstractState...) = tobasis(collect(s))
# tobasis{S<:State}(v::Array{S}) = Basis(vec(v))
# tobasis{S<:TensorState}(v::Array{S}) = TensorBasis(vec(v))
