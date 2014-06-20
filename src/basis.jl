#####################################
#Basis###############################
#####################################

abstract AbstractBasis{S<:Single} <: Dirac

immutable Basis{S<:Single} <: AbstractBasis{S}
	bsym::Symbol
	states::Vector{S}
	statemap::Dict{S,Int}
end

statemapper{S<:State}(sv::Vector{S}) = (S=>Int)[sv[i]=>i for i=1:length(sv)]

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
	statemap::Dict{Tensor{S}, Int}
end

consbasis{S<:Single,B<:Basis}(bv::Vector{B}, sv::Vector{Tensor{S}}) = TensorBasis(bv, sv, statemapper(sv))

separate(b::Basis) = [b]
separate(b::TensorBasis) = b.bases

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

# #improved type inferencing for statejoin; 
# #based on the idea that each column has its
# #own type due to the behavior of crossjoin
sjointype(sarr::Array)=typejoin(map(eltype, sarr[1,:])...)

statejoin{S<:Single}(sarr::Array{S,2}) = Tensor{S}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]
statejoin{S<:State}(sarr::Array{S,2}) = Tensor{sjointype(sarr)}[tensor(sarr[i, :]) for i=1:size(sarr, 1)]

statecross(v::Vector) = statejoin(reduce(crossjoin, v))

tensor() = error("tensor needs arguments of type T<:Basis or T<:State")
tensor(b::AbstractBasis) = error("cannot perform tensor operation on one basis")

function tensor{S<:Single}(b::AbstractBasis{S}...)
	sv = statecross([i.states for i in b])
	TensorBasis([[separate(i) for i in b]...], sv, statemapper(sv))
end

#####################################
#Misc Functions######################
#####################################

# copy{B<:Basis{Single}}(b::B) = B(copy(b.bsym), copy(b.states), copy(b.statemap))
# copy{B<:Basis{Tensor}}(b::B)= B(copy(b.bases), copy(b.states), copy(b.statemap))

# bsym(b::Basis{Single}) = b.bsym
# bsym(b::Basis{Tensor}) = map(bsym, b.bases)

# isequal{B<:Basis}(a::B, b::B) = isequal(a.states, b.states) && bsym(a)==bsym(b)
# =={B<:Basis}(a::B, b::B) = a.states==b.states && bsym(a)==bsym(b)

# for op=(:length, :endof)
# 	@eval ($op)(b::Basis) = ($op)(b.states)
# end

# dual(t::Type{KetBasis}) = BraBasis
# dual(t::Type{BraBasis}) = KetBasis
# dual{K}(t::Type{KetBasis{K}}) = BraBasis{dual(K)}
# dual{B}(t::Type{BraBasis{B}}) = KetBasis{dual(B)}
# dual(t::Type{TensorKetBasis}) = TensorBraBasis
# dual(t::Type{TensorBraBasis}) = TensorKetBasis
# dual{K}(t::Type{TensorKetBasis{K}}) = TensorBraBasis{dual(K)}
# dual{B}(t::Type{TensorBraBasis{B}}) = TensorKetBasis{dual(B)}
# dual{K,KB}(t::Type{TensorKetBasis{K,KB}}) = TensorBraBasis{dual(K), dual(KB)}
# dual{B,BB}(t::Type{TensorBraBasis{B,BB}}) = TensorKetBasis{dual(B), dual(BB)}

# ctranspose{B<:Basis{Single}}(b::B) = consbasis([ctranspose(i) for i in b.states])
# ctranspose{B<:Basis{Tensor}}(b::B) = consbasis([ctranspose(i) for i in b.bases], [ctranspose(i) for i in b.states])

# isdual{T}(a::KetBasis{Ket{T}}, b::BraBasis{Bra{T}})= length(a)==length(b) ? all(map(isdual, a.states, b.states)) : false
# isdual(a::KetBasis{Ket}, b::BraBasis{Bra}) = length(a)==length(b) ? all(map(isdual, a.states, b.states)) : false
# isdual(a::BraBasis, b::KetBasis) = isdual(b,a)
# isdual{T}(a::TensorKetBasis{Ket{T},KetBasis{Ket{T}}}, 
# 		  b::TensorBraBasis{Bra{T},BraBasis{Bra{T}}}) = length(a)==length(b) ? all(map(isdual, a.states, b.states)) : false
# isdual(a::TensorKetBasis{Ket,KetBasis}, b::TensorBraBasis{Bra,BraBasis}) = length(a)==length(b) ? all(map(isdual, a.states, b.states)) : false
# isdual(a::TensorBraBasis, b::TensorKetBasis) = isdual(b,a)
# isdual(a::Basis,b::Basis)=false 

# size(b::TensorBasis) = (length(b.states), length(b.bases))

# getindex(b::AbstractBasis, x) = b.states[x]
# setindex!(b::AbstractBasis, y, x) = setindex!(b.states, y, x)

# get{K,T}(b::Basis{K,T}, s::State{K,T}, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
# get{K}(b::TensorBasis{K}, s::TensorState{K}, notfound) = get(b.statemap, (label(s), basislabel(s)), notfound)
# get(b::AbstractBasis, s::AbstractState, notfound) = notfound
# get{K,T}(b::Basis{K,T}, s::State{K,T}) = b.statemap[(label(s), basislabel(s))]
# get{K}(b::TensorBasis{K}, s::TensorState{K}) = b.statemap[(label(s), basislabel(s))]
# get(b::AbstractBasis, s::AbstractState) = throw(KeyError(s))

# in(s::AbstractState, b::AbstractBasis)= get(b,s,"FALSE")=="FALSE" ? false : true

######################################
##Show Functions######################
######################################

# reprlabel(b::Basis{Single}) = bsym(b)
# function reprlabel(b::Basis{Tensor})
# 	labels = bsym(b)
# 	#terrible way to grow a string
# 	str = "$(labels[1])"
# 	for i=2:length(labels)
# 		str = "$str$otimes$(labels[i])"
# 	end
# 	return str
# end

# showcompact(io::IO, b::Basis) = print(io, "$(typeof(b)) $(reprlabel(b))")

# function show(io::IO, b::Basis)
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
