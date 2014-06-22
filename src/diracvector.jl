#####################################
#DiracVector#########################
#####################################

type DiracVector{S<:Single, T} <: Dirac
	coeffs::SparseMatrixCSC{T} 
	basis::AbstractBasis{S}
	function DiracVector(coeffs, basis)
		if S<:Ket
			@assert size(coeffs)==(length(basis),1) "Dimensions of coefficient array do not match basis"
			return new(coeffs, basis)
		else
			@assert size(coeffs)==(1,length(basis)) "Dimensions of coefficient array do not match basis"
			return new(coeffs, basis)
		end
	end
end

DiracVector{S,T}(coeffs::Array{T}, basis::AbstractBasis{S}) = DiracVector{S,T}(sparse(coeffs), basis)
DiracVector{S,T}(coeffs::SparseMatrixCSC{T}, basis::AbstractBasis{S}) = DiracVector{S,T}(coeffs, basis)

#####################################
#Getter-style Functions##############
#####################################
bsym(d::DiracVector) = bsym(d.basis)

#####################################
#Boolean Functions###################
#####################################
isequal(a::DiracVector, b::DiracVector) = isequal(a.coeffs, b.coeffs) && a.basis==b.basis
==(a::DiracVector, b::DiracVector) = a.coeffs==b.coeffs && a.basis==b.basis
isdual{K<:Ket,B<:Bra}(a::DiracVector{K}, b::DiracVector{B}) = isdual(a.basis, b.basis) && a.coeffs==b.coeffs'
isdual{K<:Ket,B<:Bra}(a::DiracVector{B}, b::DiracVector{K}) = isdual(b,a)
isdual(a::DiracVector, b::DiracVector) = false

#####################################
#Array-like Functions################
#####################################
copy(d::DiracVector) = DiracVector(copy(d.coeffs), copy(d.basis))

ctranspose(d::DiracVector) = DiracVector(d.coeffs', d.basis')
size(d::DiracVector, args...) = size(d.coeffs, args...)
getindex{K<:Ket}(d::DiracVector{K}, x) = d.coeffs[x,1]
getindex{B<:Bra}(d::DiracVector{B}, x) = d.coeffs[1,x]
getindex(d::DiracVector, x, y) = d.coeffs[x,y]
setindex!{K<:Ket}(d::DiracVector{K}, y, x) = setindex!(d.coeffs, y, x, 1)
setindex!{B<:Bra}(d::DiracVector{B}, y, x) = setindex!(d.coeffs, y, 1, x)
setindex!(d::DiracVector, y, x...) = setindex!(d.coeffs, y, x...)

for op=(:endof, :ndims, :eltype, :length, :find, :findn, :findnz, :nnz)
	@eval ($op)(d::DiracVector) = ($op)(d.coeffs)
end

#####################################
#Dict-like Functions#################
#####################################
getpos(d::DiracVector, s::State) = get(d.basis, s)
function get(d::DiracVector, s::State, notfound)
	try
		return d[getpos(d, s)]
	catch
		return notfound
	end
end

get(d::DiracVector, s::State) = d[getpos(d, s)]
in(s::State, d::DiracVector) = in(s, d.basis)

#####################################
#Show Functions######################
#####################################

summary{K<:Ket}(d::DiracVector{K}) = "$(size(d,1))x1 $(typeof(d))"
summary{B<:Bra}(d::DiracVector{B}) = "1x$(size(d,2)) $(typeof(d))"

function showcompact(io::IO, d::DiracVector)
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else
		tempio = IOBuffer()
		print(tempio, [" + $(d.coeffs[i])$(d.basis[i])" for i in find(d)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end
function show{S}(io::IO, d::DiracVector{S})
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else	
		println(io, summary(d))
		table = cell(length(d), 2)	
		if length(d)>=52
			for i=1:25
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
			table[26:(length(d)-25),:] = 0 # prevents access to undefined reference
			for i=(length(d)-25):length(d)
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
		else
			for i=1:length(d)
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
		end
		temp_io = IOBuffer()
		if S<:Ket
			show(temp_io, table)
		else
			show(temp_io, [transpose(table[:,2]), transpose(table[:,1])])
		end
		io_str = takebuf_string(temp_io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		print(io, io_str)
	end
end

#####################################
#Function-passing Functions##########
#####################################

find(f::Function, d::DiracVector) = find(f, d.coeffs)
findstates(f::Function, d::DiracVector) = find(f, d.basis)

map(f::Function, d::DiracVector)=DiracVector(map(f, d.coeffs), d.basis)

function map!(f::Function, d::DiracVector)
	map!(f, d.coeffs)
	return d
end

function loadcoeffs!(coeffs, d, newbasis)
	for i=1:length(newbasis)
		coeffs[i] = get(d, newbasis[i])
	end	
	return coeffs
end

function filtercons{S,T}(newbasis::AbstractBasis, d::DiracVector{S,T})
	return DiracVector(loadcoeffs!(S<:Ket ? Array(T, length(newbasis)) : Array(T, 1, length(newbasis)), d, newbasis), newbasis)
end

filterstates{S,T}(f::Function, d::DiracVector{S,T}) = filtercons(filter(f, d.basis), d)
filtercoeffs(f::Function, d::DiracVector) = filtercons(basis(d.basis[find(map(f, d.coeffs))]), d)

# qeval(f::Function, d::DiracVector) = map(x->qeval(f, x), d)

#####################################
#Arithmetic Operations###############
#####################################

for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracVector, b::DiracVector) = DiracVector(($op)(a.coeffs,b.coeffs), a.basis)
	@eval ($op)(n, d::DiracVector) = DiracVector(($op)(n,d.coeffs), d.basis)
	@eval ($op)(d::DiracVector, n) = DiracVector(($op)(d.coeffs,n), d.basis)
end

*(c::DiracCoeff, d::DiracVector) = c.*d
*(d::DiracVector, c::DiracCoeff) = c*d

ireduce(s::State, d::DiracVector) = reduce(+,[d[i]*inner(s,d.basis[i]) for i=1:length(d)])
ireduce(d::DiracVector, s::State) = reduce(+,[d[i]*inner(d.basis[i],s) for i=1:length(d)])
ireduce(a::DiracVector, b::DiracVector) = reduce(+,[a[i]*b[j]*inner(a.basis[i],b.basis[j]) for i=1:length(a), j=1:length(b)])
ireduce(a::DiracVector, b::DiracVector, eqbasis::Bool) = reduce(+,[a[i]*b[i]*inner(a.basis[i],b.basis[i]) for i=1:length(a)])

inner{B<:Bra, K<:Ket}(s::State{B}, d::DiracVector{K}) = samebasis(d,s)	? get(d, s', 0) : ireduce(s,d)
inner{B<:Bra, K<:Ket}(d::DiracVector{B}, s::State{K}) = samebasis(d,s)	? get(d, s', 0) : ireduce(d,s)
inner{B<:Bra, K<:Ket, N1<:Number, N2<:Number}(a::DiracVector{B, N1}, b::DiracVector{K, N2}) = isdual(a.basis, b.basis) ? (a.coeffs*b.coeffs)[1] : ireduce(a,b)
inner{B<:Bra, K<:Ket}(a::DiracVector{B}, b::DiracVector{K}) = isdual(a.basis, b.basis) ? ireduce(a,b,true) : ireduce(a,b)

# *(a::DiracVector{Bra}, b::DiracVector{Ket}) = inner(a,b)
# *(a::DiracVector{Bra}, b::AbstractState{Ket}) = inner(a,b)
# *(a::AbstractState{Bra}, b::DiracVector{Ket}) = inner(a,b)
# #see kron in misc.jl
# *{K}(a::DiracVector{K}, b::DiracVector{K}) = kron(a,b)
# *{K}(a::AbstractState{K}, b::DiracVector{K}) = kron(a,b)
# *{K}(a::DiracVector{K}, b::AbstractState{K}) = kron(a,b)
# *(a::DiracVector{Ket}, b::AbstractState{Bra}) = kron(a,b)
# *(a::AbstractState{Ket}, b::DiracVector{Bra}) = kron(a,b)
# *(a::DiracVector{Ket}, b::DiracVector{Bra}) = kron(a,b)


# function +{K}(d::DiracVector{K}, s::AbstractState{K})
# 	if in(s, d.basis)
# 		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
# 		res[getpos(d,s)] = 1+get(res, s)
# 		return res
# 	else
# 		@assert samebasis(d,s) "BasisMismatch"
# 		if K==Ket
# 			return DiracVector(vcat(d.coeffs, 1), basisjoin(d.basis,s))
# 		else
# 			return DiracVector(hcat(d.coeffs, 1), basisjoin(d.basis,s))
# 		end
# 	end
# end

# function +{K}(s::AbstractState{K}, d::DiracVector{K})
# 	if in(s, d.basis)
# 		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
# 		res[getpos(d,s)] = 1+get(res, s)
# 		return res
# 	else
# 		@assert samebasis(d,s) "BasisMismatch"
# 		if K==Ket
# 			return DiracVector(vcat(1, d.coeffs), basisjoin(s, d.basis))
# 		else
# 			return DiracVector(hcat(1, d.coeffs), basisjoin(s, d.basis))
# 		end
# 	end
# end

# function +{K}(a::DiracVector{K}, b::DiracVector{K})
# 	if a.basis==b.basis
# 		return DiracVector(a.coeffs+b.coeffs, a.basis)
# 	else
# 		@assert samebasis(a,b) "BasisMismatch"
# 		res = 1*a
# 		for i=1:length(b)
# 			res = res+b.basis[i]
# 			res[getpos(res, b.basis[i])] = get(res, b.basis[i]) + b[i] - 1
# 		end
# 		return res
# 	end
# end

# -{K}(d::DiracVector{K}, s::AbstractState{K}) = d+(-s)
# -{K}(s::AbstractState{K}, d::DiracVector{K}) = s+(-d)
# -{K}(a::DiracVector{K}, b::DiracVector{K}) = a+(-b)

# +(a::AbstractState{Ket}, b::AbstractState{Ket}) = a==b ? 2*a : DiracVector([1, 1], tobasis([a,b])) 
# +(a::AbstractState{Bra}, b::AbstractState{Bra}) = a==b ? 2*a : DiracVector([1 1], tobasis([a,b])) 
# -(a::AbstractState{Ket}, b::AbstractState{Ket}) = a==b ? 0 : DiracVector([1, -1], tobasis([a,b])) 
# -(a::AbstractState{Bra}, b::AbstractState{Bra}) =  a==b ? 0 : DiracVector([1 -1], tobasis([a,b])) 

# -(s::AbstractState) = -1*s
# -(d::DiracVector) = DiracVector(-1*d.coeffs, d.basis)

# norm(d::DiracVector, p::Int=2) = reduce(+,map(i->abs(i)^p, d.coeffs))^(1/p) 
# normalize(d::DiracVector) = DiracVector((1/norm(d.coeffs))*d.coeffs, d.basis)
