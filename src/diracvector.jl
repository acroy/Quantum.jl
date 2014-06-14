#####################################
#DiracVector#########################
#####################################

type DiracVector{K<:BraKet, T} <: Dirac
	coeffs::SparseMatrixCSC{T} 
	basis::AbstractBasis{K}
	function DiracVector(coeffs, basis)
		if K==Ket
			@assert size(coeffs)==(length(basis),1) "Dimensions of coefficient array do not match basis"
				new(coeffs, basis)
		else
			@assert size(coeffs)==(1,length(basis)) "Dimensions of coefficient array do not match basis"
				new(coeffs, basis)
		end
	end
end

DiracVector{K,T}(coeffs::Array{T}, basis::AbstractBasis{K}) = DiracVector{K,T}(sparse(coeffs), basis)
DiracVector{K,T}(coeffs::SparseMatrixCSC{T}, basis::AbstractBasis{K}) = DiracVector{K,T}(coeffs, basis)

#####################################
#Getter-style Functions##############
#####################################
basislabel(d::DiracVector) = label(d.basis)
kind(d::DiracVector) = kind(d.basis)

#####################################
#Boolean Functions###################
#####################################
isequal(a::DiracVector, b::DiracVector) = isequal(a.coeffs, b.coeffs) && a.basis==b.basis
==(a::DiracVector, b::DiracVector) = a.coeffs==b.coeffs && a.basis==b.basis
isdual(a::DiracVector{Ket}, b::DiracVector{Bra}) = isdual(a.basis, b.basis) && a.coeffs==b.coeffs'
isdual(a::DiracVector{Bra}, b::DiracVector{Ket}) = isdual(b,a)
isdual{K}(a::DiracVector{K}, b::DiracVector{K}) = false

#####################################
#Array-like Functions################
#####################################
copy(d::DiracVector) = DiracVector(copy(d.coeffs), copy(d.basis))

ctranspose(d::DiracVector) = DiracVector(d.coeffs', d.basis')
size(d::DiracVector, args...) = size(d.coeffs, args...)
getindex(d::DiracVector{Ket}, x) = d.coeffs[x,1]
getindex(d::DiracVector{Bra}, x) = d.coeffs[1,x]
getindex(d::DiracVector, x, y) = d.coeffs[x,y]
setindex!(d::DiracVector{Ket}, y, x) = setindex!(d.coeffs, y, x, 1)
setindex!(d::DiracVector{Bra}, y, x) = setindex!(d.coeffs, y, 1, x)
setindex!(d::DiracVector, y, x...) = setindex!(d.coeffs, y, x...)

for op=(:endof, :ndims, :eltype, :length, :find, :findn, :findnz, :nnz)
	@eval ($op)(d::DiracVector) = ($op)(d.coeffs)
end

#####################################
#Dict-like Functions#################
#####################################
getpos(d::DiracVector, s::AbstractState) = get(d.basis, s)
function get(d::DiracVector, s::AbstractState, notfound)
	try
		return d[getpos(d, s)]
	catch
		return notfound
	end
end

get(d::DiracVector, s::AbstractState) = d[getpos(d, s)]
in(s::AbstractState, d::DiracVector) = in(s, d.basis)
#####################################
#Show Functions######################
#####################################

summary(d::DiracVector{Ket}) = "$(size(d,1))x1 $(typeof(d))"
summary(d::DiracVector{Bra}) = "1x$(size(d,2)) $(typeof(d))"

function showcompact(io::IO, d::DiracVector)
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else
		tempio = IOBuffer()
		print(tempio, [" + $(d.coeffs[i])$(d.basis[i])" for i in find(d)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end
function show(io::IO, d::DiracVector)
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
		if kind(d)==Ket
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
	d.coeffs = map(f, d.coeffs)
	return d
end

function mapmatch(fstates::Function, fcoeffs::Function, d::DiracVector)
	matched = findstates(fstates, d)	
	coeffs = convert(SparseMatrixCSC{Any}, d.coeffs)
	for i in matched
		coeffs[i] = fcoeffs(coeffs[i])
	end
	#hacky concat strategy to force corect typing
	if kind(d)==Ket
		return DiracVector(vcat(coeffs...), d.basis)
	else
		return DiracVector(hcat(coeffs...), d.basis)		
	end
end

function mapmatch!(fstates::Function, fcoeffs::Function, d::DiracVector)
	matched = findstates(fstates, d)	
	for i in matched
		d[i] = fcoeffs(d[i])
	end
	return d
end

function filtercons{K,T}(newbasis::AbstractBasis, d::DiracVector{K,T})
	coeffs = K==Ket ? Array(T, length(newbasis)) : Array(T, 1, length(newbasis))
	for i=1:length(newbasis)
		coeffs[i] = get(d, newbasis[i])
	end	
	return DiracVector(coeffs, newbasis)
end

filterstates{K,T}(f::Function, d::DiracVector{K,T}) = filtercons(filter(f, d.basis), d)
filtercoeffs(f::Function, d::DiracVector) = filtercons(Basis(d.basis[find(map(f, d.coeffs))]), d)

qeval(f::Function, d::DiracVector) = map(x->qeval(f, x), d)

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

function inner(s::AbstractState{Bra}, d::DiracVector{Ket})
	if samebasis(d,s)	
		return get(d, s', 0)
	else
		return reduce(+,[d[i]*(s*d.basis[i]) for i=1:length(d)])
	end
end
function inner(d::DiracVector{Bra}, s::AbstractState{Ket})
	if samebasis(d,s)	
		return get(d, s', 0)
	else
		return reduce(+,[d[i]*(d.basis[i]*s) for i=1:length(d)])
	end
end

function inner{N1<:Number, N2<:Number}(a::DiracVector{Bra, N1}, b::DiracVector{Ket, N2})
	if a.basis==b.basis 
		return (a.coeffs*b.coeffs)[1]
	else
		return reduce(+,[a[i]*b[j]*(a.basis[i]*b.basis[j]) for i=1:length(a), j=1:length(b)])
	end
end

function inner(a::DiracVector{Bra}, b::DiracVector{Ket})
	if a.basis==b.basis 
		return reduce(+,[a[i]*b[i]*(a.basis[i]*b.basis[i]) for i=1:length(a)]) 
	else
		return reduce(+,[a[i]*b[j]*(a.basis[i]*b.basis[j]) for i=1:length(a), j=1:length(b)])
	end
end	

*(a::DiracVector{Bra}, b::DiracVector{Ket}) = inner(a,b)
*(a::DiracVector{Bra}, b::AbstractState{Ket}) = inner(a,b)
*(a::AbstractState{Bra}, b::DiracVector{Ket}) = inner(a,b)
#see kron in misc.jl
*{K}(a::DiracVector{K}, b::DiracVector{K}) = kron(a,b)
*{K}(a::AbstractState{K}, b::DiracVector{K}) = kron(a,b)
*{K}(a::DiracVector{K}, b::AbstractState{K}) = kron(a,b)
*(a::DiracVector{Ket}, b::AbstractState{Bra}) = kron(a,b)
*(a::AbstractState{Ket}, b::DiracVector{Bra}) = kron(a,b)
*(a::DiracVector{Ket}, b::DiracVector{Bra}) = kron(a,b)


function +{K}(d::DiracVector{K}, s::AbstractState{K})
	if in(s, d.basis)
		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
		res[getpos(d,s)] = 1+get(res, s)
		return res
	else
		@assert samebasis(d,s) "BasisMismatch"
		if K==Ket
			return DiracVector(vcat(d.coeffs, 1), basisjoin(d.basis,s))
		else
			return DiracVector(hcat(d.coeffs, 1), basisjoin(d.basis,s))
		end
	end
end

function +{K}(s::AbstractState{K}, d::DiracVector{K})
	if in(s, d.basis)
		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
		res[getpos(d,s)] = 1+get(res, s)
		return res
	else
		@assert samebasis(d,s) "BasisMismatch"
		if K==Ket
			return DiracVector(vcat(1, d.coeffs), basisjoin(s, d.basis))
		else
			return DiracVector(hcat(1, d.coeffs), basisjoin(s, d.basis))
		end
	end
end

function +{K}(a::DiracVector{K}, b::DiracVector{K})
	if a.basis==b.basis
		return DiracVector(a.coeffs+b.coeffs, a.basis)
	else
		@assert samebasis(a,b) "BasisMismatch"
		res = 1*a
		for i=1:length(b)
			res = res+b.basis[i]
			res[getpos(res, b.basis[i])] = get(res, b.basis[i]) + b[i] - 1
		end
		return res
	end
end

-{K}(d::DiracVector{K}, s::AbstractState{K}) = d+(-s)
-{K}(s::AbstractState{K}, d::DiracVector{K}) = s+(-d)
-{K}(a::DiracVector{K}, b::DiracVector{K}) = a+(-b)

+(a::AbstractState{Ket}, b::AbstractState{Ket}) = a==b ? 2*a : DiracVector([1, 1], tobasis([a,b])) 
+(a::AbstractState{Bra}, b::AbstractState{Bra}) = a==b ? 2*a : DiracVector([1 1], tobasis([a,b])) 
-(a::AbstractState{Ket}, b::AbstractState{Ket}) = a==b ? 0 : DiracVector([1, -1], tobasis([a,b])) 
-(a::AbstractState{Bra}, b::AbstractState{Bra}) =  a==b ? 0 : DiracVector([1 -1], tobasis([a,b])) 

-(s::AbstractState) = -1*s
-(d::DiracVector) = DiracVector(-1*d.coeffs, d.basis)

norm(d::DiracVector, p::Int=2) = reduce(+,map(i->abs(i)^p, d.coeffs))^(1/p) 
normalize(d::DiracVector) = DiracVector((1/norm(d.coeffs))*d.coeffs, d.basis)
