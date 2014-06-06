#####################################
#DiracVector#########################
#####################################

type DiracVector{T,K<:BraKet} <: Dirac
	coeffs::Array{T} 
	basis::AbstractBasis{K}
	function DiracVector(coeffs, basis)
		if samebasis("?", basis)
			error("BasisMismatch: cannot represent mixed basis object as linear combination")
		else
			if K==Ket
				if size(coeffs)==(length(basis),)
					new(coeffs, basis)
				elseif size(coeffs)==(length(basis),1)
					new(vec(coeffs), basis)
				else
					error("Dimensions of coefficient array does not match basis")
				end
			else
				if size(coeffs)==(1,length(basis))
					new(coeffs, basis)
				elseif length(coeffs)==1
			 		new(coeffs.', basis)
			 	else
			 		error("Dimensions of coefficient array does not match basis")
			 	end
			end
		end
	end
end

DiracVector{T,K}(coeffs::Array{T}, basis::AbstractBasis{K}) = DiracVector{T,K}(coeffs, basis)
copy(d::DiracVector) = DiracVector(copy(d.coeffs), copy(d.basis))
isequal(a::DiracVector, b::DiracVector) = isequal(a.coeffs, b.coeffs) && a.basis==b.basis
==(a::DiracVector, b::DiracVector) = a.coeffs==b.coeffs && a.basis==b.basis
basislabel(d::DiracVector) = label(d.basis)
isdual{A,B}(a::DiracVector{A,Ket}, b::DiracVector{B,Bra}) = isdual(a.basis, b.basis) && a.coeffs==b.coeffs'
isdual{A,B}(a::DiracVector{A,Bra}, b::DiracVector{B,Ket}) = isdual(b,a)
isdual{A,B,K}(a::DiracVector{A,K}, b::DiracVector{B,K}) = false #default to false

#####################################
#Show Functions######################
#####################################

function showcompact(io::IO, d::DiracVector)
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else
		print(io, "($(d.coeffs[1])$(d.basis[1]))",[" + ($(d.coeffs[i])$(d.basis[i]))" for i=2:length(d)]...)
	end
end
function show(io::IO, d::DiracVector)
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else	
		println("$(typeof(d)):")
		table = cell(length(d), 2)	
		if length(d.coeffs)>=52
			for i=1:25
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
			table[26:(length(d.coeffs)-25),:] = 0 # prevents access to undefined reference
			for i=(length(d.coeffs)-25):length(d.coeffs)
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
		else
			for i=1:length(d.coeffs)
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
#Array/Dict Functions################
#####################################
ctranspose(d::DiracVector) = DiracVector(d.coeffs', d.basis')
getindex(d::DiracVector, x) = d.coeffs[x]
length(d::DiracVector) = length(d.coeffs)
size(d::DiracVector, args...) = size(d.coeffs, args...)
kind(d::DiracVector) = kind(d.basis)
setindex!(d::DiracVector, y, x) = setindex!(d.coeffs, y, x)
endof(d::DiracVector) = length(d)
find(d::DiracVector) = find(d.coeffs)
find(f::Function, d::DiracVector) = find(f, d.coeffs)
findstates(f::Function, d::DiracVector) = find(f, d.basis)
getpos(d::DiracVector, s::AbstractState) = get(d.basis, s)
function get(d::DiracVector, s::AbstractState, notfound)
	try
		return d[getpos(d, s)]
	catch
		return notfound
	end
end

get(d::DiracVector, s::AbstractState) = d[getpos(d, s)]
get(d::DiracVector, label) = get(d, typeof(d.basis)<:Basis ? State(label) : TensorState(label))

#####################################
#Function-passing Functions##########
#####################################

map(f::Function, d::DiracVector)=DiracVector(map(f, d.coeffs), d.basis)

function map!(f::Function, d::DiracVector)
	d.coeffs = map(f, d.coeffs)
	return d
end

function mapmatch(fstates::Function, fcoeffs::Function, d::DiracVector)
	matched = findstates(fstates, d)	
	coeffs = convert(Array{Any}, d.coeffs)
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

function filterstates(f::Function, d::DiracVector)
	newbasis = filter(f, d.basis)
	if kind(d)==Ket
		return DiracVector(vcat([get(d, newbasis[i]) for i=1:length(newbasis)]...), newbasis)
	else 
		return DiracVector(hcat([get(d, newbasis[i]) for i=1:length(newbasis)]...), newbasis)
	end
end

function filtercoeffs(f::Function, d::DiracVector)
	matched = find(map(f, d.coeffs))
	newbasis = filter(x->in(get(d.basis,x), matched), d.basis)
	if kind(d)==Ket
		return DiracVector(vcat(getindex(d, matched)...), newbasis)
	else
		return DiracVector(hcat(getindex(d, matched)...), newbasis)
	end
end

qeval(f::Function, d::DiracVector) = map(x->qeval(f, x), d)

#####################################
#Arithmetic Functions################
#####################################

for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracVector, b::DiracVector) = DiracVector(($op)(a.coeffs,b.coeffs), a.basis)
	@eval ($op)(n, d::DiracVector) = DiracVector(($op)(n,d.coeffs), d.basis)
	@eval ($op)(d::DiracVector, n) = DiracVector(($op)(d.coeffs,n), d.basis)
end

*(c::DiracCoeff, d::DiracVector) = c.*d
*(d::DiracVector, c::DiracCoeff) = c*d

function *{T}(s::AbstractState{Bra}, d::DiracVector{T, Ket})
	if samebasis(d,s)	
		return get(d, s', 0)
	else
		return reduce(+, [d[i]*(s*d.basis[i]) for i=1:length(d)])
	end
end
function *{T}(d::DiracVector{T, Bra}, s::AbstractState{Ket})
	if samebasis(d,s)	
		return get(d, s', 0)
	else
		return reduce(+, [d[i]*(d.basis[i]*s) for i=1:length(d)])
	end
end

function *{N1<:Number, N2<:Number}(a::DiracVector{N1, Bra}, b::DiracVector{N2, Ket})
	if a.basis==b.basis 
		return (a.coeffs*b.coeffs)[1]
	else
		return reduce(+, [a[i]*b[j]*(a.basis[i]*b.basis[j]) for i=1:length(a), j=1:length(b)])
	end
end

function *{A, B}(a::DiracVector{A, Bra}, b::DiracVector{B, Ket})
	if a.basis==b.basis 
		return reduce(+, [a[i]*b[i]*(a.basis[i]*b.basis[i]) for i=1:length(a)]) 
	else
		return reduce(+, [a[i]*b[j]*(a.basis[i]*b.basis[j]) for i=1:length(a), j=1:length(b)])
	end
end	

*{A, B}(a::DiracVector{A, Ket}, b::DiracVector{B, Bra}) = DiracMatrix(a.coeffs*b.coeffs, copy(a.basis), copy(b.basis))
*{A, B, K}(a::DiracVector{A, K}, b::DiracVector{B, K}) = DiracVector(a.coeffs*b.coeffs, a.basis*b.basis)
*{T,K}(s::AbstractState{K}, d::DiracVector{T, K}) = DiracVector(d.coeffs, map(x->s*x, d.basis))
*{T, K}(d::DiracVector{T, K}, s::AbstractState{K}) = DiracVector(d.coeffs, map(x->x*s, d.basis))

function +{T,K}(d::DiracVector{T,K}, s::AbstractState{K})
	if in(s, d.basis)
		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
		res[getpos(d,s)] = 1+get(res, s)
		return res
	elseif samebasis(d,s)
		if K==Ket
			return DiracVector(vcat(d.coeffs, 1), d.basis+tobasis(s))
		else
			return DiracVector(hcat(d.coeffs, 1), d.basis+tobasis(s))
		end
	else
		error("BasisMismatch")
	end
end

function +{T,K}(s::AbstractState{K}, d::DiracVector{T,K})
	if in(s, d.basis)
		res = 1*d #forces the coeff array to accept numbers if it is InnerProduct; hacky but works
		res[getpos(d,s)] = 1+get(res, s)
		return res
	elseif samebasis(d,s)
		if K==Ket
			return DiracVector(vcat(1, d.coeffs), tobasis(s)+d.basis)
		else
			return DiracVector(hcat(1, d.coeffs), tobasis(s)+d.basis)
		end
	else
		error("BasisMismatch")
	end
end

function +{T1,T2,K}(a::DiracVector{T1,K}, b::DiracVector{T2,K})
	if a.basis==b.basis
		return DiracVector(a.coeffs+b.coeffs, a.basis)
	elseif samebasis(a,b)
		res = 1*a
		for i=1:length(b)
			res = res+b.basis[i]
			res[getpos(res, b.basis[i])] = get(res, b.basis[i]) + b[i] - 1
		end
		return res
	else
		error("BasisMismatch")

	end
end

-{T,K}(d::DiracVector{T,K}, s::AbstractState{K}) = d+(-s)
-{T,K}(s::AbstractState{K}, d::DiracVector{T,K}) = s+(-d)
-{T1,T2,K}(a::DiracVector{T1,K}, b::DiracVector{T2,K}) = a+(-b)

+{K}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1, 1], Basis([a,b])) 
-{K}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1,-1], Basis([a,b])) 
-(s::AbstractState) = -1*s
-(d::DiracVector) = DiracVector(-1*d.coeffs, d.basis)

norm(v::Vector, p::Int=2) = reduce(+, [(abs(i))^p for i in v])^(1/p)  
norm(d::DiracVector, p::Int=2) = norm(d.coeffs)

normalize(v::Vector) = (1/norm(v))*v
normalize(d::DiracVector) = DiracVector(normalize(d.coeffs), d.basis)
