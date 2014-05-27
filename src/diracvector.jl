#####################################
#DiracVector#########################
#####################################

type DiracVector{T,K<:BraKet} <: Dirac
	coeffs::Array{T} 
	basis::AbstractBasis{K}
	function DiracVector(coeffs, basis)
		if K==Ket
			if size(coeffs)==(length(basis),)
				new(coeffs, basis)
			elseif size(coeffs)==(length(basis),1)
				new(vec(coeffs), basis)
			else
				error("Dimensions of coefficient array does not match type Ket")
			end
		else
			if size(coeffs)==(1,length(basis))
				new(coeffs, basis)
			elseif length(coeffs)==1
		 		new(coeffs.', basis)
		 	else
		 		error("Dimensions of coefficient array does not match type Bra")
		 	end
		end
	end
end

DiracVector{T,K<:BraKet}(coeffs::Array{T}, basis::AbstractBasis{K}) = DiracVector{T,K}(coeffs, basis)

copy(d::DiracVector) = DiracVector(d.coeffs, d.basis)
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
size(d::DiracVector, x=nothing) = size(d.coeffs, x)
kind(d::DiracVector) = kind(d.basis)
setindex!(d::DiracVector, y, x) = setindex!(d.coeffs, y, x)
endof(d::DiracVector) = length(d)

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

#The vcat()/hcat() used below forces correct typing of the resultant coeff array, but it's sloppy. 
#I tried to define this as DiracVector(map(f, d.coeffs), d.basis), but it doesn't correctly reduce 
#the coeff array to "lowest" (i.e. most primitive) 
#common element type, and it would also yield InexactErrors for certain functions (e.g. qeval).

function map(f::Function, d::DiracVector)
	if kind(d)==Ket
		return DiracVector(vcat(map(f, d.coeffs)...), d.basis)
	else
		return DiracVector(hcat(map(f, d.coeffs)...), d.basis)
	end
end
function map!(f::Function, d::DiracVector)
	d.coeffs = kind(d)==Ket ? vcat(map(f, d.coeffs)...) : hcat(map(f, d.coeffs)...)
	return d
end

function mapmatch(fstates::Function, fcoeffs::Function, d::DiracVector)
	matched = map(x->getpos(d,x), filter(fstates, d.basis[:]))	
	#It would be more efficient to only loop over the matching
	#states, but there could be typing issues depending on 
	#what fcoeffs returns. Thus, we use the same hacky vcat method 
	#used in map(f::Function, d::DiracVector)
	if kind(d)==Ket
		return DiracVector(vcat([in(i, matched) ? fcoeffs(d[i]) : d[i] for i=1:length(d.coeffs)]...), d.basis)
	else
		return DiracVector(hcat([in(i, matched) ? fcoeffs(d[i]) : d[i] for i=1:length(d.coeffs)]...), d.basis)		
	end
end
function mapmatch!(fstates::Function, fcoeffs::Function, d::DiracVector)
	d.coeffs = mapmatch(fstates, fcoeffs, d).coeffs
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

*{T}(s::AbstractState{Bra}, d::DiracVector{T, Ket}) = get(d, s', 0)
*{T}(d::DiracVector{T, Bra}, s::AbstractState{Ket}) = get(d, s', 0)

*{N1<:Number, N2<:Number}(a::DiracVector{N1, Bra}, b::DiracVector{N2, Ket}) = (a.coeffs*b.coeffs)[1]
*{A, B}(a::DiracVector{A, Bra}, b::DiracVector{B, Ket}) = length(a)==length(b) ? reduce(+, [a[i]*b[i] for i=1:length(a)]) : throw(DimensionMismatch(""))
*{A, B}(a::DiracVector{A, Ket}, b::DiracVector{B, Bra}) = DiracMatrix(kron(a.coeffs,b.coeffs), a.basis, b.basis)
*{A, B, K<:BraKet}(a::DiracVector{A, K}, b::DiracVector{B, K}) = DiracVector(kron(a.coeffs, b.coeffs), a.basis*b.basis)


*{C<:DiracCoeff}(c::C, s::AbstractState) = DiracVector([c], statetobasis(s))
*{C<:DiracCoeff}(s::AbstractState, c::C) = *(c,s)

function +{T,K<:BraKet}(d::DiracVector{T,K}, s::AbstractState{K})
	if in(s, d.basis)
		d = 1*d #forces the coeff array to eltype DiracCoeff if it is InnerProduct; hacky but works
		d[getpos(d,s)] = 1+get(d, s)
		return d
	else
		if K==Ket
			return DiracVector(vcat(d.coeffs, 1), d.basis+statetobasis(s))
		else
			return DiracVector(hcat(d.coeffs, 1), d.basis+statetobasis(s))
		end
	end
end

function +{T,K<:BraKet}(s::AbstractState{K}, d::DiracVector{T,K})
	if in(s, d.basis)
		d = 1*d #forces the coeff array to eltype DiracCoeff if it is InnerProduct; hacky but works
		d[getpos(d,s)] = 1+get(d, s)
	else
		if K==Ket
			return DiracVector(vcat(1, d.coeffs), statetobasis(s)+d.basis)
		else
			return DiracVector(hcat(1, d.coeffs), statetobasis(s)+d.basis)
		end
	end
end

function +{T1,T2,K<:BraKet}(a::DiracVector{T1,K}, b::DiracVector{T2,K})
	if a.basis==b.basis
		return DiracVector(a.coeffs+b.coeffs, a.basis)
	else
		res = 1*copy(a)
		bdiff = setdiff(b.basis, a.basis) 
		compl = setdiff(b.basis, bdiff)
		for i in compl
			res[getpos(res, i)] = get(res, i) + get(b, i)
		end
		if K==Ket
			return DiracVector(vcat(res.coeffs, [get(b, bdiff[i]) for i=1:length(bdiff)]), res.basis+Basis(bdiff))
		else
			return DiracVector(hcat(res.coeffs, [get(b, bdiff[i]) for i=1:length(bdiff)]), res.basis+Basis(bdiff))
		end
	end
end

-{T<:DiracCoeff,K<:BraKet}(d::DiracVector{T,K}, s::AbstractState{K}) = d+(-s)
-{T,K<:BraKet}(s::AbstractState{K}, d::DiracVector{T,K}) = s+(-d)
-{T1,T2,K<:BraKet}(a::DiracVector{T1,K}, b::DiracVector{T2,K}) = a+(-b)

+{K<:BraKet}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1, 1], Basis([a,b])) 
-{K<:BraKet}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1,-1], Basis([a,b])) 
-(s::AbstractState) = -1*s
-(d::DiracVector) = DiracVector(-1*d.coeffs, d.basis)

norm(v::Vector, p::Int=2) = reduce(+, [(abs(i))^p for i in v])^(1/p)  
norm(d::DiracVector, p::Int=2) = norm(d.coeffs)

normalize(v::Vector) = (1/norm(v))*v
normalize(d::DiracVector) = DiracVector(normalize(d.coeffs), d.basis)
function normalize!(d::DiracVector) 
	d.coeffs=normalize(d.coeffs)
	return d
end
