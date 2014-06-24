#####################################
#DiracMatrix#########################
#####################################

type DiracMatrix{K<:Ket, B<:Bra, T} <: Dirac
	coeffs::SparseMatrixCSC{T, Int}
	rowb::AbstractBasis{K}
	colb::AbstractBasis{B}
	function DiracMatrix(coeffs, rowb, colb)
		@assert size(coeffs)==(length(rowb), length(colb)) "Dimensions of coefficient array do not match bases"
		return new(coeffs, rowb, colb)
	end
end
DiracMatrix{K<:Ket, B<:Bra, T}(coeffs::SparseMatrixCSC{T}, rowb::AbstractBasis{K}, colb::AbstractBasis{B}) = DiracMatrix{K,B,T}(coeffs, rowb, colb) 
DiracMatrix{K<:Ket, B<:Bra, T}(coeffs::Matrix{T}, rowb::AbstractBasis{K}, colb::AbstractBasis{B}) = DiracMatrix{K,B,T}(sparse(coeffs), rowb, colb) 
DiracMatrix{K<:Ket}(coeffs::AbstractArray, b::AbstractBasis{K}) = DiracMatrix(coeffs, b, b') 
DiracMatrix{B<:Bra}(coeffs::AbstractArray, b::AbstractBasis{B}) = DiracMatrix(coeffs, b', b) 

function constructop!{K<:Ket}(coeffs::SparseMatrixCSC, fcoeff::Function, fstate::Function, b::AbstractBasis{K})
	for i=1:length(b)
		for j=1:length(b)
			coeffs[i,j] = fcoeff(b[j])*inner(b[i]',fstate(b[j]))
		end
	end
end

function DiracMatrix{K<:Ket}(fcoeff::Function, fstate::Function, b::AbstractBasis{K}, t::DataType=typeof(fcoeff(b[1])))
	coeffs = convert(SparseMatrixCSC{t}, spzeros(length(b), length(b)))
	constructop!(coeffs, fcoeff, fstate, b)
	return DiracMatrix(coeffs, b)
end

dmat=DiracMatrix

#####################################
#Getter-style Functions##############
#####################################
bsym(op::DiracMatrix) = [bsym(op.rowb), bsym(op.colb)]

#####################################
#Boolean Functions###################
#####################################
isequal(a::DiracMatrix, b::DiracMatrix) = isequal(a.coeffs,b.coeffs) && a.rowb==b.rowb && a.colb==b.colb 
==(a::DiracMatrix, b::DiracMatrix) = a.coeffs==b.coeffs && a.rowb==b.rowb && a.colb==b.colb 
isdual(a::DiracMatrix, b::DiracMatrix) = a.coeffs'==b.coeffs && a.rowb==b.colb && b.colb==a.rowb

#####################################
#Array-like Functions################
#####################################
size(op::DiracMatrix, args...) = size(op.coeffs, args...)
for op=(:endof, :eltype, :length, :find, :findn, :findnz, :nnz,:ndims)
	@eval ($op)(d::DiracMatrix) = ($op)(d.coeffs)
end
ctranspose(op::DiracMatrix) = DiracMatrix(op.coeffs', op.colb', op.rowb')
getindex(op::DiracMatrix, x...) = op.coeffs[x...]
setindex!(op::DiracMatrix, y, x...) = setindex!(op.coeffs,y,x...)

# #####################################
# #Dict-like Functions#################
# #####################################
# getpos(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}) = (get(op.rowb, k), get(op.colb, b))
# getpos(op::DiracMatrix, o::OuterProduct) = getpos(op, o.ket, o.bra)

# get(op::DiracMatrix, s::AbstractState{Ket}) = DiracVector(op[get(op.rowb, s), :], op.colb)
# get(op::DiracMatrix, s::AbstractState{Bra}) = DiracVector(op[:, get(op.colb, s)], op.rowb)
# get(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}) = op[get(op.rowb, k), get(op.colb, b)]
# get(op::DiracMatrix, o::OuterProduct) = get(op, o.ket, o.bra)

# function get(op::DiracMatrix, s::AbstractState, notfound)
# 	try
# 		return get(op, s)
# 	catch
# 		return notfound
# 	end
# end

# function get(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}, notfound)
# 	try
# 		return get(op, k, b)
# 	catch
# 		return notfound
# 	end
# end

# #####################################
# #Show Functions######################
# #####################################
summary(op::DiracMatrix) = "$(size(op,1))x$(size(op,2)) $(typeof(op))"

function showcompact(io::IO, op::DiracMatrix)
	if length(op.coeffs)==0
		print(io, "$(typeof(op))[]")
	else
		tempio = IOBuffer()
		nz = hcat(findn(op)...)
		print(tempio, [" + $(op.coeffs[nz[i,1],nz[i,2]])$(op.rowb[nz[i,1]])$(op.colb[nz[i,2]])" for i=1:size(nz,1)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end

function show(io::IO, op::DiracMatrix)
	println(io, summary(op))
	table = cell(length(op.rowb)+1, length(op.colb)+1)	
	for i = 1:length(op.rowb)
		table[i+1,1] = op.rowb[i]
	end
	for j = 1:length(op.colb)
		table[1,j+1] = op.colb[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = full(op.coeffs)
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	print(io, io_str[searchindex(io_str, "\n")+3:end])
end

# #####################################
# #Function-passing Functions##########
# #####################################
# find(f::Function, op::DiracMatrix) = find(f, op.coeffs)
# findstates(f::Function, op::DiracMatrix) = find(f, [op.rowb[i]*op.colb[j] for i=1:length(op.rowb), j=1:length(op.colb)]) #f takes OuterProduct as argument

# map(f::Function, op::DiracMatrix) = DiracMatrix(map(f, op.coeffs), op.rowb, op.colb)

# function map!(f::Function, op::DiracMatrix)
# 	op.coeffs = map(f, op.coeffs)
# 	return op
# end

# function mapmatch(fstates::Function, fcoeffs::Function, op::DiracMatrix) #fstates takes OuterProduct as argument
# 	matched = findstates(fstates, op)
# 	coeffs = convert(Array{Any}, op.coeffs)
# 	for i in matched
# 		coeffs[i] = fcoeffs(coeffs[i])
# 	end
# 	#hacky concat strategy to force corect typing
# 	DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), op.rowb, op.colb)
# end

# function mapmatch!(fstates::Function, fcoeffs::Function, op::DiracMatrix)
# 	matched = findstates(fstates, op)
# 	for i in matched
# 		op[i] = fcoeffs(op[i])
# 	end
# 	return op
# end

# qeval(f::Function, op::DiracMatrix) = map(x->qeval(f, x), op)

# #####################################
# #Arithmetic Operations###############
# #####################################
# for op=(:.*,:.-,:.+,:./,:.^)
# 	@eval ($op)(a::DiracMatrix, b::DiracVector) = DiracMatrix(($op)(a.coeffs,b.coeffs), a.rowb, a.colb)
# 	@eval ($op)(a::DiracVector, b::DiracMatrix) = DiracMatrix(($op)(a.coeffs,b.coeffs), b.rowb, b.colb)
# 	@eval ($op)(a::DiracMatrix, b::DiracMatrix) = DiracMatrix(($op)(a.coeffs,b.coeffs), a.rowb, a.colb)
# 	@eval ($op)(n, d::DiracMatrix) = DiracMatrix(($op)(n,d.coeffs), d.rowb, d.colb)
# 	@eval ($op)(d::DiracMatrix, n) = DiracMatrix(($op)(d.coeffs,n), d.rowb, d.colb)
# end

# /(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs/d, op.rowb, op.colb)
# *(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs*d, op.rowb, op.colb)
# *(d::DiracCoeff, op::DiracMatrix) = DiracMatrix(d*op.coeffs, op.rowb, op.colb)

# function *(op::DiracMatrix, s::AbstractState{Ket}) 
# 	if in(s', op.colb) 
# 		return get(op, s') 
# 	else
# 		return reduce(+,[op.rowb[i]*reduce(+,[op[i,j]*(op.colb[j]*d) for j=1:length(op.colb)]) for i=1:length(op.rowb)])
# 	end
# end
# function *(s::AbstractState{Bra}, op::DiracMatrix)
# 	if in(s', op.rowb) 
# 		return get(op, s')
# 	else
# 		return reduce(+,[op.colb[j]*reduce(+,[op[i,j]*(d*op.rowb[i]) for j=1:length(op.rowb)]) for i=1:length(op.colb)])
# 	end
# end
# function *(op::DiracMatrix, d::DiracVector{Ket})
# 	if isdual(op.colb, d.basis)
# 		return DiracVector(op.coeffs*d.coeffs, op.rowb)
# 	else
# 		return reduce(+,[op.rowb[i]*reduce(+,[op[i,j]*(op.colb[j]*d) for j=1:length(op.colb)]) for i=1:length(op.rowb)])
# 	end
# end

# function *(d::DiracVector{Bra}, op::DiracMatrix) 
# 	if isdual(op.rowb, d.basis)
# 		return DiracVector(d.coeffs*op.coeffs, op.colb)
# 	else
# 		return reduce(+,[op.colb[j]*reduce(+,[op[i,j]*(d*op.rowb[i]) for j=1:length(op.rowb)]) for i=1:length(op.colb)])
# 	end
# end

# function *(a::DiracMatrix, b::DiracMatrix)
# 	if isdual(a.colb, b.rowb)
# 		return DiracMatrix(a.coeffs*b.coeffs, a.rowb, b.colb)
# 	else
# 		return reduce(+,[(a[i,j]*b[m,n]*(a.colb[j]*b.rowb[m]))*(a.rowb[i]*b.colb[n]) 
# 						 for n=1:size(b,2), m=1:size(b,1), j=1:size(a,2), i=1:size(a,1)])
# 	end
# end

# function +(op::DiracMatrix, o::OuterProduct)
# 	if in(o.ket, op.rowb) && in(o.bra, op.colb)
# 		res = 1*op
# 		res[getpos(op, o)...] = 1+get(op, o)
# 		return res
# 	elseif in(o.ket, op.rowb) && samebasis(o.bra, op.colb)
# 		newcol = zeros(size(op,1))
# 		newcol[get(op.rowb, o.ket)] = 1 
# 		return DiracMatrix(hcat(op.coeffs,newcol), op.rowb, basisjoin(op.colb, o.bra))
# 	elseif samebasis(o.ket, op.rowb) && in(o.bra, op.colb)
# 		newrow = zeros(1,size(op,2))
# 		newrow[get(op.colb, o.bra)] = 1 
# 		return DiracMatrix(vcat(op.coeffs,newrow), basisjoin(op.rowb, o.ket), op.colb)
# 	else
# 		@assert samebasis(op, o) "BasisMismatch"
# 		rowb = basisjoin(op.rowb, o.ket)
# 		colb = basisjoin(op.colb, o.bra)
# 		res = DiracMatrix(convert(typeof(op.coeffs), zeros(length(rowb),length(colb))), rowb, colb)
# 		res[1:size(op,1), 1:size(op,2)] = op.coeffs
# 		res[getpos(res, o)...] = 1+get(res, o)
# 		return res
# 	end
# end

# function +(o::OuterProduct, op::DiracMatrix)
# 	if in(o.ket, op.rowb) && in(o.bra, op.colb)
# 		res = 1*op
# 		res[getpos(op, o)...] = 1+get(op, o)
# 		return res
# 	elseif in(o.ket, op.rowb) && samebasis(o.bra, op.colb)
# 		newcol = zeros(size(op,1))
# 		newcol[get(op.rowb, o.ket)] = 1 
# 		return DiracMatrix(hcat(newcol, op.coeffs), op.rowb, basisjoin(o.bra,op.colb))
# 	elseif samebasis(o.ket, op.rowb) && in(o.bra, op.colb)
# 		newrow = zeros(1,size(op,2))
# 		newrow[get(op.colb, o.bra)] = 1 
# 		return DiracMatrix(vcat(newrow,op.coeffs), basisjoin(o.ket,op.rowb), op.colb)
# 	else
# 		@assert samebasis(o,op) "BasisMismatch"
# 		rowb = basisjoin(o.ket, op.rowb)
# 		colb = basisjoin(o.bra, op.colb)
# 		res = DiracMatrix(convert(typeof(op.coeffs), zeros(length(rowb),length(colb))), rowb, colb)
# 		res[(size(res,1)-size(op,1)+1):size(op,1), (size(res,2)-size(op,2)+1):size(op,2)] = op.coeffs
# 		res[getpos(res, o)...] = 1+get(res, o)
# 		return res
# 	end
# end


# function +(a::DiracMatrix, b::DiracMatrix)
# 	if a.rowb==b.rowb && a.colb==b.colb 
# 		return DiracMatrix(a.coeffs+b.coeffs, a.rowb, a.colb) 
# 	else
# 		@assert samebasis(a, b) "BasisMismatch"
# 		for i=1:size(b, 1)
# 			for j=1:size(b, 2)
# 				if b[i,j]!=0
# 					a=a+(b.rowb[i]*b.colb[j])
# 					a[getpos(a, b.rowb[i], b.colb[j])...] = get(a, b.rowb[i], b.colb[j])+b[i,j]-1
# 				end
# 			end
# 		end
# 		return a
# 	end
# end

# -(op::DiracMatrix) = -1*op
# -(a::DiracMatrix, b::DiracMatrix) = a+(-b)
# -(a::DiracMatrix, b::OuterProduct) = a+(-b)

# exp(op::DiracMatrix) = DiracMatrix(exp(op.coeffs), op.rowb, op.colb)
# ^(op::DiracMatrix, i::Integer) = DiracMatrix(^(op.coeffs, i), op.rowb, op.colb)
# ^(op::DiracMatrix, n::Number) = DiracMatrix(^(op.coeffs, n), op.rowb, op.colb)
# trace(op::DiracMatrix) = trace(op.coeffs)
# commutator(a::DiracMatrix, b::DiracMatrix) = (a*b) - (b*a)

# function ptrace(op::DiracMatrix, ind::Int)
# 	@assert isdual(op.colb, op.rowb) "BasisMismatch"
# 	trrow = btensor(vcat(separate(op.rowb)[1:ind-1], separate(op.rowb)[ind+1:end])...)
# 	trcol = trrow'
# 	len = length(trcol)
# 	coeffs = [reduce(+,[op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.rowb)[ind])]) for j=1:length(trrow), k=1:length(trcol)] 
# 	return DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), trrow, trcol) #use vcat()/hcat() trick to convert to most primitive common type
# end
