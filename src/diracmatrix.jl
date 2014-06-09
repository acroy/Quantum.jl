#####################################
#DiracMatrix#########################
#####################################

type DiracMatrix{T} <: Dirac
	coeffs::Matrix{T}
	rowbasis::AbstractBasis{Ket}
	colbasis::AbstractBasis{Bra}
	function DiracMatrix(coeffs, rowbasis, colbasis)
		if samebasis("?", rowbasis) || samebasis("?", colbasis)
			error("BasisMismatch: cannot represent mixed basis object as linear combination")
		else
			if size(coeffs)==(length(rowbasis), length(colbasis))
				new(coeffs, rowbasis, colbasis)
			else
				throw(DimensionMismatch("Bases do not match coefficient matrix"))
			end
		end
	end
end

DiracMatrix{T}(coeffs::Matrix{T}, rowbasis::AbstractBasis{Ket}, colbasis::AbstractBasis{Bra}) = DiracMatrix{T}(coeffs, rowbasis, colbasis) 
DiracMatrix(coeffs::Matrix, b::AbstractBasis{Ket}) = DiracMatrix(coeffs, b, b') 
DiracMatrix(coeffs::Matrix, b::AbstractBasis{Bra}) = DiracMatrix(coeffs, b', b) 

function DiracMatrix(fcoeff::Function, fstate::Function, b::AbstractBasis{Ket})
	coeffs = convert(Array{Any}, zeros(length(b), length(b)))
	for i=1:length(b)
		for j=1:length(b)
			coeffs[i,j] = fcoeff(b[j])*(b[i]'*fstate(b[j]))
		end
	end
	return DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), b) #use vcat()/hcat() trick to convert to most primitive common type
end

#####################################
#Misc Functions######################
#####################################

basislabel(op::DiracMatrix) = [label(op.rowbasis), label(op.colbasis)]

isequal(a::DiracMatrix, b::DiracMatrix) = isequal(a.coeffs,b.coeffs) && a.rowbasis==b.rowbasis && a.colbasis==b.colbasis 
==(a::DiracMatrix, b::DiracMatrix) = a.coeffs==b.coeffs && a.rowbasis==b.rowbasis && a.colbasis==b.colbasis 

isdual(a::DiracMatrix, b::DiracMatrix) = a'==b

ndims(op::DiracMatrix) = ndims(op.coeffs)
size(op::DiracMatrix, args...) = size(op.coeffs, args...)
length(op::DiracMatrix) = length(op.coeffs)

endof(op::DiracMatrix) = length(op)
find(op::DiracMatrix) = find(op.coeffs)
eltype(op::DiracMatrix) = eltype(op.coeffs)

ctranspose(op::DiracMatrix) = DiracMatrix(op.coeffs', op.colbasis', op.rowbasis')
getindex(op::DiracMatrix, x...) = op.coeffs[x...]
setindex!(op::DiracMatrix, y, x...) = setindex!(op.coeffs,y,x...)

getpos(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}) = (get(op.rowbasis, k), get(op.colbasis, b))
getpos(op::DiracMatrix, o::OuterProduct) = getpos(op, o.ket, o.bra)

get(op::DiracMatrix, s::AbstractState{Ket}) = DiracVector(op[get(op.rowbasis, s), :], op.colbasis)
get(op::DiracMatrix, s::AbstractState{Bra}) = DiracVector(op[:, get(op.colbasis, s)], op.rowbasis)
get(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}) = op[get(op.rowbasis, k), get(op.colbasis, b)]
get(op::DiracMatrix, o::OuterProduct) = get(op, o.ket, o.bra)

function get(op::DiracMatrix, s::AbstractState, notfound)
	try
		return get(op, s)
	catch
		return notfound
	end
end

function get(op::DiracMatrix, k::AbstractState{Ket}, b::AbstractState{Bra}, notfound)
	try
		return get(op, k, b)
	catch
		return notfound
	end
end
#####################################
#Show Functions######################
#####################################
function showcompact(io::IO, op::DiracMatrix)
	if length(op.coeffs)==0
		print(io, "$(typeof(op))[]")
	else
		tempio = IOBuffer()
		print(tempio, [" + ($(op.coeffs[i,j])$(op.rowbasis[i])$(op.colbasis[j]))" for i=1:length(op.rowbasis), j=1:length(op.colbasis)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end
function show(io::IO, op::DiracMatrix)
	println("$(typeof(op)):")
	table = cell(length(op.rowbasis)+1, length(op.colbasis)+1)	
	for i = 1:length(op.rowbasis)
		table[i+1,1] = op.rowbasis[i]
	end
	for j = 1:length(op.colbasis)
		table[1,j+1] = op.colbasis[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = op.coeffs	
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	print(io, io_str[searchindex(io_str, "\n")+3:end])
end

#####################################
#Function-passing Functions##########
#####################################

find(op::DiracMatrix) = find(op.coeffs)
find(f::Function, op::DiracMatrix) = find(f, op.coeffs)
findstates(f::Function, op::DiracMatrix) = find(f, [op.rowbasis[i]*op.colbasis[j] for i=1:length(op.rowbasis), j=1:length(op.colbasis)]) #f takes OuterProduct as argument

map(f::Function, op::DiracMatrix) = DiracMatrix(map(f, op.coeffs), op.rowbasis, op.colbasis)

function map!(f::Function, op::DiracMatrix)
	op.coeffs = map(f, op.coeffs)
	return op
end

function mapmatch(fstates::Function, fcoeffs::Function, op::DiracMatrix) #fstates takes OuterProduct as argument
	matched = findstates(fstates, op)
	coeffs = convert(Array{Any}, op.coeffs)
	for i in matched
		coeffs[i] = fcoeffs(coeffs[i])
	end
	#hacky concat strategy to force corect typing
	DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), op.rowbasis, op.colbasis)
end

function mapmatch!(fstates::Function, fcoeffs::Function, op::DiracMatrix)
	matched = findstates(fstates, op)
	for i in matched
		op[i] = fcoeffs(op[i])
	end
	return op
end

qeval(f::Function, op::DiracMatrix) = map(x->qeval(f, x), op)

#####################################
#Arithmetic Operations###############
#####################################
for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracMatrix, b::DiracVector) = DiracMatrix(($op)(a.coeffs,b.coeffs), a.rowbasis, a.colbasis)
	@eval ($op)(a::DiracVector, b::DiracMatrix) = DiracMatrix(($op)(a.coeffs,b.coeffs), b.rowbasis, b.colbasis)
	@eval ($op)(a::DiracMatrix, b::DiracMatrix) = DiracMatrix(($op)(a.coeffs,b.coeffs), a.rowbasis, a.colbasis)
	@eval ($op)(n, d::DiracMatrix) = DiracMatrix(($op)(n,d.coeffs), d.rowbasis, d.colbasis)
	@eval ($op)(d::DiracMatrix, n) = DiracMatrix(($op)(d.coeffs,n), d.rowbasis, d.colbasis)
end

/(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs/d, op.rowbasis, op.colbasis)
*(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs*d, op.rowbasis, op.colbasis)
*(d::DiracCoeff, op::DiracMatrix) = DiracMatrix(d*op.coeffs, op.rowbasis, op.colbasis)

function *(op::DiracMatrix, s::AbstractState{Ket}) 
	if in(s', op.colbasis) 
		return get(op, s') 
	else
		return sum([op.rowbasis[i]*sum([op[i,j]*(op.colbasis[j]*d) for j=1:length(op.colbasis)]) for i=1:length(op.rowbasis)])
	end
end
function *(s::AbstractState{Bra}, op::DiracMatrix)
	if in(s', op.rowbasis) 
		return get(op, s')
	else
		return sum([op.colbasis[j]*sum([op[i,j]*(d*op.rowbasis[i]) for j=1:length(op.rowbasis)]) for i=1:length(op.colbasis)])
	end
end
function *(op::DiracMatrix, d::DiracVector{Ket})
	if isdual(op.colbasis, d.basis)
		return DiracVector(op.coeffs*d.coeffs, op.rowbasis)
	else
		return sum([op.rowbasis[i]*sum([op[i,j]*(op.colbasis[j]*d) for j=1:length(op.colbasis)]) for i=1:length(op.rowbasis)])
	end
end

function *(d::DiracVector{Bra}, op::DiracMatrix) 
	if isdual(op.rowbasis, d.basis)
		return DiracVector(d.coeffs*op.coeffs, op.colbasis)
	else
		return sum([op.colbasis[j]*sum([op[i,j]*(d*op.rowbasis[i]) for j=1:length(op.rowbasis)]) for i=1:length(op.colbasis)])
	end
end

function *(a::DiracMatrix, b::DiracMatrix)
	if isdual(a.colbasis, b.rowbasis)
		return DiracMatrix(a.coeffs*b.coeffs, a.rowbasis, b.colbasis)
	else
		return sum([(a[i,j]*b[m,n]*(a.colbasis[j]*b.rowbasis[m]))*(a.rowbasis[i]*b.colbasis[n]) 
						 for n=1:size(b,2), m=1:size(b,1), j=1:size(a,2), i=1:size(a,1)])
	end
end

function +(op::DiracMatrix, o::OuterProduct)
	if in(o.bra, op.colbasis) && in(o.ket, op.rowbasis)
		res = 1*op
		res[getpos(op, o)...] = 1+get(op, o)
		return res
	elseif samebasis(op, o)
		#unecessary rehashing occurs here...
		rowb = basisjoin(op.rowbasis, o.ket)
		colb = basisjoin(op.colbasis, o.bra)
		res = DiracMatrix(convert(typeof(op.coeffs), zeros(length(rowb),length(colb))), rowb, colb)
		res[1:size(op,1), 1:size(op,2)] = op.coeffs
		res[getpos(res, o)...] = 1+get(res, o)
		return res
	else
		error("BasisMismatch")
	end
end

function +(o::OuterProduct, op::DiracMatrix)
	if in(o.bra, op.colbasis) && in(o.ket, op.rowbasis)
		res = 1*op
		res[getpos(op, o)...] = 1+get(op, o)
		return res
	elseif samebasis(o,op)
		#unecessary rehashing occurs here...
		rowb = basisjoin(o.ket, op.rowbasis)
		colb = basisjoin(o.bra, op.colbasis)
		res = DiracMatrix(convert(typeof(op.coeffs), zeros(length(rowb),length(colb))), rowb, colb)
		res[(size(res,1)-size(op,1)+1):size(op,1), (size(res,2)-size(op,2)+1):size(op,2)] = op.coeffs
		res[getpos(res, o)...] = 1+get(res, o)
		return res
	else
		error("BasisMismatch")
	end
end


function +(a::DiracMatrix, b::DiracMatrix)
	if a.rowbasis==b.rowbasis && a.colbasis==b.colbasis 
		return DiracMatrix(a.coeffs+b.coeffs, a.rowbasis, a.colbasis) 
	elseif samebasis(a, b)
		for i=1:size(b, 1)
			for j=1:size(b, 2)
				if b[i,j]!=0
					a=a+(b.rowbasis[i]*b.colbasis[j])
					a[getpos(a, b.rowbasis[i], b.colbasis[j])...] = get(a, b.rowbasis[i], b.colbasis[j])+b[i,j]-1
				end
			end
		end
		return a
	else
		error("BasisMismatch")
	end
end

-(op::DiracMatrix) = -1*op
-(a::DiracMatrix, b::DiracMatrix) = a+(-b)

exp(op::DiracMatrix) = DiracMatrix(exp(op.coeffs), op.rowbasis, op.colbasis)
expm(op::DiracMatrix) = DiracMatrix(expm(op.coeffs), op.rowbasis, op.colbasis)
^(op::DiracMatrix, i::Integer) = DiracMatrix(^(op.coeffs, i), op.rowbasis, op.colbasis)
^(op::DiracMatrix, n::Number) = DiracMatrix(^(op.coeffs, n), op.rowbasis, op.colbasis)
trace(op::DiracMatrix) = trace(op.coeffs)
commutator(a::DiracMatrix, b::DiracMatrix) = (a*b) - (b*a)

function ptrace(op::DiracMatrix, ind::Int)
	if isdual(op.colbasis, op.rowbasis)
		trrow = tensor(vcat(separate(op.rowbasis)[1:ind-1], separate(op.rowbasis)[ind+1:end])...)
		trcol = trrow'
		len = length(trcol)
		coeffs = [sum([op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.rowbasis)[ind])]) for j=1:length(trrow), k=1:length(trcol)] 
	else
		error("BasisMismatch")
	end
	return DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), trrow, trcol) #use vcat()/hcat() trick to convert to most primitive common type
end
