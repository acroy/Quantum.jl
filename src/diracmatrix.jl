#####################################
#DiracMatrix#########################
#####################################

type DiracMatrix{T} <: Dirac
	coeffs::Matrix{T}
	rowbasis::AbstractBasis{Ket}
	colbasis::AbstractBasis{Bra}
	function DiracMatrix(coeffs, rowbasis, colbasis)
		if size(coeffs)==(length(rowbasis), length(colbasis))
			new(coeffs, rowbasis, colbasis)
		else
			throw(DimensionMismatch("Bases do not match coefficient matrix"))
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

basislabel(op::DiracMatrix) = [label(op.rowbasis), label(op.colbasis)]
isdual(a::DiracMatrix, b::DiracMatrix) = a'==b

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
#Matrix/Dict Functions###############
#####################################
isequal(a::DiracMatrix, b::DiracMatrix) = isequal(a.coeffs,b.coeffs) && a.rowbasis==b.rowbasis && a.colbasis==b.colbasis 
==(a::DiracMatrix, b::DiracMatrix) = a.coeffs==b.coeffs && a.rowbasis==b.rowbasis && a.colbasis==b.colbasis 

ndims(op::DiracMatrix) = ndims(op.coeffs)
size(op::DiracMatrix, args...) = size(op.coeffs, args...)
length(op::DiracMatrix) = length(op.coeffs)

endof(op::DiracMatrix) = length(op)
find(op::DiracMatrix) = find(op.coeffs)
find(f::Function, op::DiracMatrix) = find(f, op.coeffs)
findstates(f::Function, op::DiracMatrix) = find(f, [op.rowbasis[i]*op.colbasis[j] for i=1:length(op.rowbasis), j=1:length(op.colbasis)]) #f takes outer product as argument

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
#Arithmetic Functions################
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
function *{T}(op::DiracMatrix, d::DiracVector{T, Ket})
	if isdual(op.colbasis, d.basis)
		return DiracVector(op.coeffs*d.coeffs, op.rowbasis)
	else
		return sum([op.rowbasis[i]*sum([op[i,j]*(op.colbasis[j]*d) for j=1:length(op.colbasis)]) for i=1:length(op.rowbasis)])
	end
end

function *{T}(d::DiracVector{T, Bra}, op::DiracMatrix) 
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
		terms = [(a[i,j]*b[m,n]*(a.colbasis[j]*b.rowbasis[m]))*(a.rowbasis[i]*b.colbasis[n]) 
						 for n=1:size(b,2), m=1:size(b,1), j=1:size(a,2), i=1:size(a,1)]
		terms = filter(x->x[1]!=0, terms)
		return isempty(terms) ? 0 : reduce(+, terms)
	end
end

function +(op::DiracMatrix, o::OuterProduct)
	if in(o.bra, op.colbasis) && in(o.ket, op.rowbasis)
		res = 1*op
		res[getpos(op, o)...] = 1+get(op, o)
		return res
	elseif samebasis(op, o)
		#unecessary rehashing occurs here...
		rowb = tobasis(vcat(op.rowbasis[:], o.ket))
		colb = tobasis(vcat(op.colbasis[:], o.bra))
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
		rowb = tobasis(vcat(o.ket, op.rowbasis[:]))
		colb = tobasis(vcat(o.bra, op.colbasis[:]))
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

kron(a::DiracMatrix, b::DiracMatrix) = DiracMatrix(kron(a.coeffs, b.coeffs), tensor(a.rowbasis, b.rowbasis), tensor(a.colbasis, b.colbasis)) 
kron(op::DiracMatrix, d::DiracVector{Ket}) = DiracMatrix(kron(op.coeffs, d.coeffs), tensor(op.rowbasis, d.basis), op.colbasis)
kron(op::DiracMatrix, d::DiracVector{Bra}) = DiracMatrix(kron(op.coeffs, d.coeffs), op.rowbasis, tensor(op.colbasis, d.basis))
kron(d::DiracVector{Ket}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), tensor(d.basis, op.rowbasis), op.colbasis)
kron(d::DiracVector{Bra}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), op.rowbasis, tensor(d.basis, op.colbasis))
kron{K}(a::DiracVector{K}, b::DiracVector{K}) = DiracVector(kron(a.coeffs, b.coeffs), tensor(a.basis, b.basis))
kron(a::DiracVector{Ket}, b::DiracVector{Bra}) = DiracMatrix(kron(a.coeffs, b.coeffs), a.basis, b.basis)
kron(a::DiracVector{Bra}, b::DiracVector{Ket}) = kron(b,a)

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
