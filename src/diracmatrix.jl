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

function DiracMatrix(fcoeff::Function, fstate::Function, b::AbstractBasis)
	coeffs = convert(Array{Any}, zeros(length(b), length(b)))
	for i=1:length(b)
		for j=1:length(b)
			ii = get(b, fstate(b[i]), -1)
			if ii != -1
				coeffs[ii,j] = fcoeff(b[i])*(b[j]'*b[i])
			end
		end
	end
	return DiracMatrix(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), b) #use vcat()/hcat() trick to convert to most primitive common type
end

#####################################
#Misc.Functions######################
#####################################

samebasis(a::DiracMatrix, b::DiracMatrix) = isequal(a.rowbasis,b.rowbasis) && isequal(a.colbasis, b.colbasis)

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
isequal(a::DiracMatrix, b::DiracMatrix) = isequal(a.coeffs,b.coeffs) && samebasis(a,b)
==(a::DiracMatrix, b::DiracMatrix) = a.coeffs==b.coeffs && samebasis(a,b)

ndims(op::DiracMatrix) = ndims(op.coeffs)
size(op::DiracMatrix, args...) = size(op.coeffs, args...)
length(op::DiracMatrix) = length(op.coeffs)

endof(op::DiracMatrix) = length(op)
find(op::DiracMatrix) = find(op.coeffs)
find(f::Function, op::DiracMatrix) = find(f, op.coeffs)
ctranspose(op::DiracMatrix) = DiracMatrix(op.coeffs', op.colbasis', op.rowbasis')
getindex(op::DiracMatrix, x...) = op.coeffs[x...]
setindex!(op::DiracMatrix, y, x) = setindex!(op.coeffs,y,x)

getpos(op::DiracMatrix, k::State{Ket}, b::State{Bra}) = (get(op.rowbasis, k), get(op.colbasis, b))
get(op::DiracMatrix, s::State{Ket}) = DiracVector(op[get(op.rowbasis, s), :], op.colbasis)
get(op::DiracMatrix, s::State{Bra}) = DiracVector(op[:, get(op.colbasis, s)], op.rowbasis)
get(op::DiracMatrix, k::State{Ket}, b::State{Bra}) = op[get(op.rowbasis, k), get(op.colbasis, b)]

function get(op::DiracMatrix, s::State, notfound)
	try
		return get(op, s)
	catch
		return notfound
	end
end

function get(op::DiracMatrix, k::State{Ket}, b::State{Bra}, notfound)
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
	@eval ($op)(n, d::DiracMatrix) = DiracVector(($op)(n,d.coeffs), d.rowbasis, d.colbasis)
	@eval ($op)(d::DiracMatrix, n) = DiracVector(($op)(d.coeffs,n), d.rowbasis, d.colbasis)
end

+(a::DiracMatrix, b::DiracMatrix) = samebasis(a,b) ? DiracMatrix(a.coeffs+b.coeffs, a.rowbasis, a.colbasis) : error("BasesMismatch")
-(a::DiracMatrix, b::DiracMatrix) = samebasis(a,b) ? DiracMatrix(a.coeffs-b.coeffs, a.rowbasis, a.colbasis) : error("BasesMismatch")
/(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs/d, op.rowbasis, op.colbasis)

*(a::DiracMatrix, b::DiracMatrix) = a.colbasis'==b.rowbasis ? DiracMatrix(a.coeffs*b.coeffs, a.rowbasis, b.colbasis) : error("BasesMismatch")

*(op::DiracMatrix, d::DiracCoeff) = DiracMatrix(op.coeffs*d, op.rowbasis, op.colbasis)
*(d::DiracCoeff, op::DiracMatrix) = DiracMatrix(d*op.coeffs, op.rowbasis, op.colbasis)
function *(op::DiracMatrix, s::State{Ket}) 
	if label(op.colbasis)==basislabel(s) 
		return get(op, s')
	else
		return reduce(+, [op.rowbasis[i]*((op.colbasis[j]*s)*op[i,j]) for i=1:length(op.rowbasis), j=1:length(op.colbasis)])
	end
end
function *(s::State{Bra}, op::DiracMatrix)
	if label(op.rowbasis)==basislabel(s) 
		return get(op, s')
	else
		return reduce(+, [op.colbasis[i]*((s*op.rowbasis[j])*op[i,j]) for i=1:length(op.rowbasis), j=1:length(op.colbasis)])
	end
end

*{T}(op::DiracMatrix, d::DiracVector{T, Ket}) = op.rowbasis == d.basis ? DiracVector(op.coeffs*d.coeffs, op.rowbasis) : error("BasesMismatch")
*{T}(d::DiracVector{T, Bra}, op::DiracMatrix) = op.colbasis == d.basis ? DiracVector(d.coeffs*op.coeffs, op.colbasis) : error("BasesMismatch")

trace(op::DiracMatrix) = trace(op.coeffs)

function ptrace(op::DiracMatrix, ind::Int)
	if isequal(op.colbasis, op.rowbasis')
		trrow = tensor(vcat(separate(op.rowbasis)[1:ind-1], separate(op.rowbasis)[ind+1:end])...)
		trcol = trrow'
		len = length(trcol)
		coeffs = [sum([op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.rowbasis)[ind])]) for j=1:length(trrow), k=1:length(trcol)] 
	else
		error("BasesMismatch")
	end
	return DiracMatrix(coeffs, trrow, trcol)
end
