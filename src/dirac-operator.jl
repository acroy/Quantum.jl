#####################################
#Operator############################
#####################################
type Operator <: Dirac
	coeffs::Matrix{Complex{Float64}}
	row_basis::AbstractBasis{Ket}
	col_basis::AbstractBasis{Bra}
	function Operator(coeffs::Matrix{Complex{Float64}}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra})
		if size(coeffs)==(length(row_basis), length(col_basis))
			new(coeffs, row_basis, col_basis)
		else
			throw(DimensionMismatch)
		end
	end
end

Operator{N<:Number}(coeffs::Matrix{N}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), row_basis, col_basis) 
Operator{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Ket}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), b, b') 
Operator{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Bra}) = Operator(convert(Matrix{Complex{Float64}}, coeffs), b', b) 

function Operator(label_func::Function, coeff_func::Function, b::AbstractBasis)
	coeffs = Array(Number, length(b), length(b))
	for i=1:length(b)
		for j=1:length(b)
			coeffs[i,j] = coeff_func(b[i]) * (b[j]'*label_func(b[i]))
		end
	end
	return Operator(coeffs, b)
end

#####################################
#Functions###########################
#####################################

#imported############################
isequal(a::Operator, b::Operator) = a.coeffs==b.coeffs && samebasis(a,b)

copy(op::Operator, coeffs=copy(op.coeffs)) = Operator(coeffs, op.row_basis, op.col_basis)
expm(op::Operator) = copy(op, expm(op.coeffs))

ndims(op::Operator) = 2
size(op::Operator) = size(op.coeffs)
length(op::Operator) = length(op.coeffs)
endof(op::Operator) = length(op)
find(op::Operator) = find(op.coeffs)
find(f::Function, op::Operator) = find(f::Function, op.coeffs)

size(op::Operator, i::Int) = size(op.coeffs, i)
ctranspose(op::Operator) = Operator(op.coeffs', op.col_basis',op.row_basis')
getindex(op::Operator, x...) = op.coeffs[x...]
setindex!(op::Operator, y, x) = setindex!(op.coeffs,y,x)
get(op::Operator, ket_label, bra_label) = op[get(op.row_basis, ket_label), get(op.col_basis, bra_label)]
get(op::Operator, s::State{Ket}) = op[get(op.row_basis, s), :]
get(op::Operator, s::State{Bra}) = op[:, get(op.col_basis, s)]

.+(op::Operator, n::Number) = copy(op, op.coeffs.+n)
.+(n::Number, op::Operator) = copy(op, n.+op.coeffs)
.-(op::Operator, n::Number) = copy(op, op.coeffs.-n)
.-(n::Number, op::Operator) = copy(op, n.-op.coeffs)
.^(op::Operator, n::Number) = copy(op, op.coeffs.^n)
.^(n::Number, op::Operator) = copy(op, n.^op.coeffs)

+(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs+b.coeffs) : error("BasesMismatch")
-(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs-b.coeffs) : error("BasesMismatch")
/(op::Operator, n::Number) = copy(op, op.coeffs/n)

*(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, a.coeffs*b.coeffs) : error("BasesMismatch")
*(op::Operator, n::Number) = copy(op, op.coeffs*n)
*(n::Number, op::Operator) = copy(op, n*op.coeffs)
*(op::Operator, s::State{Ket}) = StateRep(State([]), get(op, s'), op.row_basis)
*(s::State{Bra}, op::Operator) = StateRep(State([], Bra), get(op, s'), op.col_basis)
*(op::Operator, s::StateRep{Ket}) = op.row_basis == s.basis ? StateRep(s.state, op.coeffs*s.coeffs, op.row_basis) : error("BasesMismatch")
*(s::StateRep{Bra}, op::Operator) = op.col_basis == s.basis ? StateRep(s.state, s.coeffs*op.coeffs, op.col_basis) : error("BasesMismatch")
*(arr::Array, op::Operator) = copy(op, arr*op.coeffs)
*(op::Operator, arr::Array) = copy(op, op.coeffs*arr)

trace(op::Operator) = trace(op.coeffs)

function show(io::IO, op::Operator)
	println("$(typeof(op)):")
	table = cell(length(op.row_basis)+1, length(op.col_basis)+1)	
	for i = 1:length(op.row_basis)
		table[i+1,1] = op.row_basis[i]
	end
	for j = 1:length(op.col_basis)
		table[1,j+1] = op.col_basis[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = op.coeffs	
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	io_str = io_str[searchindex(io_str, "\n")+3:end]
	print(io_str)
end

#exported############################

samebasis(a::Operator, b::Operator) = isequal(a.row_basis,b.row_basis) && isequal(a.col_basis, b.col_basis)

function findstates(f::Function, op::Operator)
	numrows = size(op, 1)
	inds = find(f, op)
	result = Array(State, length(inds), 2)
	for i=1:length(inds)
		row_ind = inds[i]%numrows
		if row_ind == 0
			row_ind = numrows
		end
		col_ind = ((inds[i]-row_ind)/numrows)+1
		result[i, 2] = op.col_basis[col_ind]
		result[i, 1] = op.row_basis[row_ind]
	end
	return result
end

commutator(a::Operator, b::Operator) = samebasis(a,b) ? copy(a, (a.coeffs*b.coeffs)-(b.coeffs*a.coeffs)) : error("BasesMismatch")

function ptrace(op::Operator, ind::Int)
	if isequal(op.col_basis, op.row_basis')
		tr_row = tensor(vcat(separate(op.row_basis)[1:ind-1], separate(op.row_basis)[ind+1:end])...)
		tr_col = tr_row'
		len = length(tr_col)
		coeffs = Complex{Float64}[sum([op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.row_basis)[ind])]) for j=1:length(tr_row), k=1:length(tr_col)] 
	else
		error("BasesMismatch")
	end
	return Operator(coeffs, tr_row, tr_col)
end
