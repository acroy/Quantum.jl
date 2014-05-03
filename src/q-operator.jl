#####################################
#OperatorRep############################
#####################################
type OperatorRep <: Quantum
	coeffs::Matrix{Complex{Float64}}
	row_basis::AbstractBasis{Ket}
	col_basis::AbstractBasis{Bra}
	function OperatorRep(coeffs::Matrix{Complex{Float64}}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra})
		if size(coeffs)==(length(row_basis), length(col_basis))
			new(coeffs, row_basis, col_basis)
		else
			throw(DimensionMismatch)
		end
	end
end

OperatorRep{N<:Number}(coeffs::Matrix{N}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra}) = OperatorRep(convert(Matrix{Complex{Float64}}, coeffs), row_basis, col_basis) 
OperatorRep{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Ket}) = OperatorRep(convert(Matrix{Complex{Float64}}, coeffs), b, b') 
OperatorRep{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis{Bra}) = OperatorRep(convert(Matrix{Complex{Float64}}, coeffs), b', b) 

function OperatorRep(coeff_func::Function, label_func::Function, b::AbstractBasis)
	coeffs = zeros(Number, length(b), length(b))
	for i=1:length(b)
		for j=1:length(b)
			ii = get(b,State(label_func(b[i])), -1)
			if ii != -1
				coeffs[ii,j] = coeff_func(b[i]) * (b[j]'*b[i])
			end
		end
	end
	return OperatorRep(coeffs, b)
end

#####################################
#Functions###########################
#####################################

#imported############################
isequal(a::OperatorRep, b::OperatorRep) = a.coeffs==b.coeffs && samebasis(a,b)

copy(op::OperatorRep, coeffs=copy(op.coeffs)) = OperatorRep(coeffs, op.row_basis, op.col_basis)
expm(op::OperatorRep) = copy(op, expm(op.coeffs))

ndims(op::OperatorRep) = 2
size(op::OperatorRep) = size(op.coeffs)
length(op::OperatorRep) = length(op.coeffs)
endof(op::OperatorRep) = length(op)
find(op::OperatorRep) = find(op.coeffs)
find(f::Function, op::OperatorRep) = find(f::Function, op.coeffs)

size(op::OperatorRep, i::Integer) = size(op.coeffs, i)
ctranspose(op::OperatorRep) = OperatorRep(op.coeffs', op.col_basis',op.row_basis')
getindex(op::OperatorRep, x...) = op.coeffs[x...]
setindex!(op::OperatorRep, y, x) = setindex!(op.coeffs,y,x)
get(op::OperatorRep, ket_label, bra_label) = op[get(op.row_basis, ket_label), get(op.col_basis, bra_label)]
get(op::OperatorRep, s::State{Ket}) = op[get(op.row_basis, s), :]
get(op::OperatorRep, s::State{Bra}) = op[:, get(op.col_basis, s)]

.+(op::OperatorRep, n::Number) = copy(op, op.coeffs.+n)
.+(n::Number, op::OperatorRep) = copy(op, n.+op.coeffs)
+(arr::Array, op::OperatorRep) = copy(op, arr+op.coeffs)
+(op::OperatorRep, arr::Array) = copy(op, op.coeffs+arr)
.-(op::OperatorRep, n::Number) = copy(op, op.coeffs.-n)
.-(n::Number, op::OperatorRep) = copy(op, n.-op.coeffs)
-(arr::Array, op::OperatorRep) = copy(op, arr-op.coeffs)
-(op::OperatorRep, arr::Array) = copy(op, op.coeffs-arr)
.^(op::OperatorRep, n::Number) = copy(op, op.coeffs.^n)
.^(n::Number, op::OperatorRep) = copy(op, n.^op.coeffs)

+(a::OperatorRep, b::OperatorRep) = samebasis(a,b) ? copy(a, a.coeffs+b.coeffs) : error("BasesMismatch")
-(a::OperatorRep, b::OperatorRep) = samebasis(a,b) ? copy(a, a.coeffs-b.coeffs) : error("BasesMismatch")
/(op::OperatorRep, n::Number) = copy(op, op.coeffs/n)

*(a::OperatorRep, b::OperatorRep) = samebasis(a,b) ? copy(a, a.coeffs*b.coeffs) : error("BasesMismatch")
*(op::OperatorRep, n::Number) = copy(op, op.coeffs*n)
*(n::Number, op::OperatorRep) = copy(op, n*op.coeffs)
*(op::OperatorRep, s::State{Ket}) = StateRep(State([]), get(op, s'), op.row_basis)
*(s::State{Bra}, op::OperatorRep) = StateRep(State([], Bra), get(op, s'), op.col_basis)
*(op::OperatorRep, s::StateRep{Ket}) = op.row_basis == s.basis ? StateRep(s.state, op.coeffs*s.coeffs, op.row_basis) : error("BasesMismatch")
*(s::StateRep{Bra}, op::OperatorRep) = op.col_basis == s.basis ? StateRep(s.state, s.coeffs*op.coeffs, op.col_basis) : error("BasesMismatch")


*{N<:Number}(arr::Array{N, 2}, op::OperatorRep) = size(arr,1)==1 ? StateRep(State(), arr*op.coeffs, op.col_basis) : copy(op, arr*op.coeffs)
*{N<:Number}(op::OperatorRep, arr::Array{N, 1}) = StateRep(State([]), op.coeffs*arr, op.row_basis)
*{N<:Number}(op::OperatorRep, arr::Array{N, 2}) = size(arr,2)==1 ? StateRep(State(), op.coeffs*arr, op.row_basis) : copy(op,op.coeffs*arr)

trace(op::OperatorRep) = trace(op.coeffs)

function show(io::IO, op::OperatorRep)
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

samebasis(a::OperatorRep, b::OperatorRep) = isequal(a.row_basis,b.row_basis) && isequal(a.col_basis, b.col_basis)

function findstates(f::Function, op::OperatorRep)
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

commutator(a::OperatorRep, b::OperatorRep) = samebasis(a,b) ? copy(a, (a.coeffs*b.coeffs)-(b.coeffs*a.coeffs)) : error("BasesMismatch")

function ptrace(op::OperatorRep, ind::Int)
	if isequal(op.col_basis, op.row_basis')
		tr_row = tensor(vcat(separate(op.row_basis)[1:ind-1], separate(op.row_basis)[ind+1:end])...)
		tr_col = tr_row'
		len = length(tr_col)
		coeffs = Complex{Float64}[sum([op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.row_basis)[ind])]) for j=1:length(tr_row), k=1:length(tr_col)] 
	else
		error("BasesMismatch")
	end
	return OperatorRep(coeffs, tr_row, tr_col)
end
