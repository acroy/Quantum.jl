#####################################
#DiracMatrix#########################
#####################################

type DiracMatrix{K<:Ket, B<:Bra, T} <: AbstractOperator
	coeffs::SparseMatrixCSC{T, Int}
	rowb::AbstractBasis{K}
	colb::AbstractBasis{B}
	function DiracMatrix(coeffs, rowb, colb)
		@assert size(coeffs)==(length(rowb), length(colb)) "Dimensions of coefficient array do not match bases"
		return new(coeffs, rowb, colb)
	end
end
DiracMatrix{K<:Ket, B<:Bra, T}(coeffs::SparseMatrixCSC{T}, rowb::AbstractBasis{K}, colb::AbstractBasis{B}) = DiracMatrix{K,B,T}(coeffs, rowb, colb) 
DiracMatrix{K<:Ket, B<:Bra, T}(coeffs::Array{T}, rowb::AbstractBasis{K}, colb::AbstractBasis{B}) = DiracMatrix{K,B,T}(sparse(coeffs), rowb, colb) 
DiracMatrix{K<:Ket}(coeffs::AbstractArray, b::AbstractBasis{K}) = dmat(coeffs, b, b') 
DiracMatrix{B<:Bra}(coeffs::AbstractArray, b::AbstractBasis{B}) = dmat(coeffs, b', b) 

function constructop!(coeffs, fcoeff, fstate, b) #basically never use this function except for in the below context
	@inbounds begin
		for i=1:length(b), j=1:length(b)
			coeffs[i,j] = fcoeff(b[j])*inner(b[i]',fstate(b[j]))
		end
	end
end

function DiracMatrix{K<:Ket}(fcoeff::Function, fstate::Function, b::AbstractBasis{K}, t::DataType=typeof(fcoeff(b[1])))
	coeffs = ones(t,length(b),length(b))
	constructop!(coeffs, fcoeff, fstate, b)
	return dmat(coeffs, b)
end

dmat=DiracMatrix

copy(dm::DiracMatrix) = dmat(copy(dm.coeffs), copy(dm.rowb), copy(dm.colb))

#####################################
#Getter-style Functions##############
#####################################
bsym(dm::DiracMatrix) = [bsym(dm.rowb), bsym(dm.colb)]

#####################################
#Boolean Functions###################
#####################################
isequal(a::DiracMatrix, b::DiracMatrix) = isequal(a.coeffs,b.coeffs) && a.rowb==b.rowb && a.colb==b.colb 
==(a::DiracMatrix, b::DiracMatrix) = a.coeffs==b.coeffs && a.rowb==b.rowb && a.colb==b.colb 
isdual(a::DiracMatrix, b::DiracMatrix) = a.coeffs'==b.coeffs && isdual(a.rowb,b.colb) && isdual(b.rowb,a.colb)
#####################################
#Array-like Functions################
#####################################
size(dm::DiracMatrix, args...) = size(dm.coeffs, args...)
for op=(:endof, :eltype, :length, :find, :findn, :findnz, :nnz, :ndims, :countnz)
	@eval ($op)(dm::DiracMatrix) = ($op)(dm.coeffs)
end
ctranspose(dm::DiracMatrix) = dmat(dm.coeffs', dm.colb', dm.rowb')
getindex(dm::DiracMatrix, x...) = dm.coeffs[x...]
setindex!(dm::DiracMatrix, y, x...) = setindex!(dm.coeffs,y,x...)

#####################################
#Dict-like Functions#################
#####################################
getpos{K<:Ket}(dm::DiracMatrix{K}, k::State{K}) = get(dm.rowb, k)
getpos{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, b::State{B}) = get(dm.colb, b)
getpos{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, k::State{K}, b::State{B}) = (getpos(dm, k), getpos(dm, b))
getpos{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, o::OuterProduct{K,B}) = getpos(dm, o.ket, o.bra)
getpos(dm::DiracMatrix, arg) = throw(KeyError(arg)) 
getpos(dm::DiracMatrix, args...) = throw(KeyError(args)) 

get{K<:Ket}(dm::DiracMatrix{K}, s::State{K}) = dvec(dm[getpos(dm, s), :], dm.colb)
get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, s::State{B}) = dvec(dm[:, getpos(dm, s)], dm.rowb)
get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, k::State{K}, b::State{B}) = dm[getpos(dm, k), getpos(dm, b)]
get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, o::OuterProduct{K,B}) = get(dm, o.ket, o.bra)

get{K<:Ket}(dm::DiracMatrix{K}, s::State{K}, notfound) = in(k, dm.rowb) ? get(dm, k) : notfound
get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, s::State{B}, notfound) = in(b, dm.colb) ? get(dm, b) : notfound
function get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, k::State{K}, b::State{B}, notfound)
	if in(k, dm.rowb) && in(b, dm.colb)
		return get(dm, k, b)
	else
		return notfound
	end	
end
get{K<:Ket, B<:Bra}(dm::DiracMatrix{K,B}, o::OuterProduct{K,B}, notfound) = get(dm, o.ket, o.bra, notfound)

get(dm::DiracMatrix, arg) = throw(KeyError(arg)) 
get(dm::DiracMatrix, args...) = throw(KeyError(args)) 

#####################################
#Show Functions######################
#####################################
summary(dm::DiracMatrix) = "$(size(dm,1))x$(size(dm,2)) $(typeof(dm))"

function showcompact(io::IO, dm::DiracMatrix)
	if length(dm.coeffs)==0
		print(io, "$(typeof(dm))[]")
	else
		tempio = IOBuffer()
		nz = hcat(findn(dm)...)
		print(tempio, [" + $(dm.coeffs[nz[i,1],nz[i,2]])$(dm.rowb[nz[i,1]])$(dm.colb[nz[i,2]])" for i=1:size(nz,1)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end

function show(io::IO, dm::DiracMatrix)
	println(io, summary(dm))
	table = cell(length(dm.rowb)+1, length(dm.colb)+1)	
	for i = 1:length(dm.rowb)
		table[i+1,1] = dm.rowb[i]
	end
	for j = 1:length(dm.colb)
		table[1,j+1] = dm.colb[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = full(dm.coeffs)
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	print(io, io_str[searchindex(io_str, "\n")+3:end])
end

#####################################
#Function-passing Functions##########
#####################################
find(f::Function, dm::DiracMatrix) = find(f, dm.coeffs)
findstates(f::Function, dm::DiracMatrix) = find(f, [dm.rowb[i]*dm.colb[j] for i=1:length(dm.rowb), j=1:length(dm.colb)]) #f takes OuterProduct as argument
map(f::Function, dm::DiracMatrix) = dmat(map(f, dm.coeffs), dm.rowb, dm.colb)
qeval(f::Function, dm::DiracMatrix) = map(x->qeval(f, x), dm)

#####################################
#Arithmetic Operations###############
#####################################
for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracMatrix, b::DiracVector) = dmat(($op)(a.coeffs,b.coeffs), a.rowb, a.colb)
	@eval ($op)(a::DiracVector, b::DiracMatrix) = dmat(($op)(a.coeffs,b.coeffs), b.rowb, b.colb)
	@eval ($op)(a::DiracMatrix, b::DiracMatrix) = dmat(($op)(a.coeffs,b.coeffs), a.rowb, a.colb)
	@eval ($op)(n::DiracCoeff, d::DiracMatrix) = dmat(($op)(n,d.coeffs), d.rowb, d.colb)
	@eval ($op)(d::DiracMatrix, n::DiracCoeff) = dmat(($op)(d.coeffs,n), d.rowb, d.colb)
end

/(dm::DiracMatrix, d::DiracCoeff) = dmat(dm.coeffs/d, dm.rowb, dm.colb)
*(dm::DiracMatrix, d::DiracCoeff) = dmat(dm.coeffs*d, dm.rowb, dm.colb)
*(d::DiracCoeff, dm::DiracMatrix) = dmat(d*dm.coeffs, dm.rowb, dm.colb)

kron(c::DiracCoeff, dm::DiracMatrix) = c*dm
kron(dm::DiracMatrix, c::DiracCoeff) = dm*c

function multbystate!{K<:Ket}(dm::DiracMatrix, s::State{K}, dv::AbstractArray) 
	for i in findn(dm)[1]
		coeff = 0
		for j in findn(dm)[2]
			coeff+=dm[i,j]*(dm.colb[j]*s)
		end
		dv[i]=dv[i]+coeff
	end
	return dvec(dv, dm.rowb)
end

function multbystate!{B<:Bra}(dm::DiracMatrix, s::State{B}, dv::AbstractArray)
	for j in findn(dm)[2]
		coeff = 0
		for i in findn(dm)[1]
			coeff+=dm[i,j]*(s*dm.rowb[i])
		end
		dv[j]=dv[j]+coeff
	end
	return dvec(dv, dm.colb)
end

consmult{K<:Ket}(dm::DiracMatrix, s::State{K}) = multbystate!(dm, s, samebasis(dm.colb,s) ? zeros(length(dm.rowb)) : zeros(Any,length(dm.rowb)))
consmult{B<:Bra}(dm::DiracMatrix, s::State{B}) = multbystate!(dm, s, samebasis(dm.rowb,s) ? zeros(1,length(dm.colb)) : zeros(Any,1,length(dm.colb)))

#currently the below consmults are slow and 
#should be replaced by a something 
#more in line with the above single
#state versions that preallocate memory
#more effectively
consmult{K<:Ket}(dm::DiracMatrix, d::DiracVector{K}) = reduce(+,[dm.rowb[i]*reduce(+,[dm[i,j]*(dm.colb[j]*d) for j=1:length(dm.colb)]) for i=1:length(dm.rowb)])
consmult{B<:Bra}(dm::DiracMatrix, d::DiracVector{B}) = reduce(+,[dm.colb[j]*reduce(+,[dm[i,j]*(d*op.rowb[i]) for j=1:length(dm.rowb)]) for i=1:length(dm.colb)])
consmult(a::DiracMatrix, b::DiracMatrix) = reduce(+,[(a[i,j]*b[m,n]*(a.colb[j]*b.rowb[m]))*(a.rowb[i]*b.colb[n]) for n=1:size(b,2), m=1:size(b,1), j=1:size(a,2), i=1:size(a,1)])

*{K<:Ket}(dm::DiracMatrix, s::State{K}) = in(s', dm.colb) ? get(dm, s') : consmult(dm,s)
*{B<:Bra}(s::State{B}, dm::DiracMatrix) = in(s', dm.rowb) ? get(dm, s') : consmult(dm,s)
*{K<:Ket}(dm::DiracMatrix, dv::DiracVector{K}) = isdual(dm.colb, dv.basis) ? dvec(dm.coeffs*dv.coeffs, dm.rowb) : consmult(dm,dv)
*{B<:Bra}(dv::DiracVector{B},dm::DiracMatrix) = isdual(dm.rowb, dv.basis) ? dvec(dm.coeffs*dv.coeffs, dm.colb) : consmult(dm,dv)
*(a::DiracMatrix, b::DiracMatrix) = isdual(a.colb, b.rowb) ? dmat(a.coeffs*b.coeffs, a.rowb, b.colb) : consmult(a,b)

kron(a::DiracMatrix, b::DiracMatrix) = dmat(kron(a.coeffs, b.coeffs), tensor(a.rowb, b.rowb), tensor(a.colb, b.colb)) 
kron{K<:Ket}(dm::DiracMatrix, dv::DiracVector{K}) = dmat(kron(dm.coeffs, dv.coeffs), tensor(dm.rowb, dv.basis), dm.colb)
kron{K<:Ket}(dv::DiracVector{K}, dm::DiracMatrix) = dmat(kron(dv.coeffs, dm.coeffs), tensor(dv.basis, dm.rowb), dm.colb)
kron{B<:Bra}(dm::DiracMatrix, dv::DiracVector{B}) = dmat(kron(dm.coeffs, dv.coeffs), dm.rowb, tensor(dm.colb, dv.basis))
kron{B<:Bra}(dv::DiracVector{B}, dm::DiracMatrix) = dmat(kron(dv.coeffs, dm.coeffs), dm.rowb, tensor(dv.basis, dm.colb))
kron{K<:Ket}(s::State{K}, dm::DiracMatrix) = dmat(dm.coeffs, tensor(s, dm.rowb), dm.colb)
kron{K<:Ket}(dm::DiracMatrix, s::State{K}) = dmat(dm.coeffs, tensor(s, dm.rowb), dm.colb)
kron{B<:Bra}(s::State{B}, dm::DiracMatrix) = dmat(dm.coeffs, dm.rowb, tensor(s, dm.colb))
kron{B<:Bra}(dm::DiracMatrix, s::State{B}) = dmat(dm.coeffs, dm.rowb, tensor(s, dm.colb))
kron(o::OuterProduct, dm::DiracMatrix) = dmat(dm.coeffs, tensor(o.ket, dm.rowb), tensor(o.bra, dm.colb)) 
kron(dm::DiracMatrix, o::OuterProduct) = dmat(dm.coeffs, tensor(dm.rowb, o.ket), tensor(dm.colb, o.bra)) 


function add_out!(dm, o)
	dm[getpos(dm, o)...] = 1+get(dm, o)
	return dm
end

add_out(dm, o) = add_out!(copy(dm), o)
add_out{K,B,T<:InnerProduct}(dm::DiracMatrix{K,B,T}, o) = add_out!(1*dm, o) #converts res to DiracMatrix{K,B,ScalarExpr}

function add_bra(dm::DiracMatrix, o::OuterProduct)
	newcol = zeros(size(dm,1))
	newcol[get(dm.rowb, o.ket)] = 1 
	return dmat(hcat(dm.coeffs,newcol), dm.rowb, bcat(dm.colb, o.bra))
end

function add_bra(o::OuterProduct, dm::DiracMatrix)
	newcol = zeros(size(dm,1))
	newcol[get(dm.rowb, o.ket)] = 1 
	return dmat(hcat(newcol, dm.coeffs), dm.rowb, bcat(o.bra, dm.colb))
end

function add_ket(dm::DiracMatrix, o::OuterProduct)
	newrow = zeros(1,size(dm,2))
	newrow[get(dm.colb, o.bra)] = 1 
	return dmat(vcat(dm.coeffs,newrow), bcat(dm.rowb, o.ket), dm.colb)
end

function add_ket(o::OuterProduct, dm::DiracMatrix)
	newrow = zeros(1,size(dm,2))
	newrow[get(dm.colb, o.bra)] = 1 
	return dmat(vcat(newrow,dm.coeffs), bcat(o.ket,dm.rowb), dm.colb)
end

function add_both!(dm::DiracMatrix, o::OuterProduct)
	dm[getpos(dm, o)...] = 1+get(dm, o)
	return dm
end
function add_both!(dm::DiracMatrix, o::OuterProduct,rb,cb,res)
	res[1:size(dm,1), 1:size(dm,2)] = dm.coeffs
	return add_both!(dmat(res, rb, cb), o)
end

function add_both!(o::OuterProduct,dm::DiracMatrix,rb,cb,res)
	res[2:size(dm,1)+1, 2:size(dm,2)+1] = dm.coeffs
	return add_both!(dmat(res, rb, cb), o)
end

add_both(dm::DiracMatrix, o::OuterProduct, rb, cb) = add_both!(dm, o, rb, cb, zeros(length(rb), length(cb)))
add_both{K,B,T}(dm::DiracMatrix{K,B,T}, o::OuterProduct) = add_both(dm, o, bcat(dm.rowb, o.ket), bcat(dm.colb, o.bra))
add_both{K,B,T<:InnerProduct}(dm::DiracMatrix{K,B,T}, o::OuterProduct) = add_both(1*dm, o)

add_both(o::OuterProduct, dm::DiracMatrix, rb, cb) = add_both!(o, dm, rb, cb, zeros(length(rb), length(cb)))
add_both{K,B,T}(o::OuterProduct, dm::DiracMatrix{K,B,T}) = add_both(o, dm, bcat(o.ket, dm.rowb), bcat(o.bra, dm.colb))
add_both{K,B,T<:InnerProduct}(o::OuterProduct, dm::DiracMatrix{K,B,T}) = add_both(o,1*dm)

function +{K<:Ket, B<:Bra}(op::DiracMatrix{K,B}, o::OuterProduct{K,B})
	if in(o.ket, op.rowb)
		if in(o.bra, op.colb)
			return add_out(op, o)
		else
			return add_bra(op, o)
		end
	elseif in(o.bra, op.colb)
		return add_ket(op, o)
	else
		return add_both(op,o)
	end
end

function +{K<:Ket, B<:Bra}(o::OuterProduct{K,B}, op::DiracMatrix{K,B})
	if in(o.ket, op.rowb)
		if in(o.bra, op.colb)
			return add_out(op,o)
		else
			return add_bra(o,op)
		end
	elseif in(o.bra, op.colb)
		return add_ket(o,op)
	else
		return add_both(o,op)
	end
end

function add_dmat!(a::DiracMatrix, b::DiracMatrix)
	for i=1:findn(b)[1]
		for j=1:findn(b)[2]
			a=a+(b.rowb[i]*b.colb[j])
			a[getpos(a, b.rowb[i], b.colb[j])...] = get(a, b.rowb[i], b.colb[j])+b[i,j]-1
		end
	end
	return a
end

function +{K<:Ket, B<:Bra}(a::DiracMatrix{K,B}, b::DiracMatrix{K,B})
	if a.rowb==b.rowb && a.colb==b.colb 
		return dmat(a.coeffs+b.coeffs, a.rowb, a.colb) 
	else
		return add_dmat!(copy(a), b)
	end
end

-(op::DiracMatrix) = -1*op
-(a::DiracMatrix, b::DiracMatrix) = a+(-b)
-(a::DiracMatrix, b::OuterProduct) = a+(-b)

exp(op::DiracMatrix) = dmat(exp(op.coeffs), op.rowb, op.colb)
^(op::DiracMatrix, i::Integer) = dmat(^(op.coeffs, i), op.rowb, op.colb)
^(op::DiracMatrix, n::Number) = dmat(^(op.coeffs, n), op.rowb, op.colb)
trace(op::DiracMatrix) = trace(op.coeffs)

commutator(a::DiracMatrix, b::DiracMatrix) = (a*b) - (b*a)

function ptrace(op::DiracMatrix, ind::Int)
	@assert isdual(op.colb, op.rowb) "BasisMismatch"
	trrow = tensor(vcat(separate(op.rowb)[1:ind-1], separate(op.rowb)[ind+1:end])...)
	trcol = trrow'
	len = length(trcol)
	#the line below should be optimized by pre-allocating the array and 
	#writing an in-place function to fill it
	coeffs = [reduce(+,[op[(((i-1)*len)+k), (((i-1)*len)+j)] for i=1:length(separate(op.rowb)[ind])]) for j=1:length(trrow), k=1:length(trcol)] 
	return dmat(vcat([hcat(coeffs[i, :]...) for i=1:size(coeffs, 1)]...), trrow, trcol) #use vcat()/hcat() trick to convert to most primitive common type
end

actop{K<:Ket}(op::AbstractOperator, s::Tensor{K}, i::Integer) = kron(vcat(s[1:i-1], op*s[i], s[i+1:end])...)
actop{B<:Bra}(op::AbstractOperator, s::Tensor{B}, i::Integer) = kron(vcat(s[1:i-1], s[i]*op, s[i+1:end])...)
