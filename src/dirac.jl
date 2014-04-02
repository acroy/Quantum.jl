module m
include("rep.jl")
nfilled(S::SparseMatrixCSC) = int(S.colptr[end]-1)

const lang = "\u27E8"
const rang = "\u27E9"
const vert_ell = "\u22EE"
const horiz_ell = "\u2026"
using Base.print
import Base.show,
	   Base.getindex,
	   Base.ndims,
	   Base.size,
	   Base.length,
	   Base.slice,
	   Base.+,
	   Base.-,
	   Base.*,
	   Base.in,
	   Base.setdiff,
	   Base.get,
	   Base.!
#####Utility##########################################################
function todict(A::Array; d=Dict{eltype(A),Int}())
	sizehint(d, size(A,1))
	for i=1:size(A,1)
		d[A[i,:]...] = i
	end
	return d
end	

function crossjoin(A::Array, B::Vector)
    r1, r2 = size(A, 1), size(B, 1)
    columns = [[rep(A[:,c], 1, r2) for c=1:size(A,2)],
               [rep(B[:,c], r1, 1) for c=1:size(B,2)]]
    hcat(columns...)
end

function crossjoin(A::Array, B::Array)
	result = A;
	for i=1:length(B[1,:])
		result = crossjoin(result,B[:,i])
	end
	return result
end

tensor()=error("tensor() needs an argument of the form A::Array...")

function tensor(arr::Array...)
	if length(arr) == 2
		return crossjoin(arr[1], arr[2])
	else
		tensor(crossjoin(arr[1], arr[2]), arr[3:end]...)
	end
end

#####Basis##########################################################
abstract AbstractBasis

immutable Basis{T} <: AbstractBasis
	name::String
	bra_sym::String
	ket_sym::String
	label_map::Dict{T, Int}
	label_arr::Vector{T}
end

	#label_arr is of the form [label1,label2...]
	function Basis(name::String, label_arr::Vector; bra_sym=lang, ket_sym=rang)
		label_arr = unique(label_arr)
		label_map = todict(label_arr)
		Basis(name, bra_sym, ket_sym, label_map, label_arr)
	end

immutable TensorBasis{B<:Basis} <: AbstractBasis
	name::String
	label_map::Dict{Tuple, Int}
	label_arr::Array #memoization for future tensor products
	basis_arr::Vector{B}
	bra_sym::String
	ket_sym::String
end

	function TensorBasis{B<:Basis}(name::String, label_arr::Array, basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		label_map = todict(label_arr, d=Dict{Tuple, Int}())
		TensorBasis(name, label_map, label_arr, basis_arr, bra_sym, ket_sym)
	end

	function TensorBasis{B<:Basis}(name::String, basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		label_arr = tensor([i.label_arr for i in basis_arr]...)
		TensorBasis(name, label_arr, basis_arr, bra_sym=bra_sym, ket_sym=ket_sym)
	end

	function TensorBasis{B<:Basis}(basis_arr::Vector{B}; bra_sym=lang, ket_sym=rang)
		name = string(["_$(i.name)" for i in basis_arr]...)
		TensorBasis(name[2:end], basis_arr, bra_sym=bra_sym, ket_sym=ket_sym)
	end

	function tensor(bases::AbstractBasis...)
		TensorBasis(vcat([components(i) for i in bases]...))
	end

#####BasisFunctions##########################################################
function components(B::TensorBasis)
	return B.basis_arr
end

function components(B::Basis)
	return B
end

function extract(f::Function, A::Array)
    res = Array(eltype(A), 0)
    for i=1:size(A,1)
	    if f(A[i, :])
            push!(res, A[i,:])
        end
    end
    return vcat(res...)
end

function extract(f::Function, B::Basis; name=B.name)
	sub_arr = filter(f, B.label_arr)
	return Basis(name, sub_arr)
end

function extract(f::Function, B::TensorBasis; name=B.name)
	sub_arr = extract(f, B.label_arr)
	return TensorBasis(name, sub_arr, B.basis_arr, bra_sym=B.bra_sym, ket_sym=B.ket_sym)
end

######State##########################################################
abstract AbstractState
abstract BraKet
abstract Bra <: BraKet
abstract Ket <: BraKet

type State{B<:AbstractBasis, K<:BraKet} <: AbstractState
	kind::Type{K}
	label::String
	coeffs::SparseMatrixCSC
	basis::B
	function State(kind::Type{K}, label::String, coeffs::SparseMatrixCSC, basis)
		if length(basis)==length(coeffs)
			new(kind, label, coeffs, basis)
		elseif length(basis)>length(coeffs)
			error("coefficiemnts unspecified for $(length(basis)-length(coeffs)) basis states")
		else
			error("basis labels unspecified for $(length(coeffs)-length(basis)) coefficients")
		end	
	end	
end

State{B<:AbstractBasis,N<:Number}(label::String, coeffs::Vector{N}, basis::B)= State{B, Ket}(Ket, label, sparsevec(coeffs), basis)
State{B<:AbstractBasis,N<:Number}(label::String, coeffs::Array{N}, basis::B)= State{B, Bra}(Bra, label, sparse(coeffs), basis)

function State{B<:AbstractBasis}(label::String, coeffs::SparseMatrixCSC, basis::B)
	if size(coeffs)[2]==1
		return State{B, Ket}(Ket, label, coeffs, basis)
	else
		return State{B, Bra}(Bra, label, coeffs, basis)
	end
end

#####StateFunctions##########################################################

function magnitude(A::Number...)
	if length(A)==2
		return hypot(A[1], A[2])
	end
	return magnitude(hypot(A[1], A[2]), A[3:end]...)
end

magnitude{N<:Number}(A::Vector{N}) = magnitude(A...)
magnitude(A::SparseMatrixCSC) = magnitude(A...)
magnitude(S::AbstractState) = magnitude(S.coeffs)


function normalize{N<:Number}(A::Vector{N})
	return (1/magnitude(A))*A
end

function normalize!(S::AbstractState)
	S.coeffs = (1/magnitude(S))*S.coeffs;
	return S
end

isnorm(S::AbstractState) = magnitude(S)==1

######Operator######################################################

abstract AbstractOperator

type Operator{B1<:AbstractBasis, B2<:AbstractBasis} <:AbstractOperator
	label::String
	coeffs::SparseMatrixCSC
	row_basis::B1
	col_basis::B2
	function Operator(label::String, coeffs::SparseMatrixCSC, row_basis::B1, col_basis::B2)
		if (length(row_basis)*length(col_basis))==length(coeffs)
			new(label, coeffs, row_basis, col_basis)
		elseif length(basis)>length(coeffs)
			error("basis is larger than representation by $((length(row_basis)*length(col_basis))-length(coeffs)) labels")
		else
			error("representation is larger than basis by $(length(coeffs)-(length(row_basis)*length(col_basis))) coefficients")
		end
	end	
end



function Operator{B1<:AbstractBasis, B2<:AbstractBasis}(label::String, coeffs::SparseMatrixCSC, row_basis::B1, col_basis::B2)
	return Operator{B1, B2}(label, coeffs, row_basis, col_basis)
end

function Operator{B<:AbstractBasis}(label::String, coeffs::SparseMatrixCSC, both_basis::B)
	return Operator{B, B}(label, sparse(coeffs), both_basis, both_basis)
end

function Operator{B1<:AbstractBasis, B2<:AbstractBasis, N<:Number}(label::String, coeffs::Array{N, 2}, row_basis::B1, col_basis::B2)
	return Operator{B1, B2}(label, sparse(coeffs), row_basis, col_basis)
end

function Operator{B<:AbstractBasis, N<:Number}(label::String, coeffs::Array{N, 2}, both_basis::B)
	return Operator(label, coeffs, both_basis, both_basis)
end


#####FunctionOverloading############################################
function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) \"$(b.name)\"")
	println("$(length(b.label_arr[:,1])) Basis States:")
	for i=1:length(b.label_arr[:,1])
		println("| $(repr([b.label_arr[i, :]...])[2:end-1]) $(b.ket_sym)")
	end
end

function show(io::IO, s::AbstractState)
	if s.kind == Ket
		left = '|'
		right = s.basis.ket_sym
	else
		left = s.basis.bra_sym
		right = '|'
	end 
	println("$(typeof(s)) $left $(s.label) $right:")
	for i in find(s.coeffs)
		println("$(s.coeffs[i])	$left $(repr([s.basis.label_arr[i, :]...])[2:end-1]) $right")
	end
end

function show(io::IO, op::AbstractOperator)
	println("$(typeof(op)) \"$(op.label)\":")
	table = cell(length(op.row_basis)+1, length(op.col_basis)+1)	
	for i = 1:length(op.row_basis)
		table[i+1,1] = replace("| $(repr([op.row_basis.label_arr[i, :]...])[2:end-1]) $(op.row_basis.ket_sym)",r"\\", "")
	end
	for j = 1:length(op.col_basis)
		table[1,j+1] = replace("$(op.col_basis.bra_sym) $(repr([op.col_basis.label_arr[j, :]...])[2:end-1]) |",r"\\", "")
	end

	indent = ""
	for i=1:length(table[2])
		indent = " $indent"
	end	
	table[1,1] = indent
	table[2:end, 2:end] = full(op.coeffs)
	
	io = IOBuffer()
	show(io, table)
	io_str = takebuf_string(io)
	io_str = io_str[searchindex(io_str, "\n")+3:end]


	r_ket = "\"(?=$(op.col_basis.bra_sym))|(?<=$(op.row_basis.ket_sym))\""
	r = Regex(r_ket)
	io_str = replace(io_str,r"\"(?=\|)|\|\K\"",  " ")
	io_str = replace(io_str,r,  " ")
	io_str = replace(io_str,r"\\(?=\")",  " ")
	io_str = replace(io_str,r"\s\"\s|\s\"(?=,)",  "\" ")
	io_str = replace(io_str,r"\s\K\"(?=\s)", " ")
	print(io_str)
end

size(b::AbstractBasis, x::Int...) = size(b.label_arr, x...)
length(b::AbstractBasis) = size(b, 1)
ndims(b::AbstractBasis) = ndims(b.label_arr)

in(a, b::Basis)=in(a, collect(keys(b.label_map)))
in(a, b::TensorBasis)=in(tuple(a...), collect(keys(b.label_map)))

getindex(b::TensorBasis, x::Int) = TensorBasis("$(b.name)_$x", b.label_arr[x,:], b.basis_arr, bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::TensorBasis, x::Range1{Int}) = TensorBasis("$(b.name)_$(x[1]) to $(b.name)_$(last(x))", b.label_arr[x,:], b.basis_arr, bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::Basis, x::Int) = Basis("$(b.name)_$x", vec(b.label_arr[x,:]), bra_sym=b.bra_sym, ket_sym=b.ket_sym)
getindex(b::Basis, x::Range1{Int}) = Basis("$(b.name)_$(x[1]) to $(b.name)_$(last(x))", vec(b.label_arr[x,:]), bra_sym=b.bra_sym, ket_sym=b.ket_sym)

*(a::AbstractBasis, b::AbstractBasis) = tensor(a,b)
+(a::Basis,b::Basis) = Basis("$(a.name)+$(b.name)", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym)
+(a::TensorBasis,b::TensorBasis) = size(a,2)==size(b,2) ? TensorBasis("$(a.name)+$(b.name)", vcat(a.label_arr,b.label_arr), bra_sym=a.bra_sym, ket_sym=a.ket_sym) : error("dimension mismatch; label lengths differ")
-(a::Basis,b::Basis) = extract(x->!in(x,b), a, name="$(a.name)-$(b.name)")
-(a::TensorBasis, b::TensorBasis) = size(a,2)==size(b,2) ? extract(x->!in(x,b),a, name="$(a.name)-$(b.name)") : error("dimension mismatch; label lengths differ")
setdiff(a::AbstractBasis,b::AbstractBasis) = a-b

length(S::AbstractState) = nfilled(S.coeffs)
size(S::AbstractState) = length(S)

#The following redundancies are implemented for the sake 
#of allowing a[(key)] notation for state coefficient retrieval.
#It arises from the fact that the dicts for Basis and TensorBasis
#have keys of types T and Tuple respectively, rather than
#both having keys of type Tuple.

function get(b::AbstractBasis, key)
	v = get(b.label_map, key, "not found")
	if v=="not found"
		throw(KeyError(key))
	end
	return v
end

function get(b::AbstractBasis, key...)
	v = get(b.label_map, key, "not found")
	if v=="not found"
		throw(KeyError(key))
	end
	return v
end

get(b::TensorBasis, key::Array) = get(b, key...)
get{B<:Basis}(S::State{B}, key) = S.coeffs[get(S.basis, key)]
get{B<:TensorBasis}(S::State{B}, key) = S.coeffs[get(S.basis, key)]
get{B<:TensorBasis}(S::State{B}, key...) = S.coeffs[get(S.basis, key)]
getindex(S::State, x::Int) = S.coeffs[x]
getindex(S::State, x::Range1{Int}) = S.coeffs[x,:]
getindex{B<:Basis}(S::State{B}, x::Tuple) = get(S, x[1])
getindex{B<:TensorBasis}(S::State{B}, x::Tuple) = get(S, x)
!(K::Type{Ket}) = Bra
!(B::Type{Bra}) = Ket

function transpose(S::State)
	return State(S.label, S.coeffs.', S.basis)

end

function ctranspose(S::State)
	return State(S.label, S.coeffs', S.basis)
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::State{B1, Bra}, b::State{B2, Ket})
	return (a.coeffs*b.coeffs)[1]
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::State{B1, Ket}, b::State{B2, Ket})
	new_coeffs = tensor(full(a.coeffs), full(b.coeffs))
	return State("$(a.label),$(b.label)", new_coeffs[:,1].*new_coeffs[:,2], a.basis*b.basis)
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::State{B1, Bra}, b::State{B2, Bra})
	new_coeffs = tensor(full(a.coeffs)', full(b.coeffs)')'
	return State("$(a.label),$(b.label)", new_coeffs[1,:].*new_coeffs[2,:], a.basis*b.basis)
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::State{B1, Ket}, b::State{B2, Bra})
	return Operator("|$(a.label) $(a.basis.ket_sym) $(b.basis.bra_sym) $(b.label)|", a.coeffs*b.coeffs, a.basis, b.basis)
end

end #module
#####Tests##########################################################
using m

function main()
	A = m.Basis("A", [1:3])
	B = m.Basis("B", ["$i" for i=1:3])
	C = m.Basis("C", ["%", "\$", "#"])
	D = m.tensor(A, B, C)
	return (A,B,C,D)
end

res = main();
a = res[1]
b = res[2]
c = res[3]
d = res[4]
ds = m.State("dState", [1:27], d);
a = m.Basis("a", [1:10], ket_sym="}")
as = m.State("as", [1:10], a)
println("")

