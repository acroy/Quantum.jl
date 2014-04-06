module m
include("rep.jl")

const lang = "\u27E8"
const rang = "\u27E9"
const vert_ell = "\u22EE"
const horiz_ell = "\u2026"
import Base.show,
	   Base.getindex,
	   Base.setindex!,
	   Base.ndims,
	   Base.size,
	   Base.length,
	   Base.slice,
	   Base.(.+),
	   Base.(.^),
	   Base.(.-),
	   Base.^,
	   Base.*,
	   Base.in,
	   Base.setdiff,
	   Base.get,
	   Base.!,
	   Base.exp,
	   Base.map,
	   Base.map!,
	   Base.filter,
	   Base.isequal,
	   Base.copy
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
######StateRep##########################################################
abstract AbstractState
abstract BraKet
abstract Bra <: BraKet
abstract Ket <: BraKet

type State{K<:BraKet}
	label:Vector
	kind::K
end

State(label::Vector) = State(label,Ket)
State(label...) = State([label...],Ket)

type StateRep{B<:AbstractBasis, K<:BraKet} <: AbstractState
	state::State{K}
	coeffs::Array{Number}
	basis::B
	function StateRep(kind::Type{K}, label::String, coeffs::Array{N}, basis)
		if length(basis)==length(coeffs)
			new(kind, label, coeffs, basis)
		elseif length(basis)>length(coeffs)
			error("coefficients unspecified for $(length(basis)-length(coeffs)) basis StateReps")
		else
			error("basis labels unspecified for $(length(coeffs)-length(basis)) coefficients")
		end	
	end	
end

StateRep{B<:AbstractBasis,N<:Number}(label::String, coeffs::Vector{N}, basis::B)=StateRep{B, Ket, N}(Ket, label, coeffs, basis)

function StateRep{B<:AbstractBasis,N<:Number}(label::String, coeffs::Array{N}, basis::B)
	if size(coeffs)[2]==1
		StateRep(label, vec(coeffs), basis)
	else
		StateRep{B, Bra, N}(Bra, label, coeffs, basis)
	end
end

#####StateRepFunctions##########################################################
getcoeffs(s::StateRep) = s.coeffs
getkind(s::StateRep) = s.state.kind
getlabel(s::StateRep) = s.state.label
diracform(s::StateRep) = getkind(S)==Ket ? "| $(S.label) $(S.basis.ket_sym)" : "$(S.basis.bra_sym) $(S.label) |"
tensor{K<:BraKet}(s1::State, s2::State)

function magnitude(A::Number...)
	if length(A)==2
		return hypot(A[1], A[2])
	end
	return magnitude(hypot(A[1], A[2]), A[3:end]...)
end

magnitude{N<:Number}(A::Array{N}) = magnitude(A...)
magnitude(S::StateRep) = magnitude(S.coeffs)

normalize(S::StateRep) = copy(S, (1/magnitude(S))*S.coeffs)

isnorm(S::StateRep) = magnitude(S)==1


######OperatorRep######################################################

abstract AbstractOperator

immutable Operator <: AbstractOperator
	label
	func::Vector{Function}
end

Operator(label, func::Function...) = Operator(label, Function[func...])
Operator(label) = Operator(label, Function[])

type OperatorRep{B1<:AbstractBasis, B2<:AbstractBasis, N<:Number} <:AbstractOperator
	op::Operator
	coeffs::Matrix{N}
	row_basis::B1
	col_basis::B2
	function OperatorRep(op::Operator, coeffs::Matrix{N}, row_basis::B1, col_basis::B2)
		if (length(row_basis)*length(col_basis))==length(coeffs)
			new(label, coeffs, row_basis, col_basis)
		elseif length(basis)>length(coeffs)
			error("basis is larger than representation by $((length(row_basis)*length(col_basis))-length(coeffs)) labels")
		else
			error("representation is larger than basis by $(length(coeffs)-(length(row_basis)*length(col_basis))) coefficients")
		end
	end	
end

OperatorRep{B1<:AbstractBasis, B2<:AbstractBasis, N<:Number}(op::Operator, coeffs::Matrix{N}, row_basis::B1, col_basis::B2) = OperatorRep{B1,B2,N}(op, coeffs, row_basis, col_basis)
OperatorRep{B<:AbstractBasis, N<:Number}(op::Operator, coeffs::Matrix{N}, basis::B) = OperatorRep(op, coeffs, basis, basis)
OperatorRep{B<:AbstractBasis, N<:Number}(label, coeffs::Matrix{N}, basis::B) = OperatorRep(Operator(label), coeffs, basis)
OperatorRep{B1<:AbstractBasis, B2<:AbstractBasis, N<:Number}(label, coeffs::Matrix{N}, row_basis::B1, col_basis::B2) = OperatorRep(Operator(label), coeffs, row_basis, col_basis)

#####show()############################################
function show(io::IO, b::AbstractBasis)
	println("$(typeof(b)) \"$(b.name)\"")
	println("$(length(b.label_arr[:,1])) Basis StateReps:")
	for i=1:length(b.label_arr[:,1])
		println("| $(repr([b.label_arr[i, :]...])[2:end-1]) $(b.ket_sym)")
	end
end

function show(io::IO, s::StateRep)
	if s.kind == Ket
		left = '|'
		right = s.basis.ket_sym
	else
		left = s.basis.bra_sym
		right = '|'
	end 
	println("$(typeof(s)) $(getlabel(s)):")
	if length(s)!=0
		filled = find(s.coeffs)
		table = cell(length(filled), 2)	
		if length(filled)>=52
			for i=1:25
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= replace("$left $(repr([s.basis.label_arr[filled[i], :]...])[2:end-1]) $right",r"\\", "")
			end
			for i=(length(filled)-25):length(filled)
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= replace("$left $(repr([s.basis.label_arr[filled[i], :]...])[2:end-1]) $right",r"\\", "")
			end
		else
			for i=1:length(filled)
				table[i,1]= s.coeffs[filled[i]]
				table[i,2]= replace("$left $(repr([s.basis.label_arr[filled[i], :]...])[2:end-1]) $right",r"\\", "")
			end
		end
		io = IOBuffer()
		show(io, table)
		io_str = takebuf_string(io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		r_ket = "\"(?=$(s.basis.bra_sym))|(?<=$(s.basis.ket_sym))\""
		r = Regex(r_ket)
		io_str = replace(io_str,r"\"(?=\|)|\|\K\"",  "")
		io_str = replace(io_str,r,  "")
		io_str = replace(io_str,r"\\(?=\")",  "")
		io_str = replace(io_str,r"\s\"\s|\s\"(?=,)",  "\"")
		io_str = replace(io_str,r"\s\K\"(?=\s)", "")
		print(io_str)
	else
		println("(all coefficients are zero)")
	end
end

function show(io::IO, op::AbstractOperator)
	println("$(typeof(op)) $(op.label):")
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


isequal(a::AbstractBasis, b::AbstractBasis) = a.label_arr==b.label_arr
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

function filter(f::Function, B::Basis; name=B.name)
	sub_arr = filter(f, B.label_arr)
	return Basis(name, sub_arr)
end

function filter(f::Function, B::TensorBasis; name=B.name)
	sub_arr = extract(f, B.label_arr)
	return TensorBasis(name, sub_arr, B.basis_arr, bra_sym=B.bra_sym, ket_sym=B.ket_sym)
end

#The following redundancies are implemented for the sake 
#of allowing a[(key)] notation for StateRep coefficient retrieval.
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

isequal(a::StateRep, b::StateRep) = a.basis==b.basis && a.label==b.label && a.coeffs==b.coeffs


copy(S::StateRep) = StateRep(copy(S.label), copy(S.coeffs), S.basis)
copy{N<:Number}(S::StateRep, coeffs::Array{N}) = StateRep(copy(S.label), coeffs, S.basis)
length(S::AbstractState) = length(find(S.coeffs))
size(S::AbstractState) = length(S)

get{B<:Basis}(S::StateRep{B}, key) = S.coeffs[get(S.basis, key)]
get{B<:TensorBasis}(S::StateRep{B}, key) = S.coeffs[get(S.basis, key)]
get{B<:TensorBasis}(S::StateRep{B}, key...) = S.coeffs[get(S.basis, key)]
getindex(S::StateRep, x) = S.coeffs[x]
getindex{B<:Basis}(S::StateRep{B}, x::Tuple) = get(S, x[1])
getindex{B<:TensorBasis}(S::StateRep{B}, x::Tuple) = get(S, x)
setindex!(S::StateRep, y, x) = setindex!(S.coeffs, y, x)

function setindex!(S::StateRep, y::Int, x::Range1{Int})
	for i in x
		setindex!(S, y, i)
	end
end

function setindex!(S::StateRep, y::Array, x::Range1{Int})
	if length(y) != length(x)
		throw(BoundsError())
	else
		for i=1:length(x)
			setindex!(S, y[i], x[i])
		end
	end
end

!(K::Type{Ket}) = Bra
!(B::Type{Bra}) = Ket

function transpose(S::StateRep)
	return copy(S, S.coeffs.')
end

function ctranspose(S::StateRep)
	return copy(S, S.coeffs')
end

function map!(f::Function, s::StateRep)
	s.coeffs = map!(f, s.coeffs)
	return s
end 

map(f::Function, s::StateRep) = map!(f, copy(s))

function maplabel!(f_coeffs::Function, f_labels::Function, s::StateRep)
	labels = filter(f_labels, collect(keys(s.basis.label_map)))	
	for i in labels
		s[get(s.basis, i)] = apply(f_coeffs, get(s, i))
	end
	return s
end

maplabel(f_coeffs::Function, f_labels::Function, s::StateRep) = maplabel!(f_coeffs, f_labels, copy(s))

filter(f::Function, s::StateRep) = maplabel((x)->0, f, s)
filter!(f::Function, s::StateRep) = maplabel!((x)->0, f, s)

*(n::Number, s::StateRep) = copy(s, n*s.coeffs) 
*(s::StateRep, n::Number) = copy(s, s.coeffs*n) 
.+(s::StateRep, n::Number) = copy(s, s.coeffs.+n)
.+(n::Number, s::StateRep) = copy(s, n.+s.coeffs)
.-(s::StateRep, n::Number) = copy(s, s.coeffs.-n)
.-(n::Number, s::StateRep) = copy(s, n.-s.coeffs)
/(s::StateRep, n::Number) = copy(s, s.coeffs/n)
/(n::Number, s::StateRep) = error("cannot divide number by vector")
.^(s::StateRep, n::Number) = copy(s, s.coeffs.^n)
.^(n::Number, s::StateRep) = copy(s, n.^s.coeffs)
^(s::StateRep, n::Int, s_orig=s) = n==1 ? s : (n==2 ? s*s_orig : ^(s*s_orig, n-1, s_orig))


function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::StateRep{B1, Bra}, b::StateRep{B2, Ket})
	return (a.coeffs*b.coeffs)[1]
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::StateRep{B1, Ket}, b::StateRep{B2, Ket})
	new_coeffs = tensor(a.coeffs, b.coeffs)
	return StateRep("$(a.label),$(b.label)", new_coeffs[:,1].*new_coeffs[:,2], a.basis*b.basis)
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::StateRep{B1, Bra}, b::StateRep{B2, Bra})
	new_coeffs = tensor(a.coeffs', b.coeffs')'
	return StateRep("$(a.label),$(b.label)", new_coeffs[1,:].*new_coeffs[2,:], a.basis*b.basis)
end

function *{B1<:AbstractBasis, B2<:AbstractBasis}(a::StateRep{B1, Ket}, b::StateRep{B2, Bra})
	return OperatorRep("$(getlabel(a))$(getlabel(b))", a.coeffs*b.coeffs, a.basis, b.basis)
end

function +(a::StateRep, b::StateRep)
	if a.basis==b.basis && a.kind==b.kind
		StateRep("$(a.label) + $(b.label)", a.coeffs+b.coeffs, a.basis)
	else
		error("composite StateReps not yet supported :(")
	end
end

function -(a::StateRep, b::StateRep)
	if a.basis==b.basis && a.kind==b.kind
		StateRep("$(a.label) - $(b.label)", a.coeffs-b.coeffs, a.basis)
	else
		error("composite StateReps not yet supported :(")
	end
end

copy(op::OperatorRep) = OperatorRep(copy(op.label), copy(op.coeffs), op.row_basis, op.col_basis)
copy{N<:Number}(op::OperatorRep, coeffs::Matrix{N}) = OperatorRep(copy(op.label), coeffs, op.row_basis, op.col_basis)
expm(op::OperatorRep) = copy(op, exp(op.coeffs))
*(op::OperatorRep, n::Number) = copy(op, op.coeffs*n)
*(n::Number, op::OperatorRep) = copy(op, n*op.coeffs)
/(op::OperatorRep, n::Number) = copy(op, op.coeffs/n)
.+(op::OperatorRep, n::Number) = copy(op, op.coeffs.+n)
.+(n::Number, op::OperatorRep) = copy(op, n.+op.coeffs)
.-(op::OperatorRep, n::Number) = copy(op, op.coeffs.-n)
.-(n::Number, op::OperatorRep) = copy(op, n.-op.coeffs)
.^(op::OperatorRep, n::Number) = copy(op, op.coeffs.^n)
.^(n::Number, op::OperatorRep) = copy(op, n.^op.coeffs)


end#module
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
ds = m.StateRep("dS", [1:27], d);
a = m.Basis("a", [1:10], bra_sym="{", ket_sym="}")
as = m.StateRep("aS", [1:10], a)
println("")

