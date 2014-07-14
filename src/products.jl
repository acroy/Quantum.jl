#####################################
#InnerProduct/OuterProduct###########
#####################################

immutable OuterProduct{K<:Ket,B<:Bra} <: AbstractOperator
	ket::State{K}
	bra::State{B}
end

immutable InnerProduct <: AbstractScalar
	bra::State
	ket::State
	InnerProduct{B<:Bra,K<:Ket}(b::State{B}, k::State{K}) = new(b,k)
end

promote_rule{N<:Number, I<:InnerProduct}(::Type{I}, ::Type{N}) = ScalarExpr
promote_rule{N<:Number, I<:InnerProduct, T}(::Type{Array{I,T}}, ::Type{Array{N,T}}) = Array{ScalarExpr,T}
promote_rule{N<:Number, I<:InnerProduct, T}(::Type{SparseMatrixCSC{I,T}}, ::Type{SparseMatrixCSC{N,T}}) = SparseMatrixCSC{ScalarExpr,T}

#####################################
#Misc Functions######################
#####################################

=={O<:OuterProduct}(a::O, b::O) = a.ket==b.ket && a.bra==b.bra
isequal{O<:OuterProduct}(a::O, b::O) = isequal(a.ket,b.ket) && isequal(a.bra,b.bra)

=={I<:InnerProduct}(a::I, b::I) = a.ket==b.ket && a.bra==b.bra
isequal{I<:InnerProduct}(a::I, b::I) = isequal(a.ket,b.ket) && isequal(a.bra,b.bra)

bsym(o::OuterProduct) = [bsym(o.ket), bsym(o.bra)]
label(o::OuterProduct) = [label(o.ket), label(o.bra)]
bsym(i::InnerProduct) = [bsym(i.bra), bsym(i.ket)]
label(i::InnerProduct) = [label(i.bra), label(o.ket)]
conj(i::InnerProduct) = InnerProduct(i.ket', i.bra')
ctranspose(o::OuterProduct) = OuterProduct(o.bra', o.ket')
ctranspose(i::InnerProduct) = conj(i)

isdual(a::OuterProduct, b::OuterProduct) = isdual(a.ket, b.bra) && isdual(b.ket, a.bra)
isdual(a::InnerProduct, b::InnerProduct) = isdual(a.ket, b.bra) && isdual(b.ket, a.bra)

#####################################
#Show Functions######################
#####################################

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra))$(repr(i.ket)[2:end])");

#####################################
#Arithmetic Operations###############
#####################################

*{K<:Ket}(o::OuterProduct, s::State{K}) = dvec([inner(o.bra,s)], basis(o.ket))
*{B<:Bra}(s::State{B}, o::OuterProduct) = dvec([inner(s,o.ket)], basis(o.bra))

kron{B<:Bra}(o::OuterProduct, s::State{B}) = OuterProduct(o.ket, tensor(o.bra,s))
kron{B<:Bra}(s::State{B}, o::OuterProduct) = OuterProduct(o.ket, tensor(s,o.bra))
kron{K<:Ket}(o::OuterProduct, s::State{K}) = OuterProduct(tensor(o.ket,s), o.bra)
kron{K<:Ket}(s::State{K}, o::OuterProduct) = OuterProduct(tensor(s,o.ket), o.bra)
kron(c::DiracCoeff, o::OuterProduct) = dmat([c], basis(o.ket), basis(o.bra))
kron(o::OuterProduct, c::DiracCoeff) = kron(c,o)

kron(a::OuterProduct, b::OuterProduct) = OuterProduct(tensor(a.ket, b.ket), tensor(a.bra, b.bra))

*(a::OuterProduct, b::OuterProduct) = dmat([a.bra*b.ket], basis(a.ket), basis(b.bra))

*(c::DiracCoeff, o::OuterProduct) = kron(c,o)
*(o::OuterProduct, c::DiracCoeff) = kron(c,o)
-(o::OuterProduct) = kron(-1,o)

function +{K<:Ket, B<:Bra}(a::OuterProduct{K,B}, b::OuterProduct{K,B})
	if a==b
		return dmat([2.0], basis(a.ket), basis(b.bra))
	else
		rowb = basis(a.ket, b.ket)
		colb = basis(a.bra, b.bra)
		res = dmat(zeros(2,2)[1:length(rowb), 1:length(colb)], rowb, colb)
		res[getpos(res, a)...] = 1.0
		res[getpos(res, b)...] = 1.0
		return res
	end
end

-(a::OuterProduct, b::OuterProduct) = a+(-b)
