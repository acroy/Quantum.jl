immutable InnerProduct <: AbstractScalar
	bra::AbstractState{Bra}
	ket::AbstractState{Ket}
end

conj(i::InnerProduct) = InnerProduct(i.ket', i.bra')
bra(i::InnerProduct) = i.bra
ket(i::InnerProduct) = i.ket
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");
##########################################################################
immutable OuterProduct{N<:Number} <: Dirac
	ket::AbstractState{Ket}
	bra::AbstractState{Bra}
end

bra(o::OuterProduct) = o.bra
ket(o::OuterProduct) = o.ket

*(o::OuterProduct, s::AbstractState{Ket}) = DiracVector([(o.bra*s)], statetobasis(o.ket))
*(o::OuterProduct, s::AbstractState{Bra}) = OuterProduct(o.ket, o.bra*s)
*(s::AbstractState{Bra}, o::OuterProduct) = DiracVector([(s*o.ket)], statetobasis(o.bra))
*(s::AbstractState{Ket}, o::OuterProduct) = OuterProduct(s*o.ket, o.bra)


ctranspose(o::OuterProduct) = OuterProduct(o.bra', o.ket')

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");

##########################################################################
immutable Scalar{N<:Number} <: AbstractScalar
	coeff::N
	ivec::Vector{InnerProduct}
end

*(n::Number, s::Scalar) = Scalar(n*s.coeff, s.ivec)
*(s::Scalar, n::Number) = Scalar(s.coeff*n, s.ivec)
*(n::Number, i::InnerProduct) = Scalar(n, [i])
*(i::InnerProduct, n::Number) = *(n,i)
*(a::InnerProduct, b::InnerProduct) = Scalar(1, [a,b])
*(s::Scalar, i::InnerProduct) = Scalar(s.coeff, vcat(s.ivec, i))
*(i::InnerProduct, s::Scalar) = Scalar(s.coeff, vcat(i, s.ivec))
*(a::Scalar, b::Scalar) = Scalar(a.coeff*b.coeff, vcat(a.ivec, b.ivec))
conj(s::Scalar) = Scalar(conj(s.coeff), map(conj, s.ivec))

function show(io::IO, s::Scalar)
	print(io, s.coeff)
	for i=1:length(s.ivec)
		print(io, " * ")
		print(io, s.ivec[i])
	end
end
##########################################################################

typealias DiracCoeff Union(Number, AbstractScalar)

*{C<:DiracCoeff}(c::C, s::AbstractState) = DiracVector([c], statetobasis(s))
*{C<:DiracCoeff}(s::AbstractState, c::C) = *(c,s)
/(s::AbstractScalar, n::Number) = (1/n)*s
-(s::AbstractScalar) = -1*s


##########################################################################
immutable ScalarSum
	n::Number
	terms::Vector{Scalar}
end

+(s::AbstractScalar, n::Number) = ScalarSum(n, [s])
+(n::Number, s::AbstractScalar) = +(s, n)
-(s::AbstractScalar, n::Number) = ScalarSum(-n, [s])
-(n::Number, s::AbstractScalar) = ScalarSum(n, [-s])


function show(io::IO, s::ScalarSum)
	print(io, s.n)
	for i=1:length(s.terms)
		print(io, " + ")
		print(io, s.terms[i])
	end
end