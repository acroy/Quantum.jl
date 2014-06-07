module Quantum
	#####################################
	#constants###########################
	#####################################
	const lang = "\u27E8"
	const rang = "\u27E9"
	const otimes = "\u2297"
	const vdots ="\u205E"

	#####################################
	#imports#############################
	#####################################
	import Base:
		   show,
		   showcompact,
		   norm,
		   getindex,
		   setindex!,
		   conj,
		   copy,
		   ndims,
		   size,
		   length,
		   norm,
		   (.+),
		   (.^),
		   (.-),
		   (.*),
		   (./),
		   ^,
		   *,
		   +,
		   -,
		   /,
		   abs,
		   in,
		   setdiff,
		   get,
		   !,
		   exp,
		   expm,
		   map,
		   map!,
		   filter,
		   isequal,
		   ==,
		   endof,
		   find,
		   trace,
		   hash,
		   kron

	#####################################
	#abstract types######################
	#####################################
	abstract Dirac
	abstract BraKet <: Dirac
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Dirac
	abstract AbstractState{K<:BraKet} <: Dirac
	abstract AbstractScalar <: Dirac
	typealias DiracCoeff Union(Number, AbstractScalar)

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

	#####################################
	#includes############################
	#####################################
	include("crossjoin.jl")
	include("state.jl")
	include("basis.jl")
	include("products.jl")
	include("scalar.jl")
	include("diracvector.jl")
	include("diracmatrix.jl")

	#####################################
	#additional functions################
	#####################################

	samebasis(a::Dirac, b::Dirac)= basislabel(a)=="?" || basislabel(b)=="?" ? false : basislabel(a)==basislabel(b)
	samebasis(a::String, b::Dirac)= a=="?" ? false : a==basislabel(b)
	samebasis(a::Dirac, b::String)= samebasis(b, a)

	#additive identities
	for t=(:DiracMatrix, :DiracVector, :AbstractState, :OuterProduct)
		@eval +(n::Number, d::($t)) = n==0 ? d : error("cannot add number and $(string($t))")
		@eval +(d::($t), n::Number) = +(n,d)
		@eval -(n::Number, d::($t)) = +(n,d)
		@eval -(d::($t), n::Number) = +(n,d)
	end

	kron(a::DiracMatrix, b::DiracMatrix) = DiracMatrix(kron(a.coeffs, b.coeffs), tensor(a.rowbasis, b.rowbasis), tensor(a.colbasis, b.colbasis)) 
	kron(op::DiracMatrix, d::DiracVector{Ket}) = DiracMatrix(kron(op.coeffs, d.coeffs), tensor(op.rowbasis, d.basis), op.colbasis)
	kron(op::DiracMatrix, d::DiracVector{Bra}) = DiracMatrix(kron(op.coeffs, d.coeffs), op.rowbasis, tensor(op.colbasis, d.basis))
	kron(d::DiracVector{Ket}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), tensor(d.basis, op.rowbasis), op.colbasis)
	kron(d::DiracVector{Bra}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), op.rowbasis, tensor(d.basis, op.colbasis))
	kron{K}(a::DiracVector{K}, b::DiracVector{K}) = DiracVector(kron(a.coeffs, b.coeffs), tensor(a.basis, b.basis))
	kron(a::DiracVector{Ket}, b::DiracVector{Bra}) = DiracMatrix(kron(a.coeffs, b.coeffs), a.basis, b.basis)
	kron(a::DiracVector{Bra}, b::DiracVector{Ket}) = kron(b,a)
	kron(a::AbstractScalar, b::AbstractScalar) = a*b
	kron(a::Dirac, b::DiracCoeff) = a*b
	kron(a::DiracCoeff, b::Dirac) = b*a

	#####################################
	#exports#############################
	#####################################
	export Quantum,
		   BraKet,
		   Bra,
		   Ket,
		   AbstractBasis,
		   AbstractState,
		   AbstractScalar,
		   State,
		   TensorState,
		   Basis,
		   TensorBasis,
		   InnerProduct,
		   OuterProduct,
		   DiracVector,
		   DiracMatrix,
		   ScalarExpr,
		   qeval,
		   tobasis,
		   basislabel,
		   kind,
		   label,
		   statearr,
		   tensor,
		   inner,
		   separate,
		   normalize,
		   mapmatch,
		   mapmatch!,
		   filtercoeffs,
		   filterstates,
		   findstates,
		   getpos,
		   ptrace,
		   isdual,
		   samebasis
end
