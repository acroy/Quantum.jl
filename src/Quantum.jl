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
		   convert,	
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
		   hash

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

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

	#####################################
	#includes############################
	#####################################
	include("misc.jl")
	include("state.jl")
	include("products.jl")
	include("basis.jl")
	include("diracvector.jl")
	include("diracmatrix.jl")


	#####################################
	#additional functions################
	#####################################

	samebasis(a::AbstractBasis, b::AbstractBasis)= a==b
	samebasis(a::DiracMatrix, b::DiracMatrix) = a.rowbasis==b.rowbasis && a.colbasis==b.colbasis
	samebasis(a::DiracVector, b::DiracVector) = a.basis==b.basis
	samebasis(a::AbstractState, b::AbstractState) = basislabel(a)==basislabel(b)
	samebasis(a::InnerProduct, b::InnerProduct) = basislabel(a)==basislabel(b)
	samebasis(a::OuterProduct, b::OuterProduct) = basislabel(a)==basislabel(b)
	samebasis{T,K}(a::AbstractState{K}, b::DiracVector{T,K}) = basislabel(a)==label(b.basis)
	samebasis{T,K}(a::DiracVector{T,K}, b::AbstractState{K}) = samebasis(b,a)
	samebasis{K}(a::AbstractState{K}, b::Basis{K}) = basislabel(a)==label(b)
	samebasis{K}(a::Basis{K}, b::AbstractState{K}) = samebasis(b,a)
	
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
		   DiracCoeff,
		   ScalarExpr,
		   qeval,
		   tobasis,
		   basislabel,
		   kind,
		   label,
		   statearr,
		   tensor,
		   tensorarr,
		   statejoin,
		   inner,
		   separate,
		   normalize,
		   mapmatch,
		   mapmatch!,
		   filtercoeffs,
		   filterstates,
		   getpos,
		   ptrace,
		   isdual,
		   samebasis
end
print("")
