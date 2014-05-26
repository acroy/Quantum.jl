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
		   map,
		   map!,
		   filter,
		   isequal,
		   copy,
		   hash,
		   isequal,
		   endof,
		   start,
		   find,
		   trace

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
		   InnerProduct,
		   OuterProduct,
		   DiracVector,
		   DiracCoeff,
		   ScalarExpr,
		   qeval,
		   statetobasis,
		   basislabel,
		   kind,
		   label,
		   statearr,
		   tensor,
		   tensorarr,
		   statejoin,
		   inner,
		   separate,
		   normalize!,
		   normalize,
		   mapmatch,
		   mapmatch!,
		   filtercoeffs,
		   filterstates
end
print("")
