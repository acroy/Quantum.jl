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

	samebasis(a::Dirac, b::Dirac)= basislabel(a)=="?" || basislabel(b)=="?" ? false : basislabel(a)==basislabel(b)
	samebasis(a::String, b::Dirac)= a=="?" ? false : a==basislabel(b)
	samebasis(a::Dirac, b::String)= samebasis(b, a)
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
