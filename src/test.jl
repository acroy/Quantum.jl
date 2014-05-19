module q
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
		   slice,
		   (.+),
		   (.^),
		   (.-),
		   ^,
		   *,
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
	abstract Quantum
	abstract Dirac <: Quantum
	abstract BraKet <: Quantum
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Quantum
	abstract AbstractState{K<:BraKet} <: Dirac

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

	#####################################
	#includes############################
	#####################################
	include("misc.jl")
	include("statestuff.jl")
	include("basisstuff.jl")
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
		   State,
		   TensorState,
		   Basis,
		   InnerProduct,
		   OuterProduct,
		   DiracVector
		   basislabel,
		   kind,
		   label,
		   statearr,
		   tensor,
		   tensorarr,
		   statejoin,
		   inner
		   separate
		   statetobasis
end
using q

b= Basis([1:10])
print("")
