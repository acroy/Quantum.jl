module DiracTools
	#####################################
	#constants###########################
	#####################################
	const lang = "\u27E8"
	const rang = "\u27E9"
	const otimes = "\u2297"

	#####################################
	#imports#############################
	#####################################
	import Base:
		   show,
		   repr,
		   norm,
		   convert,	
		   getindex,
		   setindex!,
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
	abstract Dirac

	abstract BraKet <: Dirac
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Dirac

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

	#####################################
	#includes############################
	#####################################
	include("rep.jl")
	include("state.jl")
	include("basis.jl")
	include("dirac-operator.jl")
	include("dirac-expr.jl")

	#####################################
	#exports#############################
	#####################################
	export Dirac,
		   BraKet,
		   Bra,
		   Ket,
		   AbstractBasis,
		   State,
		   StateRep,
		   Basis,
		   TensorBasis,
		   Operator,
		   kind,
		   statevec,
		   tensor,
		   statejoin,
		   separate,
		   state,
		   normalize!,
		   normalize,
		   mapmatch!,
		   mapmatch,
		   samebasis,
		   findstates,
		   commutator,
		   ptrace
end
