module QuantumJL
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
		   showcompact,
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
	abstract Quantum

	abstract BraKet <: Quantum
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Quantum

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

	#####################################
	#includes############################
	#####################################
	include("rep.jl")
	include("state.jl")
	include("basis.jl")
	include("operator.jl")
	include("qexpr.jl")

	#####################################
	#exports#############################
	#####################################
	export Quantum,
		   BraKet,
		   Bra,
		   Ket,
		   AbstractBasis,
		   State,
		   StateRep,
		   Basis,
		   TensorBasis,
		   OperatorRep,
		   kind,
		   statevec,
		   tensor,
		   labeldelta,
		   statejoin,
		   separate,
		   state,
		   normalize!,
		   normalize,
		   mapmatch!,
		   mapmatch,
		   filtercoeffs, 
		   filtercoeffs!,
		   filterstates,
		   filterstates!,
		   samebasis,
		   findstates,
		   commutator,
		   ptrace,
		   qeval
end
