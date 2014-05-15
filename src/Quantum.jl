module QuantumJL
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
		   label,
		   labeldelta,
		   statejoin,
		   separate,
		   state,
		   statetobasis,
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
