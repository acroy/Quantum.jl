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
		   summary,
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
		   findn,
		   findnz,
		   trace,
		   hash,
		   kron,
		   eltype,
		   zero

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

	kind{T<:AbstractBasis{Ket}}(t::Type{T}) = Ket
	kind{T<:AbstractBasis{Bra}}(t::Type{T}) = Bra
	kind{T<:AbstractState{Ket}}(t::Type{T}) = Ket
	kind{T<:AbstractState{Bra}}(t::Type{T}) = Bra
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
	include("misc.jl")

	#####################################
	#exports#############################
	#####################################
	export Dirac,
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
		   btensor,
		   basislabel,
		   basisjoin,
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