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
			length,
			size,
			endof,
			isequal,
			==,
			copy,
			getindex,
			ctranspose,
			show,
			zero,
			one,
			conj,
			eltype,
			hash,
		    showcompact,
		    eltype,
		    get,
		    in,
		    setdiff,
		    +,
		    -,
		    *,
		    /,
		    ^,
		    (.*),
		    (.+),
		    (.-),
		    (./),
		    (.^),
		    abs,
		    exp,
		    map,
		    filter,
		    find,
		    ndims,
		    findn,
		    findnz,
		    nnz,
		    countnz,
		    ndims,
		    norm,
		    summary,
		    kron,
		    trace,
		    promote_rule,
		    convert

	#####################################
	#abstract types######################
	#####################################
	abstract type Dirac end
	abstract type AbstractScalar<:Dirac end
	abstract type AbstractOperator<:Dirac end
	abstract type State{S}<:Dirac end

	const = DiracCoeff Union(Number, AbstractScalar)

	#####################################
	#includes############################
	#####################################
	include("crossjoin.jl")
	include("state.jl")
	include("products.jl")
	include("basis.jl")
	include("scalar.jl")
	include("diracvector.jl")
	include("diracmatrix.jl")
	include("misc.jl")
	include("fock.jl")

	#####################################
	#exports#############################
	#####################################
	export Dirac,
		   State,
		   Bra,
		   Ket,
		   Tensor,
		   AbstractBasis,
		   AbstractScalar,
		   InnerProduct,
		   ScalarExpr,
		   Basis,
		   DiracVector,
		   DiracMatrix,
		   tensor,
		   basis,
		   dual,
		   label,
		   bsym,
		   separate,
		   isdual,
		   labeldelta,
		   inner,
		   samebasis,
		   kind,
		   bcat,
		   findstates,
		   filtercoeffs,
		   filterstates,
		   getpos,
		   normalize,
		   qeval,
		   dvec,
		   svec,
		   dmat,
		   bra,
		   ket,
		   fbasis,
		   fvec,
		   fcreate,
		   fdestroy,
		   fnum,
		   feye,
		   commutator,
		   ptrace,
		   scalar,
		   actop,
		   filternz,
		   showsp
end
