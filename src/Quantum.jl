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
			endof,
			isequal,
			==,
			copy,
			getindex,
			ctranspose,
			show,
			*,
			zero,
			conj,
			eltype,
			hash,
		    showcompact,
		    eltype,
		    get,
		    in,
		    setdiff,
		    +,
		    map,
		    filter,
		    find

	#####################################
	#abstract types######################
	#####################################
	abstract Dirac
	abstract AbstractScalar <: Dirac
	typealias DiracCoeff Union(Number, AbstractScalar)

	#####################################
	#includes############################
	#####################################
	include("crossjoin.jl")
	include("state.jl")
	include("products.jl")
	include("basis.jl")
	# include("scalar.jl")
	# include("diracvector.jl")
	# include("diracmatrix.jl")
	include("misc.jl")

	#####################################
	#exports#############################
	#####################################
	export Dirac,
		   State,
		   Bra,
		   Ket,
		   Tensor,
		   InnerProduct,
		   OuterProduct,
		   AbstractBasis,
		   Basis,
		   TensorBasis,
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
		   bjoin
end