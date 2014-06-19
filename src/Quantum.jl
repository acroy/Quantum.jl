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
			eltype
		   # show,
		   # summary,
		   # showcompact,
		   # norm,
		   # getindex,
		   # setindex!,
		   # conj,
		   # copy,
		   # ndims,
		   # size,
		   # length,
		   # norm,
		   # (.+),
		   # (.^),
		   # (.-),
		   # (.*),
		   # (./),
		   # ^,
		   # *,
		   # +,
		   # -,
		   # /,
		   # abs,
		   # in,
		   # setdiff,
		   # get,
		   # exp,
		   # expm,
		   # map,
		   # map!,
		   # filter,
		   # isequal,
		   # ==,
		   # endof,
		   # find,
		   # findn,
		   # findnz,
		   # trace,
		   # hash,
		   # kron,
		   # eltype,
		   # zero

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
		   Single,
		   Tensor,
		   AbstractBra,
		   AbstractKet,
		   Bra,
		   Ket,
		   TensorBra,
		   TensorKet,
		   InnerProduct,
		   OuterProduct,
		   KetBasis,
		   BraBasis,
		   TensorKetBasis,
		   TensorBraBasis,
		   tensor,
		   basis,
		   dual,
		   label,
		   bsym,
		   separate,
		   isdual,
		   labeldelta,
		   inner,
		   statearr,
		   samebasis,
		   statejoin
		   # State,
		   # TensorState,
		   # Basis,
		   # TensorBasis,
		   # InnerProduct,
		   # OuterProduct,
		   # DiracVector,
		   # DiracMatrix,
		   # ScalarExpr,
		   # qeval,
		   # tobasis,
		   # btensor,
		   # basislabel,
		   # basisjoin,
		   # kind,
		   # label,
		   # statearr,
		   # tensor,
		   # inner,
		   # separate,
		   # normalize,
		   # mapmatch,
		   # mapmatch!,
		   # filtercoeffs,
		   # filterstates,
		   # findstates,
		   # getpos,
		   # ptrace,
		   # isdual,
		   # labeldelta,
		   # samebasis,
		   # statecross,
		   # crossjoin,
		   # statemapper,
		   # statejoin,
		   # statecat
end