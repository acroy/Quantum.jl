Quantum.jl API
===
I. Type Implementations
===
1. AbstractTypes
---
The following list contains the abstract types referenced in this documentation:

	abstract Quantum
	abstract BraKet <: Quantum
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Quantum

The not operator `!` can be applied to `Bra` and `Ket` types to alternate between 
the two (since they are duals of each other):

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

2. State 
--- 
__Description__

A State is a type of an object that has only a
label (stored as a Vector) and a specification of whether it belongs to Ket-
space or Bra-space (this property is referred to as “kind”). A state is
parameterized by its kind; thus, states are either of type State{Bra} or
State{Ket}. State is a subtype of Quantum.

__Definition__

	immutable State{K<:BraKet} <: Quantum
	  label::Vector
	  kind::Type{K}
	end

__Constructors__

	State(label::Vector) #defaults to kind=Ket
	State{K<:BraKet}(label, kind::Type{K}=Ket)
	State{K<:BraKet}(label...; kind::Type{K}=Ket)

__Methods and Examples__

For the sake of this example, we'll instantiate a few `State` s thusly:

	julia> using QuantumJL

	julia> a = State(:a)
	| :a ⟩

	julia> b = State("b")
	| "b" ⟩

	julia> c = State(3)
	| 3 ⟩

	julia> abc = a*b*c
	| :a,"b",3 ⟩

You can index into a state's label just as you would index into an 
array, performing both retrieval and assignment of the label elements 
idomatically:

	julia> abc[2]
	"b"

	julia> abc[2]=:n
	:n

	julia> abc
	| :a,:n,3 ⟩

You can then separate the composite basis state bc into its component states:

	julia> separate(abc)
	3-element Array{State{Ket},1}:
	 | :a ⟩
	 | :n ⟩
	 | 3 ⟩

Several convenience functions are provided which aid in the 
contruction of states. One such function is `statevec`, which
takes in an array and produces a Vector{State{Ket}} based on
the array's contents:

	julia> arr = [1:5 6:10 11:15]
	5x3 Array{Int64,2}:
	 1   6  11
	 2   7  12
	 3   8  13
	 4   9  14
	 5  10  15

	julia> statevec(arr)
	5-element Array{State{Ket},1}:
	 | 1,6,11 ⟩
	 | 2,7,12 ⟩
	 | 3,8,13 ⟩
	 | 4,9,14 ⟩
	 | 5,10,15 ⟩

Finally, one can take the complex conjugate transpose and inner product of abstract states:

	julia> a'
	⟨ :a |

	julia> a'*a
	1

	julia> a'*b
	0

Currently, outer products are not supported in a normal sense. In the future, however, Quantum.jl will implement some form of Expr manipulation in order to
accommodate more abstract operations. 

3. Basis 
--- 
__Description__ 

A `Basis` is a mapping of labels onto the indices of a state vector
representing a given Hilbert subspace. In other words, `Basis` provides a
label for each position in a coefficient vector that represents a `State` in
that `Basis`. It is a finite dimensional subspace of the Hilbert space that
has been selected according to some arbitrary labeling or behavioral pattern,
as defined by the user.  A graph theoretic approach to thinking about 
bases would be to consider the `Basis` as a map necessary to traverse a state
representation graph, whose nodes are the states forming a given Hilbert
subspace and whose edges are weighted by coefficient values.

Taking the tensor product of `Basis` objects results in a new object of type
`TensorBasis`, which points to the original “separable” `Basis` objects in
addition to storing the new multi-system basis that resulted from taking the
tensor product. This allows for the optimization of operations on
matrices/vectors represented in a `TensorBasis` (for example, taking the partial
trace of an `OperatorRep`).

__Definition__

	immutable Basis{K<:BraKet} <: AbstractBasis{K}
		label
		states::Vector{State{K}}
		label_map::Dict{Vector, Int}
	end

	immutable TensorBasis{K<:BraKet} <: AbstractBasis{K}
		bases::Vector{Basis{K}}
		states::Vector{State{K}}
		label_map::Dict{Vector, Int}
	end

__Constructors__

	Basis{K<:BraKet}(label, states::Vector{State{K}})
	Basis{K<:BraKet}(label, states::State{K}...)
	Basis(label, label_vec::Vector)	

	TensorBasis{K<:BraKet}(bases::Vector{Basis{K}}, states::Vector{State{K}})

__Methods and Examples__

You can form a basis from a `Vector` of states, or a `Vector` of state labels:

	julia> Basis("abc", separate(abc))
	Basis{Ket} abc:
	| :a ⟩
	| :n ⟩
	| 3 ⟩

	julia> bs = Basis("bs",[1:5])
	Basis{Ket} bs:
	| 1 ⟩
	| 2 ⟩
	| 3 ⟩
	| 4 ⟩
	| 5 ⟩

It is also possible to take the tensor product of bases:

	julia> bs = bs*bs
	TensorBasis{Ket} bs ⊗ bs:
	| 1,1 ⟩
	| 1,2 ⟩
	| 1,3 ⟩
	| 1,4 ⟩
	| 1,5 ⟩
	| 2,1 ⟩
	| 2,2 ⟩
	| 2,3 ⟩
	| 2,4 ⟩
	| 2,5 ⟩
	| 3,1 ⟩
	| 3,2 ⟩
	| 3,3 ⟩
	| 3,4 ⟩
	| 3,5 ⟩
	| 4,1 ⟩
	| 4,2 ⟩
	| 4,3 ⟩
	| 4,4 ⟩
	| 4,5 ⟩
	| 5,1 ⟩
	| 5,2 ⟩
	| 5,3 ⟩
	| 5,4 ⟩
	| 5,5 ⟩


You can use `filter` to extract subspaces of your basis:

	julia> bs = filter(x->x[1]%2==0, bs)
	TensorBasis{Ket} sub_(bs ⊗ bs)_1 ⊗ sub_(bs ⊗ bs)_2:
	| 2,1 ⟩
	| 2,2 ⟩
	| 2,3 ⟩
	| 2,4 ⟩
	| 2,5 ⟩
	| 4,1 ⟩
	| 4,2 ⟩
	| 4,3 ⟩
	| 4,4 ⟩
	| 4,5 ⟩

You can separate a `TensorBasis` into component its `Basis` objects (note that
these `Basis` objects are representative of `bs` *after* the filter operation, not before it):

	julia> separate(bs)
	2-element Array{Basis{Ket},1}:
	 Basis{Ket} sub_(bs ⊗ bs)_1:
	 | 2 ⟩
	 | 4 ⟩
	 Basis{Ket} sub_(bs ⊗ bs)_2:
	 | 1 ⟩
	 | 2 ⟩
	 | 3 ⟩
	 | 4 ⟩
	 | 5 ⟩

For the representation of `State`s in the dual space, you must have
a basis in the dual space. Thus, the `ctranspose` function  

	julia> bs'
	TensorBasis{Bra} sub_(bs ⊗ bs)_1 ⊗ sub_(bs ⊗ bs)_2:
	⟨ 2,1 |
	⟨ 2,2 |
	⟨ 2,3 |
	⟨ 2,4 |
	⟨ 2,5 |
	⟨ 4,1 |
	⟨ 4,2 |
	⟨ 4,3 |
	⟨ 4,4 |
	⟨ 4,5 |

4. StateRep
---
__Description__ 

You can represent a state in a basis by providing a vector of coefficients
that corresponds to the inner product of the state with it's basis states:

	julia> sr = StateRep(a, [1:length(bs)], bs)
	StateRep{Ket} | :a ; bs ⊗ bs ⟩:
	  1.0+0.0im  | 1,1 ⟩
	  2.0+0.0im  | 1,2 ⟩
	  3.0+0.0im  | 1,3 ⟩
	  4.0+0.0im  | 1,4 ⟩
	  5.0+0.0im  | 1,5 ⟩
	  6.0+0.0im  | 2,1 ⟩
	  7.0+0.0im  | 2,2 ⟩
	  8.0+0.0im  | 2,3 ⟩
	  9.0+0.0im  | 2,4 ⟩
	 10.0+0.0im  | 2,5 ⟩
	          ⋮
	 16.0+0.0im  | 4,1 ⟩
	 17.0+0.0im  | 4,2 ⟩
	 18.0+0.0im  | 4,3 ⟩
	 19.0+0.0im  | 4,4 ⟩
	 20.0+0.0im  | 4,5 ⟩
	 21.0+0.0im  | 5,1 ⟩
	 22.0+0.0im  | 5,2 ⟩
	 23.0+0.0im  | 5,3 ⟩
	 24.0+0.0im  | 5,4 ⟩
	 25.0+0.0im  | 5,5 ⟩

	julia> normalize!(sr)
	StateRep{Ket} | :a ;åå bs ⊗ bs ⟩:
	 0.0134535+0.0im  | 1,1 ⟩
	 0.0269069+0.0im  | 1,2 ⟩
	 0.0403604+0.0im  | 1,3 ⟩
	 0.0538138+0.0im  | 1,4 ⟩
	 0.0672673+0.0im  | 1,5 ⟩
	 0.0807207+0.0im  | 2,1 ⟩
	 0.0941742+0.0im  | 2,2 ⟩
	  0.107628+0.0im  | 2,3 ⟩
	  0.121081+0.0im  | 2,4 ⟩
	  0.134535+0.0im  | 2,5 ⟩
	               ⋮
	  0.215255+0.0im  | 4,1 ⟩
	  0.228709+0.0im  | 4,2 ⟩
	  0.242162+0.0im  | 4,3 ⟩
	  0.255616+0.0im  | 4,4 ⟩
	  0.269069+0.0im  | 4,5 ⟩
	  0.282523+0.0im  | 5,1 ⟩
	  0.295976+0.0im  | 5,2 ⟩
	  0.309429+0.0im  | 5,3 ⟩
	  0.322883+0.0im  | 5,4 ⟩
	  0.336336+0.0im  | 5,5 ⟩

