Quantum.jl Documentation
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

A `State` is a type of an object that has only a
label (stored as a `Vector`) and a specification of whether it belongs to `Ket`-
space or `Bra`-space (this property is referred to as “kind”). A `State` is
parameterized by its kind; thus, states are either of type `State{Bra}` or
`State{Ket}`. `State` is a subtype of `Quantum`.

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

For the sake of this example, we'll instantiate a few states thusly:

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
idiomatically:

	julia> abc[2]
	"b"

	julia> abc[2]=:n
	:n

	julia> abc
	| :a,:n,3 ⟩

Composite states can be separated into their component states:

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

Currently, outer products of abstract states are not supported in a normal
sense. In the future, however, Quantum.jl will implement some form of Expr
manipulation in order to accommodate more abstract operations.

3. Basis 
--- 
__Description__ 

A `Basis` is a mapping of labels onto the indices of a state vector
representing a given Hilbert subspace. In other words, `Basis` provides a
label for each position in a coefficient vector that represents a `State` in
that `Basis`. A basis is a finite dimensional subspace of the Hilbert space that
has been selected according to some arbitrary labeling or behavioral pattern,
as defined by the user.  A graph theoretic approach to thinking about 
bases would be to consider objects of the type `Basis` as maps used to traverse 
a state representation graph, whose nodes are the states forming a given Hilbert
subspace and whose edges are weighted by coefficient values.

Taking the tensor product of `Basis` objects results in a new object of type
`TensorBasis`, which points to the original “separable” `Basis` objects in
addition to storing the new multi-system basis that resulted from taking the
tensor product. This allows for the optimization of operations on
states and operators represented in a `TensorBasis` (for example, taking the partial
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

	julia> qb = Basis("qb",[1,0])
	Basis{Ket} qb:
	| 1 ⟩
	| 0 ⟩

It is also possible to take the tensor product of bases:

	julia> tb = qb*qb*qb
	TensorBasis{Ket} qb ⊗ qb ⊗ qb:
	| 1,1,1 ⟩
	| 1,1,0 ⟩
	| 1,0,1 ⟩
	| 1,0,0 ⟩
	| 0,1,1 ⟩
	| 0,1,0 ⟩
	| 0,0,1 ⟩
	| 0,0,0 ⟩

Assignment and retrieval is possible in the normal manner:

	julia> tb[1]
	| 1,1,1 ⟩

You can use `filter` to extract subspaces of a basis via arbitrary selection
rules (as defined by the function passed to `filter`):

	julia> tb = filter(x->x[1]==1, tb)
	TensorBasis{Ket} sub_(qb ⊗ qb ⊗ qb)_1 ⊗ sub_(qb ⊗ qb ⊗ qb)_2 ⊗ sub_(qb ⊗ qb ⊗ qb)_3:
	| 1,1,1 ⟩
	| 1,1,0 ⟩
	| 1,0,1 ⟩
	| 1,0,0 ⟩

Note that the component bases were transformed by the filtering operation. We
can see this explicitly using `separate`:

	julia> separate(tb)
	3-element Array{Basis{Ket},1}:
	Basis{Ket} sub_(qb ⊗ qb ⊗ qb)_1:
	| 1 ⟩
	Basis{Ket} sub_(qb ⊗ qb ⊗ qb)_2:
	| 1 ⟩
	| 0 ⟩
	Basis{Ket} sub_(qb ⊗ qb ⊗ qb)_3:
	| 1 ⟩
	| 0 ⟩

For the representation of `State{Bra}`s, you must have
a basis in the dual space. Thus, the `ctranspose` function
acts accordingly:

	julia> tb'
	TensorBasis{Bra} sub_(qb ⊗ qb ⊗ qb)_1 ⊗ sub_(qb ⊗ qb ⊗ qb)_2 ⊗ sub_(qb ⊗ qb ⊗ qb)_3:
	⟨ 1,1,1 |
	⟨ 1,1,0 |
	⟨ 0,0,1 |
	⟨ 0,0,0 |

4. StateRep
---
__Description__ 

A `StateRep` is the representation of a `State` in a `Basis`. It stores
the basis alongside a corresponding vector of complex coefficients.
For example, if we want to represent `| a: ⟩` in the `qb` basis defined above, 
the coefficients are given by the theoretical operations:
	
	⟨ 1 | a: ⟩ = c_1
	⟨ 0 | a: ⟩ = c_0

In practice, the results of the inner products shown above 
are either asserted or derived from similar assertions. 

__Definition__

	type StateRep{K<:BraKet} <: Quantum
		state::State{K}
		coeffs::Array{Complex{Float64}}
		basis::AbstractBasis{K}
	end

__Constructors__

	StateRep{N<:Number, K<:BraKet}(s::State{K}, coeffs::Array{N}, basis::AbstractBasis)
	StateRep{N<:Number}(label, coeffs::Array{N}, basis::AbstractBasis)

__Methods and Examples__

Defining `sr` as the representation of `| a: ⟩` in the `qb` basis:

	julia> sr = StateRep(a, [1,1], qb)
	StateRep{Ket} | :a ; qb ⟩:
	 1.0+0.0im  | 1 ⟩
	 1.0+0.0im  | 0 ⟩

	julia> normalize!(sr)
	StateRep{Ket} | :a ; qb ⟩:
	 0.707107+0.0im  | 1 ⟩
	 0.707107+0.0im  | 0 ⟩

`sr`'s basis transforms correctly along with the representation vector:

	julia> tsr=sr*sr*sr
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	 0.353553+0.0im  | 1,1,1 ⟩
	 0.353553+0.0im  | 1,1,0 ⟩
	 0.353553+0.0im  | 1,0,1 ⟩
	 0.353553+0.0im  | 1,0,0 ⟩
	 0.353553+0.0im  | 0,1,1 ⟩
	 0.353553+0.0im  | 0,1,0 ⟩
	 0.353553+0.0im  | 0,0,1 ⟩
	 0.353553+0.0im  | 0,0,0 ⟩

	julia> tsr'
	StateRep{Bra} ⟨ :a,:a,:a ; qb ⊗ qb ⊗ qb |:
	               ⟨ 1,1,1 |                ⟨ 1,1,0 |  …                ⟨ 0,0,0 |
	 0.353553-0.0im           0.353553-0.0im              0.353553-0.0im

Assignment and retrieval refer only to the coefficients:

	julia> tsr[1]
	0.3535533905932737 + 0.0im

	julia> tsr[1]=.6+3im
	0.6 + 3.0im

	julia> tsr
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	      0.6+3.0im  | 1,1,1 ⟩
	 0.353553+0.0im  | 1,1,0 ⟩
	 0.353553+0.0im  | 1,0,1 ⟩
	 0.353553+0.0im  | 1,0,0 ⟩
	 0.353553+0.0im  | 0,1,1 ⟩
	 0.353553+0.0im  | 0,1,0 ⟩
	 0.353553+0.0im  | 0,0,1 ⟩
	 0.353553+0.0im  | 0,0,0 ⟩

You can perform inner/outer products and vector/matrix arithmetic in a normal
fashion (let `tqb=qb*qb*qb`):

	julia> tsr'*tsr
	0.9999999999999998

	julia> tqb[1]'*tsr
	0.6+3.0im

	julia> State("bob")'*tsr
	0

It is assumed when performing linear algebraic operations with normal vectors/matrices
that they are in the same basis as the state representation:

	julia> tsr+[1:8]
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	     1.6+3.0im  | 1,1,1 ⟩
	 2.35355+0.0im  | 1,1,0 ⟩
	 3.35355+0.0im  | 1,0,1 ⟩
	 4.35355+0.0im  | 1,0,0 ⟩
	 5.35355+0.0im  | 0,1,1 ⟩
	 6.35355+0.0im  | 0,1,0 ⟩
	 7.35355+0.0im  | 0,0,1 ⟩
	 8.35355+0.0im  | 0,0,0 ⟩

	julia> ([1:8]*[1:8]')*tsr
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	  12.9744+3.0im  | 1,1,1 ⟩
	  25.9487+6.0im  | 1,1,0 ⟩
	  38.9231+9.0im  | 1,0,1 ⟩
	 51.8975+12.0im  | 1,0,0 ⟩
	 64.8718+15.0im  | 0,1,1 ⟩
	 77.8462+18.0im  | 0,1,0 ⟩
	 90.8206+21.0im  | 0,0,1 ⟩
	 103.795+24.0im  | 0,0,0 ⟩

	julia> ([1:8]*[1:7]')*tsr
	ERROR: DimensionMismatch("*")

The function `filter` is implemented in two different ways; `filtercoeffs` and
`filterstates`.  As one might expect, the former checks against the
coefficients while the latter checks against the basis states. Elements that
cause the filtering function to return false cause the corresponding
coefficient to be set to zero.

A function called `mapmatch` is provided that only maps the function 
passed to it to the coefficients whose states pass a filter test:

	julia> mapmatch(c->c+100, s->s[1]==s[2], tsr)
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	  100.354+0.0im  | 1,1,1 ⟩
	  100.354+0.0im  | 1,1,0 ⟩
	 0.353553+0.0im  | 1,0,1 ⟩
	 0.353553+0.0im  | 1,0,0 ⟩
	 0.353553+0.0im  | 0,1,1 ⟩
	 0.353553+0.0im  | 0,1,0 ⟩
	  100.354+0.0im  | 0,0,1 ⟩
	  100.354+0.0im  | 0,0,0 ⟩

Finally, taking the outer product of two state representations returns 
an `OperatorRep`:

	julia> normalize!(tsr)
	StateRep{Ket} | :a,:a,:a ; qb ⊗ qb ⊗ qb ⟩:
	 0.187546+0.937729im  | 1,1,1 ⟩
	      0.110512+0.0im  | 1,1,0 ⟩
	      0.110512+0.0im  | 1,0,1 ⟩
	      0.110512+0.0im  | 1,0,0 ⟩
	      0.110512+0.0im  | 0,1,1 ⟩
	      0.110512+0.0im  | 0,1,0 ⟩
	      0.110512+0.0im  | 0,0,1 ⟩
	      0.110512+0.0im  | 0,0,0 ⟩

	julia> tsr*tsr'
	OperatorRep:
	                               ⟨ 1,1,1 |  …                      ⟨ 0,0,0 |
	  | 1,1,1 ⟩        0.914509+0.0im              0.0207261+0.103631im
	  | 1,1,0 ⟩  0.0207261-0.103631im                    0.012213+0.0im
	  | 1,0,1 ⟩  0.0207261-0.103631im                    0.012213+0.0im
	  | 1,0,0 ⟩  0.0207261-0.103631im                    0.012213+0.0im
	  | 0,1,1 ⟩  0.0207261-0.103631im           …        0.012213+0.0im
	  | 0,1,0 ⟩  0.0207261-0.103631im                    0.012213+0.0im
	  | 0,0,1 ⟩  0.0207261-0.103631im                    0.012213+0.0im
	  | 0,0,0 ⟩  0.0207261-0.103631im                    0.012213+0.0im 

5. OperatorRep 
--- 
__Description__
An `OperatorRep` is the representation of an abstract operator in a basis. It is 
a matrix of complex coefficients. 

__Definition__

	type OperatorRep <: Quantum
		coeffs::Matrix{Complex{Float64}}
		row_basis::AbstractBasis{Ket}
		col_basis::AbstractBasis{Bra}
	end

Note that it is possible for an object of type `OperatorRep` to have separate bases corresponding
to the rows and columns of a matrix. While such a construction would be quite strange, 
it would not be invalid under the rules of linear algebra, and implementing it does not
cost much.

__Constructors__

	OperatorRep{N<:Number}(coeffs::Matrix{N}, b::AbstractBasis)
	OperatorRep{N<:Number}(coeffs::Matrix{N}, row_basis::AbstractBasis{Ket}, col_basis::AbstractBasis{Bra})
	OperatorRep(coeff_func::Function, label_func::Function, b::AbstractBasis)

This last constructor allows an instance of `OperatorRep` can be constructed by defining the action 
of an abstract operator on the states of the basis (see below for examples).

__Methods and Examples__

We've already seen in the previous section that one can construct an operator representation by taking
the outer product of two state representations. Let's try instead using a functional construction of
an instance of `OperatorRep`. First, we'll define `eb` as an excitation basis:

	julia> xb = Basis("xb", [1:10])
	Basis{Ket} xb:
	| 1 ⟩
	| 2 ⟩
	| 3 ⟩
	| 4 ⟩
	| 5 ⟩
	| 6 ⟩
	| 7 ⟩
	| 8 ⟩
	| 9 ⟩
	| 10 ⟩  

Now, let's define the raising operator represented in `xb` as `r| n ⟩ = coeff_func(n) | label_func(n) ⟩ = sqrt(n+1) | n+1 ⟩`:

	julia> rrep = OperatorRep(n->sqrt(n[1]+1), n->n[1]+1, xb)

	OperatorRep:
	                     ⟨ 1 |  …               ⟨ 9 |           ⟨ 10 |
	  | 1 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 2 ⟩   1.41421+0.0im              0.0+0.0im       0.0+0.0im
	  | 3 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 4 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 5 ⟩       0.0+0.0im       …      0.0+0.0im       0.0+0.0im
	  | 6 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 7 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 8 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 9 ⟩       0.0+0.0im              0.0+0.0im       0.0+0.0im
	  | 10 ⟩      0.0+0.0im       …  3.16228+0.0im       0.0+0.0im

	julia> rrep*xb[2]
	StateRep{Ket} | #undef ; xb ⟩:
	     0.0+0.0im  | 1 ⟩
	     0.0+0.0im  | 2 ⟩
	 1.73205+0.0im  | 3 ⟩
	     0.0+0.0im  | 4 ⟩
	     0.0+0.0im  | 5 ⟩
	     0.0+0.0im  | 6 ⟩
	     0.0+0.0im  | 7 ⟩
	     0.0+0.0im  | 8 ⟩
	     0.0+0.0im  | 9 ⟩
	     0.0+0.0im  | 10 ⟩

	julia> rrep*xb[10]
	StateRep{Ket} | #undef ; xb ⟩:
	(all coefficients are zero)

In the second operation, `r| 10 ⟩ = sqrt(11) | 11 ⟩` produces an eigenvector
`| 11 ⟩` that, while an eigenstate of the abstract operator `r`, is not an
eigenvector of the matrix `rrep`. Thus, the design choice that QuantumJL be
limited to finite  dimensional subspaces of the Hilbert space is upheld in a
consistent manner.

Another feature of QuantumJL's `OperatorRep` implementation is a function 
that computes the partial trace, which is useful for performing entanglement
calculations:

	julia> q = StateRep(:q, normalize([1,1]), Basis("b", [0,1]))
	StateRep{Ket} | :q ; b ⟩:
	 0.707107+0.0im  | 0 ⟩
	 0.707107+0.0im  | 1 ⟩

	julia> qq = q*q
	StateRep{Ket} | :q,:q ; b ⊗ b ⟩:
	 0.5+0.0im  | 0,0 ⟩
	 0.5+0.0im  | 0,1 ⟩
	 0.5+0.0im  | 1,0 ⟩
	 0.5+0.0im  | 1,1 ⟩

	julia> qq[2:3] = 0
	0

	julia> normalize!(qq)
	StateRep{Ket} | :q,:q ; b ⊗ b ⟩:
	 0.707107+0.0im  | 0,0 ⟩
	      0.0+0.0im  | 0,1 ⟩
	      0.0+0.0im  | 1,0 ⟩
	 0.707107+0.0im  | 1,1 ⟩

	julia> qop = qq*qq'
	OperatorRep:
	                  ⟨ 0,0 |  …           ⟨ 1,0 |           ⟨ 1,1 |
	  | 0,0 ⟩  0.5+0.0im            0.0+0.0im         0.5+0.0im
	  | 0,1 ⟩  0.0+0.0im            0.0+0.0im         0.0+0.0im
	  | 1,0 ⟩  0.0+0.0im            0.0+0.0im         0.0+0.0im
	  | 1,1 ⟩  0.5+0.0im            0.0+0.0im         0.5+0.0im

	julia> ptrace(qop,1)
	OperatorRep:
	                ⟨ 0 |           ⟨ 1 |
	  | 0 ⟩  0.5+0.0im       0.0+0.0im
	  | 1 ⟩  0.0+0.0im       0.5+0.0im

	julia> trace(ptrace(qop,1)^2)
	0.5000000000000002 + 0.0im

6. Functions implemented in QuantumJL

The following is a list of functions that were implemented in 
QuantumJL. This list does *not* include overloaded functions 
like `filter`, `trace`, `get`, etc. 

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
	filtercoeffs, 
	filtercoeffs!,
	filterstates,
	filterstates!,
	samebasis,
	findstates,
	commutator,
	ptrace