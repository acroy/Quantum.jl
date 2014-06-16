The documentation is now out of date due to the recent redesign of the project. Once the 
new design has been more or less fully implemented, I'll come back and fix up the documentation.
<!-- Quantum.jl Documentation
===
1. Abstract Types
---
The following list contains the abstract types and type aliases referenced in this documentation:

	abstract Dirac
	abstract BraKet <: Dirac
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Dirac
	abstract AbstractState{K<:BraKet} <: Dirac
	abstract AbstractScalar <: Dirac
	typealias DiracCoeff Union(Number, AbstractScalar)


The not operator `!` can be applied to `Bra` and `Ket` types to alternate between 
the two:

	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

2. State 
--- 
__Description__

A `State` is an object that stores two labels. The first label is a unique identifier
for the state, while the second label specifies the basis for which the state is an
eigenstate. All instances of type `State` are explicit eigenstates of a basis (the 
representation of non-eigenstates is covered in the `DiracVector` section).
States are parameterized by their `kind`, which identifies whether the state
belongs to `Ket`-space or `Bra`-space.

The first label, as seen in the definition below, can be an instance of any type.

__Definition__

	immutable State{K<:BraKet, T} <: AbstractState{K}
	  label::T
	  basislabel::Symbol
	end

__Constructors__

	State{K<:BraKet,T}(label::T, basislabel::Symbol, kind::Type{K}=Ket) = State{kind, T}(label, basislabel)

__Methods and Examples__

Let's instantiate an array of eigenstates belonging to basis `:X`:

	julia> using Quantum

	julia> sv = [State(i,:X) for i=1:5]
	5-element Array{State{Ket,Int64},1}:
	 | 1:X ⟩
	 | 2:X ⟩
	 | 3:X ⟩
	 | 4:X ⟩
	 | 5:X ⟩

Taking the `ctranspose` of a state returns its dual:

	julia> sv[1]
	| 1:X ⟩

	julia> sv[1]'
	⟨ 1:X |

	julia> sv[1]'==State(1,:X,Bra)
	true

The inner product of a `Bra` state and a `Ket` state
of the same basis is a Kronecker delta function based 
on the states' labels:

	julia> sv[1]'*sv[1]
	1

	julia> sv[2]'*sv[1]
	0

The inner product of a `Bra` state and a `Ket` state
of different bases returns an `InnerProduct` object, 
which has its own section later on in this tutorial:

	julia> sv[1]'*State(1,:Y)
	⟨ 1:X |  1:Y ⟩

Similarly, the outer product of a `Ket` state and a `Bra`
state returns an `OuterProduct` object, whose behavior is
discussed in the `DiracMatrix` section:

	julia> sv[1]*State(1,:Y)'
	| 1:X ⟩⟨ 1:Y |

Finally, the product of two states with the same `kind` 
returns a `TensorState` of that `kind`:

	julia> sv[1]*State(1,:Y)
	| 1:X, 1:Y ⟩

	julia> sv[1]'*State(1,:Y)'
	⟨ 1:X, 1:Y |

For convenience, Quantum.jl provides a function called `statearr` which
takes in the following arguments and maps a `State` constructor to each 
of the elements of the Array accordingly:

	statearr{K<:BraKet, T}(arr::Array{T}, basislabel::Symbol, kind::Type{K}=Ket)

	julia> statearr([i+j for i=1:10,j=1:10], :X, Bra)
	10x10 Array{State{Bra,Int64},2}:
	 ⟨ 2:X |   ⟨ 3:X |   ⟨ 4:X |   ⟨ 5:X |   …  ⟨ 9:X |   ⟨ 10:X |  ⟨ 11:X |
	 ⟨ 3:X |   ⟨ 4:X |   ⟨ 5:X |   ⟨ 6:X |      ⟨ 10:X |  ⟨ 11:X |  ⟨ 12:X |
	 ⟨ 4:X |   ⟨ 5:X |   ⟨ 6:X |   ⟨ 7:X |      ⟨ 11:X |  ⟨ 12:X |  ⟨ 13:X |
	 ⟨ 5:X |   ⟨ 6:X |   ⟨ 7:X |   ⟨ 8:X |      ⟨ 12:X |  ⟨ 13:X |  ⟨ 14:X |
	 ⟨ 6:X |   ⟨ 7:X |   ⟨ 8:X |   ⟨ 9:X |      ⟨ 13:X |  ⟨ 14:X |  ⟨ 15:X |
	 ⟨ 7:X |   ⟨ 8:X |   ⟨ 9:X |   ⟨ 10:X |  …  ⟨ 14:X |  ⟨ 15:X |  ⟨ 16:X |
	 ⟨ 8:X |   ⟨ 9:X |   ⟨ 10:X |  ⟨ 11:X |     ⟨ 15:X |  ⟨ 16:X |  ⟨ 17:X |
	 ⟨ 9:X |   ⟨ 10:X |  ⟨ 11:X |  ⟨ 12:X |     ⟨ 16:X |  ⟨ 17:X |  ⟨ 18:X |
	 ⟨ 10:X |  ⟨ 11:X |  ⟨ 12:X |  ⟨ 13:X |     ⟨ 17:X |  ⟨ 18:X |  ⟨ 19:X |
	 ⟨ 11:X |  ⟨ 12:X |  ⟨ 13:X |  ⟨ 14:X |     ⟨ 18:X |  ⟨ 19:X |  ⟨ 20:X |

3. TensorState 
--- 
__Description__

A `TensorState` object is the result of a tensor products between `State`s.
__Definition__

	immutable TensorState{K<:BraKet} <: AbstractState{K}
	  states::Vector
	  TensorState{S<:State{K}}(v::Vector{S}) = new(v)
	  TensorState(v::Vector) = new(convert(Vector{State{K}},v))  
	end

__Constructors__

	TensorState{K<:BraKet}(labels::Vector, basislabel::Symbol, kind::Type{K}=Ket) = TensorState{kind}(statearr(labels, basislabel, kind))

__Methods and Examples__

The constructor above is provided for convenience and for it's similarity to 
the `State` constructor. Most of the time, however, one would construct `TensorState`s
using multiplication (`sv` refers to the `Vector{State{Ket,Int64}}` constructed in the 
previous section):

	julia> ts=sv[1]*State("Bob", :B)*sv[3]
	 | 1:X, "Bob":B, 3:X ⟩
	
One can take the `ctranspose` of `TensorState`s as well:

	julia> ts'
	⟨ 1:X, "Bob":B, 3:X |

The `separate` command returns the a vector of the components of a `TensorState`:

	julia> separate(ts)
	3-element Array{State{Ket,T},1}:
	 | 1:X ⟩
	 | "Bob":B ⟩
	 | 3:X ⟩

To get specific components, you can index into a `TensorState` as you would a
`Vector`:

	julia> ts[1]
	| 1:X ⟩

	julia> ts[2]
	| "Bob":B ⟩

	julia> ts[1:2]
	2-element Array{State{Ket,T},1}:
	 | 1:X ⟩
	 | "Bob":B ⟩

	julia> ts[:]
	3-element Array{State{Ket,T},1}:
	 | 1:X ⟩
	 | "Bob":B ⟩
	 | 3:X ⟩

Taking the inner product of `TensorState`s and `State`s can be ambiguous 
without specifying which component states are acting on each other. 
For example, the operation ⟨ a | a, b ⟩ is ambiguous; calculating it
as ⟨ a | a ⟩| b ⟩ or ⟨ a | b ⟩ | a ⟩ yields two different results. 

By default, Quantum.jl will always match component states to each other
by position, applying left to right. Thus, ⟨ a | a, b ⟩ is interpreted as 
⟨ a | a ⟩| b ⟩, and something like ⟨ a, c | a, b ⟩ is interpreted as 
⟨ c |(⟨ a | a, b ⟩) ->  ⟨ c | b ⟩:

	julia> ts'*ts
	1

	#sv[2]' -> ⟨ 2:X |
	julia> sv[2]'*ts
	0

	#ts' -> ⟨ 1:X, "Bob":B, 3:X |, prod(ts[1:2]) -> | 1:X, "Bob":B ⟩
	julia> ts'*prod(ts[1:2])
	⟨ 3:X |

	#prod(ts[[1,3,2]]) -> | 1:X, 3:X, "Bob":B ⟩
	julia> ts'*prod(ts[[1,3,2]])
	ScalarExpr(:(⟨ "Bob":B |  3:X ⟩ * ⟨ 3:X |  "Bob":B ⟩))

A `ScalarExpr` is an object that allows the user to transform `InnerProduct`s 
like numbers; thus they can be multiplied together, as ⟨ "Bob":B |  3:X ⟩ and 
⟨ 3:X |  "Bob":B ⟩ are above. `ScalarExpr`s are described (along with `InnerProduct`s) 
in their own section later on. 

To take inner products in a more explicit manner, one can use the `inner` function. 
The `inner` function takes a `Bra` state, `Ket` state, and an `Int` argument that
specifies index of the components on which the states act:

	#Default behavior
	julia> sv[3]'*ts #⟨ 3:X | 1:X, "Bob":B, 3:X ⟩ -> ⟨ 3:X | 1:X ⟩*| "Bob":B, 3:X ⟩ -> 0*| "Bob":B, 3:X ⟩ 
	0

	#using inner to specify which `Ket` state the `Bra` is acting on 
	julia> inner(sv[3]', ts, 3) #| 1:X, "Bob":B ⟩*⟨ 3:X | 3:X ⟩ -> | 1:X, "Bob":B ⟩*1
	| 1:X, "Bob":B ⟩ 
	
	#the integer argument always refers to an index of the TensorState being acted upon 
	julia> inner(ts', sv[3], 3) #⟨ 1:X, "Bob":B |*⟨ 3:X | 3:X ⟩ -> ⟨ 1:X, "Bob":B |*1
	⟨ 1:X, "Bob":B |

If one is taking an inner product between two `TensorState`s, the target of the indexing 
argument can be specified, and defaults to `Ket`:

	inner{K<:BraKet}(a::TensorState{Bra}, b::TensorState{Ket}, i::Int, target::Type{K}=Ket)


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

This last constructor allows an instance of `OperatorRep` to be constructed by defining the action 
of an abstract operator on the states of the basis (see below for examples).

__Methods and Examples__

We've already seen in the previous section that one can construct an operator representation by taking
the outer product of two state representations. Let's try constructing an `OperatorRep` functionally instead. First, we'll define `xb` as an excitation basis:

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
eigenvector of the matrix `rrep`. In this way, functionally constructed
instances of `OperatorRep` are consistent with the limitation that 
`Basis` objects be finite dimensional subspaces of an abstract infinite 
dimensional Hilbert space.

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

6. Functions
---

The following is a list of functions that are implemented in 
QuantumJL, not including overloaded functions 
like `filter`, `trace`, `get`, etc.:

	kind,
	statevec,
	tensor,
	statejoin,
	labeldelta,
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

The following is a list of functions imported from `Base` that are overloaded
for the types defined in QuantumJL:

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
	trace -->