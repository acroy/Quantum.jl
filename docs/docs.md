*CURRENTLY A WORK IN PROGRESS*
Quantum.jl Documentation
===

##0. Abstract Types

The abstract type `Dirac` serves as the parent type for
all quantum objects - in other words, objects represented 
with Dirac notation in Quantum.jl.

In addition, the abstract type `AbstractScalar` has been 
defined as the parent type for `Dirac` objects that behave 
like numbers, such as inner products and arithmetic expressions
involving inner products (see `ScalarExpr` and `InnerProduct` types
for more information).

##1. States and Their Operations
###1.1 Eigenkets and Eigenbras
__Description__

Quantum.jl is built on a system of eigenstates - objects which 
are either of type `Bra` or `Ket`, and are parameterized
by a primitive that explicitly specifies a basis to which 
the object belongs. Here are the type definitions and 
constructors for `Bra`s and `Ket`s:

	abstract State{S} <: Dirac

	immutable Ket{b,T} <: State{Ket{b,T}}
		label::T
	end

	immutable Bra{b,T} <: State{Bra{b,T}}
		label::T
	end

	ket{T}(b,label::T) = Ket{b,T}(label)
	bra{T}(b,label::T) = Bra{b,T}(label)

The above definition allows for each state to store a basis
identifier as a type parameter, as well as a label as a unique 
identifier within the basis it belongs to. Because the type
of the label is parameterized, `Bra` and `Ket` objects are 
as compact as the object that serves as their label!

__Examples__

It's easy to construct `Ket`s using the above constructors:

	julia> using Quantum
	
	julia> ket(:X,1)
	| 1:X ⟩

	julia> typeof(ans)
	Ket{:X,Int64} (constructor with 1 method)

	julia> ket(:S,"1")
	| "1":S ⟩

	julia> typeof(ans)
	Ket{:S,ASCIIString} (constructor with 1 method)

	julia> ket(42,["1", 1, :f])
	| {"1",1,:f}:42 ⟩

	julia> typeof(ans)
	Ket{42,Array{Any,1}} (constructor with 1 method)

As you can see, both `Symbol`s and `Int`s can be used as basis identifiers, and
anything can be used as a label. `Bra`s are constructed in a similar manner:

	julia> bra(:X,1)
	⟨ 1:X |

	julia> typeof(ans)
	Bra{:X,Int64} (constructor with 1 method)

To get the dual of an eigenstate, use the `ctranspose` function (the 
`'` operator):

	julia> bra(:X,1)'
	| 1:X ⟩

	julia> bra(:X,1)''==bra(:X,1)
	true

###1.2 Tensor Product States

__Description__

States of the same "kind" (i.e. multiple `Bra`s, or multiple `Ket`s) 
can serve as factors of a tensor product state, the definition of 
which is below:

	immutable Tensor{S<:Union(Bra, Ket)} <: State{S}
		states::Vector{S}
		Tensor{K<:Ket}(v::Vector{K}) = new(v)
		Tensor{B<:Bra}(v::Vector{B}) = new(v)
	end 

As you can see, "kind" homogeneity is enforced via the inner constructors,
and a `Tensor` object is parameterized by the element type of the `Vector` that
stores its factors. 

__Examples__

To take the tensor product of multiple states, use the `tensor` function:

	julia> k,s = ket(:K,1), ket(:S,"a")
	(| 1:K ⟩,| "a":S ⟩)

	julia> t = tensor(k,k,k)
	| 1:K, 1:K, 1:K ⟩

	julia> typeof(ans)
	Tensor{Ket{:K,Int64}} (constructor with 2 methods)

	julia> tensor(s,t)
	| "a":S, 1:K, 1:K, 1:K ⟩

	julia> tensor(t,s)
	| 1:K, 1:K, 1:K, "a":S ⟩

	julia> typeof(ans)
	Tensor{Ket{b,T}} (constructor with 2 methods)

The same thing can be done with `Bra`s:

	julia> t = tensor(k',k')
	⟨ 1:K, 1:K |

	julia> typeof(ans)
	Tensor{Bra{:K,Int64}} (constructor with 2 methods)

	julia> tensor(t, s')
	⟨ 1:K, 1:K, "a":S |

	julia> typeof(ans)
	Tensor{Bra{b,T}} (constructor with 2 methods)

Mixing `Bra`s and `Ket`s will throw an error:

	julia> tensor(k,k')
	ERROR: KindMismatch: cannot perform tensor(| 1:K ⟩, ⟨ 1:K |)
	 in error at error.jl:21
	 in tensor at /Users/jarrettrevels/data/repos/quantum/src/state.jl:46

Finally, one can obtain the dual of a `Tensor` state
by using the `ctranspose` function (`'`):

	julia> t=tensor(k,s)
	| 1:K, "a":S ⟩

	julia> t'
	⟨ 1:K, "a":S |

	julia> t''==t
	true

###1.3 Other Operations on States

Now that we know how to construct states, let's explore how they operate.
Here is a list of some binary operations on states:
	
	inner{B<:Bra, K<:Ket}(a::State{B}, b::State{K})
		computes the discrete inner product of `a` and `b`
	kron(a::State, b::State)
		computes the kronecker product of `a` and `b`
	*(a::State, b::State)
		vector multiplication between `a` and `b`

As you might imagine, there are many other operations involving 
states - some, like addition, we'll touch on later in this 
documentation. Others can be found in [Quantum.jl's API](/docs/api.md).

In general, products with `State` objects work exactly analogously 
to how they work with actual vectors - which makes sense, since 
`Bra`s and `Ket`s are just abstract row and column vectors in the
first place.

For the examples that follow, we'll construct the following states
using `svec`, a convenience function provided by Quantum.jl:

	julia> xv=svec(:X, [1:3])
	3-element Array{Ket{:X,Int64},1}:
	 | 1:X ⟩
	 | 2:X ⟩
	 | 3:X ⟩

	julia> sv=svec(:S, ["$i" for i=1:3])
	3-element Array{Ket{:S,ASCIIString},1}:
	 | "1":S ⟩
	 | "2":S ⟩
	 | "3":S ⟩

	julia> tv=Tensor{Ket}[tensor(xv[i], sv[i]) for i=1:3]
	3-element Array{Tensor{Ket{b,T}},1}:
	 | 1:X, "1":S ⟩
	 | 2:X, "2":S ⟩
	 | 3:X, "3":S ⟩

__Inner Product (inner)__

Taking the inner product between a `State{B<:Bra}` and
a `State{K<:Ket}` results in a number if both are of 
the same basis (1 if they are duals of each other, 0 if they 
are not). If not, it resolves using Quantum.jl's built-in 
`InnerProduct`, `ScalarExpr` and `DiracVector` objects, which 
are all discussed in more detail in their respective sections 
later in this documentation.

	julia> inner(sv[1]',sv[2])
	0

	julia> inner(sv[1]',sv[1])
	1

	julia> inner(sv[1]',xv[1])
	⟨ "1":S |  1:X ⟩

	julia> typeof(ans)
	InnerProduct{Bra{:S,ASCIIString},Ket{:X,Int64}} (constructor with 1 method)

Inner products involving `Tensor` states can be ambiguous 
without specifying which factor states are acting on each other. 
For example, the operation ⟨ a | a, b ⟩ is ambiguous; calculating it
as ⟨ a | a ⟩| b ⟩ or ⟨ a | b ⟩ | a ⟩ yields two different results. 

By default, Quantum.jl will apply `Bra` states to the `Ket` states
one at a time, from left to right, resolving the product as it goes. 
Thus, ⟨ a | a, b ⟩ is interpreted as ⟨ a | a ⟩| b ⟩-> 1 | b ⟩:

	julia> inner(xv[1]',tv[1]) # ⟨ 1:X | 1:X, "1":S ⟩->⟨ 1:X | 1:X ⟩| "1":S ⟩
	1x1 DiracVector{Ket{:S,ASCIIString},Int64}
	 1  | "1":S ⟩

	julia> inner(tv[1]',xv[1]) # ⟨ 1:X, "1":S | 1:X ⟩->⟨ "1":S |⟨ 1:X | 1:X ⟩
	1x1 DiracVector{Bra{:S,ASCIIString},Int64}
	  ⟨ "1":S |
	 1

`DiracVector`s are the result of multiplying a coefficient
by a state (the coefficients, here, are both ⟨ 1:X | 1:X ⟩, which have resolved to 1). 

Following this logic, the operation ⟨ a, b | c, d ⟩ resolves as 
⟨ b |⟨ a | c, d ⟩ ->  ⟨ b |⟨ a | c ⟩| d ⟩ -> ⟨ a | c ⟩⟨ b | d ⟩.
This results in the expected behavior when taking the
inner product of two `Tensor` states from the same bases: 

	julia> inner(tv[1]',tv[1]) #⟨ 1:X, "1":S | 1:X, "1":S ⟩
	1
	julia> inner(tv[1]',tv[2]) #⟨ 1:X, "1":S | 2:X, "2":S ⟩
	0

...and the following behavior for `Tensor` states of unequal 
basis:

	julia> tb = tensor(bra(:A, "a"), bra(:B,"b"))
	⟨ "a":A, "b":B |

	julia> tk = tensor(ket(:C, "c"), ket(:D, "d"))
	| "c":C, "d":D ⟩

	julia> inner(tb,tk)
	ScalarExpr(:(⟨ "a":A |  "c":C ⟩ * ⟨ "b":B |  "d":D ⟩))

The `ScalarExpr` object returned above is Quantum.jl's mechanism 
for dealing with arithmetic operations involving `InnerProduct`s, 
which represent complex numbers (once again, more information
on the `ScalarExpr` and `InnerProduct` types can be found in their
respective section below).

Often times in quantum mechanics, one would like to specify 
the index of action for states acting on tensor product states, 
such as ⟨ a_2 | a_1, b_2 ⟩->⟨ a_2 | b_2 ⟩| a_1 ⟩

To accomplish this, Quantum.jl allows the user to pass 
position arguments to `inner` in order to specify which 
states are acting on which.

For example, recalling that a scalar times a state yields a `DiracVector{typeof(state), typeof(scalar)}`:

	julia> inner(xv[1]', tv[1], 2) # ⟨ (1:X)_2 | (1:X)_1, ("1":S)_2 ⟩
	1x1 DiracVector{Ket{:X,Int64},InnerProduct{Bra{:X,Int64},Ket{:S,ASCIIString}}}
	 ⟨ 1:X |  "1":S ⟩  | 1:X ⟩

In cases where one is taking the inner product between a `Tensor` state and a single
state, the index argument is always referring to the factor of the `Tensor` state
that the single state is to act upon:

	julia> inner(tv[1]', xv[1], 2) # ⟨ (1:X)_1, ("1":S)_2 | (1:X)_2 ⟩
	1x1 DiracVector{Bra{:X,Int64},InnerProduct{Bra{:S,ASCIIString},Ket{:X,Int64}}}
	 ⟨ 1:X |
	 ⟨ "1":S |  1:X ⟩

In the case of taking the inner product of two `Tensor` states, one
can specify both the index argument, and which state the index argument
is referring to by passing `Ket` or `Bra` as the fourth argument (the 
index targets the `Ket` state by default):

	julia> inner(tb, tk, 2) # (⟨ "a":A, "b":B |)_2 | ("c":C)_1, ("d":D)_2 ⟩->⟨ "a":A, "b":B | "d":D ⟩ | "c":C ⟩
	ScalarExpr(:(⟨ "a":A | "d":D ⟩ * ⟨ "b":B | "c":C ⟩))

	julia> inner(tb, tk, 2) == inner(tb, tk, 2, Ket) #Ket is targeted by default
	true

	julia> inner(tb, tk, 2, Bra) # ⟨ ("a":A)_1, ("b":B)_2 | (| ("c":C), ("d":D) ⟩)_2->⟨ "a":A | ⟨ "b":B | "c":C, "d":D ⟩
	ScalarExpr(:(⟨ "b":B | "c":C ⟩ * ⟨ "a":A | "d":D ⟩))

__Kronecker Product (kron)__

Taking the kronecker product of two states results in 
a `Tensor` state if the states are of the same kind, and
an `OuterProduct` object if not (`OuterProduct`s will be
covered more thoroughly in the DiracMatrix section):

	julia> kron(tv[1], tv[2]')
	| 1:X, "1":S ⟩⟨ 2:X, "2":S |

	julia> kron(tv[1], xv[2]')
	| 1:X, "1":S ⟩⟨ 2:X |

	julia> kron(sv[1], xv[2]')
	| "1":S ⟩⟨ 2:X |

	julia> kron(sv[1]', xv[2])
	| 2:X ⟩⟨ "1":S |

	julia> kron(sv[1], xv[2])
	| "1":S, 2:X ⟩

	julia> kron(tv[1]', tv[2]')
	⟨ 1:X, "1":S, 2:X, "2":S |

__Vector Multiplication (*)__

The `*` operator refers to standard vector multiplication. Thus,
`*(a::State{B<:Bra}, b::State{K<:Ket})` will return an `InnerProduct`,
and `*(a::State{K<:Ket}, b::State{B<:Bra})` will return an `OuterProduct`:

	julia> xv[1]'*xv[1]
	1

	julia> xv[1]'*xv[2]
	0

	julia> xv[1]'*sv[1]
	⟨ 1:X | "1":S ⟩

	julia> xv[1]*sv[1]'
	| 1:X ⟩⟨ "1":S |

Attempting to use `*` to do a tensor product will result in an error - for
that operation, you should use `tensor` or `kron` instead:

	julia> xv[1]*sv[1]
	ERROR: Multiplication Ket{:X,Int64}*Ket{:S,ASCIIString} is undefined. Perhaps you meant to use kron(Ket{:X,Int64}, Ket{:S,ASCIIString})?
	 in error at error.jl:21
	 in * at /Users/jarrettrevels/data/repos/quantum/src/misc.jl:23

	julia> kron(xv[1], sv[1])
	| 1:X, "1":S ⟩

<!--
##2. Collections of states: Bases and TensorBases

Let's try creating an array of states.

	julia> [ket(:F,i) for i=1:5]
	5-element Array{Ket{b,Int64},1}:
	 | 1:F ⟩
	 | 2:F ⟩
	 | 3:F ⟩
	 | 4:F ⟩
	 | 5:F ⟩

As you can see, Julia's type inference system doesn't specify a basis for the element type
of the array, despite the fact that all elements share the same basis. This can be 
a performance issue for large collections of states. In that case, one would want to 
use an explicit constructor to make sure the compiler can inference types correctly:

	julia> @time [ket(1,i) for i=1:100000]; #convenience constructor ket()
	elapsed time: 0.036017889 seconds (7183712 bytes allocated)

	julia> @time [Ket{1,Int64}(i) for i=1:100000]; #explicit constructor Ket{b,T}
	elapsed time: 6.5534e-5 seconds (800048 bytes allocated)

#####################################################################

Taking the inner product of `TensorState`s and `State`s can be ambiguous 
without specifying which component states are acting on each other. 
For example, the operation ⟨ a | a, b ⟩ is ambiguous; calculating it
as ⟨ a | a ⟩| b ⟩ or ⟨ a | b ⟩ | a ⟩ yields two different results. 

By default, Quantum.jl will always match component states to each other
by position, applying left to right. Thus, ⟨ a | a, b ⟩ is interpreted as 
⟨ a | a ⟩| b ⟩, and something like ⟨ a, c | a, b ⟩ is interpreted as 
⟨ c |(⟨ a | a, b ⟩) ->  ⟨ c |⟨ a | a ⟩| b ⟩ -> ⟨ c | b ⟩:

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