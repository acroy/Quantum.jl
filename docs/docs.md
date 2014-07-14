<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [Quantum.jl Documentation](#quantumjl-documentation)
  - [0. Abstract Types](#0-abstract-types)
  - [1. States and Their Products](#1-states-and-their-products)
    - [1.1 Eigenkets and Eigenbras](#11-eigenkets-and-eigenbras)
        - [Description](#description)
        - [Examples](#examples)
    - [1.2 Tensor Product States](#12-tensor-product-states)
        - [Description](#description-1)
        - [Examples](#examples-1)
    - [1.3 Other State Products: `kron`, `inner`, and `*`](#13-other-state-products-kron-inner-and-)
        - [Inner Product (`inner`)](#inner-product-inner)
        - [Kronecker Product (`kron`)](#kronecker-product-kron)
        - [Vector Multiplication (`*`)](#vector-multiplication-)
  - [2. `Basis` and `TensorBasis`](#2-basis-and-tensorbasis)
    - [2.1 `Basis`](#21-basis)
        - [Description](#description-2)
        - [Examples](#examples-2)
    - [2.2 `TensorBasis`](#22-tensorbasis)
        - [Description](#description-3)
        - [Examples](#examples-3)
  - [3. State Representations: `DiracVector`](#3-state-representations-diracvector)
    - [3.1 Introduction](#31-introduction)
        - [Description](#description-4)
        - [Examples](#examples-4)
    - [3.2  Array-like/Dict-like Operations on `DiracVector`s](#32--array-likedict-like-operations-on-diracvectors)
  - [4. Operators and Density Matrices](#4-operators-and-density-matrices)
    - [4.1 `OuterProduct`](#41-outerproduct)
        - [Description](#description-5)
        - [Examples](#examples-5)
    - [4.2 `DiracMatrix`](#42-diracmatrix)
        - [Description](#description-6)
        - [Examples](#examples-6)
    - [4.3  Array-like/Dict-like Operations on `DiracMatrix`](#43--array-likedict-like-operations-on-diracmatrix)
    - [4.4 Density Matrices and the Partial Trace](#44-density-matrices-and-the-partial-trace)
    - [4.5 A Note on Functions as Operators](#45-a-note-on-functions-as-operators)
  - [5. Abstract Scalars and `qeval`](#5-abstract-scalars-and-qeval)
    - [5.1 `InnerProduct`](#51-innerproduct)
    - [5.2 `ScalarExpr`](#52-scalarexpr)
        - [Description](#description-7)
        - [Examples](#examples-7)
        - [Supported Arithmetic using `ScalarExpr`](#supported-arithmetic-using-scalarexpr)
        - [`qeval` with `DiracVector` and `DiracMatrix`](#qeval-with-diracvector-and-diracmatrix)
  - [6. Fock Space Functions](#6-fock-space-functions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

*CURRENTLY A WORK IN PROGRESS*
Quantum.jl Documentation
===

The purpose of this documentation is to provide summaries and examples that
describe the core parts of Quantum.jl. In effect, it's mainly focused on
giving an overview of the types the package defines, how they behave in
relation to other Julia functions and types, and how they interact with each
other.

For more explicit information on the functions provided by Quantum.jl not
covered in this file, see [Quantum.jl's API](/docs/api.md).

##0. Abstract Types

The abstract type `Dirac` serves as the parent type for
all quantum objects - in other words, objects represented 
with Dirac notation in Quantum.jl.

The abstract type `AbstractScalar` has been defined as the parent 
type for `Dirac` objects that behave like numbers, such as the 
`ScalarExpr` and `InnerProduct` types.

The abstract type `State` has been defined as the parent 
type for `Dirac` objects that behave like abstract vectors, 
such as the `Bra`, `Ket`, and `Tensor` types. 

Finally, `AbstractOperator` is the parent type for `Dirac` objects
that behave like operators, such as `DiracMatrix` and `OuterProduct`
types.

##1. States and Their Products
###1.1 Eigenkets and Eigenbras
#####Description

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

The above definition allows for each state to store a basis identifier `b` as
a type parameter, as well as a label that serves as a unique identifier
within the basis it belongs to. Because the type of the label is
parameterized, `Bra` and `Ket` objects are as compact as the object that
serves as their label!

#####Examples

It's easy to construct `Ket`s using the above constructors:

	julia> using Quantum
	
	julia> ket(:N,1) # ket(b, label) -> Ket{b,typeof(label)}(label) -> Ket{:N,Int}(1)
	| 1:N ⟩

	julia> typeof(ans)
	Ket{:N,Int64} (constructor with 1 method)

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

	julia> bra(:N,1)
	⟨ 1:N |

	julia> typeof(ans)
	Bra{:N,Int64} (constructor with 1 method)

To get the dual of an eigenstate, use the `ctranspose` function (the 
`'` operator):

	julia> bra(:N,1)'
	| 1:N ⟩

	julia> bra(:N,1)''==bra(:N,1)
	true

###1.2 Tensor Product States

#####Description

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

#####Examples

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

One can index into a `Tensor` state just like one would a vector:

	julia> t = tensor(ket(:K,1),ket(:L,"2"),ket(:M, 3.0))
	| 1:K, "2":L, 3.0:M ⟩

	julia> t[2]
	| "2":L ⟩

	julia> t[1:end]
	3-element Array{Ket{b,T},1}:
	 | 1:K ⟩
	 | "2":L ⟩
	 | 3.0:M ⟩

Finally, one can obtain the dual of a `Tensor` state
by using the `ctranspose` function (`'`):

	julia> t=tensor(k,s)
	| 1:K, "a":S ⟩

	julia> t'
	⟨ 1:K, "a":S |

	julia> t''==t
	true

###1.3 Other State Products: `kron`, `inner`, and `*`

Now that we know how to construct states and tensor product 
states, let's explore how they operate on each other.
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

For use in the examples that follow, we'll construct the following 
states using `svec`, a convenience function provided by Quantum.jl:

	julia> xv=svec(:N, [1:3])
	3-element Array{Ket{:N,Int64},1}:
	 | 1:N ⟩
	 | 2:N ⟩
	 | 3:N ⟩

	julia> sv=svec(:S, ["$i" for i=1:3])
	3-element Array{Ket{:S,ASCIIString},1}:
	 | "1":S ⟩
	 | "2":S ⟩
	 | "3":S ⟩

	julia> tv=Tensor{Ket}[tensor(xv[i], sv[i]) for i=1:3]
	3-element Array{Tensor{Ket{b,T}},1}:
	 | 1:N, "1":S ⟩
	 | 2:N, "2":S ⟩
	 | 3:N, "3":S ⟩

#####Inner Product (`inner`)

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
	⟨ "1":S |  1:N ⟩

	julia> typeof(ans)
	InnerProduct (constructor with 1 method)

Inner products involving `Tensor` states can be ambiguous 
without specifying which factor states are acting on each other. 
For example, the operation ⟨ a | a, b ⟩ is ambiguous; calculating it
as ⟨ a | a ⟩| b ⟩ or ⟨ a | b ⟩ | a ⟩ yields two different results. 

By default, Quantum.jl will apply `Bra` states to the `Ket` states
one at a time, from left to right, resolving the product as it goes. 
Thus, ⟨ a | a, b ⟩ is interpreted as ⟨ a | a ⟩| b ⟩-> 1 | b ⟩:

	julia> inner(xv[1]',tv[1]) # ⟨ 1:N | 1:N, "1":S ⟩->⟨ 1:N | 1:N ⟩| "1":S ⟩
	1x1 DiracVector{Ket{:S,ASCIIString},Int64}
	 1  | "1":S ⟩

	julia> inner(tv[1]',xv[1]) # ⟨ 1:N, "1":S | 1:N ⟩->⟨ "1":S |⟨ 1:N | 1:N ⟩
	1x1 DiracVector{Bra{:S,ASCIIString},Int64}
	  ⟨ "1":S |
	 1

`DiracVector`s are the result of multiplying a coefficient
by a state (the coefficients, here, are both ⟨ 1:N | 1:N ⟩, which have resolved to 1). 

Following this logic, the operation ⟨ a, b | c, d ⟩ resolves as 
⟨ b |⟨ a | c, d ⟩ ->  ⟨ b |⟨ a | c ⟩| d ⟩ -> ⟨ a | c ⟩⟨ b | d ⟩.
This results in the expected behavior when taking the
inner product of two `Tensor` states from the same bases: 

	julia> inner(tv[1]',tv[1]) #⟨ 1:N, "1":S | 1:N, "1":S ⟩
	1
	julia> inner(tv[1]',tv[2]) #⟨ 1:N, "1":S | 2:N, "2":S ⟩
	0

...and the following behavior for `Tensor` states of unequal 
basis:

	julia> tb = tensor(bra(:A, "a"), bra(:B,"b"))
	⟨ "a":A, "b":B |

	julia> tk = tensor(ket(:C, "c"), ket(:D, "d"))
	| "c":C, "d":D ⟩

	julia> inner(tb,tk)
	ScalarExpr(:((1 * ⟨ "a":A | "c":C ⟩) * ⟨ "b":B | "d":D ⟩))

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

For example, recalling that a scalar times a state yields a `DiracVector`:

	julia> inner(xv[1]', tv[1], 2) # ⟨ (1:N)_2 | (1:N)_1, ("1":S)_2 ⟩
	1x1 DiracVector{Ket{:N,Int64},ScalarExpr}
	 ScalarExpr(:(1 * ⟨ 1:N | "1":S ⟩))  | 1:N ⟩

In cases where one is taking the inner product between a `Tensor` state and a single
state, the index argument is always referring to the factor of the `Tensor` state
that the single state is to act upon:

	julia> inner(tv[1]', xv[1], 2) # ⟨ (1:N)_1, ("1":S)_2 | (1:N)_2 ⟩
	1x1 DiracVector{Bra{:N,Int64},ScalarExpr}
	 ⟨ 1:N |
	 ScalarExpr(:(1 * ⟨ "1":S | 1:N ⟩))

In the case of taking the inner product of two `Tensor` states, one
can specify both the index argument, and which state the index argument
is referring to by passing `Ket` or `Bra` as the fourth argument (the 
index targets the `Ket` state by default):

	julia> inner(tb, tk, 2) # (⟨ "a":A, "b":B |)_2 | ("c":C)_1, ("d":D)_2 ⟩->⟨ "a":A, "b":B | "d":D ⟩ | "c":C ⟩
	ScalarExpr(:((1 * ⟨ "a":A | "d":D ⟩) * ⟨ "b":B | "c":C ⟩))

	julia> inner(tb, tk, 2) == inner(tb, tk, 2, Ket) #Ket is targeted by default
	true

	julia> inner(tb, tk, 2, Bra) # ⟨ ("a":A)_1, ("b":B)_2 | (| ("c":C), ("d":D) ⟩)_2->⟨ "a":A | ⟨ "b":B | "c":C, "d":D ⟩
	ScalarExpr(:((1 * ⟨ "b":B | "c":C ⟩) * ⟨ "a":A | "d":D ⟩))

#####Kronecker Product (`kron`)

Taking the kronecker product of two states results in 
a `Tensor` state if the states are of the same kind, and
an `OuterProduct` object if not (`OuterProduct`s will be
covered more thoroughly in the DiracMatrix section):

	julia> kron(tv[1], tv[2]')
	| 1:N, "1":S ⟩⟨ 2:N, "2":S |

	julia> kron(tv[1], xv[2]')
	| 1:N, "1":S ⟩⟨ 2:N |

	julia> kron(sv[1], xv[2]')
	| "1":S ⟩⟨ 2:N |

	julia> kron(sv[1]', xv[2])
	| 2:N ⟩⟨ "1":S |

	julia> kron(sv[1], xv[2])
	| "1":S, 2:N ⟩

	julia> kron(tv[1]', tv[2]')
	⟨ 1:N, "1":S, 2:N, "2":S |

#####Vector Multiplication (`*`)

The `*` operator refers to standard vector multiplication. Thus,
`*(a::State{B<:Bra}, b::State{K<:Ket})` will return an `InnerProduct`,
and `*(a::State{K<:Ket}, b::State{B<:Bra})` will return an `OuterProduct`:

	julia> xv[1]'*xv[1]
	1

	julia> xv[1]'*xv[2]
	0

	julia> xv[1]'*sv[1]
	⟨ 1:N | "1":S ⟩

	julia> xv[1]*sv[1]'
	| 1:N ⟩⟨ "1":S |

Attempting to use `*` to do a tensor product will result in an error - for
that operation, you should use `tensor` or `kron` instead:

	julia> xv[1]*sv[1]
	ERROR: Multiplication Ket{:N,Int64}*Ket{:S,ASCIIString} is undefined. Perhaps you meant to use kron(Ket{:N,Int64}, Ket{:S,ASCIIString})?
	 in error at error.jl:21
	 in * at /Users/jarrettrevels/data/repos/quantum/src/misc.jl:23

	julia> kron(xv[1], sv[1])
	| 1:N, "1":S ⟩


##2. `Basis` and `TensorBasis`

###2.1 `Basis`

#####Description
The linear representations provided by Quantum.jl rely on collections
of `State` objects that form `Basis` objects. The sections on `DiracVector`
and `DiracMatrix` types will make this apparent, but for now, let's get a 
feel for working with `Basis` objects themselves. 

A `Basis` object, as one might imagine, is a collection of eigenstates
that form a truncated, finite dimensional basis. Stored in a `Basis` 
object is a map that allows the indexing of representations by eigenstates 
themselves. This map is useful for spectral analysis, state filtering, and 
more. 

A `Basis` is parameterized by the kind of states that inhabit it (`Bra{b,T}` or `Ket{b,T}`). 
In order for Julia's type inference and multiple dispatch to be fully utilized, 
*all states in a basis must be of the same concrete type.* This ensures two things:

a) Eigenstates of a `Basis` all share a basis identifier. 

b) Eigenstates of a `Basis` all share the same label type. This greatly 
improves performance, and also promotes robust design in terms of how
one selects their quantum numbers. 

The third and final constraint on the states contained in a `Basis` is 
that the states must be unique. If they are not, Quantum.jl automatically
filters out the unique states in the collection provided and uses those.

#####Examples

One can use the `basis` function to construct `Basis` objects in much
the same way as we used `svec` in section 1.3:

	julia> b=basis(:N, [0:4])
	Basis{Ket{:N,Int64}}, 5 states:
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩

You can also construct a `Basis` by feeding the `basis` function
with states, or a `Vector` of states:

	julia> basis(ket(:S,"a"), ket(:S,"b"), ket(:S,"c"))
	Basis{Ket{:S,ASCIIString}}, 3 states:
	| "a":S ⟩
	| "b":S ⟩
	| "c":S ⟩

	julia> basis([ket(:S,"a"), ket(:S,"b"), ket(:S,"c")])
	Basis{Ket{:S,ASCIIString}}, 3 states:
	| "a":S ⟩
	| "b":S ⟩
	| "c":S ⟩
	
Just like states, `Basis` objects are transformed into their duals using
the `ctranspose` operator:

	julia> b'
	Basis{Bra{:N,Int64}}, 5 states:
	⟨ 0:N |
	⟨ 1:N |
	⟨ 2:N |
	⟨ 3:N |
	⟨ 4:N |

Many functions that are defined on standard `Vector`s are defined on `Basis` objects:

	julia> b[1]
	| 0:N ⟩

	julia> b[1:3] #note that this returns a Vector instead of a Basis 
	3-element Array{Ket{:N,Int64},1}:
	 | 0:N ⟩
	 | 1:N ⟩
	 | 2:N ⟩

	julia> filter(x->label(x)%2==0, b)
	Basis{Ket{:N,Int64}}, 3 states:
	| 0:N ⟩
	| 2:N ⟩
	| 4:N ⟩

	julia> map(x->ket(:N,label(x)^2), b)
	Basis{Ket{:N,Int64}}, 5 states:
	| 0:N ⟩
	| 1:N ⟩
	| 4:N ⟩
	| 9:N ⟩
	| 16:N ⟩

One can use `bcat` to join together bases and states:

	julia> bcat(b, ket(:N, 5))
	Basis{Ket{:N,Int64}}, 6 states:
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩
	| 5:N ⟩

	julia> bcat(ket(:N, 5), b)
	Basis{Ket{:N,Int64}}, 6 states:
	| 5:N ⟩
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩

	julia> bcat(basis(:N, [-4:-1]), b)
	Basis{Ket{:N,Int64}}, 9 states:
	| -4:N ⟩
	| -3:N ⟩
	| -2:N ⟩
	| -1:N ⟩
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩

	julia> bcat(b, basis(:N, [0:7]))
	Basis{Ket{:N,Int64}}, 11 states:
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩
	| 5:N ⟩
	| 6:N ⟩
	| 7:N ⟩

Note in the above example how concatenating two `Basis` objects resulted in a `Basis` 
containing only the unique states of both. 

You can also use `Basis` objects as dicitonaries, where the keys are
states and the values are their position in the collection:

	julia> b
	Basis{Ket{:N,Int64}}, 5 states:
	| 0:N ⟩
	| 1:N ⟩
	| 2:N ⟩
	| 3:N ⟩
	| 4:N ⟩

	julia> get(b, ket(:N, 3))
	4

	julia> b[ans]
	| 3:N ⟩

	julia> get(b, ket(:S, "a"), "not there!")
	"not there!"

###2.2 `TensorBasis`

#####Description

Multi-component descriptions of quantum systems often make use of a separable
tensor product structure. We've already introduced the tensor product state
type `Tensor`; Quantum.jl's tensor product basis type is the `TensorBasis`.

Taking the tensor product of `Basis` objects results in a `TensorBasis`,
which points to the original “separable” `Basis` objects in addition to
storing the `Tensor` states that resulted from taking the tensor product.
This allows for the optimization of operations on states and operators
represented in a `TensorBasis` (for example, taking the partial trace of a
`DiracMatrix`).

`TensorBasis` objects behave similarly to `Basis` objects, and the
restrictions on `TensorBasis` objects are extensions of the restrictions of
of their component bases. Each component `Basis` of a `TensorBasis` must,
obviously, be a valid basis in and of itself, and the constraint on state
uniqueness is similarly enforced for `TensorBasis` objects, as well.

#####Examples

The easiest and most efficient way to construct a `TensorBasis` is 
to pass `Basis` objects to the `tensor` function (using the `b`
basis from the examples in section 2.1):

	julia> tb = tensor(b,b,b)
	TensorBasis{Ket{:N,Int64},Basis{Ket{:N,Int64}}}, 125 states:
	| 0:N, 0:N, 0:N ⟩
	| 0:N, 0:N, 1:N ⟩
	| 0:N, 0:N, 2:N ⟩
	| 0:N, 0:N, 3:N ⟩
	| 0:N, 0:N, 4:N ⟩
	| 0:N, 1:N, 0:N ⟩
	| 0:N, 1:N, 1:N ⟩
	| 0:N, 1:N, 2:N ⟩
	| 0:N, 1:N, 3:N ⟩
	| 0:N, 1:N, 4:N ⟩
	⁞
	| 4:N, 2:N, 4:N ⟩
	| 4:N, 3:N, 0:N ⟩
	| 4:N, 3:N, 1:N ⟩
	| 4:N, 3:N, 2:N ⟩
	| 4:N, 3:N, 3:N ⟩
	| 4:N, 3:N, 4:N ⟩
	| 4:N, 4:N, 0:N ⟩
	| 4:N, 4:N, 1:N ⟩
	| 4:N, 4:N, 2:N ⟩
	| 4:N, 4:N, 3:N ⟩
	| 4:N, 4:N, 4:N ⟩

Or, you can just pass `Tensor` states to the `basis` function:

	julia> tb_partial=basis(tb[1], tb[10], tb[100])
	TensorBasis{Ket{:N,Int64},Basis{Ket{:N,Int64}}}, 3 states:
	| 0:N, 0:N, 0:N ⟩
	| 0:N, 1:N, 4:N ⟩
	| 3:N, 4:N, 4:N ⟩

`TensorBasis` objects handle all the same functions as `Basis`
objects (`bcat`, `map`, `filter`, etc.). In addition, a `TensorBasis`
can be separated into its component bases using `separate`:

	julia> separate(tb)
	3-element Array{Basis{Ket{:N,Int64}},1}:
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 1:N ⟩,| 2:N ⟩,| 3:N ⟩,| 4:N ⟩]
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 1:N ⟩,| 2:N ⟩,| 3:N ⟩,| 4:N ⟩]
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 1:N ⟩,| 2:N ⟩,| 3:N ⟩,| 4:N ⟩]

	julia> separate(tb_partial)
	3-element Array{Basis{Ket{:N,Int64}},1}:
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 3:N ⟩]
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 1:N ⟩,| 4:N ⟩]
	 Basis{Ket{:N,Int64}}[| 0:N ⟩,| 4:N ⟩]

As you can see, the component bases for `tb_partial` had to 
be inferred from the states provided. These bases each contain 
the minimum number of states such that their tensor 
product contains the states that make up the original
`TensorBasis` (in the case of `tb_partial`: `| 0:N, 0:N, 0:N ⟩`, 
`| 0:N, 1:N, 4:N ⟩`, `| 3:N, 4:N, 4:N ⟩`).

##3. State Representations: `DiracVector`

###3.1 Introduction 
#####Description 

A `DiracVector` is the representation of a state in a
basis. It stores the basis alongside a corresponding 
sparse vector of coefficients. 

A `DiracVector{S<:State, T<:Union(Number, ScalarExpr)}` 
inherits the state parameter of its associated basis, and 
is also parameterized by the element type of its coefficient array. 
The element type of the coefficient array must be either numeric
or `ScalarExpr`.

It is enforced that a `DiracVector` with a `Ket` basis has
a column vector of coefficients, while a `DiracVector` with a
`Bra` basis has a row vector of coefficients. 

The most natural way to think of a `DiracVector` is to 
treat it like a normal `Vector` in Julia; the only 
real difference is that the vector's basis is 
explicitly stored and transformed along with the
`Vector` itself. 

#####Examples

A natural way to construct a `DiracVector` is to
simply add states together (recall from section 
1.3 that a coefficient times a state will yield
a `DiracVector`):

	julia> dv=1/sqrt(2)*ket(:N,1) + 1/sqrt(2)*ket(:N,2)
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 0.707107  | 1:N ⟩
	 0.707107  | 2:N ⟩

Note that Quantum.jl does not automatically normalize `DiracVector`s.
This is to avoid unexpected behaviors programatically, and because 
testing for normalization is often a handy test for making sure one's 
math is correct. 

Quantum.jl instead provides a `normalize` method:

	julia> ket(:N,1) + ket(:N,2)
	2x1 DiracVector{Ket{:N,Int64},Int64}
	 1  | 1:N ⟩
	 1  | 2:N ⟩

	julia> dv=normalize(ans)
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 0.707107  | 1:N ⟩
	 0.707107  | 2:N ⟩ 

All normal linear algebraic arithmetic should work 
appropriately with `DiracVector`s, as well as taking
tensor products and inner products:

	julia> dv' # ctranspose of dv
	1x2 DiracVector{Bra{:N,Int64},Float64}
	  ⟨ 1:N |   ⟨ 2:N |
	 0.707107  0.707107

	 julia> dv+dv
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 1.41421  | 1:N ⟩
	 1.41421  | 2:N ⟩

	julia> dv-dv
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 0.0  | 1:N ⟩
	 0.0  | 2:N ⟩

	julia> bra(:N,2)*dv # ⟨ 2:N | dv ⟩
	0.7071067811865475

	julia> bra(:S,"a")*dv # ⟨ "a":S | dv ⟩
	ScalarExpr(:(0.7071067811865475 * ⟨ "a":S | 1:N ⟩ + 0.7071067811865475 * ⟨ "a":S | 2:N ⟩))

	julia> kron(dv,dv)
	4x1 DiracVector{Ket{:N,Int64},Float64}
	 0.5  | 1:N, 1:N ⟩
	 0.5  | 1:N, 2:N ⟩
	 0.5  | 2:N, 1:N ⟩
	 0.5  | 2:N, 2:N ⟩

	julia> kron(dv, dv') # | dv ⟩⟨ dv |
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:N,Int64},Float64}
	          ⟨ 1:N |   ⟨ 2:N |
	  | 1:N ⟩  0.5       0.5
	  | 2:N ⟩  0.5       0.5

	julia> kron(dv,ket(:N,1))
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 0.707107  | 1:N, 1:N ⟩
	 0.707107  | 2:N, 1:N ⟩

	julia> dv'*dv # ⟨ dv | dv ⟩
	0.9999999999999998

	julia> normalize(dv+ket(:N,3))
	3x1 DiracVector{Ket{:N,Int64},Float64}
	 0.5       | 1:N ⟩
	 0.5       | 2:N ⟩
	 0.707107  | 3:N ⟩

	julia> normalize(dv-ket(:N,3))
	3x1 DiracVector{Ket{:N,Int64},Float64}
	  0.5       | 1:N ⟩
	  0.5       | 2:N ⟩
	 -0.707107  | 3:N ⟩

	julia> ket(:N,1)+ket(:N,1)+ket(:N,1)
	1x1 DiracVector{Ket{:N,Int64},Int64}
	 3  | 1:N ⟩

One cannot add states that aren't allowed
to be part of the same basis:

	julia> ket(:N,1)+ket(:S,"1")
	ERROR: cannot construct basis: all states must have same label type and same basis identifier
	 in basis at /Users/jarrettrevels/data/repos/quantum/src/basis.jl:23
	 in + at /Users/jarrettrevels/data/repos/quantum/src/diracvector.jl:278

	julia> ket(:N,1)+ket(:N,"1")
	ERROR: cannot construct basis: all states must have same label type and same basis identifier
	 in basis at /Users/jarrettrevels/data/repos/quantum/src/basis.jl:23
	 in + at /Users/jarrettrevels/data/repos/quantum/src/diracvector.jl:278

It is possible to use `InnerProduct`s and `ScalarExpr`s as 
coefficient (though the former are always converted to the 
latter for the purpose of generalizing future calculations):

	julia> (bra(:S,"1")*ket(:N, 1)) * ket(:N,1)
	1x1 DiracVector{Ket{:N,Int64},ScalarExpr}
	 ScalarExpr(:(1 * ⟨ "1":S | 1:N ⟩))  | 1:N ⟩

	julia> ans+ans
	1x1 DiracVector{Ket{:N,Int64},ScalarExpr}
	 ScalarExpr(:(1 * ⟨ "1":S | 1:N ⟩ + 1 * ⟨ "1":S | 1:N ⟩))  | 1:N ⟩

	julia> kron(ans,ans)
	1x1 DiracVector{Ket{:N,Int64},ScalarExpr}
	 ScalarExpr(:((1 * ⟨ "1":S | 1:N ⟩ + 1 * ⟨ "1":S | 1:N ⟩) * (1 * ⟨ "1":S | 1:N ⟩ + 1 * ⟨ "1":S | 1:N ⟩))) | 1:N, 1:N ⟩

###3.2  Array-like/Dict-like Operations on `DiracVector`s

As one might expect, repeatedly adding states one at a time
is generally not the most efficient way to go about things. 

The function `dvec` is provided to allow preallocation of
a `DiracVector`'s basis and coefficient array.

	julia> dv=dvec(Array(Float64, 4), basis(:N, [1:4]))
	4x1 DiracVector{Ket{:N,Int64},Float64}
	 6.94464e-310  | 1:N ⟩
	 6.94464e-310  | 2:N ⟩
	 6.94463e-310  | 3:N ⟩
	 0.0           | 4:N ⟩

The empty `Array{Float64}` instantiated above could have been any kind of 
column vector of length 4; Quantum.jl automatically converts a DiracVector's 
coefficient array to a `SparseMatrixCSC` for internal storage. 

Many array-like functions are supported on the coefficient array, 
including `getindex` and `setindex!` (for the sake of example, 
we won't bother with normalization):

	julia> dv[1:4] = [1:4];

	julia> dv
	4x1 DiracVector{Ket{:N,Int64},Float64}
	 1.0  | 1:N ⟩
	 2.0  | 2:N ⟩
	 3.0  | 3:N ⟩
	 4.0  | 4:N ⟩

Another is `map`:

	julia> map(x->x^2, dv)
	4x1 DiracVector{Ket{:N,Int64},Float64}
	  1.0  | 1:N ⟩
	  4.0  | 2:N ⟩
	  9.0  | 3:N ⟩
	 16.0  | 4:N ⟩

Three filter methods are provided as well; `filterstates`,
`filtercoeffs`, and `filternz`. The first uses the predicate
to filter by states, the second uses the predicate to 
filter by coeffs, and the third merely filters out zeros
and the states associated with those zeros:

	julia> filtercoeffs(x->x%2==0, dv)
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 2.0  | 2:N ⟩
	 4.0  | 4:N ⟩

	julia> filterstates(s->label(s)==3 || label(s)==1, dv)
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 1.0  | 1:N ⟩
	 3.0  | 3:N ⟩

	julia> dv[1:3] = 0;

	julia> dv
	4x1 DiracVector{Ket{:N,Int64},Float64}
	 0.0  | 1:N ⟩
	 0.0  | 2:N ⟩
	 0.0  | 3:N ⟩
	 4.0  | 4:N ⟩

	julia> filternz(dv)
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 4.0  | 4:N ⟩

Using `get`, you can retrieve the coefficient
associated with a state by doing the following:

	julia> get(dv, ket(:N,4))
	4.0

	julia> get(dv, ket(:N,3))
	0.0

	julia> get(dv, ket(:S,"a"))
	ERROR: key not found: | "a":S ⟩
	 in get at /Users/jarrettrevels/data/repos/quantum/src/diracvector.jl:76

	julia> get(dv, ket(:S,"a"), "not there!")
	"not there!"

##4. Operators and Density Matrices

###4.1 `OuterProduct`

#####Description

Taking the outer product of two states results in
an `OuterProduct` object, which is essentially 
Quantum.jl's projection operator. `OuterProduct`s 
form the building blocks of operators represented 
as `DiracMatrix` objects. 

#####Examples

Like states, `OuterProduct`s have a consistent
behavior for taking products, conjugate transposes,
and performing linear algebraic arithmetic:

	julia> o = ket(:N,1)*bra(:S,"a")
	| 1:N ⟩⟨ "a":S |

	julia> o'
	| "a":S ⟩⟨ 1:N |

	julia> o*ket(:S,"a") #| 1:N ⟩⟨ "a":S | "a":S ⟩
	1x1 DiracVector{Ket{:N,Int64},Int64}
	 1  | 1:N ⟩

	julia> bra(:N,1)*o #⟨ 1:N | 1:N ⟩⟨ "a":S |
	1x1 DiracVector{Bra{:S,ASCIIString},Int64}
	  ⟨ "a":S |
	 1

	julia> bra(:G,:g)*o #⟨ :g:G | 1:N ⟩⟨ "a":S |
	1x1 DiracVector{Bra{:S,ASCIIString},ScalarExpr}
	 ⟨ "a":S |
	 ScalarExpr(:(1 * ⟨ :g:G | 1:N ⟩))

	julia> o*o
	1x1 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},ScalarExpr}
	         ⟨ "a":S |
	  | 1:N ⟩  ScalarExpr(:(1 * ⟨ "a":S | 1:N ⟩))


	julia> o-ket(:N,2)*bra(:S,"b")
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |    ⟨ "b":S |
	  | 1:N ⟩  1.0          0.0
	  | 2:N ⟩  0.0         -1.0

	julia> ans-o
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |    ⟨ "b":S |
	  | 1:N ⟩  0.0          0.0
	  | 2:N ⟩  0.0         -1.0

	julia> kron(o,ket(:G,:g))
	| 1:N, :g:G ⟩⟨ "a":S |

	julia> kron(o,bra(:G,:g))
	| 1:N ⟩⟨ "a":S, :g:G |

	julia> kron(o,ket(:G,:g)*bra(:H,234.0232))
	| 1:N, :g:G ⟩⟨ "a":S, 234.0232:H |

	julia> kron(ket(:G,:g)*bra(:H,234.0232),o)
	| :g:G, 1:N ⟩⟨ 234.0232:H, "a":S |

###4.2 `DiracMatrix`

#####Description

A `DiracMatrix` serves as a representation of 
an operator in a basis of `OuterProduct`s, just 
like a `DiracVector` serves as a representation 
of a general state in a basis of eigenstates.

A `DiracMatrix{K<:Ket, B<:Bra, T<:Union(Number, ScalarExpr)}` 
explicitly stores its "ket" basis, "bra" basis, and a coefficient array. 

#####Examples

Constructing a `DiracMatrix` can be done simply 
by adding `OuterProducts`:

	julia> dm=ket(:N,1)*bra(:S,"a")+ket(:N,2)*bra(:S,"b")
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |   ⟨ "b":S |
	  | 1:N ⟩  1.0         0.0
	  | 2:N ⟩  0.0         1.0

These outer products must all map from one single basis to 
another (i.e. no basis mixing is allowed):

	julia> ket(:N,1)*bra(:S,"a")+ket(:N,2)*bra(:G,:g)
	ERROR: no method +(OuterProduct{Ket{:N,Int64},Bra{:S,ASCIIString}}, OuterProduct{Ket{:N,Int64},Bra{:G,Symbol}})

	julia> ket(:N,1)*bra(:S,"a")+ket(:N,2)*bra(:S,123)
	ERROR: no method +(OuterProduct{Ket{:N,Int64},Bra{:S,ASCIIString}}, OuterProduct{Ket{:N,Int64},Bra{:S,Int64}})

Just like `DiracVector`s, standard linear algebraic arithmetic involving
other `Dirac` objects is built-in:

	julia> dm+dm
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |   ⟨ "b":S |
	  | 1:N ⟩  2.0         0.0
	  | 2:N ⟩  0.0         2.0

	julia> dm-dm
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |   ⟨ "b":S |
	  | 1:N ⟩  0.0         0.0
	  | 2:N ⟩  0.0         0.0

	julia> dm+ket(:N,3)*bra(:S,"c")
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a":S |   ⟨ "b":S |   ⟨ "c":S |
	  | 1:N ⟩  1.0         0.0         0.0
	  | 2:N ⟩  0.0         1.0         0.0
	  | 3:N ⟩  0.0         0.0         1.0

	julia> ket(:N,0)*bra(:S,"z")-dm
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "z":S |    ⟨ "a":S |    ⟨ "b":S |
	  | 0:N ⟩  1.0          0.0          0.0
	  | 1:N ⟩  0.0         -1.0          0.0
	  | 2:N ⟩  0.0          0.0         -1.0
	  
...`kron` works as well:

	julia> kron(ket(:N,0)*bra(:S,"z"),dm)
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	               ⟨ "z":S, "a":S |   ⟨ "z":S, "b":S |
	  | 0:N, 1:N ⟩  1.0                0.0
	  | 0:N, 2:N ⟩  0.0                1.0

	julia> kron(dm,ket(:N,3)*bra(:S,"c"))
	2x2 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	               ⟨ "a":S, "c":S |   ⟨ "b":S, "c":S |
	  | 1:N, 3:N ⟩  1.0                0.0
	  | 2:N, 3:N ⟩  0.0                1.0

	julia> kron(dm,dm)
	4x4 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	               ⟨ "a":S, "a":S |   ⟨ "a":S, "b":S |  …   ⟨ "b":S, "b":S |
	  | 1:N, 1:N ⟩  1.0                0.0                   0.0
	  | 1:N, 2:N ⟩  0.0                1.0                   0.0
	  | 2:N, 1:N ⟩  0.0                0.0                   0.0
	  | 2:N, 2:N ⟩  0.0                0.0                   1.0

	julia> kron(dm, ket(:L,1.2))
	2x2 DiracMatrix{Ket{b,T},Bra{:S,ASCIIString},Float64}
	                 ⟨ "a":S |   ⟨ "b":S |
	  | 1.2:L, 1:N ⟩  1.0         0.0
	  | 1.2:L, 2:N ⟩  0.0         1.0

	julia> kron(dm, bra(:L,1.2))
	2x2 DiracMatrix{Ket{:N,Int64},Bra{b,T},Float64}
	          ⟨ 1.2:L, "a":S |   ⟨ 1.2:L, "b":S |
	  | 1:N ⟩  1.0                0.0
	  | 2:N ⟩  0.0                1.0

...as does `*`:

	julia> dm*ket(:S,"a")
	2x1 DiracVector{Ket{:N,Int64},Float64}
	 1.0  | 1:N ⟩
	 0.0  | 2:N ⟩

	julia> bra(:N,2)*dm
	1x2 DiracVector{Bra{:S,ASCIIString},Float64}
	  ⟨ "a":S |   ⟨ "b":S |
	 0.0         1.0

	julia> bra(:N,2)*dm*ket(:S,"b")
	1.0

	julia> bra(:N,1)*dm*ket(:S,"a")
	1.0

	julia> dm*(ket(:S,"a")*bra(:G,:g))
	2x1 DiracMatrix{Ket{:N,Int64},Bra{:G,Symbol},Float64}
	          ⟨ :g:G |
	  | 1:N ⟩  1.0
	  | 2:N ⟩  0.0

	julia> (ket(:G,:g)*bra(:N,2))*dm
	1x2 DiracMatrix{Ket{:G,Symbol},Bra{:S,ASCIIString},Float64}
	           ⟨ "a":S |   ⟨ "b":S |
	  | :g:G ⟩  0.0         1.0

Mixed basis products resolve themselves appropriately using `InnerProduct`s
and `ScalarExpr`s:

	julia> dm*(ket(:A,"a")*bra(:B,"b"))
	2x1 DiracMatrix{Ket{:N,Int64},Bra{:B,ASCIIString},ScalarExpr}
	         ⟨ "b":B |
	  | 1:N ⟩  ScalarExpr(:(1 * ⟨ "a":S | "a":A ⟩ + 1 * 0))
	  | 2:N ⟩  ScalarExpr(:(1 * 0 + 1 * ⟨ "b":S | "a":A ⟩))

Quantum.jl offers the function `dmat` for non-arithmetic construction
of a `DiracMatrix` (it's basically the same as `dvec` for a `DiracMatrix`):

	julia> dmat(eye(3), basis(:N, [1:3]), basis(:S, ["a,","b","c"])')
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Float64}
	          ⟨ "a,":S |   ⟨ "b":S |   ⟨ "c":S |
	  | 1:N ⟩  1.0          0.0         0.0
	  | 2:N ⟩  0.0          1.0         0.0
	  | 3:N ⟩  0.0          0.0         1.0

Note in the `dmat` call above that the first basis argument is a `Ket` basis, while 
the second basis argument is a `Bra` basis. 

When provided with a single basis argument, `dmat` works on the assumption that 
the "row" basis and "column" basis are duals of each other:

	julia> dmat(eye(3), basis(:N, [1:3]))
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:N,Int64},Float64}
	          ⟨ 1:N |   ⟨ 2:N |   ⟨ 3:N |
	  | 1:N ⟩  1.0       0.0       0.0
	  | 2:N ⟩  0.0       1.0       0.0
	  | 3:N ⟩  0.0       0.0       1.0

It is also possible to use functional definitions to construct a `DiracMatrix`
by passing function arguments to `dmat`. For example, consider the lowering operator 
`a`, which can be defined as 

	`a | n:N ⟩ = eigval(| n:N ⟩) eigstate(| n:N ⟩) = sqrt(n) | n-1:N ⟩`.

Here, we split the functions that define `a`'s action on the `N` basis by defining 
the function returning the eigeinvalue and the function returning the eigenstate 
separately. These functions separate can be passed to `dmat` as anonymous functions,
as well as the basis that you wish the operator to be defined on (currently, this
method of constructing `DiracMatrix` objects only works when the "ket" basis and 
"bra" basis are duals of each other):

	julia> a = DiracMatrix(n->sqrt(label(n)), n->ket(:N,label(n)-1), basis(:N, [1:10]))
	10x10 DiracMatrix{Ket{:N,Int64},Bra{:N,Int64},Float64}
	           ⟨ 1:N |   ⟨ 2:N |   ⟨ 3:N |  …   ⟨ 8:N |   ⟨ 9:N |   ⟨ 10:N |
	  | 1:N ⟩   0.0       1.41421   0.0          0.0       0.0       0.0
	  | 2:N ⟩   0.0       0.0       1.73205      0.0       0.0       0.0
	  | 3:N ⟩   0.0       0.0       0.0          0.0       0.0       0.0
	  | 4:N ⟩   0.0       0.0       0.0          0.0       0.0       0.0
	  | 5:N ⟩   0.0       0.0       0.0       …  0.0       0.0       0.0
	  | 6:N ⟩   0.0       0.0       0.0          0.0       0.0       0.0
	  | 7:N ⟩   0.0       0.0       0.0          2.82843   0.0       0.0
	  | 8:N ⟩   0.0       0.0       0.0          0.0       3.0       0.0
	  | 9:N ⟩   0.0       0.0       0.0          0.0       0.0       3.16228
	  | 10:N ⟩  0.0       0.0       0.0       …  0.0       0.0       0.0

	julia> filternz(a*ket(:N,2)) #using filternz to remove all zeros from output DiracVector
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 1.41421  | 1:N ⟩

	julia> filternz(a*ket(:N,3))
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 1.73205  | 2:N ⟩

	julia> filternz(a'*ket(:N,2)) # a' | 2:N ⟩
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 1.73205  | 3:N ⟩

	julia> filternz(a'*ket(:N,3))
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 2.0  | 4:N ⟩

	julia> a'*ket(:N,10) #takes us out of the basis we defined the operator on, result is zero vector
	10x1 DiracVector{Ket{:N,Int64},Float64}
	 0.0  | 1:N ⟩
	 0.0  | 2:N ⟩
	 0.0  | 3:N ⟩
	 0.0  | 4:N ⟩
	 0.0  | 5:N ⟩
	 0.0  | 6:N ⟩
	 0.0  | 7:N ⟩
	 0.0  | 8:N ⟩
	 0.0  | 9:N ⟩
	 0.0  | 10:N ⟩

###4.3  Array-like/Dict-like Operations on `DiracMatrix`

`DiracMatrix` objects have a bunch of array functions defined on them, just like `DiracVector`s.
Those functions (`map`, `getindex`, `setindex!`, etc.) all work as expected. 

The `get` function can still be used to find coefficients associated with
specific states, by passing individual `Ket`s and `Bra`s or whole
`OuterProduct`s:

	julia> dm = dmat(reshape([1:9],3,3), basis(:N, [1:3]), basis(:S, ["a,","b","c"])')
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Int64}
	          ⟨ "a,":S |   ⟨ "b":S |   ⟨ "c":S |
	  | 1:N ⟩  1            4           7
	  | 2:N ⟩  2            5           8
	  | 3:N ⟩  3            6           9

	julia> get(dm, ket(:N,2),bra(:S, "c"))
	8

	julia> get(dm, ket(:N,2),bra(:S, "c"),0)
	8

	julia> get(dm, ket(:N,2),bra(:G, :g),0)
	0

For `DiracMatrix` objects, `get` can also return `DiracVector`s:

	julia> get(dm, ket(:N,2))
	1x3 DiracVector{Bra{:S,ASCIIString},Int64}
	  ⟨ "a,":S |   ⟨ "b":S |   ⟨ "c":S |
	 2            5           8

	julia> get(dm, bra(:S,"b"))
	3x1 DiracVector{Ket{:N,Int64},Int64}
	 4  | 1:N ⟩
	 5  | 2:N ⟩
	 6  | 3:N ⟩

###4.4 Density Matrices and the Partial Trace

Here is an example of constructing a Bell state, then constructing a density matrix, 
and finally an entanglement calculation using `ptrace`:

	julia> q = 1/sqrt(2)*ket(:Q,0) + 1/sqrt(2)*ket(:Q,1)
	2x1 DiracVector{Ket{:Q,Int64},Float64}
	 0.707107  | 0:Q ⟩
	 0.707107  | 1:Q ⟩

	julia> s = 1/sqrt(2)*ket(:S,"a") + 1/sqrt(2)*ket(:S,"b")
	2x1 DiracVector{Ket{:S,ASCIIString},Float64}
	 0.707107  | "a":S ⟩
	 0.707107  | "b":S ⟩

	julia> qs = kron(q,s)
	4x1 DiracVector{Ket{b,T},Float64}
	 0.5  | 0:Q, "a":S ⟩
	 0.5  | 0:Q, "b":S ⟩
	 0.5  | 1:Q, "a":S ⟩
	 0.5  | 1:Q, "b":S ⟩

	julia> qs[2:3] = 0; qs=normalize(qs) #now we have a Bell state
	4x1 DiracVector{Ket{b,T},Float64}
	 0.707107  | 0:Q, "a":S ⟩
	 0.0       | 0:Q, "b":S ⟩
	 0.0       | 1:Q, "a":S ⟩
	 0.707107  | 1:Q, "b":S ⟩

	julia> qsdm=qs*qs'
	4x4 DiracMatrix{Ket{b,T},Bra{b,T},Float64}
	                 ⟨ 0:Q, "a":S |  …   ⟨ 1:Q, "a":S |   ⟨ 1:Q, "b":S |
	  | 0:Q, "a":S ⟩  0.5                 0.0              0.5
	  | 0:Q, "b":S ⟩  0.0                 0.0              0.0
	  | 1:Q, "a":S ⟩  0.0                 0.0              0.0
	  | 1:Q, "b":S ⟩  0.5                 0.0              0.5

	julia> ptrace(qsdm, 1) #trace out the first system
	2x2 DiracMatrix{Ket{:S,ASCIIString},Bra{:S,ASCIIString},Any}
	            ⟨ "a":S |   ⟨ "b":S |
	  | "a":S ⟩  0.5         0
	  | "b":S ⟩  0           0.5

	julia> trace(ans^2)
	0.4999999999999998

Note that I didn't need to use two separate bases here - I could have
done `kron(q,q)` just as easily. I merely chose to use the `:S` basis
to help illustrate the behavior of `ptrace`.

###4.5 A Note on Functions as Operators

As an alternative to `DiracMatrix` objects, one could just define a Julia
function on types of `State`s. This allows one to set up the input->output
definition of an operator without having to pay the storage costs of a
`DiracMatrix` for very large bases. For example:

	julia> lower(s::Ket{:N,Int}) = sqrt(label(s)) * ket(:N,label(s)-1)
	lower (generic function with 1 method)

	julia> lower(ket(:N,10))
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 3.16228  | 9:N ⟩

	julia> lower(ket(:N,1000))
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 31.6228  | 999:N ⟩

	julia> lower(ket(:N,10000000))
	1x1 DiracVector{Ket{:N,Int64},Float64}
	 3162.28  | 9999999:N ⟩

If you wanted operator-like behavior, you could overload the function for `Bra`s 
as well:

	#mimicking the action of the lowering operator on bras
	julia> lower(s::Bra{:N, Int}) = sqrt(label(s)+1) * bra(:N, label(s)+1)
	lower (generic function with 2 methods)

	julia> bra(:N,2)*a
	1x10 DiracVector{Bra{:N,Int64},Float64}
	  ⟨ 1:N |   ⟨ 2:N |   ⟨ 3:N |   ⟨ 4:N |  …   ⟨ 8:N |   ⟨ 9:N |   ⟨ 10:N |
	 0.0       0.0       1.73205   0.0          0.0       0.0       0.0

	julia> lower(bra(:N,2))
	1x1 DiracVector{Bra{:N,Int64},Float64}
	  ⟨ 3:N |
	 1.73205

The obvious downside to this is that Julia functions don't behave
like matrices/operators in the normal sense, so it's not possible to
*idiomatically* carry out multiplicative operations or operations that 
rely on a tensor product structure. In the future, Quantum.jl will 
hopefully provide an `Operator` type or something similar that
will allow such idiomatic functional behavior. 

##5. Abstract Scalars and `qeval`

This section explains two types of objects that have 
already been featured extensively in previous examples:
`InnerProduct` and `ScalarExpr`.

###5.1 `InnerProduct`

As seen in previous examples, an `InnerProduct` is the 
result of multiplying a `Bra` and `Ket` of two different bases:

	julia> i=bra(:N, 1)*ket(:S, "a")
	⟨ 1:N | "a":S ⟩

Since an `InnerProduct` is a representation of a complex number,
one can take it's complex conjugate:

	julia> i'
	⟨ "a":S | 1:N ⟩

It can also be used as a coeffcient for `DiracVector`s and `DiracMatrix`s:

	julia> i*ket(:G,:g)
	1x1 DiracVector{Ket{:G,Symbol},ScalarExpr}
	 ScalarExpr(:(1 * ⟨ 1:N | "a":S ⟩))  | :g:G ⟩

	julia> i*ket(:G,:g)*bra(:F,"f")
	1x1 DiracMatrix{Ket{:G,Symbol},Bra{:F,ASCIIString},ScalarExpr}
	          ⟨ "f":F |
	  | :g:G ⟩  ScalarExpr(:(1 * ⟨ 1:N | "a":S ⟩))

As you can see, an `InnerProduct` will automatically be converted to a `ScalarExpr` when it is
used as a coefficient of a `Dirac` object. This is mainly because all numeric values can be explicitly
converted to type `ScalarExpr`, while the same is not true of `InnerProduct`. This distinction is especially 
significant for use with sparse matrices, as `zero(::Type{ScalarExpr})` is well-defined, while 
`zero(::Type{ScalarExpr})` is not.

The function `qeval` can be used to evaluate `InnerProduct`s by passing it 
an arbitrary function of the form `f(b::(B<:Bra), k::(K<:Ket)) = (N<:Number)`:

	julia> f(b,k) = label(b) + label(k)*im
	f (generic function with 1 method)

	julia> qeval(f, bra(:N,1)*ket(:G,3))
	1 + 3im

	julia> qeval(f, bra(:N,5)*ket(:G,1.23))
	5.0 + 1.23im

###5.2 `ScalarExpr`

#####Description
`ScalarExpr` objects are used in Quantum.jl to represent calculations 
involving numbers and `InnerProduct`s. As its name implies, a `ScalarExpr` is
simply a container of an `Expr` that would evaluate to a number if all `InnerProduct`s
contained in the `Expr` are evaluated to numbers.

Evaluation of operations on a `ScalarExpr` is almost entirely lazy, with the exception of
a few idempotent or easily invertible functions (e.g. `-`, `abs`, `conj`).

#####Examples

Any numeric quantity or `InnerProduct` can be converted to a `ScalarExpr` using `scalar`:

	julia> scalar(1)
	ScalarExpr(:(1 * 1))

	julia> scalar(0)
	ScalarExpr(:(1 * 0))

	julia> scalar(3+4.2im)
	ScalarExpr(:(1 * 3.0 + 4.2im))

	julia> scalar(bra(:N,1)*ket(:G,:g))
	ScalarExpr(:(1 * ⟨ 1:N | :g:G ⟩))

`ScalarExpr` objects can be used in a variety of arithmetic 
calculations involving numbers and `InnerProduct`s:

	julia> abs(3*ans + ans^2 * ans/(bra(:K,:k)*ket(:H,"h")))
	ScalarExpr(:(abs(3 * (1 * ⟨ 1:N | :g:G ⟩) + ((1 * ⟨ 1:N | :g:G ⟩)^2 * (1 * ⟨ 1:N | :g:G ⟩)) / ⟨ :k:K | "h":H ⟩)))

`qeval` can be used to evaluate the `InnerProduct`s in the a `ScalarExpr`:

	julia> qeval((b,k)->2.341, ans) #evaluate all to the same number for trivial checking
	12.503281000000001

	julia> abs(3*2.341 + 2.341^2 * 2.341/(2.341))
	12.503281000000001

#####Supported Arithmetic using `ScalarExpr`

It is fairly easy to add `ScalarExpr` support for simple functions. 
If you desire support to be added for a function not on the 
below list, feel free to submit an issue or pull request.

List of supported arithmetic operations on `ScalarExpr`:

	^
	*
	/
	+
	-
	abs
	exp
	conj,ctranspose

#####`qeval` with `DiracVector` and `DiracMatrix`

Consider the following:
	
	julia> dm
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Int64}
	          ⟨ "a,":S |   ⟨ "b":S |   ⟨ "c":S |
	  | 1:N ⟩  1            4           7
	  | 2:N ⟩  2            5           8
	  | 3:N ⟩  3            6           9

Multiplying `dm` by itself is a mixed-basis calculation, and will result in `ScalarExpr`s 
(as you can see, formatting for long `ScalarExpr`s isn't exactly pretty):

	julia> dm*dm # | :N_i ⟩⟨ :S_j | :N_i ⟩⟨ :S_j |
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},ScalarExpr}
	         …  ⟨ "c":S |
	  | 1:N ⟩     ScalarExpr(:(((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + ((1 * 1.0 + 7 * ⟨ "a,":S | 1:N ⟩) - 1)) + 8 * ⟨ "a,":S | 2:N ⟩) - 1)) + 9 * ⟨ "a,":S | 3:N ⟩) - 1)) + 28 * ⟨ "b":S | 1:N ⟩) - 1)) + 32 * ⟨ "b":S | 2:N ⟩) - 1)) + 36 * ⟨ "b":S | 3:N ⟩) - 1)) + 49 * ⟨ "c":S | 1:N ⟩) - 1)) + 56 * ⟨ "c":S | 2:N ⟩) - 1)) + 63 * ⟨ "c":S | 3:N ⟩) - 1))
	  | 2:N ⟩     ScalarExpr(:(((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + 1 * 0.0) + 14 * ⟨ "a,":S | 1:N ⟩) - 1)) + 16 * ⟨ "a,":S | 2:N ⟩) - 1)) + 18 * ⟨ "a,":S | 3:N ⟩) - 1)) + 35 * ⟨ "b":S | 1:N ⟩) - 1)) + 40 * ⟨ "b":S | 2:N ⟩) - 1)) + 45 * ⟨ "b":S | 3:N ⟩) - 1)) + 56 * ⟨ "c":S | 1:N ⟩) - 1)) + 64 * ⟨ "c":S | 2:N ⟩) - 1)) + 72 * ⟨ "c":S | 3:N ⟩) - 1))
	  | 3:N ⟩     ScalarExpr(:(((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + (((1 + 1 * 0.0) + 21 * ⟨ "a,":S | 1:N ⟩) - 1)) + 24 * ⟨ "a,":S | 2:N ⟩) - 1)) + 27 * ⟨ "a,":S | 3:N ⟩) - 1)) + 42 * ⟨ "b":S | 1:N ⟩) - 1)) + 48 * ⟨ "b":S | 2:N ⟩) - 1)) + 54 * ⟨ "b":S | 3:N ⟩) - 1)) + 63 * ⟨ "c":S | 1:N ⟩) - 1)) + 72 * ⟨ "c":S | 2:N ⟩) - 1)) + 81 * ⟨ "c":S | 3:N ⟩) - 1))

We can use `qeval` to evaluate the above mess numerically (once again, I'm choosing to evaluate all inner products to a single number for
illustrative purposes):

	julia> qeval((b,k)->1, ans)
	3x3 DiracMatrix{Ket{:N,Int64},Bra{:S,ASCIIString},Int64}
	            ⟨ "a,":S |     ⟨ "b":S |     ⟨ "c":S |
	  | 1:N ⟩   72            180           288
	  | 2:N ⟩   90            225           360
	  | 3:N ⟩  108            270           432

The same kind of thing could be done with a `DiracVector` if one were so inclined:

	julia> dv = dvec([1:5], basis(:N, [1:5]))
	5x1 DiracVector{Ket{:N,Int64},Int64}
	 1  | 1:N ⟩
	 2  | 2:N ⟩
	 3  | 3:N ⟩
	 4  | 4:N ⟩
	 5  | 5:N ⟩

	julia> (bra(:K,"k")*dv)*dv # ⟨ "k":K | dv:N ⟩| dv:N ⟩
	5x1 DiracVector{Ket{:N,Int64},ScalarExpr}
	 ScalarExpr(:((((1 * ⟨ "k":K | 1:N ⟩ + 2 * ⟨ "k":K | 2:N ⟩) + 3 * ⟨ "k":K | 3:N ⟩) + 4 * ⟨ "k":K | 4:N ⟩) + 5 * ⟨ "k":K | 5:N ⟩))        …  | 1:N ⟩
	 ScalarExpr(:(((((1 * ⟨ "k":K | 1:N ⟩ + 2 * ⟨ "k":K | 2:N ⟩) + 3 * ⟨ "k":K | 3:N ⟩) + 4 * ⟨ "k":K | 4:N ⟩) + 5 * ⟨ "k":K | 5:N ⟩) * 2))     | 2:N ⟩
	 ScalarExpr(:(((((1 * ⟨ "k":K | 1:N ⟩ + 2 * ⟨ "k":K | 2:N ⟩) + 3 * ⟨ "k":K | 3:N ⟩) + 4 * ⟨ "k":K | 4:N ⟩) + 5 * ⟨ "k":K | 5:N ⟩) * 3))     | 3:N ⟩
	 ScalarExpr(:(((((1 * ⟨ "k":K | 1:N ⟩ + 2 * ⟨ "k":K | 2:N ⟩) + 3 * ⟨ "k":K | 3:N ⟩) + 4 * ⟨ "k":K | 4:N ⟩) + 5 * ⟨ "k":K | 5:N ⟩) * 4))     | 4:N ⟩
	 ScalarExpr(:(((((1 * ⟨ "k":K | 1:N ⟩ + 2 * ⟨ "k":K | 2:N ⟩) + 3 * ⟨ "k":K | 3:N ⟩) + 4 * ⟨ "k":K | 4:N ⟩) + 5 * ⟨ "k":K | 5:N ⟩) * 5))     | 5:N ⟩

	julia> qeval((b,k)->1, ans)
	5x1 DiracVector{Ket{:N,Int64},Int64}
	 15  | 1:N ⟩
	 30  | 2:N ⟩
	 45  | 3:N ⟩
	 60  | 4:N ⟩
	 75  | 5:N ⟩

##6. Fock Space Functions

For convenience, Quantum.jl provides a variety of functions associated with 
constructing and manipulating bases of Fock states. 

The first is `fbasis`:

	julia> fbasis(5) #5 levels, single particle
	Basis{Ket{:F,Int64}}, 5 states: 
	| 0:F ⟩
	| 1:F ⟩
	| 2:F ⟩
	| 3:F ⟩
	| 4:F ⟩

	julia> fbasis(5,5) #5 levels, 5 particles
	TensorBasis{Ket{:F,Int64},Basis{Ket{:F,Int64}}}, 3125 states:
	| 0:F, 0:F, 0:F, 0:F, 0:F ⟩
	| 0:F, 0:F, 0:F, 0:F, 1:F ⟩
	| 0:F, 0:F, 0:F, 0:F, 2:F ⟩
	| 0:F, 0:F, 0:F, 0:F, 3:F ⟩
	| 0:F, 0:F, 0:F, 0:F, 4:F ⟩
	| 0:F, 0:F, 0:F, 1:F, 0:F ⟩
	| 0:F, 0:F, 0:F, 1:F, 1:F ⟩
	| 0:F, 0:F, 0:F, 1:F, 2:F ⟩
	| 0:F, 0:F, 0:F, 1:F, 3:F ⟩
	| 0:F, 0:F, 0:F, 1:F, 4:F ⟩
	⁞
	| 4:F, 4:F, 4:F, 2:F, 4:F ⟩
	| 4:F, 4:F, 4:F, 3:F, 0:F ⟩
	| 4:F, 4:F, 4:F, 3:F, 1:F ⟩
	| 4:F, 4:F, 4:F, 3:F, 2:F ⟩
	| 4:F, 4:F, 4:F, 3:F, 3:F ⟩
	| 4:F, 4:F, 4:F, 3:F, 4:F ⟩
	| 4:F, 4:F, 4:F, 4:F, 0:F ⟩
	| 4:F, 4:F, 4:F, 4:F, 1:F ⟩
	| 4:F, 4:F, 4:F, 4:F, 2:F ⟩
	| 4:F, 4:F, 4:F, 4:F, 3:F ⟩
	| 4:F, 4:F, 4:F, 4:F, 4:F ⟩

Next is `fvec`:

	julia> fvec(5) #5 levels, | 0:F ⟩
	5x1 DiracVector{Ket{:F,Int64},Float64}
	 1.0  | 0:F ⟩
	 0.0  | 1:F ⟩
	 0.0  | 2:F ⟩
	 0.0  | 3:F ⟩
	 0.0  | 4:F ⟩

	julia> fvec(5,1) #5 levels, | 1:F ⟩
	5x1 DiracVector{Ket{:F,Int64},Float64}
	 0.0  | 0:F ⟩
	 1.0  | 1:F ⟩
	 0.0  | 2:F ⟩
	 0.0  | 3:F ⟩
	 0.0  | 4:F ⟩

	julia> fvec(5,2) #5 levels, | 2:F ⟩
	5x1 DiracVector{Ket{:F,Int64},Float64}
	 0.0  | 0:F ⟩
	 0.0  | 1:F ⟩
	 1.0  | 2:F ⟩
	 0.0  | 3:F ⟩
	 0.0  | 4:F ⟩

You can make creation and annihilation operators using
`fcreate` and `fdestroy`:

	julia> fcreate(3) #3 levels
	3x3 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	          ⟨ 0:F |   ⟨ 1:F |   ⟨ 2:F |
	  | 0:F ⟩  0.0       0.0       0.0
	  | 1:F ⟩  1.0       0.0       0.0
	  | 2:F ⟩  0.0       1.41421   0.0

	julia> fcreate(3,3,2) #3 levels, 3 particles, acting on particle 2
	27x27 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	                    ⟨ 0:F, 0:F, 0:F |  …   ⟨ 2:F, 2:F, 2:F |
	  | 0:F, 0:F, 0:F ⟩  0.0                    0.0
	  | 0:F, 0:F, 1:F ⟩  0.0                    0.0
	  | 0:F, 0:F, 2:F ⟩  0.0                    0.0
	  | 0:F, 1:F, 0:F ⟩  1.0                    0.0
	  | 0:F, 1:F, 1:F ⟩  0.0                 …  0.0
	  | 0:F, 1:F, 2:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 0:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 1:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 2:F ⟩  0.0                    0.0
	 ⋮                                       ⋱
	  | 1:F, 2:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 0:F, 0:F ⟩  0.0                    0.0
	  | 2:F, 0:F, 1:F ⟩  0.0                 …  0.0
	  | 2:F, 0:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 0:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 1:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 2:F, 0:F ⟩  0.0                 …  0.0
	  | 2:F, 2:F, 1:F ⟩  0.0                    0.0
	  | 2:F, 2:F, 2:F ⟩  0.0                    0.0

	julia> fdestroy(3)  #3 levels
	3x3 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	          ⟨ 0:F |   ⟨ 1:F |   ⟨ 2:F |
	  | 0:F ⟩  0.0       1.0       0.0
	  | 1:F ⟩  0.0       0.0       1.41421
	  | 2:F ⟩  0.0       0.0       0.0

	julia> fdestroy(3,3,3) #3 levels, 3 particles, acting on particle 3
	27x27 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	                    ⟨ 0:F, 0:F, 0:F |  …   ⟨ 2:F, 2:F, 2:F |
	  | 0:F, 0:F, 0:F ⟩  0.0                    0.0
	  | 0:F, 0:F, 1:F ⟩  0.0                    0.0
	  | 0:F, 0:F, 2:F ⟩  0.0                    0.0
	  | 0:F, 1:F, 0:F ⟩  0.0                    0.0
	  | 0:F, 1:F, 1:F ⟩  0.0                 …  0.0
	  | 0:F, 1:F, 2:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 0:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 1:F ⟩  0.0                    0.0
	  | 0:F, 2:F, 2:F ⟩  0.0                    0.0
	 ⋮                                       ⋱
	  | 1:F, 2:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 0:F, 0:F ⟩  0.0                    0.0
	  | 2:F, 0:F, 1:F ⟩  0.0                 …  0.0
	  | 2:F, 0:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 0:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 1:F ⟩  0.0                    0.0
	  | 2:F, 1:F, 2:F ⟩  0.0                    0.0
	  | 2:F, 2:F, 0:F ⟩  0.0                 …  0.0
	  | 2:F, 2:F, 1:F ⟩  0.0                    1.41421
	  | 2:F, 2:F, 2:F ⟩  0.0                    0.0

The identity and number operators are also provided:

	julia> feye(3)
	3x3 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	          ⟨ 0:F |   ⟨ 1:F |   ⟨ 2:F |
	  | 0:F ⟩  1.0       0.0       0.0
	  | 1:F ⟩  0.0       1.0       0.0
	  | 2:F ⟩  0.0       0.0       1.0

	julia> fnum(3)
	3x3 DiracMatrix{Ket{:F,Int64},Bra{:F,Int64},Float64}
	          ⟨ 0:F |   ⟨ 1:F |   ⟨ 2:F |
	  | 0:F ⟩  0.0       0.0       0.0
	  | 1:F ⟩  0.0       1.0       0.0
	  | 2:F ⟩  0.0       0.0       2.0
