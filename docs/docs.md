*CURRENTLY A WORK IN PROGRESS*
Quantum.jl Documentation
===

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

__Inner Product (`inner`)__

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

__Kronecker Product (`kron`)__

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

__Vector Multiplication (`*`)__

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


##2. `Basis` and `TensorBasis`

###2.1 `Basis`

__Description__
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

__Examples__

One can use the `basis` function to construct `Basis` objects in much
the same way as we used `svec` in section 1.3:

	julia> b=basis(:X, [0:4])
	Basis{Ket{:X,Int64}}, 5 states:
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩

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
	Basis{Bra{:X,Int64}}, 5 states:
	⟨ 0:X |
	⟨ 1:X |
	⟨ 2:X |
	⟨ 3:X |
	⟨ 4:X |

Many functions that are defined on standard `Vector`s are defined on `Basis` objects:

	julia> b[1]
	| 0:X ⟩

	julia> b[1:3] #note that this returns a Vector instead of a Basis 
	3-element Array{Ket{:X,Int64},1}:
	 | 0:X ⟩
	 | 1:X ⟩
	 | 2:X ⟩

	julia> filter(x->label(x)%2==0, b)
	Basis{Ket{:X,Int64}}, 3 states:
	| 0:X ⟩
	| 2:X ⟩
	| 4:X ⟩

	julia> map(x->ket(:X,label(x)^2), b)
	Basis{Ket{:X,Int64}}, 5 states:
	| 0:X ⟩
	| 1:X ⟩
	| 4:X ⟩
	| 9:X ⟩
	| 16:X ⟩

One can use `bcat` to join together bases and states:

	julia> bcat(b, ket(:X, 5))
	Basis{Ket{:X,Int64}}, 6 states:
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩
	| 5:X ⟩

	julia> bcat(ket(:X, 5), b)
	Basis{Ket{:X,Int64}}, 6 states:
	| 5:X ⟩
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩

	julia> bcat(basis(:X, [-4:-1]), b)
	Basis{Ket{:X,Int64}}, 9 states:
	| -4:X ⟩
	| -3:X ⟩
	| -2:X ⟩
	| -1:X ⟩
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩

	julia> bcat(b, basis(:X, [0:7]))
	Basis{Ket{:X,Int64}}, 11 states:
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩
	| 5:X ⟩
	| 6:X ⟩
	| 7:X ⟩

Note in the above example how concatenating two `Basis` objects resulted in a `Basis` 
containing only the unique states of both. 

You can also use `Basis` objects as dicitonaries, where the keys are
states and the values are their position in the collection:

	julia> b
	Basis{Ket{:X,Int64}}, 5 states:
	| 0:X ⟩
	| 1:X ⟩
	| 2:X ⟩
	| 3:X ⟩
	| 4:X ⟩

	julia> get(b, ket(:X, 3))
	4

	julia> b[ans]
	| 3:X ⟩

	julia> get(b, ket(:S, "a"), "not there!")
	"not there!"

###2.2 `TensorBasis`

__Description__

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

__Examples__

The easiest and most efficient way to construct a `TensorBasis` is 
to pass `Basis` objects to the `tensor` function (using the `b`
basis from the examples in section 2.1):

	julia> tb = tensor(b,b,b)
	TensorBasis{Ket{:X,Int64},Basis{Ket{:X,Int64}}}, 125 states:
	| 0:X, 0:X, 0:X ⟩
	| 0:X, 0:X, 1:X ⟩
	| 0:X, 0:X, 2:X ⟩
	| 0:X, 0:X, 3:X ⟩
	| 0:X, 0:X, 4:X ⟩
	| 0:X, 1:X, 0:X ⟩
	| 0:X, 1:X, 1:X ⟩
	| 0:X, 1:X, 2:X ⟩
	| 0:X, 1:X, 3:X ⟩
	| 0:X, 1:X, 4:X ⟩
	⁞
	| 4:X, 2:X, 4:X ⟩
	| 4:X, 3:X, 0:X ⟩
	| 4:X, 3:X, 1:X ⟩
	| 4:X, 3:X, 2:X ⟩
	| 4:X, 3:X, 3:X ⟩
	| 4:X, 3:X, 4:X ⟩
	| 4:X, 4:X, 0:X ⟩
	| 4:X, 4:X, 1:X ⟩
	| 4:X, 4:X, 2:X ⟩
	| 4:X, 4:X, 3:X ⟩
	| 4:X, 4:X, 4:X ⟩

Or, you can just pass `Tensor` states to the `basis` function:

	julia> tb_partial=basis(tb[1], tb[10], tb[100])
	TensorBasis{Ket{:X,Int64},Basis{Ket{:X,Int64}}}, 3 states:
	| 0:X, 0:X, 0:X ⟩
	| 0:X, 1:X, 4:X ⟩
	| 3:X, 4:X, 4:X ⟩

`TensorBasis` objects handle all the same functions as `Basis`
objects (`bcat`, `map`, `filter`, etc.). In addition, a `TensorBasis`
can be separated into its component bases using `separate`:

	julia> separate(tb)
	3-element Array{Basis{Ket{:X,Int64}},1}:
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 1:X ⟩,| 2:X ⟩,| 3:X ⟩,| 4:X ⟩]
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 1:X ⟩,| 2:X ⟩,| 3:X ⟩,| 4:X ⟩]
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 1:X ⟩,| 2:X ⟩,| 3:X ⟩,| 4:X ⟩]

	julia> separate(tb_partial)
	3-element Array{Basis{Ket{:X,Int64}},1}:
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 3:X ⟩]
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 1:X ⟩,| 4:X ⟩]
	 Basis{Ket{:X,Int64}}[| 0:X ⟩,| 4:X ⟩]

As you can see, the component bases for `tb_partial` had to 
be inferred from the states provided. 

<!--
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