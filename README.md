# Quantum.jl
*A framework for performing common quantum mechanical operations*

Julia is an up-and-coming language for scientific computing that promises a
high level of abstraction without the performance sacrifice such abstraction
usually entails. The capabilities of this new language provide an optimized
codebase for tackling various problems in physics, which are often
simultaneously conceptually difficult and computationally intensive.

Many similar technologies - Mathematica, MATLab, NumPy, to name a few - have
garnered support from the quantum physics community due to the wide variety of
third-party physics libraries available. There are usually at least one or two
projects for a given language that are dedicated to implementing quantum
mechanics operations using Dirac notation in a manner that is idiomatic to the
language.

If Julia is to properly compete with these other languages, it should have
such a library available; enter Quantum.jl.

In addition to providing the basic functionality that many quantum mechanical
implementations provide, this project provides a system for storing,
manipulating, and analyzing subspaces of the Hilbert space. In other words,
users can define their own bases, operators/states in terms of those bases, 
and perform complex combinatorial filtering and mapping on the labels of 
basis states. 

All calculations with states and operators resolve themselves using linear 
algebraic formulations whenever possible, and linear transformations act 
in accordance with user-defined bases. 

Quantum.jl currently only supports work in discrete spaces, but 
plans for implementing continuous space features are in the works 
as a post-release goal. Additionally, it has been suggested that 
Quantum.jl implement in an interface for other popular quantum 
libraries like QuTiP in order to more quickly offer up a host of 
advanced optimized features.

A discussion of the above, and elaboration on the project's plans
for the very near future, can be found [here](https://github.com/jrevels/Quantum.jl/issues/1).

Features
==========
Note: The whole package is currently undergoing a breaking refactor that
will hopefully make things even faster than before! In the meantime, 
all of the following are features that were implemented, but are in the
process of being reimplemented with the new design in mind.

The refactor is expected to be finished within the next two weeks, at which 
time the package will be registered and extensive documentation will be released.

###Current:
	- Basis, State, ScalarExpr, DiracVector, and DiracMatrix type implementations
	- Tensor product structure for Dirac objects
	- Application of arbitrary selection rules to extract subspaces from most types of 
	  Quantum.jl collection objects. Related methods:
	  	- filter
	  	- map/map!
	  	- mapmatch/mapmatch!
	  	- filtercoeffs
	 	- filterstates
	 	- find/findstates
	- The ability to index into Dirac objects using States, e.g. to find the coefficient
	  of a state in a DiracVector, one can simply call get(d::DiracVector, s::State) 
	- Arithmetic/Linear algebra operations that are logically consistent w.r.t. to user
	  defined bases. Examples: 
		- matrix/vector arithmetic for Dirac objects
		- inner/outer products implementations
		- idiomatic computation of expectation values/transition matrices
		- commutator of DiracMatrices
		- trace/partial trace, normalization, kron, and more
		- arithmetic operations between Dirac types and generic arrays
	- ScalarExpr objects for delayed calculation of mixed basis inner products. ScalarExpr
	  objects support all basic arithmetic operations. 
	- The ability to perform mixed basis calculation via the use of ScalarExprs and InnerProducts

###In The Works:
	- Refactor for optimization
	- Documentation
	- register Quantum.jl as julia package

###Post-Release Goals:
	- Interface with other similar libraries (eg. QuTiP)
	- Continuous space support
	- HDF5 support
	- Visualization system (probably via interfacing to more mature libraries)
	- incorporate units with [SIUnits.jl](https://github.com/Keno/SIUnits.jl) 
