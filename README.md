# Quantum.jl
*A framework for performing operations on user-defined quantum states*

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
such a library available.

Enter my capstone project: Quantum.jl.

In addition to providing the basic functionality that many quantum mechanical
implementations provide, my capstone will provide a system for storing,
manipulating, and analyzing subspaces of the Hilbert space. In other words,
users will be able to define their own bases, and define operators and
states in terms of those bases, and perform complex combinatorial filtering
and mapping on the labels of basis states.

Features
==========

###Current:
	- Basis, State, Scalar, DiracVector, and DiracMatrix type implementations
	- Tensor product structure for Quantum objects
	- Application of arbitrary selection rules to extract subspaces from all types of 
	  Quantum.jl objects. Related methods:
	  	- filter
	  	- map/map!
	  	- mapmatch/mapmatch!
	  	- filtercoeffs
	 	- filterstates
	- The ability to index into Quantum Objects using States, e.g. to find the coefficient
	  of a state in a DiracVector, one can simply call get(d::DiracVector, s::State) 
	- Arithmetic/Linear algebra operations that are logically consistent w.r.t. to the 
	  bases that Quantum objects are represented in. Examples: 
		- matrix/vector arithmetic with Quantum objects
		- inner/outer products of Bras and Kets
		- computing expectation values/transition matrices
		- commutator of operators
		- trace/partial trace
		- normalization
	- ScalarExpr objects for holding off calculation of mixed basis inner products. ScalarExpr
	  objects support almost all basic arithmetic operations.

###Upcoming:
	- Mixed basis operations (implement DiracSum)
	- HDF5 support
	- design a system for data visualization of common state properties/operations
	- package Quantum.jl as an open-source library that other Julia users can utilize