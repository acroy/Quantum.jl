Quantum Mechanics in Julia
===================

Julia is an up-and-coming language for scientific computing that promises a
high-level of abstraction without the performance sacrifice such abstraction
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
	- basis, state, state representation, and operator type implementations
	- tensor product structure for bases and states/state representations
	- application of arbitrary selection rules to extract subspaces from all types of 
	  Quantum.jl objects (filtering and mapping)
	- basic arithmetic operations (including linear algebra operations). Examples: 
		- normalization
		- inner/outer products of Bras and Kets
		- expectation values/transition matrices
		- commutator of operators
		- operator trace
	- partial trace of operator representations (entanglement calculations) 

###Upcoming:
	- Abstract symbolic manipulation/simplification for mixed basis operations (implement QuantumExpr)
	- HDF5 support
	- design a system for data visualization of common state properties/operations
	- package Quantum.jl as an open-source library that other Julia users can utilize