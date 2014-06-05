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

Enter Quantum.jl.

In addition to providing the basic functionality that many quantum mechanical
implementations provide, this project provides a system for storing,
manipulating, and analyzing subspaces of the Hilbert space. In other words,
users can define their own bases, operators/states in terms of those bases, 
and perform complex combinatorial filtering and mapping on the labels of 
basis states. 

All calculations with states and operators resort to using linear algebraic 
formulations whenever possible, and linear transformations act in accordance 
with user-defined bases. 

All of this will become more clear once I get a chance to write up the docs;
stay tuned for updates in the next few weeks!

Features
==========

###Current:
	- Basis, State, Scalar, DiracVector, and DiracMatrix type implementations
	- Tensor product structure for Dirac objects
	- Application of arbitrary selection rules to extract subspaces from all types of 
	  Quantum.jl objects. Related methods:
	  	- filter
	  	- map/map!
	  	- mapmatch/mapmatch!
	  	- filtercoeffs
	 	- filterstates
	- The ability to index into Dirac objects using States, e.g. to find the coefficient
	  of a state in a DiracVector, one can simply call get(d::DiracVector, s::State) 
	- Arithmetic/Linear algebra operations that are logically consistent w.r.t. to user
	  defined bases. Examples: 
		- matrix/vector arithmetic for Dirac objects
		- inner/outer products implementations
		- idiomatic computation of expectation values/transition matrices
		- commutator of DiracMatrices
		- trace/partial trace
		- normalization
	- ScalarExpr objects for delayed calculation of mixed basis inner products. ScalarExpr
	  objects support all basic arithmetic operations. 
	- The ability to perform mixed basis calculation via the use of ScalarExprs and InnerProducts

###Upcoming:
	- Optimization/code cleaning
	- Documentation
	- HDF5 support
	- design a system for data visualization of common state properties/operations
	- package Quantum.jl as an open-source library that other Julia users can utilize