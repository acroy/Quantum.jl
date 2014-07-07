# Quantum.jl

Quantum.jl is a toolset for Julia designed around the efficient and idiomatic
manipulation of abstract quantum objects, so that one can easily transfer
reasoning about quantum systems from the chalkboard to the computer without
getting too bogged down in implementation details. Almost all quantum
mechanics libraries work solely with vector/matrix representations of states,
without storing the bases in which states and operators are represented. Quantum.jl 
is designed to enable transformations of those representations to extend
naturally to the abstract quantum objects that underlie them (bases of bras,
kets, outer products, etc.), and conversely use abstract quantum objects to aid
in the manipulation and analysis of representations.

Many other libraries that offer a similarly abstract feature set 
(e.g. Quantum Mathematica) are either:

a) poorly optimized for heavy numeric work performed on matrix representations or 

b) don't preserve abstract structure when representation optimization is present. 

The goal of Quantum.jl is to be as optimized for numeric work as a library like
QuTiP, while still having the abstract feature set that other libraries have shown 
to be useful. Luckily, Julia is a language that quite naturally supports high levels 
of abstraction with little sacrifice in performance.

Quantum.jl currently only supports work in discrete spaces, but provides a foundation
that could be extended to support work in continuous spaces in the future. 

Additionally, it has been suggested that Quantum.jl implement an interface for other 
popular quantum libraries like QuTiP in order to more quickly offer up a host of advanced 
optimized features, such as continuous space support, visualization systems, 
and common algorithms used in quantum mechanics.

By the time Quantum.jl is released (hopefully within the next few weeks), it will 
have decent documentation and be registered as an actual Julia package for ease of
use.

Features
==========

###Current:
	- Basis, State, ScalarExpr, DiracVector, and DiracMatrix type implementations
	- Tensor product structure for Dirac objects
	- Application of arbitrary selection rules to extract subspaces from most types of 
	  Quantum.jl collection objects. Related methods:
	  	- filter
	  	- map/map!
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
	- Refactor for optimization; expected finish date is early July
	- Documentation
	- registration of Quantum.jl as Julia package

###Post-Release Goals/Ideas:
	- Interface with other similar libraries (eg. QuTiP)
	- Continuous space support
	- HDF5 support using HDF5.jl (https://github.com/timholy/HDF5.jl)
	- Visualization system (probably via interfacing to more mature libraries)
	- incorporate units using SIUnits.jl (https://github.com/Keno/SIUnits.jl) 
