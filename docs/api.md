*CURRENTLY A WORK IN PROGRESS*
Quantum.jl API
===

	label(s::State)
		returns the label or vector of labels identifying `s`
	labeldelta(a::State, b::State)
		returns 1 if `label(a)==label(b)`, or 0 otherwise
	inner{B<:Bra, K<:Ket}(a::State{B}, b::State{K})
		computes the inner product of `a` and `b`

	bsym(d::Dirac)
		returns the basis symbol or vector of basis symbols 
		associated with the Dirac object
	samebasis(a::Dirac, b::Dirac)
		returns true if `a` and `b` have the same basis symbol(s).
		If `a` is a state and `b` is a basis, this function does *not* 
		imply that `b` explicitly contains `a`; to check for that kind
		of relationship, use `in(a,b)`
	ctranspose(d::Dirac)
		returns the dual of a state
	isdual(a::Dirac, b::Dirac)
		checks to see whether `a` is the dual of `b`.
		Much more efficient than `a'==b` for many
		Dirac objects.
	kron(a::Dirac, b::Dirac)
		computes the kronecker product of `a` and `b`
	*(a::Dirac, b::Dirac)
		vector multiplication between `a` and `b`
