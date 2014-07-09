*CURRENTLY A WORK IN PROGRESS*
Quantum.jl API
===

	label(s::State)
		returns the label of a state
	bsym(s::State)
		returns the basis symbol associated with the state
	ctranspose(s::State)
		returns the dual of a state
	isdual(a::State, b::State)
		checks to see whether `a` is the dual of `b`
	labeldelta(a::State, b::State)
		returns 1 if label(a)==label(b), or 0 otherwise
	inner{B<:Bra, K<:Ket}(a::State{B}, b::State{K})
		computes the inner product of `a` and `b`
	kron(a::State, b::State)
		computes the kronecker product of `a` and `b`
	*(a::State, b::State)
		vector multiplication between `a` and `b`
