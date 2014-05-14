export InnerProduct
export OuterProduct

immutable InnerProduct
	bra::State{Bra}
	ket::State{Ket}
end 

show(io::IO, s::InnerProduct) = print(io, "$(repr(s.bra)) $(repr(s.ket)[2:end])");
conj(s::InnerProduct) = InnerProduct(s.ket', s.bra')
ctranspose(s::InnerProduct) = conj(s)

immutable OuterProduct <: Quantum
	bra::State{Ket}
	ket::State{Bra}
end 

show(io::IO, s::OuterProduct) = print(io, "$(repr(s.bra))$(repr(s.ket))");

*(n::Number, s::State{Ket}) = StateRep([], [n], Basis(label(s),s))
*(n::Number, s::State{Bra}) = StateRep(State([], Bra), [n]', Basis(label(s),s))
*(s::State, n::Number) = *(n,s)
