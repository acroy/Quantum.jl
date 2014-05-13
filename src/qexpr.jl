export InnerProduct
export OuterProduct

immutable InnerProduct <: Number
	bra::State{Bra}
	ket::State{Ket}
end 

show(io::IO, s::InnerProduct) = println("$(repr(s.bra)) $(repr(s.ket)[2:end])");


immutable OuterProduct <: Quantum
	bra::State{Ket}
	ket::State{Bra}
end 

show(io::IO, s::OuterProduct) = println("$(repr(s.bra))$(repr(s.ket))");

*(n::Number, s::State{Ket}) = StateRep([], [n], Basis(label(s),s))
*(n::Number, s::State{Bra}) = StateRep(State([], Bra), [n]', Basis(label(s),s))
*(s::State, n::Number) = *(n,s)
