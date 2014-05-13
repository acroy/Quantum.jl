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