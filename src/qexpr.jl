export InnerProduct
export OuterProduct
export StateProduct

immutable InnerProduct <: Quantum
	bra::State{Bra}
	ket::State{Ket}
end

conj(s::InnerProduct) = InnerProduct(s.ket', s.bra')
ctranspose(s::InnerProduct) = conj(s)

immutable OuterProduct <: Quantum
	ket::State{Ket}
	bra::State{Bra}
end 


#Assuming staterep can hold InnerProducts
#*(i::InnerProduct, o::OuterProduct) = StateRep([], [i.bra*(i.ket*o.ket)], Basis("", [o.Bra]))
#*(o::OuterProduct, i::InnerProduct) = 




show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");