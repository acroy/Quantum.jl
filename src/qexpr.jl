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

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");

*(o::OuterProduct, s::State{Bra}) = OuterProduct(o.ket, o.bra*s)
*(s::State{Ket}, o::OuterProduct) = OuterProduct(s*o.ket, o.bra)

#The following definitions assume that StateRep coefficients can be InnerProducts

# *(i::InnerProduct, o::OuterProduct)= StateRep([], [i.bra*(i.ket*o.ket)], statetobasis(o.bra))
# *(o::OuterProduct, i::InnerProduct) = StateRep([], [(o.bra*i.bra)*i.ket], statetobasis(o.ket))

# *(o::OuterProduct, s::State{Ket}) = StateRep([], [o.bra*s], statetobasis(o.ket))
# *(s::State{Bra}, o::OuterProduct) = StateRep([], [s*o.ket], statetobasis(o.bra))

# *(s::State, i::InnerProduct) = StateRep([], [i], statetobasis(s))
# *(i::InnerProduct, s::State) = *(s,i)



