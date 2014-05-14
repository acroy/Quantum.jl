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

immutable StateProduct <: Quantum
	coeff::Complex{Float64}
	inner::Vector{InnerProduct}
	outer::Vector{OuterProduct}
	StateProduct(n::Number, i, o) = new(convert(Complex{Float64}, n), i, o)
end

copy(sp::StateProduct, coeff=sp.coeff) = StateProduct(coeff, sp.inner, sp.outer) 


*(n::Number, i::InnerProduct) = StateProduct(n, [i], [])
*(i::InnerProduct, n::Number) = *(n,i)
*(n::Number, o::OuterProduct) = StateProduct(n, [], [o])
*(o::OuterProduct, n::Number) = *(n,o)

*(i::InnerProduct, sp::StateProduct) = StateProduct(sp.coeff, vcat(sp.inner, i), sp.outer) 
*(sp::StateProduct, i::InnerProduct) = *(i, sp)
*(o::OuterProduct, sp::StateProduct) = StateProduct(sp.coeff, sp.inner, vcat(sp.outer, o)) 
*(sp::StateProduct, o::OuterProduct) = *(o, sp)


*(n::Number, sp::StateProduct) = copy(sp, n*sp.coeff)
*(sp::StateProduct, n::Number) = *(n, sp)
/(sp::StateProduct, n::Number) = copy(sp, sp.coeff/n)


show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");
function show(io::IO, sp::StateProduct)
	print(io, sp.coeff)
	for i in sp.inner
		print(io, " * $(repr(i))")
	end
	for o in sp.outer
		print(io, " * $(repr(o))")
	end
end