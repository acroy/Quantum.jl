export InnerProduct
export OuterProduct
export QExpr

type InnerProduct <: Dirac
	bra::State{Bra}
	ket::State{Ket}
end

conj(s::InnerProduct) = InnerProduct(s.ket', s.bra')
ctranspose(s::InnerProduct) = conj(s)

type OuterProduct <: Dirac
	ket::State{Ket}
	bra::State{Bra}
end 

type QExpr <: Dirac
	coeff::Complex{Float64}
	ops::Vector{Function}
	args::Vector{Dirac}
	QExpr(coeff, ops, args) = coeff==0 ? 0 : new(coeff, ops, args)
end

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");
*(o::OuterProduct, s::State{Bra}) = OuterProduct(o.ket, o.bra*s)
*(s::State{Ket}, o::OuterProduct) = OuterProduct(s*o.ket, o.bra)
*(n::Number, q::QExpr) = QExpr(n*q.coeff, q.ops, q.args)
*(q::QExpr, n::Number) = *(n, q)
*(n::Number, d::Dirac) = QExpr(n, [], [d])
*(d::Dirac, n::Number) = *(n, d)

for op = (:*,:+,:-)
	@eval ($op)(q1::QExpr, q2::QExpr) = QExpr(q1.coeff, vcat(q1.ops, $op), vcat(q1.args, q2))
	@eval ($op)(d1::Dirac, d2::Dirac) = QExpr(1, [$op], [d1, d2])
	@eval ($op)(d::Dirac, q::QExpr) = QExpr(q.coeff, vcat($op, q.ops), vcat(d, q.args))
	@eval ($op)(q::QExpr, d::Dirac) = QExpr(q.coeff, vcat(q.ops, $op), vcat(q.args, d))
end

function show(io::IO, qex::QExpr) 
	if qex.coeff!=1+0im
		print(io, "$(qex.coeff) * (");
	end
	print(io, " $(qex.args[1]) ")
	for i=2:length(qex.args)
		if !isequal(qex.ops[i-1], *)
			print(io, " $(qex.ops[i-1]) ")
		end
		if typeof(qex.args[i])==QExpr
			print(io, " ($(qex.args[i])) ")
		else
			print(io, "$(qex.args[i])")
		end
	end
	if qex.coeff!=1+0im
		print(io, ")")
	end
end

#The following definitions assume that StateRep coefficients can be InnerProducts

# *(o1::OuterProduct, o2::OuterProduct) = StateRep([], [(o2.bra*o1.bra)*o2.ket], statetobasis(o1.ket))
# *(i::InnerProduct, o::OuterProduct) = StateRep([], [i.bra*(i.ket*o.ket)], statetobasis(o.bra))
# *(o::OuterProduct, i::InnerProduct) = StateRep([], [(o.bra*i.bra)*i.ket], statetobasis(o.ket))

# *(o::OuterProduct, s::State{Ket}) = StateRep([], [o.bra*s], statetobasis(o.ket))
# *(s::State{Bra}, o::OuterProduct) = StateRep([], [s*o.ket], statetobasis(o.bra))

# *(s::State, i::InnerProduct) = StateRep([], [i], statetobasis(s))
# *(i::InnerProduct, s::State) = *(s,i)