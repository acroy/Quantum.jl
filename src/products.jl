typealias DiracCoeff Union(Number, AbstractScalar)

immutable OuterProduct <: Dirac
	ket::AbstractState{Ket}
	bra::AbstractState{Bra}
end

basislabel(o::OuterProduct) = [basislabel(o.ket), basislabel(o.bra)]
label(o::OuterProduct) = [label(o.ket), label(o.bra)]
*(o::OuterProduct, s::AbstractState{Ket}) = DiracVector([(o.bra*s)], tobasis(o.ket))
*(o::OuterProduct, s::AbstractState{Bra}) = OuterProduct(o.ket, o.bra*s)
*(s::AbstractState{Bra}, o::OuterProduct) = DiracVector([(s*o.ket)], tobasis(o.bra))
*(s::AbstractState{Ket}, o::OuterProduct) = OuterProduct(s*o.ket, o.bra)
*(d::DiracCoeff, o::OuterProduct) = DiracMatrix([d]', tobasis(o.ket), tobasis(o.bra))
*(o::OuterProduct, d::DiracCoeff) = *(d,o)
-(o::OuterProduct) = -1*o

function +(a::OuterProduct, b::OuterProduct)
	if a==b
		return DiracMatrix(2.0, tobasis(a.ket), tobasis(b.bra))
	elseif basislabel(a)==basislabel(b)
		rowb = tobasis([a.ket, b.ket])
		colb = tobasis([a.bra, b.bra])
		res = DiracMatrix(zeros(2,2)[1:length(rowb), 1:length(colb)], rowb, colb)
		res[getpos(res, a)...] = 1.0
		res[getpos(res, b)...] = 1.0
		return res
	end
end

-(a::OuterProduct, b::OuterProduct) = a+(-b)

ctranspose(o::OuterProduct) = OuterProduct(o.bra', o.ket')

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");

##########################################################################

immutable InnerProduct <: AbstractScalar
	bra::AbstractState{Bra}
	ket::AbstractState{Ket}
end

basislabel(i::InnerProduct) = [basislabel(i.bra), basislabel(i.ket)]
label(i::InnerProduct) = [label(i.bra), label(o.ket)]
conj(i::InnerProduct) = InnerProduct(i.ket', i.bra')
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");
##########################################################################

##########################################################################

immutable ScalarExpr <: AbstractScalar
	ex::Expr
end

qexpr(s::ScalarExpr) = s.ex
qexpr(x) = x

function ^(d::DiracCoeff, n::Integer)
	if n==1
		return d
	elseif n==0
		return 1
	else
		ScalarExpr(:($(qexpr(d))^$(qexpr(n))))
	end
end

function ^(a::DiracCoeff, b::DiracCoeff)
	if b==1
		return a
	elseif b==0
		return 1
	else
		ScalarExpr(:($(qexpr(a))^$(qexpr(b))))
	end
end

function *(a::DiracCoeff, b::DiracCoeff)
	if a==1
		return b
	elseif b==1
		return a
	elseif a==0 || b==0
		return 0
	else
		ScalarExpr(:($(qexpr(a))*$(qexpr(b))))
	end
end

function +(a::DiracCoeff, b::DiracCoeff)
	if a==0
		return b
	elseif b==0
		return a
	else
		ScalarExpr(:($(qexpr(a))+$(qexpr(b))))
	end
end

function -(a::DiracCoeff, b::DiracCoeff)
	if a==0
		return b
	elseif b==0
		return a
	elseif a==b
		return 0
	else
		ScalarExpr(:($(qexpr(a))-$(qexpr(b))))
	end
end

-(d::DiracCoeff) = ScalarExpr(:(-$(qexpr(d))))
-(s::ScalarExpr) = length(s.ex.args)==2 && s.ex.args[1]==:- ? ScalarExpr(s.ex.args[2]) :  ScalarExpr(:(-$(qexpr(s))))

abs(i::InnerProduct) = ScalarExpr(:(abs($i)))
abs(s::ScalarExpr) = length(s.ex.args)==2 && s.ex.args[1]==:abs ? s :  ScalarExpr(:(abs($(qexpr(s)))))

exp(a::AbstractScalar) = ScalarExpr(:(exp($(qexpr(a)))))

function /(a::DiracCoeff, b::DiracCoeff)
	if a==0
		return b
	elseif b==0
		error("Zero in denominator")
	elseif a==b
		return 1
	else
		ScalarExpr(:($(qexpr(a))/$(qexpr(b))))
	end
end

conj(s::ScalarExpr)	= length(s.ex.args)==2 && s.ex.args[1]==:conj ? ScalarExpr(s.ex.args[2]) :  ScalarExpr(:(conj($(qexpr(s)))))

qeval(f::Function, s::ScalarExpr) = eval(qreduce(f, s.ex))
qeval(f::Function, i::InnerProduct) = eval(f(i.bra, i.ket))
qeval(f::Function, n::Number) = n

function qreduce(f::Function, ex::Expr)
	ex = copy(ex)
	for i=1:length(ex.args)
		if typeof(ex.args[i])==InnerProduct
			ex.args[i] = f(ex.args[i].bra, ex.args[i].ket)
		elseif typeof(ex.args[i])==Expr
			ex.args[i] = qreduce(f, ex.args[i])
		end
	end
	return ex
end

for op=(:*,:-,:+,:/,:^) #define for elementwise operators
	elop = symbol(string(:.) * string(op))
	@eval ($elop)(a::DiracCoeff, b::DiracCoeff) = ($op)(a,b)
	@eval ($elop)(a::AbstractScalar, arr::Array) = broadcast(($op), [a], arr)
	@eval ($elop)(arr::Array, a::AbstractScalar) = broadcast(($op), arr, [a])
end

*(a::AbstractScalar, arr::Array) = a.*arr
*(arr::Array, a::AbstractScalar) = arr.*a
/(arr::Array, a::AbstractScalar) = arr./a