#####################################
#ScalarExpr##########################
#####################################

immutable ScalarExpr <: AbstractScalar
	ex::Expr
end

isequal(a::ScalarExpr, b::ScalarExpr) = isequal(a.ex, b.ex)
==(a::ScalarExpr, b::ScalarExpr) = ==(a.ex, b.ex)

#####################################
#Misc Functions######################
#####################################
qexpr(s::ScalarExpr) = s.ex
qexpr(x) = x
qeval(f::Function, s::ScalarExpr) = eval(qreduce(f, s.ex))
qeval(f::Function, i::InnerProduct) = eval(f(i.bra, i.ket))
qeval(f::Function, n::Number) = n

function qreduce(f::Function, ex::Expr)
	ex = copy(ex)
	for i=1:length(ex.args)
		if typeof(ex.args[i])<:InnerProduct
			ex.args[i] = f(ex.args[i].bra, ex.args[i].ket)
		elseif typeof(ex.args[i])==Expr
			ex.args[i] = qreduce(f, ex.args[i])
		end
	end
	return ex
end

#####################################
#Arithmetic Operations###############
#####################################
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
		error("Zero in denominator of ScalarExpr")
	elseif a==b
		return 1
	else
		ScalarExpr(:($(qexpr(a))/$(qexpr(b))))
	end
end

conj(s::ScalarExpr)	= length(s.ex.args)==2 && s.ex.args[1]==:conj ? ScalarExpr(s.ex.args[2]) :  ScalarExpr(:(conj($(qexpr(s)))))

for op=(:*,:-,:+,:/,:^) #define for elementwise operators
	elop = symbol(string(:.) * string(op))
	@eval ($elop)(a::DiracCoeff, b::DiracCoeff) = ($op)(a,b)
	@eval ($elop)(a::AbstractScalar, arr::AbstractArray) = broadcast(($op), [a], arr)
	@eval ($elop)(arr::AbstractArray, a::AbstractScalar) = broadcast(($op), arr, [a])
end

*(a::AbstractScalar, arr::AbstractArray) = a.*arr
*(arr::AbstractArray, a::AbstractScalar) = arr.*a
/(arr::AbstractArray, a::AbstractScalar) = arr./a