
zero{A<:AbstractScalar}(::Type{A}) = scalar(0)
#necessary for SparseMatrixCSC{Any} to function properly;
#the fact that Any is often inferenced when ScalarExpr
#should be inferenced is a Bad Thing. As a result, in the future
#the below zeros should be removed as soon as the inferencing issue 
#is fixed
zero(::Type{Any}) = 0

orient_error(a,b) = error("Multiplication $(typeof(a))*$(typeof(b)) is undefined. Perhaps you meant to use kron($(typeof(a)), $(typeof(b)))?")

samebasis(a::Dirac, b::Dirac)= bsym(a)==bsym(b)
samebasis(a::Symbol, b::Dirac)= a==bsym(b)
samebasis(a::Dirac, b::Symbol)= samebasis(b, a)

for t=(:DiracMatrix, :DiracVector, :State, :OuterProduct)
	@eval +(n::Number, d::($t)) = error("cannot add number and $(string($t))")
	@eval +(d::($t), n::Number) = +(n,d)
	@eval -(n::Number, d::($t)) = +(n,d)
	@eval -(d::($t), n::Number) = +(n,d)
end

*(a::Dirac, b::Dirac) = orient_error(a,b)
*(a, b::Dirac) = orient_error(a,b)
*(a::Dirac, b) = orient_error(a,b)

kron(a::AbstractScalar, b::AbstractScalar) = a*b
kron(a::Dirac, b::DiracCoeff) = a*b
kron(a::DiracCoeff, b::Dirac) = b*a

convert{A<:Array}(::Type{A}, d::DiracVector) = convert(A, vec(full(d.coeffs)))
convert{A<:Array}(::Type{A}, d::DiracMatrix) = convert(A, full(d.coeffs))
convert{A<:AbstractArray}(::Type{A}, d::Dirac) = convert(A, d.coeffs)

promote_rule{A<:AbstractArray, D<:Dirac}(::Type{A}, ::Type{D}) = A
#arithmetic functrions for generic arrays
for op in (:+, :-, :*, :.+, :.-, :.*, :./, :^, :.^, :kron)
	@eval begin
	($op)(a::AbstractArray, d::Dirac) = ($op)(a, d.coeffs)
	($op)(d::Dirac, a::AbstractArray) = ($op)(d.coeffs, a)
	end
end
