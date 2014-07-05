#necessary for SparseMatrixCSC{Any} to function properly
zero(::Type{Any}) = 0
zero{D<:Dirac}(::Type{D}) = 0
zero(::Type{ScalarExpr}) = ScalarExpr(:(0+0))

samebasis(a::Dirac, b::Dirac)= bsym(a)==bsym(b)
samebasis(a::Symbol, b::Dirac)= a==bsym(b)
samebasis(a::Dirac, b::Symbol)= samebasis(b, a)

for t=(:DiracMatrix, :DiracVector, :State, :OuterProduct)
	@eval +(n::Number, d::($t)) = error("cannot add number and $(string($t))")
	@eval +(d::($t), n::Number) = +(n,d)
	@eval -(n::Number, d::($t)) = +(n,d)
	@eval -(d::($t), n::Number) = +(n,d)
end

*(a::Dirac, b::Dirac) = kron(a,b)
*(a, b::Dirac) = kron(a,b)
*(a::Dirac, b) = kron(a,b)

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
