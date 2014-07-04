#necessary for SparseMatrixCSC{Any} to function properly
zero(::Type{Any}) = 0
zero{D<:Dirac}(::Type{D}) = 0
zero(::Type{ScalarExpr}) = ScalarExpr(:(0+0))

samebasis(a::Dirac, b::Dirac)= bsym(a)==bsym(b)
samebasis(a::Symbol, b::Dirac)= a==bsym(b)
samebasis(a::Dirac, b::Symbol)= samebasis(b, a)

#additive identities
for t=(:DiracMatrix, :DiracVector, :State, :OuterProduct)
	@eval +(n::Number, d::($t)) = error("cannot add number and $(string($t))")
	@eval +(d::($t), n::Number) = +(n,d)
	@eval -(n::Number, d::($t)) = +(n,d)
	@eval -(d::($t), n::Number) = +(n,d)
end

*(a::Dirac, b::Dirac) = kron(a,b)
*(a, b::Dirac) = kron(a,b)
*(a::Dirac, b) = kron(a,b)


# kron(a::AbstractScalar, b::AbstractScalar) = a*b
# kron(a::Dirac, b::DiracCoeff) = a*b
# kron(a::DiracCoeff, b::Dirac) = b*a

# #arithmetic functrions for generic arrays
# for op in (:+, :-, :*, :.+, :.-, :.*, :./, :kron)
# 	@eval ($op){T}(a::Array{T, 1}, d::DiracVector{Ket}) = ($op)(DiracVector(a, d.basis), d)
# 	@eval ($op){T}(d::DiracVector{Ket}, a::Array{T, 1}) = ($op)(d, DiracVector(a, d.basis))
# 	@eval ($op){T}(a::Array{T, 1}, d::DiracVector{Bra}) = ($op)(DiracVector(a, d.basis'), d)
# 	@eval ($op){T}(d::DiracVector{Bra}, a::Array{T, 1}) = ($op)(d, DiracVector(a, d.basis'))
# 	@eval begin
# 			function ($op)(a::AbstractArray, d::DiracVector)
# 				if size(a,1)>1
# 					return ($op)(DiracMatrix(a, d.basis), d)
# 				elseif size(a,1)==1
# 					return ($op)(DiracVector(a, kind(d)==Ket ? d.basis' : d.basis), d)
# 				elseif size(a,2)==1
# 					return ($op)(DiracVector(a, kind(d)==Ket ? d.basis : d.basis'), d)
# 				else
# 					error("DimensionMismatch")
# 				end
# 			end
# 	end	
# 	@eval begin
# 			function ($op)(d::DiracVector, a::AbstractArray)
# 				if size(a,1)>1
# 					return ($op)(d, DiracMatrix(a, d.basis))
# 				elseif size(a,1)==1
# 					return ($op)(d, DiracVector(a, kind(d)==Ket ? d.basis' : d.basis))
# 				elseif size(a,2)==1
# 					return ($op)(d, DiracVector(a, kind(d)==Ket ? d.basis : d.basis'))
# 				else
# 					error("DimensionMismatch")
# 				end
# 			end
# 	end	
# 	@eval ($op){T}(a::Array{T, 1}, m::DiracMatrix) = ($op)(DiracVector(a, m.rowb), m)
# 	@eval ($op){T}(m::DiracMatrix, a::Array{T, 1}) = ($op)(m, DiracVector(a, m.rowb))
# 	@eval begin
# 			function ($op)(a::AbstractArray, m::DiracMatrix)
# 				if size(a,1)==size(m,1) && size(a,2)==size(m,2)
# 					return ($op)(DiracMatrix(a, m.rowb, m.colb), m)
# 				elseif size(a,1)==1
# 					return ($op)(DiracVector(a, m.colb), m)
# 				elseif size(a,2)==1
# 					return ($op)(DiracVector(a, m.rowb), m)
# 				else
# 					error("DimensionMismatch")
# 				end
# 			end
# 	end	
# 	@eval begin
# 			function ($op)(m::DiracMatrix, a::AbstractArray)
# 				if size(a,1)==size(m,1) && size(a,2)==size(m,2)
# 					return ($op)(m, DiracMatrix(a, m.rowb, m.colb))
# 				elseif size(a,1)==1
# 					return ($op)(m, DiracVector(a, m.colb))
# 				elseif size(a,2)==1
# 					return ($op)(m, DiracVector(a, m.rowb))
# 				else
# 					error("DimensionMismatch")
# 				end
# 			end
# 	end	
# end
