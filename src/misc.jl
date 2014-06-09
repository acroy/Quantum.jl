
samebasis(a::Dirac, b::Dirac)= basislabel(a)=="?" || basislabel(b)=="?" ? false : basislabel(a)==basislabel(b)
samebasis(a::String, b::Dirac)= a=="?" ? false : a==basislabel(b)
samebasis(a::Dirac, b::String)= samebasis(b, a)

#additive identities
for t=(:DiracMatrix, :DiracVector, :AbstractState, :OuterProduct)
	@eval +(n::Number, d::($t)) = n==0 ? d : error("cannot add number and $(string($t))")
	@eval +(d::($t), n::Number) = +(n,d)
	@eval -(n::Number, d::($t)) = +(n,d)
	@eval -(d::($t), n::Number) = +(n,d)
end

kron(a::DiracMatrix, b::DiracMatrix) = DiracMatrix(kron(a.coeffs, b.coeffs), tensor(a.rowbasis, b.rowbasis), tensor(a.colbasis, b.colbasis)) 
kron(op::DiracMatrix, d::DiracVector{Ket}) = DiracMatrix(kron(op.coeffs, d.coeffs), tensor(op.rowbasis, d.basis), op.colbasis)
kron(op::DiracMatrix, d::DiracVector{Bra}) = DiracMatrix(kron(op.coeffs, d.coeffs), op.rowbasis, tensor(op.colbasis, d.basis))
kron(d::DiracVector{Ket}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), tensor(d.basis, op.rowbasis), op.colbasis)
kron(d::DiracVector{Bra}, op::DiracMatrix) = DiracMatrix(kron(d.coeffs, op.coeffs), op.rowbasis, tensor(d.basis, op.colbasis))
kron{K}(a::DiracVector{K}, b::DiracVector{K}) = DiracVector(kron(a.coeffs, b.coeffs), tensor(a.basis, b.basis))
kron(a::DiracVector{Ket}, b::DiracVector{Bra}) = DiracMatrix(kron(a.coeffs, b.coeffs), a.basis, b.basis)
kron(a::DiracVector{Bra}, b::DiracVector{Ket}) = kron(b,a)
kron(a::AbstractScalar, b::AbstractScalar) = a*b
kron(a::Dirac, b::DiracCoeff) = a*b
kron(a::DiracCoeff, b::Dirac) = b*a

for op in (:+, :-, :*, :.+, :.-, :.*, :./, :kron)
	@eval ($op){T}(a::Array{T, 1}, d::DiracVector{Ket}) = ($op)(DiracVector(a, d.basis), d)
	@eval ($op){T}(d::DiracVector{Ket}, a::Array{T, 1}) = ($op)(d, DiracVector(a, d.basis))
	@eval ($op){T}(a::Array{T, 1}, d::DiracVector{Bra}) = ($op)(DiracVector(a, d.basis'), d)
	@eval ($op){T}(d::DiracVector{Bra}, a::Array{T, 1}) = ($op)(d, DiracVector(a, d.basis'))
	@eval begin
			function ($op){T}(a::Array{T, 2}, d::DiracVector)
				if size(a,1)>1
					return ($op)(DiracMatrix(a, d.basis), d)
				elseif size(a,1)==1
					return ($op)(DiracVector(a, kind(d)==Ket ? d.basis' : d.basis), d)
				elseif size(a,2)==1
					return ($op)(DiracVector(a, kind(d)==Ket ? d.basis : d.basis'), d)
				else
					error("DimensionMismatch")
				end
			end
	end	
	@eval begin
			function ($op){T}(d::DiracVector, a::Array{T, 2})
				if size(a,1)>1
					return ($op)(d, DiracMatrix(a, d.basis))
				elseif size(a,1)==1
					return ($op)(d, DiracVector(a, kind(d)==Ket ? d.basis' : d.basis))
				elseif size(a,2)==1
					return ($op)(d, DiracVector(a, kind(d)==Ket ? d.basis : d.basis'))
				else
					error("DimensionMismatch")
				end
			end
	end	
	@eval ($op){T}(a::Array{T, 1}, m::DiracMatrix) = ($op)(DiracVector(a, m.rowbasis), m)
	@eval ($op){T}(m::DiracMatrix, a::Array{T, 1}) = ($op)(m, DiracVector(a, m.rowbasis))
	@eval begin
			function ($op){T}(a::Array{T, 2}, m::DiracMatrix)
				if size(a,1)==size(m,1) && size(a,2)==size(m,2)
					return ($op)(DiracMatrix(a, m.rowbasis, m.colbasis), m)
				elseif size(a,1)==1
					return ($op)(DiracVector(a, m.colbasis), m)
				elseif size(a,2)==1
					return ($op)(DiracVector(a, m.rowbasis), m)
				else
					error("DimensionMismatch")
				end
			end
	end	
	@eval begin
			function ($op){T}(m::DiracMatrix, a::Array{T, 2})
				if size(a,1)==size(m,1) && size(a,2)==size(m,2)
					return ($op)(m, DiracMatrix(a, m.rowbasis, m.colbasis))
				elseif size(a,1)==1
					return ($op)(m, DiracVector(a, m.colbasis))
				elseif size(a,2)==1
					return ($op)(m, DiracVector(a, m.rowbasis))
				else
					error("DimensionMismatch")
				end
			end
	end	
end
