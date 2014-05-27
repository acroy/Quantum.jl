type DiracMatrix{C<:DiracCoeff} <: Dirac
	coeffs::Matrix{C}
	rowbasis::AbstractBasis{Ket}
	colbasis::AbstractBasis{Bra}
	function DiracMatrix(coeffs, rowbasis, colbasis)
		if size(coeffs)==(length(rowbasis), length(colbasis))
			new(coeffs, rowbasis, colbasis)
		else
			throw(DimensionMismatch("Bases do not match coefficient matrix"))
		end
	end
end

DiracMatrix{C<:DiracCoeff}(coeffs::Matrix{C}, rowbasis::AbstractBasis{Ket}, colbasis::AbstractBasis{Bra}) = DiracMatrix{C}(coeffs, rowbasis, colbasis) 
DiracMatrix{C<:DiracCoeff}(coeffs::Matrix{C}, b::AbstractBasis{Ket}) = DiracMatrix(coeffs, b, b') 
DiracMatrix{C<:DiracCoeff}(coeffs::Matrix{C}, b::AbstractBasis{Bra}) = DiracMatrix(coeffs, b', b) 
DiracMatrix(coeffs::Matrix, basisargs...) = DiracMatrix(convert(Matrix{DiracCoeff}, coeffs), basisargs...)

function DiracMatrix(fcoeff::Function, flabel::Function, b::AbstractBasis)
	coeffs = convert(Array{DiracCoeff}, zeros(length(b), length(b)))
	for i=1:length(b)
		for j=1:length(b)
			println("getting index of $(State(flabel(b[i])))")
			ii = get(b, State(flabel(b[i]), label(b)), -1)
			if ii != -1
				println("doin it: $(b[j]')*$(b[i])")
				coeffs[ii,j] = fcoeff(b[i])*(b[j]'*b[i])
			end
		end
	end
	return DiracMatrix(coeffs, b)
end
