#####################################
#DiracMatrix#########################
#####################################

type DiracMatrix{T} <: Dirac
	coeffs::Matrix{T}
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

DiracMatrix{T}(coeffs::Matrix{T}, rowbasis::AbstractBasis{Ket}, colbasis::AbstractBasis{Bra}) = DiracMatrix{T}(coeffs, rowbasis, colbasis) 
DiracMatrix(coeffs::Matrix, b::AbstractBasis{Ket}) = DiracMatrix(coeffs, b, b') 
DiracMatrix(coeffs::Matrix, b::AbstractBasis{Bra}) = DiracMatrix(coeffs, b', b) 

function DiracMatrix(fcoeff::Function, flabel::Function, b::AbstractBasis)
	coeffs = convert(Array{Any}, zeros(length(b), length(b)))
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

#####################################
#Show Functions######################
#####################################
function showcompact(io::IO, op::DiracMatrix)
	if length(op.coeffs)==0
		print(io, "$(typeof(op))[]")
	else
		tempio = IOBuffer()
		print(tempio, [" + ($(op.coeffs[i,j])$(op.rowbasis[i])$(op.colbasis[j]))" for i=1:length(op.rowbasis), j=1:length(op.colbasis)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end
function show(io::IO, op::DiracMatrix)
	println("$(typeof(op)):")
	table = cell(length(op.rowbasis)+1, length(op.colbasis)+1)	
	for i = 1:length(op.rowbasis)
		table[i+1,1] = op.rowbasis[i]
	end
	for j = 1:length(op.colbasis)
		table[1,j+1] = op.colbasis[j]
	end
	table[1,1] = 0
	table[2:end, 2:end] = op.coeffs	
	temp_io = IOBuffer()
	show(temp_io, table)
	io_str = takebuf_string(temp_io)
	print(io, io_str[searchindex(io_str, "\n")+3:end])
end