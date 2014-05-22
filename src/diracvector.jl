#####################################
#DiracVector#########################
#####################################
type DiracVector{C<:DiracCoeff,K<:BraKet} <: Dirac
	coeffs::Array{C} 
	basis::AbstractBasis{K}
	function DiracVector(coeffs, basis)
		if K==Ket
			if size(coeffs)==(length(basis),)
				new(coeffs, basis)
			elseif size(coeffs)==(length(basis),1)
				new(vec(coeffs), basis)
			else
				error("Dimensions of coefficient array does not match type $K")
			end
		else
			if size(coeffs)==(1,length(basis))
				new(coeffs, basis)
			else
				error("Dimensions of coefficient array does not match type $K")
			end
		end
	end
end

DiracVector{C<:DiracCoeff,K<:BraKet}(coeffs::Array{C}, basis::AbstractBasis{K}) = DiracVector{C,K}(coeffs, basis)
DiracVector(coeffs::Array, basis::AbstractBasis) = DiracVector(convert(Array{DiracCoeff}, coeffs), basis)

# #####################################
# #Functions###########################
# #####################################

showcompact(io::IO, d::DiracVector) = print(io, "($(d.coeffs[1])$(d.basis[1]))",[" + ($(d.coeffs[i])$(d.basis[i]))" for i=2:length(d)]...)
function show(io::IO, d::DiracVector)
	println("$(typeof(d)):")
	table = cell(length(d), 2)	
	if length(d.coeffs)>=52
		for i=1:25
			table[i,1]= d.coeffs[i]
			table[i,2]= d.basis[i]
		end
		table[26:(length(d.coeffs)-25),:] = 0 # prevents access to undefined reference
		for i=(length(d.coeffs)-25):length(d.coeffs)
			table[i,1]= d.coeffs[i]
			table[i,2]= d.basis[i]
		end
	else
		for i=1:length(d.coeffs)
			table[i,1]= d.coeffs[i]
			table[i,2]= d.basis[i]
		end
	end
	temp_io = IOBuffer()
	if kind(d)==Ket
		show(temp_io, table)
	else
		show(temp_io, [transpose(table[:,2]), transpose(table[:,1])])
	end
	io_str = takebuf_string(temp_io)
	io_str = io_str[searchindex(io_str, "\n")+1:end]
	print(io_str)
end

#imported############################
ctranspose(d::DiracVector) = DiracVector(d.coeffs', d.basis')
getindex(d::DiracVector, x) = d.coeffs[x]
length(d::DiracVector) = length(d.coeffs)
size(d::DiracVector, x=nothing) = size(d.coeffs, x)
kind(d::DiracVector) = kind(d.basis)
setindex!(d::DiracVector, y, x) = setindex!(d.coeffs, y, x)
endof(d::DiracVector) = length(d)

for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracVector, b::DiracVector) = DiracVector(($op)(a.coeffs,b.coeffs), a.basis)
	@eval ($op)(n, d::DiracVector) = DiracVector(($op)(n,d.coeffs), d.basis)
	@eval ($op)(d::DiracVector, n) = DiracVector(($op)(d.coeffs,n), d.basis)
end

*{N1<:Number, N2<:Number}(a::DiracVector{N1, Bra}, b::DiracVector{N2, Ket}) = (a.coeffs*b.coeffs)[1]
*{A<:DiracCoeff, B<:DiracCoeff}(a::DiracVector{A, Bra}, b::DiracVector{B, Ket}) = length(a)==length(b) ? reduce(+, [a[i]*b[i] for i=1:length(a)]) : error("DimensionMismatch")

