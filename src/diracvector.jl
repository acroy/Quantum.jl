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

getpos(d::DiracVector, s::AbstractState) = get(d.basis, s)
function get(d::DiracVector, s::AbstractState, notfound)
	try
		return d[get(d.basis, s)]
	catch
		return notfound
	end
end

get(d::DiracVector, s::AbstractState) = d[get(d.basis, s)]
get(d::DiracVector, label) = get(d, typeof(d.basis)<:Basis ? State(label) : TensorState(label))

for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracVector, b::DiracVector) = DiracVector(($op)(a.coeffs,b.coeffs), a.basis)
	@eval ($op)(n, d::DiracVector) = DiracVector(($op)(n,d.coeffs), d.basis)
	@eval ($op)(d::DiracVector, n) = DiracVector(($op)(d.coeffs,n), d.basis)
end

*(c::DiracCoeff, d::DiracVector) = c.*d
*(d::DiracVector, c::DiracCoeff) = c*d

*{C<:DiracCoeff}(s::AbstractState{Bra}, d::DiracVector{C, Ket}) = get(d, s', 0)
*{C<:DiracCoeff}(d::DiracVector{C, Bra}, s::AbstractState{Ket}) = get(d, s', 0)

*{N1<:Number, N2<:Number}(a::DiracVector{N1, Bra}, b::DiracVector{N2, Ket}) = (a.coeffs*b.coeffs)[1]
*{A<:DiracCoeff, B<:DiracCoeff}(a::DiracVector{A, Bra}, b::DiracVector{B, Ket}) = length(a)==length(b) ? reduce(+, [a[i]*b[i] for i=1:length(a)]) : throw(DimensionMismatch(""))

*{C<:DiracCoeff}(c::C, s::AbstractState) = DiracVector([c], statetobasis(s))
*{C<:DiracCoeff}(s::AbstractState, c::C) = *(c,s)

function +{C<:DiracCoeff,K<:BraKet}(d::DiracVector{C,K}, s::AbstractState{K})
	if in(s, d.states)
		d[getpos(d,s)] = 1+get(d, s)
	else
		DiracVector(vcat(d.coeffs, 1), d.basis+statetobasis(s))
	end
end

function +{C<:DiracCoeff,K<:BraKet}(s::AbstractState{K}, d::DiracVector{C,K})
	if in(s, d.states)
		d[getpos(d,s)] = 1+get(d, s)
	else
		DiracVector(vcat(1, d.coeffs), statetobasis(s)+d.basis)
	end
end

function +{C1,C2,K<:BraKet}(a::DiracVector{C1,K}, b::DiracVector{C2,K})
	if a.basis==b.basis
		return DiracVector(a.coeffs+b.coeffs, a.basis)
	else
		res = copy(a)
		bdiff = setdiff(b.basis, a.basis) 
		compl = setdiff(b.basis, bdiff)
		for i in compl
			res[getpos(res, i)] = get(res, i) + get(b, i)
		end
		res = DiracVector(vcat(res.coeffs, [get(b, bdiff[i]) for i=1:length(bdiff)]), res.basis+Basis(bdiff))
		return res
	end
end

-{C<:DiracCoeff,K<:BraKet}(d::DiracVector{C,K}, s::AbstractState{K}) = d+(-s)
-{C<:DiracCoeff,K<:BraKet}(s::AbstractState{K}, d::DiracVector{C,K}) = s+(-d)
-{C1,C2,K<:BraKet}(a::DiracVector{C1,K}, b::DiracVector{C2,K}) = a+(-b)

+{K<:BraKet}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1, 1], Basis([a,b])) 
-{K<:BraKet}(a::AbstractState{K}, b::AbstractState{K}) = DiracVector([1,-1], Basis([a,b])) 
-(s::AbstractState) = -1*s
-(d::DiracVector) = DiracVector(-1*d.coeffs, d.basis)