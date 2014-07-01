#This file contains convenience functions and structures
#as they relate to handling Fock spaces

fbasis(levels::Integer) = consbasis([Ket{Int,:F}(i) for i=0:levels-1])
fbasis(levels::Integer, systems::Integer) = tensor([fbasis(levels) for i=1:systems]...)
function fvec(levels::Integer, init::Integer=0)
	d=DiracVector(spzeros(levels,1), fbasis(levels))
	d[init+1] = 1.0
	return d
end

function raiseload!(a::Array)
	for i=2:size(a, 1) 
		a[i, i-1] = sqrt(i-1)
	end
	a
end
function lowerload!(a::Array)
	for i=2:size(a, 2) 
		a[i-1,i] = sqrt(i-1)
	end
	a
end
function numload!(a::Array)
	for i=2:size(a, 1) 
		a[i,i] = i-1	
	end
	a
end

function tensorfop(f::Function, levels::Integer, systems::Integer, act::Integer)
	@assert systems>=act && act>0 "Cannot construct operator to act on system $act out of $systems number of systems"
	id = feye(levels)
	if act==1
		res = f(levels)
	else
		res = feye(levels)
		for i=2:act-1
			res = kron(res,id)
		end
		res = kron(res, f(levels))
	end
	for i=act+1:systems
		res = kron(res,id)
	end
	return res
end

fnum(levels::Integer) = DiracMatrix(numload!(zeros(levels, levels)), fbasis(levels))
feye(levels::Integer) = DiracMatrix(speye(levels, levels), fbasis(levels))
fcreate(levels::Integer) = DiracMatrix(raiseload!(zeros(levels, levels)), fbasis(levels))
fcreate(levels::Integer, systems::Integer, act::Integer) = tensorfop(fcreate, levels, systems, act)
fdestroy(levels::Integer) = DiracMatrix(lowerload!(zeros(levels, levels)), fbasis(levels))
fdestroy(levels::Integer, systems::Integer, act::Integer) = tensorfop(fdestroy, levels, systems, act)

#TODO
# sigma/spin operators

