#This file contains convenience functions and structures
#as they relate to handling Fock spaces

fbasis(levels::Integer) = basis([0:levels-1],:F)
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

fcreate(levels::Integer) = DiracMatrix(raiseload!(zeros(levels, levels)), fbasis(levels))
fdestroy(levels::Integer) = DiracMatrix(lowerload!(zeros(levels, levels)), fbasis(levels))
fnum(levels::Integer) = DiracMatrix(numload!(zeros(levels, levels)), fbasis(levels))

#TODO
# sigma/spin operators

