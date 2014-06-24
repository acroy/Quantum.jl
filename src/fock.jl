#This file contains convenience functions and structures
#as they relate to handling Fock spaces

fbasis(levels::Integer) = basis([0:levels-1],:F)
fbasis(levels::Integer, systems::Integer) = tensor([fbasis(levels) for i=1:systems]...)
function fvec(levels::Integer, init::Integer=0)
	d=DiracVector(spzeros(levels,1), fbasis(levels))
	d[init+1] = 1.0
	return d
end

#TODO:
# Operator construction:
# lowerf()
# raisef()
# sigma/spin operators

