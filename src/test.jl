include("Quantum.jl")
using Quantum
# b = Basis("b", [1:3])
# a = Basis("a", [1:10])
# c = Basis("c", [1:10])
# ab = a[1]'*b[1]
# ba=a[1]*b[1]'
# ex = ab/ab^4
# d = DiracVector([ab, ex, 3.0, 4+im, 5:10], a)
# qd = qeval((b,k)->1, d)
# op=qd*qd'
basis1 = Basis("1", [1:10])
basis2 = Basis("2", [1:10])
rbasis = basis1*basis2

raise1 = DiracMatrix(s->sqrt(label(s[1])+1), s->State(label(s[1])+1, label(basis1))*s[2], rbasis)
raise2 = DiracMatrix(s->sqrt(label(s[2])+1), s->s[1]*State(label(s[2])+1, label(basis2)), rbasis)

raiseop1(orig1, orig2, theta) = (cos(theta)*orig1) - (sin(theta)*orig2)
raiseop2(orig1, orig2, theta) = (sin(theta)*orig1) + (cos(theta)*orig2)

function rotate(nus, raiseops)
	coeff = 1/sqrt(prod(map(x->factorial(x), nus)))
	s = prod([State(0, "$i") for i=1:length(nus)])
	for i=1:length(nus)
		s = raiseops[i]^nus[i] * s
	end
	return coeff*s
end


print("")