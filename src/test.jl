include("Quantum.jl")
using Quantum
xb = Basis([1:5], "X")
xb2 = xb*xb


rx = DiracMatrix(s->sqrt(label(s)+1), s->State(label(s)+1, label(xb)), xb)
r1 = DiracMatrix(s->sqrt(label(s[1])+1), s->State(label(s[1])+1, label(xb))*s[2], xb2)
r2 = DiracMatrix(s->sqrt(label(s[2])+1), s->s[1]*State(label(s[2])+1, label(xb)), xb2)


# raiseop1(orig1, orig2, theta) = (cos(theta)*orig1) - (sin(theta)*orig2)
# raiseop2(orig1, orig2, theta) = (sin(theta)*orig1) + (cos(theta)*orig2)

# function rotate(nus, raiseops)
# 	coeff = 1/sqrt(prod(map(x->factorial(x), nus)))
# 	s = prod([State(0, "$i") for i=1:length(nus)])
# 	for i=1:length(nus)
# 		s = raiseops[i]^nus[i] * s
# 	end
# 	return coeff*s
# end

print("")