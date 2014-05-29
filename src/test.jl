include("Quantum.jl")
using Quantum
b = Basis("b", [1:3])
a = Basis("a", [1:10])
c = Basis("c", [1:10])
ab = a[1]'*b[1]
ba=a[1]*b[1]'
ex = ab/ab^4
d = DiracVector([ab, ex, 3.0, 4+im, 5:10], a)
f = DiracVector([1:3], b)
qd = qeval((b,k)->1, d)
op=qd*qd'
raise = DiracMatrix(x->sqrt(label(x)+1), s->State(label(s)+1, label(a)), a)
print("")