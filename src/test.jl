include("Quantum.jl")
using Quantum
b = Basis("b", [1:10])
a = Basis("a", [1:10])
ab = a[1]'*b[1]
ba=a[1]*b[1]'
ex = ab/ab^4
d = DiracVector([ab, ex, 3.0, 4+im], Basis(a[1:4]))

type Foo
	arr::Array
end
print("")