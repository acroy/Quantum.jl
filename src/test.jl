include("Quantum.jl")
using QuantumJL


s1 = State(1)
s2 = State(2)
Qcoeff = Union(Complex{Float64}, InnerProduct)
a = Array(Qcoeff, 10)
a[1] = s1'*s2
a[2:end] = 1.0 + im
b = Basis(:b, [1:1000])

println("")