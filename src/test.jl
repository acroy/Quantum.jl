include("Quantum.jl")
using QuantumJL

b1 = Basis(:b1, [1:10])
s1 = StateRep(:s1, normalize([1:10]), b1)
b2 = Basis(:b2, ["a", "b", "c"])
s2 = StateRep(:s2, normalize([1:3]), b2)
println("")