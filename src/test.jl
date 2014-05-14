include("Quantum.jl")
using QuantumJL

x = State(:x)
y = State(:y)

sr = StateRep(:s, [1:2], Basis("b",[x,y]))

println("")