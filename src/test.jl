include("DiracTools.jl")
using DiracTools

s = StateRep(:psi, [1:500000], Basis("b", [1:500000]))
normalize!(s)
print("")