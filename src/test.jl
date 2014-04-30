include("DiracTools.jl")
using DiracTools

s = StateRep(:psi, [1:5], Basis("b", [1:5]))
normalize!(s)
op = s*s'
st = s*s

normalize!(s)
print("")