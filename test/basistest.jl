module BasisTest
using Base.Test 
using Quantum

xb = Basis([1:5], :X)
xbplus = Basis([2:6], :X)
xb2 = xb*xb
@test separate(xb2)==[xb,xb]
@test separate(xb)==[xb]
@test tobasis([State(i, :X) for i=1:5])==xb
@test tobasis(TensorState{Ket}[xb[j]*xb[i] for i=1:5, j=1:5])==xb2
@test isdual(xb, xb')
@test Basis([1,3,5], :X) == filter(s->label(s)%2!=0, xb)
@test tobasis(xb[1]*xb[1]) == filter(s->s==xb[1]*xb[1], xb2)
@test map(s->State(label(s)+1, :X), xb) == Basis([2:6], :X)
@test map(s->State(label(s[1])+1, :X)*s[2], xb2) == xbplus*xb
@test setdiff(xbplus, xb) == State{Ket}[xbplus[5]]
@test xb+xb == xb
@test xbplus[5]+xb == Basis([6, 1:5],:X) 
@test xb+xbplus[5] == Basis([1:6],:X) 

end