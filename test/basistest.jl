module BasisTest
using Base.Test 
using Quantum

xb = basis(:X,[1:5])
xbplus = basis(:X,[2:6])
xb2 = xb*xb
@test separate(xb2)==[xb,xb]
@test separate(xb)==[xb]
@test basis([Ket{:X,Int}(i) for i=1:5])==xb
@test basis([xb[j]*xb[i] for i=1:5, j=1:5])==xb2
@test isdual(xb, xb')
@test basis(:X, [1,3,5]) == filter(s->label(s)%2!=0, xb)
@test basis(xb[1]*xb[1]) == filter(s->s==xb[1]*xb[1], xb2)
@test map(s->ket(:X,label(s)+1), xb) == basis(:X,[2:6])
@test map(s->ket(:X,label(s[1])+1)*s[2], xb2) == xbplus*xb
@test setdiff(xbplus, xb) == [xbplus[5]]
@test xb+xb == xb
@test xbplus[5]+xb == basis(:X,[6, 1:5]) 
@test xb+xbplus[5] == basis(:X,[1:6]) 

end