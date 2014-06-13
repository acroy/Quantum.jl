module StateTest
using Base.Test 
using Quantum

xb = Basis([1:5], :X)
yb = Basis([1:5], :Y)
xyb = xb*yb

@test 0*xb[1]==0
@test 1*xb[1]==xb[1]
@test xb[1]'*xb[1]==1
@test isdual(xb[1]',xb[1])
@test sum(xb[:])==DiracVector([1,1,1,1,1], xb)
@test xb[1]'*xyb[1]==yb[1]
@test separate(xyb[1])==xyb[1][:]
@test reduce(+,[xb[i]*yb[j]' for i=1:5, j=1:5]) == DiracMatrix(ones(5,5), xb, yb')
@test inner(xyb[1]', xyb[1], 1)==1
@test inner(xyb[1]', xyb[1], 2)==xb[1]'*yb[1]*yb[1]'*xb[1]
@test inner(xyb[1]', xyb[2], 1)==0
@test inner(xyb[1]', xyb[2], 2)==xb[1]'*yb[2]*yb[1]'*xb[1]


end