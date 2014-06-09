module ScalarTest
using Base.Test 
using Quantum

xb = Basis([1:5], "X")
yb = Basis([1:5], "Y")
xyb = xb*yb

@test (1+(xb[1]'*yb[1])).ex.args==[:+, 1, InnerProduct(xb[1]', yb[1])]

end