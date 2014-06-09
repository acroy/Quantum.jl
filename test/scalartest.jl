module ScalarTest
using Base.Test 
using Quantum

xb = Basis([1:5], "X")
yb = Basis([1:5], "Y")

@test (1+(xb[1]'*yb[1])).ex.args==[:+, 1, InnerProduct(xb[1]', yb[1])]
for op=(:+, :-, :*, :/)
	@test @eval qeval((b,k)->b==yb[1]' ? 2.5 : 3.4, ($op)(xb[1]'*yb[1], yb[1]'*xb[1])) == ($op)(3.4,2.5)
end

end