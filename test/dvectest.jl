module DVeTest
using Base.Test 
using Quantum

xb = basis(:X,[1:5])
d = DiracVector([1:5], xb)
@test d==d''
@test d[:]==d.coeffs
@test d'*d == 55
@test d*d' == DiracMatrix(kron(d.coeffs, d.coeffs'), d.basis, d.basis')
@test isdual(d,d')
@test filterstates(s->label(s)%2==0, d) == DiracVector([2,4], basis(xb[2], xb[4]))
@test filtercoeffs(c->c%2==0, d) == DiracVector([2,4], basis(xb[2], xb[4]))
@test map(x->x*2, d)==d+d
@test xb[1]'*d == 1
@test d'*xb[1] == 1
@test norm(normalize(d)) < 1+1e-8 && norm(normalize(d)) > 1-1e-8
@test reduce(+,[i*xb[i]' for i=1:length(xb)]) == d'
@test (1/sqrt(2)*xb[1]) + (1/sqrt(2)*xb[2]) == dvec([1/sqrt(2), 1/sqrt(2)], basis(xb[1], xb[2]))
@test ket(:X,1)+ket(:X,1)+ket(:X,1) = dvec([3], basis(xb[1]))
@test ket(:X,1)+(bra(:S,"1")*ket(:X, 1)) * ket(:X,1)==dvec([1+1*bra(:S,"1")*ket(:X,1)], basis(ket(:X,1)))

end