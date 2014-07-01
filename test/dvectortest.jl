module DVectorTest
using Base.Test 
using Quantum

xb = Basis(:X,[1:5])
d = DiracVector([1:5], xb)
@test d==d''
@test d[:]==d.coeffs
@test d'*d == 55
@test d*d' == DiracMatrix(kron(d.coeffs, d.coeffs'), d.basis, d.basis')
@test isdual(d,d')
@test filterstates(s->label(s)%2==0, d) == DiracVector([2,4], tobasis(xb[2], xb[4]))
@test filtercoeffs(c->c%2==0, d) == DiracVector([2,4], tobasis(xb[2], xb[4]))
@test map(x->x*2, d)==d+d
@test mapmatch(s->label(s)%2==0, x->x*2, d)==DiracVector([1,4,3,8,5], xb)
@test xb[1]'*d == 1
@test d'*xb[1] == 1
@test norm(normalize(d)) < 1+1e-8 && norm(normalize(d)) > 1-1e-8
@test reduce(+,[i*xb[i]' for i=1:length(xb)]) == d'

end