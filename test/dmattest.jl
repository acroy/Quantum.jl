module DMatTest
using Base.Test 
using Quantum

xb = basis(:X,[1:3])
xbplus = basis(:X,[1:4])
xbminus = basis(:X,[0:3])

d = dmat(eye(3), xb)
@test d==d''
@test d[:,:]==d.coeffs
@test d*ket(:X,2) == dvec([0.0,1.0,0.0], xb)

@test d+ket(:X,4)*bra(:X,4) == dmat(eye(4),xbplus)
@test d+ket(:X,4)*bra(:X,3) == dmat([eye(3), [0 0 1]],xbplus, xb')
@test d+ket(:X,3)*bra(:X,4) == dmat([eye(3) [0,0,1]],xb, xbplus')
@test d+ket(:X,3)*bra(:X,3) == dmat((coeff=eye(3); coeff[3,3]=2.0; coeff),xb)

@test ket(:X,0)*bra(:X,0)+d == dmat(eye(4),xbminus)
@test ket(:X,0)*bra(:X,1)+d == dmat([[1 0 0], eye(3)],xbminus, xb')
@test ket(:X,1)*bra(:X,0)+d == dmat([[1,0,0] eye(3)],xb, xbminus')
@test ket(:X,1)*bra(:X,1)+d == dmat((coeff=eye(3); coeff[1,1]=2.0; coeff),xb)

@test d+d==2*d
@test isdual(d,d')
@test map(x->x*2, d)==d+d
@test commutator(d,d)==dmat(zeros(3,3), xb)

end