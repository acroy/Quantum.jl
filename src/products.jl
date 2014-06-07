#####################################
#InnerProduct/OuterProduct###########
#####################################

immutable OuterProduct <: Dirac
	ket::AbstractState{Ket}
	bra::AbstractState{Bra}
end

immutable InnerProduct <: AbstractScalar
	bra::AbstractState{Bra}
	ket::AbstractState{Ket}
end

#####################################
#Misc Functions######################
#####################################

basislabel(o::OuterProduct) = [basislabel(o.ket), basislabel(o.bra)]
label(o::OuterProduct) = [label(o.ket), label(o.bra)]
basislabel(i::InnerProduct) = [basislabel(i.bra), basislabel(i.ket)]
label(i::InnerProduct) = [label(i.bra), label(o.ket)]
isdual(a::InnerProduct, b::InnerProduct) = isdual(a.ket, b.bra) && isdual(b.ket, a.bra)
conj(i::InnerProduct) = InnerProduct(i.ket', i.bra')
ctranspose(o::OuterProduct) = OuterProduct(o.bra', o.ket')
isdual(a::OuterProduct, b::OuterProduct) = isdual(a.ket, b.bra) && isdual(b.ket, a.bra)

#####################################
#Show Functions######################
#####################################

show(io::IO, o::OuterProduct) = print(io, "$(repr(o.ket))$(repr(o.bra))");
show(io::IO, i::InnerProduct) = print(io, "$(repr(i.bra)) $(repr(i.ket)[2:end])");

#####################################
#Arithmetic Operations###############
#####################################
*(o::OuterProduct, s::AbstractState{Ket}) = DiracVector([(o.bra*s)], tobasis(o.ket))
*(o::OuterProduct, s::AbstractState{Bra}) = OuterProduct(o.ket, o.bra*s)
*(s::AbstractState{Bra}, o::OuterProduct) = DiracVector([(s*o.ket)], tobasis(o.bra))
*(s::AbstractState{Ket}, o::OuterProduct) = OuterProduct(s*o.ket, o.bra)
*(c::DiracCoeff,o::OuterProduct) = c==0 ? 0 : (c==1 ? o : DiracMatrix([c]', tobasis(o.ket), tobasis(o.bra)))
*(o::OuterProduct, c::DiracCoeff) = *(c,o)
-(o::OuterProduct) = -1*o

function +(a::OuterProduct, b::OuterProduct)
	if a==b
		return DiracMatrix(2.0, tobasis(a.ket), tobasis(b.bra))
	elseif samebasis(a,b)
		rowb = tobasis([a.ket, b.ket])
		colb = tobasis([a.bra, b.bra])
		res = DiracMatrix(zeros(2,2)[1:length(rowb), 1:length(colb)], rowb, colb)
		res[getpos(res, a)...] = 1.0
		res[getpos(res, b)...] = 1.0
		return res
	end
end

-(a::OuterProduct, b::OuterProduct) = a+(-b)
