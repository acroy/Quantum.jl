#####################################
#DiracVector#########################
#####################################

type DiracVector{S<:Single, T} <: Dirac
	coeffs::SparseMatrixCSC{T, Int} 
	basis::AbstractBasis{S}
	function DiracVector(coeffs, basis)
		if S<:Ket
			@assert size(coeffs)==(length(basis),1) "Dimensions of coefficient array do not match basis"
			return new(coeffs, basis)
		else
			@assert size(coeffs)==(1,length(basis)) "Dimensions of coefficient array do not match basis"
			return new(coeffs, basis)
		end
	end
end

DiracVector{S,T}(coeffs::Array{T}, basis::AbstractBasis{S}) = DiracVector{S,T}(sparse(coeffs), basis)
DiracVector{S,T}(coeffs::SparseMatrixCSC{T}, basis::AbstractBasis{S}) = DiracVector{S,T}(coeffs, basis)
DiracVector(coeffs, s::State...) = DiracVector(coeffs, basis(collect(s)))
dvec = DiracVector

#####################################
#Access Functions####################
#####################################
bsym(d::DiracVector) = bsym(d.basis)

#####################################
#Boolean Functions###################
#####################################
isequal(a::DiracVector, b::DiracVector) = isequal(a.coeffs, b.coeffs) && a.basis==b.basis
==(a::DiracVector, b::DiracVector) = a.coeffs==b.coeffs && a.basis==b.basis
isdual{K<:Ket,B<:Bra}(a::DiracVector{K}, b::DiracVector{B}) = isdual(a.basis, b.basis) && a.coeffs==b.coeffs'
isdual{K<:Ket,B<:Bra}(a::DiracVector{B}, b::DiracVector{K}) = isdual(b,a)
isdual(a::DiracVector, b::DiracVector) = false

#####################################
#Array-like Functions################
#####################################
copy(d::DiracVector) = dvec(copy(d.coeffs), copy(d.basis))

ctranspose(d::DiracVector) = dvec(d.coeffs', d.basis')
size(d::DiracVector, args...) = size(d.coeffs, args...)
getindex{K<:Ket}(d::DiracVector{K}, x) = d.coeffs[x,1]
getindex{B<:Bra}(d::DiracVector{B}, x) = d.coeffs[1,x]
getindex(d::DiracVector, x, y) = d.coeffs[x,y]
setindex!{K<:Ket}(d::DiracVector{K}, y, x) = setindex!(d.coeffs, y, x, 1)
setindex!{B<:Bra}(d::DiracVector{B}, y, x) = setindex!(d.coeffs, y, 1, x)
setindex!(d::DiracVector, y, x...) = setindex!(d.coeffs, y, x...)

for op=(:endof, :ndims, :eltype, :length, :find, :findn, :findnz, :nnz, :ndims, :countnz)
	@eval ($op)(d::DiracVector) = ($op)(d.coeffs)
end

consnz_dvec{K<:Ket}(d::DiracVector{K}, nz_info) = dvec(nz_info[3], basis(d.basis[nz_info[1]]))
consnz_dvec{B<:Bra}(d::DiracVector{B}, nz_info) = dvec(nz_info[3], basis(d.basis[nz_info[2]]))

filternz(d::DiracVector) = consnz_dvec(d, findnz(d))

#####################################
#Dict-like Functions#################
#####################################
getpos(d::DiracVector, s::State) = get(d.basis, s)
function get(d::DiracVector, s::State, notfound)
	try
		return d[getpos(d, s)]
	catch
		return notfound
	end
end

get(d::DiracVector, s::State) = d[getpos(d, s)]
in(s::State, d::DiracVector) = in(s, d.basis)

#####################################
#Show Functions######################
#####################################

summary{K<:Ket}(d::DiracVector{K}) = "$(size(d,1))x1 $(typeof(d))"
summary{B<:Bra}(d::DiracVector{B}) = "1x$(size(d,2)) $(typeof(d))"

function showcompact(io::IO, d::DiracVector)
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else
		tempio = IOBuffer()
		print(tempio, [" + $(d.coeffs[i])$(d.basis[i])" for i in find(d)]...)
		print(io, takebuf_string(tempio)[3:end])
	end
end
function show{S}(io::IO, d::DiracVector{S})
	if length(d)==0
		print(io, "$(typeof(d))[]")
	else	
		println(io, summary(d))
		table = cell(length(d), 2)	
		if length(d)>=52
			for i=1:25
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
			table[26:(length(d)-25),:] = 0 # prevents access to undefined reference
			for i=(length(d)-25):length(d)
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
		else
			for i=1:length(d)
				table[i,1]= d.coeffs[i]
				table[i,2]= d.basis[i]
			end
		end
		temp_io = IOBuffer()
		if S<:Ket
			show(temp_io, table)
		else
			show(temp_io, [transpose(table[:,2]), transpose(table[:,1])])
		end
		io_str = takebuf_string(temp_io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		print(io, io_str)
	end
end

#####################################
#Function-passing Functions##########
#####################################

find(f::Function, d::DiracVector) = find(f, d.coeffs)
findstates(f::Function, d::DiracVector) = find(f, d.basis)

map(f::Function, d::DiracVector)=dvec(map(f, d.coeffs), d.basis)

function map!(f::Function, d::DiracVector)
	map!(f, d.coeffs)
	return d
end

function loadcoeffs!(coeffs, d, newbasis)
	for i=1:length(newbasis)
		coeffs[i] = get(d, newbasis[i])
	end	
	return coeffs
end

function filtercons{S,T}(newbasis::AbstractBasis, d::DiracVector{S,T})
	return dvec(loadcoeffs!(S<:Ket ? Array(T, length(newbasis)) : Array(T, 1, length(newbasis)), d, newbasis), newbasis)
end

filterstates{S,T}(f::Function, d::DiracVector{S,T}) = filtercons(filter(f, d.basis), d)
filtercoeffs(f::Function, d::DiracVector) = filtercons(basis(d.basis[find(map(f, d.coeffs))]), d)

qeval(f::Function, d::DiracVector) = dvec(vcat([qeval(f,x) for x in d.coeffs]...), d.basis) #vcat forces correct typing, but this is hacky

#####################################
#Arithmetic Operations###############
#####################################

for op=(:.*,:.-,:.+,:./,:.^)
	@eval ($op)(a::DiracVector, b::DiracVector) = dvec(($op)(a.coeffs,b.coeffs), a.basis)
	@eval ($op)(n::DiracCoeff, d::DiracVector) = dvec(($op)(n,d.coeffs), d.basis)
	@eval ($op)(d::DiracVector, n::DiracCoeff) = dvec(($op)(d.coeffs,n), d.basis)
end

/(dv::DiracVector, c::DiracCoeff) = dvec(dv.coeffs/c, dv.basis)
*(dv::DiracVector, c::DiracCoeff) = dvec(dv.coeffs*c, dv.basis)
*(c::DiracCoeff, dv::DiracVector) = dvec(c*dv.coeffs, dv.basis)

kron(c::DiracCoeff, d::DiracVector) = c*d
kron(d::DiracVector, c::DiracCoeff) = d*c

-(d::DiracVector) = -1*d

ireduce(s::State, d::DiracVector) = reduce(+,[d[i]*inner(s,d.basis[i]) for i=1:length(d)])
ireduce(d::DiracVector, s::State) = reduce(+,[d[i]*inner(d.basis[i],s) for i=1:length(d)])
ireduce(a::DiracVector, b::DiracVector) = reduce(+,[a[i]*b[j]*inner(a.basis[i],b.basis[j]) for i=1:length(a), j=1:length(b)])
ireduce(a::DiracVector, b::DiracVector, eqbasis::Bool) = reduce(+,[a[i]*b[i]*inner(a.basis[i],b.basis[i]) for i=1:length(a)])

inner{b,T1,T2}(s::State{Bra{b,T1}}, d::DiracVector{Ket{b,T2}}) = get(d, s', 0)
inner{b,T1,T2}(d::DiracVector{Bra{b,T1}},s::State{Ket{b,T2}}) = get(d, s', 0)

inner{B<:Bra, K<:Ket}(s::State{B}, d::DiracVector{K}) = ireduce(s,d)
inner{B<:Bra, K<:Ket}(d::DiracVector{B}, s::State{K}) = ireduce(d,s)
inner{B<:Bra, K<:Ket, N1<:Number, N2<:Number}(a::DiracVector{B, N1}, b::DiracVector{K, N2}) = isdual(a.basis, b.basis) ? (a.coeffs*b.coeffs)[1] : ireduce(a,b)
inner{B<:Bra, K<:Ket}(a::DiracVector{B}, b::DiracVector{K}) = isdual(a.basis, b.basis) ? ireduce(a,b,true) : ireduce(a,b)

for t=(:Bra, :Ket)
	@eval begin
	tensor{A<:($t), B<:($t)}(a::DiracVector{A}, b::DiracVector{B}) = dvec(kron(a.coeffs, b.coeffs), tensor(a.basis, b.basis))
	tensor{A<:($t), B<:($t)}(a::DiracVector{A}, b::State{B}) = dvec(a.coeffs, tensor(a.basis, b))
	tensor{A<:($t), B<:($t)}(a::State{B}, b::DiracVector{A}) = dvec(b.coeffs, tensor(a, b.basis))
	kron{A<:($t), B<:($t)}(a::DiracVector{A}, b::DiracVector{B}) = tensor(a,b)
	kron{A<:($t), B<:($t)}(a::DiracVector{A}, b::State{B}) = tensor(a,b)
	kron{A<:($t), B<:($t)}(a::State{B}, b::DiracVector{A}) = tensor(a,b)
	*{A<:($t), B<:($t)}(a::DiracVector{A}, b::DiracVector{B}) = error("vector multiplication undefined between two $(A)s. Perhaps you meant to use elementwise multiplication (.*)?")
	*{A<:($t), B<:($t)}(a::DiracVector{A}, b::State{B}) = error("vector multiplication undefined between two $(A)s.")
	*{A<:($t), B<:($t)}(a::State{A}, b::DiracVector{B}) = error("vector multiplication undefined between two $(A)s.")
	end
end

kron{K<:Ket, B<:Bra}(a::DiracVector{K}, b::State{B}) = dmat(a.coeffs, a.basis, basis(b))
kron{K<:Ket, B<:Bra}(a::State{K}, b::DiracVector{B}) = dmat(b.coeffs, basis(a), b.basis)
kron{K<:Ket, B<:Bra}(a::DiracVector{K}, b::DiracVector{B}) = dmat(a.coeffs*b.coeffs, a.basis, b.basis)
kron{K<:Ket, B<:Bra}(a::DiracVector{B}, b::DiracVector{K}) = kron(b,a)

*{B<:Bra, K<:Ket}(a::DiracVector{B}, b::DiracVector{K}) = inner(a,b)
*{B<:Bra, K<:Ket}(a::DiracVector{B}, b::State{K}) = inner(a,b)
*{B<:Bra, K<:Ket}(a::State{B}, b::DiracVector{K}) = inner(a,b)

kron(c::DiracCoeff, s::State) = dvec([c], basis(s))
kron(s::State, c::DiracCoeff) = kron(c,s)
*(c::DiracCoeff, s::State) = kron(c,s)
*(s::State, c::DiracCoeff) = kron(c,s)
-(s::State) = kron(-1,s)

function addstate!(d,s)
	res[getpos(d,s)] = 1+get(res, s)
	return res
end

function +{K<:Ket}(d::DiracVector{K}, s::State{K})
	if in(s, d.basis)
		return addstate!(copy(d),s)
	else
		return dvec(vcat(d.coeffs, speye(1)), bcat(d.basis,s))
	end
end

function +{B<:Bra}(d::DiracVector{B}, s::State{B})
	if in(s, d.basis)
		return addstate!(copy(d),s)
	else
		return dvec(hcat(d.coeffs, speye(1)), bcat(d.basis,s))
	end
end

function +{K<:Ket}(s::State{K}, d::DiracVector{K})
	if in(s, d.basis)
		return addstate!(copy(d),s)
	else
		return dvec(vcat(speye(1), d.coeffs), bcat(s,d.basis))
	end
end

function +{B<:Bra}(s::State{B}, d::DiracVector{B})
	if in(s, d.basis)
		return addstate!(copy(d),s)
	else
		return dvec(hcat(speye(1),d.coeffs), bcat(s,d.basis))
	end
end

function add_dvec!(d::DiracVector, y::DiracVector)
	for i=1:length(y)
		d = d+y.basis[i]
		d[getpos(d, y.basis[i])] = get(d, y.basis[i]) + y[i] - 1
	end
	return res
end

for t=(:Bra,:Ket)
	@eval begin 
	function +{b,T1,T2}(x::DiracVector{($t){b,T1}}, y::DiracVector{($t){b,T2}})
		if x.basis==y.basis
			return dvec(x.coeffs+y.coeffs, x.basis)
		else
			return add_dvec!(copy(x), y)
		end
	end
	-{S1<:($t),S2<:($t)}(d::DiracVector{S1}, s::State{S2}) = d+(-s)
	-{S1<:($t),S2<:($t)}(s::State{S1}, d::DiracVector{S2}) = s+(-d)
	-{S1<:($t),S2<:($t)}(a::DiracVector{S1}, b::DiracVector{S2}) = a+(-b)
	end
end

+{K1<:Ket,K2<:Ket}(a::State{K1}, b::State{K2}) = a==b ? 2*a : dvec([1, 1], basis([a,b])) 
+{B1<:Bra,B2<:Bra}(a::State{B1}, b::State{B2}) = a==b ? 2*a : dvec([1 1], basis([a,b])) 
-{K1<:Ket,K2<:Ket}(a::State{K1}, b::State{K2}) = a==b ? dvec(spzero(1),basis(a)) : dvec([1, -1], basis([a,b])) 
-{B1<:Bra,B2<:Bra}(a::State{B1}, b::State{B2}) = a==b ? dvec(spzero(1),basis(a)) : dvec([1 -1], basis([a,b])) 

norm(d::DiracVector, p::Int=2) = reduce(+,map(i->abs(i)^p, d.coeffs))^(1/p) 
norm{S<:Single,N<:Number}(d::DiracVector{S,N}, p::Int=2) = norm(d.coeffs)
normalize(d::DiracVector) = dvec((1/norm(d))*d.coeffs, d.basis)
