#####################################
#State###############################
#####################################
immutable State{K<:BraKet} <: Dirac
  label::Vector
  kind::Type{K}
end

State(label::Vector) = State(label, Ket)
State{K<:BraKet}(label, kind::Type{K}=Ket) = State([label], kind)
State{K<:BraKet}(label...; kind::Type{K}=Ket) = State([label...], kind)

#####################################
#State Representation################
#####################################
type StateRep{K<:BraKet} <: Quantum
	state::State{K}
	coeffs::Array{Complex{Float64}}
	basis::AbstractBasis{K}
	function StateRep(s::State{K}, coeffs::Array{Complex{Float64}}, basis::AbstractBasis{K})
		if length(basis)==length(coeffs)
			new(s, coeffs, basis)
		elseif length(basis)>length(coeffs)
			error("coefficients unspecified for $(length(basis)-length(coeffs)) basis states")
		else
			error("basis states unspecified for $(length(coeffs)-length(basis)) coefficients")
		end	
	end	
end

StateRep{N<:Number, K<:BraKet}(label, coeffs::Array{N}, basis::AbstractBasis{K}) = StateRep(State(label, K), coeffs, basis)
function StateRep{N<:Number, K<:BraKet}(s::State{K}, coeffs::Array{N}, basis::AbstractBasis{K})
	if size(coeffs, 2)==1 && K==Ket
		StateRep{Ket}(s, convert(Array{Complex{Float64}},vec(coeffs)), basis)
	elseif size(coeffs, 1)==1 && K==Bra
		StateRep{Bra}(s, convert(Array{Complex{Float64}},coeffs), basis)
	else
		error("Dimensions of coefficient array does not match type $K")
	end
end

#####################################
#Functions###########################
#####################################

#imported############################
size(s::State) = size(s.label)

ndims(s::State) = 1

isequal(a::State,b::State) = a.label==b.label && a.kind==b.kind

hash(a::State) = hash(a.label)+hash(a.kind)

ctranspose(s::State) = State(s.label, !s.kind)
ctranspose(s::StateRep) = StateRep(s.state', s.coeffs', s.basis')

getindex(s::State, x) = s.label[x]
getindex(s::StateRep, x) = s.coeffs[x]

setindex!(s::StateRep, y, x) = setindex!(s.coeffs, y, x)
setindex!(s::State, y, x) = setindex!(s.label, y, x)

endof(s::State) = endof(s.label)
endof(s::StateRep) = length(s.coeffs)


label(s::State) = split(repr(s.label), ['[', ']'])[2]
pform(s::State{Bra}, extra="") = isempty(s.label) ? "$lang #undef$extra |" : "$lang $(label(s)) $extra |"

function pform(s::State{Ket}, extra="") 
	if isempty(s.label)
		return "| #undef$extra $rang"
	elseif eltype(s.label)==Any
		return "| $(label(s)) $extra $rang"
	else
		return "| $(label(s)) $extra $rang"
	end
end

function pform(s::State{Bra}, extra="") 
	if isempty(s.label)
		return "$lang #undef$extra |"
	elseif eltype(s.label)==Any
		return "$lang $(label(s)) $extra |"
	else
		return "$lang $(label(s)) $extra |"
	end
end

pform(s::StateRep) = pform(s.state, " ; $(label(s.basis))")

.*(n::Number, s::StateRep) = n*s
.*(s::StateRep, n::Number) = s*n 
.+(s::StateRep, n::Number) = copy(s, s.coeffs.+n)
.+(n::Number, s::StateRep) = copy(s, n.+s.coeffs)
+(arr::Array, s::StateRep) = copy(s, arr+s.coeffs)
+(s::StateRep, arr::Array) = copy(s, s.coeffs+arr)
+(s1::StateRep, s2::StateRep) = s1.basis==s2.basis ? StateRep("$(label(s1)) + $(label(s2))",  s1.coeffs+s2.coeffs, s1.basis) : :($s1+$s2)
.-(s::StateRep, n::Number) = copy(s, s.coeffs.-n)
.-(n::Number, s::StateRep) = copy(s, n.-s.coeffs)
-(arr::Array, s::StateRep) = copy(s, arr-s.coeffs)
-(s::StateRep, arr::Array) = copy(s, s.coeffs-arr)
-(s1::StateRep, s2::StateRep) = s1.basis==s2.basis ? StateRep("$(label(s1)) - $(label(s2))",  s1.coeffs-s2.coeffs, s1.basis) : :($s1-$s2)
./(s::StateRep, n::Number) = s/n
./(n::Number, s::StateRep) = copy(s, n./s.coeffs)
.^(n::Number, s::StateRep) = copy(s, n.^s.coeffs)
.^(s::StateRep, n::Number) = copy(s, s.coeffs.^n)

/(s::StateRep, n::Number) = copy(s, s.coeffs/n)

*(n::Number, s::StateRep) = copy(s, n*s.coeffs) 
*(s::StateRep, n::Number) = *(n, s)
*{N<:Number}(arr::Array{N, 2}, s::StateRep{Ket}) = size(arr,1)==1 ? (arr*s.coeffs)[1] : copy(s, arr*s.coeffs)
*{N<:Number}(s::StateRep{Bra}, arr::Array{N, 1}) = (s.coeffs*arr)[1]
*{N<:Number}(s::StateRep{Bra}, arr::Array{N, 2}) = size(arr,2)==1 ? (s.coeffs*arr)[1] : copy(s, s.coeffs*arr)
*{N<:Number}(arr::Array{N}, s::StateRep{Bra}) = length(arr)==length(s) ? OperatorRep(arr*s.coeffs, s.basis) : throw(DimensionMismatch)
*{N<:Number}(s::StateRep{Ket}, arr::Array{N}) = size(arr,2)==length(s) ? OperatorRep(s.coeffs*arr, s.basis) : throw(DimensionMismatch)

*(a::StateRep{Bra}, b::StateRep{Ket}) = samelabels(a.basis, b.basis) ? (a.coeffs*b.coeffs)[1] : :($a*$b)
*(a::StateRep{Ket}, b::StateRep{Ket}) = StateRep(a.state*b.state, kron(a.coeffs, b.coeffs), a.basis*b.basis)
*(a::StateRep{Bra}, b::StateRep{Bra}) = StateRep(a.state*b.state, kron(a.coeffs, b.coeffs), a.basis*b.basis)
*(a::StateRep{Ket}, b::StateRep{Bra}) = OperatorRep(a.coeffs*b.coeffs, a.basis, b.basis)
*(sr::StateRep{Bra}, s::State{Ket}) = get(sr, s')
*(s::State{Bra}, sr::StateRep{Ket}) = get(sr, s')
*{K<:BraKet}(s1::State{K}, s2::State{K}) = tensor(s1, s2)
*(s1::State{Bra}, s2::State{Ket}) = InnerProduct(s1, s2)
*(s1::State{Ket}, s2::State{Bra}) = OuterProduct(s1, s2)

copy(s::StateRep, coeffs=copy(s.coeffs)) = StateRep(s.state, coeffs, s.basis)
find(s::StateRep) = find(s.coeffs)
length(s::StateRep) = length(s.coeffs)
length(s::State) = length(s.label)

function get(s::StateRep, label, notfound) 
	ind = get(s.basis, label, 0)
	if ind==0
		return notfound
	else
		return s[ind]
	end
end
get(s::StateRep, skey::State, notfound) =  get(s, skey.label, notfound)
get(s::StateRep, key) = get(s, key, 0)
norm(s::StateRep) = norm(s.coeffs)

function map!(f::Function, s::StateRep)
	s.coeffs = map!(f, s.coeffs)
	return s
end 

map(f::Function, s::StateRep) = map!(f, copy(s))

show(io::IO, s::State) = print(io, pform(s))
showcompact(io::IO, s::StateRep) = print(io, pform(s))

function show(io::IO, s::StateRep)
	println("$(typeof(s)) $(pform(s)):")
	if any(s.coeffs.!=0)
		table = cell(length(s.coeffs), 2)	
		if length(s.coeffs)>=52
			for i=1:25
				table[i,1]= s.coeffs[i]
				table[i,2]= s.basis[i]
			end
			table[26:(length(s.coeffs)-25),:] = 0 # prevents access to undefined reference
			for i=(length(s.coeffs)-25):length(s.coeffs)
				table[i,1]= s.coeffs[i]
				table[i,2]= s.basis[i]
			end
		else
			for i=1:length(s.coeffs)
				table[i,1]= s.coeffs[i]
				table[i,2]= s.basis[i]
			end
		end
		temp_io = IOBuffer()
		if kind(s)==Ket
			show(temp_io, table)
		else
			show(temp_io, [transpose(table[:,2]), transpose(table[:,1])])
		end
		io_str = takebuf_string(temp_io)
		io_str = io_str[searchindex(io_str, "\n")+1:end]
		print(io_str)
	else
		print("(all coefficients are zero)")
	end
end

#exported############################
kind(s::State) = s.kind
kind(s::StateRep) = s.state.kind

function statevec{K<:BraKet}(v::Vector, kind::Type{K}=Ket)
	svec = Array(State{kind}, length(v))
	for i=1:length(v)
		svec[i] = State(v[i], kind)
	end
	return svec
end
function statevec{K<:BraKet}(arr::Array, kind::Type{K}=Ket)
	svec = Array(State{kind}, size(arr,1))
	for i=1:size(arr, 1)
		svec[i] = State(vec(arr[i,:]), kind)
	end
	return svec
end

labeldelta(s1::State, s2::State) = s1.label==s2.label ? 1 : 0 

tensor() = nothing
tensor{K<:BraKet}(s::State{K}...) = State(vcat([i.label for i in s]...), K)
tensor{S<:State}(state_arrs::Array{S}...) = statejoin(crossjoin(state_arrs...))

statejoin{S<:State}(v::Vector{S}) = tensor(v...)
statejoin{S<:State}(v::Vector{S}...) = broadcast(tensor, v...)
function statejoin{S<:State}(state_arr::Array{S}) 
	result = statejoin(state_arr[:,1], state_arr[:,2])
	for i=3:size(state_arr, 2)
		result = statejoin(result, state_arr[:,i])
	end
	return result
end

separate{K<:BraKet}(s::State{K}) = statevec(s.label, K)
separate{S<:State}(v::Vector{S}) = hcat(map(separate, v)...).'

state(s::StateRep) = s.state

normalize(v::Vector) = (1/norm(v))*v

function normalize!(s::StateRep) 
	s.coeffs=normalize(s.coeffs)
	return s
end

normalize(s::StateRep) = normalize!(copy(s))

function mapmatch!(f_coeffs::Function, f_states::Function, s::StateRep)
	matched_states = filter(f_states, s.basis)	
	for i=1:length(matched_states)
		s[get(s.basis, matched_states[i], nothing)] = apply(f_coeffs, get(s, matched_states[i], nothing))
	end
	return s
end
mapmatch(f_coeffs::Function, f_states::Function, s::StateRep) = mapmatch!(f_coeffs, f_states, copy(s))

filterstates(f::Function, s::StateRep) = mapmatch((x)->0, f, s)
filterstates!(f::Function, s::StateRep) = mapmatch!((x)->0, f, s)
filtercoeffs(f::Function, s::StateRep) = filtercoeffs!(f, copy(s))
filtercoeffs!(f::Function, s::StateRep) = map!(x->apply(f, x) ? x : 0, s)

