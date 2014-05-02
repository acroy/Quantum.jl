Quantum.jl API
===
I. Type Implementations
===
1. AbstractTypes
---
The following list contains the abstract types referenced in this documentation:

	abstract Dirac
	abstract BraKet <: Dirac
	abstract Bra <: BraKet
	abstract Ket <: BraKet
	abstract AbstractBasis{K<:BraKet} <: Dirac

The not operator `!` can be applied to `Bra` and `Ket` types to alternate between the two (since they are duals of each other):
	!(K::Type{Ket}) = Bra
	!(B::Type{Bra}) = Ket

2. State
---
*Description*
A State is a type of an object that has only a label (stored as a Vector) and a specification of whether it belongs to Ket-space or Bra-space (this property is referred to as “kind”). A state is parameterized by its kind; thus, states are either of type State{Bra} or State{Ket}. State is a subtype of Quantum.

*Definition*
	immutable State{K<:BraKet} <: Quantum
	  label::Vector
	  kind::Type{K}
	end

*Constructors*
	State(label::Vector) = State(label, Ket)
	State{K<:BraKet}(label, kind::Type{K}=Ket) = State([label], kind)
	State{K<:BraKet}(label...; kind::Type{K}=Ket) = State([label...], kind)