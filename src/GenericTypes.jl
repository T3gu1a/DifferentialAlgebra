import Oscar
using AbstractAlgebra


abstract type DifferentialRing <: AbstractAlgebra.Ring end
abstract type DifferentialField <: AbstractAlgebra.Field end
abstract type DifferentialRingElem <: AbstractAlgebra.RingElem end
abstract type DifferentialFieldElem <: DifferentialRingElem end 

#DiffRing --> DifferentialRing
#DiffField --> DifferentialField
#DiffRingElem --> DifferentialRingElem
#DiffFieldElem --> DifferentialFieldElem