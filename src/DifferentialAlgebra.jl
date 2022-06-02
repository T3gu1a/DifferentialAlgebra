
#module DifferentialAlgebra

import Base

using Markdown

include("GenericTypes.jl")

export DifferentialPolynomialRing


#--------

function form_derivative(varname::String, order::Integer)
    return "$(varname)^($order)"
end

function form_partial_derivative(varname::String, order::Array{Int64,1})
    return "$(varname)^($order)"
end


function list_lex_monomials(b::Array{Int64, 1})
	if length(b)==1 
		return [[j] for j in 0:b[1]]
	else
		return [append!([i],m) for i in 0:b[1] for m in list_lex_monomials(b[2:length(b)])]
	end
end

function rev_list_lex_monomials(b::Array{Int64, 1})
	if length(b)==1 
		return [[j] for j in b[1]:-1:0]
	else
		return [append!([i],m) for i in b[1]:-1:0 for m in rev_list_lex_monomials(b[2:length(b)])]
	end
end

#------------------------------------------------------------------------------

mutable struct DifferentialPolyRing <: DifferentialRing
    base_ring::AbstractAlgebra.Ring
    poly_ring::MPolyRing
    max_ord::Array{Int64,1}
    varnames::Array{String, 1}
    derivation::Array{Dict{Any, Any},1}
	
	function DifferentialPolyRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}, ranking::Symbol, ordering::Symbol, max_ord::Array{Int64,1})
        if length(max_ord)==1
			if ranking == :var_deriv
				all_varnames = [form_derivative(v, ord) for v in varnames for ord in max_ord[1]:-1:0]
			elseif ranking == :deriv_var
				all_varnames = [form_derivative(v, ord) for ord in max_ord[1]:-1:0 for v in varnames]
			end
			if ordering != :lex
				poly_ring, _ = AbstractAlgebra.PolynomialRing(R, all_varnames,ordering=ordering)
			else
				poly_ring, _ = AbstractAlgebra.PolynomialRing(R,all_varnames)
			end
			derivation = [Dict()]
			for v in varnames
				for ord in 0:(max_ord[1] - 1)
					derivation[1][str_to_var(form_derivative(v, ord), poly_ring)] = 
						str_to_var(form_derivative(v, ord + 1), poly_ring)
				end
			end
			return new(R, poly_ring, max_ord, varnames, derivation)
		else
			L = rev_list_lex_monomials(max_ord)
			all_varnames = [form_partial_derivative(v, ord) for v in varnames for ord in L]
			poly_ring, _ = AbstractAlgebra.PolynomialRing(R, all_varnames)
			derivation = [Dict() for j in 1:length(max_ord)]
			for v in varnames
				for i in 1:length(max_ord)
					for ord in L
						if ord[i]<max_ord[i] 
							updord = ord
							updord[i]+=1
							derivation[i][str_to_var(form_partial_derivative(v, ord), poly_ring)] = 
							str_to_var(form_partial_derivative(v, updord), poly_ring)
						end
					end
				end
			end
			return new(R, poly_ring, max_ord, varnames, derivation)
		end
    end
end

function DifferentialPolynomialRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}; ranking::Symbol = :var_deriv, ordering::Symbol = :lex, max_ord::Array{Int64,1}=[20])
	R = DifferentialPolyRing(R, varnames, ranking, ordering, max_ord)
    return R, Tuple([DiffIndet(R,v) for v in varnames])
end

function AbstractAlgebra.gens(R::DifferentialPolyRing)
	return [gen(R.poly_ring,(j-1)*(R.max_ord+1)+1) for j in 1:length(R.varnames)] 
    #return [DiffPoly(R, str_to_var(form_derivative(v, 0), R.poly_ring)) for v in R.varnames]
end

#-----

mutable struct DiffPoly <: DifferentialRingElem
    parent::DifferentialPolyRing
    algdata::AbstractAlgebra.MPolyElem

    function DiffPoly(R::DifferentialPolyRing, alg_poly::AbstractAlgebra.MPolyElem)
        #return new(R, parent_ring_change(alg_poly, R.poly_ring))
		return new(R, alg_poly)
    end
	
	#function DiffPoly(a::DiffIndet)
	#	return DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring))
	#end
end

#copy
mutable struct DiffIndet <: DifferentialRingElem
    parent::DifferentialPolyRing
	varname::String

    function DiffIndet(R::DifferentialPolyRing, var_name::String)
        return new(R, var_name)
    end
end

function Base.parent(a::DiffPoly)
    return a.parent
end

function Base.parent(a::DiffIndet)
    return a.parent
end

function algdata(a::DiffPoly)
    return a.algdata
end

AbstractAlgebra.elem_type(::Type{DifferentialPolyRing}) = DifferentialRingElem

AbstractAlgebra.parent_type(::Type{DifferentialRingElem}) = DifferentialPolyRing

AbstractAlgebra.parent_type(::Type{DiffPoly}) = DifferentialPolyRing

#------------------------------------------------------------------------------
#comparisons

#equality

function DiffPoly_to_Bool(a::DiffPoly)
	return convert(Bool,parse(Int64,"$(a)"))
end

function Base.:(==)(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(algdata(a) == algdata(b)))
end

function Base.:(==)(a::DiffIndet, b::DifferentialRingElem)
	return false
end

function Base.:(==)(a::DifferentialRingElem, b::DiffIndet)
	return false
end

function Base.:(==)(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) == DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
end


#lessorgreather

function Base.:<(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(algdata(a) < algdata(b)))
end

function Base.:<(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < b))
end

function Base.:>(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) > b))
end

function Base.:<(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(a < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
end

function Base.:<(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
end

function Base.:>(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return DiffPoly_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) > DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
end

#---

#-------------------------------------------------------------------------

function Base.:+(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) + algdata(b))
end

function Base.:+(a::DifferentialRingElem, b)
    return parent(a)(algdata(a) + b)
end

function Base.:+(a, b::DifferentialRingElem)
    return parent(b)(algdata(b) + a)
end

function Base.:+(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + b)
end

function Base.:+(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(a + DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:+(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:+(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return parent(a)(algdata(a) + b)
end

function Base.:+(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return parent(b)(a + algdata(b))
end

function Base.:+(a::DiffIndet, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + b)
end

function Base.:+(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
    return parent(b)(a + DiffPoly(parent(b), str_to_var(form_derivative(b.varname, 0), parent(b).poly_ring)))
end

#function Base.:+(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational, Nemo.fmpq})
#    return parent(a)(algdata(a) + b)
#end

function AbstractAlgebra.addeq!(a::DifferentialRingElem, b::DifferentialRingElem)
    a = a + b
end

#------------------------------------------------------------------------------

function Base.iszero(a::DifferentialRingElem)
    return algdata(a) == 0
end

#------------------------------------------------------------------------------

function Base.show(io::IO, R::DifferentialPolyRing)
    print(io, "Differential Polynomial Ring in " * join(R.varnames, ", ") * " over $(R.base_ring)")
end

function Base.show(io::IO, p::DifferentialRingElem)
    show(io, algdata(p))
end

function Base.show(io::IO, p::DiffIndet)
    show(io, p.varname)
end

#-----------------

function Base.zero(R::DifferentialPolyRing)
    return R(0)
end

function Base.one(R::DifferentialPolyRing)
    return R(1)
end

#------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

#------------------------------------------------------------------------------

function d_aux(p::MPolyElem, der::Dict{Any, Any})
    result = zero(parent(p))
    for v in vars(p)
        if !(v in keys(der))
            throw(DomainError("No derivative defined for $v. Most likely you have exceeded the maximal order."))
        end
        result += der[v] * derivative(p, v)
    end
    return result
end

function d(a::DiffPoly)
    return DiffPoly(parent(a), d_aux(algdata(a), parent(a).derivation[1]))
end

function d(a::DiffIndet)
    return DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 1), parent(a).poly_ring))
end

#- today
function d(a::Union{AbstractFloat, Integer, Rational})
    return 0
end
#---

#-- 
function d(a::DifferentialRingElem, ord::Integer)
    if ord == 0
        return a
    end
    return d(d(a), ord - 1)
end

function d(a::Union{Integer,Rational}, ord::Integer)
    return 0
end

#------------------------------------------------------------------------------

function Base.hash(a::DifferentialRingElem)
    return hash(algdata(a))
end

#------------------------------------------------------------------------------

function (R::DifferentialPolyRing)(b)
    return DiffPoly(R, R.poly_ring(b))
end

function (R::DifferentialPolyRing)(b::DiffPoly)
    return DiffPoly(R, algdata(b))
end

function (R::DifferentialPolyRing)()
    return zero(R)
end

#------------------------------------------------------------------------------

function Base.:*(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) * algdata(b))
end

function Base.:*(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * algdata(b))
end

function Base.:*(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
    return parent(a)(algdata(a) * DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:*(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:*(a::DiffIndet, b::Union{AbstractFloat, Integer, Rational})
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * b)
end

function Base.:*(a::Union{AbstractFloat, Integer, Rational}, b::DiffIndet)
    return parent(b)(a * DiffPoly(parent(b), str_to_var(form_derivative(b.varname, 0), parent(b).poly_ring)))
end

function Base.:*(a::RingElem, b::DifferentialRingElem)
    if typeof(a) <: DifferentialRingElem
        return parent(a)(algdata(a) * algdata(b))
    end
    return parent(b)(a * algdata(b))
end

function Base.:*(a::DifferentialRingElem, b::RingElem)
    if typeof(b) <: DifferentialRingElem
        return parent(b)(algdata(a) * algdata(b))
    end
    return parent(a)(algdata(a) * b)
end

#function Base.:*(a, b::DifferentialRingElem)
#    return parent(b)(algdata(b) * a)
#end

function Base.:*(a::Union{AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return parent(b)(a * algdata(b))
end

function AbstractAlgebra.mul!(a::DifferentialRingElem, b::DifferentialRingElem, c::DifferentialRingElem)
    a = b * c
end

#------------------------------------------------------------------------------

function Base.:-(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
    return parent(a)(algdata(a) - algdata(b))
end

function Base.:-(a::DifferentialRingElem, b)
    return parent(a)(algdata(a) - b)
end

function Base.:-(a, b::DifferentialRingElem)
    return parent(b)(-algdata(b) + a)
end

function Base.:-(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return parent(a)(algdata(a) - b)
end

function Base.:-(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return parent(b)(a - algdata(b))
end

function Base.:-(a::DiffIndet, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - b)
end

function Base.:-(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
    return parent(b)(a - DiffPoly(parent(b), str_to_var(form_derivative(b.varname, 0), parent(b).poly_ring)))
end

function Base.:-(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:-(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - b)
end

function Base.:-(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(a - DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

#------------------------------------------------------------------------------

function Base.:-(a::DifferentialRingElem)
    return parent(a)(-algdata(a))
end

#-today

function Base.:-(a::DiffIndet)
    return (-1)*a
end

function Base.://(a::DiffPoly, b::Union{Integer, Rational})
	return a*(1//b)
end

function Base.://(a::DiffIndet, b::Union{Integer, Rational})
	return a*(1//b)
end

function Base.:/(a::DiffPoly, b::Union{AbstractFloat,Integer,Rational})
	return a*(1/b)
end

function Base.:/(a::DiffIndet, b::Union{AbstractFloat,Integer,Rational})
	return a*(1/b)
end

#-----------------------------------------------------------------------------

function Base.isless(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end


#------------------------------------------------------------------------------

#end # module