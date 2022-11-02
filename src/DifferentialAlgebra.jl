
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

#not used
function list_lex_monomials(b::Array{Int64, 1})
	if length(b)==1 
		return [[j] for j in 0:b[1]]
	else
		return [append!([i],m) for i in 0:b[1] for m in list_lex_monomials(b[2:length(b)])]
	end
end

function list_degrevlex_monomials(b::Array{Int64, 1})
	if length(b)==1 
		return [[j] for j in 0:b[1]]
	else 
		return sort!(rev_list_lex_monomials(b),lt=degrevlex_isless,rev=true)
	end
end

function rev_list_lex_monomials(b::Array{Int64, 1})
	if length(b)==1 
		return [[j] for j in b[1]:-1:0]
	else
		return [append!([i],m) for i in b[1]:-1:0 for m in rev_list_lex_monomials(b[2:length(b)])]
	end
end

function degrevlex_isless(v1::Array{Int64, 1},v2::Array{Int64, 1})
	s1 = sum(v1)
	s2 = sum(v2)
	if s1 == s2
		return v1 < v2
	else
		return s1 < s2
	end
end


#------------------------------------------------------------------------------

mutable struct DifferentialPolyRing <: DifferentialRing
    base_ring::AbstractAlgebra.Ring
    poly_ring::MPolyRing
    max_ord::Array{Int64,1}
    varnames::Array{String, 1}
    derivation::Array{Dict{Any, Any},1}
	
	function DifferentialPolyRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}, ranking_dependent::Symbol, ranking_independent::Symbol, ordering::Symbol, max_ord::Array{Int64,1})
        if length(max_ord)==1
			if ranking_dependent == :var_deriv
				all_varnames = [form_derivative(v, ord) for v in varnames for ord in max_ord[1]:-1:0]
			elseif ranking_dependent == :deriv_var
				all_varnames = [form_derivative(v, ord) for ord in max_ord[1]:-1:0 for v in varnames]
			else
				throw(DomainError("Please use the Julia symbol :var_deriv or :deriv_var for ranking_dependent"))
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
			if ranking_independent == :lex
				L = rev_list_lex_monomials(max_ord)
			elseif ranking_independent == :degrevlex
				L = list_degrevlex_monomials(max_ord)
			end
			if ranking_dependent == :var_deriv
				all_varnames = [form_partial_derivative(v, ord) for v in varnames for ord in L]
			else
				all_varnames = [form_partial_derivative(v, ord) for ord in L for v in varnames]
			end
			if ordering != :lex
				poly_ring, _ = AbstractAlgebra.PolynomialRing(R, all_varnames,ordering=ordering)
			else
				poly_ring, _ = AbstractAlgebra.PolynomialRing(R,all_varnames)
			end
			derivation = [Dict() for j in 1:length(max_ord)]
			for v in varnames
				for i in 1:length(max_ord)
					for ord in L
						if ord[i]<max_ord[i] 
							updord = copy(ord)
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

function DifferentialPolynomialRing(R::AbstractAlgebra.Ring, varnames::Array{String, 1}; ranking_dependent::Symbol = :var_deriv, ranking_independent::Symbol = :lex, ordering::Symbol = :lex, max_ord::Array{Int64,1}=[20])
	R = DifferentialPolyRing(R, varnames, ranking_dependent, ranking_independent, ordering, max_ord)
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

#the comparison returns a boolean integer (1 or 0)
#which the following code converts to true or false
function Binary_to_Bool(a::DiffPoly)
	return convert(Bool,parse(Int64,"$(a)"))
end

function Base.:(==)(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
	return Binary_to_Bool(parent(a)(algdata(a) == algdata(b)))
end

function Base.:(==)(a::DiffIndet, b::DifferentialRingElem)
	return false
end

function Base.:(==)(a::DifferentialRingElem, b::DiffIndet)
	return false
end

function Base.:(==)(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, Z), parent(a).poly_ring)) == DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, Z), parent(a).poly_ring))))
	else
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) == DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
	end
end


#lessorgreather

function Base.:<(a::DifferentialRingElem, b::DifferentialRingElem)
    check_parent(a, b)
	return Binary_to_Bool(parent(a)(algdata(a) < algdata(b)))
end

function Base.:<(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, Z), parent(a).poly_ring)) < b))
	else
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < b))
	end
end

function Base.:>(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, Z), parent(a).poly_ring)) > b))
	else
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) > b))
	end
end

function Base.:<(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(a < DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, Z), parent(a).poly_ring))))
	else
		return Binary_to_Bool(parent(a)(a < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
	end
end

function Base.:<(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, Z), parent(a).poly_ring)) < DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, Z), parent(a).poly_ring))))
	else
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
	end
end

function Base.:>(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		Z=zeros(Int64,length(parent(a).derivation))
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, Z), parent(a).poly_ring)) > DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, Z), parent(a).poly_ring))))
	else
		return Binary_to_Bool(parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) > DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring))))
	end
end

#--------------------------------------------------------------------------------------

function lex_leader_index(p::DiffPoly)
	#Assuming that the terms (monomials) in p are sorted with respect to lex
    return findfirst(x -> x!=0, leading_exponent_vector(p.algdata))
end

function lex_leader(p::DiffPoly)
    return gens(parent(p).poly_ring)[lex_leader_index(p)]
end

function lex_leader_index_degree(p::DiffPoly)
	#Assuming that the terms (monomials) in p are sorted with respect to lex
    v=leading_exponent_vector(p.algdata)
	i= findfirst(x -> x!=0, v)
	return i,v[i]
end

#--

#------------ Independently of the chosen monomial ordering -----------------------

#compute the leader as the maximum of the maximums of the variables in each term
#with respect to the chosen ranking
function leader(P::Union{DiffPoly,DiffIndet})
	p=P+0 #working with DiffIndet as DiffPoly
	V = filter!(v->v!=[],map(vars,terms(p.algdata)))
	if V==[]
		return 1
	else
		return  parent(p)(maximum(map(v->maximum(v),V)))
	end
end

#for internal use
function lead_er(P::Union{DiffPoly,DiffIndet})
	p=P+0
	V = filter!(v->v!=[],map(vars,terms(p.algdata)))
	if V==[]
		return 1
	else
		return  maximum(map(v->maximum(v),V))
	end
end

#we have a very easy way to compute the separant :)
#derivative is a command from the AbstractAlgebra package
function separant(P::Union{DiffPoly,DiffIndet})
	p=P+0
	return parent(p)(derivative(p.algdata,lead_er(p)))
end

#for internal use
function sep_arant(p::DiffPoly)
	return derivative(p.algdata,lead_er(p))
end

function leader_degree_initial(P::Union{DiffPoly,DiffIndet})
	p=P+0
	ld = lead_er(p)
	
	#hasleader is a function that select terms having the leader
	hasleader(t) = (ld in vars(t)) ? true : false
	
	#select all terms with the leader (collect is important here in order to use findall)
	pdata=collect(terms(p.algdata))
	term_with_ld = map(i->pdata[i],findall(hasleader, pdata))
	
	#index of the leader degree in the exponent vector
	index_ld_deg=findfirst(x->x!=0,degrees(ld))
	
	#degrees of the leader for all terms having the leader
	deg_term_with_ld=map(t->degrees(t)[index_ld_deg], term_with_ld)
	
	#we have the leader degree and the index of the term containing the initial
	ld_deg, init_index = findmax(deg_term_with_ld) 
	
	#the initial is then computed by removing the leader from the term with index init_index
	_, init = remove(term_with_ld[init_index],ld)
	
	#return the results
	return ld_deg, parent(p)(init)
end

#extracting the order from an indeterminate
function indet_order(v::DiffPoly)
	v_str="$v"
	ord=v_str[findfirst('(',v_str)+1:findfirst(')',v_str)-1]
	if length(parent(v).derivation)>1
		return parse.(Int64, split(chop(ord; head=1, tail=1), ','))
	else
		return parse(Int64, ord)
	end
end

#funny: check_parent(x,x^1) is (was) not true!!
function Base.:^(a::Union{DiffPoly,DifferentialRingElem}, i::Integer)
	return parent(a)(algdata(a)^i)
end

function Base.:^(a::DiffIndet, i::Integer)
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring))^i)
	else
		return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring))^i)
	end
end

function leader_isgreater(l1::Union{DiffPoly,Integer},l2::Union{DiffPoly,Integer})
	if typeof(l1)==DiffPoly && typeof(l2)==DiffPoly
		return l1>l2
	elseif typeof(l1)==DiffPoly
		return true
	else
		return false
	end
end

#Alrogithm for one differential indeterminate
function diffreduction(p::Union{DiffPoly,DiffIndet}, q::Union{DiffPoly,DiffIndet})
    """
    Performs a differential reduction of g with respect to f or vice verca
	when there is only one differential indeterminate
    """
    check_parent(p, q)
	if length(parent(p).varnames)>1
		#head reduction (to be implemented)
        throw(DomainError("More than one differential indeterminate. To be defined..."))
    end
	leadg = leader(p)
	leadf = leader(q)
    g = (leadg > leadf) ? p+0 : q+0
	f = (leadf < leadg) ? q+0 : p+0
	while leader_isgreater(leadg, leadf)
		deg_g, init_g = leader_degree_initial(g)
		dord = indet_order(leadg)-indet_order(leadf)
		g = separant(f)*g - init_g*(leadg^(deg_g-1))*d(f,dord)
		leadg = leader(g)
	end
	return parent(p)(divrem(g.algdata,f.algdata)[2])
end

#--------------------------------------------------------------------------------

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
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) + b)
	else
		return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + b)
	end
end

function Base.:+(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(a + DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
	return parent(a)(a + DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:+(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) + DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:+(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return parent(a)(algdata(a) + b)
end

function Base.:+(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return parent(b)(a + algdata(b))
end

function Base.:+(a::DiffIndet, b::Union{RingElem, AbstractFloat, Integer, Rational})
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) + b)
	end
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) + b)
end

function Base.:+(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
	if length(parent(b).derivation)>1
		return parent(b)(a + DiffPoly(parent(b), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(b).poly_ring)))
	end
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
    print(io, p.varname)
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
	if length(parent(a).derivation)>1
		#if in the partial case then ERROR
		throw(DomainError("Missing partial orders for the derivation"))
	end
    return DiffPoly(parent(a), d_aux(algdata(a), parent(a).derivation[1]))
end

function d(a::DiffIndet)
	if length(parent(a).derivation)>1
		#if in the partial case then ERROR
		throw(DomainError("Missing partial orders for the derivation"))
	end
    return DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 1), parent(a).poly_ring))
end

#-
function d(a::Union{AbstractFloat, Integer, Rational})
	if length(parent(a).derivation)>1
		#if in the partial case then ERROR
		throw(DomainError("Missing partial orders for the derivation"))
	end
    return 0
end
#---

#-- 
function d(a::DifferentialRingElem, ord::Integer)
	if length(parent(a).derivation)>1
		#if in the partial case then ERROR
		throw(DomainError("Missing partial orders for the derivation"))
	end
    if ord == 0
        return a
    end
    return d(d(a), ord - 1)
end

function d(a::Union{Integer,Rational}, ord::Integer)
	if length(parent(a).derivation)>1
		#if in the partial case then ERROR
		throw(DomainError("Missing partial orders for the derivation"))
	end
    return 0
end
#--------------------------------- extending the bound ---------

function embedDiffPoly(p::DiffPoly, R::DifferentialPolyRing)
	h = Oscar.hom(parent(p), R, gens(parent(p)))
	return h(p)
end

#----------partial d -----------------------------------

function d(a::DiffPoly,r::Array{Int64, 1})
	result=algdata(a)
	for i in 1:length(r)
		for j in 1:r[i]
			result = d_aux(result, parent(a).derivation[i])
		end
	end
	return DiffPoly(parent(a), result)
end

function d(a::DiffIndet,r::Array{Int64, 1})
	return DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, r), parent(a).poly_ring))
end

#-
function d(a::Union{AbstractFloat, Integer, Rational},r::Array{Int64, 1})
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
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) * algdata(b))
	end
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * algdata(b))
end

function Base.:*(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(algdata(a) * DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
    return parent(a)(algdata(a) * DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:*(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) * DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:*(a::DiffIndet, b::Union{AbstractFloat, Integer, Rational})
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) * b)
	end
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) * b)
end

function Base.:*(a::Union{AbstractFloat, Integer, Rational}, b::DiffIndet)
	if length(parent(b).derivation)>1
		return parent(b)(a * DiffPoly(parent(b), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(b).poly_ring)))
	end
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
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) - b)
	end
    return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - b)
end

function Base.:-(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
	if length(parent(b).derivation)>1
		return parent(b)(a - DiffPoly(parent(b), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(b).poly_ring)))
	end
    return parent(b)(a - DiffPoly(parent(b), str_to_var(form_derivative(b.varname, 0), parent(b).poly_ring)))
end

function Base.:-(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) - DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

function Base.:-(a::DiffIndet, b::DifferentialRingElem)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(DiffPoly(parent(a), str_to_var(form_partial_derivative(a.varname, zeros(Int64,length(parent(a).derivation))), parent(a).poly_ring)) - b)
	end
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) - b)
end

function Base.:-(a::DifferentialRingElem, b::DiffIndet)
    check_parent(a, b)
	if length(parent(a).derivation)>1
		return parent(a)(a - DiffPoly(parent(a), str_to_var(form_partial_derivative(b.varname, zeros(Int64,length(parent(b).derivation))), parent(a).poly_ring)))
	end
	return parent(a)(a - DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end

#------------------------------------------------------------------------------

function Base.:-(a::DifferentialRingElem)
    return parent(a)(-algdata(a))
end

#-------------------------------------

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

#---------------------------------------------------------------------------------------

function Base.isless(a::DiffIndet, b::DiffIndet)
    check_parent(a, b)
	return parent(a)(DiffPoly(parent(a), str_to_var(form_derivative(a.varname, 0), parent(a).poly_ring)) < DiffPoly(parent(a), str_to_var(form_derivative(b.varname, 0), parent(a).poly_ring)))
end


#------------------------------------------------------------------------------

#end # module
