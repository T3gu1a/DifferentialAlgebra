function Base.:>(a::DiffIndet, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return true
end

function Base.:<(a::DiffIndet, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return false
end

function Base.:<(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
    return true
end

function Base.:>(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DiffIndet)
    return false
end

function Base.:>(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return true
end

function Base.:<(a::DifferentialRingElem, b::Union{RingElem, AbstractFloat, Integer, Rational})
    return false
end

function Base.:<(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return true
end

function Base.:>(a::Union{RingElem, AbstractFloat, Integer, Rational}, b::DifferentialRingElem)
    return false
end

#the comparison returns a boolean integer (1 or 0)
#which the following code converts to true or false
function Binary_to_Bool(a::DiffPoly)
	return convert(Bool,parse(Int64,"$(a)"))
end