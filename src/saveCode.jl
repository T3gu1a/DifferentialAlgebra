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