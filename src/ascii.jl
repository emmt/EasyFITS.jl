# Implement ASCII encoding for FITS cards.

@inline function Base.thisind(s::FitsCardPart{:ascii}, i::Int)
    i == 0 && return 0
    n = ncodeunits(s)
    i == n + 1 && return i
    @boundscheck between(i, 1, n) || throw(BoundsError(s, i))
    return i
end

@inline function Base.nextind(s::FitsCardPart{:ascii}, i::Int)
    i == 0 && return 1
    n = ncodeunits(s)
    @boundscheck between(i, 1, n) || throw(BoundsError(s, i))
    return i + 1
end

@inline function Base.iterate(s::FitsCardPart{:ascii}, i::Int = firstindex(s))
    (i % UInt) - 1 < ncodeunits(s) || return nothing
    b = @inbounds codeunit(s, i)
    return Char(b), i + 1
end

@propagate_inbounds Base.getindex(s::FitsCardPart{:ascii}, i::Int) = Char(codeunit(s, i))

Base.length(A::FitsCardPart{:ascii}) = ncodeunits(A)
function Base.length(A::FitsCardPart{:ascii}, i::Int, j::Int)
    @boundscheck begin
        np1 = ncodeunits(s) + 1
        0 < i ≤ np1 || throw(BoundsError(s, i))
        0 ≤ j < np1 || throw(BoundsError(s, j))
    end
    return ifelse(j < i, 0, j - i + 1)
end

function Base.isvalid(s::FitsCardPart{:ascii}, i::Int)
    @boundscheck return checkbounds(Bool, s, i)
    return true
end
