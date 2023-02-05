module SmallVectors

export SmallVector

"""
    SmallVector{N,T}(undef) -> A
    SmallVector{N,T}(vals) -> A
    SmallVector{N}(vals) -> A
    SmallVector(vals) -> A

yields a vector storing `N` elements of type `T`. The initial values `vals` may
be specified as a vector or as an `N`-tuple . If `vals` is a vector, the type
parameter `N` must be specified (for type-stability) and must be equal to
`length(vals)`. Initial values may also be specified as individual arguments.

Small vectors are similar to `MVector` in `StaticArrays`. However they can be
used as `Cstring` buffers.

"""
mutable struct SmallVector{N,T} <: AbstractVector{T}
    # NOTE: The structure must be mutable to have writable elements and element
    #       type must be plain data type (see Base.pointer(::SmallVector)).
    vals::NTuple{N,T}
    SmallVector{N,T}(::UndefInitializer) where {N,T} = new{check_N(N),check_T(T)}()
    SmallVector{N,T}(vals::NTuple{N,<:Any}) where {T,N} = new{check_N(N),check_T(T)}(vals)
end

# Constructors of small vectors from abstract vectors.
SmallVector{N}(vals::AbstractVector{T}) where {N,T} = SmallVector{N,T}(vals)
function SmallVector{N,T}(vals::AbstractVector) where {N,T}
    check_length(vals, check_N(N))
    A = SmallVector{N,T}(undef)
    off = firstindex(vals) - firstindex(A)
    @inbounds @simd for i in keys(A)
        A[i] = vals[i+off]
    end
    return A
end

# Constructors of small vectors from tuples.
SmallVector(vals::NTuple{N}) where {N} = SmallVector{N}(vals)
function SmallVector{N}(vals::NTuple{N}) where {N}
    T = promote_type(map(typeof, vals)...)
    return SmallVector{N,T}(vals)
end
SmallVector{N}(vals::Tuple) where {N} = bad_length(vals, N)
SmallVector{N,T}(vals::Tuple) where {N,T} = bad_length(vals, N)

# Constructors of small vectors from list of values.
SmallVector(a, b...) = SmallVector((a, b...))
SmallVector{N}(a, b...) where {N} = SmallVector{N}((a, b...))
SmallVector{N,T}(a, b...) where {N,T} = SmallVector{N,T}((a, b...))

Base.Tuple(A::SmallVector) = getfield(A, :vals)

# NOTE: pointer_from_objref requires mutable object and unsafe_convert to
#       Ptr{T} requires that T be are plain data type.
Base.unsafe_convert(::Type{Ptr{Nothing}}, A::SmallVector) =
    pointer_from_objref(A)
Base.unsafe_convert(::Type{Ptr{T}}, A::SmallVector{N,T}) where {N,T} =
    Ptr{T}(pointer_from_objref(A))
Base.pointer(A::SmallVector{N,T}) where {N,T} =
    Ptr{T}(pointer_from_objref(A))

Base.cconvert(::Type{Cstring}, A::SmallVector{N,UInt8}) where {N} = A
Base.unsafe_convert(::Type{Cstring}, A::SmallVector{N,UInt8}) where {N} =
    Cstring(Base.unsafe_convert(Ptr{Cchar}, A))

# Implement abstract array API.
Base.length(A::SmallVector{N,T}) where {N,T} = N
Base.IndexStyle(::Type{<:SmallVector}) = IndexLinear()
Base.size(A::SmallVector) = (length(A),)
Base.axes(A::SmallVector) = (keys(A),)
Base.keys(A::SmallVector) = Base.OneTo(length(A))
Base.firstindex(A::SmallVector) = 1
Base.lastindex(A::SmallVector) = length(A)

@inline function Base.getindex(A::SmallVector, i::Int)
    @boundscheck checkbounds(A, i)
    return GC.@preserve A unsafe_load(pointer(A), i)
end

@inline function Base.setindex!(A::SmallVector{N,T}, x, i::Int) where {N,T}
    @boundscheck checkbounds(A, i)
    GC.@preserve A unsafe_store!(pointer(A), convert(T, x), i)
    return A
end

check_N(N::Int) = N ≥ 0 ? N : error("type parameter N must be a nonnegative")
check_N(N) = error("type parameter N must be an `Int`")

check_T(T) = isbitstype(T) ? T : not_bits_type(T)
@noinline not_bits_type(T::Type) = error("type `$T` is not a plain data type")

check_length(A::Union{AbstractVector,Tuple},N::Integer) =
    length(A) == N || bad_length(A, N)
@noinline bad_length(A::AbstractArray, N::Integer) =
    error("array of $(length(A)) values while $N required")
@noinline bad_length(A::Tuple, N::Integer) =
    error("tuple of $(length(A)) values while $N required")

end # module
