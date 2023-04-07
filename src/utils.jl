# Yield error message.
function errmsg(status::Integer)
    buf = Array{UInt8}(undef, CFITSIO.FLEN_ERRMSG)
    ptr = pointer(buf)
    GC.@preserve buf CFITSIO.fits_get_errstatus(status, ptr)
    return unsafe_string(ptr)
end

function Base.show(io::IO, err::FitsError)
    print(io, "FitsError(")
    print(io, err.code)
    print(io, ')')
end

function Base.show(io::IO, ::MIME"text/plain", err::FitsError)
    show(io, err)
    print(io, ": \"")
    print(io, errmsg(err.code))
    print(io, "\"")
end

Base.showerror(io::IO, err::FitsError) = show(io, MIME("text/plain"), err)

# Yield whether a pointer is null.
isnull(ptr::Ptr) = (ptr === null(ptr))

# Yield a null pointer of a given type.
null(ptr::Ptr) = null(typeof(ptr))
null(::Type{Ptr{T}}) where {T} = Ptr{T}(0)

# Check argument whether it is a status code returned by one the CFITSIO
# library function or a pointer to a FITS file.
check(status::Ref{Status}) = check(status[])
check(status::Status) = status == 0 ? nothing : throw(FitsError(status))
check(ptr::Ptr{CFITSIO.fitsfile}) = isnull(ptr) ? bad_argument("FITS file has been closed") : ptr

bad_argument(str::ArgumentError.types[1]) = throw(ArgumentError(str))
@noinline bad_argument(args...) = bad_argument(string(args...))

@noinline invalid_property(obj, key::Symbol) =
    error("objects of type $(typeof(obj)) have no property `$key`")

@noinline readonly_property(obj, key::Symbol) =
    error("attempt to write read-only property `$key` for an object of type $(typeof(obj))")

@noinline throw_file_already_exists(args::AbstractString...) =
    error(file_already_exists(args...))

file_already_exists(path::AbstractString) =
    string("file \"", path, "\" already exists")

file_already_exists(path::AbstractString, usage::AbstractString) =
    string(file_already_exists(path), ", ", usage)

function Base.isequal(::FitsLogic, s1::AbstractString, s2::AbstractString)
    i1, last1 = firstindex(s1), lastindex(s1)
    i2, last2 = firstindex(s2), lastindex(s2)
    @inbounds while (i1 ≤ last1)&(i2 ≤ last2)
        uppercase(s1[i1]) == uppercase(s2[i2]) || return false
        i1 = nextind(s1, i1)
        i2 = nextind(s2, i2)
    end
    @inbounds while i1 ≤ last1
        isspace(s1[i1]) || return false
        i1 = nextind(s1, i1)
    end
    @inbounds while i2 ≤ last2
        isspace(s1[i2]) || return false
        i2 = nextind(s2, i2)
    end
    return true
end
Base.isequal(::FitsLogic, x) = y -> isequal(FitsLogic(), x, y)

"""
    EasyFITS.get_field_index(T::Type, name::Union{AbstractSring,Symbol}) -> i::Int

yields the index `i` of the field `name` of an object of type `T`. This
functions is used for introspection and meta-programming, and is not meant to
be fast.

"""
function get_field_index(::Type{T}, name::Symbol) where {T}
    for i in 1:fieldcount(T)
        fieldname(T, i) === name && return i
    end
    error("structure $T have no field `$name`")
end
get_field_index(::Type{T}, name::AbstractString) where {T} = get_field_index(T, Symbol(name))

"""
    EasyFITS.type_to_bitpix(T)

yields the code for FITS image pixels of type `T`. Argument can also be an
array instance or type.

Basic numeric types are recognized by this method which may be extended by
other packages to yield the CFITSIO codes equivalent to their own types. The
CFITSIO constants (to be prefixed by `EasyFITS.CFITSIO.`) and their
corresponding Julia types and standard BITPIX code are:

| CFITSIO Constant | Julia Type | `BITPIX` |
|:-----------------|:-----------|---------:|
| `BYTE_IMG`       | `UInt8`    |        8 |
| `SBYTE_IMG`      | `Int8`     |          |
| `SHORT_IMG`      | `Int16`    |       16 |
| `USHORT_IMG`     | `UInt16`   |          |
| `LONG_IMG`       | `Int32`    |       32 |
| `ULONG_IMG`      | `UInt32`   |          |
| `LONGLONG_IMG`   | `Int64`    |       64 |
| `ULONGLONG_IMG`  | `UInt64`   |          |
| `FLOAT_IMG`      | `Float32`  |      -32 |
| `DOUBLE_IMG`     | `Float64`  |      -64 |

Note that CFITSIO can read/write non-standard pixel types (those without a
`BITPIX` value above) by setting keywords `BSCALE` and `BZERO` with special
values as explicitely allowed by the FITS Standard (version 4).

"""
type_to_bitpix(arr::AbstractArray) = type_to_bitpix(typeof(arr))
type_to_bitpix(::Type{<:AbstractArray{T}}) where {T} = type_to_bitpix(T)
type_to_bitpix(::Type{Bool}) = type_to_bitpix(UInt8)

"""
    EasyFITS.type_from_bitpix(b) -> T

yields the Julia type `T` corresponding to FITS "BITPIX" `b`.

"""
type_from_bitpix(b::Integer) = type_from_bitpix(Int(b)::Int)

let expr = :(bad_argument("invalid BITPIX value"))
    for (type, bitpix) in ((UInt8,   CFITSIO.BYTE_IMG),
                           (Int8,    CFITSIO.SBYTE_IMG),
                           (UInt16,  CFITSIO.USHORT_IMG),
                           (Int16,   CFITSIO.SHORT_IMG),
                           (UInt32,  CFITSIO.ULONG_IMG),
                           (Int32,   CFITSIO.LONG_IMG),
                           (UInt64,  CFITSIO.ULONGLONG_IMG),
                           (Int64,   CFITSIO.LONGLONG_IMG),
                           (Float32, CFITSIO.FLOAT_IMG),
                           (Float64, CFITSIO.DOUBLE_IMG))
        bitpix = Int(bitpix)::Int
        @eval type_to_bitpix(::Type{$type}) = $bitpix
        expr = :(b == $bitpix ? $type : $expr)
    end
    @eval type_from_bitpix(b::Int) = $expr
end

# A few type assertions.
CFITSIO.LONGLONG === Clonglong || throw(AssertionError(
    "CFITSIO.LONGLONG = `$(CFITSIO.LONGLONG)` while Clonglong = `$Clonglong`"))
CFITSIO.ULONGLONG === Culonglong || throw(AssertionError(
    "CFITSIO.LONGLONG = `$(CFITSIO.ULONGLONG)` while Culonglong = `$Culonglong`"))

"""
    EasyFITS.type_to_code(T)

yields the CFITSIO type code for a keyword value or table cells of type `T`.
Argument can also be an array instance or type.

Basic numeric types and string types are recognized by this method which may be
extended by other packages to yield the CFITSIO codes equivalent to their own
types. The CFITSIO type constants (to be prefixed by `EasyFITS.CFITSIO.`) and
their corresponding C and Julia types are:

| CFITSIO Constant | C Types              | Julia Types        |
|:-----------------|:---------------------|:-------------------|
| `TLOGICAL`       | `char`               | `Cchar`, `Bool`    |
| `TBYTE`          | `unsigned char`      | `UInt8`, `Bool`    |
| `TSBYTE`         | `signed char`        | `Int8`             |
| `TUSHORT`        | `unsigned short`     | `Cushort`          |
| `TSHORT`         | `short`              | `Cshort`           |
| `TUINT`          | `unsigned int`       | `Cuint`            |
| `TINT`           | `int`                | `Cint`             |
| `TULONG`         | `unsigned long`      | `Culong`           |
| `TLONG`          | `long`               | `Clong`            |
| `TULONGLONG`     | `unsigned long long` | `Culonglong`       |
| `TLONGLONG`      | `long long`          | `Clonglong`        |
| `TFLOAT`         | `float`              | `Cfloat`           |
| `TDOUBLE`        | `double`             | `Cdouble`          |
| `TCOMPLEX`       | `float complex`      | `Complex{Cfloat}`  |
| `TDBLCOMPLEX`    | `double complex`     | `Complex{Cdouble}` |
| `TSTRING`        | `char*`              |                    |
| `TBIT`           |                      |                    |

"""
type_to_code(arr::AbstractArray) = type_to_code(typeof(arr))
type_to_code(::Type{<:AbstractArray{T}}) where {T} = type_to_code(T)

"""
    EasyFITS.pixeltype_to_code(T)

yields the CFITSIO type code for an image extension storing an array of
elements type `T`. Argument can also be an array instance or type.

"""
pixeltype_to_code(arr::AbstractArray) = pixeltype_to_code(typeof(arr))
pixeltype_to_code(::Type{<:AbstractArray{T}}) where {T} = pixeltype_to_code(T)
pixeltype_to_code(::Type{T}) where {T} = type_to_code(T)
pixeltype_to_code(::Type{Bool}) =
    # NOTE: For FITS image extension Booleans must be considered as bytes in
    # CFITSIO library.
    type_to_code(UInt8)


"""
    EasyFITS.type_from_code(c) -> T

yields the Julia type `T` corresponding to CFITSIO type code `c`.

"""
type_from_code(c::Integer) = type_from_code(Int(c)::Int)

# Provide default value, if symbol constant is not defined in C FITSIO library
# (e.g. MS-Windows). Otherwise, check that assumed and actual values agree.
# NOTE: This "hack" is needed to cope with some symbols under MS-Windows. If
# the assumed constant is wrong, the call will clash when tested under, say,
# Linux.
function _const(key::Symbol, val)
    if isdefined(CFITSIO, key)
        cval = getfield(CFITSIO, key)
        val == cval || error("CFITSIO constant $sym is $cval, not $val")
    end
    return val
end

# Data type code as assumed by the CFITSIO library. NOTE: Sets are used to
# avoid duplicates.
let types = Set{DataType}(),
    codes = Set{Int}(),
    expr = :(bad_argument("invalid type code"))
    for (type, code) in ((String,           CFITSIO.TSTRING),
                         (Bool,             CFITSIO.TLOGICAL),
                         (Int8,             CFITSIO.TSBYTE),
                         (UInt8,            _const(:TBYTE, 11)),
                         (Cshort,           CFITSIO.TSHORT),
                         (Cushort,          CFITSIO.TUSHORT),
                         (Cint,             CFITSIO.TINT),
                         (Cuint,            CFITSIO.TUINT),
                         (Clong,            CFITSIO.TLONG),
                         (Culong,           CFITSIO.TULONG),
                         (Clonglong,        CFITSIO.TLONGLONG),
                         (Culonglong,       CFITSIO.TULONGLONG),
                         (Cfloat,           CFITSIO.TFLOAT),
                         (Cdouble,          CFITSIO.TDOUBLE),
                         (Complex{Cfloat},  CFITSIO.TCOMPLEX),
                         (Complex{Cdouble}, CFITSIO.TDBLCOMPLEX),
                         (Bit,              CFITSIO.TBIT))
        code = Int(code)::Int
        if ! (type ∈ types)
            push!(types, type)
            if type === String
                @eval type_to_code(::Type{<:AbstractString}) = $code
            else
                @eval type_to_code(::Type{$type}) = $code
            end
        end
        if ! (code ∈ codes)
            push!(codes, code)
            expr = :(c == $code ? $type : $expr)
        end
    end
    @eval type_from_code(c::Int) = $expr
end

"""
    EasyFITS.type_to_letter(T)

yields the letter of the `TFORMn` keyword representing table cells of type `T`
in FITS. Argument can also be an array instance or type.

"""
type_to_letter(arr::AbstractArray) = type_to_letter(typeof(arr))
type_to_letter(::Type{<:AbstractArray{T}}) where {T} = type_to_letter(T)

"""
    EasyFITS.type_from_letter(c) -> T

yields the Julia type `T` corresponding to CFITSIO column type letter `c` as
assumed for the `TFORMn` keywords.

"""
type_from_letter(c::Integer) = type_from_letter(Char(c)::Char)

# FITS table extension column types as specified by the TFORMn keywords.
let expr = :(bad_argument("invalid TFORM type letter"))
    for (type, letter) in ((String,     'A'),
                           (Bool,       'L'),
                           (UInt8,      'B'),
                           (Int16,      'I'),
                           (Int32,      'J'),
                           (Int64,      'K'),
                           (Float32,    'E'),
                           (Float64,    'D'),
                           (ComplexF32, 'C'),
                           (ComplexF64, 'M'),
                           (Bit,        'X'),
                           # The followings are specific to CFITSIO.
                           (Int8,       'S'),
                           (UInt16,     'U'),
                           (UInt32,     'V'),
                           (UInt64,     'W'))
        if type === String
            @eval type_to_letter(::Type{<:AbstractString}) = $letter
        else
            @eval type_to_letter(::Type{$type}) = $letter
        end
        expr = :(c == $letter ? $type : $expr)
    end
    @eval type_from_letter(c::Char) = $expr
end

# Named tuples to associate short and long suffixes to Julia types.
const TYPE_SUFFIXES =
    (Int8            => (:sb, :sbyt),   UInt8            => (:b,   :byt),
     Cshort          => (:i,  :sht),    Cushort          => (:ui,  :usht),
     Cint            => (:k,  :int),    Cuint            => (:uk,  :uint),
     Clong           => (:j,  :lng),    Culong           => (:uj,  :ulng),
     Clonglong       => (:jj, :lnglng), Culonglong       => (:ujj, :ulnglng),
     Cfloat          => (:e,  :flt),    Cdouble          => (:d,   :dbl),
     Complex{Cfloat} => (:c,  :cmp),    Complex{Cdouble} => (:m,   :dblcmp),
     Bool            => (:l,  :log),    AbstractString   => (:s,   :str),
     Bit             => (:x,  :bit),    Nothing          => (:u,   :null),
     Missing         => (:u,  :null),   UndefInitializer => (:u,   :null))

let types = Set{DataType}()
    for (T, sfx) in TYPE_SUFFIXES
        T ∈ types && continue
        push!(types, T)
        S = isconcretetype(T) ? T : :(<:$T)
        @eval begin
            short_suffix(::Type{$S}) = $(String(sfx[1]))
            long_suffix( ::Type{$S}) = $(String(sfx[2]))
        end
    end
    @eval const NUMERIC_TYPES = $((filter(T -> T <: Number, types)...,))
    @eval const INTEGER_TYPES = $((filter(T -> T <: Integer, types)...,))
    @eval const   FLOAT_TYPES = $((filter(T -> T <: AbstractFloat, types)...,))
    @eval const COMPLEX_TYPES = $((filter(T -> T <: Complex, types)...,))
end

# Numeric types supported by the library for reading/writing.
const NumericTypes = Union{NUMERIC_TYPES...,}
const IntegerTypes = Union{INTEGER_TYPES...,}
const   FloatTypes = Union{  FLOAT_TYPES...,}
const ComplexTypes = Union{COMPLEX_TYPES...,}

"""
    EasyFITS.cfunc(pfx::Union{AbstractString,Symbol}, T::Type) -> sym

yields the symbolic name of the function in the CFITSIO library whose prefix is
`pfx` and whose suffix is deduced from the type `T`. Long/short function names
are supported and automatically detected depending on whether `pfx` ends with
an underscore character.

"""
cfunc(func::Symbol, ::Type{T}) where {T} = cfunc(String(func), T)
cfunc(func::AbstractString, ::Type{T}) where {T} =
    Symbol(func, endswith(func, '_') ? long_suffix(T) : short_suffix(T))

"""
    EasyFITS.ctype(T) -> T′

yields the C type equivalent to `T` in CFITSIO library. This is mostly for
Booleans which are treated as C `char` in the library, the other basic
numerical types being unchanged. An error is thrown if the returned type has a
different size than that of `T`. This is to make sure that Julia arrays with
elements of type `T` can safely be used to store values of type `T′`.

!!! warning
    This only applies to element type of arrays. For boolean scalars, a `Cint`
    is the correct type.

"""
function ctype(::Type{Bool})
    T = Cchar
    sizeof(T) == sizeof(Bool) && return T
    throw(AssertionError("in CFITSIO, logical values are of type `$T` but sizeof(Bool) != sizeof($T)"))
end

function ctype(::Type{T}) where {T<:Union{Int8,UInt8,Int16,UInt16,
                                          Int32,UInt32,Int64,UInt64,
                                          Float32,ComplexF32,
                                          Float64,ComplexF64}}
    return T
end

"""
    EasyFITS.cpointer(arr::AbstractArray) -> ptr::Ptr{ctype(eltype(arr))}

yields a pointer to the elements of array `arr` that can be used in calls to
functions of the CFITSIO library. Compared to `Ptr{Cvoid}(pointer(arr))`, this
function yields a typed pointer which prevents using arguments of the wrong
type.

!!! warning
    Do not forget to use `GC.@protect arr ...` to avoid `arr` being garbage
    collected while its address is in use.

"""
cpointer(arr::AbstractArray) = convert(Ptr{ctype(eltype(arr))}, pointer(arr))

"""
    EasyFITS.fix_booleans!(arr::AbstractArray{Bool}) -> arr

makes sure that the values in array of booleans `arr` read by a function of the
CFITSIO library are only `true` or `false`.

"""
function fix_booleans!(arr::AbstractArray{Bool})
    T = ctype(Bool)
    @inbounds @simd for i in eachindex(arr)
        arr[i] = !iszero(reinterpret(T, arr[i]))
    end
    return arr
end

"""
    EasyFITS.map_recursively(f, r, a1, a2; missing=missing)

yields the result of applying function `f` recursively to update result `r` for
each pair of values of tuples `a1` and `a2`. The recursion stops when the
longest of `a1` and `a2` is exhausted. In this process, missing values are
specified by the keyword `missing`. This is equivalent to:

    for i in 1:max(length(a1), length(a2))
        r = f(r, (i ≤ length(a1) ? a1[i] : missing),
                 (i ≤ length(a2) ? a2[i] : missing))
    end
    return r

except that in-lining is used to unroll the loop.

"""
@inline map_recursively(f::Function, r, a1::Tuple{}, a2::Tuple{}; missing=missing) = r
@inline map_recursively(f::Function, r, a1::Tuple, a2::Tuple{}; missing=missing) =
    map_recursively(f, f(r, first(a1), missing), Base.tail(a1), (); missing=missing)
@inline map_recursively(f::Function, r, a1::Tuple{}, a2::Tuple; missing=missing) =
    map_recursively(f, f(r, missing, first(a2)), (), Base.tail(a2); missing=missing)
@inline map_recursively(f::Function, r, a1::Tuple, a2::Tuple; missing=missing) =
    map_recursively(f, f(r, first(a1), first(a2)), Base.tail(a1), Base.tail(a2); missing)

"""
    EasyFITS.subarray_params(dims, inds) -> (d, f, s, l)

yields resulting dimensions `d`, first indices `f`, steps `s`, and last indices
`l` when sub-indexing an array of size `dims` with sub-indices `inds`. Tuples
`f`, `s`, and `l` have the same number of values as `dims`; `d` may have fewer
values.

"""
subarray_params(dims::Dims, inds::SubArrayIndices) =
    map_recursively(_update_subarray_params, ((),(),(),()), dims, inds)

# More dimensions than sub-indices is forbidden.
_update_subarray_params(p, d::Integer, i::Missing) = too_few_subindices()
@noinline too_few_subindices() = throw(DimensionMismatch("too few sub-indices specified"))

# Update parameters for sub-indexing a given dimension (`d` is an integer) or
# beyond the number of dimensions (`d` is `missing`).
_update_subarray_params(p, d::Union{Integer,Missing}, i::Integer) = begin
    verify_subindex(d, i)
    if d === missing
        return p
    else
        return (p[1],         # no additional sub-array dimension
                (p[2]..., i), # first index
                (p[3]..., 1), # step
                (p[4]..., i)) # last index
    end
end

_update_subarray_params(p, d::Union{Integer,Missing}, ::Colon) =
    ((p[1]..., (d === missing ? 1 : d)), # sub-array dimension
     (p[2]..., 1),                       # first index
     (p[3]..., 1),                       # step
     (p[4]..., (d === missing ? 1 : d))) # last index

_update_subarray_params(p, d::Union{Integer,Missing}, r::SubArrayIndex) = begin
    verify_subindex(d, r)
    return ((p[1]..., length(r)), # sub-array dimension
            (p[2]..., first(r)),  # first index
            (p[3]..., step(r)),   # step
            (p[4]..., last(r)))   # last index
end

is_valid_subindex(d::Missing, i::Integer) = isone(i)
is_valid_subindex(d::Integer, i::Integer) = (i ≥ one(i)) & (i ≤ d)

is_valid_subindex(d::Union{Integer,Missing}, ::Colon) = true

is_valid_subindex(d::Missing, r::IndexRange) = isone(first(r)) & isone(length(r))
is_valid_subindex(d::Integer, r::IndexRange) = begin
    i_min, i_max = extrema(r)
    return (i_min ≤ i_max) & (i_min ≥ one(i_min)) & (i_max ≤ d)
end

function verify_subindex(d::Union{Integer,Missing}, i::SubArrayIndex)
    is_valid_subindex(d, i) || out_of_bounds_subindex()
    return nothing
end
@noinline out_of_bounds_subindex() = throw(ArgumentError("out of bound sub-index"))

"""
    EasyFITS.new_array(T, dims...) -> arr

yields a new array with element type `T` and dimensions `dims...` which may be
specified as for the `Array` constructor or as a vector of integers. In this
latter case, the result is not type-stable.

If the number of dimensions, say `N`, is known, a type-stable result is
returned by:

    EasyFITS.new_array(T, Val(N), dims...) -> arr

"""
new_array(::Type{T}, dims::Integer...) where {T} = Array{T}(undef, dims)
new_array(::Type{T}, dims::NTuple{N,Integer}) where {T,N} = Array{T,N}(undef, dims)
new_array(::Type{T}, dims::AbstractVector{<:Integer}) where {T} =
    new_array(T, Val(length(dims)), dims)
new_array(::Type{T}, n::Val{N}, dims::Integer...) where {T,N} = new_array(T, n, dims)
@generated function new_array(::Type{T}, ::Val{N}, dims) where {T,N}
    d = [:(dims[$k]) for k in 1:N]
    return quote
        length(dims) == N || throw(DimensionMismatch("incompatible number of dimensions"))
        return Array{$T,$N}(undef, $(d...))
    end
end

"""
    EasyFITS.dense_array(arr)

yields a dense array, that is an array that stores its elements contiguously,
with same dimensions and values as `arr`. If `arr` is already a dense array, it
is returned.

"""
dense_array(arr::DenseArray) = arr
dense_array(arr::AbstractArray{T,N}) where {T,N} = convert(Array{T,N}, arr)

function to_string!(buf::Vector{UInt8})
    @inbounds for i in 1:length(buf)
        if iszero(buf[i])
            resize!(buf, i - 1)
            break
        end
    end
    return String(buf)
end

"""
    EasyFITS.string_length(str)

yields the length of the string `str` not counting non-significant trailing
spaces. As assumed in the CFITSIO library and by FITS standard, if there are
only spaces in `str`, the first space is considered to be significant. This is
intended to distinguish null (empty) and non-null strings.

"""
function string_length(str::AbstractString)
    num = 0 # number of characters
    len = 0 # number of character without non-significant trailing spaces
    @inbounds for c in str
        num += 1
        if !is_space(c)
            len = num
        end
    end
    return ((num > 0)&(len < 1)) ? 1 : len
end

function maximum_length(A::AbstractArray{<:AbstractString})
    maxlen = 0
    @inbounds for str in A
        maxlen = max(maxlen, string_length(str))
    end
    return maxlen
end
