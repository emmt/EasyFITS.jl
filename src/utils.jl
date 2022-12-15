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
check(status::Status) = (status == 0 || throw(FitsError(status)))
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
get_field_index(T::Type, name::AbstractString) = get_field_index(T, Symbol(name))

"""
    EasyFITS.type_to_bitpix(T)

yields the code for FITS image pixels of type `T`. Argument can also be an
array instance or type.

"""
type_to_bitpix(arr::AbstractArray) = type_to_bitpix(typeof(arr))
type_to_bitpix(::Type{AbstractArray{T}}) where {T} = type_to_bitpix(T)

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

"""
    EasyFITS.type_to_code(T)

yields the CFITSIO type code for a keyword value or table cells of type `T`.
Argument can also be an array instance or type.

"""
type_to_code(arr::AbstractArray) = type_to_code(typeof(arr))
type_to_code(::Type{AbstractArray{T}}) where {T} = type_to_code(T)

"""
    EasyFITS.type_from_code(c) -> T

yields the Julia type `T` corresponding to CFITSIO type code `c`.

"""
type_from_code(c::Integer) = type_from_code(Int(c)::Int)

# Data type code as assumed by the CFITSIO library. NOTE: Sets are used to
# avoid duplicates.
let types = Set{DataType}(),
    codes = Set{Int}(),
    expr = :(bad_argument("invalid type code"))
    for (type, code) in ((String,           CFITSIO.TSTRING),
                         (Bool,             CFITSIO.TLOGICAL),
                         (Int8,             CFITSIO.TSBYTE),
                         (UInt8,            CFITSIO.TBYTE),
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
type_to_letter(::Type{AbstractArray{T}}) where {T} = type_to_letter(T)

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

# Named tuples to associate suffixes to types.
const TYPE_SHORT_SUFFIXES = (sb  = Int8,            b   = UInt8,
                             i   = Cshort,          ui  = Cushort,
                             k   = Cint,            uk  = Cuint,
                             j   = Clong,           uj  = Culong,
                             jj  = Clonglong,       ujj = Culonglong,
                             e   = Cfloat,          d   = Cdouble,
                             c   = Complex{Cfloat}, m   = Complex{Cdouble},
                             l   = Bool,            s   = String,
                             x   = Bit)

const TYPE_LONG_SUFFIXES = (sbyt    = Int8,            byt      = UInt8,
                            sht     = Cshort,          usht     = Cushort,
                            int     = Cint,            uint     = Cuint,
                            lng     = Clong,           ulng     = Culong,
                            lnglng  = Clonglong,       ulnglng  = Culonglong,
                            flt     = Cfloat,          dbl      = Cdouble,
                            cmp     = Complex{Cfloat}, dblcmp   = Complex{Cdouble},
                            log     = Bool,            str      = String,
                            bit     = Bit)

let S = Set{DataType}()
    for (sym, T) in pairs(TYPE_SHORT_SUFFIXES)
        T ∈ S && continue
        push!(S, T)
        if T === String
            @eval short_suffix(::Type{<:AbstractString}) = $(String(sym))
        else
            @eval short_suffix(::Type{$T}) = $(String(sym))
        end
    end
end

let S = Set{DataType}()
    for (sym, T) in pairs(TYPE_LONG_SUFFIXES)
        T ∈ S && continue
        push!(S, T)
        if T === String
            @eval long_suffix(::Type{<:AbstractString}) = $(String(sym))
        else
            @eval long_suffix(::Type{$T}) = $(String(sym))
        end
    end
end

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
        length(dims) == N || bad_argument("incompatible number of dimensions")
        return Array{$T,$N}(undef, $(d...))
    end
end

function to_string!(buf::Vector{UInt8})
    @inbounds for i in 1:length(buf)
        if iszero(buf[i])
            resize!(buf, i - 1)
            break
        end
    end
    return String(buf)
end

to_type(::Type{T}, x::T) where {T} = x
to_type(::Type{T}, x) where {T} = convert(T, x)::T
