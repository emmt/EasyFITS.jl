#
# EasyFITS.jl -
#
# Utilities to simplify reading/writing FITS files.
#
#-------------------------------------------------------------------------------
#
# This file if part of EasyFITS software (https://github.com/emmt/EasyFITS.jl)
# licensed under the MIT license.
#
# Copyright (C) 2018-2019, Éric Thiébaut.
#

module EasyFITS

export
    FitsArray,
    FitsComment,
    FitsFile,
    FitsHDU,
    FitsHeader,
    FitsIO,
    FitsImage,
    FitsImageHDU,
    FitsTableHDU,
    createfits!,     # FIXME: deprecated
    createfits,      # FIXME: deprecated
    exists,
    getfitskey,      # FIXME: deprecated
    openfits,        # FIXME: deprecated
    readfits,        # FIXME: deprecated
    setfitskey!,     # FIXME: deprecated
    tryreadfitskey,  # FIXME: deprecated
    tryreadfitskeys, # FIXME: deprecated
    write!,
    writefits!,      # FIXME: deprecated
    writefits        # FIXME: deprecated

using FITSIO
using FITSIO.Libcfitsio
using FITSIO.Libcfitsio: libcfitsio

using Base: getfield, elsize, tail, OneTo, throw_boundserror, @propagate_inbounds
import Base: get, getindex, setindex!, keys, haskey, getkey,
    read, write, open, close

#
# Predefine all types.
#

struct FitsIO
    fh::FITS
end

struct FitsHDU{T<:HDU}
    hdu::T
end

const FitsImageHDU = FitsHDU{ImageHDU}
const FitsTableHDU = FitsHDU{TableHDU}

struct FitsHeader
    hdr::FITSHeader
end

struct FitsImage{T,N} <: DenseArray{T,N}
    arr::Array{T,N}
    hdr::FitsHeader
end

# Singleton type for writing a FITS file.
struct FitsFile end

# Abstract type for retrieving the array part of an Image HDU in a FITS file.
abstract type FitsArray{T,N} <: AbstractArray{T,N} end

# Singleton type for retrieving a keyword comment.
struct FitsComment end

# Singleton type for marking a missing keyword.
struct Missing end

# Annotated objects have a FITS header and implements indexation by keywords,
# they can also be directly read from a FITS file.
const Annotated = Union{FitsImage,FitsHeader}

# Readable objects are those which can be directly read from a FITS file.
const Readable = Union{Annotated,FitsArray}

# Allowed types for keywords values.
const KeywordValues = Union{Bool,Int,Float64,String}

# Valid types to specify a FITS extension.
const Extension = Union{Integer,AbstractString}

"""

```julia
exists(path) -> boolean
```

yields whether file `path` exists.  Argument can be a file name or an instance
of `Base.Filesystem.StatStruct` as returned by the `stat` method.

"""
exists(st::Base.Filesystem.StatStruct) = (st.nlink != 0)
exists(path) = exists(stat(path))

"""

```julia
EasyFITS.hduname(T) -> name, vers=1
```

yields the name and revision number of FITS Header Data Unit for storing data
of type `T`.  Other packages are encouraged to extend this method with their
types.

""" hduname
hduname(::T) where {T} = hduname(T)
@noinline hduname(::Type{T}) where {T} =
    error(string("method `EasyFITS.hduname` has not been specialized ",
                 "for type `", T, "`"))

# Make sure that a (name, vers) tuple is returned.
fixhduname(x) = _fixhduname(hduname(x))
_fixhduname(str::AbstractString) = _fixhduname(str, 1)
_fixhduname(arg::Tuple{AbstractString,Integer}) = _fixhduname(arg[1], Int(arg[2]))
_fixhduname(arg::Tuple{AbstractString,Int}) = _fixhduname(String(arg[1]), arg[2])
_fixhduname(arg::Tuple{String,Int}) = arg

"""

```julia
FitsIO(path, mode="r") -> io
```

opens of creates a FITS file whose name is `path` and returns an instance of
`FitsIO` which is an opaque handle to a FITS file.  Argument `mode` can be:

- `"r"` or `"r+"` to read or append to the contents of and existing FITS file.
  If file `path` does not exist but `path` does not end with the `".gz"`
  extension and `"\$path.gz"` does exist, then the compressed file
  `"\$path.gz"` is open instead.

- `"w"` to create a new FITS file named `path` that must not already exist.  An
  error is thrown if the file already exists.

- `"w!"` to open FITS file named `path` for writing.  If the file already
  exists, it is (silently) overwritten.

The do-block syntax is supported to automatically close the FITS file:

```julia
FitsIO(filename, mode="r") do io
    # use FITS handle io
    ...
end
```

An instance of `FitsIO` is a collection of *Header Data Units* (HDU) and
implements indexation and iteration.  Assuming `io` is a `FitsIO` object, then:

- `io[i]` yields the `i`-th HDU.

- `io[name]` or `io[name,vers]` yields the HDU whose `EXTNAME` (or `HDUNAME`)
  keyword is equal to `name` (a string) and, optionally, whose `EXTVER` (or
  `HDUVER`) keyword is equal to `vers` (an integer).

- You can do `for hdu in io; ...; end` to iterate through all HDU's of `io`.

- Methods `findfirst(p,io)`, `findlast(p,io)`, `findnext(p,io,i)` and
  `findprev(p,io,i)` can be used on `FitsIO` object `io` to search for a
  specific HDU.  These methods test each HDU (starting at initial index `i` for
  `findnext` and `findprev`) with the predicate function `p` (called with a
  `FitsHDU` argument) and return the index (an integer) of the first HDU for
  which the predicate yields `true` or `nothing` if this never occurs.

Also see: [`EasyFITS.find`](@ref).

""" FitsIO

FitsIO(path::AbstractString, mode::AbstractString="r") = begin
    local handle::FITS
    if mode == "r" || mode == "r+"
        handle = FITS((isfile(path) || endswith(path, ".gz") ||
                       !isfile(path*".gz") ? path : path*".gz"), mode)
    else
        if mode == "w"
            if exists(path)
                throw_file_already_exists(path, "use mode \"w!\" to overwrite")
            end
        elseif mode != "w!"
            throw(ArgumentError("mode must be \"r\", \"r+\", \"w\" or \"w!\""))
        end
        handle = FITS(path, "w")
    end
    return FitsIO(handle)
end

FitsIO(func::Function, args::AbstractString...) = begin
    io = FitsIO(args...)
    try
        func(io)
    finally
        close(io)
    end
end

Base.get(::Type{FITS}, obj::FitsIO) = getfield(obj, :fh)
Base.get(::Type{FitsIO}, obj::FitsIO) = obj

Base.convert(::Type{FitsIO}, obj::FITS) = FitsIO(obj)
Base.convert(::Type{FITS}, obj::FitsIO) = get(FITS, obj)

Base.close(io::FitsIO) = close(get(FITS, io))
Base.length(io::FitsIO) = length(get(FITS, io))
Base.getindex(io::FitsIO, args...) = FitsHDU(getindex(get(FITS, io), args...))
Base.iterate(io::FitsIO, args...) = begin
    result = iterate(get(FITS, io), args...)
    if result === nothing
        return nothing
    end
    return (FitsHDU(result[1]), result[2])
end
Base.lastindex(io::FitsIO) = lastindex(get(FITS, io))
Base.show(io::IO, obj::FitsIO) = show(io, get(FITS, obj))

Base.findfirst(pred::Function, io::FitsIO) = findnext(pred, io, 1)
Base.findlast(pred::Function, io::FitsIO) = findprev(pred, io, length(io))
Base.findnext(pred::Function, io::FitsIO, i::Integer) = findnext(pred, io, Int(i))
Base.findnext(pred::Function, io::FitsIO, i::Int=1) = find(pred, io, i:length(io))
Base.findprev(pred::Function, io::FitsIO, i::Integer) = findprev(pred, io, Int(i))
Base.findprev(pred::Function, io::FitsIO, i::Int=1) = find(pred, io, i:-1:1)

"""

```julia
find(pred, io, I=1:length(io))
```

yields the HDU first index `i ∈ I` of FITS instance `io` for which the
predicate `pred(io[i])` returns `true`, or [`nothing`](@ref) if `pred(io[i])`
returns `false` for all `i ∈ I`.

For instance:

```julia
i = find(hdu -> get(String,hdu,"EXTNAME","") == "CALIBRATION", io)
if i === nothing
    # not found
    ...
else
    # found at index i
    hdu = io[i]
end
```

"""
function find(pred::Function, io::FitsIO,
              I::AbstractVector{Int} = 1:length(io)) :: Union{Nothing,Int}
    for i in I
        if pred(io[i])
            return i
        end
    end
    return nothing
end

"""

An instance of type `FitsHDU{T}` is an opaque handle to a FITS *Header Data
Unit* (HDU).  Parameter `T` is to distinguish between the different kinds of
HDU (*Image* or *Binary Table*).

""" FitsHDU

Base.get(::Type{HDU}, obj::FitsHDU) = getfield(obj, :hdu)
Base.get(::Type{ImageHDU}, obj::FitsImageHDU) = get(HDU, obj)
Base.get(::Type{TableHDU}, obj::FitsTableHDU) = get(HDU, obj)
Base.get(::Type{T}, obj::T) where {T<:FitsHDU} = obj

Base.convert(::Type{FitsHDU}, obj::HDU) = FitsHDU(obj)
Base.convert(::Type{FitsImageHDU}, obj::ImageHDU) = FitsHDU(obj)
Base.convert(::Type{FitsTableHDU}, obj::TableHDU) = FitsHDU(obj)

Base.convert(::Type{HDU}, obj::FitsHDU) = get(HDU, obj)
Base.convert(::Type{ImageHDU}, obj::FitsImageHDU) = get(HDU, obj)
Base.convert(::Type{TableHDU}, obj::FitsTableHDU) = get(HDU, obj)

Base.show(io::IO, obj::FitsHDU) = show(io, get(HDU, obj))
Base.size(obj::FitsImageHDU, args...) = size(get(HDU, obj), args...)
Base.ndims(obj::FitsImageHDU) = ndims(get(HDU, obj))
Base.length(obj::FitsImageHDU) = length(get(HDU, obj))
Base.lastindex(obj::FitsImageHDU) = lastindex(get(HDU, obj))

"""

An object of type `FitsHeader` stores the header part of a FITS *Header Data
Unit* (HDU).  There are several ways to build an instance of `FitsHeader`.

```julia
FitsHeader()
```

yields an empty instance of `FitsHeader`.

```julia
FitsHeader(io, ext=1)
```

yields an instance of `FitsHeader` built over FITS header read from extension
`ext` in FITS file `io`.

```julia
FitsHeader(; key1=val1, key2=val2, ...)
```

yields an instance of `FitsHeader` whose contents is set by keywords (values
can be a tuple with a value and a comment).  This is can also be done by
specifying key-value pairs:

```julia
FitsHeader("key1" => val1, "key2" => val2, ...)
```

To avoid ambiguities the two styles cannot be mixed.

An object of type `FitsHeader`, say `obj`, implements indexation by keywords,
as with `obj[key]`, and `obj.key` syntax.


## Implementation

Type `FitsHeader` is a simple wrapper over `FITSIO.FITSHeader` (beware of the
different spellings) to implement indexation by keywords and `obj.key` syntax.
We have to define a new type for that because implementing the extended
interface over the existing `FITSIO.FITSHeader` would constitute a *type
piracy*.  An instance of `FitsHeader` can be simply built over an object `hdr`
of type `FITSIO.FITSHeader` by calling:

```julia
FitsHeader(hdr)
```

""" FitsHeader

FitsHeader(hdu::FitsHDU) = FitsHeader(read_header(get(HDU, hdu)))

FitsHeader(io::FitsIO, ext::Extension = 1) = FitsHeader(io[ext])

function FitsHeader(; kwds...)
    hdr = FitsHeader(FITSHeader(String[], [], String[]))
    for (key, val) in kwds
        setproperty!(hdr, key, val)
    end
    return hdr
end

function FitsHeader(args::Pair{<:AbstractString}...)
    hdr = FitsHeader(FITSHeader(String[], [], String[]))
    for (key, val) in args
        hdr[key] = val
    end
    return hdr
end

Base.get(::Type{Array}, obj::FitsHeader) = nothing # FIXME: ???
Base.get(::Type{FitsHeader}, obj::FitsHeader) = obj
Base.get(::Type{FITSHeader}, obj::FitsHeader) = getfield(obj, :hdr)

Base.convert(::Type{FitsHeader}, hdr::FITSHeader) = FitsHeader(hdr)
Base.convert(::Type{FITSHeader}, hdr::FitsHeader) = get(FITSHeader, hdr)

Base.length(obj::FitsHeader) = length(get(FITSHeader, obj))
Base.show(io::IO, obj::FitsHeader) = show(io, get(FITSHeader, obj))

"""

Type `FitsComment` is a singleton used as a marker to indicate that the comment
of a FITS keyword is to be returned by the `get` method.

""" FitsComment

"""

```julia
FitsImage(arr, hdr)
```

yields a FITS Image that behaves like a Julia array when indexed by integers of
Cartesian indices and like a dictionary of FITS keywords when indexed by
strings.  Argument `arr` specifies the array contents of the object and `hdr`
specifies the keywords both are shared by the returned object.

The initial contents of the header may be specified by keywords.  For instance:

```julia
FitsImage(arr; KEY1 = val1, KEY2 = (val2, com2), ...)
```

The embedded array may be creted by the constructor:

```julia
FitsImage{T}(init, dims...; kwds...)
```

yields a FITS Image with an empty header and array-like data of dimensions
`dims` and with elements of type `T` initially set by `init`.  Specifying the
number of dimensions is also supported:

```julia
FitsImage{T,N}(init, dims...; kwds...)
```

""" FitsImage

FitsImage(arr::AbstractArray; kwds...) = FitsImage(arr, FitsHeader(; kwds...))
FitsImage(arr::AbstractArray{T,N}, hdr::FitsHeader) where {T,N} =
    FitsImage{T,N}(convert(Array{T,N}, arr), hdr)
FitsImage(arr::AbstractArray, hdr::FITSHeader) =
    FitsImage(arr, FitsHeader(hdr))
FitsImage{T}(init, dims::Integer...; kwds...) where {T} =
    FitsImage{T}(init, dims; kwds...)
FitsImage{T}(init, dims::NTuple{N,Integer}; kwds...) where {T,N} =
    FitsImage{T,N}(init, dims; kwds...)
FitsImage{T,N}(init, dims::Integer...; kwds...) where {T,N} =
    FitsImage{T,N}(init, dims; kwds...)
FitsImage{T,N}(init, dims::NTuple{N,Integer}; kwds...) where {T,N} =
    FitsImage{T,N}(Array{T,N}(init, dims), FitsHeader(; kwds...))

Base.get(::Type{Array}, obj::FitsImage) = getfield(obj, :arr)
Base.get(::Type{FitsHeader}, obj::FitsImage) = getfield(obj, :hdr)
Base.get(::Type{FITSHeader}, obj::FitsImage) = get(FITSHeader, get(FitsHeader, obj))
Base.get(::Type{FitsImage}, obj::FitsImage) = obj

Base.convert(::Type{FitsImage}, arr::AbstractArray) = FitsImage(arr)
Base.convert(::Type{FitsHeader}, obj::FitsImage) = get(FitsHeader, obj)
Base.convert(::Type{FITSHeader}, obj::FitsImage) = get(FITSHeader, obj)
Base.convert(::Type{Array}, obj::FitsImage) = get(Array, obj)
Base.convert(T::Type{<:AbstractArray}, obj::FitsImage) =
    convert(T, get(Array, obj))

# Make FitsImage instances behave like arrays (indexing is considered later).
#Base.eltype(::FitsImage{T,N}) where {T,N} = T # FIXME: not needed
#Base.ndims(::FitsImage{T,N}) where {T,N} = N # FIXME: not needed
Base.length(A::FitsImage) = length(get(Array, A))
Base.size(A::FitsImage) = size(get(Array, A))
Base.size(A::FitsImage, d) = size(get(Array, A), d)
Base.axes(A::FitsImage) = axes(get(Array, A))
Base.axes(A::FitsImage, d) = axes(get(Array, A), d)
@inline Base.axes1(A::FitsImage) = axes1(get(Array, A))
Base.IndexStyle(::Type{<:FitsImage}) = IndexLinear()
Base.similar(::Type{FitsImage{T}}, dims::NTuple{N,Int}) where {T,N} =
    FitsImage{T,N}(undef, dims)
Base.elsize(::Type{FitsImage{T,N}}) where {T,N} = elsize(Array{T,N})
Base.sizeof(A::FitsImage) = sizeof(get(Array, A))

# Make FitsImage's efficient iterators.
@inline Base.iterate(A::FitsImage, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

# FIXME: copyto!, convert, unsafe_convert, pointer, etc.


"""

```julia
get([T,] obj, key[, def])
```

yields the value of the FITS keyword `key` in object `obj`.  Argument `T` may
be used to specify the type of the result (i.e. `Int`, `Bool`, `Float64` or
`String`) or `FitsComment` to retrieve the comment associated with the FITS
keyword `key`.  Argument `def` may be used to specify the value to return when
the keyword is not found (`def` is not converted to type `T` if this argument
is specified).  A `KeyError` exception is thrown if the keyword is not found
and no default value `def` is specifed.  An error is thrown if the keyword is
found but its value cannot be converted to `T` if it is specifed.

Multiple keywords can be specified if `obj` is a FITS HDU:

```julia
get(T, obj, keys, def)
```

yields the value of the first FITS keyword out of `keys` found in object `obj`
with type `T` or `def` if none of the keywords is found.  Argument `keys` can
be a single keyword or a tuple of keywords.  An `ErrorException` is thrown if
the value of the first matching keyword cannot be converted to type `T`.

The syntax `get(T,obj)` is also extended to retrieve various things from object
`obj` of type `FitsHeader`, `FitsHDU` or `FitsImage`.  Examples:

```julia
get(Array, obj)            -> arr # get array associated to a `FitsImage`
get(FitsHeader, obj)       -> hdr # get FITS header of object `obj`
get(FitsComment, obj, key) -> str # yields comment for FITS keyword `key` in object `obj`
```

"""
get(obj::Annotated, key::AbstractString, def) = begin
    hdr = get(FITSHeader, obj)
    return (haskey(hdr, key) ? getindex(hdr, key) : def)
end

get(obj::Annotated, key::AbstractString) = begin
    hdr = get(FITSHeader, obj)
    haskey(hdr, key) || missing_keyword(key)
    return getindex(hdr, key)
end

get(::Type{FitsComment}, obj::Annotated, key::AbstractString, def) = begin
    hdr = get(FITSHeader, obj)
    return (haskey(hdr, key) ? FITSIO.get_comment(hdr, key) : def)
end

get(::Type{FitsComment}, obj::Annotated, key::AbstractString) = begin
    hdr = get(FITSHeader, obj)
    haskey(hdr, key) || missing_keyword(key)
    return get_comment(hdr, key)
end

for S in (AbstractFloat, Integer, AbstractString, Bool)
    @eval begin
        function get(::Type{T}, obj::Annotated,
                          key::AbstractString, def) where {T<:$S}
            hdr = get(FITSHeader, obj)
            return (haskey(hdr, key) ?
                    checkvalue(T, $S, getindex(hdr, key), key) : def)
        end

        function get(::Type{T}, obj::Annotated,
                          key::AbstractString) where {T<:$S}
            hdr = get(FITSHeader, obj)
            haskey(hdr, key) || missing_keyword(key)
            return checkvalue(T, $S, getindex(hdr, key), key)
        end
    end
end

function get(::Type{T}, obj::Union{FitsIO,FitsHDU},
             key::AbstractString, def) where {T<:KeywordValues}
    return get(T, obj, (key,), def)
end

function get(::Type{T}, obj::Union{FitsIO,FitsHDU},
             key::AbstractString) :: T where {T<:KeywordValues}
    value = get(T, obj, key, Missing())
    isa(value, T) || missing_keyword(key)
    return value
end

function get(::Type{T}, obj::Union{FitsIO,FitsHDU},
             keys::Tuple{Vararg{AbstractString}},
             def) where {T<:KeywordValues}
    # FIXME: This is a modified version of `FITSIO.fits_try_read_keys`.
    file = getfile(obj)
    status = Ref{Cint}()
    value = Vector{UInt8}(undef, 71)
    for key in keys
        status[] = 0
        ccall((:ffgkey, libcfitsio), Cint,
              (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint}),
              file.ptr, key, value, C_NULL, status)

        # If the key is found, return it. If there was some other error
        # besides key not found, throw an error.
        if status[] == 0
            value = FITSIO.try_parse_hdrval(T, unsafe_string(pointer(value)))
            if isa(value, T)
                return value
            end
            error(string("value of FITS keyword \"", key,
                         "\" cannot be converted to type `", T, "`"))
        elseif status[] != 202
            error(FITSIO.fits_get_errstatus(status[1]))
        end
    end
    return def
end

#
# Override `getindex` and `setindex!` for efficient indexation by array
# indices.
#

@inline @propagate_inbounds getindex(A::FitsImage, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds getindex(get(Array, A), i)
end

@inline @propagate_inbounds setindex!(A::FitsImage, x, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds setindex!(get(Array, A), x, i)
end

@inline Base.checkbounds(::Type{Bool}, A::FitsImage, i::Int) =
    1 ≤ i ≤ length(A)

#
# Override `getindex` and `setindex!` to implement `obj[key]`, `obj[key] = val`
# and `obj[key] = (val, com)` syntaxes.
#

@inline getindex(obj::Annotated, key::AbstractString) =
   getindex(get(FITSHeader, obj), key)

setindex!(obj::Annotated, val, key::AbstractString) =
    setindex!(get(FITSHeader, obj), val, key)

setindex!(obj::Annotated, val::Tuple{Any,AbstractString}, key::AbstractString) = begin
    hdr = get(FITSHeader, obj)
    result = setindex!(hdr, val[1], key)
    set_comment!(hdr, key, val[2])
    return result
end

setindex!(obj::Annotated, val::Tuple{Any,Nothing}, key::AbstractString) =
    setindex!(obj, (val[1], ""), key)

setindex!(obj::Annotated, val::Tuple{Any}, key::AbstractString) =
    setindex!(obj, val[1], key)

# Capture other calls as an error.
setindex!(obj::Annotated, val::Tuple, key::AbstractString) =
    throw(ArgumentError(string("unexpected value type `", typeof(val),
                               "` for FITS keyword")))

#
# Override `getproperty` and `setproperty!` to implement `obj.key` syntax.
#

@inline Base.getproperty(obj::T, sym::Symbol) where {T<:Annotated} =
    getindex(obj, propertyname(T, sym))

@inline Base.setproperty!(obj::T, sym::Symbol, val) where {T<:Annotated} =
    setindex!(obj, val, propertyname(T, sym))

"""

```julia
EasyFITS.propertyname(T, sym)
```

converts symbol `sym` to a suitable key for object of type `T` throwing an
error if this conversion is not supported.

"""
@inline propertyname(::Type{<:Annotated}, sym::Symbol) = String(sym)
@noinline propertyname(::Type{T}, sym::Symbol) where {T} =
    error(string("Converting symbolic key for objects of type `", T, "` ",
                 "is not implemented. As a result, syntax `obj.key` is not ",
                 "supported for such objects."))

keys(obj::Annotated) = keys(get(FITSHeader, obj))
nkeys(obj::Annotated) = nkeys(get(FITSHeader, obj))
nkeys(obj::Union{FITSHeader,AbstractDict}) = length(obj)
haskey(obj::Annotated, key) = haskey(get(FITSHeader, obj), key)
getkey(obj::Annotated, key, def) = getkey(get(FITSHeader, obj), key, def)

function checkvalue(::Type{T}, ::Type{S}, val,
                    key::AbstractString)::T where {S,T}
    isa(val, S) || bad_type(key)
    return T(val)
end

@noinline missing_keyword(key::AbstractString) = throw(KeyError(key))

@noinline bad_type(key::AbstractString) =
    error("bad type for FITS keyword \"$key\"")

"""

FIXME:

```julia
read(T::Union{FitsImage,FitsHeader}, src, ext=1) -> A
```

read an object of type `T` from FITS extension `ext` in source `src`.  Argument
`src` can be the name of a FITS file or an object of type `FitsIO` or
`FitsHDU`.

If `T` is `FitsImage`, the returned value is a pseudo-array `A` with the
contents of the FITS extension `ext` in `src`.

The optional HDU number, the first one by
default, must correspond to a FITS *Image* extension.

The result is indexable.  Using string index yields the value of the
corresponding FITS keyword in the header part of the HDU.  Any other indices
are used to access the contents of data part of the HDU (as a regular Julia
array).

Argument `T` can be used to specify the type of the result. `T` can be `Array`
or `Array{E}` or `Array{E,N}` to retrieve only the data part of the FITS HDU as
a regular array with, if these parameters are specified, element type `E` and
`N` dimensions.  `T` can also be `FitsImage` (the default if `T` is not
specified) or `FitsImage{E}` or `FitsImage{E,N}` to retrieve a pseudo-array,
possibly converted to suitable element type.  Finally, `T` can be `FitsHeader`
to only retrieve the header part.

Examples:

```julia
using EasyFITS
A = read(FitsImage, "image.fits")  # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
A.BITPIX                           # idem
getfitscomment(A, "BITPIX")        # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
A["STUFF"] = (1, "Blah")           # idem with a value-comment pair
A.STUFF = 2                        # set value of FITS keyword STUFF
A.STUFF = (3, "Blah")              # idem with a comment
arr = get(Array, A)                # get the data part (a regular Julia array)
hdr = get(FitsHeader, A)           # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
length(hdr)                        # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
```

See also: [`FitsIO`](@ref).

"""
read(::Type{T}, path::AbstractString, args...; kwds...) where {T<:Readable} =
    FitsIO(path, "r") do io
        return read(T, io, args...; kwds...)
    end

function read(::Type{T}, io::FitsIO,
              ext::Extension = 1) where {T<:Union{Readable,Array}}
    read(T, io[ext])
end

# Read header from FITS ImageHDU.
read(::Type{FitsHeader}, hdu::FitsHDU) =
    FitsHeader(read_header(get(HDU,hdu)))

# Read header and data from FITS Image HDU.
read(::Type{FitsImage{T,N}}, hdu::FitsHDU) where {T,N} =
    FitsImage(read(Array{T,N}, hdu), read(FitsHeader, hdu))
read(::Type{FitsImage}, hdu::FitsHDU) =
    FitsImage(read(Array, hdu), read(FitsHeader, hdu))
read(::Type{FitsImage{T}}, hdu::FitsHDU) where {T} =
    FitsImage(read(Array{T}, hdu), read(FitsHeader, hdu))

# Read array from FITS Image HDU.
read(hdu::FitsImageHDU, args...) = read(get(HDU, hdu), args...)
read(::Type{Array}, hdu::FitsImageHDU, args...) = read(hdu, args...)
read(::Type{Array{T}}, hdu::FitsImageHDU) where {T} =
    convert(Array{T}, read(hdu))
read(::Type{Array{T,N}}, hdu::FitsImageHDU) where {T,N} =
    convert(Array{T,N}, read(hdu))

read(::Type{FitsArray}, hdu::FitsHDU, args...) =
    read(Array, hdu, args...)
read(::Type{FitsArray{T}}, hdu::FitsHDU, args...) where {T} =
    read(Array{T}, hdu, args...)
read(::Type{FitsArray{T,N}}, hdu::FitsHDU, args...) where {T,N} =
    read(Array{T,N}, hdu, args...)

"""

```julia
write(FitsFile, path, args...; overwrite=false, kwds...)
```

creates a new FITS file `path` with contents built from the provided arguments
`args...` and keywords `kwds...`.  Unless `overwrite` is `true`, the file must
not already exist.  Instead of setting the `overwrite` keyword, simply call:

```julia
write!(FitsFile, path, args...; kwds...)
```

to silently overwrites file `path` if it already exits.

Typical examples are:

```julia
write(FitsFile, path, arr; KEY1 = val1, KEY2 = (val2, com2), ...)
write(FitsFile, path, hdr, arr)
write(FitsFile, path, (hdr1, arr1), arr2, (hdr3, arr3), ...)
```

with `arr*` arrays, `KEY*` keywords, `val*` keyword values, `com*` keyword
comments and `hdr*` objects of type `FitsHeader`.  In the last example, a FITS
file with 3 (or more) HDU's is created.

"""
function write!(::Type{FitsFile}, path::AbstractString, args...; kwds...)
    FitsIO(path, "w!") do io
        write(io, args...; kwds...)
    end
    nothing
end

function write(::Type{FitsFile}, path::AbstractString, args...;
               overwrite::Bool=false, kwds...)
    (overwrite == false && exists(path)) &&
        throw_file_already_exists(path, "try with `overwrite=true`")
    FitsIO(path, (overwrite ? "w!" : "w")) do io
        write(io, args...; kwds...)
    end
    nothing
end

write(io::FitsIO, arr::AbstractArray{<:Real}; kwds...) =
    write(io, FitsHeader(; kwds...), arr)

write(io::FitsIO, arg::Tuple{FitsHeader,AbstractArray{<:Real}}) = write(io, arg...)
write(io::FitsIO, arg::Tuple{AbstractArray{<:Real},FitsHeader}) = write(io, arg...)
write(io::FitsIO, arr::AbstractArray{<:Real}, hdr::FitsHeader) = write(io, hdr, arr)
write(io::FitsIO, hdr::FitsHeader, arr::AbstractArray{<:Real}) =
    write(get(FITS, io), arr, header = get(FITSHeader, hdr))
write(io::FitsIO, obj::FitsImage) =
    write(io, get(FitsHeader, obj), get(Array, obj))
write(io::FitsIO, args...) =
    for arg in args
        write(io, arg)
    end

"""

```julia
EasyFITS.getfile(obj) -> file
```

yeilds the FITS file associated with object `obj` if it is open; otherwise
throws an error.  Object `obj` can be an instance of `FitsIO`, `FitsHDU`,
`FITSIO.FITS`, `FITSIO.HDU` or `FITSIO.FITSFile`.

```julia
EasyFITS.getfile(obj, ext) -> file
```

throws an error if the FITS file associated with `obj` is not open; otherwise
moves to HDU number `ext` and returns the FITS file.  Object `obj` can be an
instance of `FitsIO`, `FITSIO.FITS`, or `FITSIO.FITSFile`.

This method is considered as *low-level* and is not exported.

"""
function getfile(file::FITSFile)
    fits_assert_open(file)
    return file
end

function getfile(file::FITSFile, ext::Integer)
    fits_assert_open(file)
    fits_movabs_hdu(file, ext)
    return file
end

getfile(fh::FITS) = getfile(fh.fitsfile)
getfile(fh::FITS, ext::Integer) = getfile(fh.fitsfile, ext)

getfile(hdu::HDU) = getfile(hdu.fitsfile, hdu.ext)

getfile(io::FitsIO) = getfile(get(FITS, io))
getfile(io::FitsIO, ext::Integer) = getfile(get(FITS, io), ext)

getfile(hdu::FitsHDU) = getfile(get(HDU, hdu))


@noinline throw_file_already_exists(args::AbstractString...) =
    error(file_already_exists(args...))

file_already_exists(path::AbstractString) =
    string("file \"", path, "\" already exists")

file_already_exists(path::AbstractString, usage::AbstractString) =
    string(file_already_exists(path), ", ", usage)

include("deprecate.jl")

end # module
