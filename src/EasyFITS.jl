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
# Copyright (C) 2018-2020, Éric Thiébaut.
#

module EasyFITS

export
    FitsArray,
    FitsBitpix,
    FitsComment,
    FitsFile,
    FitsHDU,
    FitsHeader,
    FitsIO,
    FitsImage,
    FitsImageHDU,
    FitsTableHDU,
    exists,
    readfits,
    writefits,
    writefits!,
    write!

using FITSIO
using CFITSIO
using CFITSIO: libcfitsio

using Base: getfield, elsize, tail, OneTo,
    throw_boundserror, @propagate_inbounds
import Base: get, getindex, setindex!, keys, haskey, getkey,
    read, write, open, close

# Apply patches.
include("patches.jl")

#
# Predefine all types.
#

struct FitsIO
    fh::FITS
end

struct FitsBitpix{N} end

struct FitsHDU{T<:HDU}
    hdu::T
end

const FitsImageHDU{T<:ImageHDU} = FitsHDU{T}
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

# This decoration is used to specialize `tryparse` for FITS keyword values.
struct FitsUnparsedValue
    str::String
end

# Annotated objects have a FITS header and implements indexation by keywords,
# they can also be directly read from a FITS file.
const Annotated = Union{FitsImage,FitsHeader}

# Readable objects are those which can be directly read from a FITS file.
const Readable = Union{Annotated,FitsArray}

# Allowed types for keywords values.
const KeywordValues = Union{Bool,Int,Float64,String}

# Valid types to specify a FITS extension.
const Extension = Union{Integer,AbstractString}

# Valid types to specify a FITS header.
const HeaderLike = Union{FitsHeader,NamedTuple,
                         Tuple{Vararg{Pair{<:AbstractString}}}}

readfits(filename::AbstractString; ext::Extension = 1) =
    read(FitsImage, filename; ext=ext)

readfits(::Type{Array}, filename::AbstractString; ext::Extension = 1) =
    read(FitsArray, filename; ext=ext)

function readfits(::Type{Array{T}}, filename::AbstractString;
                  ext::Extension = 1) where {T}
    read(FitsArray{T}, filename; ext=ext)
end

function readfits(::Type{Array{T,N}}, filename::AbstractString;
                  ext::Extension = 1) where {T,N}
    read(FitsArray{T,N}, filename; ext=ext)
end

writefits(filename::AbstractString, args...; kwds...) =
    write(FitsFile, filename, args...; kwds...)

writefits!(filename::AbstractString, args...; kwds...) =
    write(FitsFile, filename, args...; kwds...)

"""
    exists(path) -> boolean

yields whether file `path` exists.  Argument can be a file name or an instance
of `Base.Filesystem.StatStruct` as returned by the `stat` method.

"""
exists(st::Base.Filesystem.StatStruct) = (st.nlink != 0)
exists(path) = exists(stat(path))

"""
    EasyFITS.hduname(T) -> name, vers=1

yields the name and revision number of FITS Header Data Unit for storing data
of type `T`.  Other packages are encouraged to extend this method with their
types.

"""
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
    FitsIO(path, mode="r") -> io

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

Call `close(io)` to close the FITS file associated with the `FitsIO` instance
`io`.  Call `isopen(io)` to check whether the FITS file associated with the
`FitsIO` instance `io` is open.  Closing the FITS file is automatically done,
if needed, when the instance is garbage collected.

The do-block syntax is supported to automatically close the FITS file:

    FitsIO(filename, mode="r") do io
        # use FITS handle io
        ...
    end

An instance of `FitsIO` is a collection of *Header Data Units* (HDU) and
implements indexation and iteration.  Assuming `io` is a `FitsIO` object, then:

- `io[i]` yields the `i`-th HDU.

- `length(io)` yields the number of HDUs.

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

FitsIO(func::Function, path::AbstractString, mode::AbstractString="r") = begin
    io = FitsIO(path, mode)
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
    EasyFITS.find(pred, io, I=1:length(io))

yields the first index HDU `i ∈ I` of FITS instance `io` for which the
predicate `pred(io[i])` returns `true`, or [`nothing`](@ref) if `pred(io[i])`
returns `false` for all `i ∈ I`.

For instance:

    i = EasyFITS.find(hdu -> get(String,hdu,"EXTNAME","") == "CALIBRATION", io)
    if i === nothing
        # not found
        ...
    else
        # found at index i
        hdu = io[i]
    end

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
    isprimary(hdu)

yields whether `hdu` is a primary FIT header data units (HDU).

"""
isprimary(obj::FitsHDU) = isprimary(get(HDU, obj))
isprimary(obj::HDU) = (getfield(obj, :ext) == 1)

"""

An object of type `FitsHeader` stores the header part of a FITS *Header Data
Unit* (HDU).  There are several ways to build an instance of `FitsHeader`.

    FitsHeader()

yields an empty instance of `FitsHeader`.

    FitsHeader(io, ext=1)

yields an instance of `FitsHeader` built over FITS header read from extension
`ext` in FITS file `io`.

    FitsHeader(; key1=val1, key2=val2, ...)

yields an instance of `FitsHeader` whose contents is set by keywords (values
can be a tuple with a value and a comment).  This is can also be done by
specifying key-value pairs:

    FitsHeader("key1" => val1, "key2" => val2, ...)

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

    FitsHeader(hdr)

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
    FitsBitpix(arg)

yields a singleton type which encapsulates the FITS *bitpix* (for
*bits-per-pixel*) code identifying array element type according to `arg`.
Argument `arg` can be a Julia type, an integer (interpreted as a FITS bitpix
code), an array, a FITS HDU/image/header.

Call `eltype(bpx)` to convert FITS bitpix `bpx` into a Julia type.

"""
FitsBitpix(n::Integer) = FitsBitpix{Int(n)}()
FitsBitpix(::Type{T}) where {T} = FitsBitpix(bitpix_from_type(T))
FitsBitpix(A::AbstractArray) = FitsBitpix(typeof(A))
FitsBitpix(::Type{<:AbstractArray{T}}) where {T} = FitsBitpix(T)
FitsBitpix(hdr::Union{<:FitsHeader,<:FITSHeader}) = FitsBitpix(hdr["BITPIX"])
FitsBitpix(hdu::Union{<:FitsImageHDU,<:ImageHDU}) =
    FitsBitpix(fits_get_img_equivtype(getfile(hdu)))

Base.eltype(::FitsBitpix{N}) where {N} = type_from_bitpix(Val(Cint(N)))

# FIXME: Base.eltype could also be extended to FITSHeader and ImageHDU.
Base.eltype(hdr::FitsHeader) =
    type_from_bitpix(Val(Cint(hdr["BITPIX"])))
Base.eltype(hdu::FitsImageHDU) =
    type_from_bitpix(Val(Cint(fits_get_img_equivtype(getfile(hdu)))))

"""
    FitsImage(arr, hdr)

yields a FITS Image that behaves like a Julia array when indexed by integers of
Cartesian indices and like a dictionary of FITS keywords when indexed by
strings.  Argument `arr` specifies the array contents of the object and `hdr`
specifies the keywords both are shared by the returned object.

The initial contents of the header may be specified by keywords.  For instance:

    FitsImage(arr; KEY1 = val1, KEY2 = (val2, com2), ...)

The embedded array may be creted by the constructor:

    FitsImage{T}(init, dims...; kwds...)

yields a FITS Image with an empty header and array-like data of dimensions
`dims` and with elements of type `T` initially set by `init`.  Specifying the
number of dimensions is also supported:

    FitsImage{T,N}(init, dims...; kwds...)

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
Base.similar(A::FitsImage, T::Type=eltype(A), dims::NTuple{N,Int}=size(A)) where {N} =
    FitsImage{T,N}(undef, dims)
Base.elsize(::Type{FitsImage{T,N}}) where {T,N} = elsize(Array{T,N})
Base.sizeof(A::FitsImage) = sizeof(get(Array, A))

# For AbstractArray types, a more specific version of `Base.show` must be extended.
Base.show(io::IO, ::MIME"text/plain", A::FitsImage) = show(io, A)
function Base.show(io::IO, A::FitsImage{T,N}) where {T,N}
    println(io, "FitsImage{$T,$N} of size $(size(A)) and keywords:")
    show(io, FitsHeader(A))
end

# Make FitsImage's efficient iterators.
@inline Base.iterate(A::FitsImage, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

# FIXME: copyto!, convert, unsafe_convert, pointer, etc.


# Maximum size of the value field including the final '\0'.
const VALUE_SIZE = 72

"""
    findkeys(obj, keys [, buf])

finds a FITS card in `obj` with one of the keywords in `keys` and returns a
string with the unparsed value of the keyword if found and `nothing` otherwise.

Optional argument `buf` is a small scratch buffer of bytes (`UInt8`) used to
temporarily store the unparsed value before it is converted into a Julia
string, its size is augmented if it is smaller than `VALUE_SIZE`.

"""
function findkeys(obj::Union{FitsIO,FitsHDU},
                  keys::NTuple{N,AbstractString},
                  buf::Vector{UInt8}=Vector{UInt8}(undef, VALUE_SIZE)) where {N}
    file = getfile(obj)
    status = Ref{Cint}()
    length(buf) ≥ VALUE_SIZE || resize!(buf, VALUE_SIZE)
    for key in keys
        status[] = 0
        ccall((:ffgkey, libcfitsio), Cint,
              (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint}),
              file.ptr, key, buf, C_NULL, status)

        # If the key is found, return it. If there was some other error
        # besides key not found, throw an error.
        if status[] == 0
            return unsafe_string(pointer(buf))
        elseif status[] != 202
            error(FITSIO.fits_get_errstatus(status[1]))
        end
    end
    return nothing
end

@doc @doc(trygetkeys) VALUE_SIZE

"""
    get([T,] obj, key[, def])

yields the value of the FITS keyword `key` in object `obj`.  Argument `T` may
be used to specify the type of the result (i.e. `Int`, `Bool`, `Float64` or
`String`) or `FitsComment` to retrieve the comment associated with the FITS
keyword `key`.  Argument `def` may be used to specify the value to return when
the keyword is not found (`def` is not converted to type `T` if this argument
is specified).  A `KeyError` exception is thrown if the keyword is not found
and no default value `def` is specifed.  An error is thrown if the keyword is
found but its value cannot be converted to `T` if it is specifed.

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

function get(::Type{T}, obj::Union{FitsIO,FitsHDU}, key::AbstractString,
             def = Missing()) where {T<:KeywordValues}
    raw = findkeys(obj, (key,))
    if raw === nothing
        def === Missing() && missing_keyword(key)
        return def
    end
    if iscomment(key)
        if T == String
            return raw
        end
    else
        val = tryparse(T, FitsUnparsedValue(raw))
        if isa(val, T)
            return val
        end
    end
    error("value of FITS keyword \"", key,
          "\" cannot be converted to type `", T, "`")
end

function get(obj::Union{FitsIO,FitsHDU}, key::AbstractString, def = Missing())
    raw = findkeys(obj, (key,))
    if raw === nothing
        def === Missing() && missing_keyword(key)
        return def
    end
    if iscomment(key)
        return raw
    end
    c = first(raw)
    if c == 'T' || c == 'F'
        if length(raw) == 1
            return (c == 'T')
        end
    elseif c == Char(39) # single quote, using '\'' break Emacs fontify
        sval = tryparse(String, FitsUnparsedValue(raw))
        if isa(sval, String)
            return sval
        end
    else
        ival = tryparse(Int, FitsUnparsedValue(raw))
        if isa(ival, Int)
            return ival
        end
        rval = tryparse(Float64, FitsUnparsedValue(raw))
        if isa(rval, Float64)
            return rval
        end
    end
    error("unknown type of value for FITS keyword \"", key, "\"")
end

iscomment(key::AbstractString) =
    (key == "HISTORY" || key == "COMMENT")

Base.tryparse(::Type{T}, val::FitsUnparsedValue) where {T<:Real} =
    tryparse(T, val.str)

function Base.tryparse(::Type{Bool}, val::FitsUnparsedValue)
    if length(val.str) == 1
        if val.str[1] == 'T'
            return true
        elseif val.str[1] == 'F'
            return false
        end
    end
    return nothing
end

function Base.tryparse(::Type{String}, val::FitsUnparsedValue)
    # Check that the string in enclosed by single quotes.
    QUOTE = Char(39) # a single quote, avoing breaking Emacs' indentation
    str = val.str
    i = firstindex(str)
    j = lastindex(str)
    if length(str) < 2 || str[i] != QUOTE || str[j] != QUOTE
        return nothing
    end

    # Replace pairs of quotes by a single quote and discard trailing spaces.
    buf = Vector{UInt8}(undef, length(str) + 1)
    spaces = 0  # number of unwritten trailing spaces
    esc = false # previous character was a quote?
    n = 0
    for k in (i+1):(j-1)
        c = str[k]
        if esc
            if c != QUOTE
                return nothing
            end
            esc = false
        elseif c == QUOTE
            esc = true
            continue
        end
        if c == ' ' # FIXME: use `isspace`?
            spaces += 1
        else
            # Write any unwritten spaces and then the new non-space character.
            while spaces > 0
                spaces -= 1
                n += 1
                buf[n] = ' '
            end
            n += 1
            buf[n] = c
        end
    end
    if esc
        return nothing
    end
    buf[n+1] = 0
    return unsafe_string(pointer(buf))
end

#
# Override `getindex` and `setindex!` for efficient indexation by array
# indices.
#

@inline getindex(A::FitsImage, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds val = get(Array, A)[i]
    return val
end

@inline setindex!(A::FitsImage, val, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds get(Array, A)[i] = val
    return A
end

@inline Base.checkbounds(::Type{Bool}, A::FitsImage, i::Int) =
    1 ≤ i ≤ length(A)

#
# Override `getindex` and `setindex!` to implement `obj[key]`, `obj[key] = val`
# and `obj[key] = (val, com)` syntaxes.
#

@inline getindex(obj::Annotated, key::AbstractString) =
   getindex(get(FITSHeader, obj), key)

setindex!(obj::Annotated, val, key::AbstractString) = begin
    setindex!(get(FITSHeader, obj), val, key)
    return obj
end

setindex!(obj::Annotated, val::Tuple{Any,AbstractString}, key::AbstractString) = begin
    hdr = get(FITSHeader, obj)
    setindex!(hdr, val[1], key)
    set_comment!(hdr, key, val[2])
    return obj
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
    EasyFITS.propertyname(T, sym)

converts symbol `sym` to a suitable key for object of type `T` throwing an
error if this conversion is not supported.

"""
@inline propertyname(::Type{<:Annotated}, sym::Symbol) = String(sym)
@noinline propertyname(::Type{T}, sym::Symbol) where {T} =
    error(string("Converting symbolic key for objects of type `", T, "` ",
                 "is not implemented. As a result, syntax `obj.key` is not ",
                 "supported for such objects."))

keys(obj::Annotated) = keys(get(FITSHeader, obj))
nkeys(obj::Annotated) = length(get(FITSHeader, obj))
nkeys(obj::Union{FITSHeader,AbstractDict}) = length(obj)
haskey(obj::Annotated, key) = haskey(get(FITSHeader, obj), key)
getkey(obj::Annotated, key, def) = getkey(get(FITSHeader, obj), key, def)

Base.delete!(obj::Annotated, i::Integer) = delete!(obj, Int(i))
Base.delete!(obj::Annotated, key::String) = delete!(obj, first_card(obj, key))
Base.delete!(obj::Annotated, idx::Int) = begin
    delete_card!(get(FITSHeader, obj), idx)
    return obj
end

Base.pop!(obj::Annotated, idx::Integer) = pop!(obj, Int(idx))
Base.pop!(obj::Annotated, idx::Integer, def) = pop!(obj, Int(idx), def)
Base.pop!(obj::Annotated, key::String) = pop!(obj, first_card(obj, key))
Base.pop!(obj::Annotated, key::String, def) = begin
    idx = first_card(obj, key, -1)
    return (idx == 0 ? def : pop!(obj, idx))
end
Base.pop!(obj::Annotated, idx::Int) = begin
    1 ≤ idx ≤ nkeys(obj) || error("out of bound FITS card index $idx")
    return unsafe_pop!(obj, idx)
end
Base.pop!(obj::Annotated, idx::Int, def) = begin
    1 ≤ idx ≤ nkeys(obj) || return def
    return unsafe_pop!(obj, idx)
end

unsafe_pop!(obj::Annotated, idx::Int) = begin
    hdr = get(FITSHeader, obj)
    val = hdr.values[idx]
    delete_card!(hdr, idx)
    return val
end

# Delete FITS card in header.  NOTE: Do not extend `delete!` to avoid
# type-piracy.  FIXME: Deleting a FITS card is expensive.
delete_card!(hdr::FITSHeader, i::Int) = begin
    if 1 ≤ i ≤ length(hdr)
        # Remove the entry from hdr.keys, hdr.values and hdr.comments and
        # rebuild the keyword map.
        delete_entry!(hdr.keys,     i)
        delete_entry!(hdr.values,   i)
        delete_entry!(hdr.comments, i)
        empty!(hdr.map)
        for j in 1:length(hdr.keys)
            key = hdr.keys[j]
            if !haskey(hdr.map, key)
                hdr.map[key] = j
            end
        end
    end
    nothing
end

# Delete an entry in a vector.
delete_entry!(A::Vector, i::Int) = begin
    n = length(A)
    if 1 ≤ i ≤ n
        for j in i:n-1
            A[j] = A[j+1]
        end
        resize!(A, n - 1)
    end
    nothing
end

# Yields index of first FITS card with given keyword.
first_card(obj::Annotated, key::String, def) =
    get(get(FITSHeader, obj).map, key, def)
first_card(obj::Annotated, key::String) = begin
    idx = first_card(obj, key, 0)
    idx == 0 && throw(KeyError(key))
    return idx
end

function checkvalue(::Type{T}, ::Type{S}, val,
                    key::AbstractString)::T where {S,T}
    isa(val, S) || bad_type(key)
    return T(val)
end

@noinline missing_keyword(key::AbstractString) = throw(KeyError(key))

@noinline bad_type(key::AbstractString) =
    error("bad type for FITS keyword \"$key\"")

"""
    read(::Type{T}, src) -> A

read an object of type `T<:Union{FitsImage,FitsHeader,FitsArray}` from FITS
extension `ext` in source `src` which can be the name of a FITS file or an
instance of `FitsIO` or `FitsHDU`.

If `src` is an instance of `FitsIO` or the name of a FITS file, keyword `ext`
may be used to specify the number or the name of the HDU number to read (the
first one by default).

Argument `T` specifies the type of the result.  Use `T<:Array` or
`T<:FitsArray` to retrieve the data part of the FITS HDU as a regular array.
Use `T<:FitsImage` to get a pseudo-array combining the header and data parts of
the FITS HDU.  Finally, use `T<:FitsHeader` to only retrieve the header part.

For `T<:Array` `T<:FitsArray`, or `T<:FitsImage`, the element type and number
of dimensions of the result can be specified, e.g. with `FitsImage{Float32,2}`

To avoid *type-piracy*, `T<:FitsArray` must be used instead of `T<:Array` when
`src` is a file name.

The result is indexable.  For `T<:FitsImage` and `T<:FitsHeader`, using a
string index yields the value of the corresponding FITS keyword in the header
part of the HDU.  Any other indices are used to access the contents of data
part of the HDU (as a regular Julia array).

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

function read(::Type{T}, io::FitsIO, args...;
              ext::Extension = 1) where {T<:Union{Readable,Array}}
    read(T, io[ext], args...)
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
read(::Type{Array{T}}, hdu::FitsImageHDU, args...) where {T} =
    convert(Array{T}, read(hdu, args...))
read(::Type{Array{T,N}}, hdu::FitsImageHDU, args...) where {T,N} =
    convert(Array{T,N}, read(hdu, args...))

read(::Type{FitsArray}, hdu::FitsHDU, args...) =
    read(Array, hdu, args...)
read(::Type{FitsArray{T}}, hdu::FitsHDU, args...) where {T} =
    read(Array{T}, hdu, args...)
read(::Type{FitsArray{T,N}}, hdu::FitsHDU, args...) where {T,N} =
    read(Array{T,N}, hdu, args...)

function Base.read!(hdu::FitsHDU, dst, args...)
    read!(get(HDU, hdu), dst, args...)
    return dst
end

"""
    write(FitsFile, path, args...; overwrite=false, kwds...)

creates a new FITS file `path` with contents built from the provided arguments
`args...` and keywords `kwds...`.  Unless `overwrite` is `true`, the file must
not already exist.  Instead of setting the `overwrite` keyword, simply call:

    write!(FitsFile, path, args...; kwds...)

to silently overwrites file `path` if it already exits.

Typical examples are:

    write(FitsFile, path, arr; KEY1 = val1, KEY2 = (val2, com2), ...)
    write(FitsFile, path, hdr, arr)
    write(FitsFile, path, (hdr1, arr1), arr2, (hdr3, arr3), ...)

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

write!(path::AbstractString, A::FitsImage) =
    write!(FitsFile, path, A)

write(path::AbstractString, A::FitsImage; overwrite::Bool=false) =
    write(FitsFile, path, A; overwrite=overwrite)

"""
    write(io::FitsIO, args...; kwds...)

writes data `args...` in FITS output `io` with keywords `kwds...`.  This method
can be specialized by foreign packages based on the types of the arguments
`args...`.

"""
write(io::FitsIO, arr::AbstractArray{<:Real}; kwds...) =
    write(io, FitsHeader(; kwds...), arr)
write(io::FitsIO, arr::AbstractArray{<:Real}, args::Pair{<:AbstractString}...) =
    write(io, FitsHeader(args...), arr)
write(io::FitsIO, hdr::FitsHeader, arr::AbstractArray{<:Real}) =
    write(get(FITS, io), arr, header = get(FITSHeader, hdr))
write(io::FitsIO, obj::FitsImage) =
    write(io, get(FitsHeader, obj), get(Array, obj))
write(io::FitsIO, args...) =
    for arg in args
        write(io, arg)
    end

# The following versions can build the header on the fly.
write(io::FitsIO, arg::Tuple{HeaderLike,AbstractArray{<:Real}}) =
        write(io, arg...)
write(io::FitsIO, arg::Tuple{AbstractArray{<:Real},HeaderLike}) =
    write(io, arg...)
write(io::FitsIO, arr::AbstractArray{<:Real}, hdr::HeaderLike) =
    write(io, hdr, arr)
function write(io::FitsIO,
               hdr::Tuple{Vararg{Pair{<:AbstractString}}},
               arr::AbstractArray)
    write(io::FitsIO, FitsHeader(hdr...), arr)
end
function write(io::FitsIO,
               hdr::NamedTuple,
               arr::AbstractArray)
    write(io::FitsIO, FitsHeader(; hdr...), arr)
end

"""
    EasyFITS.getfile(obj) -> file

yields the FITS file associated with object `obj` if it is open; otherwise
throws an error.  Object `obj` can be an instance of `FitsIO`, `FitsHDU`,
`FITSIO.FITS`, `FITSIO.HDU` or `FITSIO.FITSFile`.

    EasyFITS.getfile(obj, ext) -> file

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

get(::Type{FITSFile}, io::FitsIO) = getfield(get(FITS, io), :fitsfile)
get(::Type{FITSFile}, hdu::FitsHDU) = getfield(get(HDU, hdu), :fitsfile)
Base.isopen(obj::Union{FitsIO,<:FitsHDU}) = isopen(get(FITSFile, obj))
Base.close(obj::Union{FitsIO,<:FitsHDU}) = close(get(FITSFile, obj))

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

end # module
