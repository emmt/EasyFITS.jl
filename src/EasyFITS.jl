#
# EasyFITS.jl -
#
# Utility routines to simplify reading/writing FITS files.
#
#-------------------------------------------------------------------------------
#
# This file if part of the TAO software (https://github.com/emmt/TAO) licensed
# under the MIT license.
#
# Copyright (C) 2018-2019, Éric Thiébaut.
#

module EasyFITS

export
    FitsComment,
    FitsHeader,
    FitsImage,
    createfits!,
    createfits,
    exists,
    getfitskey,
    openfits,
    readfits,
    setfitskey!,
    tryreadfitskey,
    tryreadfitskeys,
    writefits!,
    writefits

using FITSIO
using FITSIO.Libcfitsio
using FITSIO: fits_try_read_keys

using Base: elsize, tail, OneTo, throw_boundserror, @propagate_inbounds
import Base: getindex, setindex!, keys, haskey, getkey

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

Type `FitsHeader` is a simple wrapper over `FITSIO.FITSHeader` (beware of the
different spellings) to implement indexation by keywords and `obj.key` syntax.
We have to define a new type for that because implementing the extended
interface over the existing `FITSIO.FITSHeader` would be considered as *type
piracy*.

There are several ways to built an instance of `FitsHeader`:

```julia
FitsHeader(hdr)
```

yields an instance of `FitsHeader` built over FITS header object `hdr`.

```julia
FitsHeader(fh, ext=1)
```

yields an instance of `FitsHeader` built over FITS header read from
extension `ext` in FITS file `fh`.

```julia
FitsHeader()
```

yields an empty instance of `FitsHeader`.

```julia
FitsHeader(; key1=val1, key2=val2, ...)
```

yields an instance of `FitsHeader` whose contents is set by keywords (values
can be a tuple with a value and a comment).  This is can also be done by
specifying key-value pairs

```julia
FitsHeader("key1" => val1, "key2" => val2, ...)
```

To avoid ambiguities the two styles cannot be mixed.

"""
struct FitsHeader
    hdr::FITSHeader
end

FitsHeader(hdu::HDU) = FitsHeader(read_header(hdu))
FitsHeader(fh::FITS, ext::Union{Integer,AbstractString} = 1) =
    FitsHeader(fh[ext])

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

Base.convert(::Type{FitsHeader}, hdr::FITSHeader) = FitsHeader(hdr)
Base.convert(::Type{FITSHeader}, hdr::FitsHeader) = get(FITSHeader, hdr)

"""

Type `FitsComment` is a singleton used as a marker to indicate that the coment
of a FITS keyword is to be returned by the `get` method.

"""
struct FitsComment end

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

"""
struct FitsImage{T,N} <: DenseArray{T,N}
    arr::Array{T,N}
    hdr::FitsHeader
end

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

# Make FitsImage instances behave like arrays (indexing is considered later).
#Base.eltype(::FitsImage{T,N}) where {T,N} = T # FIXME: not needed
#Base.ndims(::FitsImage{T,N}) where {T,N} = N # FIXME: not needed
Base.length(obj::FitsHeader) = length(get(FITSHeader, obj))
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


# Annotated objects have a FITS header and implements indexation by keywords.
const Annotated = Union{FitsImage,FitsHeader}

"""

The syntax `get(T,obj)` is extended to retrieve various things from FITS object `obj`:

```julia
get(Array, obj)            -> arr # get array associated a `FitsImage`
get(FitsHeader, obj)       -> hdr # get FITS header of object `obj`
get(FitsComment, obj, key) -> str # yields comment for FITS keyword `key` in object `obj`
```

"""
@inline Base.get(::Type{Array}, obj::FitsImage) = Base.getfield(obj, :arr)
@inline Base.get(::Type{Array}, obj::FitsHeader) = nothing
@inline Base.get(::Type{FitsHeader}, obj::FitsImage) = Base.getfield(obj, :hdr)
@inline Base.get(::Type{FITSHeader}, obj::FitsImage) = get(FITSHeader, get(FitsHeader, obj))
@inline Base.get(::Type{FitsHeader}, obj::FitsHeader) = obj
@inline Base.get(::Type{FITSHeader}, obj::FitsHeader) = Base.getfield(obj, :hdr)
@inline Base.get(::Type{FitsComment}, obj::Annotated, key) = FITSIO.get_comment(get(FITSHeader, obj), key)
FITSIO.get_comment(obj::Annotated, key) = getfitscomment(obj, key)

#
# Override `getindex` and `setindex!` for indexation by array indices or by
# keywords and override `getproperty` and `setproperty!` to implement `obj.key`
# syntax.
#
@inline @propagate_inbounds getindex(A::FitsImage, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds r = getindex(get(Array, A), i)
    return r
end

@inline @propagate_inbounds setindex!(A::FitsImage, x, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds r = setindex!(get(Array, A), x, i)
    return r
end

@inline getindex(obj::Annotated, key::AbstractString) =
    getfitskey(obj, key)

@inline setindex!(obj::Annotated, val, key::AbstractString) =
    setfitskey!(obj, key, val)

@inline Base.getproperty(obj::T, sym::Symbol) where {T<:Annotated} =
    getindex(obj, propertyname(T, sym))

@inline Base.setproperty!(obj::T, sym::Symbol, val) where {T<:Annotated} =
    setindex!(obj, val, propertyname(T, sym))

@inline Base.checkbounds(A::FitsImage, i::Int) =
    1 ≤ i ≤ length(A) || throw_boundserror(A, i)

keys(obj::Annotated) = keys(get(FITSHeader, obj))
nkeys(obj::Annotated) = nkeys(get(FITSHeader, obj))
nkeys(hdr::FITSHeader) = length(hdr)
haskey(obj::Annotated, key) = haskey(get(FITSHeader, obj), key)
getkey(obj::Annotated, key, def) = getkey(get(FITSHeader, obj), key, def)

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

"""

```julia
getfitskey([T,] obj, key[, def]) -> val :: T
```

yields the value of FITS keyword `key` in `obj`.  Argument `obj` can be an
instance of `FitsImage`, `FitsHeader` or `FITSIO.FITSHeader`.  The returned
value has type `T`.  Optional argument `def` is to specify a default value if
keyword `key` is not found in `hdr`.  An error is thrown if keyword `key` is
not found in `hdr` and no default value is provided or if the value to be
returned has not a type that can be converted to `T`.

Call `getfistcomment` to retrieve the comment assiciated with a FITS keyword.

See also: [`setfitskey!`](@ref), [`getfitscomment`](@ref), [`readfits`](@ref).

"""
getfitskey(obj::Annotated, key::AbstractString) =
    getindex(get(FITSHeader, obj), key)

getfitskey(::Type{T}, obj::Annotated, args...) where {T} =
    getfitskey(T, get(FITSHeader, obj), args...)

for S in (AbstractFloat, Integer, AbstractString, Bool)
    @eval begin
        function getfitskey(::Type{T}, hdr::FITSHeader,
                            key::AbstractString, def) where {T<:$S}
            val = haskey(hdr, key) ? getindex(hdr, key) : def
            return checkvalue(T, $S, val, key)
        end

        function getfitskey(::Type{T}, hdr::FITSHeader,
                            key::AbstractString) where {T<:$S}
            haskey(hdr, key) || missing_keyword(key)
            return checkvalue(T, $S, getindex(hdr, key), key)
        end
    end
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

```julia
setfitskey!(obj, key, val[, com])
```

set FITS keyword `key` in `obj` to the value `val`.  Argument `obj` can be an
instance of `FitsImage`, `FitsHeader` or `FITSIO.FITSHeader`.  Argument `com`
is an optional comment; if it is not specified and `key` already exists in
`dst`, its comment is preserved.  Argument `com` can also be `nothing` to erase
the comment.

Arguments `val` and `com` can also be specified as a tuple.

See also: [`getfitskey`](@ref), [`readfits`](@ref).

"""
setfitskey!(obj::Annotated, key::AbstractString, args...) =
    setfitskey!(get(FITSHeader, obj), key, args...)

setfitskey!(hdr::FITSHeader, key::AbstractString, val) = begin
    setindex!(hdr, val, key)
    return nothing
end

setfitskey!(hdr::FITSHeader, key::AbstractString, val, com::AbstractString) = begin
    setindex!(hdr, val, key)
    set_comment!(hdr, key, com)
    return nothing
end

setfitskey!(hdr::FITSHeader, key::AbstractString, val, com::Nothing) =
    setfitskey!(hdr, key, val, "")

setfitskey!(hdr::FITSHeader, key::AbstractString, val::Tuple{<:Any}) =
    setfitskey!(hdr, key, val[1], nothing)

function setfitskey!(hdr::FITSHeader, key::AbstractString,
                     val::Tuple{<:Any,<:Union{AbstractString,Nothing}})
    setfitskey!(hdr, key, val[1], val[2])
end

"""

```julia
openfits(filename) -> fh
```

yields a `FITS` handle to read the contents of the existing FITS file whose
name is `filename`.  If `filename` does not exist but `filename` does not
end with the `".gz"` extension and `"\$filename.gz"` does exist, then the
compressed file `"\$filename.gz"` is open instead.

The do-block syntax is supported to automatically close the FITS file:

```julia
openfits(filename) do fh
   # use FITS handle fh
   ...
end
```

See also: [`readfits`](@ref), [`createfits`](@ref).

"""
openfits(filename::AbstractString) =
    FITS(exists(filename) || endswith(filename, ".gz") ||
         !exists(filename*".gz") ? filename : filename*".gz")

openfits(func::Function, filename::AbstractString) = begin
    io = openfits(filename)
    try
        func(io)
    finally
        #println("closing FITS file...")
        close(io)
    end
end

"""

```julia
createfits!(filename) -> fh
```

creates a new FITS file named `fillename` and return its handle.  If the
file already exists, it is (silently) overwritten.

The `createfits` method can be used instead to avoid that.  The syntax is
identical:

```julia
createfits!(filename) -> fh
```

The do-block syntax is supported to automatically close the FITS file:

```julia
createfits!(filename) do fh
   # use FITS handle fh
   ...
end
```

or, if you want to make sure that `filename` does not exists:

```julia
createfits(filename) do fh
   # use FITS handle fh
   ...
end
```

See also: [`readfits`](@ref), [`createfits`](@ref).

"""
createfits!(filename::AbstractString) = FITS(filename, "w")
createfits!(func::Function, filename::AbstractString) = begin
    io = createfits!(filename)
    try
        func(io)
    finally
        #println("closing FITS file...")
        close(io)
    end
end

createfits(filename::AbstractString; overwrite::Bool=false) = begin
    if !overwrite && exists(filename)
        error("file \"$filename\" already exists (consider using `createfits!`)")
    end
    return createfits!(filename)
end
createfits(func::Function, filename::AbstractString; kwds...) = begin
    io = createfits(filename; kwds...)
    try
        func(io)
    finally
        #println("closing FITS file...")
        close(io)
    end
end

@doc @doc(createfits!) createfits

"""

```julia
readfits([T,] arg, ext=1) -> A
```

If `T` is unspecified, the above statement yields a pseudo-array `A` with the
contents of the FITS extension `ext` in `arg`.  Argument `arg` can be the name
of a FITS file or a FITS handle.  The optional HDU number, the first one by
default, must correspond to a FITS *FitsImage* extension.  The result is
indexable.  Using string index yields the value of the corresponding FITS
keyword in the header part of the HDU.  Any other indices are used to access
the contents of data part of the HDU (as a regular Julia array).

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
A = readfits("image.fits")         # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
A.BITPIX                           # idem
getfitscomment(A, "BITPIX")        # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
A.STUFF = 1                        # idem
setfitskey!(A, "STUFF", 3, "Blah") # idem with a comment
A["STUFF"] = (3, "Blah")           # idem with value-comment pair
A.STUFF = (3, "Blah")              # idem
arr = get(Array, A)               # get the data part (a regular Julia array)
hdr = get(FitsHeader, A)           # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
```

See also: [`openfits`](@ref).

"""
readfits(filename::AbstractString, args...; kwds...) =
    readfits(FitsImage, filename, args...; kwds...)

readfits(::Type{T}, filename::AbstractString, args...; kwds...) where {T} =
    openfits(filename) do io
        return readfits(T, io, args...; kwds...)
    end

readfits(fh::FITS, args...; kwds...) = readfits(FitsImage, fh, args...; kwds...)
readfits(::Type{T}, fh::FITS, ext::Union{Integer,AbstractString} = 1) where {T} =
    readfits(T, fh[ext])

# Read header from FITS ImageHDU.
readfits(::Type{FITSHeader}, hdu::HDU) = read_header(hdu)
readfits(::Type{FitsHeader}, hdu::HDU) = FitsHeader(read_header(hdu))

# Read header and data from FITS ImageHDU.
readfits(::Type{FitsImage}, hdu::HDU) =
    FitsImage(readfits(Array, hdu), readfits(FitsHeader, hdu))
readfits(::Type{FitsImage{T}}, hdu::HDU) where {T} =
    FitsImage(readfits(Array{T}, hdu), readfits(FitsHeader, hdu))
readfits(::Type{FitsImage{T,N}}, hdu::HDU) where {T,N} =
    FitsImage(readfits(Array{T,N}, hdu), readfits(FitsHeader, hdu))

# Read array from FITS ImageHDU.
readfits(::Type{Array}, hdu::ImageHDU) =
    read(hdu)
readfits(::Type{Array{T}}, hdu::ImageHDU) where {T} =
    convert(Array{T}, read(hdu))
readfits(::Type{Array{T,N}}, hdu::ImageHDU) where {T,N} =
    convert(Array{T,N}, read(hdu))

"""

```julia
writefits!(filename, args...; kwds...)
```

creates a new FITS file named `filename` with contents built from the provided
arguments and keywords.  If the file already exists, it is silently
overwritten.

Examples:

```julia
writefits!(filename, arr; KEY1 = VAL1, KEY2 = (VAL2, COM2), ...)
writefits!(filename, arr, hdr)
writefits!(filename, obj)
writefits!(filename, obj, (arr, hdr), ...)
```

"""
writefits!(filename::AbstractString, args...; kwds...) =
    createfits!(filename) do io
        writefits(io, args...; kwds...)
    end

"""

```julia
writefits(filename, args...; overwrite=false, kwds...)
```

behaves like `writefits!` but throws an error if `filename` already exits and
keyword `overwrite` is `false`.

This method can also be used to append a new HDU to a FITS file openned
for writing.  For instance:

```julia
writefits(fh, arr; KEY1 = VAL1, KEY2 = (VAL2, COM2), ...)
writefits(fh, arr, hdr)
writefits(fh, (arr, hdr))
writefits(fh, obj)
writefits(fh, obj, (arr, hdr), ...)
```

where `fh` is an instance of `FITSIO.FITS`.

"""
writefits(filename::AbstractString, args...; overwrite::Bool=false, kwds...) = begin
    if !overwrite && exists(filename)
        error("file \"$filename\" already exists (consider using `writefits!`)")
    end
    writefits!(filename, args...; kwds...)
end

writefits(io::FITS, arr::AbstractArray; kwds...) = write(io, arr, FitsHeader(; kwds...))
writefits(io::FITS, arg::Tuple{<:AbstractArray,<:FitsHeader}) = write(io, arg...)
writefits(io::FITS, arr::AbstractArray, hdr::FitsHeader) = write(io, arr, hdr)
writefits(io::FITS, obj::FitsImage) = write(io, obj)
writefits(io::FITS, args...) =
    for arg in args
        writefits(io, arg)
    end

Base.write(io::FITS, arr::AbstractArray, hdr::FitsHeader) =
    write(io, arr, header = get(FITSHeader, hdr))

Base.write(io::FITS, obj::FitsImage) =
    write(io, get(Array, obj), get(FitsHeader, obj))


# FIXME: The following constitute *type piracy*.
Base.findfirst(pred::Function, fh::FITS) = findnext(pred, fh, 1)
Base.findlast(pred::Function, fh::FITS) = findprev(pred, fh, length(fh))
Base.findnext(pred::Function, fh::FITS, i::Integer) = findnext(pred, fh, Int(i))
Base.findnext(pred::Function, fh::FITS, i::Int=1) = find(pred, fh, i:length(fh))
Base.findprev(pred::Function, fh::FITS, i::Integer) = findprev(pred, fh, Int(i))
Base.findprev(pred::Function, fh::FITS, i::Int=1) = find(pred, fh, i:-1:1)

"""

```julia
EasyFITS.find(pred, fh[, I = 1:length(fh)])
```

yields the HDU first index `i ∈ I` of FITS instance `fh` for which the
predicate `pred(fh[i])` returns `true`, or [`nothing`](@ref) if `pred(fh[i])`
returns `false` for all `i ∈ I`.

For instance:

```julia
i = EasyFITS.find(hdu -> EasyFITS.extname(hdu) == "CALIBRATION", fh)
if i === nothing
    # not found
    ...
else
    # found at index i
    hdu = fh[i]
end
```

"""
function find(pred::Function, fh::FITS,
              I::AbstractVector{Int} = 1:length(fh)) :: Union{Nothing,Int}
    for i in I
        if pred(fh[i])
            return i
        end
    end
    return nothing
end

"""

```julia
EasyFITS.getfile(arg) -> file
```

throws an error if the FITS file associated with `arg` is not open;
otherwise returns it.

```julia
EasyFITS.getfile(arg, ext) -> file
```

throws an error if the FITS file associated with `arg` is not open;
otherwise moves to HDU number `ext` and returns the FITS file.

Argument `arg` can be an instance of `FITSIO.FITS`, `FITSIO.HDU` or
`FITSIO.FITSFile`.

This method is considered as *low-level* and is not exported.

See also: [`EasyFITS.tryreadfitskeys`](@ref).

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

getfile(fh::FITS, ext::Integer) = getfile(fh.fitsfile, ext)

getfile(hdu::HDU) = getfile(hdu.fitsfile, hdu.ext)

"""

```julia
tryreadfitskey(src, T, key) -> val :: Union{T,Nothing}
```

attempts to read the value of FITS keywrod `key` keys in source `src` and
returns a value of type `T` or `nothing` if keyword is not found or if
parsing its value is unsuccessful.

See also: [`tryreadfitskeys`](@ref).

"""
tryreadfitskey(src::Union{FITSFile,FITS,HDU}, ::Type{T}, key::AbstractString) where T =
    tryreadfitskeys(src, T, (key,))

"""

```julia
tryreadfitskeys(src, T, keys) -> val :: Union{T,Nothing}
```

attempts to read the raw FITS keywords `keys` in given order in source
`src` and returns a value of type `T` or `nothing` if no keywords are found
or if parsing the value of the first matching keyword is unsuccessful.

See also: [`tryreadfitskey`](@ref).

"""
function tryreadfitskeys(src::Union{FITSFile,FITS,HDU}, ::Type{T},
                         keys::Union{AbstractArray{<:AbstractString},
                                     Tuple{Vararg{AbstractString}}}) where T
    return fits_try_read_keys(getfile(src), T, keys)
end

extname(hdu::HDU) =
    tryreadfitskey(getfile(hdu), String, "EXTNAME")

extversion(hdu::HDU) :: Int = begin
    version = tryreadfitskey(getfile(hdu), Int, "EXTVER")
    if version === nothing
        return 1
    end
    return version
end

hduname(hdu::HDU) =
    tryreadfitskey(getfile(hdu), String, "HDUNAME")

hduversion(hdu::HDU) :: Int = begin
    version = tryreadfitskey(getfile(hdu), Int, "HDUVER")
    if version === nothing
        return 1
    end
    return version
end

include("deprecate.jl")

end # module
