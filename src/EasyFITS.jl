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
    createfits!,
    createfits,
    exists,
    getfitscomment,
    getfitsdata,
    getfitsheader,
    getfitskey,
    openfits,
    readfits,
    setfitskey!,
    tryreadfitskey,
    tryreadfitskeys

using FITSIO
using FITSIO.Libcfitsio
#import FITSIO.Libcfitsio: libcfitsio
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

```julia
EasyFITS.Image(arr, hdr)
```

yields a FITS Image that behaves like a Julia array when indexed by
integers of Cartesian indices and like a dictionary of FITS keywords when
indexed by strings.  Argument `arr` specifies the array contents of the
object and `hdr` specifies the keywords.

"""
struct Image{T,N} <: DenseArray{T,N}
    arr::Array{T,N}
    hdr::FITSHeader
end

Image(arr::AbstractArray) = Image(arr, header())
Image(arr::AbstractArray{T,N}, hdr::FITSHeader) where {T,N} =
    Image{T,N}(convert(Array{T,N}, arr), hdr)
Image{T}(init, dims::Integer...) where {T} =
    Image{T}(init, dims)
Image{T}(init, dims::NTuple{N,Integer}) where {T,N} =
    Image{T,N}(init, dims)
Image{T,N}(init, dims::Integer...) where {T,N} =
    Image{T,N}(init, dims)
Image{T,N}(init, dims::NTuple{N,Integer}) where {T,N} =
    Image{T,N}(Array{T,N}(init, dims), header())

# Make Image instances behave like arrays (indexing is considered later).
#Base.eltype(::Image{T,N}) where {T,N} = T # FIXME: not needed
#Base.ndims(::Image{T,N}) where {T,N} = N # FIXME: not needed
Base.length(A::Image) = length(getfitsdata(A))
Base.size(A::Image) = size(getfitsdata(A))
Base.size(A::Image, d) = size(getfitsdata(A), d)
Base.axes(A::Image) = axes(getfitsdata(A))
Base.axes(A::Image, d) = axes(getfitsdata(A), d)
@inline Base.axes1(A::Image) = axes1(getfitsdata(A))
Base.IndexStyle(::Type{<:Image}) = IndexLinear()
Base.similar(::Type{Image{T}}, dims::NTuple{N,Int}) where {T,N} =
    Image{T,N}(undef, dims)
Base.elsize(::Type{Image{T,N}}) where {T,N} = elsize(Array{T,N})
Base.sizeof(A::Image) = sizeof(getfitsdata(A))

# Make Image's efficient iterators.
@inline Base.iterate(A::Image, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

# FIXME: copyto!, convert, unsafe_convert, pointer, etc.

getfitsdata(A::Image) = A.arr
getfitsheader(dat::Image) = dat.hdr
getfitscomment(dat::Image, k) = FITSIO.get_comment(getfitsheader(dat), k)

FITSIO.get_comment(dat::Image, k) = getfitscomment(dat, k)

getindex(dat::Image, key::AbstractString) = getindex(getfitsheader(dat), key)
#getindex(dat::Image, inds...) = getindex(getfitsdata(dat), inds...)
@inline @propagate_inbounds getindex(A::Image, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds r = getindex(getfitsdata(A), i)
    return r
end

setindex!(dat::Image, val, key::AbstractString) =
    setfitskey!(getfitsheader(dat), key, val)
#setindex!(dat::Image, val, inds...) = setindex!(getfitsdata(dat), val, inds...)
@inline @propagate_inbounds setindex!(A::Image, x, i::Int) = begin
    @boundscheck checkbounds(A, i)
    @inbounds r = setindex!(getfitsdata(A), x, i)
    return r
end

@inline Base.checkbounds(A::Image, i::Int) =
    1 ≤ i ≤ length(A) || throw_boundserror(A, i)

keys(dat::Image) = keys(getfitsheader(dat))
nkeys(dat::Image) = nkeys(getfitsheader(dat))
nkeys(hdr::FITSHeader) = length(hdr)
haskey(dat::Image, key) = haskey(getfitsheader(dat), key)
getkey(dat::Image, key, def) = getkey(getfitsheader(dat), key, def)

"""

```julia
getfitskey(T, obj, key[, def]) -> val :: T
```

yields the value of FITS keyword `key` in `obj`.  Argument `obj` can be an
instance of `EasyFITS.Image` or `FITSIO.FITSHeader`. The returned value has
type `T`.  Optional argument `def` is to specify a default value if keyword
`key` is not found in `hdr`.  An error is thrown if keyword `key` is not
found in `hdr` and no default value is provided or if the value to be
returned has not a type that can be converted to `T`.

Call `getfistcomment` to retrieve the comment assiciated with a FITS
keyword.

"""
getfitskey(::Type{T}, dat::Image, args...) where {T} =
    getfitskey(T, getfitsheader(dat), args...)

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
setfitskey!(dst, key, val[, com])
```

set FITS keyword `key` in FITS header `dst` to the value `val`.  Argument `com`
is an optional comment; if it is not specified and `ky` already exists in
`dst`, its comment is preserved.

Argument `dst` can be an instance of `FITSIO.FITSHeader` or of `EasyFITS.Image`.

See also: [`readfits`](@ref).

"""
setfitskey!(dat::Image, args...) = setfitskey!(getfitsheader(dat), args...)
setfitskey!(hdr::FITSHeader, key::AbstractString, val) = begin
    setindex!(hdr, val, key)
    return nothing
end
setfitskey!(hdr::FITSHeader, key::AbstractString, val, com::AbstractString) = begin
    setindex!(hdr, val, key)
    set_comment!(hdr, key, com)
    return nothing
end

"""

```julia
EasyFITS.header()
```

yields an empty `FITSIO.FITSHeader`.

"""
header() = FITSHeader(String[], [], String[])
# FIXME: should be Base.empty(::Type{FITSHeader}) = ...

"""

```julia
openfits(filename) -> fh
```

yields a `FITS` handle to read the contents of the existing FITS file whose
name is `filename`.  If `filename` does not exist but `filename` does not end
with the `".gz"` extension and `"\$filename.gz"` does exist, then the
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

creates a new FITS file named `fillename` and return its handle.  If the file
already exists, it is (silently) overwritten.

The `createfits` method can be used instead to avoid that.  The syntax is identical:

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
readfits(arg, hdu=1) -> A
```

yields a pseudo-array `A` with the contents of the FITS HDU (*header data
unit*) `hdu` in `arg`.  Argument `arg` can be the name of a FITS file or a FITS
handle.  The optional HDU number, the first one by default, must correspond to
a FITS *Image* extension.  The result is indexable.  Using string index yields
the value of the corresponding FITS keyword in the header part of the HDU.  Any
other indices are used to access the contents of data part of the HDU (as a
regular Julia array).

Examples:

```julia
using EasyFITS
A = readfits("image.fits")         # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
getfitscomment(A, "BITPIX")        # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
setfitskey!(A, "STUFF", 3, "Blah") # idem with a comment
arr = getfitsdata(A)               # get the data part (a regular Julia array)
hdr = getfitsheader(A)             # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
```

See also: [`openfits`](@ref).

"""
readfits(filename::AbstractString, hdu::Integer=1) =
    openfits(filename) do io
        return readfits(io, hdu)
    end

readfits(fh::FITS, hdu::Integer=1) =
    Image(read(fh[hdu]), read_header(fh[hdu]))

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
