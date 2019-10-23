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
    loadfits,
    openfits,
    setkey!,
    tryreadkey,
    tryreadkeys

using FITSIO
using FITSIO.Libcfitsio
#import FITSIO.Libcfitsio: libcfitsio
using FITSIO: fits_try_read_keys

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

struct Data{T,N}
    hdr::FITSHeader
    arr::Array{T,N}
end

getheader(dat::Data) = dat.hdr
getdata(dat::Data) = dat.arr
getcomment(dat::Data, k) = FITSIO.get_comment(getheader(dat), k)

FITSIO.get_comment(dat::Data, k) = getcomment(dat, k)

getindex(dat::Data, key::AbstractString) = getindex(getheader(dat), key)
getindex(dat::Data, inds...) = getindex(getdata(dat), inds...)
setindex!(dat::Data, val, inds...) = setindex(getdata(dat), val, inds...)
setindex!(dat::Data, val, key::AbstractString) =
    setkey!(getheader(dat), key, val)
keys(dat::Data) = keys(getheader(dat))
nkeys(dat::Data) = nkeys(getheader(dat))
nkeys(hdr::FITSHeader) = length(hdr)
haskey(dat::Data, key) = haskey(getheader(dat), key)
getkey(dat::Data, key, def) = getkey(getheader(dat), key, def)

"""

```julia
getkey(T, dat, key[, def]) -> val :: T
```

yields the value of keyword `key` in FITS data `dat`.  The returned value
has type `T`.  Optional argument `def` is to specify a default value if
keyword `key` is not found in `hdr`.  An error is thrown if keyword `key`
is not found in `hdr` and no default value is provided or if the value to
be returned has not a type that can be converted to `T`.

"""
getkey(::Type{T}, dat::Data, args...) where {T} =
    _getkey(T, getheader(dat), args...)

"""

```julia
_getkey(T, hdr, key[, def]) -> val :: T
```

yields the value of keyword `key` in FITS header `hdr`.  The returned value
has type `T`.  Optional argument `def` is to specify a default value if
keyword `key` is not found in `hdr`.  An error is thrown if keyword `key`
is not found in `hdr` and no default value is provided or if the value to
be returned has not a type that can be converted to `T`.


!!! note
    This method is named `_getkey` because extending `getkey` for a type
    (`FITSHeader`) that is not defined in this package would be considered as
    *type piracy*.

""" _getkey

for S in (AbstractFloat, Integer, AbstractString, Bool)
    @eval begin
        function _getkey(::Type{T}, hdr::FITSHeader,
                         key::AbstractString, def) where {T<:$S}
            val = haskey(hdr, key) ? getindex(hdr, key) : def
            return _checkvalue(T, $S, val, key)
        end

        function _getkey(::Type{T}, hdr::FITSHeader,
                         key::AbstractString) where {T<:$S}
            haskey(hdr, key) || missing_keyword(key)
            return _checkvalue(T, $S, getindex(hdr, key), key)
        end
    end
end

function _checkvalue(::Type{T}, ::Type{S}, val,
                     key::AbstractString)::T where {S,T}
    isa(val, S) || bad_type(key)
    return T(val)
end

@noinline missing_keyword(key::AbstractString) = throw(KeyError(key))

@noinline bad_type(key::AbstractString) =
    error("bad type for FITS keyword \"$key\"")

"""

```julia
setkey!(dst, key, val[, com])
```

set FITS keyword `key` in FITS header `dst` to the value `val`.  Argument `com`
is an optional comment; if it is not specified and `ky` already exists in
`dst`, its comment is preserved.

Argument `dst` can be an instance of `FITSIO.FITSHeader` or of `EasyFITS.Data`.

See also: [`loadfits`](@ref).

"""
setkey!(dat::Data, args...) = setkey!(getheader(dat), args...)
setkey!(hdr::FITSHeader, key::AbstractString, val) = begin
    setindex!(hdr, val, key)
    return nothing
end
setkey!(hdr::FITSHeader, key::AbstractString, val, com::AbstractString) = begin
    setindex!(hdr, val, key)
    set_comment!(hdr, key, com)
    return nothing
end

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

See also: [`loadfits`](@ref), [`createfits`](@ref).

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

See also: [`loadfits`](@ref), [`createfits`](@ref).

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
loadfits(arg, hdu=1) -> A
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
A = loadfits("image.fits")       # load the first HDU
A[2,3]                           # get value of data at indices (2,3)
A["BITPIX"]                      # get FITS bits per pixel
EasyFITS.getcomment(A, "BITPIX") # get the associated comment
A["STUFF"] = 1                   # set value of FITS keyword STUFF
setkey!(A, "STUFF", 3, "Blah")   # idem with a comment
arr = EasyFITS.getdata(A)        # get the data part (a regular Julia array)
hdr = EasyFITS.getheader(A)      # get the header part
EasyFITS.nkeys(A)                # get the number of keywords
EasyFITS.nkeys(hdr)              # get the number of keywords
keys(A)                          # get the list of keywords
keys(hdr)                        # get the list of keywords
```

See also: [`openfits`](@ref).

"""
loadfits(filename::AbstractString, hdu::Integer=1) =
    openfits(filename) do io
        return loadfits(io, hdu)
    end

loadfits(fh::FITS, hdu::Integer=1) =
    Data(read_header(fh[hdu]), read(fh[hdu]))

# FIXME: The following constitute *type piracy*.
Base.findfirst(pred::Function, fh::FITS) = findnext(pred, fh, 1)
Base.findlast(pred::Function, fh::FITS) = findprev(pred, fh, length(fh))
Base.findnext(pred::Function, fh::FITS, i::Integer) = findnext(pred, fh, Int(i))
Base.findnext(pred::Function, fh::FITS, i::Int=1) = find(pred, fh, i:length(fh))
Base.findprev(pred::Function, fh::FITS, i::Integer) = findprev(pred, fh, Int(i))
Base.findprev(pred::Function, fh::FITS, i::Int=1) = find(pred, fh, length(fh):-1:i)

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

See also: [`EasyFITS.tryreadkeys`](@ref).

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
EasyFITS.tryreadkey(src, T, key) -> val :: Union{T,Nothing}
```

attempts to read the value of FITS keywrod `key` keys in source `src` and
returns a value of type `T` or `nothing` if keyword is not found or if
parsing its value is unsuccessful.

See also: [`EasyFITS.tryreadkeys`](@ref).

"""
tryreadkey(src::Union{FITSFile,FITS,HDU}, ::Type{T}, key::AbstractString) where T =
    tryreadkeys(src, T, (key,))

"""

```julia
EasyFITS.tryreadkeys(src, T, keys) -> val :: Union{T,Nothing}
```

attempts to read the raw FITS keywords `keys` in given order in source
`src` and returns a value of type `T` or `nothing` if no keywords are found
or if parsing the value of the first matching keyword is unsuccessful.

See also: [`EasyFITS.tryreadkey`](@ref).

"""
function tryreadkeys(src::Union{FITSFile,FITS,HDU}, ::Type{T},
                     keys::Union{AbstractArray{<:AbstractString},
                                 Tuple{Vararg{AbstractString}}}) where T
    return fits_try_read_keys(getfile(src), T, keys)
end

extname(hdu::HDU) =
    tryreadkey(getfile(hdu), String, "EXTNAME")

extversion(hdu::HDU) :: Int = begin
    version = tryreadkey(getfile(hdu), Int, "EXTVER")
    if version === nothing
        return 1
    end
    return version
end

hduname(hdu::HDU) =
    tryreadkey(getfile(hdu), String, "HDUNAME")

hduversion(hdu::HDU) :: Int = begin
    version = tryreadkey(getfile(hdu), Int, "HDUVER")
    if version === nothing
        return 1
    end
    return version
end

end # module
