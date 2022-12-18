Base.open(path::FitsFile, mode::AbstractString = "r"; kwds...) =
    FitsIO(path, mode; kwds...)
Base.open(func::Function, path::FitsFile, mode::AbstractString = "r"; kwds...) =
    FitsIO(func, path, mode; kwds...)

openfits(path::AbstractString, mode::AbstractString = "r"; kwds...) =
    open(FitsFile(path), mode; kwds...)
openfits(func::Function, path::AbstractString, mode::AbstractString = "r"; kwds...) =
    open(func, FitsFile(path), mode; kwds...)

readfits(path::AbstractString, ext::Union{Integer,AbstractString} = 1; kwds...) =
    read(FitsFile(path), ext; kwds...)

writefits(path::AbstractString; kwds...) = write(FitsFile(path), mode; kwds...)
writefits!(path::AbstractString; kwds...) = write!(FitsFile(path), mode; kwds...)

# Implement do-block syntax.
function FitsIO(f::Function, path::AbstractString, mode::AbstractString = "r"; kwds...)
    io = FitsIO(path, mode; kwds...)
    try
        return f(io)
    finally
        close(io)
    end
end

"""
    write!(path::FitsFile, args...; kwds...)

writes FITS file whose name is given by `path`. If the file already exists, it
is silently overwritten. This method is a shortcut for:

    write(path, args...; overwrite=true, kwds...)

"""
write!(path::FitsFile, args...; kwds...) =  write(path, args...; overwrite=true, kwds...)

"""
    write(path::FitsFile, arr, hdr=nothing)
    write(path::FitsFile, hdr, arr)

create a new FITS file whose name is given by `path` and writes array `arr`
whith header `hdr` into this file. If the file already exists, an error is
thrown if keyword `overwrite` is false; otherwise the file is silently
overwritten. Order of arguments `arr` and `hdr` is irrelevant.

Specify keyword `extended = true` to use CFITSIO extended filename syntax.

"""
function write(path::FitsFile,
               hdr::Union{Nothing,Header},
               A::AbstractArray;
               #bitpix::Integer = type_to_bitpix(eltype(A)),
               overwrite::Bool = false,
               kwds...)
    open(path, overwrite ? "w!" : "w"; kwds...) do io
        write(io, A, hdr)
    end
end

function write(path::FitsFile,
               A::AbstractArray,
               hdr::Union{Nothing,Header} = nothing;
               kwds...)
    return write(path, hdr, A; kwds...)
end

"""
    read(R::Type{<:Array}=Array, path::FitsFile, ext=1; extended=false) -> arr::R

reads image extension `ext` (a header data unit number or name) in FITS file
`path` and returns its contents as an array of type `R`.

    read(R::Type{<:Array}=Array, path::FitsFile, ext, col; extended=false) -> arr::R

reads column `col` in table extension `ext` (a header data unit number or name)
of FITS file `path` and returns its contents as an array of type `R`.

Array type parameters may be specified in `R`. For example, specify `R =
Array{Float32}` to ensure that the result be a single precision floating-point
array.

"""
function read(path::FitsFile,
              ext::Union{AbstractString,Integer} = 1; kwds...)
    return read(Array, path, ext; kwds...)
end

function read(path::FitsFile,
              ext::Union{AbstractString,Integer},
              col::Union{AbstractString,Integer}; kwds...)
    return read(Array, path, ext, col; kwds...)
end

function read(R::Type{<:Array}, path::FitsFile,
              ext::Union{AbstractString,Integer} = 1; kwds...)
    open(path, "r"; kwds...) do io
        hdu = io[ext]
        hdu isa FitsImageHDU || error("not a FITS image extension")
        return read(R, hdu)
    end
end

function read(R::Type{<:Array}, path::FitsFile,
              ext::Union{AbstractString,Integer},
              col::Union{AbstractString,Integer}; kwds...)
    open(path, "r"; kwds...) do io
        hdu = io[ext]
        hdu isa FitsTableHDU || error("not a FITS table extension")
        return read(R, hdu, col)
    end
end

function read!(arr::DenseArray, path::FitsFile,
               ext::Union{AbstractString,Integer} = 1; kwds...)
    open(path, "r"; kwds...) do io
        hdu = io[ext]
        hdu isa FitsImageHDU || error("not a FITS image extension")
        return read!(arr, hdu)
    end
end

function read!(arr::DenseArray, path::FitsFile,
               ext::Union{AbstractString,Integer},
               col::Union{AbstractString,Integer}; kwds...)
    open(path, "r"; kwds...) do io
        hdu = io[ext]
        hdu isa FitsTableHDU || error("not a FITS table extension")
        return read!(arr, hdu, col)
    end
end

"""
    pointer(io::FitsIO)

yields the pointer to the FITS file for `io`. It is the caller responsibility
to insure that the pointer is valid, for example, by checking that the file is
open with `isopen(io)`.

"""
Base.pointer(io::FitsIO) = getfield(io, :handle)

# Extend unsafe_convert to automatically extract and check the FITS file handle
# from a FitsIO object. This secures and simplifies calls to functions of the
# CFITSIO library.
Base.unsafe_convert(::Type{Ptr{CFITSIO.fitsfile}}, io::FitsIO) = check(pointer(io))

"""
    isopen(io::FitsIO)

returns whether `io` is open.

"""
Base.isopen(io::FitsIO) = !isnull(pointer(io))

"""
    close(io::FitsIO)

closes the file associated with `io`.

"""
Base.close(io::FitsIO) = (check(_close(io)); nothing)

function _close(io::FitsIO)
    status = Ref{Status}(0)
    ptr = getfield(io, :handle)
    if ! isnull(ptr)
        CFITSIO.fits_close_file(ptr, status)
        setfield!(io, :handle, null(ptr))
    end
    return status[]
end

"""
    pathof(io::FitsIO) -> str

yields the name of the FITS file associated with `io`.

"""
Base.pathof(io::FitsIO) = getfield(io, :path)

"""
    filemode(io::FitsIO)

yields `:r`, `:rw`, or `:w` depending whether `io` is open for reading, reading
and writing, or writing.

"""
Base.filemode(io::FitsIO) = getfield(io, :mode)

"""
    isreadable(io::FitsIO)

returns whether `io` is readable.

"""
Base.isreadable(io::FitsIO) = (filemode(io) !== :w) && isopen(io)

"""
    isreadonly(io::FitsIO)

returns whether `io` is read-only.

"""
Base.isreadonly(io::FitsIO) = (filemode(io) === :r) && isopen(io)

"""
    iswritable(io::FitsIO)

returns whether `io` is writable.

"""
Base.iswritable(io::FitsIO) = (filemode(io) !== :r) && isopen(io)

"""
    seek(io::FitsIO, n) -> type

moves to `n`-th HDU of FITS file `io` and returns an integer identifying the
type of the HDU:

* `FITS_IMAGE_HDU` if the `n`-th HDU contains an image.

* `FITS_BINARY_TABLE_HDU` if the `n`-th HDU contains a binary table.

* `FITS_ASCII_TABLE_HDU` if the `n`-th HDU contains an ASCII table.

* `FITS_ANY_HDU` if the `n`-th HDU is undefined.

An error is thrown if the file has been closed.

See also [`seekstart(::FitsIO)`](@ref), [`seekend(::FitsIO)`](@ref), and
[`position(::FitsIO)`](@ref).

"""
function Base.seek(io::FitsIO, i::Integer)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(io, i, type, Ref{Status}(0)))
    return Int(type[])
end

"""
    seekstart(io::FitsIO) -> type

moves to the first HDU of FITS file `io` and returns an integer identifying the
type of the HDU. See [`seek(::FitsIO)`](@ref).

"""
Base.seekstart(io::FitsIO) = seek(io, firstindex(io))

"""
    seekend(io::FitsIO) -> type

moves to the last HDU of FITS file `io` and returns an integer identifying the
type of the HDU. See [`seek(::FitsIO)`](@ref).

"""
Base.seekend(io::FitsIO) = seek(io, lastindex(io))

"""
    position(io::FitsIO) -> n

yields the current HDU number of FITS file `io`. An error is thrown if the file
has been closed. See [`seek(::FitsIO)`](@ref).

"""
function Base.position(io::FitsIO)
    num = Ref{Cint}()
    return Int(CFITSIO.fits_get_hdu_num(io, num))
end

"""
    flush(f::Union{FitsIO,FitsHDU})

flushes the internal data buffers of `f` to the associated output FITS file.

"""
Base.flush(f::Union{FitsIO,FitsHDU}) =
    check(CFITSIO.fits_flush_buffer(out, Ref{Status}(0)))

# Implement abstract array API for FitsIO objects.
function Base.length(io::FitsIO)
    ptr = pointer(io)
    isnull(ptr) && return 0
    num = Ref{Cint}()
    check(CFITSIO.fits_get_num_hdus(ptr, num, Ref{Status}(0)))
    return Int(num[])
end
Base.size(io::FitsIO) = (length(io),)
Base.axes(io::FitsIO) = (keys(io),)
Base.IndexStyle(::Type{FitsIO}) = IndexLinear()
Base.firstindex(::FitsIO) = 1
Base.lastindex(io::FitsIO) = length(io)
Base.keys(io::FitsIO) = Base.OneTo(length(io))
Base.getindex(io::FitsIO, i::Int) = FitsHDU(io, i)
function Base.getindex(io::FitsIO, s::AbstractString)
    hdu = findfirst(s, io)
    hdu === nothing && error("no FITS header data unit named \"$s\"")
    return hdu
end

function Base.get(io::FitsIO, s::AbstractString, default)
    hdu = findfirst(s, io)
    hdu === nothing && return default
    return hdu
end

"""
    FitsFile(path)
    fits"path"

yield a decorated string that represents a FITS filename. This is useful to
have methods `open`, `read`, `write`, and `write!` behave specifically for FITS
files.

"""
FitsFile(x::FitsFile) = x

macro fits_str(str)
    :(FitsFile($str))
end

Base.String(x::FitsFile) = getfield(x, :path)
Base.string(x::FitsFile) = String(x)
Base.string(io::IO, x::FitsFile) = string(io, string(x))
Base.convert(::Type{String}, x::FitsFile) = String(x)
Base.show(io::IO, ::MIME"text/plain", x::FitsFile) = show(io, x)
Base.show(io::IO, x::FitsFile) = begin
    print(io, "fits")
    show(io, String(x))
    print(io, '"')
end

# Extend AbstractString API (see start of base/strings/basic.jl) plus methods
# specialized for String (see base/strings/string.jl).
Base.ncodeunits(s::FitsFile) = ncodeunits(String(s))
Base.codeunit(s::FitsFile) = codeunit(String(s))
Base.codeunit(s::FitsFile, i::Integer) = codeunit(String(s), i)
Base.isvalid(s::FitsFile, i::Integer) = isvalid(String(s), i)
Base.length(s::FitsFile) = length(String(s))
Base.pointer(s::FitsFile) = pointer(String(s))
Base.pointer(s::FitsFile, i::Integer) = pointer(String(s), i)
Base.iterate(s::FitsFile) = iterate(String(s))
Base.iterate(s::FitsFile, i::Int) = iterate(String(s), i)
@inline @propagate_inbounds Base.thisind(s::FitsFile, i::Int) = thisind(String(s), i)
@inline @propagate_inbounds Base.nextind(s::FitsFile, i::Int) = nextind(String(s), i)
@inline @propagate_inbounds Base.getindex(x::FitsFile, i::Int) = getindex(String(x), i)
Base.getindex(s::FitsFile, r::AbstractUnitRange{<:Integer}) = getindex(String(s), r)
Base.isascii(s::FitsFile) = isascii(String(s))
Base.cmp(a::FitsFile, b::AbstractString) = cmp(String(a), b)
Base.cmp(a::AbstractString, b::FitsFile) = cmp(a, String(b))
Base.cmp(a::FitsFile, b::FitsFile) = cmp(String(a), String(b))
Base.:(==)(a::FitsFile, b::AbstractString) = String(a) == b
Base.:(==)(a::AbstractString, b::FitsFile) = a == String(b)
Base.:(==)(a::FitsFile, b::FitsFile) = String(a) == String(b)
