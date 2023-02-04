#------------------------------------------------------------------------------
# Open FITS files.

"""
    openfits(filename, mode="r"; kwds...) -> file

opens FITS file named `filename` with access `mode`. See [`FitsFile`](@ref) for
the different modes and keywords.

"""
function openfits(filename::AbstractString, mode::AbstractString = "r"; kwds...)
    return FitsFile(filename, mode; kwds...)
end

"""
    open(FitsFile, filename, mode="r"; kwds...) -> file

opens FITS file named `filename` with access `mode`. See [`FitsFile`](@ref) and
[`openfits`](@ref) for the different modes and keywords.

"""
function Base.open(::Type{FitsFile}, filename::AbstractString,
                   mode::AbstractString = "r"; kwds...)
    return openfits(filename, mode; kwds...)
end

"""
    FitsFile(filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block.

"""
function FitsFile(func::Function, filename::AbstractString,
                  mode::AbstractString = "r"; kwds...)
    file = FitsFile(filename, mode; kwds...)
    try
        func(file)
    finally
        close(file)
    end
end

"""
    openfits(filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block. See [`FitsFile`](@ref) for the different modes and keywords.

"""
function openfits(func::Function, filename::AbstractString,
                  mode::AbstractString = "r"; kwds...)
    return FitsFile(func, filename, mode; kwds...)
end

"""
    open(FitsFile, filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block. See [`FitsFile`](@ref) and [`openfits`](@ref) for the different
modes and keywords.

"""
function Base.open(func::Function, ::Type{FitsFile}, filename::AbstractString,
                   mode::AbstractString = "r"; kwds...)
    return openfits(func, filename, mode; kwds...)
end

#------------------------------------------------------------------------------
# Read FITS files.

"""
    readfits(R::Type=Array, filename; ext=1, col=nothing, extended=false) -> data::R

reads some data in extension `ext` (a Header Data Unit number or a name) in
FITS file `filename`. The data is returned as an object of type `R`. Array type
parameters may be specified in `R`. For example, specify `R = Array{Float32}`
to ensure that the result be a single precision floating-point array.

If the extension is a table, keyword `col` may be used to specify which
column(s) to read.

Specify keyword `extended = true` to use CFITSIO extended filename syntax.

"""
function readfits(filename::AbstractString, args...; kwds...)
    return readfits(Array, filename, args...; kwds...)
end

function readfits(R::Type, filename::AbstractString;
                  ext::Union{AbstractString,Integer} = 1,
                  col = nothing, kwds...)
    openfits(filename, "r"; kwds...) do file
        hdu = file[ext]
        if col === nothing
            return read(R, hdu)
        else
            hdu isa FitsTableHDU || error(
                "column(s) may only be specified for a FITS table extension")
            return read(R, hdu, col)
        end
    end
end

"""
    read(R::Type=Array, FitsFile, filename, args...; kwds...) -> data::R

reads some data in FITS file `filename`. See [`readfits`](@ref) for the meaning
of arguments and for possible keywords.

"""
function read(::Type{FitsFile}, filename::AbstractString, args...; kwds...)
    return readfits(filename, args...; kwds...)
end

function read(R::Type, ::Type{FitsFile}, filename::AbstractString, args...; kwds...)
    return readfits(R, filename, args...; kwds...)
end

"""
    readfits!(dest, filename, args...; kwds...) -> dest

overwrites destination `dest` with some data read from FITS file named
`filename`. This is more efficient but is similar to:

    copyto!(dest, readfits(typeof(dest), filename, args...; kwds...))

See [`readfits`](@ref) for the meaning of arguments and for possible keywords.

"""
function readfits!(dest, filename::AbstractString;
                   ext::Union{AbstractString,Integer} = 1,
                   col = nothing, kwds...)
    openfits(filename, "r"; kwds...) do file
        hdu = file[ext]
        if col === nothing
            read!(dest, hdu)
        else
            hdu isa FitsTableHDU || error(
                "column(s) may only be specified for a FITS table extension")
            read!(dest, hdu, col)
        end
    end
    return dest
end

"""
    read!(dest, FitsFile, filename, args...; kwds...) -> dest

overwrites destination `dest` with some data read from FITS file named
`filename`. See [`readfits!`](@ref) for the meaning of arguments and for
possible keywords.

"""
function read!(dest, ::Type{FitsFile}, filename::AbstractString, args...; kwds...)
    return readfits!(dest, filename, args...; kwds...)
end

#------------------------------------------------------------------------------
# Write FITS files.

"""
    writefits(filename, hdr, dat, args...; overwrite = false, kwds...)

creates a new FITS file named `filename` whose contents is specified by `hdr`,
`dat`, and `args...`. If the file already exists, the method fails unless
keyword `overwrite` is `true`. See [`FitsFile`](@ref) for other keywords that
may be specified when opening the file.

Arguments `hdr` and `dat` are the header and the data of a 1st Header Data Unit
(HDU) to write. Trailing arguments `args...` are headers and data of optional
additional HDUs.

See also [`writefits!`](@ref) and [`FitsFile`](@ref).

"""
function writefits(filename::AbstractString, args...;
                   overwrite::Bool = false, kwds...)
    openfits(filename, overwrite ? "w!" : "w"; kwds...) do file
        write(file, args...)
    end
    nothing
end

"""
    write(FitsFile, filename, args...; overwrite = false, kwds...)

creates a new FITS file named `filename` whose contents is specified by
`args...`. If the file already exists, the method fails unless keyword
`overwrite` is `true`. This method is equivalent to:

    writefits(filename, args...; overwrite = overwrite, kwds...)

See [`writefits`](@ref) for the meaning of `args...` and [`FitsFile`](@ref) for
other keywords that may be specified when opening the file.

"""
function write(::Type{FitsFile}, filename::AbstractString, args...; kwds...)
    return writefits(filename, args...; kwds...)
end

"""
    writefits!(filename, args...; kwds...)

creates a new FITS file named `filename` whose contents is specified by
`args...`. If the file already exists, it is (silently) overwritten. This
method is equivalent to:

    writefits(filename, args...; overwrite = true, kwds...)

See [`writefits`](@ref) for the meaning of `args...` and [`FitsFile`](@ref) for
other keywords that may be specified when opening the file.

"""
function writefits!(filename::AbstractString, args...; kwds...)
    return writefits(filename, args...; overwrite = true, kwds...)
end

"""
    write!(FitsFile, filename, args...; kwds...)

creates a new FITS file named `filename` whose contents is specified by
`args...`. If the file already exists, it is (silently) overwritten. This
method is equivalent to:

    writefits(filename, args...; overwrite = true, kwds...)

See [`writefits`](@ref) for the meaning of `args...` and [`FitsFile`](@ref) for
other keywords that may be specified when opening the file.

"""
function write!(::Type{FitsFile}, filename::AbstractString, args...; kwds...)
    return writefits(filename, args...; overwrite = true, kwds...)
end

# NOTE: It is assumed that:
#
#    write(file::FitsFile, header, data)
#
# is implemented elsewhere for different data types (image or table).

function write(file::FitsFile)
    # Nothing to do.
    return file
end

function write(file::FitsFile, hdr::Union{Nothing,Header},
               A::Union{ImageData,TableData}, args...)
    write(file, hdr, A)
    write(file, args...)
    return file
end

#------------------------------------------------------------------------------
# Interface to FITS files.

"""
    EasyFITS.get_handle(file::FitsFile)

yields the pointer to the opaque FITS file structure for `file`. It is the
caller responsibility to insure that the pointer is and remains valid as long
as it is needed.

!!! warning
    This function should never be directly called. When calling a function of
    the CFITSIO library (with `ccall` or equivalent), directly pass the
    `FitsFile` object so that (1) the validity of the pointer is checked and
    (2) the `FitsFile` object is preserved to not be garbage collected before
    the C function be called thus eliminating the risk of the file being closed
    and the pointer becoming invalid. `EasyFITS` simply achieves this by
    properly extending `Base.cconvert` and `Base.unsafe_convert`. In fact there
    are only 2 functions in `EasyFITS` which calls `get_handle`: `Base.isopen`
    which amounts to just checking whether the pointer is not null and, of
    course, `Base.unsafe_convert`.

"""
get_handle(file::FitsFile) = getfield(file, :handle)

# Extend Base.unsafe_convert to automatically extract and check the FITS file
# handle from a FitsFile object. This secures and simplifies calls to functions
# of the CFITSIO library. See `EasyFITS.get_handle` doc.
Base.unsafe_convert(::Type{Ptr{CFITSIO.fitsfile}}, file::FitsFile) =
    check(get_handle(file))

"""
    isopen(file::FitsFile)

returns whether `file` is open.

"""
Base.isopen(file::FitsFile) = !isnull(get_handle(file))

"""
    close(file::FitsFile)

closes the file associated with `file`.

"""
function Base.close(file::FitsFile)
    check(close_handle(file))
    nothing
end

# The following method is used to finalize or to close the object.
function close_handle(file::FitsFile)
    status = Ref{Status}(0)
    if isopen(file)
        CFITSIO.fits_close_file(file, status)
        setfield!(file, :handle, Ptr{CFITSIO.fitsfile}(0))
    end
    return status[]
end

"""
    pathof(file::FitsFile) -> str

yields the name of the FITS file associated with `file`.

"""
Base.pathof(file::FitsFile) = getfield(file, :path)

"""
    filemode(file::FitsFile)

yields `:r`, `:rw`, or `:w` depending whether `file` is open for reading, reading
and writing, or writing.

"""
Base.filemode(file::FitsFile) = getfield(file, :mode)

"""
    isreadable(file::FitsFile)

returns whether `file` is readable.

"""
Base.isreadable(file::FitsFile) = (filemode(file) !== :w) && isopen(file)

"""
    isreadonly(file::FitsFile)

returns whether `file` is read-only.

"""
Base.isreadonly(file::FitsFile) = (filemode(file) === :r) && isopen(file)

"""
    iswritable(file::FitsFile)

returns whether `file` is writable.

"""
Base.iswritable(file::FitsFile) = (filemode(file) !== :r) && isopen(file)

"""
    seek(file::FitsFile, n) -> type

moves to `n`-th HDU of FITS file `file` and returns an integer identifying the
type of the HDU:

* `FITS_IMAGE_HDU` if the `n`-th HDU contains an image.

* `FITS_BINARY_TABLE_HDU` if the `n`-th HDU contains a binary table.

* `FITS_ASCII_TABLE_HDU` if the `n`-th HDU contains an ASCII table.

* `FITS_ANY_HDU` if the `n`-th HDU is undefined.

An error is thrown if the file has been closed.

See also [`seekstart(::FitsFile)`](@ref), [`seekend(::FitsFile)`](@ref), and
[`position(::FitsFile)`](@ref).

"""
function Base.seek(file::FitsFile, i::Integer)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(file, i, type, Ref{Status}(0)))
    return Int(type[])
end

"""
    seekstart(file::FitsFile) -> type

moves to the first HDU of FITS file `file` and returns an integer identifying
the type of the HDU. See [`seek(::FitsFile)`](@ref).

"""
Base.seekstart(file::FitsFile) = seek(file, firstindex(file))

"""
    seekend(file::FitsFile) -> type

moves to the last HDU of FITS file `file` and returns an integer identifying
the type of the HDU. See [`seek(::FitsFile)`](@ref).

"""
Base.seekend(file::FitsFile) = seek(file, lastindex(file))

"""
    position(file::FitsFile) -> n

yields the current HDU number of FITS file `file`. An error is thrown if the
file has been closed. See [`seek(::FitsFile)`](@ref).

"""
function Base.position(file::FitsFile)
    num = Ref{Cint}()
    return Int(CFITSIO.fits_get_hdu_num(file, num))
end

function get_nhdus(file::FitsFile)
    num = Ref{Cint}()
    check(CFITSIO.fits_get_num_hdus(file, num, Ref{Status}(0)))
    return Int(num[])
end

"""
    flush(f::Union{FitsFile,FitsHDU})

flushes the internal data buffers of `f` to the associated output FITS file.

"""
Base.flush(f::Union{FitsFile,FitsHDU}) =
    check(CFITSIO.fits_flush_buffer(f, 0, Ref{Status}(0)))

# Implement abstract array API for FitsFile objects.
Base.length(file::FitsFile) = getfield(file, :nhdus)
Base.size(file::FitsFile) = (length(file),)
Base.axes(file::FitsFile) = (keys(file),)
Base.IndexStyle(::Type{FitsFile}) = IndexLinear()
Base.firstindex(::FitsFile) = 1
Base.lastindex(file::FitsFile) = length(file)
Base.keys(file::FitsFile) = Base.OneTo(length(file))
Base.getindex(file::FitsFile, i::Int) = FitsHDU(file, i)
function Base.getindex(file::FitsFile, str::AbstractString)
    hdu = findfirst(str, file)
    hdu === nothing && error("no FITS Header Data Unit named \"$str\"")
    return hdu
end

function Base.get(file::FitsFile, str::AbstractString, default)
    hdu = findfirst(str, file)
    hdu === nothing && return default
    return hdu
end
