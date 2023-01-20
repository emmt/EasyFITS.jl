"""
    FITSFile(filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block.

"""
function FITSFile(func::Function, filename::AbstractString,
                  mode::AbstractString = "r"; kwds...)
    file = FITSFile(filename, mode; kwds...)
    try
        func(file)
    finally
        close(file)
    end
end

"""
    open(FITSFile, filename, mode="r"; kwds...)

opens FITS file named `filename` with access `mode`. See [`FITSFile`](@ref) for
the different modes and keywords.

"""
function Base.open(::Type{FITSFile}, filename::AbstractString,
                   mode::AbstractString = "r"; kwds...)
    FITSFile(filename, mode; kwds...)
end

"""
    open(FITSFile, filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block. See [`FITSFile`](@ref) for the different modes and keywords.

"""
function Base.open(func::Function, ::Type{FITSFile}, filename::AbstractString,
                   mode::AbstractString = "r"; kwds...)
    FITSFile(func, filename, mode; kwds...)
end

"""
    openfits(filename, mode="r"; kwds...)

opens FITS file named `filename` with access `mode`. See [`FITSFile`](@ref) for
the different modes and keywords.

"""
function openfits(filename::AbstractString, mode::AbstractString = "r"; kwds...)
    open(FITSFile, filename, mode; kwds...)
end

"""
    openfits(filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of
the block. See [`FITSFile`](@ref) for the different modes and keywords.

"""
function openfits(func::Function, filename::AbstractString,
                  mode::AbstractString = "r"; kwds...)
    open(func, FITSFile, filename, mode; kwds...)
end

function readfits(filename::AbstractString, args...; kwds...)
    read(FITSFile, filename, args...; kwds...)
end

"""
    writefits(filename, hdr, dat, args...; overwrite = false, kwds...)

creates a new FITS file named `filename` and with contents specified by `hdr`,
`dat`, and `args...`. If the file already exists, the method fails unless
keyword `overwrite` is `true`. See [`FITSFile`](@ref) for other keywords that
may be specified when opening the file.

Arguments `hdr` and `dat` are the header and the data of a 1st Header Data Unit
(HDU) to write. Arguments `args...` denotes headers and data for optional
additional HDUs.

See also [`writefits!`](@ref) and [`FITSFile`](@ref).

"""
function writefits(filename::AbstractString, args...;
                   overwrite::Bool = false, kwds...)
    openfits(filename, overwrite ? "w!" : "w"; kwds...) do file
        write(file, args...)
    end
    nothing
end

function write(file::FITSFile)
    # Nothing to do.
    nothing # FIXME: return file?
end

function write(file::FITSFile, hdr::Union{Nothing,Header}, A::ImageData)
    # Write a FITS image extension.
    # FIXME: write(file, hdr, A)
end

function write(file::FITSFile, hdr::Union{Nothing,Header}, A::TableData)
    # Write a FITS binary table extension.
    # FIXME: write(file, hdr, A)
end

function write(file::FITSFile, hdr::Union{Nothing,Header},
               A::Union{ImageData,TableData}, args...)
    write(file, hdr, A)
    write(file, args...)
end

"""
    writefits!(filename, args...; kwds...)

creates a new FITS file named `filename` and with contents specified by
`args...`. If the file already exists, it is silently overwritten. This method
is equivalent to:

    writefits(filename, args...; overwrite = true, kwds...)

See also [`writefits`](@ref) and [`FITSFile`](@ref).

"""
function writefits!(filename::AbstractString, args...; kwds...)
    writefits(filename, args...; overwrite = true, kwds...)
end

"""
    write(FITSFile, filename, args...; overwrite = false, kwds...)

creates a new FITS file named `filename` and with contents specified `args...`.
If the file already exists, the method fails unless keyword `overwrite` is
`true`. See [`FITSFile`](@ref) for other keywords that may be specified when
opening the file.

Arguments `hdr` and `dat` are the header and the data of a 1st Header Data Unit
(HDU) to write. Arguments `args...` denotes headers and data for optional
additional HDUs.

See also [`writefits!`](@ref) and [`FITSFile`](@ref).

    write(FITSFile, filename, args...; kwds...)

creates a new FITS file named `filename` and with contents specified by
`args...`. If the file already exists, it is silently overwritten. This method
is equivalent to:

    writefits!(filename, args...; kwds...)

See also [`writefits!`](@ref) and [`FITSFile`](@ref).

"""
function write(::Type{FITSFile}, filename::AbstractString,
               A::AbstractArray,
               hdr::Union{Nothing,Header} = nothing;
               kwds...)
    return write(filename, hdr, A; kwds...)
end

"""
    write!(FITSFile, filename, args...; kwds...)

creates a new FITS file named `filename` and with contents specified by
`args...`. If the file already exists, it is silently overwritten. This method
is equivalent to:

    writefits!(filename, args...; kwds...)

See also [`writefits!`](@ref) and [`FITSFile`](@ref).

"""
function write!(::Type{FITSFile}, filename::AbstractString, args...; kwds...)
    writefits!(filename, args...; kwds...)
end

"""
    write!(FITSFile, filename, args...; kwds...)

creates a new FITS file named `filename` and with contents specified by
`args...`. This method is equivalent to:

    writefits(filename, args...; kwds...)

See also: [`writefits`](@ref) and [`FITSFile`](@ref).

    write(FITSFile, filename, arr, hdr=nothing)
    write(FITSFile, filename, hdr, arr)

create a new FITS file whose name is `filename` and writes array `arr` whith
header `hdr` into this file. If the file already exists, an error is thrown if
keyword `overwrite` is false; otherwise the file is silently overwritten. Order
of arguments `arr` and `hdr` is irrelevant.

Specify keyword `extended = true` to use CFITSIO extended filename syntax.

"""
function write(::Type{FITSFile}, filename::AbstractString, args...; kwds...)
    writefits(filename, args...; kwds...)
end

"""
    read(R::Type{<:Array}=Array, FITSFile, filename, ext=1; extended=false) -> arr::R

reads image extension `ext` (a header data unit number or name) in FITS file
`filename` and returns its contents as an array of type `R`.

    read(R::Type{<:Array}=Array, FITSFile, filename, ext, col; extended=false) -> arr::R

reads column `col` in table extension `ext` (a header data unit number or name)
of FITS file `filename` and returns its contents as an array of type `R`.

Array type parameters may be specified in `R`. For example, specify `R =
Array{Float32}` to ensure that the result be a single precision floating-point
array.

"""
function read(::Type{FITSFile}, filename::AbstractString,
              ext::Union{AbstractString,Integer} = 1; kwds...)
    return read(Array, filename, ext; kwds...)
end

function read(::Type{FITSFile}, filename::AbstractString,
              ext::Union{AbstractString,Integer},
              col::Union{AbstractString,Integer}; kwds...)
    return read(Array, filename, ext, col; kwds...)
end

function read(R::Type{<:Array}, ::Type{FITSFile}, filename::AbstractString,
              ext::Union{AbstractString,Integer} = 1; kwds...)
    open(filename, "r"; kwds...) do file
        hdu = file[ext]
        hdu isa FITSImageHDU || error("not a FITS image extension")
        return read(R, hdu)
    end
end

function read(R::Type{<:Array}, ::Type{FITSFile}, filename::AbstractString,
              ext::Union{AbstractString,Integer},
              col::Union{AbstractString,Integer}; kwds...)
    open(filename, "r"; kwds...) do file
        hdu = file[ext]
        hdu isa FITSTableHDU || error("not a FITS table extension")
        return read(R, hdu, col)
    end
end

function read!(arr::DenseArray, ::Type{FITSFile}, filename::AbstractString,
               ext::Union{AbstractString,Integer} = 1; kwds...)
    open(filename, "r"; kwds...) do file
        hdu = file[ext]
        hdu isa FITSImageHDU || error("not a FITS image extension")
        return read!(arr, hdu)
    end
end

function read!(arr::DenseArray, ::Type{FITSFile}, filename::AbstractString,
               ext::Union{AbstractString,Integer},
               col::Union{AbstractString,Integer}; kwds...)
    open(filename, "r"; kwds...) do file
        hdu = file[ext]
        hdu isa FITSTableHDU || error("not a FITS table extension")
        return read!(arr, hdu, col)
    end
end

"""
    EasyFITS.get_handle(file::FITSFile)

yields the pointer to the opaque FITS file structure for `file`. It is the
caller responsibility to insure that the pointer is and remains valid as long
as it is needed.

"""
get_handle(file::FITSFile) = getfield(file, :handle)

# Extend unsafe_convert to automatically extract and check the FITS file handle
# from a FITSFile object. This secures and simplifies calls to functions of the
# CFITSIO library.
Base.unsafe_convert(::Type{Ptr{CFITSIO.fitsfile}}, file::FITSFile) =
    check(get_handle(file))

"""
    isopen(file::FITSFile)

returns whether `file` is open.

"""
Base.isopen(file::FITSFile) = !isnull(get_handle(file))

"""
    close(file::FITSFile)

closes the file associated with `file`.

"""
function Base.close(file::FITSFile)
    check(close_handle(file))
    nothing
end

# The following method is used to finalize or to close the object.
function close_handle(file::FITSFile)
    status = Ref{Status}(0)
    ptr = get_handle(file)
    if ! isnull(ptr)
        CFITSIO.fits_close_file(ptr, status)
        setfield!(file, :handle, null(ptr))
    end
    return status[]
end

"""
    pathof(file::FITSFile) -> str

yields the name of the FITS file associated with `file`.

"""
Base.pathof(file::FITSFile) = getfield(file, :path)

"""
    filemode(file::FITSFile)

yields `:r`, `:rw`, or `:w` depending whether `file` is open for reading, reading
and writing, or writing.

"""
Base.filemode(file::FITSFile) = getfield(file, :mode)

"""
    isreadable(file::FITSFile)

returns whether `file` is readable.

"""
Base.isreadable(file::FITSFile) = (filemode(file) !== :w) && isopen(file)

"""
    isreadonly(file::FITSFile)

returns whether `file` is read-only.

"""
Base.isreadonly(file::FITSFile) = (filemode(file) === :r) && isopen(file)

"""
    iswritable(file::FITSFile)

returns whether `file` is writable.

"""
Base.iswritable(file::FITSFile) = (filemode(file) !== :r) && isopen(file)

"""
    seek(file::FITSFile, n) -> type

moves to `n`-th HDU of FITS file `file` and returns an integer identifying the
type of the HDU:

* `FITS_IMAGE_HDU` if the `n`-th HDU contains an image.

* `FITS_BINARY_TABLE_HDU` if the `n`-th HDU contains a binary table.

* `FITS_ASCII_TABLE_HDU` if the `n`-th HDU contains an ASCII table.

* `FITS_ANY_HDU` if the `n`-th HDU is undefined.

An error is thrown if the file has been closed.

See also [`seekstart(::FITSFile)`](@ref), [`seekend(::FITSFile)`](@ref), and
[`position(::FITSFile)`](@ref).

"""
function Base.seek(file::FITSFile, i::Integer)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(file, i, type, Ref{Status}(0)))
    return Int(type[])
end

"""
    seekstart(file::FITSFile) -> type

moves to the first HDU of FITS file `file` and returns an integer identifying
the type of the HDU. See [`seek(::FITSFile)`](@ref).

"""
Base.seekstart(file::FITSFile) = seek(file, firstindex(file))

"""
    seekend(file::FITSFile) -> type

moves to the last HDU of FITS file `file` and returns an integer identifying
the type of the HDU. See [`seek(::FITSFile)`](@ref).

"""
Base.seekend(file::FITSFile) = seek(file, lastindex(file))

"""
    position(file::FITSFile) -> n

yields the current HDU number of FITS file `file`. An error is thrown if the
file has been closed. See [`seek(::FITSFile)`](@ref).

"""
function Base.position(file::FITSFile)
    num = Ref{Cint}()
    return Int(CFITSIO.fits_get_hdu_num(file, num))
end

"""
    flush(f::Union{FITSFile,FITSHDU})

flushes the internal data buffers of `f` to the associated output FITS file.

"""
Base.flush(f::Union{FITSFile,FITSHDU}) =
    check(CFITSIO.fits_flush_buffer(f, 0, Ref{Status}(0)))

# Implement abstract array API for FITSFile objects.
Base.length(file::FITSFile) = getfield(file, :nhdus)
Base.size(file::FITSFile) = (length(file),)
Base.axes(file::FITSFile) = (keys(file),)
Base.IndexStyle(::Type{FITSFile}) = IndexLinear()
Base.firstindex(::FITSFile) = 1
Base.lastindex(file::FITSFile) = length(file)
Base.keys(file::FITSFile) = Base.OneTo(length(file))
Base.getindex(file::FITSFile, i::Int) = FITSHDU(file, i)
function Base.getindex(file::FITSFile, str::AbstractString)
    hdu = findfirst(str, file)
    hdu === nothing && error("no FITS Header Data Unit named \"$str\"")
    return hdu
end

function get_num_hdus(file::FITSFile)
    num = Ref{Cint}()
    check(CFITSIO.fits_get_num_hdus(file, num, Ref{Status}(0)))
    return Int(num[])
end

function Base.get(file::FITSFile, str::AbstractString, default)
    hdu = findfirst(str, file)
    hdu === nothing && return default
    return hdu
end
