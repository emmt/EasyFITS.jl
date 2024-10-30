#----------------------------------------------------------------------------------------
# Open FITS files.

"""
    openfits(filename, mode="r"; kwds...) -> file

opens FITS file named `filename` with access `mode`. See [`FitsFile`](@ref) for the
different modes and keywords.

"""
function openfits(filename::AbstractString, mode::AbstractString = "r"; kwds...)
    return FitsFile(filename, mode; kwds...)
end

"""
    FitsFile(filename, mode="r"; kwds...) do file
        ... # use file
    end

do-block syntax to open a FITS file which is automatically closed at the end of the block.

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

do-block syntax to open a FITS file which is automatically closed at the end of the block.
See [`FitsFile`](@ref) for the different modes and keywords.

"""
function openfits(func::Function, filename::AbstractString,
                  mode::AbstractString = "r"; kwds...)
    return FitsFile(func, filename, mode; kwds...)
end

#----------------------------------------------------------------------------------------
# Read FITS files.

"""
    read(FitsHeader, filename; ext=1, kwds...) -> hdr::FitsHeader

yields the header of the `ext` extension of the FITS file `filename`. See
[`FitsFile`](@ref) for the possible keywords `kwds...`.

"""
function read(::Type{FitsHeader}, filename::AbstractString;
              ext::Union{AbstractString,Integer} = 1, kwds...)
    return FitsFile(filename, "r"; kwds...) do file; FitsHeader(file[ext]) end
end

"""
    readfits([R::Type,] filename, args...; ext=1, extended=false, kwds...) -> data

reads some data in extension `ext` (a Header Data Unit number or a name) in FITS file
`filename`. Specify keyword `extended = true` to use CFITSIO extended filename syntax.

If `R` is specified, the data is returned as an object of type `R`. Array type parameters
may be specified in `R`. For example, specify `R = Array{Float32}` to ensure that the
result be a single precision floating-point array.

If the extension is an image, `args...` specifies the ranges of pixels to read along the
dimensions. The default is to read all pixels.

If the extension is a table, `args...` consist in up to 2 arguments `cols` and `rows` to
select a subset of columns and of rows respectively. The default is to read all columns
and rows.

"""
function readfits(filename::AbstractString, args...; extended::Bool = false,
                  ext::Union{AbstractString,Integer} = 1, kwds...)
    openfits(filename, "r"; extended) do file
        return read(file[ext], args...; kwds...)
    end
end

function readfits(::Type{R}, filename::AbstractString, args...; extended::Bool = false,
                  ext::Union{AbstractString,Integer} = 1, kwds...)  where {R}
    openfits(filename, "r"; extended) do file
        return read(R, file[ext], args...; kwds...)
    end
end

"""
    readfits!(dest, filename, args...; kwds...) -> dest

overwrites destination `dest` with some data read from FITS file named `filename`. This is
more efficient but is similar to:

    copyto!(dest, readfits(typeof(dest), filename, args...; kwds...))

See [`readfits`](@ref) for the meaning of arguments and for possible keywords.

"""
function readfits!(dest, filename::AbstractString, args...; extended::Bool = false,
                   ext::Union{AbstractString,Integer} = 1, kwds...)
    openfits(filename, "r"; extended) do file
        return read!(dest, file[ext], args...; kwds...)
    end
end

#----------------------------------------------------------------------------------------
# Write FITS files.

"""
    writefits(filename, hdr, dat, args...; overwrite = false, kwds...)

creates a new FITS file named `filename` whose contents is specified by `hdr`, `dat`, and
`args...`. If the file already exists, the method fails unless keyword `overwrite` is
`true`. See [`FitsFile`](@ref) for other keywords that may be specified when opening the
file.

Arguments `hdr` and `dat` are the header and the data of a 1st Header Data Unit (HDU) to
write. Trailing arguments `args...` are headers and data of optional additional HDUs.

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
    writefits!(filename, args...; kwds...)

creates a new FITS file named `filename` whose contents is specified by `args...`. If the
file already exists, it is (silently) overwritten. This method is equivalent to:

    writefits(filename, args...; overwrite = true, kwds...)

See [`writefits`](@ref) for the meaning of `args...` and [`FitsFile`](@ref) for other
keywords that may be specified when opening the file.

"""
function writefits!(filename::AbstractString, args...; kwds...)
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

function write(file::FitsFile, header::OptionalHeader,
               data::Union{ImageData,TableData}, args...)
    write(file, header, data)
    write(file, args...)
    return file
end

# catch errors
@noinline function write(file::FitsFile,
                         @nospecialize(header::OptionalHeader),
                         @nospecialize(data::Any))
    error("no method to write FITS extension for data of type $(typeof(data))")
end

#----------------------------------------------------------------------------------------
# Interface to FITS files.

"""
    EasyFITS.get_handle(file::FitsFile)

yields the pointer to the opaque FITS file structure for `file`. It is the caller
responsibility to insure that the pointer is and remains valid as long as it is needed.

!!! warning
    This function should never be directly called. When calling a function of the CFITSIO
    library (with `ccall` or equivalent), directly pass the `FitsFile` object so that (1)
    the validity of the pointer is checked and (2) the `FitsFile` object is preserved to
    not be garbage collected before the C function be called thus eliminating the risk of
    the file being closed and the pointer becoming invalid. `EasyFITS` simply achieves
    this by properly extending `Base.cconvert` and `Base.unsafe_convert`. In fact there
    are only 2 functions in `EasyFITS` which calls `get_handle`: `Base.isopen` which
    amounts to just checking whether the pointer is not null and, of course,
    `Base.unsafe_convert`.

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

yields `:r`, `:rw`, or `:w` depending whether `file` is open for reading, reading and
writing, or writing.

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

moves to `n`-th HDU of FITS file `file` and returns an integer identifying the type of the
HDU:

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

moves to the first HDU of FITS file `file` and returns an integer identifying the type of
the HDU.

See also [`seek(::FitsFile)`](@ref).

"""
Base.seekstart(file::FitsFile) = seek(file, firstindex(file))

"""
    seekend(file::FitsFile) -> type

moves to the last HDU of FITS file `file` and returns an integer identifying the type of
the HDU.

See also [`seek(::FitsFile)`](@ref).

"""
Base.seekend(file::FitsFile) = seek(file, lastindex(file))

"""
    position(file::FitsFile) -> n

yields the current HDU number of FITS file `file`. An error is thrown if the file has been
closed.

See also [`seek(::FitsFile)`](@ref).

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
Base.length(file::FitsFile) = isopen(file) ? getfield(file, :nhdus) : 0
Base.size(file::FitsFile) = (length(file),)
Base.axes(file::FitsFile) = (keys(file),)
Base.IndexStyle(::Type{FitsFile}) = IndexLinear()
Base.firstindex(::FitsFile) = 1
Base.lastindex(file::FitsFile) = length(file)
Base.keys(file::FitsFile) = Base.OneTo(length(file))

Base.show(io::IO, file::FitsFile) = show(io, MIME"text/plain"(), file)
function Base.show(io::IO, mime::MIME"text/plain", file::FitsFile)
    if isopen(file)
        # Mimics `show` for abstract vectors.
        len = length(file)
        print(io, len, "-element FitsFile")
        if len > zero(len)
            println(io, ":")
            for i in 1:len
                print(io, " ")
                show(io, mime, file[i])
                i < len && println(io)
            end
        end
    else
        print(io, "closed FitsFile")
    end
    nothing
end

"""
    getindex(file::FitsFile, ext) -> hdu
    file[ext] -> hdu

yield the FITS Header Data Unit (HDU) of the FITS file `file` at `ext`, the FITS extension
number or name.

The returned object has the following read-only properties:

    hdu.file     # associated FITS file
    hdu.number   # HDU number (the index `i` above)
    hdu.type     # same as FitsHDUType(hdu)
    hdu.xtension # value of the XTENSION card (never nothing)
    hdu.extname  # value of the EXTNAME card or nothing
    hdu.hduname  # value of the HDUNAME card or nothing

"""
function Base.getindex(file::FitsFile, i::Integer)
    status = Ref{Status}(0)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(file, i, type, status))
    type = FitsHDUType(type[])
    if type == FITS_ASCII_TABLE_HDU
        return FitsTableHDU(BareBuild(), file, i, true)
    elseif type == FITS_BINARY_TABLE_HDU
        return FitsTableHDU(BareBuild(), file, i, false)
    elseif type == FITS_IMAGE_HDU
        bitpix = Ref{Cint}()
        check(CFITSIO.fits_get_img_equivtype(file, bitpix, status))
        ndims = Ref{Cint}()
        check(CFITSIO.fits_get_img_dim(file, ndims, status))
        N = as(Int, ndims[])
        T = type_from_bitpix(bitpix[])
        return FitsImageHDU{T,N}(BareBuild(), file, i)
    else
        return FitsAnyHDU(BareBuild(), file, i)
    end
end

function Base.getindex(file::FitsFile, str::AbstractString)
    i = findfirst(str, file)
    i === nothing && error("no FITS Header Data Unit named \"$str\"")
    return file[i]
end

function Base.get(file::FitsFile, i::Integer, def)
    i = as(keytype(file), i)
    return checkbounds(Bool, file, i) ? file[i] : def
end

function Base.get(file::FitsFile, str::AbstractString, def)
    i = findfirst(str, file)
    return i === nothing ? def : file[i]
end

"""
    nameof(hdu::FitsHDU) -> str

yields the name of the FITS header data unit `hdu`. The result is the value of the first
keyword of `"EXTNAME"` or `"HDUNAME"` which exists and has a string value. If none of
these keywords exist, the result is `hdu.xtension` which is the name of the FITS extension
of `hdu`, that is `"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"` depending on whether
`hdu` is an image, an ASCII table, a binary table, or anything else.

"""
function Base.nameof(hdu::FitsHDU)
    (str = hdu.hduname) === nothing || return str
    (str = hdu.extname) === nothing || return str
    return hdu.xtension
end

"""
    EasyFITS.is_named(hdu, pat) -> bool

yields whether pattern `pat` is equal to (in the FITS sense if `pat` is a string) or
matches (if `pat` is a regular expression) the extension of the FITS header data unit
`hdu`, or to the value of one of its `"EXTNAME"` or `"HDUNAME"` keywords. These are
respectively given by `hdu.xtension`, `hdu.extname`, or `hdu.hduname`.

This method is used as a predicate for the search methods `findfirst`, `findlast`,
`findnext`, and `findprev`.

The extension `hdu.xtension` is `"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"` depending
on whether `hdu` is an image, an ASCII table, a binary table, or anything else.

"""
is_named(hdu::FitsHDU, pat::Union{AbstractString,Regex}) =
    # Since a match only fails if no matching name is found, the order of the
    # tests is irrelevant. We therefore start with the costless ones.
    same_name(hdu.xtension, pat) ||
    same_name(hdu.hduname, pat) ||
    same_name(hdu.extname, pat)
is_named(pat::Union{AbstractString,Regex}) = Base.Fix2(is_named, pat)

# Compare HDU name with some pattern , HDU name may be `nothing` which can
# never be considered as a success.
same_name(name::Nothing, pat::Union{AbstractString,Regex}) = false
same_name(name::AbstractString, pat::AbstractString) = isequal(FitsLogic(), name, pat)
same_name(name::AbstractString, pat::Regex) = match(pat, name) !== nothing

for func in (:findfirst, :findlast)
    @eval Base.$func(pat::Union{AbstractString,Regex}, file::FitsFile) =
        $func(is_named(pat), file)
end

for func in (:findnext, :findprev)
    @eval Base.$func(pat::Union{AbstractString,Regex}, file::FitsFile, start::Integer) =
        $func(is_named(pat), file, start)
end

function Base.findfirst(f::Function, file::FitsFile)
    for i ∈ keys(file)
        f(file[i]) && return i
    end
    return nothing
end

function Base.findlast(f::Function, file::FitsFile)
    for i ∈ reverse(keys(file))
        f(file[i]) && return i
    end
    return nothing
end

function Base.findnext(f::Function, file::FitsFile, start::Integer)
    start = as(keytype(file), start)
    start < firstindex(file) && throw(BoundsError(file, start))
    for i ∈ start:lastindex(file)
        f(file[i]) && return i
    end
    return nothing
end

function Base.findprev(f::Function, file::FitsFile, start::Integer)
    start = as(keytype(file), start)
    start > lastindex(file) && throw(BoundsError(file, start))
    for i ∈ start:-1:firstindex(file)
        f(file[i]) && return i
    end
    return nothing
end

Base.haskey(file::FitsFile, ext::Integer) = 1 ≤ ext ≤ length(file)
Base.haskey(file::FitsFile, ext::AbstractString) = findfirst(ext, file) !== nothing

"""
    eachmatch(pat, file::FitsFile)

yields an iterator over the Header Data Units (HDUs) of FITS `file` matching pattern
`pat`. Pattern `pat` can be a string or a regular expression to be matched against the
name of the HDUs of `file` or a predicate function taking a HDU as argument and returning
whether it matches.

For example:

    for hdu in eachmatch(pat, file)
        ... # do something
    end

is a shortcut for:

    i = findfirst(pat, file)
    while i !== nothing
        hdu = file[i]
        ... # do something
        i = findnext(pat, file, i+1)
    end

while:

    for hdu in reverse(eachmatch(pat, file))
        ... # do something
    end

is equivalent to:

    i = findlast(pat, file)
    while i !== nothing
        hdu = file[i]
        ... # do something
        i = findprev(pat, file, i-1)
    end

"""
Base.eachmatch(pat, file::FitsFile) = FileIterator(pat, file)

struct FileIterator{O<:Ordering,P}
    pattern::P
    file::FitsFile
    FileIterator(ord::O, pat::P, file::FitsFile) where {O,P} =
        new{O,P}(pat, file)
end
FileIterator(pat, file::FitsFile) = FileIterator(Forward, pat, file)
FileIterator(ord::Ordering, pat::AbstractString, file::FitsFile) =
    FileIterator(ord, is_named(pat), file)
FileIterator(ord::Ordering, pat::Regex, file::FitsFile) =
    FileIterator(ord, is_named(pat), file)

Base.IteratorEltype(::Type{<:FileIterator}) = Base.HasEltype()
Base.eltype(::Type{<:FileIterator}) = FitsCard

Base.IteratorSize(::Type{<:FileIterator}) = Base.HasLength()
Base.length(iter::FileIterator) = length(iter.file)

Base.reverse(iter::FileIterator{typeof(Forward)}) =
    FileIterator(Reverse, iter.pattern, iter.file)
Base.reverse(iter::FileIterator{typeof(Reverse)}) =
    FileIterator(Forward, iter.pattern, iter.file)

# Iterate over entries in forward order.
function Base.iterate(iter::FileIterator{typeof(Forward)})
    j = findfirst(iter.pattern, iter.file)
    j === nothing ? nothing : ((@inbounds iter.file[j]), j+1)
end
function Base.iterate(iter::FileIterator{typeof(Forward)}, i::Int)
    j = findnext(iter.pattern, iter.file, i)
    j === nothing ? nothing : ((@inbounds iter.file[j]), j+1)
end

# Iterate over entries in reverse order.
function Base.iterate(iter::FileIterator{typeof(Reverse)})
    j = findlast(iter.pattern, iter.file)
    j === nothing ? nothing : ((@inbounds iter.file[j]), j-1)
end
function Base.iterate(iter::FileIterator{typeof(Reverse)}, i::Int)
    j = findprev(iter.pattern, iter.file, i)
    j === nothing ? nothing : ((@inbounds iter.file[j]), j-1)
end
