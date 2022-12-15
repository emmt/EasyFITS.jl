module EasyFITS

export
    # FITS filename.
    @fits_str,
    FitsFile,

    # FITS header data units.
    FitsHDU,
    FitsHDUType,
    FitsImageHDU,
    FitsTableHDU,
    FITS_ANY_HDU,
    FITS_ASCII_TABLE_HDU,
    FITS_BINARY_TABLE_HDU,
    FITS_IMAGE_HDU,

    # FITS header cards.
    FitsCard,
    FitsCardType,
    FITS_UNKNOWN,
    FITS_EMPTY,
    FITS_LOGICAL,
    FITS_INTEGER,
    FITS_FLOAT,
    FITS_COMPLEX,
    FITS_STRING,
    FITS_COMMENT,

    # FITS exception, etc.
    FitsError,
    FitsLogic,

    # FITS i/o.
    FitsIO,
    openfits,
    readfits,
    write!,
    writefits,
    writefits!

using Base: @propagate_inbounds, string_index_err
import Base: open, read, read!, write

include("../deps/deps.jl")
include("types.jl")
include("utils.jl")
include("files.jl")
include("hdus.jl")
include("cards.jl")
include("images.jl")
include("tables.jl")
include("init.jl")

function Base.rm(io::FitsIO)
    status = Ref{Status}(0)
    ptr = check(getfield(io, :handle)) # FIXME: check pointer
    CFITSIO.fits_delete_file(ptr, status)
    setfield!(io, :handle, null(ptr)) # FIXME: close instead
    check(status)
    nothing
end

function read_keyword(io::FitsIO)
    status = Ref{Status}(0)
    existing = Ref{Cint}()
    remaining = Ref{Cint}()
    check(CFITSIO.fits_get_keyname(io, existing, remaining, status))
    return (Int(existing[]), Int(remaining[]))
end

function get_card_name(card::Union{AbstractString,Vector{UInt8}})
    key = Vector{UInt8}(undef, CFITSIO.FLEN_KEYWORD)
    len = Ref{Cint}()
    check(CFITSIO.fits_get_keyname(card, pointer(key), len, Ref{Status}(0)))
    return key
    return String(resize!(key, len[])) # NOTE: a bit faster
    key[len[] + 1] = 0
    return unsafe_string(pointer(key))
end

function get_card_value(card::Union{AbstractString,Vector{UInt8}})
    val = Vector{UInt8}(undef, CFITSIO.FLEN_VALUE)
    com = Vector{UInt8}(undef, CFITSIO.FLEN_COMMENT)
    check(CFITSIO.fits_parse_value(card, pointer(val), pointer(com), Ref{Status}(0)))
    return (unsafe_string(pointer(val)), unsafe_string(pointer(com)))
end

# Yields type of value as one of: 'C', 'L', 'I', 'F', or 'X', for character
# string, logical, integer, floating point, or complex, respectively.
function get_value_type(value::Union{AbstractString,Vector{UInt8}})
    type = Ref{Cchar}()
    check(CFITSIO.fits_get_keytype(value, type, Ref{Status}(0)))
    return Char(type[])
end

end # module
