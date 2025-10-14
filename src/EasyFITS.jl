module EasyFITS

export
    # Re-export from FITSHeaders:
    @Fits_str,
    FitsKey,
    FitsCard,
    FitsCardType,
    FitsHeader,
    FITS_LOGICAL,
    FITS_INTEGER,
    FITS_FLOAT,
    FITS_STRING,
    FITS_COMPLEX,
    FITS_COMMENT,
    FITS_UNDEFINED,
    FITS_END,

    # FITS header data units.
    FitsHDU,
    FitsHDUType,
    FitsImageHDU,
    FitsTableHDU,
    FITS_ANY_HDU,
    FITS_ASCII_TABLE_HDU,
    FITS_BINARY_TABLE_HDU,
    FITS_IMAGE_HDU,

    # FITS exception, etc.
    FitsError,
    FitsLogic,

    # FITS file.
    FitsFile,
    readfits,
    readfits!,
    openfits,
    write!,
    writefits,
    writefits!

using TypeUtils

using FITSHeaders
using FITSHeaders:
    CardComment,
    CardName,
    CardPair,
    CardValue,
    FitsComplex,
    FitsFloat,
    FitsInteger,
    Undefined,
    is_comment,
    is_end,
    is_naxis,
    is_structural

using Base: @propagate_inbounds, string_index_err
using Base.Order: Ordering, Forward, Reverse
import Base: open, read, read!, write

if !isdefined(Base, :Memory)
    const Memory{T} = Vector{T}
end

include("CFITSIO.jl")
include("types.jl")
include("utils.jl")
include("files.jl")
include("hdus.jl")
include("images.jl")
include("tables.jl")
include("init.jl")
include("deprecated.jl")

end # module
