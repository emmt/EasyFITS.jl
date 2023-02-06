module EasyFITS

export
    # Re-export from FITSBase:
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
    openfits,
    readfits,
    write!,
    writefits,
    writefits!

using FITSBase
using FITSBase:
    FitsComplex,
    FitsInteger,
    FitsFloat

using Base: @propagate_inbounds, string_index_err
using Base.Order: Ordering, Forward, Reverse
import Base: open, read, read!, write

include("../deps/deps.jl")
include("types.jl")
include("SmallVectors.jl")
using .SmallVectors
include("utils.jl")
include("files.jl")
include("hdus.jl")
include("images.jl")
include("tables.jl")
include("init.jl")

end # module
