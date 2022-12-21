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
    FITS_UNDEFINED,
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

end # module
