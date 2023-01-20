module EasyFITS

export
    # Re-export from FITSCards:
    @FITS_str,
    FITSKey,
    FITSCard,
    FITSCardType,
    FITS_LOGICAL,
    FITS_INTEGER,
    FITS_FLOAT,
    FITS_STRING,
    FITS_COMPLEX,
    FITS_COMMENT,
    FITS_UNDEFINED,
    FITS_END,

    # FITS header data units.
    FITSHDU,
    FITSHDUType,
    FITSImageHDU,
    FITSTableHDU,
    FITS_ANY_HDU,
    FITS_ASCII_TABLE_HDU,
    FITS_BINARY_TABLE_HDU,
    FITS_IMAGE_HDU,

    # FITS exception, etc.
    FITSError,
    FITSLogic,

    # FITS file.
    FITSFile,
    openfits,
    readfits,
    write!,
    writefits,
    writefits!

using FITSCards

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
