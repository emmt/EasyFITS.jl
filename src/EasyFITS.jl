module EasyFITS

export
    # Re-export from BaseFITS:
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

using TypeUtils

using BaseFITS
using BaseFITS:
    CardName,
    CardValue,
    CardComment,
    CardPair,
    FitsComplex,
    FitsInteger,
    FitsFloat,
    Undefined

using Base: @propagate_inbounds, string_index_err
using Base.Order: Ordering, Forward, Reverse
import Base: open, read, read!, write, merge, merge!, push!

let file = joinpath(@__DIR__, "..", "deps", "deps.jl")
    if !isfile(file)
        error("File \"$file\" does not exists. You may may generate it by:\n",
              "    using Pkg\n",
              "    Pkg.build(\"$(@__MODULE__)\")")
    end
    include(file)
end
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
