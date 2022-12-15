using Requires

"""
    EasyFits.library_version()

yields the version number of the CFITSIO library.

"""
library_version() = LIBRARY_VERSION[]

const LIBRARY_VERSION = Ref{VersionNumber}()

function __init__()
    # Setup version.
    version = round(Int, CFITSIO.fits_get_version(Ref{Float32}())*10_000.0)
    major = div(version, 10_000)
    minor, micro = divrem(version - 10_000*major, 100)
    LIBRARY_VERSION[] = VersionNumber(major, minor, micro)

    # Extend DataFrames when this package is loaded.
    @require DataFrames="a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
        """
            nrow(hdu::FitsTableHDU)

        yields the number of rows of the FITS table extension in `hdu`.

        """
        DataFrames.nrow(hdu::FitsTableHDU) = hdu.nrows

        """
            ncol(hdu::FitsTableHDU)

        yields the number of columns of the FITS table extension in `hdu`.

    """
        DataFrames.ncol(hdu::FitsTableHDU) = hdu.ncols
    end
end
