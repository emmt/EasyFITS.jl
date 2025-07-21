using Requires

"""
    EasyFITS.CFITSIO_VERSION

is the version number of the CFITSIO library for which `EasyFITS` has been built. When
`EasyFITS` is loaded, it is checked that the version of the CFITSIO library does match this
version.

"""
const CFITSIO_VERSION = VersionNumber(CFITSIO.CFITSIO_MAJOR,
                                      CFITSIO.CFITSIO_MINOR,
                                      CFITSIO.CFITSIO_MICRO)

function __init__()
    # Check version.
    version = round(Int, CFITSIO.fits_get_version(Ref{Float32}())*10_000.0)
    major = div(version, 10_000)
    minor, micro = divrem(version - 10_000*major, 100)
    lib_version = VersionNumber(major, minor, micro)

    CFITSIO_VERSION == lib_version || error(
        "`$(@__MODULE__)` has been built for CFITSIO version $CFITSIO_VERSION ",
        "while the version of the dynamic library is $lib_version. \n",
        "`$(@__MODULE__)` must be rebuilt, for example by:\n",
              "    using Pkg\n",
              "    Pkg.build(\"$(@__MODULE__)\")")

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
