# Retrieve version number from the library.
function fits_get_version()
    version = round(Int, CFITSIO.fits_get_version(Ref{Cfloat}())*10_000.0)
    major = div(version, 10_000)
    minor, micro = divrem(version - 10_000*major, 100)
    return VersionNumber(major, minor, micro)
end

function __init__()
    # Version number according to header files used to generate this code.
    gen_version = VersionNumber(CFITSIO_MAJOR,CFITSIO_MINOR,CFITSIO_MICRO)

    # Version number of the library provided by the artifact.
    lib_version = fits_get_version()

    # Check for compatibility.
    (lib_version.major == gen_version.major && lib_version.minor ≥ gen_version.minor) || @warn """
`CFITSIO_jll` library has version $(lib_version) while the headers used to generate the code
of `CFITSIO.jl` have version $(gen_version). You should regenerate code in `CFITSIO.jl`
following instructions in `$(@__DIR__)/README.md`.
"""

end
