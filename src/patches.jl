#
# patches.jl -
#
# Methods that should be part of FITSIO or CFITSIO.
#

using Pkg

package_version(name) =
    filter(x-> x.second.name == name, Pkg.dependencies()) |> x -> first(x)[2].version

if package_version("FITSIO") < v"0.17.0"
    # Should be in FITSIO.jl/src/header.jl
    """
        get(hdr, key, def)

    yields value of keyword `key` in FITS header `hdr` or `def` if not found.

    """
    Base.get(hdr::FITSHeader, key::String, def) =
        ((i = get(hdr.map, key, -1)) > 0 ? hdr.values[i] : def)
end

# Should be in CFITSIO.jl/src/CFITSIO.jl
Base.isopen(f::FITSFile) = (f.ptr != C_NULL)
