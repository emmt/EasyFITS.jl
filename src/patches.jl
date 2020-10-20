#
# patches.jl -
#
# Methods that should be part of FITSIO or CFITSIO.
#

# Should be in FITSIO.jl/src/header.jl
"""
    get(hdr, key, def)

yields value of keyword `key` in FITS header `hdr` or `def` if not found.

"""
Base.get(hdr::FITSHeader, key::String, def) =
    ((i = get(hdr.map, key, -1)) > 0 ? hdr.values[i] : def)

# Should be in CFITSIO.jl/src/CFITSIO.jl
Base.isopen(f::FITSFile) = (f.ptr != C_NULL)

# Should be in CFITSIO.jl/src/CFITSIO.jl, this is to avoid an infinite loop
# with something like type_from_bitpix(32) (because Int may not be Cint).
CFITSIO.type_from_bitpix(code::Integer) = type_from_bitpix(Val(Cint(code)))
