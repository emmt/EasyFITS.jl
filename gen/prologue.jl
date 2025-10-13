using CFITSIO_jll: libcfitsio

struct Status
    code::Cint
end

fits_open_file(A, B, C, D) = ffopentest(CFITSIO_SONAME, A, B, C, D)
