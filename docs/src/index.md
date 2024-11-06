# Introduction

`EasyFITS` is a [Julia](https://julialang.org/) package for reading and writing files in
[FITS](https://fits.gsfc.nasa.gov/fits_standard.html) format, a *Flexible Image Transport
System*, widely used in astronomy.

The source code of `EasyFITS` is on [GitHub](https://github.com/emmt/EasyFITS.jl).


## Related software

- [`FITSHeaders`](https://github.com/emmt/FITSHeaders.jl) is used by `EasyFITS` to
  efficiently manage FITS header cards. You may read `FITSHeaders` documentation to learn
  how to deal with FITS header cards when using `EasyFITS`.

- [`CFITSIO`](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html) is the
  standard C library for managing FITS files and is used by `EasyFITS`.

- [`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) is another Julia package for
  reading/writing with FITS files. `FITSIO.jl` is also based on the `CFITSIO` C library.
  Compared to `FITSIO.jl`, `EasyFITS` attempts to be more intuitive and more type-stable.

- [`CFITSIO_jll`](https://github.com/JuliaBinaryWrappers/CFITSIO_jll.jl) is the Julia
   artifact providing the `CFITSIO` C library used by the `EasyFITS` and `FITSIO.jl`
   packages.


## Table of contents

```@contents
Pages = ["structure.md", "reading.md", "writing.md", "files.md","images.md", "tables.md", "library.md", "links.md"]
```
