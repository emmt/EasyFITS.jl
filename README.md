# Using FITS files made easier for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package is to facilitate the use of FITS files (widely used in
Astronomy).  It uses the great [FITSIO][fitsio-url] package which provides
a [Julia][julia-url] interface to the [CFITSIO][cfitsio-url] library.


## Usage

To use the EasyFITS package:

```julia
using EasyFITS
```

### High level methods

With EasyFITS, you can load a FITS Image (that is a multi-dimensional
array) and its header as an array-like object which can be indexed by
integers or Cartesian indices to get/set array values of by strings to
get/set keyword values.

The call

```julia
loadfits(arg, hdu=1) -> A
```

yields a pseudo-array `A` with the contents of the FITS HDU (*header data
unit*) `hdu` in `arg`.  Argument `arg` can be the name of a FITS file or a FITS
handle.  The optional HDU number, the first one by default, must correspond to
a FITS *Image* extension.  The result is indexable.  Using string index yields
the value of the corresponding FITS keyword in the header part of the HDU.  Any
other indices are used to access the contents of data part of the HDU (as a
regular Julia array).

Examples:

```julia
using EasyFITS
A = loadfits("image.fits")       # load the first HDU
A[2,3]                           # get value of data at indices (2,3)
A["BITPIX"]                      # get FITS bits per pixel
EasyFITS.getcomment(A, "BITPIX") # get the associated comment
A["STUFF"] = 1                   # set value of FITS keyword STUFF
setkey!(A, "STUFF", 3, "Blah")   # idem with a comment
arr = EasyFITS.getdata(A)        # get the data part (a regular Julia array)
hdr = EasyFITS.getheader(A)      # get the header part
EasyFITS.nkeys(A)                # get the number of keywords
EasyFITS.nkeys(hdr)              # get the number of keywords
keys(A)                          # get the list of keywords
keys(hdr)                        # get the list of keywords
```


### Lower level methods

To check whether a file `path` already exists in the file system, call:

```julia
exists(path) -> bool
```

To open an existing FITS file for reading (or updating), call:

```julia
openfits(path) -> fh
```

which yields a FITS handle `fh` to read/update the file contents.

To create a new FITS file for writing, call:

```julia
createfits(path; overwrite=false) -> fh
```

which yields a FITS handle `fh` to write the file contents.  If keyword
`overwrite` is `true`, the file is (silently) overwritten if it already
exists; otherwise (the default), an error is thrown if if the file already
exists.  A shortcut for creating the file even though it may already exists is to
call:

```julia
createfits!(path) -> fh
```

The do-block syntax is supported by `openfits`, `createfits` and
`createfits!` to automatically close the FITS file.  For instance:

```julia
openfits(path) do io

end
```

Also:

```julia
loadfits(path, hdu=1)

setkey!(dst, key, val[, com])
getkey(T, dat, key[, def]) -> val :: T

tryreadkey(src, T, key)
tryreadkeys(src, T, keys)

EasyFITS.getfile(arg [, ext])
EasyFITS.find(pred, )

```



[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/EasyFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/EasyFITS.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/EasyFITS.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/EasyFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/EasyFITS-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/EasyFITS.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/EasyFITS.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/EasyFITS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/EasyFITS.jl?branch=master

[fitsio-url]: https://github.com/JuliaAstro/FITSIO.jl
[julia-url]: http://julialang.org/
[cfitsio-url]: http://heasarc.gsfc.nasa.gov/fitsio/
