# Using FITS files made easier for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package is to facilitate the use of FITS files (widely used in
Astronomy) in [Julia][julia-url].  EasyFITS uses the great
[FITSIO][fitsio-url] package which provides a [Julia][julia-url] interface
to the [CFITSIO][cfitsio-url] library.


## Usage

To use the EasyFITS package:

```julia
using EasyFITS
```

### High level methods

EasyFITS provides objects of type `EasyFITS.Image` that, like FITS Image,
combine a data part which is a multi-dimensional array and a header part.
Such objects can be indexed by integers or Cartesian indices to get/set
array values of by strings to get/set keyword values.  Accessing these
objects as arrays, they should be as fast as oridnary Julia arrays.

To create an `EasyFITS.Image` from an existing array `arr`, call:

```julia
EasyFITS.Image(arr, hdr=EasyFITS.header()) -> A
```

where optional argument `hdr` is a FITS header.  By default, an empty
header is used.  If array `arr` is an `Array` instance, its contents is
shared by `A`; otherwise, `arr` is converted to an `Array` instance.  If a
header is provided its contents is also shared by `A`.  If you do not want
to share the contents of `arr`, just make a copy:

```julia
A = EasyFITS.Image(copy(arr))
```

It is also possible to create an `EasyFITS.Image` with a new data part of
type `T`, dimensions `dims` and an, initially, empty header:


```julia
A = EasyFITS.Image{T}(undef, dims)
```

To load a FITS Image extension as an instance of `EasyFITS.Image`, call:

```julia
readfits(arg, hdu=1) -> A
```

which yields a pseudo-array `A` with the contents of the FITS HDU (*header
data unit*) `hdu` in `arg`.  Argument `arg` can be the name of a FITS file
or a FITS handle.  The optional HDU number, the first one by default, must
correspond to a FITS *Image* extension.

Examples:

```julia
using EasyFITS
A = readfits("image.fits")         # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
getfitscomment(A, "BITPIX")        # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
setfitskey!(A, "STUFF", 3, "Blah") # idem with a comment
arr = getfitsdata(A)               # get the data part (a regular Julia array)
hdr = getfitsheader(A)             # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
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

which yields a `FITSIO.FITS` handle `fh` to read/update the contents of the
FITS file.  This is basically the same as calling `FITSIO.FITS(path)`
except that if `path` does not exist but `path` does not end with the
`".gz"` extension and `"$path.gz"` does exist, then the compressed file
`"$path.gz"` is open instead.

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
    # Read data from FITS handle io
    ...
end
```

Also:

```julia
setfitskey!(dst, key, val[, com])
getfitskey(T, dat, key[, def]) -> val :: T

tryreadfitskey(src, T, key)
tryreadfitskeys(src, T, keys)

EasyFITS.getfile(arg [, ext])
EasyFITS.find(pred, )
```

## Naming conventions

To avoid conflicts such as *type piracy*, all exported methods but `exists`
have the word **fits** embedded in their name.

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/EasyFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/EasyFITS.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/EasyFITS.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/EasyFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/EasyFITS-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/github/emmt/EasyFITS.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/emmt/EasyFITS.jl?branch=master

[codecov-img]: https://codecov.io/gh/emmt/EasyFITS.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/emmt/EasyFITS.jl

[fitsio-url]: https://github.com/JuliaAstro/FITSIO.jl
[julia-url]: http://julialang.org/
[cfitsio-url]: http://heasarc.gsfc.nasa.gov/fitsio/
