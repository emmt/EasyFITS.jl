# Using FITS files made easier for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package is to facilitate the use of FITS files (widely used in Astronomy)
in [Julia][julia-url].  EasyFITS uses the great [FITSIO][fitsio-url] package
which provides a [Julia][julia-url] interface to the [CFITSIO][cfitsio-url]
library.


## Usage

To use the EasyFITS package:

```julia
using EasyFITS
```

### High level methods

EasyFITS provides objects of type `FitsImage` that, like FITS Image,
combine a data part which is a multi-dimensional array and a header part.  Such
objects can be indexed by integers or Cartesian indices to get/set array values
of by strings to get/set keyword values.  Accessing these objects as arrays,
they should be as fast as oridnary Julia arrays.

To create an `FitsImage` from an existing array `arr`, call:

```julia
FitsImage(arr, hdr=EasyFITS.header()) -> A
```

where optional argument `hdr` is a FITS header.  By default, an empty header is
used.  If array `arr` is an `Array` instance, its contents is shared by `A`;
otherwise, `arr` is converted to an `Array` instance.  If a header is provided
its contents is also shared by `A`.  If you do not want to share the contents
of `arr`, just make a copy:

```julia
A = FitsImage(copy(arr))
```

It is also possible to create an `FitsImage` with a new data part of type
`T`, dimensions `dims` and an, initially, empty header:


```julia
A = FitsImage{T}(undef, dims)
```

To load a FITS Image extension as an instance of `FitsImage`, call:

```julia
readfits(arg, ext=1) -> A
```

which yields a pseudo-array `A` with the contents of the FITS extension `ext`
in `arg`.  Argument `arg` can be the name of a FITS file or a FITS handle.  The
optional extension can be a HDU number (the first one by default), or an
extension name.  It must correspond to a FITS *Image* extension.

Examples:

```julia
using EasyFITS
A = readfits("image.fits")         # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
A.BITPIX                           # idem
get(FitsComment, A, "BITPIX")      # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
A.STUFF = 1                        # idem
setfitskey!(A, "STUFF", 3, "Blah") # idem with a comment
A["STUFF"] = (3, "Blah")           # idem with value-comment pair
A.STUFF = (3, "Blah")              # idem
arr = get(Array, A)                # get the data part (a regular Julia array)
hdr = get(FitsHeader, A)           # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
```

It is also possible to specify a type as the first argument of `readfits`
to constrain the type of the result.  For instance:

```julia
using EasyFITS
readfits("data.fits")                 # load the first array and header
readfits(FitsImage, "data.fits")      # idem
readfits(FitsHeader, "data.fits")     # reads only the header part
readfits(Array, "data.fits")          # only load the array part (as a regular array)
readfits(FitsImage{T}, "data.fits")   # yields pseudo-array with elements of type T
readfits(Array{T}, "data.fits")       # yields regular array with elements of type T
readfits(FitsImage{T,N}, "data.fits") # yields N-dimensional pseudo-array with elements of type T
readfits(Array{T,N}, "data.fits")     # yields N-dimensional regular array with elements of type T
```

Note that the result of `readfits(FitsHeader,"data.fits")` can be indexed by
strings to access FITS keywords and implements the `obj.key` syntax.


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
FITS file.  This is basically the same as calling `FITSIO.FITS(path)` except
that if `path` does not exist but `path` does not end with the `".gz"`
extension and `"$path.gz"` does exist, then the compressed file `"$path.gz"` is
open instead.

To create a new FITS file for writing, call:

```julia
createfits(path; overwrite=false) -> fh
```

which yields a FITS handle `fh` to write the file contents.  If keyword
`overwrite` is `true`, the file is (silently) overwritten if it already exists;
otherwise (the default), an error is thrown if if the file already exists.  A
shortcut for creating the file even though it may already exists is to call:

```julia
createfits!(path) -> fh
```

The do-block syntax is supported by `openfits`, `createfits` and `createfits!`
to automatically close the FITS file.  For instance:

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
