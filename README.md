# Using FITS files made easier for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package is to facilitate the use of FITS files (widely used in Astronomy)
in [Julia][julia-url].  EasyFITS uses [FITSIO.jl][fitsio-url]
[CFITSIO.jl][fitsio-url] packages which provide a [Julia][julia-url] interface
to the [CFITSIO][cfitsio-url] library.


## Usage

To use the EasyFITS package:

```julia
using EasyFITS
```

This imports several data types (all prefixed with `Fits`) and a few methods
(`exists` and `write!`).  In principle, with EasyFITS, you do not need to
import `FITSIO` (except to deal with FITS tables).


## EasyFITS objects

### FITS Image pseudo-arrays

EasyFITS provides objects of type `FitsImage` that, like FITS Image extensions,
combine a data part which is a multi-dimensional array and a header part.  Such
objects can be indexed by integers or Cartesian indices to get/set array values
of by strings to get/set keyword values.  Accessing these objects as arrays,
they should be as fast as ordinary Julia arrays.

To create an `FitsImage` from an existing array `arr`, call:

```julia
FitsImage(arr, hdr=FitsHeader()) -> A
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

### FITS header objects

A FITS header instance is created by one of:

```julia
hdr = FitsHeader(; key1=val1, key2=val2, ...)
hdr = FitsHeader("key1" => val1, "key2" => val2, ...)
```

showing that the initial header contents can be specified by keywords or
key-value pairs.  To avoid ambiguities the two styles cannot be mixed.  Using
either of the above syntax is a matter of taste.  Julia keyword syntax is
lighter, but FITS keywords with spaces (like the `"HIERARCH ..."` ones) are
easier to specify with key-value pairs.  Remember that, by convention, FITS
keywords are in uppercase letters.  Values `val2`, `val1` etc. can be 2-tuples
providing the value and the comment of the FITS keyword.

An instance of `FitsHeader`, say `hdr`, implements indexation by keywords, as
with `hdr[key]`, and `hdr.key` syntax.  The two styles are usable to retrieve
or to set the value of a keyword.  To retrieve just the comment part:

```julia
get(FitsComment, hdr, key)
```

This also works for an instance of `FitsImage`.

To retrieve the whole header from a given HDU of a `FitsIO` instance:

```julia
FitsHeader(io, ext=1)
```


### FITS Input/Output object

An instance of `FitsIO` is a wrapper around a FITS file open for reading or for
writing:

```julia
io = FitsIO(path, mode)
```

where `path` is the name of the FITS file and `mode` can be:

- `"r"` or `"r+"` to read or append to the contents of and existing FITS file.
  If file `path` does not exist but `path` does not end with the `".gz"`
  extension and `"\$path.gz"` does exist, then the compressed file
  `"$path.gz"` is open instead.

- `"w"` to create a new FITS file named `path` that must not already exist.  An
  error is thrown if the file already exists.

- `"w!"` to open FITS file named `path` for writing.  If the file already
  exists, it is (silently) overwritten.

Call `close(io)` to close the FITS file associated with the `FitsIO` instance
`io`.  Call `isopen(io)` to check whether the FITS file associated with the
`FitsIO` instance `io` is open.  Closing the FITS file is automatically done,
if needed, when the instance is garbage collected.

The do-block syntax is supported to automatically close the FITS file:

```julia
FitsIO(filename, mode="r") do io
    # use FITS handle io
    ...
end
```

An instance of `FitsIO` is a collection of *Header Data Units* (HDU) and
implements indexation and iteration.  Assuming `io` is a `FitsIO` object, then:

- `io[i]` yields the `i`-th HDU.

- `length(io)` yields the number of HDUs.

- `io[name]` or `io[name,vers]` yields the HDU whose `EXTNAME` (or `HDUNAME`)
  keyword is equal to `name` (a string) and, optionally, whose `EXTVER` (or
  `HDUVER`) keyword is equal to `vers` (an integer).

- You can do `for hdu in io; ...; end` to iterate through all HDU's of `io`.

- Methods `findfirst(p,io)`, `findlast(p,io)`, `findnext(p,io,i)` and
  `findprev(p,io,i)` can be used on `FitsIO` object `io` to search for a
  specific HDU.  These methods test each HDU (starting at initial index `i` for
  `findnext` and `findprev`) with the predicate function `p` (called with a
  `FitsHDU` argument) and return the index (an integer) of the first HDU for
  which the predicate yields `true` or `nothing` if this never occurs.
  These methods are build upon `EasyFITS.find` which may be directly
  called.


## Reading data from FITS files

To load a FITS Image extension as an instance of `FitsImage`, call:

```julia
read(FitsImage, src, args...; ext=1) -> A
```

which yields a pseudo-array `A` with the contents of the extension `ext` read
from FITS source `src`.  Argument `src` can be the name of a FITS file or a
FITS handle (an instance of `FitsIO`).  The keyword `ext` is to specify the
name or the number of the HDU to consider (the first one by default).  It must
correspond to a FITS *Image* extension.

Examples:

```julia
using EasyFITS
A = read(FitsImage, "image.fits")  # load the first HDU
A[2,3]                             # get value of data at indices (2,3)
A["BITPIX"]                        # get FITS bits per pixel
A.BITPIX                           # idem
get(FitsComment, A, "BITPIX")      # get the associated comment
A["STUFF"] = 1                     # set value of FITS keyword STUFF
A["STUFF"] = (1, "Some value")     # idem with value-comment pair
A.STUFF = 3                        # set value
A.STUFF = (3, "Some other value")  # idem
arr = get(Array, A)                # get the data part (a regular Julia array)
hdr = get(FitsHeader, A)           # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
```

It is also possible to specify other `Fits*` types as the first argument of `read`
to constrain the type of the result.  For instance:

```julia
using EasyFITS
read(FitsImage, "data.fits")          # load the first array and header
read(FitsHeader, "data.fits")         # reads only the header part
read(FitsImage{T}, "data.fits")       # yields pseudo-array with elements of type T
read(FitsImage{T,N}, "data.fits")     # yields N-dimensional pseudo-array with elements of type T
read(FitsArray, "data.fits")          # only load the array part (as a regular array)
read(FitsArray{T}, "data.fits")       # yields regular array with elements of type T
read(FitsArray{T,N}, "data.fits")     # yields N-dimensional regular array with elements of type T
```

Array slicing is also possible by specifying additional arguments:

```julia
using EasyFITS
read(FitsImage, "data.fits", :, 3:4)  # load slice A[:,3:4]
```

where it has been assumed that `A` is the full contents of the first HDU of
`"data.fits"`.

Note that the result of `read(FitsHeader,"data.fits")` can be indexed by
strings to access FITS keywords and implements the `obj.key` syntax.


## Writing FITS files

With **EasyFITS**, you can write rather complex FITS files in very lines of
code.

### Concise writing of FITS files

The `do ... end` construction is supported for instances of `FitsIO`.  For
instance something like:

```julia
FitsIO(path, "w!") do io
    write(io, dat1, hdr1)
    write(io, hdr2, dat2)
    write(io, img3)
end
```

can be used to create a new FITS file `path` with 3 extensions: `dat1` and
`hdr1` are the data and header parts of the 1st extension, `dat2` and `hdr2`
are those of the 2nd extension and `img3` is stored in the 3rd extension.  If
`img3` is an instance of `FitsImage`, both its data and header part will be
used.  Note that the order of the header and data parts is irrelevant, they are
recognized by their types.  The header part is optional, in the above example,
if `img3` is another kind of Julia array than `FitsImage`, a minimal FITS
header will be written.  At the end of the `do ... end` block, the FITS file is
closed so you do not have to worry about this.

The above example can be shortened as:

```julia
FitsIO(path, "w!") do io
    write(io, (dat1, hdr1), (hdr2, dat2), img3)
end
```

where each argument after `io` specifies the contents of a given extension.  If
data and header parts are two distinct instances, they have to be made into a
tuple so that there are no ambiguities.

Finally an even shorter equivalent of this example is:

```julia
write(FitFile, path, (dat1, hdr1), (hdr2, dat2), img3; overwrite=true)
```

Note the use of the `overwrite=true` keyword which indicates that if `path`
already exists, it can be overwritten without complaining.  This is the
counterpart of the `"w!"` mode in the call to `FitsIO`.  If you do not want to
overwrite an existing file, call `FitsIO` with `"w"` or `write(FitFile,...)`
without `overwrite=true` or with `overwrite=false` (which is the default).

As you can guess from these examples, the simplest way to save a Julia
array into a FITS file (as an *image* extension) is:

```julia
write(FitsFile, path, arr)
```

and to read it back, one of:

```julia
img = read(FitsImage, path)
arr = read(FitsArray, path)
```

to retrieve the contents of the file (both header and data parts in `img`, only
the data part in `arr` which is a regular Julia array).


### Building header information on the fly

In all previous example, the header parts can be built on the fly in three
different ways: using the `FitsHeader(...)` constructor, specifying a tuple of
`"key" => value` pairs, or specifying a named tuple `(key = value, ...)`.  For
instance, assuming `path` is a file name and `arr` a Julia array, the following
statements are equivalent:

```julia
write(FitsFile, path, arr,
      FitsHeader("HDUNAME" => ("Custom Extension", "Some comment"),
                 "VERSION" => 1.2))
write(FitsFile, path, arr,
      FitsHeader(HDUNAME = ("Custom Extension", "Some comment"),
                 VERSION = 1.2))
write(FitsFile, path, arr, ("HDUNAME" => ("Custom Extension", "Some comment"),
                            "VERSION" => 1.2))
write(FitsFile, path, arr, (HDUNAME = ("Custom Extension", "Some comment"),
                            VERSION = 1.2))
```

Note that we followed the convention that FITS keywords are in capital letters
and that an optional comment can be given for each FITS keyword by specifying
its value as a 2-tuple `(value,comment)`.  Above, the header specifications
have been wrapped into a single object: a `FitsHeader` or a tuple.  Using a
named tuple or a tuple of `"key" => value` pairs can also be used for writing
multiple extensions.  When a single extension is written, the tuple can be
avoided:

```julia
write(FitsFile, path, arr, "HDUNAME" => ("Custom Extension", "Some comment"),
      "VERSION" => 1.2)
write(FitsFile, path, arr, HDUNAME = ("Custom Extension", "Some comment"),
      VERSION = 1.2)
```

To avoid ambiguities, these two styles cannot be mixed.


### Extending for other kinds of data types

The API of **EasyFITS** is designed to be easily extensible for other data
types.   For instance, to benefit from the available methods, it may be possible
to just extend:

```julia
Base.write(io::FitsIO, dat::CustomDataType) = ...
```

to save `dat` appropriately in a FITS file and so that the end user can just call:

```julia
write(FitsFile, path, dat)
```

to save `dat` (and deal with compression and overwriting issues).

allowing for extraneous keywords to save in the header requires to extend
more variants of `write(io,dat,...)`.


## Automatic (de)compression

With **EasyFITS**, file names ending with `.gz` are automatically recognized
for compressed files.  When reading a FITS file, say `path`, if no file named
`path` exists but the is a file named `"$(path).gz"` this latter file will be
automatically open.  When creating a FITS file, the file is automatically
compressed if its name ends with `.gz`.


## Lower level methods

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
[fitsio-url]: https://github.com/JuliaAstro/CFITSIO.jl
[julia-url]: http://julialang.org/
[libcfitsio-url]: http://heasarc.gsfc.nasa.gov/fitsio/
