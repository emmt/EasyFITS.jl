# Easy reading/writing of FITS files in Julia

[![License][license-img]][license-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Build Status][appveyor-img]][appveyor-url]
[![Coverage][codecov-img]][codecov-url]

This package is to facilitate reading/writing
[FITS](https://fits.gsfc.nasa.gov/fits_standard.html) files (widely used in
Astronomy) in [Julia][julia-url]. `EasyFITS` uses [Clang.jl][clang-url] to call
the functions of the [CFITSIO][cfitsio-url] library.


## Usage

To use the `EasyFITS` package:

```julia
using EasyFITS
```

This imports a bunch of types (all prefixed by `Fits...`), some constants (all
prefixed by `FITS_...`), a few methods (`write!`, `readfits`, `openfits`,
`writefits`, and `writefits!`), and the macro `@Fits_str`.


## FITS files

Open FITS files are represented by objects of type `FitsFile`. To open an
existing FITS file there are several possibilities:

```julia
file = FitsFile(filename)
file = openfits(filename)
file = open(FitsFile, filename)
```

where `filename` is the name of the FITS file.

The methods to open a FITS file all take an optional argument `mode` after
`filename` which can be:

- `"r"` (the default mode) to open an existing file for reading only.

- `"r+"` to read an existing file for reading and appending to its contents.

- `"w"` to create a new FITS file that must not already exist. An error is
  thrown if the file already exists.

- `"w!"` to open the file for writing. If the file already exists, it is
  (silently) overwritten.

The keyword `extended` can be used to specify whether to use the [extended file
name
syntax](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html)
implemented by the `FITSIO` library.

It is not mandatory to call `close(file)` to close the FITS file, this is
automatically done when the `file` object is garbage collected. Calling
`close(file)` is however needed in order to ensure that a FITS file open for
writing is up to date before `file` be garbage collected. To automatically
close a FITS file, use the do-block syntax:

``` julia
open(FitsFile, "test.fits.gz", "w!") do file
    write(file, ...)
    ...
end
```

FITS files consist in a concatenation of header data units (HDUs) which each
have a header part followed by a data part. Open FITS files can be indexed to
walk along their HDUs:

``` julia
first(file)             # yields the first HDU
file[firstindex(file)]  # idem
file[1]                 # idem
file[2]                 # yields the second HDU
last(file)              # yields the last HDU
file[lastindex(file)]   # idem
file[end]               # idem
file[length(file)]      # idem
file["name"]            # next HDU matching given name or nothing
```

In other words, FITS files behave as vectors of HDUs with 1-based integer
indices but can also be indexed by strings.

Note that indexation by name yields the next matching HDU. To rewind the search
of HDUs by their names, just call `seekstart(file)`. This makes easy to travel
through all HDUs of a given name:

``` julia
seekstart(file)
while true
    hdu = file["BINTABLE"]
    hdu === nothing && break
    # do something with this HDU
    ...
end
```

Call `seek(file,n)` to move to the `n`-th HDU without retrieving the HDU object.
`seekstart(file)` and `seekend(file)` are similar but move to the first/last HDU.
Call `position(file)` to figure out the number of the current HDU.


## `EasyFITS` objects

### FITS file objects

An instance of `FitsFile` is a wrapper around a FITS file open for reading or for
writing:

```julia
file = FitsFile(filename, mode)
```

where `filename` is the name of the FITS file and `mode` can be:

- `"r"` or `"r+"` to read or append to the contents of and existing FITS file.
  If file `filename` does not exist but `filename` does not end with the `".gz"`
  extension and `"\$filename.gz"` does exist, then the compressed file `"$filename.gz"`
  is open instead.

- `"w"` to create a new FITS file named `filename` that must not already exist. An
  error is thrown if the file already exists.

- `"w!"` to open FITS file named `filename` for writing. If the file already
  exists, it is (silently) overwritten.

An alternative is to call `open` as:

```julia
file = open(FitsFile, filename, mode)
file = openfits(filename, mode)
```


Call `close(file)` to close the FITS file associated with the `FitsFile` instance
`file`. Call `isopen(file)` to check whether the FITS file associated with the
`FitsFile` instance `file` is open. Closing the FITS file is automatically done, if
needed, when the instance is garbage collected.

The do-block syntax is supported to automatically close the FITS file:

```julia
FitsFile(filename, mode="r") do file
    # use FITS handle file
    ...
end
```

An instance of `FitsFile` is a collection of *Header Data Units* (HDU) and
implements indexation and iteration. Assuming `file` is a `FitsFile` object, then:

- `file[i]` yields the `i`-th HDU.

- `length(file)` yields the number of HDUs.

- `file[name]` or `file[name,vers]` yields the HDU whose `EXTNAME` (or `HDUNAME`)
  keyword is equal to `name` (a string) and, optionally, whose `EXTVER` (or
  `HDUVER`) keyword is equal to `vers` (an integer).

- You can do `for hdu in file; ...; end` to iterate through all HDU's of `file`.

- Methods `findfirst(p,file)`, `findlast(p,file)`, `findnext(p,file,i)` and
  `findprev(p,file,i)` can be used on `FitsFile` object `file` to search for a
  specific HDU. These methods test each HDU (starting at initial index `i` for
  `findnext` and `findprev`) with the predicate function `p` (called with a
  `FitsHDU` argument) and return the index (an integer) of the first HDU for
  which the predicate yields `true` or `nothing` if this never occurs. These
  methods are build upon `EasyFITS.find` which may be directly called.


## Reading data from FITS files

To load a FITS Image extension as an instance of `FitsImage`, call:

```julia
read(FitsImage, src, args...; ext=1) -> A
```

which yields a pseudo-array `A` with the contents of the extension `ext` read
from FITS source `src`.  Argument `src` can be the name of a FITS file or a
FITS handle (an instance of `FitsFile`).  The keyword `ext` is to specify the
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
arr = convert(Array, A)            # get the data part (a regular Julia array)
hdr = get(FitsHeader, A)           # get the header part
EasyFITS.nkeys(A)                  # get the number of keywords
EasyFITS.nkeys(hdr)                # get the number of keywords
keys(A)                            # get the list of keywords
keys(hdr)                          # get the list of keywords
delete!(hdr, key)                  # delete FITS keyword from header
delete!(A, key)                    # delete FITS keyword from annotated array
pop!(hdr, key[, def])              # pop FITS keyword out of header
pop!(A, key[, def])                # pop FITS keyword out of annotated array
```

It is also possible to specify other `FITS*` types as the first argument of
`read` to constrain the type of the result.  For instance:

```julia
using EasyFITS
read(FitsImage, "data.fits")          # load the first array and header
read(FitsHeader, "data.fits")         # reads only the header part
read(FitsImage{T}, "data.fits")       # yields pseudo-array with elements of type T
read(FitsImage{T,N}, "data.fits")     # yields N-dimensional pseudo-array with elements of type T
read(Array, fits"data.fits")          # only load the array part (as a regular array)
read(Array{T}, fits"data.fits")       # yields regular array with elements of type T
read(Array{T,N}, fits"data.fits")     # yields N-dimensional regular array with elements of type T
```

Above, the `fits` string prefix was needed to indicate that the file was a FITS
file and avoid type-piracy. The very same result is obtained by using the
`readfits` function instead:

```julia
readfits(Array, "data.fits")          # only load the array part (as a regular array)
readfits(Array{T}, "data.fits")       # yields regular array with elements of type T
readfits(Array{T,N}, "data.fits")     # yields N-dimensional regular array with elements of type T
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

With `EasyFITS`, you can write rather complex FITS files in very lines of
code.

### Writing a FITS Image

If `A` is an instance of `FitsImage` (see above), then saving its contents
(data and header parts) as a FITS file is as easy as:

```julia
write(filename, A)
```

with `filename` the name of the output file.  Keyword `overwrite` may be used to
specify whether overwriting an existing file is allowed.  Another possibility
to force overwriting is to call:

```julia
write!(filename, A)
```


### Concise writing of FITS files

The `do ... end` construction is supported for instances of `FitsFile`.  For
instance something like:

```julia
FitsFile(filename, "w!") do file
    write(file, dat1, hdr1)
    write(file, hdr2, dat2)
    write(file, img3)
end
```

can be used to create a new FITS file `filename` with 3 extensions: `dat1` and
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
FitsFile(filename, "w!") do file
    write(file, (dat1, hdr1), (hdr2, dat2), img3)
end
```

where each argument after `file` specifies the contents of a given extension.  If
data and header parts are two distinct instances, they have to be made into a
tuple so that there are no ambiguities.

Finally an even shorter equivalent of this example is:

```julia
write(FitFile, filename, (dat1, hdr1), (hdr2, dat2), img3; overwrite=true)
```

Note the use of the `overwrite=true` keyword which indicates that if `filename`
already exists, it can be overwritten without complaining.  This is the
counterpart of the `"w!"` mode in the call to `FitsFile`.  If you do not want to
overwrite an existing file, call `FitsFile` with `"w"` or `write(FitFile,...)`
without `overwrite=true` or with `overwrite=false` (which is the default).

As you can guess from these examples, the simplest way to save a Julia
array into a FITS file (as an *image* extension) is:

```julia
write(FitsFile, filename, arr)
```

and to read it back, one of:

```julia
img = read(FitsImage, filename)
arr = read(FitsArray, filename)
```

to retrieve the contents of the file (both header and data parts in `img`, only
the data part in `arr` which is a regular Julia array).


### Building header information on the fly

In all previous example, the header parts can be built on the fly in three
different ways: using the `FitsHeader(...)` constructor, specifying a tuple of
`"key" => value` pairs, or specifying a named tuple `(key = value, ...)`.  For
instance, assuming `filename` is a file name and `arr` a Julia array, the following
statements are equivalent:

```julia
write(FitsFile, filename, arr,
      FitsHeader("HDUNAME" => ("Custom Extension", "Some comment"),
                 "VERSION" => 1.2))
write(FitsFile, filename, arr,
      FitsHeader(HDUNAME = ("Custom Extension", "Some comment"),
                 VERSION = 1.2))
write(FitsFile, filename, arr, ("HDUNAME" => ("Custom Extension", "Some comment"),
                            "VERSION" => 1.2))
write(FitsFile, filename, arr, (HDUNAME = ("Custom Extension", "Some comment"),
                            VERSION = 1.2))
```

Note that we followed the convention that FITS keywords are in capital letters
and that an optional comment can be given for each FITS keyword by specifying
its value as a 2-tuple `(value,comment)`. Above, the header specifications have
been wrapped into a single object: a `FitsHeader` or a tuple. Using a named
tuple or a tuple of `"key" => value` pairs can also be used for writing
multiple extensions. When a single extension is written, the tuple can be
avoided:

```julia
write(FitsFile, filename, arr, "HDUNAME" => ("Custom Extension", "Some comment"),
      "VERSION" => 1.2)
write(FitsFile, filename, arr, HDUNAME = ("Custom Extension", "Some comment"),
      VERSION = 1.2)
```

To avoid ambiguities, these two styles cannot be mixed.


### Extending for other kinds of data types

The API of `EasyFITS` is designed to be easily extensible for other data
types.   For instance, to benefit from the available methods, it may be possible
to just extend:

```julia
Base.write(file::FitsFile, dat::CustomDataType) = ...
```

to save `dat` appropriately in a FITS file and so that the end user can just call:

```julia
write(FitsFile, filename, dat)
```

to save `dat` (and deal with compression and overwriting issues).

allowing for extraneous keywords to save in the header requires to extend
more variants of `write(file,dat,...)`.


## Automatic (de)compression

With `EasyFITS`, file names ending with `.gz` are automatically recognized for
compressed files. When reading a FITS file, say `filename`, if no file named `filename`
exists but there is a file named `"$(filename).gz"` this latter file will be
automatically open. When creating a FITS file, the file is automatically
compressed if its name ends with `.gz`.


## Lower level methods, types

To check whether a file `filename` already exists in the file system, call:

```julia
isfile(filename) -> bool
```

### FITS bits per pixel (BITPIX)

Call

```julia
FitsBitpix(arg)
```

to get a singleton type which encapsulates the FITS *bitpix* (for
*bits-per-pixel*) code identifying the array element type relevant for `arg`.
Argument `arg` can be a Julia type, an integer (interpreted as a FITS bitpix
code), an array, a FITS HDU/image/header.

Conversely call `eltype(bpx)` to convert FITS bitpix `bpx` into a Julia type.


## Naming conventions

To avoid conflicts such as *type piracy*, all exported methods but `exists` and
`write!` have their names prefixed by `FITS*`.

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/EasyFITS.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/EasyFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[github-ci-img]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml?query=branch%3Amaster

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/EasyFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/EasyFITS-jl/branch/master

[codecov-img]: http://codecov.io/github/emmt/EasyFITS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/EasyFITS.jl?branch=master

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/

[fitsbase-url]: https://github.com/emmt/FITSIO.jl
[fitsio-url]: https://github.com/JuliaAstro/FITSIO.jl
[cfitsio-url]: https://github.com/JuliaAstro/CFITSIO.jl
[julia-url]: http://julialang.org/
[libcfitsio-url]: http://heasarc.gsfc.nasa.gov/fitsio/
[clang-url]:https://github.com/JuliaInterop/Clang.jl
