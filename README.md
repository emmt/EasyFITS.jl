# Easy reading/writing of FITS files in Julia

[![Doc][doc-dev-img]][doc-dev-url]
[![License][license-img]][license-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Coverage][codecov-img]][codecov-url]

`EasyFITS` is a [Julia][julia-url] package designed to make it easier to read and write
data in [FITS](https://fits.gsfc.nasa.gov/fits_standard.html) format without sacrificing
performances, flexibility, or readability.


## A few examples

The full documentation is available [on-line][doc-dev-url].

The *Flexible Image Transport System* (or
[FITS](https://fits.gsfc.nasa.gov/fits_standard.html) for short) is a file
format widely used in Astronomy to store many kinds of data (images, tables,
etc.) and metadata. FITS files consist in a concatenation of Header Data Units
(HDUs) which each have a header part followed by a data part.

The following example demonstrates how to write a FITS file with 3 HDUs, an
*Image Extension* and two *Table Extensions*:

```julia
using Dates, EasyFITS
filename = "/tmp/test.fits";
arr = rand(Float32, (3,4,5));
nrows = 20;
inds = 1:nrows;
speed = rand(Float64, nrows);
mass = rand(Float32, nrows);
position = rand(Float32, 3, nrows);
phase = (1:7) .// 3;
amplitude = exp.(-1:-1:-7);
x = amplitude.*cos.(phase);
y = amplitude.*sin.(phase);
writefits(filename,
          #-----------------------------------------------------------------
          # First HDU must be a FITS "image", but data may be empty.
          #
          # Header part as a vector of `key=>val` or `key=>(val,com)` pairs:
          ["DATE"    => (now(), "date of creation"),
           "HISTORY" => "This file has been produced by EasyFITS",
           "USER"    => ENV["USER"]],
          # Data part as an array:
          arr,
          #-----------------------------------------------------------------
          # Second HDU, here a FITS "table".
          #
          # Header part of 2nd HDU as a tuple of pairs:
          ("EXTNAME" => ("MY-EXTENSION", "Name of this extension"),
           "EXTVER"  => (1, "Version of this extension")),
          # Data part is a table in the form of a named tuple:
          (Speed    = (speed, "km/s"),  # this column has units
           Indices  = inds,             # not this one
           Mass     = (mass, "kg"),
           Position = (position, "cm")),
          #-----------------------------------------------------------------
          # Third HDU, another FITS "table".
          #
          # Header part of 3rd HDU as a named tuple (note that keywords must
          # be in uppercase letters):
          (EXTNAME = ("MY-OTHER-EXTENSION", "Name of this other extension"),
           EXTVER  = (1, "Version of this other extension"),
           COMMENT = "This is an interesting comment"),
          # Data part is a table in the form of a vector of pairs (column names
          # can be strings or symbols but not a mixture):
          [:phase => ((180/π).*phase, "deg"),
           :amplitude => (amplitude, "V"),
           :xy => (hcat(x,y)', "V")])
```

Each HDU has a header part (the metadata) and a data part which is reflected by the pairs
of arguments after`filename`, the name of the file, in the above call to `writefits`. The
headers are provided by collections (a vector for the 1st one, a tuple for the 2nd one) of
pairs or by a named tuples (3rd one) associating a keyword with a value and a comment
(both optional). In a FITS *Image Extension*, the data can be any real-valued Julia array.
In a FITS *Table Extension* , the data part is provided by a collection of column names
associated with columns values and optional units. The columns in a FITS table must have
the same trailing dimension (interpreted as the rows of the table) but may have different
leading dimensions corresponding to the sizes of the column cells. In the above example,
the `"Position"` column has 3 values per cell (presumably the 3D coordinates), while other
columns have a single value per cell. Note that the 1st HDU, so-called *Primary HDU*, of a
FITS file must be a *Fits Image* (possibly empty), not a *FITS Table*.

To read the headers of the 1st and 2nd HDU of the file:

```julia
hdr1 = read(FitsHeader, filename)
hdr2 = read(FitsHeader, filename, ext=2)
```

yield two instance of `FitsHeader`. Reading the data parts is very easy:

```julia
dat1 = readfits(filename)
dat2 = readfits(filename, ext=2)
```

will yield an array `dat1` equal to `arr` and a dictionary `dat2` indexed by the column
names (in uppercase letters by default). For example:

``` julia
dat2["SPEED"] == speed
```

should hold.


## Installation

The easiest way to install `EasyFITS` is via Julia registry
[`EmmtRegistry`](https://github.com/emmt/EmmtRegistry):

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
pkg"add EasyFITS"
```

Adding the `General` registry (2nd line of the above example) is mandatory to have access
to the official Julia packages if you never have used the package manager before.


## Related projects

The [FITSIO](https://github.com/JuliaAstro/FITSIO.jl) package is another alternative to
read/write FITS files. `EasyFITS` is no longer based on `FITSIO` and now exploits
[Clang.jl][clang-url] to directly call the functions of the [CFITSIO][cfitsio-url] library
and [`FITSHeaders`](https://github.com/emmt/FITSHeaders.jl) to parse metadata (FITS header
cards).


[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/EasyFITS.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/EasyFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[github-ci-img]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml?query=branch%3Amaster

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
