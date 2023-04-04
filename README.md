# Easy reading/writing of FITS files in Julia

[![Doc][doc-dev-img]][doc-dev-url]
[![License][license-img]][license-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Build Status][appveyor-img]][appveyor-url]
[![Coverage][codecov-img]][codecov-url]

`EasyFITS` is a [Julia][julia-url] package designed to make it easier to read
and write data in [FITS (for)](https://fits.gsfc.nasa.gov/fits_standard.html) format without
sacrificing performances, flexibility, or readability.


## A few examples

The full documentation is available [on-line][doc-dev-url].

 The *Flexible Image Transport System* (or
 [FITS](https://fits.gsfc.nasa.gov/fits_standard.html) for short) is a file
 format widely used in Astronomy to store many kinds of data (images, tables,
 etc.) and metadata. FITS files consist in a concatenation of Header Data Units
 (HDUs) which each have a header part followed by a data part.

The following example demonstrates how to write a FITS file with 2 HDUs, an
*Image Extension* and a *Table Extension*:

```julia
using Dates, EasyFITS
arr = rand(Float32, (3,4,5));
nrows = 20;
inds = 1:nrows;
speed = rand(Float64, nrows);
mass = rand(Float32, nrows);
position = rand(Float32, 3, nrows);
writefits(filename,
          # Header part of 1st HDU.
          ["DATE"    => (now(), "date of creation"),
           "HISTORY" => "This file has been produced by EasyFITS",
           "USER"    => ENV["USER"]],
          # Data part (here, a FITS "image") of 1st HDU .
          arr,
          # Header part of 2nd HDU.
          ("EXTNAME" => ("MY-EXTENSION", "Name of this extension"),
           "EXTVER"  => (1, "Version of this extension")),
          # Data part (here, a FITS "table") of 2nd HDU .
          (Speed    = (speed, "km/s"),  # this column has units
           Indices  = inds,             # not this one
           Mass     = (mass, "kg"),
           Position = (position, "cm")))
```

Each HDU has a the header part (the metadata) and a data part which is
reflected by the pairs of arguments after the name of the file `filename` in
the above call to `writefits`. The two headers are provided by collections (a
vector for the 1st one, a tuple for the 2nd) of pairs associating a keyword
with a value and a comment (both optional). The data in a FITS *Image
Extension* is any real-valued Julia array. The data part in a FITS *Table
Extension* is provided by a collection of column names associated with columns
values and optional units. The columns in a FITS table must have the same
trailing dimension (interpreted as the rows of the table) but may have
different leading dimensions corresponding to the sizes of the column cells. In
the above example, the `"Position"` column has 3 values per cell (presumably
the 3D coordinates), while other columns have a single value per cell.

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

will yield an array `dat1` equal to `arr` and a dictionary `dat2` indexed by the
column names (in uppercase letters by default).  For example:

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

Adding the `General` registry (2nd line of the above example) is mandatory to
have access to the official Julia packages if you never have used the package
manager before.


## Related projects

The [FITSIO](https://github.com/JuliaAstro/FITSIO.jl) package is another
alternative to read/write FITS files. `EasyFITS` is no longer based on `FITSIO`
and now exploits [Clang.jl][clang-url] to directly call the functions of the
[CFITSIO][cfitsio-url] library and
[`BaseFITS`](https://github.com/emmt/BaseFITS.jl) to parse metadata (FITS
header cards).


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
