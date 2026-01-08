# Reading FITS files

## Direct reading of FITS data

The simplest way to read some data in a FITS file is to call [`readfits`](@ref):

```julia
data = readfits(filename, args...; kwds...)
```

where `filename` is the name of the file while `args...` and `kwds...` are optional
arguments and keywords to specify which HDU and which part of the data to read.

By default, the returned data are read in the first FITS extension of `filename`. Keyword
`ext` may however be set with a number, a name, or a predicate function to select another
extension. Another possibility is to specify the keyword `extended = true` to open the file
using the [extended file name
syntax](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html) implemented
by the `CFITSIO` library.

If no optional arguments `args...` are specified, all the data part of the selected FITS HDU
is read and returned. Otherwise, arguments `args...` indicate which part the HDU data to
read:

- For a **FITS Image HDU**, `args...` specifies the ranges of pixels to read along the
  dimensions with a syntax similar to that of sub-arrays or views in Julia. For example:

  ```julia
  data = readfits(filename, :, :, 5)
  ```

  would read the 5th slice in a 3-dimensional image. This is equivalent but more efficient
  than:

  ```julia
  data = readfits(filename)[:, :, 5]
  ```

  which amounts to reading all the data and then only keep the 5th slice.

- For a **FITS Table HDU**, `args...` may be up to 2 arguments, say `cols` and `rows`, to
  respectively select a subset of columns and rows (if the latter is not specified, all rows
  of the table are read). For example:

  ```julia
  A = readfits(filename, ("Speed", "Height"))
  B = readfits(filename, :, 11:40)
  ```

  respectively yield the columns named `Speed` and `Height` of the table and all the columns
  of the table but only for rows in the range `11:40`. Keyword `case` may be used to
  indicate whether letter case does matter in the column names.

!!! note
    In `EasyFITS`, the *rows* of a table correspond to the last dimension of arrays. This is
    to have the same storage order in memory and in the FITS file. Method `permutedims` can
    be used is this convention does not suit you.

The type of the object returned by [`readfits`](@ref) depends on the kind of the FITS
extension and may also depend on the optional arguments `args...`:

- For a **FITS Image HDU**, the read data is returned as a Julia `Array`.

- For a **FITS Table HDU**, if `cols` is a single column name or number, an array of the
  columns values is returned, otherwise, a dictionary indexed by the column names is
  returned. Note that, a column range like `4:4` would yield a dictionary with a single
  column (the 4th one). Keywords `case` and `rename` are available to indicate how to search
  the columns by name in the table and how to translate these names into dictionary keys.

To avoid ambiguities or for improved type-stability, an optional leading type argument can
be specified in [`readfits`](@ref) to indicate the expected type of the returned data:

```julia
data = readfits(R::Type, filename, args...; kwds...)
```

which warrants that `data isa R` holds. For a FITS Image extension, array type parameters
such as the element type and the number of dimensions may be specified in `R`. For example:

```julia
data = readfits(Array{Float32,3}, filename)
```

ensures that `data` is a single precision floating-point 3-dimensional array; while:

```julia
data = readfits(Dict, filename, cols; ext=2)
```

ensures that the table in 2nd FITS Header Data Unit is returned as a dictionary even though
`cols` specifies a single column.


## Direct reading of a FITS header

For efficiency, FITS headers are returned as instances of `FitsHeader` provided by the
[`FITSHeaders`](https://github.com/emmt/FITSHeaders.jl) package. To avoid type piracy, there
is as yet no direct method of reading a FITS header. The most compact way to read the header
of a known FITS HDU is to do:

```julia
hdr = FitsHeader(FitsFile(filename)[ext])
```

where `filename` is the name of the file and `ext` is an HDU number, name, or a predicate
function to read the first HDU for which the predicate function yields true. The above
one-liner exploits the fact that FITS files are automatically closed when garbage collected.


## Advanced reading of FITS files

Direct reading of a FITS file with [`readfits`](@ref) is suitable to read the data of a
single FITS HDU at a known location (HDU number or name). More advanced operations such as
reading several pieces of data, searching the file, checking the contents, etc., can be
performed by the following steps:

1. Open the FITS file for reading with the [`FitsFile`](@ref) constructor or the
   [`openfits`](@ref) function.

2. Move to the FITS HDU of interest (see [Indexing and searching HDUs](@ref)).

3. Read some HDU content: [read the header](@ref FitsHeader(::FitsHDU)), [read image
   *pixels*](@ref read(::FitsImageHDU)), [read a table column](@ref read(::FitsTableHDU,
   ::String)), or [read several table columns](@ref read(::FitsTableHDU)).

4. Eventually close the FITS file with `close`. Closing the FITS file is automatically done
   when the object representing the file is no longer used and garbage collected, closing
   the file is therefore optional.

Remarks:

- Steps 2 and 3 can be repeated as many times as necessary to read more data and/or other
  headers.

- In-place reading is also possible with `read!(dest, hdu, ...)`.


For example, reading the data and the header of the HDU at extension `ext` in FITS file
named `filename` is done by:

```julia
file = FitsFile(filename)
hdu = file[ext]
hdr = FitsHeader(hdu)
data = read(hdu)
close(file)
```

or, using the `do`-block syntax:

```julia
hdr, data = FitsFile(filename) do file
    hdu = file[ext]
    FitsHeader(hdu), read(hdu)
end
```

In these examples, the [`FitsFile`](@ref) constructor can be replaced by the
[`openfits`](@ref) function to suit your taste.
