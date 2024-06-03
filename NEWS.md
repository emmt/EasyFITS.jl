# User visible changes for EasyFITS

## Version 0.5.16

- `FITSHeaders` is an official package.

## Version 0.5.15

- Update compatibility for `TypeUtils`.

## Version 0.5.14

- Add property `hdu.data_axes` for FITS Image HDUs.

## Version 0.5.13

- Can read columns with given number of dimensions.

- Column data dimensions are considered in a more flexible way: the last array
  dimension is the row index, the other leading dimensions are the cell
  dimensions. For columns of strings, the first cell dimension is not smaller
  than the maximal number of bytes (ASCII characters) of any string in that
  column. The cell dimensions are considered to be the same if the leading
  dimensions are equal and the extra trailing dimensions, if any, are all ones.

- Package `FITSBase` has been renamed as `FITSHeaders`.

## Version 0.5.12

- Columns of strings in FITS tables:
  - An error is raised when attempting to write strings longer than the maximal
    length in a column of strings. The previous behavior was to truncate the
    strings and display a warning.
  - When cell dimensions are not specified in a column definition, it is
    assumed that a single value is stored by each cell of that column. This is
    inappropriate for strings for which the default amounts to storing a single
    character per cell. Cell dimensions are therefore now mandatory for
    defining columns of strings in a table HDU. Note that this is the simplest
    way to specify the maximum length of the strings.
  - Reading/writing a column of strings of length 1 has been fixed.

## Version 0.5.11

- Extend `haskey` for `FitsHDU` instances.

## Version 0.5.10

- Fix compatibility.

## Version 0.5.9

- New non-exported constant `EasyFITS.OptionalHeader` to match `nothing` or
  anything that can represent a FITS header.

- Other packages may simply extend
  `Base.write(file::FitsFile,header::EasyFITS.OptionalHeader,data::CustomType)`
  and/or `Base.read(::Type{<:CustomType},hdu::FitsHDU)` for their own type of
  data `CustomType` to specify how to save and/or load such kind of data.

## Version 0.5.8

- Import predicate functions `is_comment`, `is_end`, `is_naxis`, and
  `is_structural` from `BaseFITS`.

## Version 0.5.7

- Package `AsType` is now [`TypeUtils`](https://github.com/emmt/TypeUtils.jl).

## Version 0.5.6

This version introduces many changes mostly for reading/writing FITS table
extensions:

- The `rename` keyword can be used to specify a function to rename column names
  when reading a FITS table in a dictionary.

- Keywords `first` and `last` to specify the range of rows to read in a FITS
  table with the `read` method have been replaced by a `rows` argument which
  may be a single row index, a unit range of row numbers, or a colon `:` to
  read all rows (the default). The `read!` method keeps its `first` keyword to
  specify the first row to read.

- Reading a single row (e.g. by specifying `rows = 3` but not `rows = 3:3`),
  for a table yield a 0-dimensional result along the corresponding dimension.

- `read!(dict,hdu[,cols[,rows]])` replaces the contents of dictionary `dict`
  with columns `cols` of the FITS table `hdu` while
  `merge!(dict,hdu[,cols[,rows]])` merges columns `cols` of the FITS table
  `hdu` to the contents the dictionary `dict`.

- Reading column(s) from a FITS table can yield the column(s) values or the
  column(s) values *and* their units.

- FITS table extensions can be written by `write(file,hdr,cols)` with `hdr`
  specifying additional header records, columns `cols` specified by a
  collection of pairs like `key => vals` or `key => (vals, units)` with `key`
  the (symbolic) name of the column, `vals` its values, and `units` its
  optional units. The collection `cols` can be a dictionary, a named tuple, a
  vector of pairs, or a tuple of pairs.

- To avoid ambiguities, when writing complete FITS extensions in a single
  `write` call, two arguments must be supplied for each extension: one for the
  header (possibly `nothing`) and one for the data.

- **Columns of strings:**

  - Columns with string values and, possibly, multi-dimensional cells can be
    read/written as strings or as raw bytes. Although specified by FITS
    standard, multi-dimensional cells of strings may not be correctly supported
    by all software (notably not by
    [`fv`](https://heasarc.gsfc.nasa.gov/docs/software/ftools/fv/)). CFITSIO
    itself implements its own mechanisms for multi-dimensional cells of strings
    which are not part of FITS standard and does not understand the FITS
    standard rules. For this reason in `EasyFITS`, raw bytes are always used as
    an intermediate by to represent strings so as to shortcut the handling of
    strings by CFITSIO.

  - As assumed in the CFITSIO library and by the FITS standard, trailing spaces
    are not significant (and discarded) unless the string only consists in
    spaces if which case the first space is considered as significant (and
    kept). This is intended to distinguish null (empty) and non-null strings.

There are also some changes not related to FITS table extensions:

- The constant `EasyFITS.CFITSIO_VERSION` gives the version of the CFITSIO
  library for which the package has been built. At load time, it is checked
  that this version matches that of the dynamic library.

## Version 0.5.5

- `read(FitsHeader,filename;ext=...)` can be used to read the header of FITS
  extension `ext` in file `filename`.

- Empty FITS image extensions can be written.

## Version 0.5.4

- Use `AsType` package.

- Save arrays of Booleans as bytes in FITS image extensions.

## Version 0.5.3

- Fix directly writing a FITS image extension in a file with a given header and
  image array.

## Version 0.5.2

- Call `write(file::FitsFile,FitsImageHDU{T},dims...)` to create a new FITS image HDU
  with elements of type `T` and dimensions `dims...`.

- Fix writing image extensions with values of type `Bool`. In CFITSIO library,
  logical values must be considered as bytes in FITS image extensions.

## Version 0.5.1

- The function `hduname` imported by some other packages is back.

## Version 0.5.0

This version introduces lots of changes in the API. The syntax should be more
intuitive and consistent. To reduce the number of exported functions and avoid
using non-exported ones, *properties* are used extensively.

### General changes

- Type-stability has been improved in many places.

- The number of exported functions has been reduced to `write!`, `readfits`,
  `openfits`, `writefits`, and `writefits!`. All exported types are prefixed by
  `Fits` and all exported constants are prefixed by `FITS_`.

- String prefix `fits` can be used to mark names of FITS files when calling
  standard functions like `open`, `read`, and `write`. Decoration
  `FitsFile(filename)` can also be used if the file name is bound to a variable
  (also works for literal strings).

- `EasyFITS` no longer wraps over `FITSIO` but directly use `CFITSIO_jll`.

- `exists(path)` is no longer provided, call `ispath(path)` instead.


### FITS headers

- Access to header cards has been completely rewritten for speed and type
  stability.

- Headers are now specified as vectors of `key => val` pairs with `key` the
  keyword string and `val` the associated value and/or comment.

- Properties, that is the `card.field` syntax, are used to retrieve the
  different parts of a FITS header card (name, value, and comment), to parse
  the card value, and to get the card units if any.


### FITS table extensions

- When creating table extensions, columns are defined by vectors of `col =>
  def` pairs with column name `col` and column definition `def` that must
  include the type of the column values and may also include their units and
  the dimensions of the column cells.

- Writing columns is done by `write(hdu,col=>vals,...)` where `hdu` is the
  header data unit object of the FITS table extension, `col` is the column
  name/number, and `vals` are the values to write. The ellipsis `...` may be
  any other column data to write. Automatic type conversion is performed.
  Partial writing is possible by using the keyword `firstrow` to specify the
  first row to write.

- Reading column values is done by `read(hdu,col)` where `hdu` is the header
  data unit object of the FITS table extension and `col` is the column
  name/number. Partial reading is possible by using the keywords `first` and
  `last` to specify the first and last rows to read.


### FITS image extensions


## Versions 0.2.4, 0.2.5

Fix type piracy for `get` method when `FITSIO` version is greater or equal
0.17.

## Version 0.2.0

This version requires [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) and
[CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl).

## Version 0.1.0

This release is intended to work with
[FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) prior to the introduction
of [CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl) for the low-level
interface to the C library.
