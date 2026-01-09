# User visible changes for `EasyFITS`

This page describes the most important changes in `EasyFITS`. The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic
Versioning](https://semver.org).

## Unreleased

### Fixed

- `writefits!` always overwrites even though `overwrite=false` is specified.

- Doc. has been updated and fixed.

## Version 0.7.2 [2026-01-08]

### Fixed

- When reading a column of a FITS Table, it is always a `FitsError` that is thrown if the
  error occurs in a call to a C function of the `CFITSIO` library. Formerly, `EasyFITS`
  could throw a detailed `ErrorException` if the error was due to a non-existing column
  (with status `EasyFITS.CFITSIO.COL_NOT_FOUND`). This should fix a bug in the `OIFITS`
  package.

## Version 0.7.1 [2025-10-16]

### Fixed

- `EasyFITS` no longer crashes on 32-bit machines. There were two potential issues: (i)
  having the status returned as `struct Status; code::Cint; end` results in stack overflows
  (it is now again a simple `Cint`) and (ii) `Ptr{Ptr{fitsfile}}` had to be replaced by
  `Ptr{Cvoid}` in `@ccall` to avoid segmentation faults.

### Changed

- It should no longer be necessary to re-build `EasyFITS` for each new versions of
  `CFITSIO_jll`. The generated wrapper code is shipped with `EasyFITS` and
  [Clang.jl](https://github.com/JuliaInterop/Clang.jl) is only needed to regenerate this
  code when the API of the CFITSIO library has breaking changes.

- The absolute path of the FITS file is saved in the `FitsFile` structure and returned by
  the `pathof` function.

### Added

- Non-exported public type `EasyFITS.OutputCstring` for output C string arguments in
  `ccall`.

- Non-exported public function `EasyFITS.cfitsio_errmsg` to retrieve CFITSION error
  messages.

## Version 0.7.0 [2025-07-21]

### Changed

- [Clang.jl](https://github.com/JuliaInterop/Clang.jl) is no longer used to generate the
  wrapper code. The new generator is a simple parser of CFITSIO header files, it is faster
  and has no dependencies other than `CFITSIO_jll`.

- Pass all [`Aqua.jl`](https://github.com/JuliaTesting/Aqua.jl) tests.


### Added

- `hdu = file[pred]` can be used to retrieve the first HDU of the FITS `file` for which
  predicate function `pred` is true.


### Breaking changes

- `read!(dict::AbstractDict,hdu::FitsTableHDU,...)` merges some columns of the FITS Table
  in `hdu` with the contents of `dict` but no longer deletes existing contents of `dict`.
  Call `read!(empty!(dict),hdu,...)` for that.

- `merge!(dict::AbstractDict,hdu::FitsTableHDU,...)` to merge some columns of the FITS
  Table in `hdu` with the contents of `dict` is no longer supported as it is inconsistent
  with the usual meaning of `merge!`. Call `read!(dict,hdu,...)` instead.

- `push!(vec::AbstractVector,hdu::FitsTableHDU,...)` to append some columns of the FITS
  Table in `hdu` to the vector `vec` is no longer supported as it is inconsistent with the
  usual meaning of `push!`. Call `append!(vec, read(Vector,hdu,...))` instead.

- To avoid type-piracy, `read(FitsHeader, filename)` is replaced by `readfits(FitsHeader,
  filename)`.


### Deprecations

Deprecate some methods that are either redundant with shorter calls or which do not
implement the usual behavior of a base method:

- Deprecate `hdu = write(T::Type{<:FitsHDU},file::FitsFile,...)` in favor of `hdu =
  T(file,...)` because (i) it is shorter to directly call the constructor, (ii) creating
  the HDU does not actually write something to the file (this is deferred until data is
  written or another HDU is created), and (iii) is is unusual in Julia to store the value
  returned by a `write` call in a variable for further use.

- Deprecate `file = open(FitsFile,filename,...)` in favor of `file =
  openfits(filename,...)` or `file = FitsFile(filename,...)`.

- Deprecate `read([R::Type,]FitsFile,filename,...)` in favor of
  `readfits([R::Type,],filename,...)`.

- Deprecate `read!(dest,FitsFile,filename,...)` in favor of
  `readfits!(dest,filename,...)`.

- Deprecate `write(FitsFile,filename,...)` and `write!(FitsFile,filename,...)` in favor of
  `writefits(filename,...)` and `writefits!(filename,...)`.


## Version 0.6.1

- Avoid an ambiguity: union `EasyFITS.ColumnIdent` is for specifying a single column while
  union `EasyFITS.Columns` is for specifying several columns.

- Fix `Base.show` and `Base.length` for closed `FitsFile` (solves issue #10).

- Extend `Base.haskey` for `FitsFile`.

- Restrict version of `CFITSIO_jll` to solve CFITSIO bug related to non-US locales (solves
  issue #7).

## Version 0.6.0

- HDUs are created using their constructor, `FitsImageHDU` or `FitsTableHDU`. Users should
  replace calls like `write(file,FitsImageHDU{T},dims)` by `FitsImageHDU{T}(file,dims)`
  and calls like `write(file,FitsTableHDU,cols)` by `FitsTableHDU{T}(file,cols)`.

- Calling an HDU constructor to retrieve the HDU of a specific extension in a FITS file is
  no longer supported. Users should use the `file[ext]` syntax or call the `getindex` base
  method, possibly, with a type assertion. For example write
  `file[ext]::FitsImageHDU{Float32,3}` instead of `FitsImageHDU{Float32,3}(file,ext)`.

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

- Columns of strings in FITS Tables:
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

This version introduces many changes mostly for reading/writing FITS Table
extensions:

- The `rename` keyword can be used to specify a function to rename column names
  when reading a FITS Table in a dictionary.

- Keywords `first` and `last` to specify the range of rows to read in a FITS
  Table with the `read` method have been replaced by a `rows` argument which
  may be a single row index, a unit range of row numbers, or a colon `:` to
  read all rows (the default). The `read!` method keeps its `first` keyword to
  specify the first row to read.

- Reading a single row (e.g. by specifying `rows = 3` but not `rows = 3:3`),
  for a table yield a 0-dimensional result along the corresponding dimension.

- `read!(dict,hdu[,cols[,rows]])` replaces the contents of dictionary `dict`
  with columns `cols` of the FITS Table `hdu` while
  `merge!(dict,hdu[,cols[,rows]])` merges columns `cols` of the FITS Table
  `hdu` to the contents the dictionary `dict`.

- Reading column(s) from a FITS Table can yield the column(s) values or the
  column(s) values *and* their units.

- FITS Table extensions can be written by `write(file,hdr,cols)` with `hdr`
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

There are also some changes not related to FITS Table extensions:

- The constant `EasyFITS.CFITSIO_VERSION` gives the version of the CFITSIO
  library for which the package has been built. At load time, it is checked
  that this version matches that of the dynamic library.

## Version 0.5.5

- `read(FitsHeader,filename;ext=...)` can be used to read the header of FITS
  extension `ext` in file `filename`.

- Empty FITS Image extensions can be written.

## Version 0.5.4

- Use `AsType` package.

- Save arrays of Booleans as bytes in FITS Image extensions.

## Version 0.5.3

- Fix directly writing a FITS Image extension in a file with a given header and
  image array.

## Version 0.5.2

- Call `write(file::FitsFile,FitsImageHDU{T},dims...)` to create a new FITS Image HDU
  with elements of type `T` and dimensions `dims...`.

- Fix writing image extensions with values of type `Bool`. In CFITSIO library,
  logical values must be considered as bytes in FITS Image extensions.

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


### FITS Table extensions

- When creating table extensions, columns are defined by vectors of `col =>
  def` pairs with column name `col` and column definition `def` that must
  include the type of the column values and may also include their units and
  the dimensions of the column cells.

- Writing columns is done by `write(hdu,col=>vals,...)` where `hdu` is the
  header data unit object of the FITS Table extension, `col` is the column
  name/number, and `vals` are the values to write. The ellipsis `...` may be
  any other column data to write. Automatic type conversion is performed.
  Partial writing is possible by using the keyword `firstrow` to specify the
  first row to write.

- Reading column values is done by `read(hdu,col)` where `hdu` is the header
  data unit object of the FITS Table extension and `col` is the column
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
