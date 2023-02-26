# User visible changes for EasyFITS

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
