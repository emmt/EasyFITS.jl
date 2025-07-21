# Things to do in EasyFITS

- Document how to:
  - Copy a HDU from a file to another.
  - Change column name.
  - Change keyword name.
  - Edit keyword comment.

- Add properties `hdu.column_eltype`, `hdu.cell_size`, `hdu.cell_ndims`,
  `hdu.cell_eltype` or alike.

- Fix `@inferred read(Array{UInt8}, hdu, :col1)` for which inferred type is
  `Any` not `Array{UInt8}`.

- Writing of a table columns with a dictionary.

- Check performances of the different ways to dispatch on the pairs specifying
  columns when writing a FITS table.

- Reading a table as a vector of columns is not tested. Methods to push a table
  column in a vector are ambiguous and should probably be removed.

- Make `EasyFITS` understand the rules of CFITSIO for reading multi-dimensional
  cells of strings. Not for writing though.

- Deal with long strings (`CONTINUE` FITS keyword).

- Deal with columns of bits.

- Deal with variable length columns.

- Use `@inbounds` to optimize some more loops.

- Abstract type `FitsHDU` and `FitsHeader(hdu::FitsHDU)` should be defined in
  `BaseFITS`.

- In `utils.jl` use a more elegant and secure way to deal with FITS Booleans
  which are implemented as `Cchar` in CFITSIO and thus arrays of Booleans are
  unfortunately thought as `CString` by Julia code wrapper.

- Extension `ext` may be specified as a predicate function.

- FITS cards are (restricted) ASCII strings. The FITS standard states that FITS
  header cards exclusively consist in the characters whose hexadecimal values
  are in the range `0x20` through `0x7E`. The ASCII control characters with
  hexadecimal values less than `0x20` (including the null, tab, carriage
  return, and line-feed characters), and the delete character (hexadecimal
  value `0x7F`) must not appear anywhere within a keyword record.

  For keywords, all digits `0` through `9` (hexadecimal codes `0x30` to `0x39`)
  and upper case Latin alphabetic characters `A` through `Z` (hexadecimal codes
  `0x41` to `0x5A`) are permitted; lower-case characters shall not be used. The
  underscore (`_`, hexadecimal code `0x5F`) and hyphen (`-`, hexadecimal code
  `0x2D`) are also permitted. Space ` ` (hexadecimal code `0x20`) may appear in
  `HIERARCH` keywords. Equal sign `=` (hexadecimal code `0x3D`) followed by a
  space is used to indicate a keyword value. A slash `/` (hexadecimal code
  `0x2F`) is used to separate value and comment.


## Examples

Create an empty Image HDU with a header:

```julia
f = openfits(filename, "w!);
hdu = FitsImageHDU(f); # same as: hdu = FitsImageHDU{UInt8,0}(f,());
hdu[key_1] = (value_1, "comment_1");
hdu[key_2] = (value_2, "comment_2");
# etc.
hdu = ... # Create anothe HDU
Base.write(f::FitsFile, hdr::HeaderLike, nothing)
```
