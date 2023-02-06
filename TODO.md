# Things to do in EasyFITS

- Deal with long strings (`CONTINUE` FITS keyword).

- Deal with columns of bits.

- Deal with variable length columns.

- Use `@inbounds` to optimize some more loops.

- Abstract type `FitsHDU` and `FitsHeader(hdu::FitsHDU)` should be defined in
  `BaseFITS`.

- In `utils.jl` use a more elegant and secure way to deal with FITS booleans
  which are implemented as `Cchar` in CFITSIO and thus arrays of booleans are
  unfortunately thought as `CString` by Julia code wrapper.

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
