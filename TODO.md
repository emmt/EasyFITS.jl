* Deprecate `getheader` in favor of `header`.
* Deprecate `getdata` in favor of `data` or `array`.
* Deprecate `getcomment` in favor of `comment`.
* Provide `comment!` to set the comment part.
* HDU can be specified by name ("EXTNAME" or "HDUNAME") when loading a FITS
  Image.
* `readfits(...)` can also be `read(FitsImage, ...)` and should perhaps be
  renamed as `loadfitsimage` or `readfitsimage` (to prepare for
  `readfitstable`).
* Provide `FitsTable` for FITS Table extensions.
* When an "IMAGE" (resp. "TABLE") extension is expected without other
  constraint, find the first one.
* Filter keyword values (e.g., `Irrational`) to prevent problems.
* Deprecate `exists(path)` in favor of `isfile(path)`.
* `"r+"` mode.
* Implement `eltype(hdu)`.
