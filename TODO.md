* Deprecate `getheader` in favor of `header`.
* Deprecate `getdata` in favor of `data` or `array`.
* Deprecate `getcomment` in favor of `comment`.
* Provide `comment!` to set the comment part.
* HDU can be specified by name ("EXTNAME" or "HDUNAME") when loading a FITS
  Image.
* `loadfits(...)` can be `read(EasyFITS.Image, ...)` and should perhaps be
  renamed as `loadfitsimage` or `readfitsimage` (to prepare for
  `readfitstable`).
* Provide means to write an `EasyFITS.Image`.
* Provide `FitsTable` for FITS Table extensions.
* When an "IMAGE" (resp. "TABLE") extension is expected without other
  constraint, find the first one.
* Make `FitsHeader` not an alias to `FITSIO.FITSHeader` but an object indexable
  as the `FitsImage`.
