# Things to do in EasyFITS

* Provide `FitsTable` for FITS Table extensions.

* Use only low-level `CFITSIO` module.
* Cleanup `get(...)` extensions.
* Extend `Array(...)` to yield the array contents of a FITS image.
* Provide `comment!` to set the comment part.
* HDU can be specified by name ("EXTNAME" or "HDUNAME") when loading a FITS
  Image.
* When an "IMAGE" (resp. "TABLE") extension is expected without other
  constraint, find the first one.
* Filter keyword values (e.g., `Irrational`) to prevent problems.
* Deprecate `exists(path)` in favor of `isfile(path)`.
* `"r+"` mode.
