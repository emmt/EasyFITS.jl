# The *Grammar* of `EasyFITS`

FITS file may have quite complex structure and `EasyFITS` implement a *grammar*
that is intended to make clear the intention when reading the code and to help
to guess the correct syntax when writing code. Most of the assumed rules should
be familiar to Julia users.

The general syntax to read some data in a FITS file follows the following
pattern:

```julia
data = read([type,] src[, what...])
```

where:

- `type` is an optional Julia type to specify the expected type of the result.

- `src` is the source to read: a filename name, a `FitsFile` object
  representing an open FITS file, or a `FitsHDU` object representing a FITS
  Header Data Unit in an open FITS file.

- `what...` denotes optional arguments to specify which part(s) to read: index
  ranges of an array if reading a FITS image extension, column(s) and rows if
  reading a FITS table extension.

Overwriting the contents of an existing destination object `dest` may be done
by:

```julia
read!(dest, src, what...)
```

or by:

```julia
merge!(dest, src, what...)
```

when the intention is to preserve (part of) the contents of `dest`, or also by:

```julia
push!(dest, src, what...)
```

to append the read data to the contents of `dest`.


Create a new HDU in `file`:

```julia
hdu = write(file, hdutype, defs...)
```

Write some data in a FITS file:

```julia
write(file, header, data) -> file
```

which yields `file` so that chaining calls is possible to write several HDUs in
a row:

```julia
write(write(write(file, header1, data1), header2, data2), header3, data3)
```

which is equivalent to the more readable form:

```julia
write(file, header1, data1, header2, data2, header3, data3)
```
