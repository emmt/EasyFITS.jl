# Navigating in FITS files

This section explains how to open and navigate into a FITS file, that is move to a given
HDU.

## Open a FITS file

In `EasyFITS`, a FITS file is represented by an instance of [`FitsFile`](@ref). To create a
new FITS file or to open an existing FITS file, simply call the [`FitsFile`](@ref)
constructor or the [`openfits`](@ref) function:

```julia
file = FitsFile(filename, mode="r"; extended=false)
```

with `filename` the name of the FITS file and one of the following modes:

- `"r"` (the default) to open an existing FITS file for reading only.

- `"r+"` to open an existing FITS file for reading and appending to its contents.

- `"w"` to create a new FITS file, throwing an error if the file already exists.

- `"w!"` to create a new FITS file, silently overwriting the file if it already exists.

The `extended` keyword specifies whether to use the [extended file name
syntax](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html)
implemented by the `CFITSIO` library.

It is not mandatory to call `close(file)` to close the FITS file, this is automatically
done when the `file` object is garbage collected. Calling `close(file)` is however needed
to ensure that a FITS file open for writing is up to date before `file` be garbage
collected. To make sure the file contents is up to date without explicitly closing it,
call `flush(file)` instead. To automatically close a FITS file, use the `do`-block syntax:

``` julia
FitsFile("test.fits.gz", "w!") do file
    do_something_with(file)
end
```

For all these above examples, the [`openfits`](@ref) function can be used in-place of the
[`FitsFile`](@ref) constructor.

Standard Julia methods for file-like objects are available for an open FITS file:

```julia
isopen(file)     # yields whether file is open
close(file)      # close the file
pathof(file)     # yields the path to the file
filemode(file)   # yields the mode, one of `:r`, `:rw`, or `:w`
isreadable(file) # yields whether file is readable
isreadonly(file) # yields whether file is readable and not writable
iswritable(file) # yields whether file is writable
```

## Indexing and searching HDUs

An instance, say `file`, of `FitsFile` can be thought as an abstract vector of FITS *Header
Data Units* (HDUs) indexed by their number: `file[1]` is the primary HDU, `file[2]` is the
second HDU, `file[end]` is the last HDU, and so on.

The index may also be the HDU name (given by the `HDUNAME` or `EXTNAME` keywords)
or a predicate function. For example:

```julia
hdu = file[x -> nameof(x) == "SOME_NAME"]
hdu = file["SOME_NAME"]
```

both yield the first HDU in `file` whose name is `"SOME_NAME"`. In the latter case, the
comparison is done ignoring letter case as assumed by FITS standard for comparing strings.

To find the index of the first (resp. last) HDU matching `pat` call:

```julia
i = findfirst(pat, file)
i = findlast(pat, file)
```

where `pat` is a string, a regular expression or a predicate function. On return, `i` is
`nothing` if no match is found and an integer index otherwise. Then, to find the next (resp.
previous) HDU matching `pat` after (resp. before) `start`, call:

```julia
i = findnext(pat, file, start)
i = findprev(pat, file, start)
```

which also return `nothing` or an integer. The method `eachmatch` may be used to execute
some code on each HDU matching `pat`:

```julia
for hdu in eachmatch(pat, file)
    ... # do something
end
```

To check whether `ext` (a string or an integer) is an extension is `file`, call:

```julia
bool = haskey(file::FitsFile, ext::Union{AbstractString,Integer})
```

To retrieve the corresponding HDU (or `def` if the extension is not found):

```julia
hdu = get(file::FitsFile, ext::Union{AbstractString,Integer}, def)
```

## Stream-like operations

A `FitsFile` instance `file` may also be thought as a stream of HDUs:

```julia
hdutype = seek(file, i)   # move to `i`-th HDU of file
hdutype = seekstart(file) # move to first HDU of file
hdutype = seekend(file)   # move to first HDU of file
hdunum = position(file)   # get index of current HDU in file
flush(file)               # flush unwritten content to file
```

where `hdutype` is one of the following enumeration value:

* `FITS_IMAGE_HDU` if the current HDU contains an image.

* `FITS_BINARY_TABLE_HDU` if the current HDU contains a binary table.

* `FITS_ASCII_TABLE_HDU` if the current HDU contains an ASCII table.

* `FITS_ANY_HDU` if the current HDU has undefined content.

Note that `flush` may also be applied to a FITS HDU.
