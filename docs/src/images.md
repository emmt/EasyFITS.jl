# FITS image extensions

FITS image extensions store multi-dimensional arrays with numerical values
exactly as regular Julia arrays. In `EasyFITS`, a FITS image extension in an
open file is represented by an object of type `FitsImageDHU{T,N}` with `T` the
element type and `N` the number of dimensions.


## Image HDU properties

An image HDU has the following properties:

| Property      | Description                     |
|:--------------|:--------------------------------|
| `data_size`   | Array dimensions                |
| `data_axes`   | Array axes                      |
| `data_ndims`  | Number of dimensions `N`        |
| `data_eltype` | Element type `T`                |
| `extname`     | Extension name                  |
| `hduname`     | HDU name                        |
| `file`        | Associated i/o FITS file        |
| `number`      | HDU number                      |
| `type`        | HDU type, i.e. `FITS_IMAGE_HDU` |
| `xtension`    | Extension, i.e. `"IMAGE"`       |


## Reading a FITS image

To read the data stored by the *Header Data Unit* (HDU) object `hdu` of type
`FitsImageDHU` HDU as an array `arr`, call [`read`](@ref) as:

``` julia
arr = read(hdu)
```

The result is of type `Array{T,N}` if `hdu isa FitsImageDHU{T,N}` holds. To
choose another element type, say, `S`, just do:

``` julia
arr = read(Array{S}, hdu)
```

This is similar but more efficient than any of:

``` julia
arr = S.(read(hdu))
arr = convert(Array{S}, read(hdu))
```

For type-stability, the expected number of dimension, say `N`, may also be
specified. For example, any of:

``` julia
arr = read(Array{Float32,2}, hdu)
arr = read(Matrix{Float32}, hdu)
```

ensure that `arr` will be a 2-dimensional image with pixels of type `Float32`.

Call [`read!`](@ref) to overwrite the elements of an existing array with the
contents of the FITS image extension. For example:

``` julia
read!(arr, hdu)
```

An hyper-rectangular sub-image can be read using the same syntax as for a Julia
`view` by specifying indices and/or index ranges after the `hdu` argument.
Index ranges may have non-unit steps but steps must all be positive. For
example:

``` julia
arr = read(hdu, :, :, 2)
```

yields the 2nd slice in a 3-dimensional FITS image.

The result is similar to:

``` julia
arr = read(R, hdu)[:,:,2]
```

but should be more efficient as no array other than the result is allocated and
fewer values are read.

Call `read!` instead of `read` to overwrite the contents of an existing array.
Following the previous example, reading the next slice could be done by:

``` julia
read!(arr, hdu, :, :, 3)
```


## Creating an image extension

To create a new FITS image extension in an open FITS file

``` julia
write(file::FitsFile, FitsImageHDU{T}, dims...=()) -> hdu
write(file::FitsFile, FitsImageHDU, T::Type=UInt8, dims...=()) -> hdu
write(file::FitsFile, FitsImageHDU, bitpix::Integer, dims=()) -> hdu
```

create a new primary array or image extension in FITS file `file` with a
specified pixel type `T` and size `dims...`. If the FITS file is currently
empty then a primary array is created, otherwise a new image extension is
appended to the file. Pixel type can be specified as a numeric type `T` or as
an integer BITPIX code `bitpix`.

An object to manage the new extension is returned which can be used to push
header cards and then to write the data.

For example:

``` julia
hdu = write(file, FitsImageHDU, eltype(arr), size(arr))
hdu["KEY1"] = val1             # add a 1st header record
hdu["KEY2"] = (val2, str2)     # add a 2nd header record
hdu["KEY3"] = (nothing, str3)  # add a 3rd header record
write(hdu, arr)                # write data
```

will create a new Header Data Unit (HDU) storing array `arr` with 3 additional
header records: one named `"KEY1"` with value `val1` and no comments, another
named `"KEY2"` with value `val2` and comment string `str2`, and yet another one
named `"KEY3"` with no value and with comment string `str3`. Note that special
names `"COMMENT"`, `"HISTORY"`, and `""` indicating commentary entries have no
associated, only a comment string, say `str` which can be specified as `str` or
as `(nothing,str)`.
