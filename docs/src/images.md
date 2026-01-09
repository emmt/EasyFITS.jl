# FITS Image HDUs

FITS Image HDUs store multi-dimensional arrays with numerical values exactly as regular
Julia arrays (of type `Array`). In `EasyFITS`, a FITS Image HDU is represented by an object
of type [`FitsImageDHU{T,N}`](@ref FitsImageHDU) with `T` the element type and `N` the
number of dimensions.


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
| `xtension`    | Extension name, i.e. `"IMAGE"`  |


## Reading a FITS Image

To read the data stored by the *Header Data Unit* (HDU) object `hdu` of type `FitsImageDHU`
HDU as an array `arr`, call [`read`](@ref read(::FitsImageHDU)) as:

``` julia
arr = read(hdu)
```

The result is type-stable: `arr` is of type `Array{T,N}` if `hdu isa FitsImageDHU{T,N}`
holds. To choose another element type, say `S`, just do:

``` julia
arr = read(Array{S}, hdu)
```

This is similar but more efficient than any of:

``` julia
arr = S.(read(hdu))
arr = convert(Array{S}, read(hdu))
```

The expected number of dimension, say `N`, may also be specified. For example:

``` julia
arr = read(Array{Float32,2}, hdu)
arr = read(Matrix{Float32}, hdu)
```

both warrant that `arr` will be a 2-dimensional image with pixels of type `Float32`.

Call [`read!`](@ref) to overwrite the elements of an existing array with the contents of
the FITS Image HDU. For example:

``` julia
read!(arr, hdu)
```

An hyper-rectangular sub-image can be read using the same syntax as for a Julia `view` by
specifying indices, index ranges, or colons after the `hdu` argument. Index ranges may have
non-unit steps but steps must all be positive. For example:

``` julia
arr = read(hdu, :, :, 2)
```

yields the 2nd slice in a 3-dimensional FITS Image.

The result is similar to:

``` julia
arr = read(hdu)[:,:,2]
```

but should be more efficient as no array other than the result is allocated and fewer values
are read.

Call `read!` instead of `read` to overwrite the contents of an existing array. Following the
previous example, reading the next slice could be done by:

``` julia
read!(arr, hdu, :, :, 3)
```


## Creating an image HDU

To start a new FITS Image HDU in an open FITS `file`, there are several possibilities:

``` julia
hdu = FitsImageHDU(file, dims...; bitpix=...)
hdu = FitsImageHDU{T}(file, dims...)
hdu = FitsImageHDU{T,N}(file, dims...)
```

with pixel type specified by its Julia data-type via the parameter `T` or by its FITS BITPIX
code via the keyword `bitpix` and image size given by `dims...`. The number of dimensions
`N` can be inferred from the image size but may be explicitly specified for assertion.

If the FITS file is currently empty then a primary array is created, otherwise a new image
HDU is appended to the file.

The returned `hdu` is an object to manage the new HDU, it can be used to push header cards
and then to write the data.

If the array `arr` to be written is available, the element type and dimensions can be
inferred from `arr` itself:

``` julia
hdu = FitsImageHDU(file, arr)
hdu = FitsImageHDU{T}(file, arr)
hdu = FitsImageHDU{T,N}(file, arr)
```

where `T = eltype(arr)` and `N = ndims(arr)` are assumed if these parameters are not
explicitly specified. Specifying a different pixel type than `eltype(arr)` is possible,
conversion will automatically be performed when writing the array values.

After starting a new HDU, the non structural keywords of its header part of the HDU can be
written and then the data part of the HDU must be written (unless empty).

For example:

``` julia
hdu = FitsImageHDU(file, arr)
hdu["KEY1"] = val1             # add a 1st header record
hdu["KEY2"] = (val2, com2)     # add a 2nd header record
hdu["KEY3"] = (nothing, com3)  # add a 3rd header record
write(hdu, arr)                # write data
```

will create a new Header Data Unit (HDU) storing array `arr` with 3 additional header
records: one named `"KEY1"` with value `val1` and no comments, another named `"KEY2"` with
value `val2` and comment string `com2`, and yet another one named `"KEY3"` with no value and
with comment string `com3`. Note that special names `"COMMENT"`, `"HISTORY"`, and `""`
indicating commentary entries have no associated, only a comment string, say `com` which can
be specified as `com` or as `(nothing,com)`.
