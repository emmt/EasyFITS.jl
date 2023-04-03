# Access to FITS files

## Direct reading of data

The simplest way to read some data in a FITS file is to call
[`readfits`](@ref):

```julia
data = readfits(filename, args...; kwds...)
```

with `filename` the name of the file.

By default, the data of the first FITS extension in `filename` is read. Keyword
`ext` may however be set with a number, a name, or a symbol to select another
extension. Another possibility is to specify the keyword `extended = true` to
open the file using the [extended file name
syntax](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html)
implemented by the `CFITSIO` library.

The default is to read all the data part of the selected FITS extension but
optional arguments `args...` and keywords `kwds...` may be specified to
restrict the read contents:

- If the FITS extension is an image, `args...` specifies the ranges of pixels
  to read along the dimensions. For example:

  ```julia
  data = readfits(filename, :, :, 5)
  ```

  would read the 5th slice in a 3-dimensional image. This is equivalent but
  more efficient than:

  ```julia
  data = readfits(filename)[:, :, 5]
  ```

  which amounts to reading all the data and then only keep the 5th slice.

- If the FITS extension is a table, `args...` may be up to 2 arguments, `cols`
  and `rows`, to respectively select a subset of columns and rows. For example:

  ```julia
  A = readfits(filename, ("Speed", "Height"))
  B = readfits(filename, :, 11:40)
  ```

  respectively yield the columns named `Speed` and `Height` of the table and
  all the columns of the table but only for rows in the range `11:40`. Keyword
  `case` may be used to indicate whether letter case does matter in the column
  names. Note that, in `EasyFITS` the rows of a table correspond to the last
  dimension of arrays. This is to have the same storage order in memory and in
  the FITS file. Method `permutedims` can be used is this convention does not
  suit you.

The type of the object returned by [`readfits`](@ref) depends on the kind of
the FITS extension and may also depend on the optional arguments `args...` and
on the keywords `kwds...`:

- If the FITS extension is an image, the data is read as a Julia `Array`.

- If the FITS extension is a table, then, if `cols` is a single column name or
  number, an array of the columns values is returned, otherwise, a dictionary
  indexed by the column names is returned. Note that, a column range like `4:4`
  would yield a dictionary with a single column (the 4th one) in that context.
  Keywords `case` and `rename` are available to indicate how to search the
  columns by name in the table and how to translate these names into dictionary
  keys.

To avoid ambiguities or for improved type-stability, a leading type argument
can be specified in [`readfits`](@ref) to indicate the expected type for the
returned data:

```julia
data = readfits(R::Type, filename, args...; kwds...)::R
```

Array type parameters such as element type and number of dimensions may be
specified in `R`. For example:

```julia
data = readfits(Array{Float32,3}, filename)
```

ensures that the result be a single precision floating-point 3-dimensional
array; while:

```julia
data = readfits(Dict, filename, cols; ext=2)
```

ensures that the table in 2nd FITS Header Data Unit is returned as a dictionary
even though `cols` specifies a single column.


## Direct writing of data

In principle, directly writing data in a FITS file is as simple as direct
reading of data.
