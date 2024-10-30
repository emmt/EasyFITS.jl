# FITS Table extensions

In a given table column, cells may contain more than a single value.

!!! note
    The elements of a column of a FITS table extension are stored in the same memory order
    as in an ordinary Julia array. As a result, the **rows** of a column in a FITS table
    extension correspond to the **last** index in the equivalent Julia array. Method
    `permutedims` can be used is this convention does not suit you.


## Reading a single table column

To read a single column in a FITS table extension as an array `arr`, call `read` as:

``` julia
arr = read(hdu, col)
```

with `hdu` the *Header Data Unit* containing the table and `col` the column name or
number. The names of the columns are given by `hdu.column_names` and valid column numbers
are in the range `hdu.columns` (see [Table HDU Properties](#Table-HDU-Properties) for a
list of all properties). If `col` is a string or a symbol, keyword `case` can be used to
specify whether uppercase/lowercase matters (`case` is false by default).

By default, all rows are read but it is possible to specify which rows to read:

``` julia
arr = read(hdu, col, rows)
```

where `rows` can be an integer to read a single row, a unit range of integers to read
these rows, or a colon `:` to read all the rows (the default).
[Properties](#Table-HDU-Properties) `hdu.first_row` and `hdu.last_row` can be used to
retrieve the first and last row numbers of `hdu`.

The `read` method automatically guess the number of dimensions and the type of the
elements of the Julia array representing the column of a table. To improve type stability,
the `read` method takes an optional first argument to specify the type of the result, from
the least specific `Array` to `Array{T}` or `Array{T,N}`. For instance:

``` julia
arr = read(Vector{Float32}, hdu, col)
arr = read(Vector{Float32}, hdu, col, rows)
```

are both guaranteed to yield an instance of `Vector{Float32}` (an exception is thrown if
the elements of the column cannot be converted to the given type or if the column cells
are not scalar numbers).

In fact, the calls `read(hdu,col)` and `read(hdu,col,rows)` are respectively equivalent to
`read(Array,hdu,col)` and `read(Array,hdu,col,rows)`.

The `read!` method can be called to overwrite the contents of an existing array (it must
have contiguous entries) with the cells of a column:

``` julia
read!(arr, hdu, col)
```

where keyword `first` may be used to specify the first row to read (the first one,
`hdu.first_row`, by default), the number of rows to read being deduced from the size of
`arr`. Other keywords are `anynull` and `null` to deal with undefined values.


## Reading several table columns

In the most simple form, reading several table columns is done by one of:

``` julia
dict = read(hdu, cols)
dict = read(hdu, cols, rows)
```

where, as for [reading a single column](#Reading-a-single-table-column), `hdu` is the
*Header Data Unit* containing the table and `rows` may be specified to select the rows to
read but `cols` is a vector or a tuple of column names/numbers, or a colon `:` to read all
columns. Specify keyword `case=true` if columns are specified by their (symbolic) names
and the case of these names matters.

The result of reading several columns is a dictionary whose keys are the names of the
columns and whose values are the cells of the columns. The keys of the dictionary are the
names of the table columns unconverted if keyword `case=true` and converted to upper-case
letters otherwise. If this is not suitable, use the keyword `rename` to specify a function
to modify the names of the columns. For example, with `rename=identity` the names of the
columns will be left unchanged whatever the keyword `case` while with `rename=lowercase`
the names of the columns will be converted to lower-case letters.

In order to retrieve the units of the columns, specify keyword `units=String` and get a
dictionary whose values are 2-tuples of the form `(data,units)` with `data` the column
cells and `units` the column units as a string. The default behavior, that is to get a
dictionary whose values are simply be the cells of the columns, corresponds to
`units=nothing`.

To have a dictionary when `cols` is a single column name/number, just specify the type of
the expected result as the leading argument. For example:

``` julia
dict = read(Dict, hdu)
dict = read(Dict, hdu, cols)
dict = read(Dict, hdu, cols, rows)
```

Of course, this also works if the `cols` argument represents several columns without
ambiguities.

The leading type argument may be more specific. For example:

``` julia
dict = read(Dict{String,Array}, hdu, cols)
```

is the same as with `Dict`, while:

``` julia
dict = read(Dict{String,Vector}, hdu, cols, rows)
```

can only be successful if all read columns have 0-dimensional cells.

If the leading argument is more specific than `Dict`, the `units=String` keyword is not
allowed and retrieving the units is done by specifying that the values are tuples of an
array and a string. For example:

``` julia
dict = read(Dict{String,Tuple{Array,String}}, hdu, cols)
```

Using the leading type argument, it is also possible to retrieve a vector of the table
columns (in the same order as specified in `cols`). For example:

``` julia
vect = read(Vector, hdu, cols)
vect = read(Vector, hdu, cols, rows)
vect = read(Vector{Array}, hdu, cols, rows)
```

It is also possible to retrieve the units at the same time:

``` julia
vect = read(Vector, hdu, cols; units=String)
vect = read(Vector, hdu, cols, rows; units=String)
vect = read(Vector{Tuple{Array,String}}, hdu, cols, rows)
```

With the `read!` method, the contents of an existing dictionary may be replaced by columns
from a FITS table:

``` julia
read!(dict, hdu)
read!(dict, hdu, cols)
read!(dict, hdu, cols, rows)
```

where, as for the `read` method, keywords `case` and `rename` may be supplied. Call the
`merge!` method instead to preserve pre-existing contents:

``` julia
merge!(dict, hdu)
merge!(dict, hdu, cols)
merge!(dict, hdu, cols, rows)
```

With the `push!` method, a vector of columns can be augmented with additional columns:

``` julia
push!(vect, hdu)
push!(vect, hdu, cols)
push!(vect, hdu, cols, rows)
```


## Creating a table extension

To create a new FITS table extension in the FITS file `file`, call:

``` julia
hdu = write(file, FitsTableHDU, cols)
```

with `cols` defining the columns of the table. Each column definition is a pair `name =>
format` where `name` is the column name while `format` specifies the type of the column
values and, optionally, their units and the size of the column cells. The following
definitions are possible:

``` julia
name => type
name => (type,)
name => (type, units)
name => (type, dims)
name => (type, units, dims)
name => (type, dims, units)
```

where `type` is either a Julia type (`Number` or `String`) or a letter (see
table below), `dims` is an integer or a tuple of integers, and `units` is a
string. By default, `dims = (1,)` and `units = ""`.

| Type               | Letter | Remarks          |
|:-------------------|:-------|:-----------------|
| `AbstractString`   | `'A'`  | ASCII string     |
| `Bool`             | `'L'`  |                  |
| `Int8`             | `'S'`  | CFITSIO specific |
| `UInt8`            | `'B'`  |                  |
| `Int16`            | `'I'`  |                  |
| `UInt16`           | `'U'`  | CFITSIO specific |
| `Int32`            | `'J'`  |                  |
| `UInt32`           | `'V'`  | CFITSIO specific |
| `Int64`            | `'K'`  |                  |
| `UInt64`           | `'W'`  | CFITSIO specific |
| `Float32`          | `'E'`  |                  |
| `Float64`          | `'D'`  |                  |
| `Complex{Float32}` | `'C'`  |                  |
| `Complex{Float64}` | `'M'`  |                  |
| `FitsBit`          | `'X'`  |  bits

The returned object, `hdu`, can be used to add FITS keywords to the header of
the table and, then, to write column data. Typically:

```julia
hdu = write(file, FitsTableHDU, cols)
push!(hdu, key1 => val1)         # add a first header keyword
push!(hdu, key2 => (val2, com2)) # add a second header keyword with a comment
...                              # add other header keywords
write(hdu, col1 => arr1)         # write a first column
write(hdu, col2 => arr2)         # write a second column
...                              # write other columns
```

where `key1 => val1`, `key2 => (val2,com2)`, etc. specify header cards, while `col1 =>
arr1`, `col2 => arr2`, etc. specify columns names and associated data. Such a table may be
created in a single call:

```julia
write(file,
      [key1 => val1, key2 => (val2, com2), ...],
      [col1 => arr1, col2 => arr2, ...])
```

which follows the `write(dest, header, data)` convention in `EasyFITS` with `dest` the
destination, `header` the header, and `data` the data.

The header may be specified in [any forms](#Header_forms) accepted by the `EasyFITS`
methods. Similarly, the columns may be specified in various forms as explained below.


## Writing table columns

To write a single column into the FITS table extension `hdu`:

```julia
write(hdu, col => arr, ...; first=hdu.first_row, case=false, null=nothing) -> hdu
```

where `col` is the column name/number and `arr` is an array of column values. Column
values are converted as needed and are written starting at the row specified by `first`.
The leading dimensions of `arr` should be the same as those specified by the corresponding
`TDIMn` keyword (with `n` the column number) and the remaining last dimension, if any,
corresponds to the *row* index of the table.

If `col` is a column name, keyword `case` may be used to indicate whether case of letters
matters (default is `false`).

Keyword `null` may be used to specify the value of undefined elements in `arr`.

Any number of columns may be specified as subsequent arguments. The same keywords apply to
all columns.


## Table HDU Properties

The following table lists all properties of a FITS table HDU.

| Property        | Description                          |
|:----------------|:-------------------------------------|
| `nrows`         | Number of rows                       |
| `rows`          | Index range of rows                  |
| `first_row`     | Index of first row                   |
| `last_row`      | Index of last row                    |
| `ncols`         | Number of columns                    |
| `columns`       | Index range of columns               |
| `first_column`  | Index of first column                |
| `last_column`   | Index of last column                 |
| `column_name`   | Column name accessor                 |
| `column_names`  | Column names                         |
| `column_number` | Column number accessor               |
| `column_units`  | Column units accessor                |
| `data_size`     | Table dimensions                     |
| `data_ndims`    | Number of table dimensions           |
| `data_axes`     | Indices along table dimensions       |
| `extname`       | Extension name                       |
| `hduname`       | HDU name                             |
| `file`          | Associated FITS file                 |
| `num`           | HDU number                           |
| `type`          | HDU type: `FITS_BINARY_TABLE_HDU`    |
| `xtension`      | Extension: `"BINTABLE"` or `"TABLE"` |

To retrieve the units, the number, or the name of column `col` in the FITS
table `hdu` object, use the properties `column_units`, `column_number`, or
`column_name` properties as follows:

```julia
hdu.column_units(col; case=false) -> units::String
hdu.column_name(col; case=false) -> name::String
hdu.column_number(col; case=false) -> number::Int
```

Keyword `case` specifies whether the case of letters does matters when `col` is a
(symbolic) name. The result of `hdu.column_units(col)` is always a string, possibly empty.
