# Notes for developers

## The FITS standard

The FITS format is described in the [*FITS Standard
Document*](https://fits.gsfc.nasa.gov/fits_standard.html).


## The `CFITSIO` sub-module

`EasyFITS` makes use of the [`Clang`](https://github.com/JuliaInterop/Clang.jl)
Julia package to automatically build file `deps/deps.jl` with constants, types,
and low level functions to call the functions of the CFITSIO library with
arguments of the correct type. All these are available in the
`EasyFITS.CFITSIO` sub-module.

## Calls to functions in the CFITSIO library

When calling functions of the CFITSIO library, there are several things to take
care of:

- Passing correct arguments. This is partially ensured by the type assertions
  in the `@ccall` macro. It is also necessary to check whether a pointer to
  some opaque structure in the library is valid.

- Preserving objects from being destroyed while being in use. Of course, this
  is automatically done by Julia for Julia objects, but must be handled for
  references or pointers to objects provided by the CFITSIO library.

When an object `obj` is specified for an argument of type `Ptr{T}` to be passed
to a C function, Julia `ccall` does something like:

``` julia
ref = Base.cconvert(Ptr{T}, obj)
ptr = Base.unsafe_convert(Ptr{T}, ref)
result = GC.@preserve ref call_some_c_function(..., ptr, ...)
```

Here `ref` is an object (by default, `Base.cconvert(Ptr{T},obj)` yields `obj`
itself) to be used with `Base.unsafe_convert(Ptr{T},ref)` to get the pointer
and to be preserved from being garbage collected in order to warrant that the
pointer remains valid.

Thanks to this mechanism, it is quite simple to ensure that valid pointers to
opaque structures of the CFITSIO library be passed to a function of this
library. For example, the code below is how is handled a pointer to a
`fitsfile` C structure in our code:

``` julia
isnull(ptr::Ptr{T}) where {T} = ptr === Ptr{T}(0)
check(ptr::Ptr) = isnull(ptr) ? error("invalid NULL pointer") : ptr
get_handle(file::FitsFile) = getfield(obj, :handle)
get_file(hdu::AbstractHDU) = getfield(hdu, :file)
Base.unsafe_convert(Ptr{CFITSIO.fitsfile}, obj::FitsFile) = check(get_handle(obj))
Base.cconvert(Ptr{CFITSIO.fitsfile}, hdu::AbstractHDU) = get_file(hdu)
```

Private methods `isnull` and `check` are introduced for readability. Private
methods `get_handle` and `get_file` are two of the private accessors introduced
to hide the fields of a FITS file object and let some other *public* properties
be implemented. The former yields the pointer to the `fitsfile` C structure
that is managed by a `FitsFile` object, while the latter yields the `FitsFile`
object storing the HDU as it is the one that must be used and preserved when
calling a C function requiring a pointer to a `fitsfile` C structure.

Note that `FitsFile` objects have a finalizer that automatically releases
resources such as the associated `fitsfile` C structure when the object is
garbage collected.

## Helper functions

The following non-exported functions are provided for meta-programming and for
dealing with the types of arguments in calls to the functions of the CFITSIO
library.

```@docs
EasyFITS.cfunc
EasyFITS.ctype
EasyFITS.cpointer
```

## Pixel types

The following table lists conventions used by `CFITSIO` for pixel types, that
is the `BITPIX` keyword in FITS image extensions.

| Type Code       | Julia Type | `BITPIX` |
|:----------------|:-----------|---------:|
| `BYTE_IMG`      | `UInt8`    |        8 |
| `SBYTE_IMG`     | `Int8`     |          |
| `SHORT_IMG`     | `Int16`    |       16 |
| `USHORT_IMG`    | `UInt16`   |          |
| `LONG_IMG`      | `Int32`    |       32 |
| `ULONG_IMG`     | `UInt32`   |          |
| `LONGLONG_IMG`  | `Int64`    |       64 |
| `ULONGLONG_IMG` | `UInt64`   |          |
| `FLOAT_IMG`     | `Float32`  |      -32 |
| `DOUBLE_IMG`    | `Float64`  |      -64 |

Types without a value in the `BITPIX` column are converted by the CFITSIO
library into the other signed/unsigned type using special values of the
`BSCALE` and `BZERO` keywords to allow for the reciprocal conversion. This is
explicitly allowed by the FITS Standard (version 4.0).

The above equivalence rules are implemented by the following two non-exported
functions.

```@docs
EasyFITS.type_to_bitpix
EasyFITS.type_from_bitpix
```

## Array data types

The following table lists conventions used by `CFITSIO` for array element
types.

| Type Code     | C Type               | Short Suffix | Long Suffix |
|:--------------|:---------------------|:-------------|:------------|
| `TLOGICAL`    | `char`               | `l`          | `_log`      |
| `TBYTE`       | `unsigned char`      | `b`          | `_byt`      |
| `TSBYTE`      | `signed char`        | `sb`         | `_sbyt`     |
| `TUSHORT`     | `unsigned short`     | `ui`         | `_usht`     |
| `TSHORT`      | `short`              | `i`          | `_sht`      |
| `TUINT`       | `unsigned int`       | `uk`         | `_uint`     |
| `TINT`        | `int`                | `k`          | `_int`      |
| `TULONG`      | `unsigned long`      | `uj`         | `_ulng`     |
| `TLONG`       | `long`               | `j`          | `_lng`      |
| `TULONGLONG`  | `unsigned long long` | `ujj`        | `_ulnglng`  |
| `TLONGLONG`   | `long long`          | `jj`         | `_lnglng`   |
| `TFLOAT`      | `float`              | `e`          | `_flt`      |
| `TDOUBLE`     | `double`             | `d`          | `_dbl`      |
| `TCOMPLEX`    | `float complex`      | `c`          | `_cmp`      |
| `TDBLCOMPLEX` | `double complex`     | `m`          | `_dblcmp`   |
| `TSTRING`     | `char*`              | `s`          | `_str`      |
| `TBIT`        |                      | `x`          | `_bit`      |
|               |                      | `u`          | `_null`     |

Complex types `float complex` and `double complex` are stored as pairs of
single/double precision floating-point values (this is not guaranteed by C99
standard so strict equivalence does not hold here).

The above equivalence rules are implemented by the following two non-exported
functions.

```@docs
EasyFITS.type_to_code
EasyFITS.type_from_code
```


## Column data types

The following table lists the correspondences between the `TFORMn` letter in
FITS table extensions and the column data type.

| Type Code     | Julia Type   | `TFORM` | Description                          |
|:--------------|:-------------|:--------|:-------------------------------------|
| `TLOGICAL`    | `Bool`       | `’L’`   | Logical (1 byte)                     |
| `TBIT`        |              | `’X’`   | Bit (special)                        |
| `TBYTE`       | `UInt8`      | `’B’`   | 8-bit unsigned integer               |
| `TSHORT`      | `Int16`      | `’I’`   | 16-bit signed integer                |
| `TLONG`       | `Int32`      | `’J’`   | 32-bit signed integer                |
| `TLONGLONG`   | `Int64`      | `’K’`   | 64-bit signed integer                |
| `TSTRING`     | `String`     | `’A’`   | Character (1 byte, used for strings) |
| `TFLOAT`      | `Float32`    | `’E’`   | 32-bit floating point                |
| `TDOUBLE`     | `Float64`    | `’D’`   | 64-bit floating point                |
| `TCOMPLEX`    | `ComplexF32` | `’C’`   | 64-bit complex                       |
| `TDBLCOMPLEX` | `ComplexF64` | `’M’`   | 128-bit complex                      |
|               |              | `’P’`   | 32-bit array descriptor              |
|               |              | `’Q’`   | 64-bit array descriptor              |

A few non-standard `TFORM` letters are allowed by CFITSIO. These are converted
by the library into other types using `TSCALE` and `TZERO` keywords to allow
for the reciprocal conversion following the same principles as for the `BITPIX`
code and the `BSCALE` and `BZERO` keywords.

| Type Code    | Julia Type | `TFORM` | Description             |
|:-------------|:-----------|:--------|:------------------------|
| `TSBYTE`     | `Int8`     | `’S’`   | 8-bit signed integer    |
| `TUSHORT`    | `UInt16`   | `’U’`   | 16-bit unsigned integer |
| `TULONG`     | `UInt32`   | `’V’`   | 32-bit unsigned integer |
| `TULONGLONG` | `UInt64`   | `’W’`   | 64-bit unsigned integer |

 The *Type Code* column indicates the code used by CFITSIO (it is not always
consistent with the C types as defined in the above table, so my guess is that
this code is only used to keep track of the column data type internally).

The above equivalence rules are implemented by the following two non-exported
functions.

```@docs
EasyFITS.type_to_letter
EasyFITS.type_from_letter
```
