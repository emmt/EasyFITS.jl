# Writing FITS files

This section explains how to directly write a FITS file in a single function call whose
arguments readily reflects the structure of the FITS file or how to create and write the
contents of a FITS file *piece-by-piece*.


## Direct writing of a FITS file

A complex (i.e. with multiple HDUs) FITS file may be directly written in a single call
to [`writefits`](@ref) or [`writefits!`](@ref) as:

```julia
writefits(filename, hdr1, dat1, hdr2, dat2, ...; kwds...)
```

where `filename` is the name of the FITS file, `hdr1` and `dat1` specify the header and
data parts of the first HDU, `hdr2` and `dat2` specify the header and data parts of the
second HDU and so on. Here `kwds...` represent keywords that are passed to
[`FitsFile`](@ref). The call `writefits!(args...; kwds...)` is a shortcut for
`writefits(args...; overwrite = true, kwds...)` which overwrites `filename` if it already
exists.

In `EasyFITS`, there are many different possible ways to specify a HDU header and data:

* A header may be `nothing` if there are no additional keywords other than the *structural
  keywords* describing the data part. Otherwise, the possible types for a header are given
  by the union [`EasyFITS.Header`](@ref): a header may be a collection of `FitsCard`
  instances (i.e. an instance of `FitsHeader`, a tuple, or a vector of FITS cards), a
  named tuple, a tuple or a vector of pairs like `key => val` or `key => (val, com)` with
  `key` the keyword name, `val` its value, and `com` a comment.

* For a FITS Image HDU, the data part is specified as a numerical Julia array.

* For a FITS Table HDU, the data part may be specified by a dictionary (whose keys are the
  column names and whose values are the corresponding column data), a vector of column
  data, a named tuple, a tuple or a vector of pairs like `col => val` with `col` the
  column name and `val` the column data. The data of each column of a table is a Julia
  array, the number of rows of the table is the last dimension of these arrays which must
  be the same for all columns.

Calling [`writefits`](@ref) or [`writefits!`](@ref), the structure of the resulting FITS
file can be made obvious by the code as in the following complex example with array `arr`
saved in the primary HDU and two tables in the next HDUs:

```julia
using Dates, EasyFITS
filename = "/tmp/test.fits";
arr = rand(Float32, (3,4,5));
nrows = 20;
inds = 1:nrows;
speed = rand(Float64, nrows);
mass = rand(Float32, nrows);
position = rand(Float32, 3, nrows);
phase = (1:7) .// 3;
amplitude = exp.(-1:-1:-7);
x = amplitude.*cos.(phase);
y = amplitude.*sin.(phase);
writefits(filename,
          #-----------------------------------------------------------------
          # First HDU must be a FITS "image", but data may be empty.
          #
          # Header part as a vector of `key=>val` or `key=>(val,com)` pairs:
          ["DATE"    => (now(), "date of creation"),
           "HISTORY" => "This file has been produced by EasyFITS",
           "USER"    => ENV["USER"]],
          # Data part as an array:
          arr,
          #-----------------------------------------------------------------
          # Second HDU, here a FITS "table".
          #
          # Header part of 2nd HDU as a tuple of pairs:
          ("EXTNAME" => ("MY-EXTENSION", "Name of this extension"),
           "EXTVER"  => (1, "Version of this extension")),
          # Data part is a table in the form of a named tuple:
          (Speed    = (speed, "km/s"),  # this column has units
           Indices  = inds,             # not this one
           Mass     = (mass, "kg"),
           Position = (position, "cm")),
          #-----------------------------------------------------------------
          # Third HDU, another FITS "table".
          #
          # Header part of 3rd HDU as a named tuple (note that keywords must
          # be in uppercase letters):
          (EXTNAME = ("MY-OTHER-EXTENSION", "Name of this other extension"),
           EXTVER  = (1, "Version of this other extension"),
           COMMENT = "This is an interesting comment"),
          # Data part is a table in the form of a vector of pairs (column names
          # can be strings or symbols but not a mixture):
          [:phase => ((180/π).*phase, "deg"),
           :amplitude => (amplitude, "V"),
           :xy => (hcat(x,y)', "V")])
```


## Advanced writing of a FITS file

Direct writing of a FITS file with [`writefits`](@ref) or [`writefits!`](@ref) is designed
to write all the content of a FITS file in a single function call and with a syntax that
readily shows the structure of the file. This is not suitable when not all content is
immediately available or when it is more convenient to write the FITS file piece by piece.
The latter is typically done by the following steps:

1. Create a new empty FITS file for writing with the [`FitsFile`](@ref) constructor.

2. Append a new HDU to with the [`FitsImageHDU`](@ref) constructor or the
   [`FitsTableHDU`](@ref) constructor.

3. Instantiate the FITS keywords of the header part of the HDU.

4. Write the data part of the HDU. This may be done by chunks and this
   depend on whether the HDU is an image or a table.

4. Eventually close the FITS file with `close`. Closing the FITS file is automatically done
   when the object representing the file is no longer used and garbage collected, closing
   the file is therefore optional.

Remarks:

- Steps 2 to 4 can be repeated as many times as necessary to create more than one HDU.


For example, with `arr` a Julia array to be written in the data part of the HDU:

```julia
# Create the file:
file = FitsFile(filename, "w")
# Appends a new HDU to store array data `arr`:
hdu = FitsImageHDU(file, arr)
# Set some keywords in the header part of the HDU:
hdu["HDUNAME"] = ("SOME_NAME", "Name of this HDU")
hdu["COMMENT"] = "Some comment."
# Write the data part of the HDU:
write(hdu, arr)
# Close the file:
close(file)
```
