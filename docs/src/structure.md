# Structure of FITS files

A FITS file is a collection of *Header Data Units* (HDUs). The first HDU of a FITS file is
called the *Primary HDU*, other HDUs are called *extenstions*. Each HDU consists in two
parts: the header part is a list of FITS keywords which have a name and, usually, a value
and a comment; the data part may be a multi-dimensional array or a table. The *primary
HDU* and so-called *FITS Image extensions* are HDUs storing a multi-dimensional array in
their data part, while *FITS Table extensions* are HDUs storing a table in their data
part. The HDUs of a FITS file are indexed by their number (starting at 1 for the *primary
HDU*) or by their name (specified by the value of the `HDUNAME` or `EXTNAME` keyword of
their header part). The first HDU which must be an image may be empty. Arrays stored in
FITS Image extensions have elements of the same numeric type and may have 1 to 999
dimensions. Tables stored in FITS Image extensions may have a variable number of columns
identified by their names, their entries (at a given row and column) are called *cells*
and may have multiple dimensions. Generally, all cells of a given column have the same
data type and dimensions. FITS headers are stored in textual form encoded with a
restricted subset of the ASCII characters. Although FITS supports ASCII Tables, the data
part of FITS HDUs is generally in binary format (in big endian byte order and following
IEEE stabdard for floating point values).

References and related resources:

* The [IAU FITS Working Group](http://fits.gsfc.nasa.gov/iaufwg/) maintains the FITS
  standard and endorses its evolution.

* The [NASA FITS Standard Document](https://fits.gsfc.nasa.gov/fits_standard.html)
  hosts the offical reference document that defines the FITS format.
