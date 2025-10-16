# Aliases used for sub-indexing.
const IndexRange = OrdinalRange{<:Integer,<:Integer}
const SubArrayIndex = Union{Colon,Integer,IndexRange}
const SubArrayIndices{N} = NTuple{N,SubArrayIndex}

# Alias for specifying array size. C functions directly take vectors of integers (if of
# the correct type) and tuples of integers (if passed by reference).
const DimsLike = Union{Tuple{Vararg{Integer}},AbstractVector{<:Integer}}

if !isdefined(Base, :Memory)
    const Memory{T} = Vector{T}
end

"""
    EasyFITS.Header

is the (union of) type(s) that are accepted to specify a FITS header in
`EasyFITS` package.

A header may be a vector of pairs like `key => val`, `key => (val,com)`, or `key => com`
with `key` the keyword name, `val` its value, and `com` its comment. The keyword name `key`
is a string or a symbol which is automatically converted to uppercase letters and trailing
spaces discarded. The syntax `key => com`, with `com` a string, is only allowed for
commentary keywords `COMMENT` or `HISTORY`. For other keywords, the value is mandatory but
the comment is optional, not specifying the comment is like specifying `nothing` for the
comment; otherwise, the comment must be a string. The value `val` may be `missing` or
`undef` to indicate that it is undefined. If the comment is too long, it is automatically
split across multiple records for commentary keywords and it is truncated for other
keywords. A non-commentary keyword may have units specified in square brackets at the
beginning of the associated comment. Commentary keywords may appear more than once, other
keywords are unique.

For example:

    ["VERIFIED" => true,
     "COUNT" => (42, "Fundamental number"),
     "SPEED" => (2.1, "[km/s] Speed of gizmo"),
     "USER" => "Julia",
     "AGE" => (undef, "Some undefined value."),
     "JOB" => (missing, "Another undefined value."),
     "HISTORY" => "Some historical information.",
     "COMMENT" => "Some comment.",
     "COMMENT" => "Some other comment.",
     "HISTORY" => "Some other historical information."]

defines a possible FITS header with several records: a keyword `VERIFIED` having a logical
value and no comments, a keyword `COUNT` having an integer value and a comment, a keyword
`SPEED` having a floating-point value and a comment with units, a keyword `USER` having a
string value, keywords `AGE` and `JOB` having comments but undefined values, and a few
additional commentary keywords.

A header may also be specified as a named tuple with entries `key = val`, `key = (val,com)`,
or `key = com`. The same rules apply as above except that `key` must be allowed as a
variable symbolic name (no embedded hyphen `'-'`).

Finally, most methods assume that `nothing` can be used to indicate an empty header.

!!! note
    Specifying a FITS header as a dictionary is purposely not implemented because, to a
    certain extend, the order of keywords in a FITS header is relevant and because some
    keywords (`COMMENT`, `HISTORY`, and `CONTINUE`) may appear more than once.

"""
const Header = Union{
    # Collections of FITS cards:
    FitsHeader, Tuple{Vararg{FitsCard}}, AbstractVector{FitsCard},
    # Collection of pairs:
    NamedTuple, Tuple{Vararg{CardPair}}, AbstractVector{<:CardPair}}

"""
    EasyFITS.OptionalHeader

is the union of `Nothing` and of type(s) that are accepted to specify a FITS header in
`EasyFITS` package.

"""
const OptionalHeader = Union{Nothing,Header}

"""
    EasyFITS.ImageData{T,N}

is the possible type(s) for the pixels of a `N`-dimensional FITS image extension with
elements of type `T`.

"""
const ImageData{T<:Number,N} = AbstractArray{T,N}

"""
    EasyFITS.ColumnName

is the union of possible types for specifying the name of a column in a FITS table
extension.

"""
const ColumnName = Union{AbstractString,Symbol}

"""
    EasyFITS.ColumnIdent

is the union of possible types for identifying a single column in a FITS table extension.

"""
const ColumnIdent = Union{ColumnName,Integer}

"""
    EasyFITS.Columns

is the union of possible types for specifying several columns in a FITS table extension.

The method:

    EasyFITS.columns_to_read(hdu::FitsTableHDU, cols::Columns)

yields an iterable object over the column indices to read in table.

"""
const Columns = Union{Colon,
                      OrdinalRange{<:Integer,<:Integer},
                      AbstractVector{<:ColumnIdent},
                      # If specified as a tuple, don't mix indices and names:
                      Tuple{Vararg{Integer}},
                      Tuple{Vararg{ColumnName}}}

"""
    EasyFITS.Rows

is the union of possible types for specifying one or several rows in a FITS table
extension.

The following methods are provided:

    EasyFITS.rows_to_read(hdu::FitsTableHDU, rows::Rows)
    EasyFITS.first_row_to_read(hdu::FitsTableHDU, rows::Rows)
    EasyFITS.last_row_to_read(hdu::FitsTableHDU, rows::Rows)

to yield a iterable object over the row indices to read in table and the first/last such row
indices.

"""
const Rows = Union{Colon,Integer,AbstractUnitRange{<:Integer}}

"""
    EasyFITS.ColumnData{T,N}

is the possible type(s) for the cells of a `N`-dimensional column with values of type `T` in
a FITS table extension.

"""
const ColumnData{T,N} = AbstractArray{T,N}

"""
    EasyFITS.ColumnUnits

is the possible type(s) for the units of a column in a FITS table extension.

"""
const ColumnUnits = AbstractString

"""
    EasyFITS.ColumnEltype

is the possible type(s) for specifying the types of the values in a column of a FITS table
extension.

"""
const ColumnEltype = Union{Type,Char}

"""
    EasyFITS.ColumnDims

is the possible type(s) for specifying the cell dimensions in a column of a FITS table
extension.

"""
const ColumnDims = Union{Integer,Tuple{Vararg{Integer}}}

"""
    EasyFITS.ColumnSpec

is the possible type(s) for specifying the type of values, cell dimensions, and units of a
column of a FITS table extension. The element type is mandatory and must be specified first.

"""
const ColumnSpec = Union{ColumnEltype,
                         Tuple{ColumnEltype},
                         Tuple{ColumnEltype,ColumnUnits},
                         Tuple{ColumnEltype,ColumnDims},
                         Tuple{ColumnEltype,ColumnDims,ColumnUnits},
                         Tuple{ColumnEltype,ColumnUnits,ColumnDims}}

"""
    EasyFITS.ColumnNameSpecPair

is the possible type(s) for specifying a column when creating a table. Instances of this
kind are pairs like `col => type` or `col => (type,dims,units)` with `col` the column name
number, `type` the type of the column values, `dims` the cell dimensions, and `units` the
units of the values. `dims` and `units` are optional and may appear in any order after
`type` which is mandatory.

"""
const ColumnNameSpecPair = Pair{<:ColumnName,<:ColumnSpec}

"""
    EasyFITS.ColumnIdentDataPair

is the possible type(s) for specifying a column with its data and, optionally, its units to
be written in a FITS table extension. Instances of this kind are pairs like `col => vals` or
`col => (vals, units)` with `col` the column name or number, `vals` the column values, and
`units` optional units.

[`EasyFITS.ColumnNameDataPair`](@ref) is similar except that `col` cannot be a column
number.

"""
const ColumnIdentDataPair = Pair{<:ColumnIdent,
                                 <:Union{ColumnData,Tuple{ColumnData,ColumnUnits}}}

"""
    EasyFITS.ColumnNameDataPair

is the possible type(s) for specifying a column with its data and, optionally, its units to
be written in a FITS table extension. Instances of this kind are pairs like `col => vals` or
`col => (vals, units)` with `col` the column name, `vals` the column values, and `units`
optional units.

[`EasyFITS.ColumnIdentDataPair`](@ref) is similar except that `col` can be a column number.

"""
const ColumnNameDataPair = Pair{<:ColumnName,
                                <:Union{ColumnData,Tuple{ColumnData,ColumnUnits}}}

"""
    EasyFITS.TableData

is the union of types that can possibly be that of FITS table data. Instances of this kind
are collections of `key => vals` or `key => (vals, units)` pairs with `key` the column name,
`vals` the column values, and `units` the optional units of these values. Such collections
can be dictionaries, named tuples, vectors, or tuples.

For table data specified by dictionaries or vectors, the names of the columns must all be
of the same type.

Owing to the variety of possibilities for representing column values with optional units,
`EasyFITS.TableData` cannot be specific for the values of the pairs in the collection. The
package therefore rely on *error catcher* methods to detect column with invalid associated
data.

Another consequence is that there is a non-empty intersection between `EasyFITS.TableData`
and `EasyFITS.Header` which imposes to rely on position of arguments to distinguish them.

"""
const TableData = Union{AbstractDict{<:ColumnName,<:Any},
                        AbstractVector{<:Pair{<:ColumnName,<:Any}},
                        Tuple{Vararg{Pair{<:ColumnName,<:Any}}},
                        NamedTuple}

"""
    FitsHDU

is the abstract type of FITS Header Data Units which consists in a header and a data parts.
Concrete instances of `FitsHDU` behave as vectors whose elements are FITS header records,
a.k.a. FITS cards, and which can be indexed by integers or by names.

For faster access to the records of a header, consider creating a FITS header object from a
HDU object:

    hdr = FitsHeader(hdu::FitsHDU)

"""
abstract type FitsHDU <: AbstractVector{FitsCard} end

# Enumeration of HDU type identifiers.
@enum FitsHDUType::Cint begin
    FITS_IMAGE_HDU = CFITSIO.IMAGE_HDU
    FITS_BINARY_TABLE_HDU = CFITSIO.BINARY_TBL
    FITS_ASCII_TABLE_HDU = CFITSIO.ASCII_TBL
    FITS_ANY_HDU = CFITSIO.ANY_HDU
end

mutable struct FitsFile <: AbstractVector{FitsHDU}
    handle::Ptr{CFITSIO.fitsfile}
    mode::Symbol
    path::String
    nhdus::Int
    function FitsFile(filename::AbstractString, mode::AbstractString = "r";
                      extended::Bool=false)
        access = mode == "r" ? :r :
            mode == "r+" ? :rw :
            mode == "w" || mode == "w!" ? :w :
            throw(ArgumentError("mode must be \"r\", \"r+\", \"w\" or \"w!\""))
        status = Ref{Cint}(0)
        handle = Ref{Ptr{CFITSIO.fitsfile}}(0)
        num = Ref{Cint}(0)
        if access === :w
            if extended
                prefix = (mode == "w!" && ! startswith(filename, '!') ? "!" : "")
                check(CFITSIO.fits_create_file(handle, prefix*filename, status))
            else
                if isfile(filename)
                    mode == "w!" || throw_file_already_exists(filename, "use mode \"w!\" to overwrite")
                    rm(filename; force=true)
                end
                check(CFITSIO.fits_create_diskfile(handle, filename, status))
            end
        else
            iomode = (access === :rw ? CFITSIO.READWRITE : CFITSIO.READONLY)
            if extended
                check(CFITSIO.fits_open_file(handle, filename, iomode, status))
            else
                check(CFITSIO.fits_open_diskfile(handle, filename, iomode, status))
            end
            if !iszero(CFITSIO.fits_get_num_hdus(handle[], num, status))
                errcode = status[]
                status[] = 0
                CFITSIO.fits_close_file(ptr, status)
                throw(FitsError(errcode))
            end
        end
        return finalizer(close_handle, new(handle[], access, filename, num[]))
    end
end

"""
    EasyFITS.Invalid

is the singleton type of the object used to indicate invalid arguments while sparing
throwing an exception.

"""
struct Invalid end

# Singleton type to indicate that the inner constructor should be called.
struct BareBuild end

# Any other FITS extension than Image and Table (who knows...).
struct FitsAnyHDU <: FitsHDU
    file::FitsFile
    num::Int
    FitsAnyHDU(::BareBuild, file::FitsFile, num::Integer) = new(file, num)
end

# FITS Image extension (Array for Julia).
struct FitsImageHDU{T,N} <: FitsHDU
    file::FitsFile
    num::Int
    function FitsImageHDU{T,N}(::BareBuild, file::FitsFile, num::Integer) where {T,N}
        isbitstype(T) || bad_argument("parameter `T=$T` is not a plain type")
        isa(N, Int) && N ≥ 0 || bad_argument("parameter `N=$N` must be a nonnegative `Int`")
        return new{T,N}(file, num)
    end
end

# FITS Table extension.
struct FitsTableHDU <: FitsHDU
    file::FitsFile
    num::Int
    ascii::Bool
    FitsTableHDU(::BareBuild, file::FitsFile, num::Integer, ascii::Bool) =
        new(file, num, ascii)
end

"""
    FitsLogic()

yields a singleton object to indicate that FITS rules should be applied for some logical
operation. For example:

    isequal(FitsLogic(), s1, s2)

compares strings `s1` and `s2` according to FITS rules, that is case of letters and trailing
spaces are irrelevant.

    isequal(FitsLogic(), x) -> f

yields a predicate function `f` such that `f(y)` yields `isequal(FitsLogic(),x,y)`.

"""
struct FitsLogic end

struct Bit end

"""
    FitsError <: Exception

is the type of exceptions thrown when an error occurs in the FITSIO library used to manage
FITS files.

"""
struct FitsError <: Exception
    code::Cint
end
