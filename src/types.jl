const Status = Cint

const OptionalString = Union{AbstractString,Nothing}

# Types indicating an undefined card value, Nothing is used for the value of
# commentary cards.
const UndefinedValue = Union{Missing,UndefInitializer}

# FITS header card name, a.k.a. FITS keyword.
const CardName = Union{Symbol,AbstractString}

# Acceptable types for the value of a FITS header card.
const CardValue = Union{Nothing,UndefinedValue,Real,Complex,AbstractString}

# The different ways to represent a FITS header card value + comment.
const CardData = Union{CardValue,
                       Tuple{CardValue},
                       Tuple{CardValue,OptionalString}}

const CardPair{K<:CardName,V<:CardData} = Pair{K,V}

# NOTE: The type of the pair values cannot be restricted, it must be `<:Any`,
#       not even `V where {V<:Any}`.
const VectorOfCardPairs{K<:CardName} = AbstractVector{<:Pair{K,<:Any}}

# Aliases used for sub-indexing.
const IndexRange = OrdinalRange{<:Integer,<:Integer}
const SubArrayIndex = Union{Colon,Integer,IndexRange}
const SubArrayIndices{N} = NTuple{N,SubArrayIndex}

const ColumnName = Union{AbstractString,Symbol}

"""
    EasyFITS.Header

is the (union of) type(s) that are accepted to specify a FITS header in
`EasyFITS` package.

A header may be a vector of pairs like `key => val`, `key => (val,com)`, or
`key => com` with `key` the keyword name, `val` its value, and `com` its
comment. The keyword name `key` is a string or a symbol which is automatically
converted to uppercase letters and trailing spaces discarded. The syntax `key
=> com`, with `com` a string, is only allowed for commentary keywords `COMMENT`
or `HISTORY`. For other keywords, the value is mandatory but the comment is
optional, not specifying the comment is like specifying `nothing` for the
comment; otherwise, the comment must be a string. The value `val` may be
`missing` or `undef` to indicate that it is undefined. If the comment is too
long, it is automatically split across multiple records for commentary keywords
and it is truncated for other keywords. A non-commentary keyword may have units
specified in square brackets at the beginning of the associated comment.
Commentary keywords may appear more than once, other keywords are unique.

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

defines a possible FITS header with several records: a keyword `VERIFIED`
having a logical value and no comments, a keyword `COUNT` having an integer
value and a comment, a keyword `SPEED` having a floating-point value and a
comment with units, a keyword `USER` having a string value, keywords `AGE` and
`JOB` having comments but undefined values, and a few additional commentary
keywords.

A header may also be specified as a named tuple with entries `key = val`, `key
= (val,com)`, or `key = com`. The same rules apply as above except that `key`
must be allowed as a variable symbolic name (no embedded hyphen `'-'`).

Finally, most methods assume that `nothing` can be used to indicate an empty
header.

!!! note
    Specifying a FITS header as a dictionary is purposely not implemented
    because, to a certain extend, the order of keywords in a FITS header is
    relevant and because some keywords (`COMMENT`, `HISTORY`, and `CONTINUE`)
    may appear more than once.

"""
const Header = Union{FitsHeader,NamedTuple,VectorOfCardPairs}

"""
    EasyFITS.ImageData{T,N}

"""
const ImageData{T,N} = AbstractArray{T,N}

"""
    EasyFITS.TableData

is the union of types of arguments that can be specified to represent data in a
FITS table extension.

A loop through all columns of `dat::TableData` writes:

    for k in keys(dat)
        name, vals = EasyFITS.get_column(dat, k)
    end

where non-exported method `get_column(dat,k)` yields a pair `name => vals` of
the name and the values of the `k`-th column in dat.

"""
const TableData = Union{AbstractVector{<:Pair{<:ColumnName,<:AbstractArray}},
                        Tuple{Vararg{<:Pair{<:ColumnName,<:AbstractArray}}},
                        NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}},
                        AbstractDict{<:ColumnName,<:AbstractArray}}

ncols(dat::TableData) = length(dat)

function get_column(dat::AbstractVector{<:Pair{<:ColumnName,<:AbstractArray}}, k::Int)
    key, vals = dat[k]
    String(key) => vals
end
function get_column(dat::Tuple{Vararg{<:Pair{<:ColumnName,<:AbstractArray}}}, k::Int)
    key, vals = dat[k]
    String(key) => vals
end
function get_column(dat::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}, key::Symbol)
    vals = dat[key]
    String(key) => vals
end
function get_column(dat::AbstractDict{K,<:AbstractArray}, key::K) where {K<:ColumnName}
    vals = dat[key]
    String(key) => vals
end

"""
    FitsHDU

is the abstract type of FITS Header Data Units which consists in a header and a
data parts. Concrete instances of `FitsHDU` behave as vectors whose elements
are FITS header records, a.k.a. FITS cards, and which can be indexed by
integers or by names.

For faster access to the records of a header, consider creating a FITS header
object from a HDU object:

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

"""
    FitsFile(filename, mode="r"; extended=false) -> file

opens FITS file `filename` for reading if `mode` is `"r"`, for reading and
writing if mode is "r+", or creates a new file if mode is `"w"` or `"w!"`. File
must not exists if mode is `"w"`. File is overwritten if it exists and mode is
`"w!"`. The file is automatically closed when the `file` object is finalized so
it is not necessary to call `close(file)`.

Keyword `extended` specifies whether to use extended file name syntax featured
by the CFITSIO library.

"""
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
        status = Ref{Status}(0)
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

is the singleton type of the object used to indicate invalid arguments while
sparing throwing an exception.

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
        isbitstype(T) || bad_argument("parameter T=$T is not a plain type")
        isa(N, Int) && N ≥ 0 || bad_argument("parameter N=$N must be a nonnegative `Int`")
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

yields a singleton object used to indicate that FITS rules should be applied for some
logical operation.  For example:

    isequal(FitsLogic(), s1, s2)

compares strings `s1` and `s2` according to FITS rules, that is case of letters
and trailing spaces are irrelevant.

    isequal(FitsLogic(), x) -> f

yields a predicate function `f` such that `f(y)` yields
`isequal(FitsLogic(),x,y)`.

"""
struct FitsLogic end

struct Bit end

struct FitsError <: Exception
    code::Status
end
