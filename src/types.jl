# According to FITS standard, the comment separator (a space followed by a '/'
# character) is strongly recommended but not mandatory.
const MANDATORY_COMMENT_SEPARATOR = false

# Encoding for FITS cards should be ASCII but UTF8 is also supported.
const CARD_ENCODING = :utf8

const Status = Cint
const OptionalString = Union{AbstractString,Nothing}

# Types indicating an undefined card value, Nothing is used for the value of
# commentary cards.
const UndefinedValue = Union{Missing,UndefInitializer}

# FITS header card name, a.k.a. FITS keyword.
const CardName = Union{Symbol,AbstractString}

# Acceptable types for the value of a FITS header card.
const CardValue = Union{Nothing,UndefinedValue,Real,Complex,AbstractString}

# The different ways to represent a FITS header card data.
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
     "UNDEF" => (undef, "Some undefined value."),
     "MISSING" => (missing, "Another undefined value."),
     "HISTORY" => "Some historical information.",
     "COMMENT" => "Some comment.",
     "COMMENT" => "Some other comment.",
     "HISTORY" => "Some other historical information."]

defines a possible FITS header with several records: a keyword `VERIFIED`
having a logical value and no comments, a keyword `COUNT` having an integer
value and a comment, a keyword `SPEED` having a floating-point value and a
comment with units, a keyword `USER` having a string value, keywords `UNDEF`
and `MISSING` having comments but undefined values, and a few additional
commentary keywords.

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
const Header = Union{NamedTuple,VectorOfCardPairs}

# String decoration to represent a FITS file name.
struct FitsFile <: AbstractString
    path::String
end

abstract type FitsHDU end

# Enumeration of HDU type identifiers.
@enum FitsHDUType::Cint begin
    FITS_IMAGE_HDU = CFITSIO.IMAGE_HDU
    FITS_BINARY_TABLE_HDU = CFITSIO.BINARY_TBL
    FITS_ASCII_TABLE_HDU = CFITSIO.ASCII_TBL
    FITS_ANY_HDU = CFITSIO.ANY_HDU
end

# Enumeration of keyword value type identifiers.
@enum FitsCardType::Cint begin
    FITS_UNKNOWN   = -1
    FITS_UNDEFINED =  0 # no value given
    FITS_LOGICAL   =  1
    FITS_INTEGER   =  2
    FITS_FLOAT     =  3
    FITS_COMPLEX   =  4
    FITS_STRING    =  5
    FITS_COMMENT   =  6
end

"""
    FitsIO(path, mode="r"; extended=false) -> io

opens FITS file `path` for reading if `mode` is `"r"`, for reading and writing
if mode is "r+", or creates a new file if mode is `"w"` or `"w!"`. File must
not exists if mode is `"w"`. File is overwritten if it exists and mode is
`"w!"`. The file is automatically closed when the `io` object is finalized so
it is not necessary to call `close(io)`.

Keyword `extended` specifies whether to use extended file name syntax featured
by CFITSIO library.

"""
mutable struct FitsIO <: AbstractVector{FitsHDU}
    handle::Ptr{CFITSIO.fitsfile}
    mode::Symbol
    path::String
    nhdus::Int
    function FitsIO(path::AbstractString, mode::AbstractString = "r"; extended::Bool=false)
        access = mode == "r" ? :r :
            mode == "r+" ? :rw :
            mode == "w" || mode == "w!" ? :w :
            throw(ArgumentError("mode must be \"r\", \"r+\", \"w\" or \"w!\""))
        status = Ref{Status}(0)
        handle = Ref{Ptr{CFITSIO.fitsfile}}(0)
        num = Ref{Cint}(0)
        if access === :w
            if extended
                prefix = (mode == "w!" && ! startswith(path, '!') ? "!" : "")
                check(CFITSIO.fits_create_file(handle, prefix*path, status))
            else
                if isfile(path)
                    mode == "w!" || throw_file_already_exists(path, "use mode \"w!\" to overwrite")
                    rm(path; force=true)
                end
                check(CFITSIO.fits_create_diskfile(handle, path, status))
            end
        else
            iomode = (access === :rw ? CFITSIO.READWRITE : CFITSIO.READONLY)
            if extended
                check(CFITSIO.fits_open_file(handle, path, iomode, status))
            else
                check(CFITSIO.fits_open_diskfile(handle, path, iomode, status))
            end
            if !iszero(CFITSIO.fits_get_num_hdus(handle[], num, status))
                errcode = status[]
                status[] = 0
                CFITSIO.fits_close_file(ptr, status)
                throw(FitsError(errcode))
            end
        end
        return finalizer(_close, new(handle[], access, path, num[]))
    end
end

"""
    EasyFits.Invalid

is a singleton used to indicate invalid arguments while sparing throwing an
exception.

"""
struct Invalid end

# Singleton type to indicate that the inner constructor should be called.
struct BareBuild end

struct FitsAnyHDU <: FitsHDU
    io::FitsIO
    num::Int
    FitsAnyHDU(::BareBuild, io::FitsIO, num::Integer) = new(io, num)
end

# FITS Image (Array for Julia).
struct FitsImageHDU{T,N} <: FitsHDU
    io::FitsIO
    num::Int
    function FitsImageHDU{T,N}(::BareBuild, io::FitsIO, num::Integer) where {T,N}
        isbitstype(T) || bad_argument("parameter T=$T is not a plain type")
        isa(N, Int) && N ≥ 0 || bad_argument("parameter N=$N is a nonnegative Int")
        return new{T,N}(io, num)
    end
end

# FITS table.
struct FitsTableHDU <: FitsHDU
    io::FitsIO
    num::Int
    ascii::Bool
    FitsTableHDU(::BareBuild, io::FitsIO, num::Integer, ascii::Bool) =
        new(io, num, ascii)
end

# Object to store a single FITS header card. Must be mutable to have access to
# its address.
mutable struct FitsCard <: AbstractVector{UInt8}
    data::NTuple{CFITSIO.FLEN_CARD,UInt8}
    type::FitsCardType
    name_offset::Int    # offset to the name part
    name_length::Int    # length of the name part
    value_offset::Int   # offset to the value part
    value_length::Int   # length of the value part
    comment_offset::Int # offset to the comment part
    comment_length::Int # length of the comment part
    function FitsCard()
        card = new()
        card[firstindex(card)] = zero(eltype(card))
        set_type!(card, FITS_UNKNOWN)
        set_name_offset!(card, 0)
        set_name_length!(card, 0)
        set_value_offset!(card, 0)
        set_value_length!(card, 0)
        set_comment_offset!(card, 0)
        set_comment_length!(card, 0)
        return card
    end
end

# Any FITS card part (name, value, or comment) is behaving as an abstract
# string whose encoding is specified by the type parameter E.
abstract type FitsCardPart{E} <: AbstractString end

# The 3 possible parts of a FITS card are simple decorations behaving as
# abstract strings.
struct FitsCardName    <: FitsCardPart{CARD_ENCODING}; parent::FitsCard; end
struct FitsCardValue   <: FitsCardPart{CARD_ENCODING}; parent::FitsCard; end
struct FitsCardComment <: FitsCardPart{CARD_ENCODING}; parent::FitsCard; end

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
