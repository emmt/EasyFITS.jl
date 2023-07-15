# Management of FITS header data units (HDUs) and of FITS header cards.

"""
    FitsHDU(file::FitsFile, i) -> hdu
    getindex(file::FitsFile, i) -> hdu
    file[i] -> hdu

yield the `i`-th FITS header data unit of FITS file `file`.

The returned object has the following read-only properties:

    hdu.file     # associated FITS file
    hdu.number   # HDU number (the index `i` above)
    hdu.type     # same as FitsHDUType(hdu)
    hdu.xtension # value of the XTENSION card (never nothing)
    hdu.extname  # value of the EXTNAME card or nothing
    hdu.hduname  # value of the HDUNAME card or nothing

"""
function FitsHDU(file::FitsFile, i::Integer)
    status = Ref{Status}(0)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(file, i, type, status))
    type = FitsHDUType(type[])
    if type == FITS_ASCII_TABLE_HDU
        return FitsTableHDU(BareBuild(), file, i, true)
    elseif type == FITS_BINARY_TABLE_HDU
        return FitsTableHDU(BareBuild(), file, i, false)
    elseif type == FITS_IMAGE_HDU
        bitpix = Ref{Cint}()
        check(CFITSIO.fits_get_img_equivtype(file, bitpix, status))
        ndims = Ref{Cint}()
        check(CFITSIO.fits_get_img_dim(file, ndims, status))
        N = Int(ndims[])::Int
        T = type_from_bitpix(bitpix[])
        return FitsImageHDU{T,N}(BareBuild(), file, i)
    else
        return FitsAnyHDU(BareBuild(), file, i)
    end
end

FitsHDU(file::FitsFile, s::AbstractString) = file[s]

FitsTableHDU(file::FitsFile, i::Union{Integer,AbstractString}) =
    FitsHDU(file, i)::FitsTableHDU

FitsImageHDU(file::FitsFile, i::Union{Integer,AbstractString}) =
    FitsHDU(file, i)::FitsImageHDU
FitsImageHDU{T}(file::FitsFile, i::Union{Integer,AbstractString}) where {T} =
    FitsHDU(file, i)::FitsImageHDU{T}
FitsImageHDU{T,N}(file::FitsFile, i::Union{Integer,AbstractString}) where {T,N} =
    FitsHDU(file, i)::FitsImageHDU{T,N}

Base.propertynames(::FitsHDU) = (:extname, :hduname, :file, :number, :type, :xtension)
Base.getproperty(hdu::FitsHDU, sym::Symbol) = getproperty(hdu, Val(sym))

Base.getproperty(hdu::FitsHDU, ::Val{:extname})  = get_extname(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:hduname})  = get_hduname(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:file})     = get_file(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:number})   = get_number(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:type})     = FitsHDUType(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:xtension}) = get_xtension(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{sym}) where {sym} = invalid_property(hdu, sym)

Base.setproperty!(hdu::FitsHDU, sym::Symbol, x) =
    sym ∈ propertynames(hdu) ? readonly_property(hdu, sym) : invalid_property(hdu, sym)

"""
    hdu.hduname
    EasyFITS.hduname(hdu::FitsHDU)

yield the value of the keyword `HDUNAME` in the FITS Header Data Unit `hdu` or
`nothing` if no such keyword exists.

It is suggested that foreign packages extend this function so that:

    EasyFITS.hduname(T::Type) -> (name, vers)

yields the name and the version of the FITS Header Data Unit when saving an
object of type `T` in a FITS file. This is also useful to search for such
contents in a FITS file.

"""
hduname(hdu::FitsHDU) = get_hduname(hdu)

"""
    EasyFITS.get_file_at(hdu::FitsHDU) -> file

yields the FITS file associated with FITS Header Data Unit `hdu` checking that
the file is still open and moving the file position of `file` to that of `hdu`.

"""
function get_file_at(hdu::FitsHDU)
    file = get_file(hdu)
    check(CFITSIO.fits_movabs_hdu(file, get_number(hdu),
                                  Ptr{Cint}(0), Ref{Status}(0)))
    return file
end

# Yields the object to preserve and to use to retrieve the pointer to the
# opaque C structure.
Base.cconvert(::Type{Ptr{CFITSIO.fitsfile}}, hdu::FitsHDU) = get_file_at(hdu)

get_file(hdu::FitsHDU) = getfield(hdu, :file)
get_number(hdu::FitsHDU) = getfield(hdu, :num)

get_xtension(hdu::FitsHDU) = "ANY"
get_xtension(hdu::FitsImageHDU) = "IMAGE"
get_xtension(hdu::FitsTableHDU) = isascii(hdu) ? "TABLE" : "BINTABLE"

for (func, key) in ((:get_extname, "EXTNAME"),
                    (:get_hduname, "HDUNAME"))
    @eval function $func(hdu::FitsHDU)
        card = get(hdu, $key, nothing)
        card === nothing && return nothing
        card.type == FITS_STRING || return nothing
        return card.string
    end
end

"""
    FitsHDUType(hdu)
    hdu.type

yield the type code of the FITS header data unit `hdu`, one of:

    FITS_ANY_HDU           # HDU type unknown or not yet defined
    FITS_IMAGE_HDU         # FITS Image, i.e. multi-dimensional array
    FITS_ASCII_TABLE_HDU   # FITS ASCII Table
    FITS_BINARY_TABLE_HDU  # FITS Binary Table

"""
FitsHDUType(hdu::FitsHDU)      = FITS_ANY_HDU
FitsHDUType(hdu::FitsImageHDU) = FITS_IMAGE_HDU
FitsHDUType(hdu::FitsTableHDU) = isascii(hdu) ? FITS_ASCII_TABLE_HDU : FITS_BINARY_TABLE_HDU

"""
    isascii(hdu::FitsHDU) -> bool

yields whether FITS HDU `hdu` is an ASCII table.

"""
Base.isascii(hdu::FitsHDU) = false
Base.isascii(hdu::FitsTableHDU) = getfield(hdu, :ascii)

"""
    getindex(hdu::FitsHDU, key) -> card::FitsCard
    hdu[key] -> card

yield the FITS header card at index `key` (can be an integer or a string) in
header data unit `hdu`. The result is a `FitsCard` object with the following
properties:

    card.type                  # type of card: FITS_LOGICAL, FITS_INTEGER, etc.
    card.key                   # quick key of card: Fits"BITPIX", Fits"HIERARCH", etc.
    card.name                  # name of card
    card.value                 # callable object representing the card value
    card.comment               # comment of card
    card.units                 # units of card value
    card.unitless              # comment of card without the units part if any
    card.logical :: Bool       # alias for card.value(Bool)
    card.integer :: $FitsInteger      # alias for card.value(Integer)
    card.float   :: $FitsFloat    # alias for card.value(Real)
    card.complex :: $FitsComplex # alias for card.value(Complex)
    card.string  :: String     # alias for card.value(String)

Examples:

    ndims = hdu["NAXIS"].value(Int)
    dims = ntuple(i -> hdu["NAXIS\$i"].value(Int), Val(ndims))

"""
function Base.getindex(hdu::FitsHDU, key::Union{CardName,Integer})
    buf = SmallVector{CFITSIO.FLEN_CARD,UInt8}(undef)
    status = @inbounds try_read!(buf, hdu, key)
    if iszero(status)
        return parse_cstring(FitsCard, buf)
    elseif status == (key isa Integer ? CFITSIO.KEY_OUT_BOUNDS : CFITSIO.KEY_NO_EXIST)
        key isa Integer ? throw(BoundsError(hdu, key)) : throw(KeyError(key))
    else
        throw(FitsError(status))
    end
end

function Base.get(hdu::FitsHDU, key::Union{CardName,Integer}, def)
    buf = SmallVector{CFITSIO.FLEN_CARD,UInt8}(undef)
    status = @inbounds try_read!(buf, hdu, key)
    if iszero(status)
        return parse_cstring(FitsCard, buf)
    elseif status == (key isa Integer ? CFITSIO.KEY_OUT_BOUNDS : CFITSIO.KEY_NO_EXIST)
        return def
    else
        throw(FitsError(status))
    end
end

Base.haskey(hdu::FitsHDU, key::Integer) = key ∈ keys(hdu)
function Base.haskey(hdu::FitsHDU, key::CardName)
    buf = SmallVector{CFITSIO.FLEN_CARD,UInt8}(undef)
    status = @inbounds try_read!(buf, hdu, key)
    iszero(status) && return true
    status == CFITSIO.KEY_NO_EXIST && return false
    throw(FitsError(status))
end

@inline function try_read!(buf::AbstractVector{UInt8}, hdu::FitsHDU,
                           key::Union{CardName,Integer})
    length(buf) ≥ CFITSIO.FLEN_CARD || error("buffer is too small")
    if !(key isa Integer)
        return CFITSIO.fits_read_card(hdu, key, buf, Ref{Status}(0))
    elseif key ≥ firstindex(hdu)
        return CFITSIO.fits_read_record(hdu, key, buf, Ref{Status}(0))
    else
        return CFITSIO.KEY_OUT_BOUNDS
    end
end

# This function is needed to truncate C-string at 1st null, we take the
# opportunity of this filtering to strip trailing spaces.
# FIXME: This maybe done elsewhere?
function parse_cstring(::Type{FitsCard}, buf::AbstractVector{UInt8})
    first = firstindex(buf)
    last = first - 1
    @inbounds for i ∈ eachindex(buf)
        b = buf[i]
        if b != 0x20 # not a space
            if b == 0x00
                break
            end
            last = i
        end
    end
    return FitsCard(@inbounds view(buf, first:last))
end

"""
    FitsHeader(hdu::FitsHDU)

reads all records of the header of `hdu`.

"""
function BaseFITS.FitsHeader(hdu::FitsHDU)
    file = get_file_at(hdu)
    len = length(hdu)
    hdr = sizehint!(FitsHeader(), len)
    buf = SmallVector{CFITSIO.FLEN_CARD,UInt8}(undef)
    status = Ref{Status}(0)
    @inbounds for i in 1:len
        check(CFITSIO.fits_read_record(file, i, buf, status))
        push!(hdr, parse_cstring(FitsCard, buf))
    end
    return hdr
end

"""
    reset(hdu::FitsHDU) -> hdu

reset the search by wild card characters in FITS header of `hdu` to the first
record.

"""
function Base.reset(hdu::FitsHDU)
    check(CFITSIO.fits_read_record(hdu, 0, Ptr{Cchar}(0), Ref{Status}(0)))
    return hdu
end

"""
    setindex!(hdu::FitsHDU, x, key) -> hdu
    hdu[key] = x
    push!(hdu, key => x)

updates or appends a record associating the keyword `key` with `x` in the
header of the FITS header data unit `hdu` (see `push!`). Depending on the type
of the FITS keyword `key`, the argument `x` can be `val`, `com`, or `(val,com)`
with `val` and `com` the value and the comment of the FITS card.

"""
Base.setindex!(hdu::FitsHDU, x, key::CardName) = push!(hdu, key => x)

"""
    push!(hdu::FitsHDU, rec; append=false) -> hdu

updates or appends header record `rec` to FITS Header Data Units `hdu`. If the
name of `rec` does not yet exist in the header part of `hdu` or if it is a
commentary or continuation FITS keyword (`"COMMENT"`, `"HISTORY"`, `""`, or
`"CONTINUE"`), a new record is appended to the header part of `hdu`; otherwise,
the existing record in `hdu` is updated. If keyword `append` is set true, the
record is appended whether another record with the same name already exists or
not. Forcing append is not recommended as it may result in an invalid header.

Argument `rec` may be a FITS card (of type `FitsCard`) or anything that can be
converted into a FITS card. This includes a pair `key => val`, `key => com`, or
`key => (val,com)` with `key` the keyword of the record, `val` its value and
`com` its comment.

To push more than one record, call `merge!` instead of `push!`.

"""
function push!(hdu::FitsHDU, card::FitsCard; append::Bool = false)
    # Private method.
    function set_key(hdu::FitsHDU, key::String, val, com::String; append::Bool = false)
        if append
            write_key(hdu, key, val, com)
        else
            update_key(hdu, key, val, com)
        end
        nothing
    end

    if card.type === FITS_LOGICAL
        set_key(hdu, card.name, card.logical, card.comment; append)
    elseif card.type === FITS_INTEGER
        set_key(hdu, card.name, card.integer, card.comment; append)
    elseif card.type === FITS_FLOAT
        set_key(hdu, card.name, card.float, card.comment; append)
    elseif card.type === FITS_STRING
        set_key(hdu, card.name, card.string, card.comment; append)
    elseif card.type === FITS_COMPLEX
        set_key(hdu, card.name, card.complex, card.comment; append)
    elseif card.type === FITS_COMMENT
        if card.key == Fits"COMMENT"
            write_comment(hdu, card.comment)
        elseif card.key == Fits"HISTORY"
            write_history(hdu, card.comment)
        elseif card.key == Fits"" || card.key == Fits"CONTINUE"
            write_key(hdu, card.name, nothing, card.comment)
        else
            update_key(hdu, card.name, nothing, card.comment)
        end
    elseif card.type === FITS_UNDEFINED
        set_key(hdu, card.name, missing, card.comment; append)
    elseif card.type !== FITS_END
        error("unexpected FITS card type")
    end
    return hdu
end

push!(hdu::FitsHDU, rec; kwds...) = push!(hdu, FitsCard(rec); kwds...)

"""
    merge!(hdu::FitsHDU, recs) -> hdu

pushes all FITS header cards in `recs` into FITS Header Data Units `hdu` and
returns it. Examples:

    merge!(hdu::FitsHDU, ["key1" => dat1, "key2" => dat2, ...]) -> hdu
    merge!(hdu::FitsHDU, (key1 = dat1, key2 = dat2, ...)) -> hdu

In most cases, calling `merge!` is a shortcut to:

    for rec in recs
        push!(hdu, rec)
    end

"""
function merge!(hdu::FitsHDU, hdr::Header; append::Bool = false)
    # By default, assume an iterable object producing cards or equivalent.
    for rec in hdr
        push!(hdu, rec; append)
    end
    return hdu
end

merge!(hdu::FitsHDU, recs::Nothing; append::Bool = false) = hdu

function merge!(hdu::FitsHDU, hdr::NamedTuple; append::Bool = false)
    for key in keys(hdr)
        push!(hdu, key => hdr[key]; append)
    end
    return hdu
end

"""
    delete!(hdu::FitsHDU, key) -> hdu

deletes from FITS header data unit `hdu` the header card identified by the name
or number `key`.

"""
function Base.delete!(hdu::FitsHDU, key::CardName)
    check(CFITSIO.fits_delete_key(hdu, key, Ref{Status}(0)))
    return hdu
end

function Base.delete!(hdu::FitsHDU, key::Integer)
    check(CFITSIO.fits_delete_record(hdu, key, Ref{Status}(0)))
    return hdu
end

"""
    EasyFITS.write_key(dst, key, val, com=nothing) -> dst

appends a new FITS header card in `dst` associating value `val` and comment
`com` to the keyword `key`.

""" write_key

"""
    EasyFITS.update_key(dst, key, val, com=nothing) -> dst

updates or appends a FITS header card in `dst` associating value `val` and
comment `com` to the keyword `key`.

The card is considered to have an undefined value if `val` is `missing` or
`undef`.

If `val` is `nothing` and `com` is a string, the comment of the FITS header
card is updated to be `com`.

"""
update_key(hdu::FitsHDU, key::CardName, val::Nothing, com::Nothing) = hdu
function update_key(hdu::FitsHDU, key::CardName, val::Nothing, com::AbstractString)
    # BUG: When modifying the comment of an existing keyword which has an
    #      undefined value, the keyword becomes a commentary keyword.
    check(CFITSIO.fits_modify_comment(hdu, key, com, Ref{Status}(0)))
    return hdu
end

unsafe_optional_string(s::AbstractString) = s
unsafe_optional_string(s::Nothing) = Ptr{Cchar}(0)

for func in (:update_key, :write_key),
    (V, T) in ((AbstractString, String),
               (Bool,           Bool),
               (Integer,        Int),
               (Real,           Cdouble),
               (Complex,        Complex{Cdouble}),
               (Undefined,      Missing))
    if T === Missing
        # Commentary (nothing) card or undefined value (undef or missing).
        @eval function $func(dst, key::CardName, val::$V, com::CardComment=nothing)
            check(CFITSIO.$(Symbol("fits_",func,"_null"))(
                dst, key, unsafe_optional_string(com), Ref{Status}(0)))
            return dst
        end
    elseif T === String
        @eval function $func(dst, key::CardName, val::$V, com::CardComment=nothing)
            check(CFITSIO.$(Symbol("fits_",func,"_str"))(
                dst, key, val, unsafe_optional_string(com), Ref{Status}(0)))
            return dst
        end
    elseif T === Bool
        @eval function $func(dst, key::CardName, val::$V, com::CardComment=nothing)
            check(CFITSIO.$(Symbol("fits_",func,"_log"))(
                dst, key, val, unsafe_optional_string(com), Ref{Status}(0)))
            return dst
        end
    else # FIXME: use more specialized CFITSIO routines?
        @eval function $func(dst, key::CardName, val::$V, com::CardComment=nothing)
            check(CFITSIO.$(Symbol("fits_",func))(
                dst, type_to_code($T), key, Ref{$T}(val),
                unsafe_optional_string(com), Ref{Status}(0)))
            return dst
        end
    end
end

# Implement abstract vector API for HDUs.
Base.length(hdu::FitsHDU) = get_hdrspace(hdu)[1]
Base.size(hdu::FitsHDU) = (length(hdu),)
Base.axes(hdu::FitsHDU) = (keys(hdu),)
Base.firstindex(hdu::FitsHDU) = 1
Base.lastindex(hdu::FitsHDU) = length(hdu)
Base.keys(hdu::FitsHDU) = Base.OneTo(length(hdu))

# Yield number of existing and remaining undefined keys in the current HDU.
function get_hdrspace(hdu::FitsHDU)
    existing = Ref{Cint}()
    remaining = Ref{Cint}()
    check(CFITSIO.fits_get_hdrspace(hdu, existing, remaining, Ref{Status}(0)))
    return (Int(existing[]), Int(remaining[]))
end

"""
    EasyFITS.write_comment(dst, str) -> dst

appends a new FITS comment record in `dst` with text string `str`. The comment
record will be continued over multiple cards if `str` is longer than 70
characters.

""" write_comment

"""
    EasyFITS.write_history(dst, str) -> dst

appends a new FITS history record in `dst` with text string `str`. The history
record will be continued over multiple cards if `str` is longer than 70
characters.

""" write_history

for func in (:write_comment, :write_history)
    @eval begin
        function $func(dst, str::AbstractString)
            check(CFITSIO.$(Symbol("fits_",func))(dst, str, Ref{Status}(0)))
            return dst
        end
    end
end

"""
    EasyFITS.write_date(hdu::FitsHDU) -> hdu

creates or updates a FITS comment record of `hdu` with the current date.

"""
function write_date(hdu::FitsHDU)
    check(CFITSIO.fits_write_date(hdu, Ref{Status}(0)))
    return hdu
end

function write_comment(hdu::FitsHDU, str::AbstractString)
    check(CFITSIO.fits_write_comment(hdu, str, Ref{Status}(0)))
    return hdu
end

function write_history(hdu::FitsHDU, str::AbstractString)
    check(CFITSIO.fits_write_history(hdu, str, Ref{Status}(0)))
    return hdu
end

function write_record(f::Union{FitsFile,FitsHDU}, card::FitsCard)
    check(CFITSIO.fits_write_record(f, card, Ref{Status}(0)))
    return f
end

function update_record(f::Union{FitsFile,FitsHDU}, key::CardName, card::FitsCard)
    check(CFITSIO.fits_update_card(f, key, card, Ref{Status}(0)))
    return f
end

Base.show(io::IO, hdu::FitsImageHDU{T,N}) where {T,N} = begin
    print(io, "FitsImageHDU{")
    print(io, T)
    print(io, ',')
    print(io, N)
    print(io, '}')
end
Base.show(io::IO, hdu::FitsTableHDU) = print(io, "FitsTableHDU")
Base.show(io::IO, hdu::FitsAnyHDU) = print(io, "FitsAnyHDU")

function Base.show(io::IO, mime::MIME"text/plain", hdu::FitsHDU)
    show(io, hdu)
    if hdu isa FitsTableHDU
        print(io, isascii(hdu) ? " (ASCII)" : " (Binary)")
    end
    print(io, ": num = ")
    print(io, hdu.number)
    if hdu isa FitsImageHDU || hdu isa FitsTableHDU
        for name in ("HDUNAME", "EXTNAME")
            card = get(hdu, name, nothing)
            if card !== nothing
                print(io, ", ")
                print(io, name)
                print(io, " = ")
                show(io, mime, card.value)
            end
        end
    end
    if hdu isa FitsImageHDU
        naxis = hdu["NAXIS"].integer
        if naxis > 0
            print(io, ", size = ")
            for i in 1:naxis
                i > 1 && print(io, '×')
                print(io, hdu["NAXIS$i"].integer)
            end
        end
    end
end
