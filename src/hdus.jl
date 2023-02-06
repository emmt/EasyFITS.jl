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

@inline function try_read!(buf::AbstractVector{UInt8}, hdu::FitsHDU, key::CardName)
    @boundscheck length(buf) ≥ CFITSIO.FLEN_CARD || error("buffer is too small")
    return CFITSIO.fits_read_card(hdu, key, buf, Ref{Status}(0))
end

@inline function try_read!(buf::AbstractVector{UInt8}, hdu::FitsHDU, key::Integer)
    @boundscheck length(buf) ≥ CFITSIO.FLEN_CARD || error("buffer is too small")
    return key ≤ 0 ? CFITSIO.KEY_OUT_BOUNDS :
        CFITSIO.fits_read_record(hdu, key, buf, Ref{Status}(0))
end

# This function is needed to truncate C-string at 1st null, we take the
# opportunity of thsi filtering to strip trailing spaces.
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
function FITSBase.FitsHeader(hdu::FitsHDU)
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
    setindex!(hdu::FitsHDU, dat, key) -> hdu
    hdu[key] = dat

updates or appends a record associating the keyword `key` with the data `dat`
in the header of the FITS header data unit `hdu`. See [`EasyFITS.Header`](@ref)
for the possible forms of `dat`.

"""
Base.setindex!(hdu::FitsHDU, dat::CardData, key::CardName) =
    set_key(hdu, key => dat; update=true)

"""
    push!(hdu::FitsHDU, key => dat) -> hdu
    hdu[key] = dat

appends a new record associating the keyword `key` with the data `dat` in the
header of the FITS header data unit `hdu`. See [`EasyFITS.Header`](@ref) for
the possible forms of such pairs.

A vector of pairs or a named tuple may be specified to push more than one
record in a single call:

    push!(hdu::FitsHDU, ["key1" => dat1, "key2" => dat2, ...]) -> hdu
    push!(hdu::FitsHDU, (key1 = dat1, key2 = dat2, ...)) -> hdu

"""
Base.push!(hdu::FitsHDU, pair::Pair{<:CardName}) =
    set_key(hdu, pair; update=false)

Base.push!(hdu::FitsHDU, ::Nothing) = hdu

function Base.push!(hdu::FitsHDU, cards::VectorOfCardPairs)
    for card in cards
        push!(hdu, card)
    end
    return hdu
end

function Base.push!(hdu::FitsHDU, cards::NamedTuple)
    for key in keys(cards)
        push!(hdu, key => cards[key])
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
   EasyFITS.is_comment_keyword(key)

yields whether `key` is the name of a commentary FITS keyword.

"""
function is_comment_keyword(key::AbstractString)
    len = length(key)
    len == 0 && return true
    c = uppercase(first(key))
    c == 'C' && return isequal(FitsLogic(), key, "COMMENT")
    c == 'H' && return isequal(FitsLogic(), key, "HISTORY")
    c == ' ' && return isequal(FitsLogic(), key, "")
    return false
end
is_comment_keyword(key::Symbol) =
    key === :COMMENT || key === :HISTORY ||
    key === :comment || key === :history

@inline function set_key(hdu::FitsHDU, pair::CardPair; update::Bool)
    key = first(pair)
    dat = last(pair)
    type = keyword_type(key)
    if type === :COMMENT || type === :HISTORY
        com = dat isa AbstractString ? dat :
            dat isa Nothing ? "" : invalid_card_value(key, dat, true)
        if type === :HISTORY
            write_history(hdu, com)
        else
            write_comment(hdu, com)
        end
    else
        val, com = normalize_card_data(dat)
        val === Invalid() && invalid_card_value(key, dat, false)
        if update
            update_key(hdu, key, val, com)
        else
            write_key(hdu, key, val, com)
        end
    end
    return hdu
end

"""
    @eq b c
    @eq c b

yields whether expression/symbol `b` resulting in a byte value (of type
`UInt8`) is equal to char `c`. The comparison is case insensitive and assumes
ASCII characters.

"""
macro eq(b::Union{Expr,Symbol}, c::Char)
    esc(_eq(b, c))
end

macro eq(c::Char, b::Union{Expr,Symbol})
    esc(_eq(b, c))
end

_eq(b::UInt8, c::UInt8) = b == c
_eq(b::UInt8, c1::UInt8, c2::UInt8) = ((b == c1)|(b == c2))
_eq(b::Union{Expr,Symbol}, c::Char) =
    'a' ≤ c ≤ 'z' ? :(_eq($b, $(UInt8(c) & ~0x20), $(UInt8(c)))) :
    'A' ≤ c ≤ 'Z' ? :(_eq($b, $(UInt8(c)), $(UInt8(c) | 0x20))) :
    :(_eq($b, $(UInt8(c))))

"""
    EasyFITS.keyword_type(key)

yields the type of FITS keyword `key` (a symbol or a string) as one of:
`:COMMENT`, `:HISTORY`, `:CONTINUE`, or `:HIERARCH`. The comparison is
insensitive to case of letters and to trailing spaces. The method is meant to
be fast.

"""
function keyword_type(str::Union{String,Symbol})
    # NOTE: This single-pass version takes at most 7ns for all keywords on my
    # laptop dating from 2015.
    #
    # NOTE: `Base.unsafe_convert(Ptr{UInt8},str::String)` and
    # `Base.unsafe_convert(Ptr{UInt8},sym::Symbol)` both yield a pointer to a
    # null terminated array of bytes, thus suitable for direct byte-by-byte
    # comparison between ASCII C strings. A `String` may have embedded nulls
    # (not a `Symbol`), however these nulls will just stop the comparison which
    # is what we want. This save us from checking for nulls as would be the
    # case if we used `Cstring` instead of `Ptr{UInt8}`.
    obj = Base.cconvert(Ptr{UInt8}, str)
    ptr = Base.unsafe_convert(Ptr{UInt8}, obj)
    GC.@preserve obj begin
        b1 = unsafe_load(ptr, 1)
        if @eq('C', b1)
            # May be COMMENT or CONTINUE.
            if @eq('O', unsafe_load(ptr, 2))
                b3 = unsafe_load(ptr, 3)
                if @eq('M', b3) &&
                    @eq('M', unsafe_load(ptr, 4)) &&
                    @eq('E', unsafe_load(ptr, 5)) &&
                    @eq('N', unsafe_load(ptr, 6)) &&
                    @eq('T', unsafe_load(ptr, 7))
                    unsafe_only_spaces(ptr, 8) && return :COMMENT
                elseif @eq('N', b3) &&
                    @eq('T', unsafe_load(ptr, 4)) &&
                    @eq('I', unsafe_load(ptr, 5)) &&
                    @eq('N', unsafe_load(ptr, 6)) &&
                    @eq('U', unsafe_load(ptr, 7)) &&
                    @eq('E', unsafe_load(ptr, 8))
                    unsafe_only_spaces(ptr, 9) && return :CONTINUE
                end
            end
        elseif @eq('H', b1)
            # May be HISTORY or HIERARCH.
            if @eq('I', unsafe_load(ptr, 2))
                b3 = unsafe_load(ptr, 3)
                if @eq('S', b3) &&
                    @eq('T', unsafe_load(ptr, 4)) &&
                    @eq('O', unsafe_load(ptr, 5)) &&
                    @eq('R', unsafe_load(ptr, 6)) &&
                    @eq('Y', unsafe_load(ptr, 7))
                    unsafe_only_spaces(ptr, 8) && return :HISTORY
                elseif @eq('E', b3) &&
                    @eq('R', unsafe_load(ptr, 4)) &&
                    @eq('A', unsafe_load(ptr, 5)) &&
                    @eq('R', unsafe_load(ptr, 6)) &&
                    @eq('C', unsafe_load(ptr, 7)) &&
                    @eq('H', unsafe_load(ptr, 8)) &&
                    @eq(' ', unsafe_load(ptr, 9))
                    unsafe_only_spaces(ptr, 10) || return :HIERARCH
                end
            end
        end
    end
    return :OTHER
end

# NOTE: This version takes less than 7ns for any keyword for all keywords
# on my laptop dating from 2015.
function keyword_type(str::AbstractString)
   n = ncodeunits(str)
    if n ≥ 7
        b1 = _load(str, 1)
        if @eq('C', b1)
            # Can be COMMENT or CONTINUE.
            if @eq('O', _load(str, 2))
                b3 = _load(str, 3)
                if @eq('M', b3) &&
                    @eq('M', _load(str, 4)) &&
                    @eq('E', _load(str, 5)) &&
                    @eq('N', _load(str, 6)) &&
                    @eq('T', _load(str, 7))
                    while n > 7 && @eq(' ', _load(str, n))
                        n -= 1
                    end
                    n == 7 && return :COMMENT
                elseif n ≥ 8 &&
                    @eq('N', b3) &&
                    @eq('T', _load(str, 4)) &&
                    @eq('I', _load(str, 5)) &&
                    @eq('N', _load(str, 6)) &&
                    @eq('U', _load(str, 7)) &&
                    @eq('E', _load(str, 8))
                    while n > 8 && @eq(' ', _load(str, n))
                        n -= 1
                    end
                    n == 8 && return :CONTINUE
                end
            end
        elseif @eq('H', b1)
            # Can be HISTORY or HIERARCH.
            if @eq('I', _load(str, 2))
                b3 = _load(str, 3)
                if @eq('S', b3) &&
                    @eq('T', _load(str, 4)) &&
                    @eq('O', _load(str, 5)) &&
                    @eq('R', _load(str, 6)) &&
                    @eq('Y', _load(str, 7))
                    while n > 7 && @eq(' ', _load(str, n))
                        n -= 1
                    end
                    n == 7 && return :HISTORY
                elseif n > 9 && @eq('E', b3) &&
                    @eq('R', _load(str, 4)) &&
                    @eq('A', _load(str, 5)) &&
                    @eq('R', _load(str, 6)) &&
                    @eq('C', _load(str, 7)) &&
                    @eq('H', _load(str, 8)) &&
                    @eq(' ', _load(str, 9))
                    for i in 10:n
                        @eq(' ', _load(str, i)) || return :HIERARCH
                    end
                end
            end
        end
    end
    return :OTHER
end

# Private method to mimic `unsafe_load` for other type of arguments.
_load(str::AbstractString, i::Int) = @inbounds codeunit(str, i)

@inline function unsafe_only_spaces(ptr::Ptr{UInt8}, i::Int=1)
    while true
        b = unsafe_load(ptr, i)
        iszero(b) && return true
        @eq(' ', b) || return false
        i += 1
    end
end

# Convert card data (e.g., from a key=>dat pair) into a 2-tuple (val,com).
# NOTE: An empty comment "" is the same as unspecified comment (nothing).
normalize_card_data(val::Any) = (normalize_card_value(val), nothing)
normalize_card_data(val::Tuple{Any}) = (normalize_card_value(val[1]), nothing)
normalize_card_data(val::Tuple{Any,OptionalString}) = (normalize_card_value(val[1]), val[2])

"""
   EasyFITS.normalize_card_value(val)

converts `val` to a suitable keyword value, yielding `Invalid()` for invalid
value type. Argument `val` shall be the bare value of a non-commentary FITS
keyword without the comment part.

"""
normalize_card_value(val::UndefinedValue) = missing
normalize_card_value(val::Nothing) = nothing
normalize_card_value(val::Bool) = val
normalize_card_value(val::Integer) = to_type(Int, val)
normalize_card_value(val::Real) = to_type(Cdouble, val)
normalize_card_value(val::Complex) = to_type(Complex{Cdouble}, val)
normalize_card_value(val::AbstractString) = val
normalize_card_value(val::Any) = Invalid() # means error

unsafe_optional_string(s::AbstractString) = s
unsafe_optional_string(s::Nothing) = Ptr{Cchar}(0)

@noinline invalid_card_value(key::CardName, val::Any, commentary::Bool=false) = bad_argument(
    "invalid value of type ", typeof(val), " for ", (commentary ? "commentary " : ""),
    "FITS keyword \"", normalize_card_name(key), "\"")

normalize_card_name(key::Symbol) = normalize_card_name(String(key))
normalize_card_name(key::AbstractString) = uppercase(rstrip(key))

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

for func in (:update_key, :write_key),
    (V, T) in ((AbstractString, String),
               (Bool,           Bool),
               (Integer,        Int),
               (Real,           Cdouble),
               (Complex,        Complex{Cdouble}),
               (UndefinedValue, Missing))
    if T === Missing
        # Commentary (nothing) card or undefined value (undef or missing).
        @eval function $func(dst, key::CardName, val::$V, com::OptionalString=nothing)
            _com = unsafe_optional_string(com)
            check(CFITSIO.$(Symbol("fits_",func,"_null"))(
                dst, key, _com, Ref{Status}(0)))
            return dst
        end
    elseif T === String
        @eval function $func(dst, key::CardName, val::$V, com::OptionalString=nothing)
            _com = unsafe_optional_string(com)
            check(CFITSIO.$(Symbol("fits_",func,"_str"))(
                dst, key, val, _com, Ref{Status}(0)))
            return dst
        end
    elseif T === Bool
        @eval function $func(dst, key::CardName, val::$V, com::OptionalString=nothing)
            _com = unsafe_optional_string(com)
            check(CFITSIO.$(Symbol("fits_",func,"_log"))(
                dst, key, val, _com, Ref{Status}(0)))
            return dst
        end
    else
        @eval function $func(dst, key::CardName, val::$V, com::OptionalString=nothing)
            _com = unsafe_optional_string(com)
            _val = Ref{$T}(val)
            check(CFITSIO.$(Symbol("fits_",func))(
                dst, type_to_code($T), key, _val, _com, Ref{Status}(0)))
            return dst
        end
    end
end

# HDUs are similar to ordered sets (lists) of header cards. NOTE: eltype,
# ndims, and size are used for the data part not the header part.
Base.length(hdu::FitsHDU) = get_hdrspace(hdu)[1]
Base.firstindex(hdu::FitsHDU) = 1
Base.lastindex(hdu::FitsHDU) = length(hdu)
Base.keys(hdu::FitsHDU) = Base.OneTo(length(hdu))
Base.iterate(hdu::FitsHDU, state::Int = firstindex(hdu)) =
    state ≤ length(hdu) ? (hdu[state], state + 1) : nothing
Base.IteratorSize(::Type{<:FitsHDU}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FitsHDU}) = Base.HasEltype() # FIXME: not Image HDU
Base.eltype(::Type{<:FitsHDU}) = FitsCard # FIXME: not Image HDU

# FIXME: implement iterator API?

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

"""
    nameof(hdu::FitsHDU) -> str

yields the name of the FITS header data unit `hdu`. The name is the value of
the first keyword of `"EXTNAME"` or `"HDUNAME"` which exists and has a string
value. Otherwise, the name is that of the FITS extension of `hdu`, that is
`"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"` depending on whether `hdu` is an
image, an ASCII table, a binary table, or anything else.

"""
function Base.nameof(hdu::FitsHDU)
    (str = hdu.hduname) === nothing || return str
    (str = hdu.extname) === nothing || return str
    return hdu.xtension
end

# Compare HDU name with some pattern , HDU name may be `nothing` which can
# never be considered as a success.
same_name(name::Nothing, pat::Union{AbstractString,Regex}) = false
same_name(name::AbstractString, pat::AbstractString) = isequal(FitsLogic(), name, pat)
same_name(name::AbstractString, pat::Regex) = match(pat, name) !== nothing

"""
    EasyFITS.is_named(hdu, pat) -> bool

yields whether pattern `pat` is equal to (in the FITS sense if `pat` is a
string) or matches (if `pat` is a regular expression) the extension of the FITS
header data unit `hdu`, or to the value of one of its `"EXTNAME"` or
`"HDUNAME"` keywords. These are respectively given by `hdu.xtension`,
`hdu.extname`, or `hdu.hduname`.

This method is used as a predicate for the search methods `findfirst`,
`findlast`, `findnext`, and `findprev`.

The extension `hdu.xtension` is `"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"`
depending on whether `hdu` is an image, an ASCII table, a binary table, or
anything else.

"""
is_named(hdu::FitsHDU, pat::Union{AbstractString,Regex}) =
    # Since a match only fails if no matching name is found, the order of the
    # tests is irrelevant. We therefore start with the costless ones.
    same_name(hdu.xtension, pat) ||
    same_name(hdu.hduname, pat) ||
    same_name(hdu.extname, pat)
is_named(pat::Union{AbstractString,Regex}) = Base.Fix2(is_named, pat)

for func in (:findfirst, :findlast)
    @eval Base.$func(pat::Union{AbstractString,Regex}, file::FitsFile) =
        $func(is_named(pat), file)
end

for func in (:findnext, :findprev)
    @eval Base.$func(pat::Union{AbstractString,Regex}, file::FitsFile, start::Integer) =
        $func(is_named(pat), file, start)
end

function Base.findfirst(f::Function, file::FitsFile)
    @inbounds for i ∈ keys(file)
        f(file[i]) && return i
    end
    return nothing
end

function Base.findlast(f::Function, file::FitsFile)
    @inbounds for i ∈ reverse(keys(file))
        f(file[i]) && return i
    end
    return nothing
end

function Base.findnext(f::Function, file::FitsFile, start::Integer)
    start = to_type(keytype(file), start)
    start < firstindex(file) && throw(BoundsError(file, start))
    @inbounds for i ∈ start:lastindex(file)
        f(file[i]) && return i
    end
    return nothing
end

function Base.findprev(f::Function, file::FitsFile, start::Integer)
    start = to_type(keytype(file), start)
    start > lastindex(file) && throw(BoundsError(file, start))
    @inbounds for i ∈ start:-1:firstindex(file)
        f(file[i]) && return i
    end
    return nothing
end
