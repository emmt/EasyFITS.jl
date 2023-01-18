# Management of FITS cards.

function Base.get(hdu::FitsHDU, key::Union{CardName,Integer}, def)
    card = FitsCard()
    status = try_read!(card, hdu, key)
    iszero(status) && return initialize!(card)
    status == (key isa Integer ?
        CFITSIO.KEY_OUT_BOUNDS : CFITSIO.KEY_NO_EXIST) || throw(FitsError(status))
    return def
end

function FitsCard(hdu::FitsHDU, key::Union{CardName,Integer})
    card = FitsCard()
    status = try_read!(card, hdu, key)
    iszero(status) && return initialize!(card)
    status == (key isa Integer ?
        CFITSIO.KEY_OUT_BOUNDS : CFITSIO.KEY_NO_EXIST) || throw(FitsError(status))
    throw(KeyError(key))
end

# Unchecked read of a FITS header card.
try_read!(card::FitsCard, hdu::FitsHDU, key::CardName) =
    CFITSIO.fits_read_card(hdu, key, pointer(card), Ref{Status}(0))
function try_read!(card::FitsCard, hdu::FitsHDU, key::Integer)
    key > 0 || return Status(CFITSIO.KEY_OUT_BOUNDS)
    return CFITSIO.fits_read_record(hdu, key, pointer(card), Ref{Status}(0))
end

let offset = fieldoffset(FitsCard, get_field_index(FitsCard, :data))
    @eval Base.pointer(A::FitsCard) = Ptr{UInt8}(pointer_from_objref(A)) + $offset
end

# Implement abstract array API for FITS cards as a vector of bytes.
Base.firstindex(A::FitsCard) = 1
Base.lastindex(A::FitsCard) = 80
Base.length(A::FitsCard) = 80
Base.IndexStyle(::Type{<:FitsCard}) = IndexLinear()
Base.size(A::FitsCard) = (length(A),)
Base.axes(A::FitsCard) = (keys(A),)
Base.keys(A::FitsCard) = Base.OneTo(length(A))
@inline function Base.getindex(A::FitsCard, i::Int)
    @boundscheck checkbounds(A, i)
    return GC.@preserve A unsafe_load(pointer(A), i)
end
@inline function Base.setindex!(A::FitsCard, x, i::Int)
    @boundscheck checkbounds(A, i)
    GC.@preserve A unsafe_store!(pointer(A), convert(eltype(A), x), i)
    return A
end

# Private accessors.  The end-user only sees the properties.
get_type(A::FitsCard) = getfield(A, :type)
get_type(A::FitsCardValue) = get_type(parent(A))
set_type!(A::FitsCard, val::FitsCardType) = setfield!(A, :type, val)
for (part, child) in (("name",    FitsCardName),
                      ("value",   FitsCardValue),
                      ("comment", FitsCardComment))
    for field in ("offset", "length")
        part_field = part*"_"*field
        i = get_field_index(FitsCard, part_field)
        @eval begin
            $(Symbol("get_"*field))(A::$child) = getfield(parent(A), $i)
            $(Symbol("get_"*part_field))(A::FitsCard) = getfield(A, $i)
            $(Symbol("set_"*part_field*"!"))(A::FitsCard, x::Integer) =
                setfield!(A, $i, Int(x))
        end
    end
end

Base.show(io::IO, A::FitsCard) = print_struct(print, io, A)

function Base.show(io::IO, ::MIME"text/plain", A::FitsCard)
    print(io, '"')
    print_struct(print_esc, io, A)
    print(io, '"')
end

# Print structure contents in a human readable form with an auxiliary function
# to print parts.
function print_struct(print_part::Function, io::IO, A::FitsCard)
    print(io, A.name)
    len = length(A.name)
    while len < 8
        print(io, ' ')
        len += 1
    end
    if A.type == FITS_COMMENT
        print(io, ' ')
        print_part(io, A.comment)
    elseif A.type ∈ (FITS_LOGICAL, FITS_INTEGER, FITS_FLOAT, FITS_COMPLEX,
                     FITS_STRING, FITS_UNDEFINED)
        print(io, "= ")
        print_part(io, A.value)
        if length(A.comment) > 0
            print(io, " / ")
            print_part(io, A.comment)
        end
    end
end

# Like print but escape some characters.
function print_esc(io::IO, s::FitsCardPart)
    for c in s
        if (c == '"') | (c == '\\')
            print(io, '\\')
        end
        print(io, c)
    end
    return io
end

function Base.show(io::IO, mime::MIME"text/plain", s::FitsCardValue)
    type = s.type
    if type == FITS_LOGICAL
        show(io, mime, s.logical)
    elseif type == FITS_INTEGER
        show(io, mime, s.integer)
    elseif type == FITS_FLOAT
        show(io, mime, s.float)
    elseif type == FITS_COMPLEX
        show(io, mime, s.complex)
    elseif type == FITS_STRING
        show(io, mime, s.string)
    elseif type == FITS_COMMENT
        show(io, mime, nothing)
    elseif type == FITS_UNDEFINED
        show(io, mime, missing)
    end
end

#------------------------------------------------------------------------------
# PARSING OF FITS CARD.
#
# Even though FITS standard imposes that only ASCII characters are allowed in
# FITS header cards, we may allow for UTF8 encoding. Since the structure of a
# FITS card is entirely determined by the occurrence of ordinary ASCII
# characters, the FITS card can be accessed as an array of bytes during its
# parsing. Things are a bit more complicated for long strings and comments that
# may be split across several cards.
#
# For these reasons a `FitsCard` is an abstract array of bytes while its
# different parts (its name, value, and comment or respective types
# `FitsCardName`, `FitsCardValue`, and `FitsCardComment`) are abstract UTF8
# strings. For efficiency, these parts are represented by simple decorations
# that wrap their parent FITS card.

"""
    @starts_with buf str

yields unrolled code that checks whether vector of bytes bound to variable
`buf` starts with the characters of the literal string `str`.

!!! warning:
    Element type of `buf` must be `UInt8` and its first index must be 1, none
    of these can be checked.

"""
macro starts_with(buf::Symbol, str::AbstractString)
    # NOTE: This assumes that buf has first index 1.
    return esc(_starts_with(buf, 1, str, firstindex(str), lastindex(str)))
end

function _starts_with(buf::Symbol, i::Int, str::AbstractString, j::Int, last::Int)
    ex = :($buf[$i] == $(UInt8(str[j])))
    if j < last
        return Expr(:(&&), ex, _starts_with(buf, i + 1, str, nextind(str, j), last))
    else
        return ex
    end
end

# A few methods to compare characters and avoid conversion of bytes to Char in
# ASCII strings.
between(b::UInt8, lo::Char, hi::Char) = between(b, UInt8(lo), UInt8(hi))
between(x::T, lo::T, hi::T) where {T<:Integer} = (x ≥ lo) & (x ≤ hi)
equal(b::UInt8, c::Char) = (b == UInt8(c))

# NOTE: Trimming in a FITS card only concerns ordinary ASCII spaces so that the
#       string can be accessed as an ASCII string.

@inline @propagate_inbounds function trim_left(card::FitsCard, r::AbstractUnitRange{Int})
    i_first, i_last = first(r), last(r)
    while i_first ≤ i_last && equal(card[i_first], ' ')
        i_first += 1
    end
    return i_first : i_last
end

@inline @propagate_inbounds function trim_right(card::FitsCard, r::AbstractUnitRange{Int})
    i_first, i_last = first(r), last(r)
    while i_last ≥ i_first && equal(card[i_last], ' ')
        i_last -= 1
    end
    return i_first : i_last
end

@inline @propagate_inbounds function trim(card::FitsCard, r::AbstractUnitRange{Int})
    i_first, i_last = first(r), last(r)
    while i_first ≤ i_last && equal(card[i_first], ' ')
        i_first += 1
    end
    while i_last ≥ i_first && equal(card[i_last], ' ')
        i_last -= 1
    end
    return i_first : i_last
end

# NOTE: In the 2 following methods, it is assumed that first(r) is the 1st
#       possible index for the value part while last(r) is the position of the
#       last non-space character of the card.

@propagate_inbounds function find_value_and_comment(card::FitsCard, r::AbstractUnitRange{Int})
    # Find value, then find beginning of comment.
    type, value = find_value(card, r)
    flag = false # to remember whether a slash was already found
    i = last(value) + 1
    while i ≤ last(r)
        c = card[i]
        if ! equal(c, ' ')
            if equal(c, '/')
                flag && break
                flag = true
            else
                break
            end
        end
        i += 1
    end
    comment = i : last(r)
    return type, value, comment
end

@propagate_inbounds function find_value(card::FitsCard, r::AbstractUnitRange{Int})
    for i in r
        c = card[i]
        if equal(c, ' ')
            # Skip leading spaces.
            continue
        elseif equal(c, '\'')
            # Quoted string value. Find the closing quote. NOTE: The search
            # loop cannot be a for-loop because index j is incremented twice
            # when an escaped quote is encountered.
            j = i
            while j < last(r)
                j += 1
                if equal(card[j], '\'')
                    if j < last(r) && equal(card[j+1], '\'')
                        # Skip escaped quote.
                        j += 1
                    else
                        # Current character is the closing quote.
                        return FITS_STRING, i:j
                    end
                end
            end
            error("no closing quote in string value")
        elseif equal(c, '(')
            # Complex value.
            for j in i+1:last(r)
                if equal(card[j], ')')
                    return FITS_COMPLEX, i:j
                end
            end
            error("no closing parenthesis in complex value")
        elseif equal(c, 'F') | equal(c, 'T')
            return FITS_LOGICAL, i:i
        elseif equal(c, '+') | equal(c, '-') | equal(c, '.') | between(c, '0', '9')
            # Integer or float value.
            type = equal(c, '.') ? FITS_FLOAT : FITS_INTEGER
            for j in i + 1 : last(r)
                c = card[j]
                if equal(c, ' ') | equal(c, '/')
                    return type, i : j - 1
                elseif ! between(c, '0', '9')
                    type = FITS_FLOAT
                end
            end
            return type, i : last(r)
        elseif equal(c, '/')
            # Comment marker found before value, the value is undefined.
            #
            # NOTE: The position of the last character n the value part is used later for
            # finding the comment part, so we set the first index of the value part so as to have an empty range.
            return FITS_UNDEFINED, i : i - 1
        elseif MANDATORY_COMMENT_SEPARATOR
            error("unexpected leading character '$c' for value part")
        else
            break
        end
    end
    return FITS_UNDEFINED, last(r) + 1 : last(r)
end

function initialize!(card::FitsCard)
    # Find end of card, trimming trailing spaces. NOTE: This is assumed by the
    # helper methods find_value and find_value_and_comment.
    r = firstindex(card) : lastindex(card)
    mark = first(r) - 1 # position of last non-space
    for i in r
        c = card[i]
        if ! equal(c, ' ')
            iszero(c) && break
            mark = i # remember position of last non-space
        end
    end
    r = first(r) : mark # range of allowed indices

    # Find the different parts of the card: its name, value, and comment.
    # According to FITS standard:
    #
    # * Card name consists in upper case Latin letters (A:Z), digits (0:9),
    #   hyphen (_) or underscore (_) characters, it is left justified in
    #   columns 1:8, right padded with ordinary spaces. Exception: the HIERARCH
    #   convention.
    #
    # * If value indicator "= " is present in columns 9:10, the card may have a
    #   value starting at column 11, followed by an optional comment. Value and
    #   comment are separated by a slash /. Otherwise, the card is a comment.
    #
    # The procedure shall give the same result as `ffpsvc` (a.k.a.
    # `fits_parse_value`) in `fitscore.c` of the CFITSIO library but extracting
    # the card name as well and sparing the needs to allocate buffers for the
    # value and comment strings.
    #
    # NOTE: To speed-up the procedure, the card is assumed to be syntactically
    #       correct and a minimal number of checks are performed (similar
    #       assumptions are made in ffpsvc).
    len = length(r)
    if len ≥ 9 && @starts_with(card, "HIERARCH ")
        # HIERARCH convention. Skip leading spaces and search for the value
        # marker, a '=' character, keeping track of the position of the last
        # non-space character.
        mark = last(r) + 1
        for i in first(r) + 9 : last(r)
            if equal(card[i], '=')
                mark = i
                break
            end
        end
        if mark ≤ last(r)
            name = trim(card, first(r) + 9 : mark - 1)
            if isempty(name)
                name = first(r) : first(r) + 7
            end
            type, value, comment = find_value_and_comment(card, mark + 1 : last(r))
        else
            # No value marker found. Assume a commentary card whose name is
            # "HIERARCH".
            type = FITS_COMMENT
            name = first(r) : first(r) + 7
            value = first(r) : first(r) - 1 # empty value
            comment = first(r) + 8 : last(r)
        end
    else
        # Ordinary key name in columns 1:8. Check for value indicator in columns 9:10.
        name = trim_right(card, first(r) : min(first(r) + 7, last(r)))
        if len < 10 || ! equal(card[first(r) + 8], '=') || ! equal(card[first(r) + 9], ' ')
            # Commentary card starting in column 9 and without trailing spaces.
            comment = first(r) + 8 : last(r)
            type = FITS_COMMENT
            value = first(r) : first(r) - 1 # empty value
        else
            type, value, comment = find_value_and_comment(card, first(r) + 10 : last(r))
        end
    end

    # Instantiate card fields.
    set_type!(card, type)
    set_name_offset!(card, first(name) - first(r))
    set_name_length!(card, length(name))
    set_value_offset!(card, first(value) - first(r))
    set_value_length!(card, length(value))
    set_comment_offset!(card, first(comment) - first(r))
    set_comment_length!(card, length(comment))
    return card
end

#------------------------------------------------------------------------------
# Methods for FITS card parts.

Base.parent(s::FitsCardPart) = getfield(s, :parent)
Base.pointer(s::FitsCardPart) = pointer(parent(s)) + get_offset(s)

# Implement abstract string API for FITS card parts.
Base.ncodeunits(s::FitsCardPart) = get_length(s)
Base.codeunit(s::FitsCardPart) = UInt8
@inline function Base.codeunit(s::FitsCardPart, i::Int)
    @boundscheck (i % UInt) - 1 < ncodeunits(s) || throw(BoundsError(s, i))
    return @inbounds getindex(parent(s), i + get_offset(s))
end
function Base.isascii(s::FitsCardPart)
    @inbounds for i = 1:ncodeunits(s)
        codeunit(s, i) >= 0x80 && return false
    end
    return true
end

# Code for FITS cards as ASCII and UTF8 strings.
include("ascii.jl")
include("utf8.jl")

Base.getindex(s::FitsCardPart, r::AbstractUnitRange{<:Integer}) = s[Int(first(r)):Int(last(r))]

@inline function Base.getindex(s::FitsCardPart, r::UnitRange{Int})
    isempty(r) && return ""
    i, j = first(r), last(r)
    @boundscheck begin
        checkbounds(s, r)
        @inbounds isvalid(s, i) || string_index_err(s, i)
        @inbounds isvalid(s, j) || string_index_err(s, j)
    end
    off = i - firstindex(s) # offset in bytes
    siz = nextind(s, j) - i # number in bytes
    return GC.@preserve s unsafe_string(pointer(s) + off, len)
end

Base.propertynames(::FitsCard) = (:comment, :complex, :name, :pair, :parsed,
                                  :string, :type, :value, :units)
@inline Base.getproperty(card::FitsCard, key::Symbol) = getproperty(card, Val(key))

Base.getproperty(card::FitsCard, ::Val{:pair}) = Pair(card)
Base.getproperty(card::FitsCard, ::Val{:parsed}) = parsed_value(card)
Base.getproperty(card::FitsCard, ::Val{:type}) = get_type(card)
Base.getproperty(card::FitsCard, ::Val{:name}) = FitsCardName(card)
Base.getproperty(card::FitsCard, ::Val{:value}) = FitsCardValue(card)
Base.getproperty(card::FitsCard, ::Val{:comment}) = FitsCardComment(card)
for key in (:logical, :integer, :float, :complex, :string)
    @eval @inline Base.getproperty(card::FitsCard, ::$(Val{key})) = card.value.$key
end
for key in (:units, :unitless)
    @eval @inline Base.getproperty(card::FitsCard, ::$(Val{key})) = card.comment.$key
end

Base.propertynames(::FitsCardComment) = (:units, :unitless)

@inline Base.getproperty(A::FitsCardComment, key::Symbol) =
    key === :units ? units(A) :
    key === :unitless ? unitless(A) :
    invalid_property(A, key)

function find_units_marks(A::FitsCardComment)
    i_first, i_last = firstindex(A), lastindex(A)
    i = i_first
    # FIXME: @inbounds
    while i ≤ i_last
        c = A[i]
        j = nextind(A, i)
        if c == '['
            while j ≤ i_last
                A[j] == ']' && return i:j
                j = nextind(A, j)
            end
        elseif c != ' '
            break
        end
        i = j
    end
    return i_first : prevind(A, i_first)
end

function empty_substring(s::AbstractString)
    i = firstindex(s)
    return SubString(s, i, prevind(s, i))
end

function units(A::FitsCardComment)
    r = find_units_marks(A)
    if !isempty(r)
        i = nextind(A, first(r)) # skip opening [
        j = prevind(A, last(r))  # skip closing ]
        # FIXME: @inbounds
        while i ≤ j
            A[j] == ' ' || return SubString(A, i, j)
            j = prevind(A, j)
        end
    end
    return empty_substring(A)
end

function unitless(A::FitsCardComment)
    r = find_units_marks(A)
    i_last = lastindex(A)
    isempty(r) && return SubString(A, firstindex(A), i_last)
    i = nextind(A, last(r))
    # FIXME: @inbounds
    while i ≤ i_last
        A[i] == ' ' || return SubString(A, i, i_last)
        i = nextind(A, i)
    end
    return empty_substring(A)
end

Base.convert(::Type{T}, card::FitsCard) where {T<:Pair} = T(card)

Base.Pair(card::FitsCard) = card.name => (parsed_value(card), card.comment)
Base.Pair{K}(card::FitsCard) where {K<:AbstractString} =
    convert(K, card.name) => (parsed_value(card), card.comment)
Base.Pair{K,Tuple{T,S}}(card::FitsCard) where {K<:AbstractString,T,S} =
    convert(K, card.name) => (convert(T, card.value),
                              convert_comment(S, card.comment))

convert_comment(::Type{T}, comment::T) where {T} = comment
convert_comment(::Type{T}, comment::AbstractString) where {T<:AbstractString} = convert(T, comment)
convert_comment(::Type{T}, comment::Nothing) where {T<:AbstractString} = convert(T, "")
convert_comment(::Type{T}, comment::OptionalString) where {T<:Nothing} = nothing

Base.propertynames(::FitsCardValue) = (:complex, :float, :integer, :logical, :parsed, :string, :type)

@inline Base.getproperty(A::FitsCardValue, key::Symbol) =
    key === :logical ? convert(Bool,             A) :
    key === :integer ? convert(Int,              A) :
    key === :float   ? convert(Cdouble,          A) :
    key === :complex ? convert(Complex{Cdouble}, A) :
    key === :parsed  ? parsed_value(             A) :
    key === :string  ? convert(String,           A) :
    key === :type    ? get_type(A)                  :
    invalid_property(A, key)

@noinline convert_value_error(::Type{T}, A::FitsCardValue) where {T} =
    error("cannot parse value of $(type_name(A)) FITS card as $T")

"""
    EasyFITS.parsed_value(A::Union{FitsCard,FitsCardValue})

yields the parsed value of `A`, a FITS card or a FITS card value.

!!! warning
    This function is not type-stable.

"""
parsed_value(card::FitsCard) = parsed_value(card.value)
parsed_value(value::FitsCardValue) =
    value.type == FITS_LOGICAL   ? value.logical :
    value.type == FITS_INTEGER   ? value.integer :
    value.type == FITS_FLOAT     ? value.float   :
    value.type == FITS_COMPLEX   ? value.complex :
    value.type == FITS_STRING    ? value.string  :
    value.type == FITS_COMMENT   ? nothing       :
    value.type == FITS_UNDEFINED ? missing       : error("unknown type of FITS card value")

type_name(A::Union{FitsCard,FitsCardValue}) = type_name(get_type(A))
type_name(type::FitsCardType) =
    type == FITS_LOGICAL   ? "LOGICAL"   :
    type == FITS_INTEGER   ? "INTEGER"   :
    type == FITS_FLOAT     ? "FLOAT"     :
    type == FITS_COMPLEX   ? "COMPLEX"   :
    type == FITS_STRING    ? "STRING"    :
    type == FITS_COMMENT   ? "COMMENT"   :
    type == FITS_UNDEFINED ? "UNDEFINED" : "UNKNOWN"

for T in (:Bool, :UInt8, :Int8, :UInt16, :Int16, :UInt32, :Int32, :UInt64, :Int64,
          :UInt128, :Int128, :BigInt, :Float32, :Float64, :BigFloat, :String)
    if isdefined(Base, T)
        @eval begin
            Base.$T(A::FitsCardValue) = convert($T, A)
            Base.$T(A::FitsCard) = convert($T, FitsCardValue(A))
        end
    end
end

function Base.convert(::Type{T}, A::FitsCardValue) where {T<:Number}
    type = A.type
    if type == FITS_LOGICAL
        val = _try_parse_logical(A)
        val !== nothing && return convert(T, val)
    elseif type == FITS_INTEGER
        val = _try_parse_integer(A)
        val !== nothing && return convert(T, val)
    elseif type == FITS_FLOAT
        val = _try_parse_float(A)
        val !== nothing && return convert(T, val)
    elseif type == FITS_COMPLEX
        val = _try_parse_complex(A)
        val !== nothing && return convert(T, val)
    end
    convert_value_error(T, A)
end

function Base.convert(::Type{T}, A::FitsCardValue) where {T<:String}
    if A.type == FITS_STRING
        val = _try_parse_string(A)
        val !== nothing && return val
    end
    convert_value_error(T, A)
end

function Base.convert(::Type{T}, A::FitsCardValue) where {T<:Nothing}
    A.type == FITS_COMMENT && return T()
    convert_value_error(T, A)
end

function Base.convert(::Type{T}, A::FitsCardValue) where {T<:UndefinedValue}
    A.type == FITS_UNDEFINED && return T()
    convert_value_error(T, A)
end

# Private functions to convert the value of a FITS card. It is the caller's
# responsibility to check that the card type is correct.
_try_parse_logical(A::FitsCardValue) = first(A) == 'T'

_try_parse_integer(A::FitsCardValue) = tryparse(Int, A)

function _try_parse_float(A::FitsCardValue)
    io = IOBuffer()
    # FIXME: @inbounds
    for c in A
        write(io, ifelse((c == 'D') | (c == 'd'), oftype(c, 'E'), c))
    end
    return tryparse(Cdouble, String(take!(io)))
end

function _try_parse_complex(A::FitsCardValue)
    i = firstindex(A)
    i_last = lastindex(A)
    state = 0
    @inbounds while i ≤ i_last
        c = A[i]
        if c != ' '
            if c == '('
                state = 1
            else
                break
            end
        end
        i = nextind(A, i)
    end
    state == 1 || return nothing
    io = IOBuffer()
    @inbounds while i ≤ i_last
        c = A[i]
        if (c == 'D') | (c == 'd')
            write(io, 'E')
        elseif c == ','
            state = 2
            i = nextind(A, i)
            break
        else
            write(io, c)
        end
        i = nextind(A, i)
    end
    state == 2 || return nothing
    re = tryparse(Cdouble, String(take!(io)))
    re === nothing && return nothing
    @inbounds while i ≤ i_last
        c = A[i]
        if (c == 'D') | (c == 'd')
            write(io, 'E')
        elseif c == ')'
            @inbounds while i < i_last
                # Check that there are only trailing spaces.
                i = nextind(A, i)
                A[i] == ' ' || return nothing
            end
            state = 3
            break
        else
            write(io, c)
        end
        i = nextind(A, i)
    end
    state == 3 || return nothing
    im = tryparse(Cdouble, String(take!(io)))
    im === nothing && return nothing
    return complex(re, im)
end

function _try_parse_string(A::FitsCardValue)
    i_first = firstindex(A)
    length(A) > 0 && A[i_first] == '\'' || return nothing
    i_last = lastindex(A)
    io = IOBuffer()
    i = i_first
    nspaces = 0
    while i < i_last
        i = nextind(A, i)
        c = A[i] # FIXME: @inbounds
        if c == ' '
            nspaces += 1
        else
            if c == '\''
                if i < i_last
                    i = nextind(A, i)
                    c = A[i] # FIXME: @inbounds
                    c == '\'' || break
                    while nspaces > 0
                        write(io, ' ')
                        nspaces -= 1
                    end
                    write(io, c)
                else
                    return String(take!(io))
                end
            else
                while nspaces > 0
                    write(io, ' ')
                    nspaces -= 1
                end
                write(io, c)
            end
        end
    end
    return nothing
end

"""
    FitsCardType(T)

yields the FITS header card type code corresponding to Julia type `T`. The
returned value is one of: `FITS_UNDEFINED`, `FITS_LOGICAL`, `FITS_INTEGER`,
`FITS_FLOAT`, `FITS_COMPLEX`, `FITS_STRING`, `FITS_COMMENT`, or `FITS_UNKNOWN`.
The latter indicates unsupported type of card value.

---

    FitsCardType(x::Union{FitsCard,FitsCardValue})
    x.type

yield the type code of the FITS header card or value `x`.

"""
FitsCardType(::Type)                   = FITS_UNKNOWN
FitsCardType(::Type{<:UndefinedValue}) = FITS_UNDEFINED
FitsCardType(::Type{<:Bool})           = FITS_LOGICAL
FitsCardType(::Type{<:Integer})        = FITS_INTEGER
FitsCardType(::Type{<:AbstractFloat})  = FITS_FLOAT
FitsCardType(::Type{<:Complex})        = FITS_COMPLEX
FitsCardType(::Type{<:AbstractString}) = FITS_STRING
FitsCardType(::Type{<:Nothing})        = FITS_COMMENT
FitsCardType(x::Union{FitsCard,FitsCardValue}) = get_type(x)
