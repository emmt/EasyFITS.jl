# Management of FITS header data units (HDUs) and of FITS header cards.

"""
    FitsHDU(io::FitsIO, i) -> hdu
    getindex(io::FitsIO, i) -> hdu
    io[i] -> hdu

yield the `i`-th FITS header data unit of FITS file `io`.

The returned object has the following read-only properties:

    hdu.io       # same as FitsIO(hdu)
    hdu.num      # yields the HDU number (the index `i` above)
    hdu.type     # same as FITSHDUType(hdu)
    hdu.xtension # value of the XTENSION card (never nothing)
    hdu.extname  # value of the EXTNAME card or nothing
    hdu.hduname  # value of the HDUNAME card or nothing

"""
function FitsHDU(io::FitsIO, i::Integer)
    ptr = check(pointer(io))
    status = Ref{Status}(0)
    type = Ref{Cint}()
    check(CFITSIO.fits_movabs_hdu(ptr, i, type, status))
    type = FitsHDUType(type[])
    if type == FITS_ASCII_TABLE_HDU
        return FitsTableHDU(BareBuild(), io, i, true)
    elseif type == FITS_BINARY_TABLE_HDU
        return FitsTableHDU(BareBuild(), io, i, false)
    elseif type == FITS_IMAGE_HDU
        bitpix = Ref{Cint}()
        check(CFITSIO.fits_get_img_equivtype(ptr, bitpix, status))
        ndims = Ref{Cint}()
        check(CFITSIO.fits_get_img_dim(ptr, ndims, status))
        N = Int(ndims[])::Int
        T = type_from_bitpix(bitpix[])
        return FitsImageHDU{T,N}(BareBuild(), io, i)
    else
        return FitsAnyHDU(BareBuild(), io, i)
    end
end

FitsHDU(io::FitsIO, s::AbstractString) = io[s]

FitsTableHDU(io::FitsIO, i::Union{Integer,AbstractString}) =
    FitsHDU(io, i)::FitsTableHDU

FitsImageHDU(io::FitsIO, i::Union{Integer,AbstractString}) =
    FitsHDU(io, i)::FitsImageHDU
FitsImageHDU{T}(io::FitsIO, i::Union{Integer,AbstractString}) where {T} =
    FitsHDU(io, i)::FitsImageHDU{T}
FitsImageHDU{T,N}(io::FitsIO, i::Union{Integer,AbstractString}) where {T,N} =
    FitsHDU(io, i)::FitsImageHDU{T,N}

Base.propertynames(::FitsHDU) = (:extname, :hduname, :io, :num, :type, :xtension)
Base.getproperty(hdu::FitsHDU, sym::Symbol) = getproperty(hdu, Val(sym))

Base.getproperty(hdu::FitsHDU, ::Val{:extname})  = get_extname(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:hduname})  = get_hduname(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:io})       = FitsIO(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:num})      = getfield(hdu, :num)
Base.getproperty(hdu::FitsHDU, ::Val{:type})     = FitsHDUType(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{:xtension}) = get_xtension(hdu)
Base.getproperty(hdu::FitsHDU, ::Val{sym}) where {sym} = invalid_property(hdu, sym)

Base.setproperty!(hdu::FitsHDU, sym::Symbol, x) =
    sym ∈ propertynames(hdu) ? readonly_property(hdu, sym) : invalid_property(hdu, sym)

get_xtension(hdu::FitsHDU) = "ANY"
get_xtension(hdu::FitsImageHDU) = "IMAGE"
get_xtension(hdu::FitsTableHDU) = isascii(hdu) ? "TABLE" : "BINTABLE"

for (func, key) in ((:get_extname, "EXTNAME"),
                    (:get_hduname, "HDUNAME"))
    @eval function $func(hdu::FitsHDU)
        card = get(hdu, $key, nothing)
        card === nothing && return nothing
        card.type == FITS_STRING || return nothing
        return card.value.string
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
    FitsIO(hdu::FitsHDU) -> io
    hdu.io -> io

yield the FITS file associated with FITS Header Data Unit `hdu` moving the file
position of `io` to that of `hdu` and checking that the file is still open.

"""
function FitsIO(hdu::FitsHDU)
    io = getfield(hdu, :io)
    check(CFITSIO.fits_movabs_hdu(io, getfield(hdu, :num), Ptr{Cint}(0), Ref{Status}(0)))
    return io
end

"""
    pointer(hdu::FitsHDU)

yields the pointer to the FITS file associated with FITS Header Data Unit `hdu`
moving the file position to that of `hdu` and checking that the file is still
open. It is the caller's responsibility to make sure that the pointer remains
valid, e.g., by wrapping the code in a `GC@preserve` block.

"""
Base.pointer(hdu::FitsHDU) = pointer(FitsIO(hdu))

# Extend unsafe_convert to automatically extract and check the FITS file handle
# from a FitsHDU object. This secures and simplifies calls to
# functions of the CFITSIO library.
Base.unsafe_convert(::Type{Ptr{CFITSIO.fitsfile}}, hdu::FitsHDU) = pointer(hdu)

"""
    getindex(hdu::FitsHDU, key) -> card
    hdu[key] -> card

yield the FITS header card at index `key` (can be an integer or a string) in
header data unit `hdu`. The result is a small object with the following
properties:

    card.type      # card type (an enumeration `FitsCardType`)
    card.name      # card name as an abstract string
    card.value     # card unparsed value as an abstract string
    card.parsed    # parsed value (FIXME: not type-stable)
    card.pair      # an equivalent key=>(val,com) pair (FIXME: not type-stable)
    card.integer   # card value as an integer
    card.logical   # card value as a boolean
    card.float     # card value as a floating-point
    card.complex   # card value as a complex
    card.string    # card value as an abstract string
    card.comment   # card comment as an abstract string (leading and trailing spaces stripped)
    card.units     # card units as an abstract string
    card.unitless  # card comment without units part if any

Note that the card parts `card.name`, `card.value`, and `card.comment` are simple
objects. The 2 latter have properties:

    card.value.integer     # card value as an integer
    card.value.logical     # card value as a boolean
    card.value.float       # card value as a floating-point
    card.value.complex     # card value as a complex
    card.value.string      # card value as an abstract string
    card.value.parsed      # parsed card value (FIXME: not type-stable)

    card.comment.units     # card units as an abstract string
    card.comment.unitless  # card comment without units part if any

To retrieve the card value with type `T`, call `convert(T,card)` or
`convert(T,card.value)`. For the most common types, just call `T(card)` or
`T(card.value)`. If `T` is `Bool`, `Int`, `Float64`, `ComplexF64`, or `String`,
then `convert(T,x)` and `T(x)` with `x=card` or `x=card.value` are respectively
equivalent to `x.logical`, `x.integer`, `x.float`, `x.complex`, or `x.string`.
For example:

    convert(Int, hdu["NAXIS"])
    convert(Int, hdu["NAXIS"].value)
    Int(hdu["NAXIS"])
    Int(hdu["NAXIS"].value)
    hdu["NAXIS"].integer
    hdu["NAXIS"].value.integer

are all equivalent. Which one is to be preferred is a matter of taste and
readability.

"""
Base.getindex(hdu::FitsHDU, key::Union{CardName,Integer}) = FitsCard(hdu, key)

"""
    reset(hdu::FitsHDU) -> hdu

reset the search by wild card character in FITS header of `hdu` to the first
record.

"""
function Base.reset(hdu::FitsHDU)
    check(CFITSIO.fits_read_record(hdu, 0, Ptr{Cchar}(0), Ref{Status}(0)))
    return hdu
end

@inline function set_key(f::Function, hdu::FitsHDU, key::CardName, dat::CardData)
    if ! isa(dat, Tuple) && isa(dat, AbstractString) && is_comment_keyword(key)
        # Specific case of a commentary card.
        f(hdu, key, nothing, dat)
    else
        val, com = normalize_card_data(dat)
        f(hdu, key, val, com)
    end
    return hdu
end

"""
    setindex!(hdu::FitsHDU, dat, key) -> hdu
    hdu[key] = dat

updates or appends a record associating the keyword `key` with the data `dat`
in the header of the FITS header data unit `hdu`. See [`EasyFITS.Header`](@ref)
for the possible forms of `dat`.

"""
Base.setindex!(hdu::FitsHDU, dat::CardData, key::CardName) = set_key(update_key, hdu, key, dat)

"""
    push!(hdu::FitsHDU, key => dat) -> hdu
    hdu[key] = dat

appends a new record associating the keyword `key` with the data `dat` in the
header of the FITS header data unit `hdu`. See [`EasyFITS.Header`](@ref) for
the possible forms of such pairs.

A vector of pairs may be specified to push more than one record in a single
call:

    push!(hdu::FitsHDU, [key1 => dat1, key2 => val2, ...]) -> hdu

"""
Base.push!(hdu::FitsHDU, pair::Pair{<:CardName,<:CardData}) = set_key(write_key, hdu, pair...)
Base.push!(hdu::FitsHDU, ::Nothing) = hdu
function Base.push!(hdu, pairs::Header)
    for pair in pairs
        push!(hdu, pair)
    end
    return hdu
end

"""
    delete!(hdu::FitsHDU, key) -> hdu

deletes from FITS header data unit `hdu` the header card identified by `key`
the card name or number.

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

"""
   EasyFITS.normalize_card_value(x)

converts `x` to a suitable keyword value.

"""
normalize_card_value(x::Nothing) = x
normalize_card_value(x::Bool) = x
normalize_card_value(x::Integer) = to_type(Int, x)
normalize_card_value(x::Real) = to_type(Cdouble, x)
normalize_card_value(x::Complex) = to_type(Complex{Cdouble}, x)
normalize_card_value(x::AbstractString) = x

# Convert card data (e.g., from a key=>dat pair) into a 2-tuple (val,com).
normalize_card_data(dat::CardValue) = (normalize_card_value(dat), nothing)
normalize_card_data(dat::Tuple{CardValue}) = (normalize_card_value(dat[1]), nothing)
normalize_card_data(dat::Tuple{CardValue,OptionalString}) = (normalize_card_value(dat[1]), dat[2])

unsafe_optional_string(s::AbstractString) = s
unsafe_optional_string(s::Nothing) = Ptr{Cchar}(0)

"""
    EasyFITS.write_key(dst, key, val, com=nothing) -> dst

appends a new FITS header card in `dst` associating value `val` and comment
`com` to the keyword `key`.

""" write_key

"""
    EasyFITS.update_key(dst, key, val, com=nothing) -> dst

updates or appends a FITS header card in `dst` associating value `val` and
comment `com` to the keyword `key`.

""" update_key

for func in (:update_key, :write_key),
    (V, T) in ((AbstractString, String),
               (Bool, Bool),
               (Integer, Int),
               (Real, Cdouble),
               (Complex, Complex{Cdouble}),
               (Union{Nothing,UndefInitializer,Missing}, Nothing))
    local f = cfunc(func === :update_key ? "ffuky" : "ffpky", T)
    if T === Nothing
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

"""
    EasyFITS.modify_comment(hdu::FitsHDU, key, com) -> hdu

modifies the comment of the FITS header card of `hdu` whose name is `key`. The
new comment is specified by `com`.

"""
function modify_comment(hdu::FitsHDU, key::CardName, com::AbstractString)
    check(CFITSIO.fits_modify_comment(hdu, key, com, Ref{Status}(0)))
    return hdu
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
    print(io, hdu.num)
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

function Base.isequal(::FitsLogic, s1::AbstractString, s2::AbstractString)
    i1, last1 = firstindex(s1), lastindex(s1)
    i2, last2 = firstindex(s2), lastindex(s2)
    @inbounds while (i1 ≤ last1)&(i2 ≤ last2)
        uppercase(s1[i1]) == uppercase(s2[i2]) || return false
        i1 = nextind(s1, i1)
        i2 = nextind(s2, i2)
    end
    @inbounds while i1 ≤ last1
        isspace(s1[i1]) || return false
        i1 = nextind(s1, i1)
    end
    @inbounds while i2 ≤ last2
        isspace(s1[i2]) || return false
        i2 = nextind(s2, i2)
    end
    return true
end
Base.isequal(::FitsLogic, x) = y -> isequal(FitsLogic(), x, y)

"""
    nameof(hdu::FitsHDU) -> str

yields the name of the FITS header data unit `hdu`. The name is the value of
the first card of `"EXTNAME"` or `"HDUNAME"` which exists and has a string
value. Otherwise, the name is that of the FITS extension of `hdu`, that is
`"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"` depending on whether `hdu` is an
image, an ASCII table, a binary table, or anything else.

"""
function Base.nameof(hdu::FitsHDU)
    (str = hdu.hduname) === nothing || return str
    (str = hdu.extname) === nothing || return str
    return hdu.xtension
end

# Compare 2 HDU names, `nothing` can never be considered as a success.
same_name(s1::Nothing, s2::Nothing) = false
same_name(s1::AbstractString, s2::Nothing) = false
same_name(s1::Nothing, s2::AbstractString) = false
same_name(s1::AbstractString, s2::AbstractString) = isequal(FitsLogic(), s1, s2)

"""
    EasyFITS.is_named(hdu, str) -> bool

yields whether string `str` is equal (in the FITS sense) to the extension of
the FITS header data unit `hdu`, or to the value of one of its `"EXTNAME"` or
`"HDUNAME"` keywords. These are respectively given by `hdu.xtension`,
`hdu.extname`, or `hdu.hduname`.

This method is used as a predicate for the search methods `findfirst`,
`findlast`, `findnext`, and `findprev`.

The extension `hdu.xtension` is `"IMAGE"`, `"TABLE"`, `"BINTABLE"`, or `"ANY"`
depending on whether `hdu` is an image, an ASCII table, a binary table, or
anything else.

"""
function is_named(hdu::FitsHDU, str::AbstractString)
    # Since a match only fails if no matching name is found, the order of the
    # tests is irrelevant. We therefore start with the costless ones.
    same_name(hdu.xtension, str) && return true
    same_name(hdu.hduname, str) && return true
    same_name(hdu.extname, str) && return true
    return false
end
is_named(str::AbstractString) = x -> is_named(x, str)

for func in (:findfirst, :findlast)
    @eval begin
        Base.$func(str::AbstractString, hdu::FitsHDU) = $func(str, hdu.io)
        Base.$func(f::Function, hdu::FitsHDU) = $func(f, hdu.io)
        Base.$func(str::AbstractString, io::FitsIO) = $func(is_named(str), io)
    end
end

for func in (:findnext, :findprev)
    @eval begin
        Base.$func(str::AbstractString, hdu::FitsHDU) = $func(str, hdu.io, hdu)
        Base.$func(f::Function, hdu::FitsHDU) = $func(f, hdu.io, hdu)
        Base.$func(str::AbstractString, io::FitsIO, hdu::FitsHDU) =
            $func(is_named(str), io, hdu)
    end
end

Base.findfirst(f::Function, io::FitsIO) = find(f, io, firstindex(io):lastindex(io))
Base.findlast(f::Function, io::FitsIO) = find(f, io, lastindex(io):-1:firstindex(io))
Base.findnext(f::Function, io::FitsIO, hdu::FitsHDU) = find(f, io, hdu.num+1:lastindex(io))
Base.findprev(f::Function, io::FitsIO, hdu::FitsHDU) = find(f, io, hdu.num-1:-1:firstindex(io))

function find(f::Function, io::FitsIO, r::OrdinalRange{<:Integer})
    for num in r
        hdu = io[num]
        f(hdu) && return hdu
    end
    return nothing
end
