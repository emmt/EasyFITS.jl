#------------------------------------------------------------------------------
# FITS TABLES PROPERTIES

Base.propertynames(::FitsTableHDU) = (:nrows, :rows, :first_row, :last_row,
                                      :ncols, :columns, :first_column, :last_column,
                                      :column_name, :column_names, :column_number, :column_units,
                                      :data_size, :data_ndims, :data_axes,
                                      :extname, :hduname, :file, :num, :type, :xtension)

Base.getproperty(hdu::FitsTableHDU, ::Val{:data_ndims}) = 2
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_size}) = (get_nrows(hdu),
                                                          get_ncols(hdu))
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_axes}) = (Base.OneTo(get_nrows(hdu)),
                                                          Base.OneTo(get_ncols(hdu)))

Base.getproperty(hdu::FitsTableHDU, ::Val{:rows}) = Base.OneTo(get_nrows(hdu))
Base.getproperty(hdu::FitsTableHDU, ::Val{:first_row}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_row}) = get_nrows(hdu)
Base.getproperty(hdu::FitsTableHDU, ::Val{:nrows}) = get_nrows(hdu)
function get_nrows(f::Union{FitsFile,FitsTableHDU})
    nrows = Ref{Clonglong}()
    check(CFITSIO.fits_get_num_rowsll(f, nrows, Ref{Status}(0)))
    return as(Int, nrows[])
end

Base.getproperty(hdu::FitsTableHDU, ::Val{:columns}) = Base.OneTo(get_ncols(hdu))
Base.getproperty(hdu::FitsTableHDU, ::Val{:first_column}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_column}) = get_ncols(hdu)
Base.getproperty(hdu::FitsTableHDU, ::Val{:ncols}) = get_ncols(hdu)
function get_ncols(f::Union{FitsFile,FitsTableHDU})
    ncols = Ref{Cint}()
    check(CFITSIO.fits_get_num_cols(f, ncols, Ref{Status}(0)))
    return as(Int, ncols[])
end

struct UnitsGetter
    hdu::FitsTableHDU
end
(obj::UnitsGetter)(col::Column; case::Bool = false) = get_units(obj.hdu, col; case)
Base.getproperty(hdu::FitsTableHDU, ::Val{:column_units}) = UnitsGetter(hdu)

struct ColumnNumberGetter
    hdu::FitsTableHDU
end
(obj::ColumnNumberGetter)(col::Column; case::Bool = false) = get_colnum(obj.hdu, col; case)
Base.getproperty(hdu::FitsTableHDU, ::Val{:column_number}) = ColumnNumberGetter(hdu)

struct ColumnNameGetter
    hdu::FitsTableHDU
end
(obj::ColumnNameGetter)(col::Column; case::Bool = false) = get_colname(obj.hdu, col; case)[1]
Base.getproperty(hdu::FitsTableHDU, ::Val{:column_name}) = ColumnNameGetter(hdu)

"""
    names(hdu::FitsTableHDU)
    hdu.column_names

yield a vector of the column names of the FITS table extension in `hdu`.

"""
Base.names(hdu::FitsTableHDU) = get_colnames(hdu)
Base.getproperty(hdu::FitsTableHDU, ::Val{:column_names}) = get_colnames(hdu)
function get_colnames(hdu::FitsTableHDU)
    ncols = hdu.ncols
    names = Array{String}(undef, ncols)
    for col in 1:ncols
        names[col] = hdu["TTYPE$col"].string
    end
    return names
end

#------------------------------------------------------------------------------
# READING FITS TABLES

columns_to_read(hdu::FitsTableHDU, cols::Columns) = cols
columns_to_read(hdu::FitsTableHDU, cols::Colon) = hdu.columns

rows_to_read(hdu::FitsTableHDU, rows::Rows) = rows
rows_to_read(hdu::FitsTableHDU, rows::Colon) = hdu.rows

first_row_to_read(hdu::FitsTableHDU, rows::Rows) = as(Int, first(rows))
first_row_to_read(hdu::FitsTableHDU, rows::Colon) = hdu.first_row

last_row_to_read(hdu::FitsTableHDU, rows::Rows) = as(Int, last(rows))
last_row_to_read(hdu::FitsTableHDU, rows::Colon) = hdu.last_row

"""
    EasyFITS.get_units(hdu, col; case=false) -> str
    hdu.column_units(col; case=false) -> str

yield the units for the column `col` of FITS table `hdu`. Keyword `case`
specifies whether the case of letters matters when `col` is a (symbolic) name.
The result is always a string, possibly empty.

"""
function get_units(hdu::FitsTableHDU, col::Column; case::Bool = false)
    num = get_colnum(hdu, col; case)
    card = get(hdu, "TUNIT$num", nothing)
    return (card === nothing || card.type != FITS_STRING ? "" : card.string)
end

"""
    EasyFITS.get_colnum(hdu::FitsTableHDU, col; case=false) -> num
    hdu.column_number(col; case=false) -> num

yields the column number of column matching `col` (a string, a symbol, or an
integer) in FITS table extension of `hdu`. Keyword `case` specifies whether the
case of letters matters if `col` is a string or a symbol.

"""
@inline function get_colnum(hdu::FitsTableHDU, col::Integer; case::Bool = false)
    @boundscheck checkbounds(hdu, col)
    return as(Int, col)
end

function get_colnum(hdu::FitsTableHDU, col::ColumnName; case::Bool = false)
    colnum = Ref{Cint}()
    check(CFITSIO.fits_get_colnum(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                  col, colnum, Ref{Status}(0)))
    return as(Int, colnum[])
end

Base.checkbounds(hdu::FitsTableHDU, col::Integer) =
    col ∈ hdu.columns || bad_argument("out of range column index")

"""
    EasyFITS.get_colname(hdu::FitsTableHDU, col; case=false) -> (str, num)
    hdu.column_name(col; case=false) -> str

yields the column name and number of column matching `col` (a string, a symbol,
or an integer) in FITS table extension of `hdu`. Keyword `case` specifies
whether the case of letters matters if `col` is a string or a symbol. If `case`
is false, the column name `str` is converted to upper-case letters.

"""
function get_colname(hdu::FitsTableHDU, col::Integer; case::Bool = false)
    checkbounds(hdu, col)
    card = get(hdu, "TTYPE$col", nothing) # FIXME: optimize?
    (card === nothing || card.type != FITS_STRING) && return ("COL#$col", as(Int, col))
    return (case ? card.string : uppercase(card.string), as(Int, col))
end

function get_colname(hdu::FitsTableHDU, col::ColumnName; case::Bool = false)
    colnum = Ref{Cint}()
    buffer = Vector{UInt8}(undef, CFITSIO.FLEN_VALUE)
    check(CFITSIO.fits_get_colname(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                   col, pointer(buffer), colnum, Ref{Status}(0)))
    colname = to_string!(buffer)
    return (case ? colname : uppercase(colname), as(Int, colnum[]))
end

for func in (:get_coltype, :get_eqcoltype)
    local cfunc = Symbol("fits_",func,"ll")
    @eval function $func(f::Union{FitsFile,FitsTableHDU}, col::ColumnName;
                         case::Bool = false)
        return $func(f, get_colnum(f, col; case))
    end
    @eval function $func(f::Union{FitsFile,FitsTableHDU}, col::Integer;
                         case::Bool = false)
        type = Ref{Cint}()
        repeat = Ref{Clonglong}()
        width = Ref{Clonglong}()
        check(CFITSIO.$cfunc(f, col, type, repeat, width, Ref{Status}(0)))
        return (type[], as(Int, repeat[]), as(Int, width[]))
    end
end

read_tdim(f::Union{FitsFile,FitsTableHDU}, col::ColumnName; case::Bool=false) =
    read_tdim(f, get_colnum(f, col; case))

function read_tdim(f::Union{FitsFile,FitsTableHDU}, col::Integer)
    1 ≤ col ≤ 9999 || throw(FitsError(CFITSIO.BAD_COL_NUM))
    ndims = Ref{Cint}()
    dims = Vector{Clonglong}(undef, 5)
    again = true # recall function?
    while again
        check(CFITSIO.fits_read_tdimll(
            f, col, length(dims), ndims, dims, Ref{Status}(0)))
        again = ndims[] > length(dims)
        ndims[] != length(dims) && resize!(dims, ndims[])
    end
    return dims
end

"""
    read([R=Array,] hdu::FitsTableHDU, col[, rows]; kwds...) -> vals::R

reads column `col` of the FITS table extension in `hdu` and returns its values
and, possibly, its units.

The column `col` may be specified by its name or by its number. If `col` is a
string or a symbol, keyword `case` specifies whether the case of letters
matters (`case = false` by default).

Optional leading argument `R` is to specify the type of the result which can be
an array type, to only retrieve the column values, or a tuple of array and string
types, to retrieve the column values and their units.

Optional argument `rows` is to specify which rows to read. It can be an integer
to read a single row, a unit range of integers to read these rows, or a colon
`:` to read all rows (the default). Use `hdu.first_row` and `hdu.last_row`
to retrieve the first and last row numbers.

See [`read!(::DenseArray,::FitsTableHDU,::ColumName)`](@ref) for the other
possible keywords.

"""
function read(hdu::FitsTableHDU, col::Column, rows::Rows = Colon(); kwds...)
    return read(Array, hdu, col, rows; kwds...)
end

# Read column values and units.
function read(::Type{Tuple{A,String}}, hdu::FitsTableHDU, col::Column,
              rows::Rows = Colon(); case::Bool = false, kwds...) where {A<:Array}
    num = get_colnum(hdu, col; case)
    return (read(A, hdu, num, rows; kwds...), get_units(hdu, num))
end

# Convert column name to number.
function read(::Type{A}, hdu::FitsTableHDU, col::ColumnName,
              rows::Rows = Colon(); case::Bool = false, kwds...) where {A<:Array}
    return read(A, hdu, get_colnum(hdu, col; case), rows; kwds...)
end

function read(::Type{A}, hdu::FitsTableHDU, col::Integer, rows::Rows = Colon();
              kwds...) where {A<:Array}
    # Get equivalent type stored by the column.
    type, repeat, width = get_eqcoltype(hdu, col)
    type > zero(type) || error("reading variable-length array is not implemented")

    # Get dimension of array to read. Only non-string can be shrinked.
    dims = size_to_read(hdu, col, rows, abs(type) != CFITSIO.TSTRING)

    # Allocate array for storing the result, read the data, and eventually
    # convert its element type.
    return read!(args_for_read!(type_from_code(type), dims, A)...,
                 hdu, col; first = first_row_to_read(hdu, rows), kwds...)
end

"""
    EasyFITS.args_for_read!(S::Type, dims, A::Type) -> (arr, T)

given the element type `S` and the dimensions `dims` of the column data in the
file, and the array type `A` requested by the caller, this private function
yields an array `arr` to store the column values with `read!` and the type `T`
of the elements of the array to return to the caller.

If `eltype(arr) == T` holds, `arr` can be directly returned to the caller;
otherwise, a conversion by [`convert_eltype`](@ref) is needed.

If `S` is `String`, then:

- if `A` does not specify an element type, `T` is returned as `String`;
  otherwise, `T` is returned as `eltype(A)` which must be `UInt8` or `String`;

- the number of dimensions `length(dims)` includes the character index and data
  is read as bytes, so `eltype(arr) == UInt8` holds and, if `A` specifies a
  number of dimensions `N`, then:

      N == (T == String ? length(dims) - 1 : length(dims))

  must hold.

Otherwise `S` is a numeric type, and

- if `A` does not specify an element type, `T` is returned as `S`;
  otherwise, `T` is returned as `eltype(A)` and `T <: Real` must hold;

- if `A` specifies a number of dimensions, say `N`, then:

      N == length(dims)

  must hold.

"""
args_for_read!(::Type{S}, dims::Vector{<:Integer}, ::Type{<:Array}) where {S} =
    # Caller specifies no element type and no number of dimensions.
    new_array(eltype_to_read(S), dims), S
args_for_read!(::Type{S}, dims::Vector{<:Integer}, ::Type{<:Array{<:Any,N}}) where {S,N} =
    # Caller specifies the number of dimensions but no element type.
    new_array(eltype_to_read(S), ndims_to_read(dims, Array{S,N}), dims), S
args_for_read!(::Type{S}, dims::Vector{<:Integer}, ::Type{<:Array{T}}) where {S,T} =
    # Caller specifies the element type but no number of dimensions.
    new_array(eltype_to_read(S, T), dims), T
args_for_read!(::Type{S}, dims::Vector{<:Integer}, ::Type{<:Array{T,N}}) where {S,T,N} =
    # Caller specifies the element type and the number of dimensions.
    new_array(eltype_to_read(S, T), ndims_to_read(dims, Array{T,N}), dims), T

"""
    EasyFITS.ndims_to_read(dims, Array{T,N}) -> Val(length(dims))

yields the number of dimensions of the array to read in a column and checks the
number of dimensions of the expected result of type `Array{T,N}` after possible
conversion. This method is type-stable as it yields the following predictable
result:

    Val(N+1)  if T <: String
    Val(N)    otherwise

"""
ndims_to_read(dims::Vector{<:Integer}, ::Type{Array{String,N}}) where {N} =
    (n = length(dims)) == N+1 ? Val(N+1) : bad_column_ndims(N, n-1)

ndims_to_read(dims::Vector{<:Integer}, ::Type{Array{T,N}}) where {T,N} =
    (n = length(dims)) == N ? Val(N) : bad_column_ndims(N, n)

@noinline bad_column_ndims(m::Int, n::Int) = throw(DimensionMismatch(
    "bad number of dimensions for column, got $m, should be $n"))

"""
    EasyFITS.eltype_to_read(S) -> R
    EasyFITS.eltype_to_read(S, T) -> R

yield the type `R` of the elements to read for a column with elements of type
`S`. Optional argument `T` is the element type requested by the caller. An
exception is thrown if `T` is invalid.

"""
eltype_to_read(::Type{String}) = UInt8
eltype_to_read(::Type{S}) where {S<:NumericTypes} = S

# Special case of strings:
eltype_to_read(::Type{String}, ::Type{T}) where {T<:Union{UInt8,String}} = UInt8

# NOTE: Some conversions like float->integer are not allowed.
const AllowedTypesFromComplex = Complex{<:Union{AbstractFloat,Rational}}
const AllowedTypesFromFloat   = Union{AbstractFloat,Rational,AllowedTypesFromComplex}
const AllowedTypesFromInteger = Union{Integer,AllowedTypesFromFloat}

# Conversions handled by CFITSIO:
eltype_to_read(::Type{S}, ::Type{T}) where {S<:IntegerTypes,T<:NumericTypes} = T
eltype_to_read(::Type{S}, ::Type{T}) where {S<:FloatTypes,T<:Union{FloatTypes,ComplexTypes}} = T
eltype_to_read(::Type{S}, ::Type{T}) where {S<:ComplexTypes,T<:ComplexTypes} = T
# Conversions not handled by CFITSIO:
eltype_to_read(::Type{S}, ::Type{T}) where {S<:IntegerTypes,T<:AllowedTypesFromInteger} = S
eltype_to_read(::Type{S}, ::Type{T}) where {S<:FloatTypes,  T<:AllowedTypesFromFloat  } = S
eltype_to_read(::Type{S}, ::Type{T}) where {S<:ComplexTypes,T<:AllowedTypesFromComplex} = S

# Error catchers:
@noinline eltype_to_read(::Type{S}) where {S} = throw(ArgumentError(
    "reading column values of type `$S` is not supported"))
@noinline eltype_to_read(::Type{String}, ::Type{T}) where {T} = throw(ArgumentError(
    "reading column of strings as values of type `$T` is not supported"))
@noinline eltype_to_read(::Type{S}, ::Type{T}) where {S,T} = throw(ArgumentError(
    "reading column values of type `$S` as values of type `$T` is not supported"))

# Yields size of column data as read. If reading numbers and cell has size
# (1,), an empty size is assumed.
function size_to_read(hdu::FitsTableHDU, col::Integer, rows::R,
                      shrink::Bool) where {R<:Rows}
    all_rows = hdu.rows
    if R <: Colon
        nrows = length(all_rows)
    else
        nrows = length(rows)
        nrows == 0 || ((first(rows) ≥ first(all_rows)) & (last(rows) ≤ last(all_rows))) ||
            bad_argument("out of bounds row(s) to read")
    end
    dims = read_tdim(hdu, col)
    shrink &= (length(dims) == 1 && first(dims) == 1) # only shrink if cells are scalars
    if R <: Integer
        # A scalar row is to be read, there is no row index dimension.
        shrink && empty!(dims)
    else
        # Set the last dimension to be the row index.
        if shrink
            dims[firstindex(dims)] = nrows
        else
            push!(dims, nrows)
        end
    end
    return dims
end

"""
    read!(wrk::AbstractArray, T::Type, hdu::FitsTableHDU, col) -> arr

overwrites the elements of workspace array `wrk` with values of the column
`col` read in the FITS table extension `hdu` and returns it with elements of
type `T`. If `eltype(wrk) <: T` holds, `wrk` is returned as `arr`; otherwise
`arr` is a freshly allocated array.

If `T` is `String` and `eltype(wrk)` is `UInt8` the result is an array of
strings built from the bytes along the first dimension of `wrk` for all other
dimensions of `wrk` and with trailing spaces removed.

"""
function read!(wrk::AbstractArray, ::Type{T}, hdu::FitsTableHDU,
               col::Column; kwds...) where {T}
    return convert_eltype(T, read!(wrk, hdu, col; kwds...))
end

convert_eltype(::Type{T}, A::AbstractArray{<:T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = copyto!(similar(A, T), A)

# Convert an array of bytes to an array of strings. The leading dimension
# indexes the characters of the strings.
#
# NOTE: Trailing spaces are stripped, a single space is however significant as
#       assumed in CFITSIO and FITS standard to distinguish null (empty) and
#       non-null strings.
function convert_eltype(::Type{String}, A::AbstractArray{UInt8,N}) where {N}
    B = Array{String}(undef, size(A)[2:N])
    inds = axes(A)
    I = inds[1]
    J = CartesianIndices(inds[2:N])
    buf = Array{UInt8}(undef, length(I))
    guard = Base.cconvert(Ptr{UInt8}, buf)
    ptr = Base.unsafe_convert(Ptr{UInt8}, guard)
    l = firstindex(B) - 1 # fast linear index in B
    nbads = 0 # number of bad characters
    @inbounds for j in J
        flag = false # any character found before the null?
        k = firstindex(buf) - 1 # index of last written character in buf
        last = k # index of last non-space character in buf
        for i in I
            c = A[i,j]
            is_null(c) && break
            if ! is_ascii(c)
                # FIXME: throw(ArgumentError("invalid character `$c`"))
                nbads += 1
                continue
            end
            flag = true
            buf[k += 1] = c
            if ! is_space(c)
                last = k
            end
        end
        len = last - (firstindex(buf) - 1)
        if flag & (len < 1)
            # See NOTE above.
            len = 1
        end
        B[l += 1] = GC.@preserve guard unsafe_string(ptr, len)
    end
    nbads == 0 || @warn "$nbads non-ASCII characters have been filtered"
    return B
end

# Convert an array of strings into an array of bytes. The leading dimension in
# the result indexes the characters of the strings.
#
# NOTE: Trailing spaces are stripped, a single space is however significant as
#       assumed in CFITSIO and FITS standard to distinguish null (empty) and
#       non-null strings.
function convert_eltype(::Type{UInt8}, A::AbstractArray{<:AbstractString,N};
                        firstdim::Integer = maximum_length(A),
                        fillchar::Char = '\0') where {N}
    fillchar ∈ (' ', '\0') || throw(ArgumentError("invalid fill character `$fillchar`"))
    fillchar = UInt8(fillchar)
    B = Array{UInt8}(undef, firstdim, size(A)...)
    maxlen = 0 # maximum length of strings (without non-significant trailing spaces)
    off = firstindex(B) - 1 # offset index in B
    nbads = 0 # number of bad characters
    @inbounds for str in A
        num = 0 # number of input characters
        len = 0 # length of output string without non-significant trailing spaces
        for c in str
            if is_null(c) | (! is_ascii(c))
                # FIXME: throw(ArgumentError("invalid character `$c`"))
                nbads += 1
                continue
            end
            num += 1
            if !is_space(c)
                len = num
            end
            if num ≤ firstdim
                B[off + num] = c
            end
        end
        if (num > 0)&(len < 1)
            # See NOTE above.
            len = 1
        end
        maxlen = max(maxlen, len)
        for k in len+1:firstdim
            B[off + k] = fillchar
        end
        off += firstdim
    end
    firstdim ≥ maxlen || throw(ArgumentError(
        "Cannot add a string of $maxlen bytes, the column maximum is set to $firstdim"))
    nbads == 0 || @warn "$nbads non-ASCII or null characters have been filtered"
    return B
end

# NOTE: `isascii` only takes character argument.
is_ascii(c::UInt8) = c < 0x80
is_ascii(c::Char) = isascii(c)

# NOTE: We are only interested in ASCII spaces.
is_space(c::UInt8) = c == 0x20
is_space(c::Char) = c == ' '

is_null(c::UInt8) = iszero(c)
is_null(c::Char) = c == '\0'

"""
    read([Dict,] hdu::FitsTableHDU[, cols[, rows]]) -> dict

reads some columns of the FITS table extension in `hdu` as a dictionary indexed
by the column names. The columns to read can be specified by `cols` which may
be a single column name/index, a tuple/range/vector of column names/numbers, or
a colon `:` to read all columns (the default). Column names may be strings or
symbols (not a mixture of these). The rows to read can be specified by `rows`
as a single row index, a unit range of row numbers, or a colon `:` to read all
rows (the default).

Keyword `rename` is to specify a function to change column names. If
unspecified, the colmun names are left unchanged if keyword `case` is true and
converted to uppercase letters otherwise.

Keyword `units` can be used to indicate whether to retrieve the units of the
columns. If `units` is `String`, the values of the dictionary will be 2-tuples
`(data,units)` with `data` the column data and `units` the column units as a
string. Otherwise, if `units` is `nothing` (the default), the values of the
dictionary will just be the columns data.

To avoid the `units` keyword, the following methods are provided:

    read(Dict{String,Array},               hdu[, cols[, rows]])
    read(Dict{String,Tuple{Array,String}}, hdu[, cols[, rows]])

to yield the same result as `read(hdu,...)` with respectively `units=nothing`
and `units=String`.

"""
function read(::Type{Dict}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon();
              units = nothing, kwds...)
    if units === String
        return read(Dict{String,Tuple{Array,String}}, hdu, cols, rows; kwds...)
    elseif units === nothing
        return read(Dict{String,Array}, hdu, cols, rows; kwds...)
    else
        throw(ArgumentError("invalid value for keyword `units`"))
    end
end

function read(::Type{D}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon();
              kwds...) where {D<:Union{Dict{String,<:AbstractArray},
                                       Dict{String,<:Tuple{<:AbstractArray,String}}}}
    return merge!(D(), hdu, cols, rows; kwds...)
end

# The default is to read table columns as a dictionary.
function read(hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon(); kwds...)
    return read(Dict, hdu, cols, rows; kwds...)
end

"""
    read!(dict, hdu::FitsTableHDU[, cols[, rows]]) -> dict

overwrites the contents of the dictionary `dict` with the colum(s) `cols` read
from the FITS table extension in `hdu` and returns the dictionary. Any previous
contents is erased, call `merge!(dict,hdu,...)` to preserve contents.

"""
function read!(dict::AbstractDict, hdu::FitsTableHDU,
               cols::Columns = Colon(), rows::Rows = Colon(); kwds...)
    return merge!(empty!(dict), hdu, cols, rows; kwds...)
end

"""
    merge!(dict, hdu::FitsTableHDU[, cols[, rows]]) -> dict

merges the contents of the dictionary `dict` with the colum(s) `cols` read from
the FITS table extension in `hdu` and returns the dictionary.

"""
function merge!(dict::AbstractDict{K,V}, hdu::FitsTableHDU,
                cols::Columns = Colon(), rows::Rows = Colon();
                case::Bool = false, rename::Function = (case ? identity : uppercase),
                kwds...) where {K<:String,V<:Union{AbstractArray,
                                                   Tuple{AbstractArray,String}}}
    names = hdu.column_names
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col; case)
        key = rename(names[num]) # FIXME: as(K, ...) and relax type of parameter K
        push!(dict, key => read(V, hdu, num, rows; kwds...))
    end
    return dict
end

"""
    read(Vector, hdu::FitsTableHDU[, cols[, rows]]) -> vec::Vector

reads some columns of the FITS table extension in `hdu` as a vector. The
columns to read can be specified by `cols` which may be a single column
name/index, a tuple/range/vector of column names/numbers, or a colon `:` to
read all columns (the default). Column names may be strings or symbols (not a
mixture of these). The rows to read can be specified by `rows` as a single row
index, a unit range of row numbers, or a colon `:` to read all rows (the
default). `V` is the type of the result.

Keyword `units` can be used to indicate whether to retrieve the units of the
columns. If `units` is `String`, the values of the dictionary will be 2-tuples
`(data,units)` with `data` the column data and `units` the column units as a
string. Otherwise, if `units=nothing` (the default), the values of the
dictionary will just be the columns data.

To avoid the `units` keyword and allow more control on the type of the result,
the following 2 methods are provided:

    read(Vector{<:Array},               hdu[, cols[, rows]])
    read(Vector{Tuple{<:Array,String}}, hdu[, cols[, rows]])

to yield the same result as `read(hdu,...)` with respectively `units=nothing`
and `units=String`.

"""
function read(::Type{Vector}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon();
              units = nothing, kwds...)
    if units === String
        return read(Vector{Tuple{Array,String}}, hdu, cols, rows; kwds...)
    elseif units === nothing
        return read(Vector{Array}, hdu, cols, rows; kwds...)
    else
        throw(ArgumentError("invalid value for keyword `units`"))
    end
end

# Read selection of columns with or without their units.
function read(::Type{V}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon(); case::Bool = false,
              kwds...) where {T<:Union{AbstractArray,Tuple{AbstractArray,String}},
                              V<:AbstractVector{T}}
    cols = columns_to_read(hdu, cols)
    vec = V(undef, length(cols))
    i = firstindex(vec) - 1
    for col in cols
        vec[i += 1] = read(T, hdu, col, rows; kwds...)
    end
    return vec
end

"""
    push!(vec, hdu::FitsTableHDU[, cols[, rows]]) -> vec

appends rows `rows` of columns `cols` read from FITS table extension `hdu` to
the vector `vec` and returns it.

"""
function push!(vec::AbstractVector{A},
               hdu::FitsTableHDU, cols::Columns = Colon(), rows::Rows = Colon();
               case::Bool = false, kwds...) where {A<:AbstractArray}
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col; case)
        push!(vec, read(A, hdu, num, rows; kwds...))
    end
    return vec
end

function push!(vec::AbstractVector{<:Tuple{A,String}},
               hdu::FitsTableHDU, cols::Columns = Colon(), rows::Rows = Colon();
               case::Bool = false, kwds...) where {A<:AbstractArray}
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col; case)
        push!(vec, (read(A, hdu, num, rows; kwds...), get_units(hdu, num)))
    end
    return vec
end

"""
    read!(arr::DenseArray, hdu::FitsTableHDU, col) -> arr

overwrites the elements of array `arr` with values of the column `col` of the
FITS table extension in `hdu` and returns `arr`.

The column `col` may be specified by its name or by its number. If `col` is a
string or a symbol, keyword `case` indicates whether the case of letters
matters (`case = false` by default).

Keyword `first` may be specified with the index of the first row to read. By
default, `first = hdu.first_row` and reading starts at the first row of the
table.

Keyword `anynull` may be specified with a reference to a boolean
(`Ref{Bool}()`) to retrieve whether any of the read values is undefined.

Keyword `null` may be specified with a reference to a value of the same type as
the elements of the destination `arr` (`Ref{eltype(arr)}()`) to retrieve the
value of undefined values. Keyword `null` may also be set with an array of
`Bool` of same size as `arr` and which will be set to `true` for undefined
values and to `false` elsewhere.

Output arrays `arr` and `null` must have contiguous elements, in other words,
they must be *dense arrays*.

"""
function read!(arr::DenseArray, hdu::FitsTableHDU, col::ColumnName;
               case::Bool = false, kwds...)
    return read!(arr, hdu, get_colnum(hdu, col; case); kwds...)
end

@noinline function read!(arr::DenseArray{T,N}, hdu::FitsTableHDU, col::Integer;
                         kwds...) where {T,N}
    error("reading column of values of type `$T` is not supported")
end

function read!(arr::DenseArray{T,N}, hdu::FitsTableHDU, col::Integer;
               first::Integer = hdu.first_row,
               null::Union{DenseArray{Bool,N},Number,Nothing} = nothing,
               anynull::AnyNull = nothing) where {T<:NumericTypes,N,
                                                  AnyNull<:Union{Ref{Bool},Nothing}}
    temp = AnyNull <: Ref{Bool} ? Ref(zero(Cint)) : C_NULL
    unsafe_read_col!(hdu, col, first, 1, length(arr), arr, null, temp)
    AnyNull <: Ref{Bool} && (anynull[] = !iszero(temp[]))
    return arr
end

for T in NUMERIC_TYPES

    @eval function unsafe_read_col!(hdu::FitsTableHDU, col::Integer, row::Integer,
                                    elem::Integer, nelem::Integer, arr::DenseArray{$T,N},
                                    null::Null, anynull) where {N, Null<:Union{DenseArray{Bool,N},
                                                                               Number,Nothing}}
        if Null <: AbstractArray
            axes(null) == axes(arr) || error("incompatible array axes")
            if nelem > 0
                GC.@preserve arr null check(CFITSIO.$(cfunc("fits_read_colnull_", T))(
                    hdu, col, row, elem, nelem, cpointer(null), cpointer(arr), anynull, Ref{Status}(0)))
                fix_booleans!(null)
            end
        else
            null = (Null <: Number ? as($T, null) : zero($T))
            if nelem > 0
                GC.@preserve arr check(CFITSIO.$(cfunc("fits_read_col_", T))(
                    hdu, col, row, elem, nelem, null, cpointer(arr), anynull, Ref{Status}(0)))
            end
        end
        return nothing
    end

    @eval function unsafe_write_col(hdu::FitsTableHDU, col::Integer, row::Integer,
                                    elem::Integer, nelem::Integer, arr::DenseArray{$T},
                                    null::Null = nothing) where {Null<:Union{Number,Nothing}}
        if nelem > 0
            if Null <: Nothing
                GC.@preserve arr check(CFITSIO.$(cfunc("fits_write_col_", T))(
                    hdu, col, row, elem, nelem, cpointer(arr), Ref{Status}(0)))
            else
                GC.@preserve arr check(CFITSIO.$(cfunc("fits_write_colnull_", T))(
                    hdu, col, row, elem, nelem, cpointer(arr), Ref{$T}(null), Ref{Status}(0)))
            end
        end
    end
end

#------------------------------------------------------------------------------
# WRITING FITS TABLES

# Private functions to parse column definitions. Possible types for `def` are
# given by the union `ColumnDefinition`.

column_eltype(def::ColumnEltype) = def
column_eltype(def::Tuple{ColumnEltype}) = def[1]
column_eltype(def::Tuple{ColumnEltype,ColumnDims}) = def[1]
column_eltype(def::Tuple{ColumnEltype,ColumnUnits}) = def[1]
column_eltype(def::Tuple{ColumnEltype,ColumnDims,ColumnUnits}) = def[1]
column_eltype(def::Tuple{ColumnEltype,ColumnUnits,ColumnDims}) = def[1]

column_dims(def::ColumnEltype) = (1,)
column_dims(def::Tuple{ColumnEltype}) = (1,)
column_dims(def::Tuple{ColumnEltype,ColumnUnits}) = (1,)
column_dims(def::Tuple{ColumnEltype,ColumnDims}) = Tuple(def[2])
column_dims(def::Tuple{ColumnEltype,ColumnDims,ColumnUnits}) = Tuple(def[2])
column_dims(def::Tuple{ColumnEltype,ColumnUnits,ColumnDims}) = Tuple(def[3])

column_units(def::ColumnEltype) = ""
column_units(def::Tuple{ColumnEltype}) = ""
column_units(def::Tuple{ColumnEltype,ColumnDims}) = ""
column_units(def::Tuple{ColumnEltype,ColumnUnits}) = def[2]
column_units(def::Tuple{ColumnEltype,ColumnUnits,ColumnDims}) = def[2]
column_units(def::Tuple{ColumnEltype,ColumnDims,ColumnUnits}) = def[3]

# Yields TFORM given column definition.
function column_tform(def::ColumnDefinition)
    T = column_eltype(def)
    letter = T isa Char ? T : type_to_letter(T)
    repeat = prod(column_dims(def))
    if repeat > 1
        return string(repeat, letter)
    else
        return string(letter)
    end
end

# Yields TDIM given column definition.
function column_tdim(def::ColumnDefinition)
    dims = column_dims(def)
    if length(dims) > 1
        io = IOBuffer()
        write(io, '(')
        join(io, dims, ',')
        write(io, ')')
        return String(take!(io))
    else
        return ""
    end
end

# Yields TUNIT given column definition.
column_tunit(def::ColumnDefinition) = column_units(def)

# Write TDIM header card.
function write_tdim(f::Union{FitsFile,FitsTableHDU}, colnum::Integer,
                    dims::Tuple{Vararg{Integer}})
    if Clong === Clonglong || maximum(dims) ≤ typemax(Clong)
        write_tdim(f, colnum, convert(Tuple{Vararg{Clong}}, dims))
    else
        write_tdim(f, colnum, convert(Tuple{Vararg{Clonglong}}, dims))
    end
end
function write_tdim(f::Union{FitsFile,FitsTableHDU}, colnum::Integer,
                    dims::AbstractVector{<:Integer})
    if Clong === Clonglong || maximum(dims) ≤ typemax(Clong)
        write_tdim(f, colnum, convert(Vector{Clong}, dims))
    else
        write_tdim(f, colnum, convert(Vector{Clonglong}, dims))
    end
end
for T in (Clong === Clonglong ? (Clong,) : (Clong, Clonglong))
    local cfunc = (T === Clong ? :fits_write_tdim : :fits_write_tdimll)
    @eval begin
        function write_tdim(f::Union{FitsFile,FitsTableHDU}, colnum::Integer,
                            dims::DenseVector{$T})
            check(CFITSIO.$cfunc(f, colnum, length(dims), dims, Ref{Status}(0)))
        end
        function write_tdim(f::Union{FitsFile,FitsTableHDU}, colnum::Integer,
                            dims::NTuple{N,$T}) where {N}
            check(CFITSIO.$cfunc(f, colnum, N, Ref(dims), Ref{Status}(0)))
        end
    end
end

function write_tdim(dest::Union{FitsFile,FitsTableHDU},
                    def::ColumnDefinition, k::Integer)
    dims = column_dims(def)
    if length(dims) > 1
        write_tdim(dest, k, dims)
    end
    nothing
end

"""
    write(file, FitsTableHDU, cols...) -> hdu

creates a new FITS table extension in FITS file `file` with columns defined by
`cols...`. Each column definition is a pair `name => format` where `name` is
the column name while `format` specifies the type of the column values and,
optionally, their units and the size of the column cells. The following
definitions are possible:

    name => type
    name => (type,)
    name => (type, units)
    name => (type, dims)
    name => (type, units, dims)
    name => (type, dims, units)

where `type` is either a Julia type (`Number` or `String`) or a letter (see
table below), `dims` is an integer or a tuple of integers, and `units` is a
string. By default, `dims = (1,)` and `units = ""`.

| Type               | Letter | Remarks          |
|:-------------------|:-------|:-----------------|
| `AbstractString`   | `'A'`  | ASCII string     |
| `Bool`             | `'L'`  |                  |
| `Int8`             | `'S'`  | CFITSIO specific |
| `UInt8`            | `'B'`  |                  |
| `Int16`            | `'I'`  |                  |
| `UInt16`           | `'U'`  | CFITSIO specific |
| `Int32`            | `'J'`  |                  |
| `UInt32`           | `'V'`  | CFITSIO specific |
| `Int64`            | `'K'`  |                  |
| `UInt64`           | `'W'`  | CFITSIO specific |
| `Float32`          | `'E'`  |                  |
| `Float64`          | `'D'`  |                  |
| `Complex{Float32}` | `'C'`  |                  |
| `Complex{Float64}` | `'M'`  |                  |
| `FitsBit`          | `'X'`  |  bits

The returned object can be used to add FITS keywords to the header of the table
and, then, to write column data. Typically:

    hdu = write(file, FitsTableHDU, cols)
    push!(hdu, key1 => val1) # add a first header keyword
    ...                      # add other header keywords
    write(hdu, col1 => arr1) # write a first column
    ...                      # write other columns

where `key1 => val1`, `key2 => val2`, etc. specify header cards, while `col1 =>
arr1`, `col2 => arr2`, etc. specify columns names and associated data. Such a
table may be created in a single call:

    write(file, [key1 => val1, key2 => val2, ...], [col1 => arr1, col2 => arr2, ...])

"""
function write(file::FitsFile, ::Type{FitsTableHDU},
               cols::Pair{<:ColumnName,<:Any}...; kwds...)
    return write(file, FitsTableHDU, cols; kwds...)
end

function write(file::FitsFile, ::Type{FitsTableHDU},
               cols::Union{AbstractVector{<:Pair{<:ColumnName,<:Any}},
                           Tuple{Vararg{Pair{<:ColumnName,<:Any}}},
                           NamedTuple};
               ascii::Bool=false, nrows::Integer=0)

    function parse_column_definition!(ttype::AbstractVector{String},
                                      tform::AbstractVector{String},
                                      tunit::AbstractVector{String},
                                      name::String,
                                      def::ColumnDefinition,
                                      k::Integer)
        # NOTE: see column_tform, column_tunit, and column_tdim
        T = column_eltype(def)
        letter = T isa Char ? T : type_to_letter(T)
        letter === 'A' && def isa DimensionlessColumnDefinition && throw(ArgumentError(
            "dimensions must be specified for column of strings \"$name\""))
        repeat = prod(column_dims(def))
        ttype[k] = name
        tform[k] = repeat > 1 ? string(repeat, letter) : string(letter)
        tunit[k] = column_units(def)
        nothing
    end

    # Error catcher.
    function parse_column_definition!(@nospecialize(ttype::AbstractVector{String}),
                                      @nospecialize(tform::AbstractVector{String}),
                                      @nospecialize(tunit::AbstractVector{String}),
                                      name::String,
                                      @nospecialize(def),
                                      @nospecialize(k::Integer))
        bad_argument("invalid definition for column \"$name\"")
    end

    ascii && error("only creating binary tables is supported")
    tbl = ascii ? CFITSIO.ASCII_TBL : CFITSIO.BINARY_TBL

    # Write column defintions.
    ncols = length(cols)
    ttype = Vector{String}(undef, ncols)
    tform = Vector{String}(undef, ncols)
    tunit = Vector{String}(undef, ncols)
    k = 0
    for (key, def) in pairs_producer(cols)
        parse_column_definition!(ttype, tform, tunit, String(key), def, k += 1)
    end
    check(CFITSIO.fits_create_tbl(file, tbl, nrows, ncols, ttype, tform, tunit,
                                  C_NULL, Ref{Status}(0)))

    # Write cell dimensions if needed.
    k = 0
    for (key, def) in pairs_producer(cols)
        write_tdim(file, def, k += 1)
    end

    # The number of HDUs as returned by fits_get_num_hdus is only incremented
    # after writing data.
    n = position(file)
    if length(file) < n
        setfield!(file, :nhdus, n)
    end
    return FitsTableHDU(BareBuild(), file, n, ascii)
end

# This function yields an iterable object that produces pairs.
pairs_producer(x::AbstractDict) = x
pairs_producer(x::AbstractVector{<:Pair}) = x
pairs_producer(x::Tuple{Vararg{Pair}}) = x
pairs_producer(x::NamedTuple) = pairs(x)
pairs_producer(x::Any) = throw(ArgumentError("argument is not a collection of pairs"))

"""
    write(hdu::FitsTableHDU, cols...;
          first=hdu.first_row, case=false, null=nothing) -> hdu

writes columns `cols...` into the FITS table extension of `hdu`. Columns
`cols...` are specified as pairs like `col => vals` or `col => (vals, units)`
with `col` the column name/number, `vals` an array of column values, and
`units` optional units. Columns `cols...` can also be named tuples, vectors or
tuples of pairs.

For a given column, if `col` is a column name, keyword `case` indicates whether
the case of letters matters (default is `false`) and, if units are specified,
they must match the value of the FITS keyword `"TUNIT\$n"` in the header of
`hdu` with `n` the column number.

Column values are converted as needed and are written starting at the row
specified by `first`. The leading dimensions of `vals` should be the same as
those specified by the FITS keyword `"TDIM\$n"` (with `n` the column number) in
the header of `hdu` and the remaining last dimension, if any, corresponds to
the *row* index of the table.

Keyword `null` may be used to specify the value of undefined elements in `arr`.

The same keywords apply to all columns.

"""
function write(hdu::FitsTableHDU,
               pair::ColumnDataPair;
               case::Bool = false,
               null::Union{Number,AbstractString,Nothing} = nothing,
               first::Integer = hdu.first_row)
    # Private function to retrieve the values and the units of a column.
    data_units(col::Pair{<:Column,<:ColumnData}) = (pair.second, nothing)
    data_units(col::Pair{<:Column,<:Tuple{ColumnData,ColumnUnits}}) = pair.second

    # Get column number.
    col = pair.first
    num = get_colnum(hdu, col; case)

    # Get column equivalent type.
    type, repeat, width = get_eqcoltype(hdu, num)
    type > zero(type) || error("reading variable length array in column `$col` is not implemented")

    # Retrieve column values and, if units are specified, check that they are
    # the same as those already set.
    vals, units = data_units(pair)
    if units !== nothing
        card = get(hdu, "TUNIT$num", nothing)
        (card === nothing ? isempty(units) :
            card.type == FITS_STRING && units == card.string) || bad_argument(
                "invalid units for column `$col`")
    end

    # Check string case.
    if abs(type) == CFITSIO.TSTRING
        eltype(vals) <: Union{AbstractString,UInt8} || bad_argument(
            "columns of strings can only be written as bytes or strings in column `$col`")
    else
        eltype(vals) <: AbstractString && bad_argument(
            "strings can only be written to columns of strings in column `$col`")
    end

    # Check dimensions. The number of rows is adjusted automatically.
    cell_dims = read_tdim(hdu, num)
    prod(cell_dims) == repeat || error(
        "unexpected column repeat count ($repeat) and product of cell dimensions ($(prod(cell_dims))) for column `$col`")
    cell_ndims = length(cell_dims)
    if !(eltype(vals) <: AbstractString) && cell_ndims == 1 && Base.first(cell_dims) == 1
        # For columns of values, not strings, TDIM = (1,) is the same as no TDIM.
        cell_ndims = 0
    end
    vals_dims = size(vals)
    vals_ndims = ndims(vals)
    if eltype(vals) <: AbstractString
        off = 1 # the first dimension is to index characters and does not count here
    else
        off = 0 # all leading dimensions must be identical
    end
    n = cell_ndims - off # number of dimension to compare
    # a single table row | several table rows
    ((vals_ndims == n) | (vals_ndims == n+1)) || throw(DimensionMismatch(
        "bad number of dimensions for column `$col`, got $vals_ndims instead of $n or $(n+1)"))
    for i in 1:n
        vals_dims[i] == cell_dims[i+off] || throw(DimensionMismatch(
            "incompatible dimensions for column `$col`"))
    end

    # Write column values.
    if eltype(vals) <: AbstractString
        null isa Union{AbstractString,Nothing} || bad_argument(
            "invalid `null` type for column `$col`")
        fillchar =
            null === nothing || isempty(null) ? '\0' :
            null == " " ? ' ' : bad_argument("invalid `null` value for column `$col`")
        firstdim = Base.first(cell_dims) # number of character per string as stored in this column
        temp = convert_eltype(UInt8, vals; firstdim, fillchar)
        unsafe_write_col(hdu, num, first, 1, length(temp), temp, nothing)
    else
        null isa AbstractString && bad_argument("invalid `null` type for column `$col`")
        unsafe_write_col(hdu, num, first, 1, length(vals), dense_array(vals), null)
    end
    return hdu
end

write(hdu::FitsTableHDU; kdws...) = hdu # nothing to write
@inline function write(hdu::FitsTableHDU, col::ColumnDataPair,
                       cols::ColumnDataPair...; kwds...)
    write(hdu, col; kwds...)
    return write(hdu, cols...; kwds...)
end

write(hdu::FitsTableHDU, ::Tuple{}; kdws...) = hdu # no more columns to write
@inline function write(hdu::FitsTableHDU, cols::Tuple{Vararg{ColumnDataPair}};
                       kwds...)
    write(hdu, first(cols); kwds...)
    return write(hdu, Base.tail(cols)...; kwds...)
end

function write(hdu::FitsTableHDU, cols::NamedTuple; kwds...)
    for col in pairs(cols)
        write(hdu, col; kwds...)
    end
    return hdu
end

function write(hdu::FitsTableHDU,
               # NOTE: For a vector of input column data, we cannot be very
               # specific for the key and value types of the pairs if we want
               # to allow for a mixture of these types. The "error catcher"
               # version of the `write` method will be eventually called.
               cols::AbstractVector{<:Pair};
               kwds...)
    for col in cols
        write(hdu, col)
    end
    return hdu
end

# Error catcher.
@noinline write(hdu::FitsTableHDU, col::Pair{K,V}, args...; kwds...) where {K,V} =
    throw(ArgumentError("invalid column-data pair of type Pair{$K,$V}"))

"""
    write(file::FitsFile, header, cols; ascii=false) -> file

creates a FITS table extension in `file` with additional keywords given by
`header` and columns `cols...`. The columns `cols` are specified as a
collection of pairs like `key => vals` or `key => (vals, units)` with `key` the
(symbolic) name of the column, `vals` its values, and `units` its optional
units. The collection can be a dictionary, a named tuple, a vector of pairs, or
a tuple of pairs.

"""
function write(file::FitsFile, header::OptionalHeader,
               cols::TableData; ascii::Bool=false)

    @noinline bad_column_data(key, val) = throw(ArgumentError("invalid data for column $key"))

    column_definition(key, val::ColumnData{T,1}) where {T} = (type_to_letter(T), ())
    column_definition(key, val::ColumnData{T,N}) where {T,N} = (type_to_letter(T), size(val)[1:N-1])
    column_definition(key, val::ColumnData{<:AbstractString,1})  = ('A', (maximum_length(val),))
    column_definition(key, val::ColumnData{<:AbstractString,N}) where {N} = ('A', (maximum_length(val),
                                                                                   size(val)[1:N-1]...))
    column_definition(key, val::Tuple{ColumnData,ColumnUnits}) = (column_definition(key, val[1])..., val[2])
    column_definition(key, val::Any) = bad_column_data(key, val)

    column_values(key, val::AbstractArray) = val
    column_values(key, val::Tuple{AbstractArray,String}) = val[1]
    column_values(key, val::Any) = bad_column_data(key, val)

    # Create a new FITS table HDU with column definitions.
    defs = sizehint!(Vector{Pair{String,Any}}(undef, 0), length(cols))
    for (key, val) in pairs_producer(cols)
        push!(defs, String(key) => column_definition(key, val))
    end
    hdu = write(file, FitsTableHDU, defs; ascii)

    # Merge the keywords provided by the caller.
    merge!(hdu, header)

    # Write values of columns.
    for (key, val) in pairs_producer(cols)
        write(hdu, key => column_values(key, val))
    end
    return file
end

#------------------------------------------------------------------------------
