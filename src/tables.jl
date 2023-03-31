#------------------------------------------------------------------------------
# FITS TABLES PROPERTIES

Base.propertynames(::FitsTableHDU) = (:nrows, :ncols, :column_names,
                                      :first_row, :last_row,
                                      :first_column, :last_column,
                                      :data_size, :data_ndims, :data_axes,
                                      :extname, :hduname, :file, :num, :type, :xtension)

Base.getproperty(hdu::FitsTableHDU, ::Val{:data_ndims}) = 2
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_size}) = (get_nrows(hdu),
                                                          get_ncols(hdu))
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_axes}) = (Base.OneTo(get_nrows(hdu)),
                                                          Base.OneTo(get_ncols(hdu)))

Base.getproperty(hdu::FitsTableHDU, ::Val{:first_row}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_row}) = get_nrows(hdu)
Base.getproperty(hdu::FitsTableHDU, ::Val{:nrows}) = get_nrows(hdu)

Base.getproperty(hdu::FitsTableHDU, ::Val{:first_column}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_column}) = get_ncols(hdu)
Base.getproperty(hdu::FitsTableHDU, ::Val{:ncols}) = get_ncols(hdu)

function get_nrows(f::Union{FitsFile,FitsTableHDU})
    nrows = Ref{Clonglong}()
    check(CFITSIO.fits_get_num_rowsll(f, nrows, Ref{Status}(0)))
    return as(Int, nrows[])
end

function get_ncols(f::Union{FitsFile,FitsTableHDU})
    ncols = Ref{Cint}()
    check(CFITSIO.fits_get_num_cols(f, ncols, Ref{Status}(0)))
    return as(Int, ncols[])
end

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

# The different possibilities to specify a column and a list of columns to read.
const Column = Union{ColumnName,Integer}
const Columns = Union{Colon,Integer,
                      AbstractVector{<:Integer},
                      Tuple{Vararg{Integer}},
                      OrdinalRange{<:Integer,<:Integer},
                      AbstractString,
                      AbstractVector{<:AbstractString},
                      Tuple{Vararg{AbstractString}},
                      Symbol,
                      AbstractVector{<:Symbol},
                      Tuple{Vararg{Symbol}}}

# There different possibilities to specify a list of rows to read.
const Rows = Union{Colon,Integer,AbstractUnitRange{<:Integer}}

columns_to_read(hdu::FitsTableHDU, cols::Columns) = cols
columns_to_read(hdu::FitsTableHDU, cols::Colon) = hdu.first_column:hdu.last_column

rows_to_read(hdu::FitsTableHDU, rows::Rows) = rows
rows_to_read(hdu::FitsTableHDU, rows::Colon) =
    first_row_to_read(hdu,rows):last_row_to_read(hdu,rows)

first_row_to_read(hdu::FitsTableHDU, rows::Rows) = as(Int, first(rows))
first_row_to_read(hdu::FitsTableHDU, rows::Colon) = hdu.first_row

last_row_to_read(hdu::FitsTableHDU, rows::Rows) = as(Int, last(rows))
last_row_to_read(hdu::FitsTableHDU, rows::Colon) = hdu.last_row

"""
    EasyFITS.get_units(hdu, col, case=false) -> str

yields the units for the column `col` of FITS table `hdu`. If `case` is true
and `col` is a (symbolic) name, the case of the letters is taken into account.

"""
function get_units(hdu::FitsTableHDU, col::ColumnName, case::Bool = false)
    return get_units(hdu, get_colnum(hdu, col, case), false)
end

function get_units(hdu::FitsTableHDU, col::Integer, case::Bool = false)
    card = get(hdu, "TUNIT$col", nothing)
    if card !== nothing && card.type == FITS_STRING
        return card.string
    else
        return ""
    end
end

"""
    EasyFITS.get_colnum(hdu::FitsTableHDU, col, case=false) -> num

yields the column number of column matching `col` (a string, a symbol, or an
integer) in FITS table extension of `hdu`. Optional argument `case` specify
whether upper-/lower-case matters if `col` is a string or a symbol.

"""
@inline function get_colnum(hdu::FitsTableHDU, col::Integer, case::Bool = false)
    @boundscheck check_colnum(hdu, col)
    return as(Int, col)
end

function get_colnum(hdu::FitsTableHDU, col::ColumnName, case::Bool = false)
    colnum = Ref{Cint}()
    check(CFITSIO.fits_get_colnum(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                  col, colnum, Ref{Status}(0)))
    return as(Int, colnum[])
end

check_colnum(hdu::FitsTableHDU, col::Integer) =
    ((hdu.first_column ≤ col) & (col ≤ hdu.last_column)) || bad_argument("out of range column index")

"""
    EasyFITS.get_colname(hdu::FitsTableHDU, col, case=false) -> (str, num)

yields the column name and number of column matching `col` (a string, a symbol,
or an integer) in FITS table extension of `hdu`. Optional argument `case`
specify whether upper-/lower-case matters if `col` is a string or a symbol. If
`case` is false, the column name `str` is converted to upper-case letters.

"""
function get_colname(hdu::FitsTableHDU, col::Integer, case::Bool = false)
    check_colnum(hdu, col)
    card = get(hdu, "TTYPE$col", nothing) # FIXME: optimize?
    (card === nothing || card.type != FITS_STRING) && return ("COL#$col", as(Int, col))
    return (case ? card.string : uppercase(card.string), as(Int, col))
end

function get_colname(hdu::FitsTableHDU, col::ColumnName, case::Bool = false)
    colnum = Ref{Cint}()
    colname = Vector{UInt8}(undef, CFITSIO.FLEN_VALUE)
    check(CFITSIO.fits_get_colname(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                   col, pointer(colname), colnum, Ref{Status}(0)))
    return (to_string!(colname), as(Int, colnum[]))
end

for func in (:get_coltype, :get_eqcoltype)
    local cfunc = Symbol("fits_",func,"ll")
    @eval function $func(f::Union{FitsFile,FitsTableHDU}, col::ColumnName,
                         case::Bool = false)
        return $func(f, get_colnum(f, col, case))
    end
    @eval function $func(f::Union{FitsFile,FitsTableHDU}, col::Integer)
        type = Ref{Cint}()
        repeat = Ref{Clonglong}()
        width = Ref{Clonglong}()
        check(CFITSIO.$cfunc(f, col, type, repeat, width, Ref{Status}(0)))
        return (type[], as(Int, repeat[]), as(Int, width[]))
    end
end

read_tdim(f::Union{FitsFile,FitsTableHDU}, col::ColumnName, case::Bool=false) =
    read_tdim(f, get_colnum(f, col, case))

function read_tdim(f::Union{FitsFile,FitsTableHDU}, col::Integer)
    1 ≤ col ≤ 999 || error("invalid column number")
    naxis = Ref{Cint}()
    naxes = Vector{Clonglong}(undef, 5)
    again = true # recall function?
    while again
        check(CFITSIO.fits_read_tdimll(f, col, length(naxes), naxis, naxes,
                                       Ref{Status}(0)))
        again = naxis[] > length(naxes)
        if naxis[] != length(naxes)
            resize!(naxes, naxis[])
        end
    end
    return as(Vector{Int}, naxes)
end

"""
    read([Array,] hdu::FitsTableHDU, col[, rows]; kwds...) -> vals::Array

reads column `col` of the FITS table extension in `hdu` and returns its values
as an array `vals`.

The column `col` may be specified by its name or by its number. If `col` is a
string or a symbol, keyword `case` specifies whether uppercase/lowercase
matters (`case = false` by default).

Optional argument `rows` is to specify which rows to read. It can be an integer
to read a single row, a unit range of integers to read these rows, or a colon
`:` to read all rows (the default). Properties `hdu.first_row` and
`hdu.last_row` can be used to retrieve the first and last row numbers.

See [`read!(::DenseArray,::FitsTableHDU,::Integer)`](@ref) for the other
possible keywords.

"""
function read(hdu::FitsTableHDU, col::Column, rows::Rows = Colon(); kwds...)
    return read(Array, hdu, col, rows; kwds...)
end

function read(::Type{A}, hdu::FitsTableHDU, col::ColumnName,
              rows::Rows = Colon(); case::Bool = false, kwds...) where {A<:Array}
    return read(A, hdu, get_colnum(hdu, col, case), rows; kwds...)
end

function read(::Type{Array}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...)
    type, repeat, width = get_eqcoltype(hdu, col)
    if type == CFITSIO.TSTRING
        return read(Array{String}, hdu, col, rows; kwds...)
    else
        dims = size_to_read(hdu, col, rows, true)
        return read!(new_array(type_from_code(type), dims), hdu, col;
                     first = first_row_to_read(hdu, rows), kwds...)
    end
end

function read(::Type{Array{T}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {T<:Number}
    dims = size_to_read(hdu, col, rows, true)
    return read!(new_array(T, dims), hdu, col;
                 first = first_row_to_read(hdu, rows), kwds...)
end

function read(::Type{Array{T,N}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {T<:Number,N}
    dims = size_to_read(hdu, col, rows, true)
    length(dims) == N || error("invalid number of dimensions")
    return read!(new_array(T, Val(N), dims), hdu, col;
                 first = first_row_to_read(hdu, rows), kwds...)
end

# Read array of strings as bytes. FIXME: null and anynull keywords

function read(::Type{Array{String}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...)
    dims = size_to_read(hdu, col, rows, false)
    return bytes_to_strings(read!(new_array(UInt8, dims), hdu, col;
                                  first = first_row_to_read(hdu, rows), kwds...))
end

function read(::Type{Array{String,N}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {N}
    dims = size_to_read(hdu, col, rows, false)
    length(dims) == N+1 || error("invalid number of dimensions")
    return bytes_to_strings(read!(new_array(UInt8, Val(N+1), dims), hdu, col;
                                  first = first_row_to_read(hdu, rows), kwds...))
end

# Yields size of column data as read.
function size_to_read(hdu::FitsTableHDU, col::Integer, rows::Rows, compress::Bool = true)
    first = first_row_to_read(hdu, rows)
    last = last_row_to_read(hdu, rows)
    nrows = max(last + 1 - first, 0)
    if nrows > 0 && (first < hdu.first_row || last > hdu.last_row)
        bad_argument("out of bounds first row to read")
    end
    dims = read_tdim(hdu, col)
    if compress && length(dims) == 1 && Base.first(dims) == 1
        dims[firstindex(dims)] = nrows
    elseif !(rows isa Integer)
        push!(dims, nrows)
    end
    return dims
end

function bytes_to_strings(A::AbstractArray{UInt8,N}) where {N}
    I = axes(A)[1]
    J = CartesianIndices(axes(A)[2:N])
    B = Array{String}(undef, size(J))
    buf = Array{UInt8}(undef, length(I))
    space = UInt8(' ')
    k = 0
    @inbounds for j in J
        len = 0
        l = 0
        for i in I
            c = A[i,j]
            iszero(c) && break
            l += 1
            buf[l] = c
            if c != space
                len = l
            end
        end
        k += 1
        B[k] = GC.@preserve buf unsafe_string(pointer(buf), len)
    end
    return B
end

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
string. Otherwise, if `units=nothing` (the default), the values of the
dictionary will just be the columns data.

To avoid the `units` keyword, the following 2 calls are provided:

   read(Dict{String,Array},               hdu[, cols[, rows]])
   read(Dict{String,Tuple{Array,String}}, hdu[, cols[, rows]])

to respectively yields the same result as `read(hdu,...)` with
`units=nothing` and as with `units=String`.

"""
function read(hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon(); kwds...)
    return read(Dict, hdu, cols, rows; kwds...)
end

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
function merge!(dict::AbstractDict{String,<:Array}, hdu::FitsTableHDU,
                cols::Columns = Colon(), rows::Rows = Colon();
                case::Bool = false, rename::Function = (case ? identity : uppercase), kwds...)
    names = hdu.column_names
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col, case)
        key = rename(names[num])
        push!(dict, key => read(Array, hdu, num, rows; kwds...))
    end
    return dict
end

function merge!(dict::AbstractDict{String,Tuple{<:Array,String}}, hdu::FitsTableHDU,
                cols::Columns = Colon(), rows::Rows = Colon();
                case::Bool = false, rename::Function = (case ? identity : uppercase), kwds...)
    names = hdu.column_names
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col, case)
        key = rename(names[num])
        push!(dict, key => (read(Array, hdu, num, rows; kwds...),
                            get_units(hdu, num)))
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
the following 2 calls are provided:

   read(Vector{<:Array},               hdu[, cols[, rows]])
   read(Vector{Tuple{<:Array,String}}, hdu[, cols[, rows]])

to respectively yields the same result as `read(hdu,...)` with
`units=nothing` and as with `units=String`.

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

# Read selection of columns without their units.
function read(::Type{V}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon(); case::Bool = false,
              kwds...) where {A<:AbstractArray, V<:AbstractVector{A}}
    cols = columns_to_read(hdu, cols)
    vec = V(undef, length(cols))
    i = firstindex(vec) - 1
    for col in cols
        num = get_colnum(hdu, col, case)
        vec[i += 1] = read(A, hdu, num, rows; kwds...)
    end
    return vec
end

# Read selection of columns with their units.
function read(::Type{V}, hdu::FitsTableHDU,
              cols::Columns = Colon(), rows::Rows = Colon(); case::Bool = false,
              kwds...) where {A<:AbstractArray, V<:AbstractVector{Tuple{A,String}}}
    cols = columns_to_read(hdu, cols)
    vec = V(undef, length(cols))
    i = firstindex(vec) - 1
    for col in cols
        num = get_colnum(hdu, col, case)
        vec[i += 1] = (read(A, hdu, num, rows; kwds...), get_units(hdu, num))
    end
    return vec
end

function push!(vec::Type{V}, hdu::FitsTableHDU,
               cols::Columns = Colon(), rows::Rows = Colon(); case::Bool = false,
               kwds...) where {A<:AbstractArray, V<:AbstractVector{A}}
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col, case)
        push!(vec, read(A, hdu, num, rows; kwds...))
    end
    return vec
end

function push!(vec::Type{V}, hdu::FitsTableHDU,
               cols::Columns = Colon(), rows::Rows = Colon(); case::Bool = false,
               kwds...) where {A<:AbstractArray, V<:AbstractVector{Tuple{A,String}}}
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col, case)
        push!(vec, (read(A, hdu, num, rows; kwds...), get_units(hdu, num)))
    end
    return vec
end

"""
    read!(arr::DenseArray, hdu::FitsTableHDU, col) -> arr

overwrites the elements of array `arr` with values of the column `col` of the
FITS table extension in `hdu` and returns `arr`.

The column `col` may be specified by its name or by its number. If `col` is a
string or a symbol, keyword `case` indicates whether uppercase/lowercase
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
function read!(arr::DenseArray, hdu::FitsTableHDU,
               col::AbstractString; case::Bool = false, kwds...)
    return read!(arr, hdu, get_colnum(hdu, col, case); kwds...)
end

function read!(arr::DenseArray{String,N}, hdu::FitsTableHDU, col::Integer;
               kwds...) where {N}
    error("reading column of strings not yet implemented")
end

function read!(arr::DenseArray{T,N}, hdu::FitsTableHDU, col::Integer;
               null::Union{DenseArray{Bool,N},Number,Nothing} = nothing,
               anynull::Union{Ref{Bool},Nothing} = nothing,
               first::Integer = hdu.first_row) where {T<:Number,N}
    if null === nothing
        _null = zero(T)
    elseif null isa Number
        _null = as(T, null)
    else
        _null = null
    end
    if anynull isa Ref{Bool}
        _anynull = Ref(zero(Cint))
    else
        _anynull = C_NULL
    end
    read_col!(hdu, col, first, 1, arr, _null, _anynull)
    if anynull isa Ref{Bool}
        anynull[] = !iszero(_anynull[])
    end
    return arr
end

for T in (Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64,
          Float32, Float64, Complex{Float32}, Complex{Float64}, Bool)

    @eval function read_col!(hdu::FitsTableHDU,
                             col::Integer, row::Integer, elem::Integer,
                             arr::DenseArray{$T},
                             null::$T, anynull)
        len = length(arr)
        if len > 0
            GC.@preserve arr begin
                check(CFITSIO.$(cfunc("fits_read_col_", T))(
                    hdu, col, row, elem, len, null, cpointer(arr), anynull,
                    Ref{Status}(0)))
            end
        end
        nothing
    end

    @eval function read_col!(hdu::FitsTableHDU,
                             col::Integer, row::Integer, elem::Integer,
                             arr::DenseArray{$T,N},
                             null::DenseArray{Bool,N}, anynull) where {N}
        len = length(arr)
        axes(null) == axes(arr) || error("incompatible array axes")
        if len > 0
            GC.@preserve arr null begin
                check(CFITSIO.$(cfunc("fits_read_colnull_", T))(
                    hdu, col, row, elem, len, cpointer(null), cpointer(arr), anynull,
                    Ref{Status}(0)))
            end
            fix_booleans!(null)
        end
        nothing
    end

    @eval function write_col(hdu, col::Integer, row::Integer, elem::Integer,
                             arr::DenseArray{$T}, null::Nothing = nothing)
        len = length(arr)
        if len > 0
            GC.@preserve arr begin
                check(CFITSIO.$(cfunc("fits_write_col_", T))(
                    hdu, col, row, elem, len, cpointer(arr), Ref{Status}(0)))
            end
        end
        nothing
    end

    @eval function write_col(hdu, col::Integer, row::Integer, elem::Integer,
                             arr::DenseArray{$T,N}, null::Number) where {N}

        len = length(arr)
        if len > 0
            GC.@preserve arr begin
                check(CFITSIO.$(cfunc("fits_write_colnull_", T))(
                    hdu, col, row, elem, len, cpointer(arr), Ref{T}(null),
                    Ref{Status}(0)))
            end
        end
        nothing
    end

end

#------------------------------------------------------------------------------
# WRITING FITS TABLES

# Private definitions and functions to decode column definitions.
const ColumnType = Union{Type,Char}
const ColumnDims = Union{Integer,Tuple{Vararg{Integer}}}
const ColumnUnits = AbstractString
const ColumnDefinition = Union{ColumnType,
                               Tuple{ColumnType},
                               Tuple{ColumnType,ColumnDims},
                               Tuple{ColumnType,ColumnUnits},
                               Tuple{ColumnType,ColumnUnits,ColumnDims}}
#
column_type(def::ColumnType) = def
column_type(def::Tuple{ColumnType}) = def[1]
column_type(def::Tuple{ColumnType,ColumnDims}) = def[1]
column_type(def::Tuple{ColumnType,ColumnUnits}) = def[1]
column_type(def::Tuple{ColumnType,ColumnDims,ColumnUnits}) = def[1]
column_type(def::Tuple{ColumnType,ColumnUnits,ColumnDims}) = def[1]
#
column_dims(def::ColumnType) = (1,)
column_dims(def::Tuple{ColumnType}) = (1,)
column_dims(def::Tuple{ColumnType,ColumnUnits}) = (1,)
column_dims(def::Tuple{ColumnType,ColumnDims}) = Tuple(def[2])
column_dims(def::Tuple{ColumnType,ColumnDims,ColumnUnits}) = Tuple(def[2])
column_dims(def::Tuple{ColumnType,ColumnUnits,ColumnDims}) = Tuple(def[3])
#
column_units(def::ColumnType) = ""
column_units(def::Tuple{ColumnType}) = ""
column_units(def::Tuple{ColumnType,ColumnDims}) = ""
column_units(def::Tuple{ColumnType,ColumnUnits}) = def[2]
column_units(def::Tuple{ColumnType,ColumnUnits,ColumnDims}) = def[2]
column_units(def::Tuple{ColumnType,ColumnDims,ColumnUnits}) = def[3]

# Yields TFORM given column definition.
function column_tform(def::ColumnDefinition)
    T = column_type(def)
    if T isa Char
        letter = T
    else
        letter = type_to_letter(T)
    end
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

"""
    write(file, FitsTableHDU, cols) -> hdu

creates a new FITS table extension in FITS file `file` with columns defined by
`cols`. Each column definition is a pair `name => format` where `name` is the
column name while `format` specifies the type of the column values and,
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

Such a table may be created in a single call:

    write(file, [key1 => val1, key2 => val2, ...], [col1 => arr1, col2 => arr2, ...])

where `key1 => val1`, `key2 => val2`, etc. specify header cards, while `col1 =>
arr1`, `col2 => arr2`, etc. specify columns names and associated data.

"""
function write(file::FitsFile, ::Type{FitsTableHDU},
               cols::AbstractVector{Pair{K,V}};
               ascii::Bool=false, nrows::Integer=0) where {K<:AbstractString,V}
    ascii && error("only creating binary tables is supported")
    tbl = ascii ? CFITSIO.ASCII_TBL : CFITSIO.BINARY_TBL
    ncols = length(cols)
    ttype = Vector{String}(undef, ncols)
    tform = Vector{String}(undef, ncols)
    tunit = Vector{String}(undef, ncols)
    k = 0
    for (name, def) in cols
        k += 1
        def isa ColumnDefinition || bad_argument("invalid definition for column \"$name\"")
        ttype[k] = name
        tform[k] = column_tform(def)
        tunit[k] = column_tunit(def)
    end
    check(CFITSIO.fits_create_tbl(file, tbl, nrows, ncols, ttype, tform, tunit,
                                  C_NULL, Ref{Status}(0)))
    k = 0
    for (name, def) in cols
        k += 1
        dims = column_dims(def)
        if length(dims) > 1
            write_tdim(file, k, dims)
        end
    end
    # The number of HDUs as returned by fits_get_num_hdus is only incremented
    # after writing data.
    n = position(file)
    if length(file) < n
        setfield!(file, :nhdus, n)
    end
    return FitsTableHDU(BareBuild(), file, n, ascii)
end

"""
    write(hdu::FitsTableHDU, col => arr, ...; first=hdu.first_row, case=false, null=nothing) -> hdu

writes column data into FITS table extension of `hdu`. `col` is the column
name/number, `arr` is an array of column values. Column values are converted as
needed and are written starting at the row specified by `first`. The leading
dimensions of `arr` should be the same as those specified by the corresponding
`TDIMn` keyword (with `n` the column number) and the remaining last dimension,
if any, corresponds to the *row* index of the table.

If `col` is a column name, keyword `case` may be used to indicate whether case
of letters matters (default is `false`).

Keyword `null` may be used to specify the value of undefined elements in `arr`.

Any number of columns may be specified as subsequent arguments. The same
keywords apply to all columns.

"""
function write(hdu::FitsTableHDU,
               pair::Pair{<:Integer,<:AbstractArray};
               case::Bool = false,
               null::Union{Number,Nothing} = nothing,
               first::Integer = hdu.first_row)
    col = pair.first
    arr = dense_array(pair.second)
    write_col(hdu, col, first, 1, arr, null)
    return hdu
end

function write(hdu::FitsTableHDU,
               pair::Pair{<:AbstractString,<:AbstractArray};
               case::Bool = false, kwds...)
    col = get_colnum(hdu, pair.first, case)
    return write(hdu, col => pair.second; kwds...)
end

# Union of possible types for specifying column data to write.
const ColumnInputData = Union{Pair{<:Integer,<:AbstractArray},
                              Pair{<:AbstractString,<:AbstractArray}}

# Write more than one column at a time.
@inline function write(hdu::FitsTableHDU, data::ColumnInputData,
                       args::ColumnInputData...; kwds...)
    write(hdu, data; kwds...)
    return write(hdu, args...; kwds...)
end

@inline function write(hdu::FitsTableHDU,
                       data::Tuple{Vararg{ColumnInputData}};
                       kwds...)
    write(hdu, Base.first(data); kwds...)
    return write(hdu, Base.tail(data)...; kwds...)
end

function write(hdu::FitsTableHDU,
               # NOTE: For a vector of input column data, it is not possible to
               # be more specific for the key type if we want to allow for a
               # mixture of key types.
               pairs::AbstractVector{<:Pair{<:Any,<:AbstractArray}};
               kwds...)
    for pair in pairs
        write(hdu, pair)
    end
    return hdu
end

# No more columns to write.
write(hdu::FitsTableHDU, ::Tuple{}; kdws...) = hdu
write(hdu::FitsTableHDU; kdws...) = hdu

#------------------------------------------------------------------------------
