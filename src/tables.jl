#------------------------------------------------------------------------------
# FITS TABLES PROPERTIES

Base.propertynames(::FitsTableHDU) = (:nrows, :rows, :first_row, :last_row,
                                      :ncols, :columns, :first_column, :last_column,
                                      :column_name, :column_names, :column_number, :column_units,
                                      :data_size, :data_ndims, :data_axes,
                                      :extname, :hduname, :file, :num, :type, :xtension)

Base.getproperty(hdu::FitsTableHDU, ::Val{:data_ndims}) = 2
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_size}) = (hdu.nrows, hdu.ncols)
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_axes}) = (hdu.rows, hdu.columns)

Base.getproperty(hdu::FitsTableHDU, ::Val{:rows}) = Base.OneTo(hdu.nrows)
Base.getproperty(hdu::FitsTableHDU, ::Val{:first_row}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_row}) = hdu.nrows
Base.getproperty(hdu::FitsTableHDU, ::Val{:nrows}) = get_nrows(hdu)
function get_nrows(f::Union{FitsFile,FitsTableHDU})
    nrows = Ref{Clonglong}()
    check(CFITSIO.fits_get_num_rowsll(f, nrows, Ref{Status}(0)))
    return as(Int, nrows[])
end

Base.getproperty(hdu::FitsTableHDU, ::Val{:columns}) = Base.OneTo(hdu.ncols)
Base.getproperty(hdu::FitsTableHDU, ::Val{:first_column}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_column}) = hdu.ncols
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
              rows::Rows = Colon(); case::Bool = false,
              kwds...) where {A<:AbstractArray}
    num = get_colnum(hdu, col; case)
    return (read(A, hdu, num, rows; kwds...), get_units(hdu, num))
end

function read(::Type{A}, hdu::FitsTableHDU, col::ColumnName,
              rows::Rows = Colon(); case::Bool = false,
              kwds...) where {A<:AbstractArray}
    return read(A, hdu, get_colnum(hdu, col; case), rows; kwds...)
end

function read(::Type{Array}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...)
    type, repeat, width = get_eqcoltype(hdu, col)
    if type == CFITSIO.TSTRING
        return read(Array{String}, hdu, col, rows; kwds...)
    else
        dims = size_to_read(Number, hdu, col, rows)
        return read!(new_array(type_from_code(type), dims), hdu, col;
                     first = first_row_to_read(hdu, rows), kwds...)
    end
end

function read(::Type{Array{T}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {T<:Number}
    dims = size_to_read(Number, hdu, col, rows)
    return read!(new_array(T, dims), hdu, col;
                 first = first_row_to_read(hdu, rows), kwds...)
end

function read(::Type{Array{T,N}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {T<:Number,N}
    dims = size_to_read(Number, hdu, col, rows)
    length(dims) == N || error("invalid number of dimensions")
    return read!(new_array(T, Val(N), dims), hdu, col;
                 first = first_row_to_read(hdu, rows), kwds...)
end

# Read array of strings as bytes. FIXME: null and anynull keywords

function read(::Type{Array{String}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...)
    dims = size_to_read(String, hdu, col, rows)
    return bytes_to_strings(read!(new_array(UInt8, dims), hdu, col;
                                  first = first_row_to_read(hdu, rows), kwds...))
end

function read(::Type{Array{String,N}}, hdu::FitsTableHDU, col::Integer,
              rows::Rows = Colon(); kwds...) where {N}
    dims = size_to_read(String, hdu, col, rows)
    length(dims) == N+1 || error("invalid number of dimensions")
    return bytes_to_strings(read!(new_array(UInt8, Val(N+1), dims), hdu, col;
                                  first = first_row_to_read(hdu, rows), kwds...))
end

# Yields size of column data as read. If reading numbers and cell has size
# (1,), an empty size is assumed.
function size_to_read(::Type{T}, hdu::FitsTableHDU,
                      col::Integer, rows::Rows) where {T<:Union{String,Number}}
    all_rows = hdu.rows
    if rows isa Colon
        nrows = length(all_rows)
    else
        nrows = length(rows)
        nrows == 0 || (first(rows) ≥ first(all_rows) && last(rows) ≤ last(all_rows)) ||
            bad_argument("out of bounds row(s) to read")
    end
    dims = read_tdim(hdu, col)
    if T !== String && length(dims) == 1 && first(dims) == 1
        # Cells of this column have no dimensions and do not have string
        # values.
        if rows isa Integer
            empty!(dims)
        else
            dims[firstindex(dims)] = nrows
        end
    elseif !(rows isa Integer)
        push!(dims, nrows)
    end
    return dims
end

# Convert an array of bytes to an arry of strings. The leading dimension
# indexes the characters of the strings. Trailing spaces are stripped.
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
            buf[l += 1] = c
            if c != space
                len = l
            end
        end
        B[k += 1] = GC.@preserve buf unsafe_string(pointer(buf), len)
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
function push!(vec::AbstractVector{<:AbstractArray},
               hdu::FitsTableHDU, cols::Columns = Colon(), rows::Rows = Colon();
               case::Bool = false, kwds...)
    for col in columns_to_read(hdu, cols)
        num = get_colnum(hdu, col; case)
        push!(vec, read(A, hdu, num, rows; kwds...))
    end
    return vec
end

function push!(vec::AbstractVector{<:Tuple{<:AbstractArray,String}},
               hdu::FitsTableHDU, cols::Columns = Colon(), rows::Rows = Colon();
               case::Bool = false, kwds...)
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
    return write(file, FITSTableHDU, cols; kwds...)
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
        ttype[k] = name
        tform[k] = column_tform(def)
        tunit[k] = column_tunit(def)
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
               null::Union{Number,Nothing} = nothing,
               first::Integer = hdu.first_row)
    # Private function to retrieve the values and the units of a column.
    data_units(col::Pair{<:Column,<:ColumnData}) = (pair.second, nothing)
    data_units(col::Pair{<:Column,<:Tuple{ColumnData,ColumnUnits}}) = pair.second

    # Get column number.
    num = get_colnum(hdu, pair.first; case)

    # Retrieve column values and, if units are specified, check that they are
    # the same as those already set.
    vals, units = data_units(pair)
    if units !== nothing
        card = get(hdu, "TUNIT$num", nothing)
        (card === nothing ? isempty(units) :
            cart.type == FITS_STRING && units == card.string) || bad_argument("invalid column units")
    end

    # Write column values.
    unsafe_write_col(hdu, num, first, 1, length(vals), dense_array(vals), null)
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
function write(file::FitsFile, header::Union{Nothing,Header},
               cols::TableData; ascii::Bool=false)

    @noinline bad_column_data(key, val) = throw(ArgumentError("invalid data for column $key"))

    function max_length(A::AbstractArray{<:AbstractString})
        maxlen = 0
        for str in A
            maxlen = max(maxlen, length(str))
        end
        return maxlen
    end

    column_definition(key, val::ColumnData{T,1}) where {T} = (type_to_letter(T), ())
    column_definition(key, val::ColumnData{T,N}) where {T,N} = (type_to_letter(T), size(val)[1:N-1])
    column_definition(key, val::ColumnData{<:AbstractString,1})  = ('A', (max_length(val),))
    column_definition(key, val::ColumnData{<:AbstractString,N}) where {N} = ('A', (max_length(val),
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
