#------------------------------------------------------------------------------
# FITS TABLES PROPERTIES

Base.propertynames(::FitsTableHDU) = (:nrows, :ncols, :column_names,
                                      :first_row, :last_row,
                                      :first_column, :last_column,
                                      :data_size, :data_ndims, :data_axes,
                                      :extname, :hduname, :io, :num, :type, :xtension)

Base.getproperty(hdu::FitsTableHDU, ::Val{:data_ndims}) = 2
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_size}) = (get_num_rows(hdu),
                                                          get_num_cols(hdu))
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_axes}) = (Base.OneTo(get_num_rows(hdu)),
                                                          Base.OneTo(get_num_cols(hdu)))

Base.getproperty(hdu::FitsTableHDU, ::Val{:first_row}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_row}) = get_num_rows(hdu)

Base.getproperty(hdu::FitsTableHDU, ::Val{:first_column}) = 1
Base.getproperty(hdu::FitsTableHDU, ::Val{:last_column}) = get_num_cols(hdu)

Base.getproperty(hdu::FitsTableHDU, ::Val{:nrows}) = get_num_rows(hdu)

function get_num_rows(f::Union{FitsIO,FitsTableHDU})
    nrows = Ref{Clonglong}()
    check(CFITSIO.fits_get_num_rowsll(f, nrows, Ref{Status}(0)))
    return to_type(Int, nrows[])
end

Base.getproperty(hdu::FitsTableHDU, ::Val{:ncols}) = get_num_cols(hdu)
function get_num_cols(f::Union{FitsIO,FitsTableHDU})
    ncols = Ref{Cint}()
    check(CFITSIO.fits_get_num_cols(f, ncols, Ref{Status}(0)))
    return to_type(Int, ncols[])
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
        names[col] = hdu["TTYPE$col"].value.string
    end
    return names
end

#------------------------------------------------------------------------------
# READING FITS TABLES

const ColumnName = Union{AbstractString,Symbol}

"""
    EasyFITS.get_colnum(hdu::FitsTableHDU, col, case=false) -> num

yields the column number of column matching `col` (a string, a symbol, or an
integer) in FITS table extension of `hdu`. Optional argument `case` specify
whether upper-/lower-case matters if `col` is a string or a symbol.

"""
function get_colnum(hdu::FitsTableHDU, col::Integer, case::Bool = false)
    @boundscheck check_colnum(hdu, col)
    return to_type(Int, col)
end

function get_colnum(hdu::FitsTableHDU, col::ColumnName, case::Bool = false)
    colnum = Ref{Cint}()
    check(CFITSIO.fits_get_colnum(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                  col, colnum, Ref{Status}(0)))
    return to_type(Int, colnum[])
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
    (card === nothing || card.type != FITS_STRING) && return ("COL#$col", to_type(Int, col))
    return (case ? card.value.string : uppercase(card.value.string), to_type(Int, col))
end

function get_colname(hdu::FitsTableHDU, col::ColumnName, case::Bool = false)
    colnum = Ref{Cint}()
    colname = Vector{UInt8}(undef, CFITSIO.FLEN_VALUE)
    check(CFITSIO.fits_get_colname(hdu, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                   col, pointer(colname), colnum, Ref{Status}(0)))
    return (to_string!(colname), to_type(Int, colnum[]))
end

for func in (:get_coltype, :get_eqcoltype)
    local cfunc = Symbol("fits_",func,"ll")
    @eval function $func(f::Union{FitsIO,FitsTableHDU}, col::ColumnName,
                         case::Bool = false)
        return $func(f, get_colnum(f, col, case))
    end
    @eval function $func(f::Union{FitsIO,FitsTableHDU}, col::Integer)
        type = Ref{Cint}()
        repeat = Ref{Clonglong}()
        width = Ref{Clonglong}()
        check(CFITSIO.$cfunc(f, col, type, repeat, width, Ref{Status}(0)))
        return (type[], to_type(Int, repeat[]), to_type(Int, width[]))
    end
end

read_tdim(f::Union{FitsIO,FitsTableHDU}, col::ColumnName, case::Bool=false) =
    read_tdim(f, get_colnum(f, col, case))

function read_tdim(f::Union{FitsIO,FitsTableHDU}, col::Integer)
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
    return to_type(Vector{Int}, naxes)
end

"""
    read(R::Type{<:Array}=Array, hdu::FitsTableHDU, col; kwds...) -> arr::R

reads column `col` of the FITS table extension in `hdu` and returns its values
as an array `arr` of type `R`.

The column `col` may be specified by its name or by its number. If `col` is a
string or a symbol, keyword `case` specifies whether uppercase/lowercase
matters (`case = false` by default).

Keywords `first` and/or `last` may be specified with the index of the
first/last row to read. By default, `first = hdu.first_row`, that is reading
starts at the first row of the table, and `last = hdu.last_row`, that is
reading stops at the last row of the table. An alternative is to specify the
row range after the column name/index:

    read([R,] hdu, col, first:last)

See [`read!(::DenseArray,::FitsTableHDU,::Integer)`](@ref) for the other
possible keywords.

"""
Base.read(hdu::FitsTableHDU, col::Union{Integer,AbstractString}; kwds...) =
    read(Array, hdu, col; kwds...)

function Base.read(hdu::FitsTableHDU, col::Union{Integer,AbstractString},
                   rng::AbstractUnitRange{<:Integer}; kwds...)
    return read(Array, hdu, col;
                first = Int(Base.first(rng)),
                last = Int(Base.last(rng)), kwds...)
end

function Base.read(hdu::FitsTableHDU, col::Union{Integer,AbstractString},
                   ::Colon; kwds...)
    return read(Array, hdu, col;
                first = hdu.first_row,
                last = hadu.last_row, kwds...)
end

function Base.read(::Type{R}, hdu::FitsTableHDU, col::AbstractString;
                   case::Bool = false, kwds...) where {R<:Array}
    return read(R, hdu, get_colnum(hdu, col, case); kwds...)
end

function Base.read(::Type{Array}, hdu::FitsTableHDU, col::Integer;
                   first::Integer = hdu.first_row,
                   last::Integer = hdu.last_row,
                   kwds...)
    type, repeat, width = get_eqcoltype(hdu, col)
    if type == CFITSIO.TSTRING
        return read(Array{String}, hdu, col;
                    first=first, last=last, kwds...)
    else
        dims = size_to_read(hdu, col, first, last, true)
        return read!(new_array(type_from_code(type), dims), hdu, col;
                     first=first, kwds...)
    end
end

function Base.read(::Type{Array{T}}, hdu::FitsTableHDU, col::Integer;
                   first::Integer = hdu.first_row,
                   last::Integer = hdu.last_row,
                   kwds...) where {T<:Number}
    dims = size_to_read(hdu, col, first, last, true)
    return read!(new_array(T, dims), hdu, col; first=first, kwds...)
end

function Base.read(::Type{Array{T,N}}, hdu::FitsTableHDU, col::Integer;
                   first::Integer = hdu.first_row,
                   last::Integer = hdu.last_row,
                   kwds...) where {T<:Number,N}
    dims = size_to_read(hdu, col, first, last, true)
    length(dims) == N || error("invalid number of dimensions")
    return read!(new_array(T, Val(N), dims), hdu, col; first=first, kwds...)
end

# Read array of strings as bytes. FIXME: null and anynull keywords

function Base.read(::Type{Array{String}}, hdu::FitsTableHDU, col::Integer;
                   first::Integer = hdu.first_row,
                   last::Integer = hdu.last_row,
                   kwds...)
    dims = size_to_read(hdu, col, first, last, false)
    return bytes_to_strings(read!(new_array(UInt8, dims), hdu, col;
                                  first=first, kwds...))
end

function Base.read(::Type{Array{String,N}}, hdu::FitsTableHDU, col::Integer;
                   first::Integer = hdu.first_row,
                   last::Integer = hdu.last_row,
                   kwds...) where {N}
    dims = size_to_read(hdu, col, first, last, false)
    length(dims) == N+1 || error("invalid number of dimensions")
    return bytes_to_strings(read!(new_array(T, Val(N+1), dims), hdu, col;
                                  first=first, kwds...))
end

# Yields size of column data as read.
function size_to_read(hdu::FitsTableHDU, col::Integer,
                      first::Integer = hdu.first_row,
                      last::Integer = hdu.last_row,
                      compress::Bool = true)
    nrows = max(Int(last + 1 - first)::Int, 0)
    if nrows > 0 && (first < hdu.first_row || last > hdu.last_row)
        bad_argument("out of bounds first row to read")
    end
    dims = read_tdim(hdu, col)
    if compress && length(dims) == 1 && Base.first(dims) == 1
        dims[firstindex(dims)] = nrows
    else
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
    read(R::Type{<:Dict} = Dict{String,Array}, hdu::FitsTableHDU[, cols[, rows]])

reads some columns of the FITS table extension in `hdu` as a dictionary indexed
by the column names. The columns to read can be specified by `cols` which may
be a single column name/index or a tuple/range/vector of column namas/indices.
Column names may be strings or keywords (not a mixture of these). The rows to
read can be specified by `rows` as a single row index or a unit range of row
indices. Keywords `first` and `last` may also be used to specify the range of
rows to read. By default, all columns and all rows are read.

"""
read(hdu::FitsTableHDU, args...; kwds...) =
    read(Dict{String,Array}, hdu, args...; kwds...)

read(::Type{Dict}, hdu::FitsTableHDU, args...; kwds...) =
    read(Dict{String,Array}, hdu, args...; kwds...)

read(::Type{Dict{String}}, hdu::FitsTableHDU, args...; kwds...) =
    read(Dict{String,Array}, hdu, args...; kwds...)

# Get rid of the selection of rows.
function read(::Type{Dict{String,Array}}, hdu::FitsTableHDU, cols,
              rows::Union{Integer,AbstractUnitRange{<:Integer}}; kwds...)
    return read(Dict{String,Array}, hdu, cols;
                first = Int(Base.first(rows))::Int,
                last = Int(Base.last(rows))::Int, kwds...)
end

# Read all rows.
function read(::Type{Dict{String,Array}}, hdu::FitsTableHDU, cols,
              rows::Colon; kwds...)
    return read(Dict{String,Array}, hdu, cols;
                first = hdu.first_row,
                last = hdu.last_row, kwds...)
end

# Read selection of columns.
function read(::Type{Dict{String,Array}},
              hdu::FitsTableHDU,
              cols::Union{Integer,
                          AbstractVector{<:Integer},
                          Tuple{Vararg{Integer}},
                          OrdinalRange{<:Integer,<:Integer},
                          AbstractString,
                          AbstractVector{<:AbstractString},
                          Tuple{Vararg{AbstractString}},
                          Symbol,
                          AbstractVector{<:Symbol},
                          Tuple{Vararg{Symbol}}} = hdu.first_column:hdu.last_column;
              first::Integer = hdu.first_row,
              last::Integer = hdu.last_row,
              case::Bool = false,
              kwds...)
    dict = Dict{String,Array}()
    names = hdu.column_names
    if !case
        for i in eachindex(names)
            names[i] = uppercase(names[i])
        end
    end
    for col in cols
        if col isa Integer
            num = to_type(Int, col)
        elseif col isa AbstractString || col isa Symbol
            num = get_colnum(hdu, col, case)
        end
        key = names[num]
        push!(dict, key => read(Array, hdu, num; first = first, last = last, kwds...))
    end
    return dict
end

"""
    read(R::Type{<:Union{Vector,Vector{Array}}}, hdu::FitsTableHDU[, cols[, rows]])

reads some columns of the FITS table extension in `hdu` as a vector. The
columns to read can be specified by `cols` which may be a single column
name/index or a tuple/range/vector of column names/indices. Column names may be
strings or keywords (not a mixture of these). The rows to read can be specified
by `rows` as a single row index or a unit range of row indices. Keywords
`first` and `last` may also be used to specify the range of rows to read. By
default, all columns and all rows are read.

"""
read(::Type{Vector}, hdu::FitsTableHDU, args...; kwds...) =
    read(Vector{Array}, hdu, args...; kwds...)

# Get rid of the selection of rows.
function read(::Type{Vector{Array}}, hdu::FitsTableHDU, cols,
              rows::Union{Integer,AbstractUnitRange{<:Integer}}; kwds...)
    return read(Vector{Array}, hdu, cols;
                first = Int(Base.first(rows))::Int,
                last = Int(Base.last(rows))::Int, kwds...)
end

# Read all rows.
function read(::Type{Vector{Array}}, hdu::FitsTableHDU, cols,
              rows::Colon; kwds...)
    return read(Vector{Array}, hdu, cols;
                first = hdu.first_row,
                last = hdu.last_row, kwds...)
    return read(Vector{Array}, hdu, cols; kwds...)
end

# Read selection of columns.
function read(::Type{Vector{Array}},
              hdu::FitsTableHDU,
              cols::Union{Integer,
                          AbstractVector{<:Integer},
                          Tuple{Vararg{Integer}},
                          OrdinalRange{<:Integer,<:Integer},
                          AbstractString,
                          AbstractVector{<:AbstractString},
                          Tuple{Vararg{AbstractString}},
                          Symbol,
                          AbstractVector{<:Symbol},
                          Tuple{Vararg{Symbol}}} = hdu.first_column:hdu.last_column;
              first::Integer = hdu.first_row,
              last::Integer = hdu.last_row,
              case::Bool = false,
              kwds...)
    vect = Vector{Array}(undef, length(cols))
    i = firstindex(vect)
    for col in cols
        if col isa Integer
            num = to_type(Int, col)
        elseif col isa AbstractString || col isa Symbol
            num = get_colnum(hdu, col, case)
        end
        vect[i] = read(Array, hdu, num; first = first, last = last, kwds...)
        i += 1
    end
    return vect
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
function Base.read!(arr::DenseArray, hdu::FitsTableHDU,
                    col::AbstractString; case::Bool = false, kwds...)
    return read!(arr, hdu, get_colnum(hdu, col, case); kwds...)
end

function Base.read!(arr::DenseArray{String,N}, hdu::FitsTableHDU, col::Integer;
                    kwds...) where {N}
    error("reading column of strings not yet implemented")
end

function Base.read!(arr::DenseArray{T,N}, hdu::FitsTableHDU, col::Integer;
                    null::Union{DenseArray{Bool,N},Number,Nothing} = nothing,
                    anynull::Union{Ref{Bool},Nothing} = nothing,
                    first::Integer = hdu.first_row) where {T<:Number,N}
    if null === nothing
        _null = zero(T)
    elseif null isa Number
        _null = to_type(T, null)
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
function write_tdim(f::Union{FitsIO,FitsTableHDU}, colnum::Integer,
                    dims::Tuple{Vararg{Integer}})
    if Clong === Clonglong || maximum(dims) ≤ typemax(Clong)
        write_tdim(f, colnum, convert(Tuple{Vararg{Clong}}, dims))
    else
        write_tdim(f, colnum, convert(Tuple{Vararg{Clonglong}}, dims))
    end
end
function write_tdim(f::Union{FitsIO,FitsTableHDU}, colnum::Integer,
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
        function write_tdim(f::Union{FitsIO,FitsTableHDU}, colnum::Integer,
                            dims::DenseVector{$T})
            check(CFITSIO.$cfunc(f, colnum, length(dims), dims, Ref{Status}(0)))
        end
        function write_tdim(f::Union{FitsIO,FitsTableHDU}, colnum::Integer,
                            dims::NTuple{N,$T}) where {N}
            check(CFITSIO.$cfunc(f, colnum, N, Ref(dims), Ref{Status}(0)))
        end
    end
end

"""
    write(io, FitsTableHDU, cols) -> hdu

creates a new FITS table extension in FITS file `io` with columns defined by
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
string. By default, `dims = 1` and `units = ""`.

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

    hdu = write(io, FitsTableHDU, cols)
    push!(hdu, key1 => val1) # add a first header keyword
    ...                      # add other header keywords
    write(hdu, col1 => arr1) # write a first column
    ...                      # write other columns

Such a table may be created in a single call:

    write(io, [key1 => val1, key2 => val2, ...], [col1 => arr1, col2 => arr2, ...])

where `key1 => val1`, `key2 => val2`, etc. specify header cards, while `col1 =>
arr1`, `col2 => arr2`, etc. specify columns names and associated data.

"""
function Base.write(io::FitsIO, ::Type{FitsTableHDU},
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
    check(CFITSIO.fits_create_tbl(io, tbl, nrows, ncols, ttype, tform, tunit,
                                  C_NULL, Ref{Status}(0)))
    k = 0
    for (name, def) in cols
        k += 1
        dims = column_dims(def)
        if length(dims) > 1
            write_tdim(io, k, dims)
        end
    end
    # The number of HDUs as returned by fits_get_num_hdus is only incremented
    # after writing data.
    n = position(io)
    if length(io) < n
        setfield!(io, :nhdus, n)
    end
    return FitsTableHDU(BareBuild(), io, n, ascii)
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
function Base.write(hdu::FitsTableHDU,
                    pair::Pair{<:Integer,<:AbstractArray};
                    case::Bool = false,
                    null::Union{Number,Nothing} = nothing,
                    first::Integer = hdu.first_row)
    col = pair.first
    arr = dense_array(pair.second)
    write_col(hdu, col, first, 1, arr, null)
    return hdu
end

function Base.write(hdu::FitsTableHDU,
                    pair::Pair{<:AbstractString,<:AbstractArray};
                    case::Bool = false, kwds...)
    col = get_colnum(hdu, pair.first, case)
    return write(hdu, col => pair.second; kwds...)
end

# Union of possible types for specifying column data to write.
const ColumnInputData = Union{Pair{<:Integer,<:AbstractArray},
                              Pair{<:AbstractString,<:AbstractArray}}

# Write more than one column at a time.
@inline function Base.write(hdu::FitsTableHDU, data::ColumnInputData,
                            args::ColumnInputData...; kwds...)
    write(hdu, data; kwds...)
    return write(hdu, args...; kwds...)
end

@inline function Base.write(hdu::FitsTableHDU,
                            data::Tuple{Vararg{ColumnInputData}};
                            kwds...)
    write(hdu, Base.first(data); kwds...)
    return write(hdu, Base.tail(data)...; kwds...)
end

function Base.write(hdu::FitsTableHDU,
                    # NOTE: For a vector of input column data, it is not
                    # possible to be more specific for the key type if we want
                    # to allow for a mixture of key types.
                    pairs::AbstractVector{<:Pair{<:Any,<:AbstractArray}};
                    kwds...)
    for pair in pairs
        write(hdu, pair)
    end
    return hdu
end

# No more columns to write.
Base.write(hdu::FitsTableHDU, ::Tuple{}; kdws...) = hdu
Base.write(hdu::FitsTableHDU; kdws...) = hdu

#------------------------------------------------------------------------------
