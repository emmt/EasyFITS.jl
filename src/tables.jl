#------------------------------------------------------------------------------
# FITS TABLES PROPERTIES

Base.propertynames(::FitsTableHDU) = (:nrows, :ncols, :column_names,
                                      :first_column, :last_column,
                                      :data_size, :data_ndims,
                                      :extname, :hduname, :io, :num, :type, :xtension)

Base.getproperty(hdu::FitsTableHDU, ::Val{:data_ndims}) = 2
Base.getproperty(hdu::FitsTableHDU, ::Val{:data_size}) = (get_num_rows(hdu),
                                                          get_num_cols(hdu))

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

function get_colnum(f::Union{FitsIO,FitsTableHDU}, col::AbstractString,
                    case::Bool = false)
    colnum = Ref{Cint}()
    check(CFITSIO.fits_get_colnum(f, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                  col, colnum, Ref{Status}(0)))
    return to_type(Int, colnum[])
end

function get_colname(f::Union{FitsIO,FitsTableHDU}, col::AbstractString,
                     case::Bool = false)
    colnum = Ref{Cint}()
    colname = Vector{UInt8}(undef, CFITSIO.FLEN_VALUE)
    check(CFITSIO.fits_get_colname(f, (case ? CFITSIO.CASESEN : CFITSIO.CASEINSEN),
                                   col, pointer(colname), colnum, Ref{Status}(0)))
    return (to_string!(colname), to_type(Int, colnum[]))
end

for func in (:get_coltype, :get_eqcoltype)
    local cfunc = Symbol("fits_",func,"ll")
    @eval function $func(f::Union{FitsIO,FitsTableHDU}, col::Integer)
        type = Ref{Cint}()
        repeat = Ref{Clonglong}()
        width = Ref{Clonglong}()
        check(CFITSIO.$cfunc(f, col, type, repeat, width, Ref{Status}(0)))
        return (type[], to_type(Int, repeat[]), to_type(Int, width[]))
    end
end

read_tdim(f::Union{FitsIO,FitsTableHDU}, col::AbstractString, case::Bool=false) =
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
string, keyword `case` specifies whether uppercase/lowercase matters (`case =
false` by default).

Keyword `firstrow` may be specified with the index of the first row to read. By
default, `firstrow = 1` and reading starts at the first row of the table.

Keyword `nrows` may be specified with the number of rows to read. By default,
all the rows at and after row specified by keyword `firstrow` are read.

See [`read!(::DenseArray,::FitsTableHDU,::Integer)`](@ref) for the other
possible keywords.

"""
Base.read(hdu::FitsTableHDU, col::Union{Integer,AbstractString}; kwds...) =
    read(Array, hdu, col; kwds...)

function Base.read(::Type{R}, hdu::FitsTableHDU, col::AbstractString;
                   case::Bool = false, kwds...) where {R<:Array}
    return read(R, hdu, get_colnum(hdu, col, case); kwds...)
end

function Base.read(::Type{Array}, hdu::FitsTableHDU, col::Integer;
                   firstrow::Integer = 1,
                   nrows::Union{Integer,Nothing} = nothing,
                   kwds...)
    type, repeat, width = get_eqcoltype(hdu, col)
    if type == CFITSIO.TSTRING
        return read(Array{String}, hdu, col;
                    firstrow=firstrow, nrows=nrows, kwds...)
    else
        dims = size_to_read(hdu, col, firstrow, nrows, true)
        return read!(new_array(type_from_code(type), dims), hdu, col;
                     firstrow=firstrow, kwds...)
    end
end

function Base.read(::Type{Array{T}}, hdu::FitsTableHDU, col::Integer;
                   firstrow::Integer = 1,
                   nrows::Union{Integer,Nothing} = nothing,
                   kwds...) where {T<:Number}
    dims = size_to_read(hdu, col, firstrow, nrows, true)
    return read!(new_array(T, dims), hdu, col; firstrow=firstrow, kwds...)
end

function Base.read(::Type{Array{T,N}}, hdu::FitsTableHDU, col::Integer;
                   firstrow::Integer = 1,
                   nrows::Union{Integer,Nothing} = nothing,
                   kwds...) where {T<:Number,N}
    dims = size_to_read(hdu, col, firstrow, nrows, true)
    length(dims) == N || error("invalid number of dimensions")
    return read!(new_array(T, Val(N), dims), hdu, col; firstrow=firstrow, kwds...)
end

# Read array of strings as bytes. FIXME: null and anynull keywords

function Base.read(::Type{Array{String}}, hdu::FitsTableHDU, col::Integer;
                   firstrow::Integer = 1,
                   nrows::Union{Integer,Nothing} = nothing,
                   kwds...)
    dims = size_to_read(hdu, col, firstrow, nrows, false)
    return bytes_to_strings(read!(new_array(UInt8, dims), hdu, col;
                                  firstrow=firstrow, kwds...))
end

function Base.read(::Type{Array{String,N}}, hdu::FitsTableHDU, col::Integer;
                   firstrow::Integer = 1,
                   nrows::Union{Integer,Nothing} = nothing,
                   kwds...) where {N}
    dims = size_to_read(hdu, col, firstrow, nrows, false)
    length(dims) == N+1 || error("invalid number of dimensions")
    return bytes_to_strings(read!(new_array(T, Val(N+1), dims), hdu, col;
                                  firstrow=firstrow, kwds...))
end

function nrows_to_read(hdu::FitsTableHDU, firstrow::Integer, nrows::Nothing = nothing)
    1 ≤ firstrow || bad_argument("out of range value for `firstrow` keyword")
    return max(0, Int(hdu.nrows + 1 - firstrow)::Int)
end

function nrows_to_read(hdu::FitsTableHDU, firstrow::Integer, nrows::Integer)
    nrows ≤ nrows_to_read(hdu, firstrow) || bad_argument("too many rows to read")
    return Int(nrows)::Int
end

# Yields size of column data as read.
function size_to_read(hdu::FitsTableHDU, col::Integer,
                      firstrow::Integer = 1,
                      nrows::Union{Integer,Nothing} = nothing,
                      compress::Bool = true)
    _nrows = nrows_to_read(hdu, firstrow, nrows)
    dims = read_tdim(hdu, col)
    if compress && length(dims) == 1 && first(dims) == 1
        dims[firstindex(dims)] = _nrows
    else
        push!(dims, _nrows)
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
    read!(arr::DenseArray, hdu::FitsTableHDU, col) -> arr

overwrites the elements of array `arr` with values of the column `col` of the
FITS table extension in `hdu` and returns `arr`.

The column `col` may be specified by its name or by its number. If `col` is a
string, keyword `case` specifies whether uppercase/lowercase matters (`case =
false` by default).

Keyword `firstrow` may be specified with the index of the first row to read. By
default, `firstrow = 1` and reading starts at the first row of the table.

Keyword `anynull` may be specified with a reference to a boolean
(`Ref{Bool}()`) to retrieve whether any of the read values is undefined.

Keyword `null` may be specified with a reference to a value of the same type as
the elements of the destination `arr` (`Ref{eltype(arr)}()`) to retrieve the
value of undefined values. Keyword `null` may also be set with an array of
`Bool` of same size as `arr` and which will be set to `1` for undefined values
and to `0` elsewhere.

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
                    firstrow::Integer = 1) where {T<:Number,N}
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
    read_col!(hdu, col, firstrow, 1, arr, _null, _anynull)
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
    write(hdu::FitsTableHDU, col => arr, ...; first=1, case=false, null=nothing) -> hdu

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
                    first::Integer = 1)
    col = pair.first
    arr = dense_array(pair.second)
    write_col(hdu, col, first, 1, arr, null)
    return hdu
end

function Base.write(hdu::FitsTableHDU,
                    pair::Pair{<:AbstractString,<:AbstractArray};
                    case::Bool = false, kwds...)
    col = get_colnum(hdu, first(pair), case)
    return write(hdu, col => last(pair); kwds...)
end

function Base.write(hdu::FitsTableHDU,
                    pairs::Union{Pair{<:Integer,<:AbstractArray},
                                 Pair{<:AbstractString,<:AbstractArray}}...;
                    kwds...)
    for pair in pairs
        write(hdu, pair; kwds...)
    end
    return hdu
end

#------------------------------------------------------------------------------
