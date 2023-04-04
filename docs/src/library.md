# Reference

## FITS Tables Extensions

```@docs
names(::FitsTableHDU)
read(::FitsTableHDU, ::Column)
read(::Type{Dict}, ::FitsTableHDU, ::Columns, ::Rows)
read(::Type{Vector}, ::FitsTableHDU, ::Columns, ::Rows)
read!(::DenseArray, ::FitsTableHDU, ::ColumnName)
read!(::AbstractDict, ::FitsTableHDU, ::Columns, ::Rows)
merge!(::AbstractDict{String,<:Array}, ::FitsTableHDU, ::Columns, ::Rows)
push!(::AbstractVector{<:AbstractArray}, ::FitsTableHDU, ::Columns, ::Rows)
write(::FitsFile, ::Type{FitsTableHDU}, ::Pair{<:ColumnName,<:Any}..)
write(::FitsTableHDU, ::ColumnDataPair)
write(f::FitsFile, ::Union{Nothing,Header}, ::TableData)
```

## FITS Image Extensions

```@docs
read(::FitsImageHDU)
read!(::DenseArray{<:Number}, ::FitsImageHDU, ::SubArrayIndex...)
read!(::DenseArray, ::FitsImageHDU)
write(::FitsFile, ::Type{FitsImageHDU{T}}, ::NTuple{N,Integer}) where {T,N}
write(::FitsFile, ::Union{Header,Nothing}, ::AbstractArray{T,N}) where {T<:Number,N}
write(::FitsImageHDU{<:Any,N}, ::AbstractArray{T,L}) where {T<:Number,L,N}
```
