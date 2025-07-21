import Base: open, read, read!, write
@deprecate open(::Type{FitsFile}, filename::AbstractString, mode::AbstractString = "r"; kwds...) openfits(filename, mode; kwds...) false
@deprecate open(func::Function, ::Type{FitsFile}, filename::AbstractString, mode::AbstractString = "r"; kwds...) openfits(func, filename, mode; kwds...) false
@deprecate read!(dest, ::Type{FitsFile}, filename::AbstractString, args...; kwds...) readfits!(dest, filename, args...; kwds...) false
@deprecate read(::Type{FitsFile}, filename::AbstractString, args...; kwds...) readfits(filename, args...; kwds...) false
@deprecate read(::Type{R}, ::Type{FitsFile}, filename::AbstractString, args...; kwds...) where {R} readfits(R, filename, args...; kwds...) false
@deprecate write!(::Type{FitsFile}, filename::AbstractString, args...; kwds...) writefits(filename, args...; overwrite = true, kwds...) false
@deprecate write(::Type{FitsFile}, filename::AbstractString, args...; kwds...) writefits(filename, args...; kwds...) false
@deprecate write(::Type{T}, file::FitsFile, args...; kwds...) where {T<:FitsHDU} T(file, args...; kwds...) false
