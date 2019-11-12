@deprecate(loadfits(args...; kwds...),
           readfits(args...; kwds...))
@deprecate(getcomment(args...; kwds...),
           getfitscomment(args...; kwds...))
@deprecate(setkey!(args...; kwds...),
           setfitskey!(args...; kwds...))
@deprecate(tryreadkey(args...; kwds...),
           tryreadfitskey(args...; kwds...))
@deprecate(tryreadkeys(args...; kwds...),
           tryreadfitskeys(args...; kwds...))
@deprecate(getheader(args...; kwds...),
           getfitsheader(args...; kwds...))
@deprecate(getdata(args...; kwds...),
           getfitsdata(args...; kwds...))
@deprecate(parent(A::FitsImage),
           getfitsdata(A))
@deprecate(Image, FitsImage)
@deprecate(header(),
           FitsHeader())
@deprecate(getheader(obj),
           get(FITSHeader, obj))
@deprecate(getfitsheader(obj),
           get(FitsHeader, obj))
@deprecate(getfitsdata(obj),
           get(Array, obj))
@deprecate(getfitscomment(obj),
           get(FitsComment, obj))
@deprecate(openfits(path::AbstractString),
           FitsIO(path, "r"))
@deprecate(openfits(func::Function, path::AbstractString),
           FitsIO(func, path, "r"))
@deprecate(createfits(path::AbstractString; overwrite::Bool=false),
           FitsIO(path, (overwrite ? "w!" : "w")))
@deprecate(createfits(func::Function, path::AbstractString; overwrite::Bool=false),
           FitsIO(func, path, (overwrite ? "w!" : "w")))
@deprecate(createfits!(path::AbstractString),
           FitsIO(path, "w!"))
@deprecate(createfits!(func::Function, path::AbstractString),
           FitsIO(func, path, "w!"))
@deprecate(readfits(T::Type, fh::FITS, ext::Union{Integer,AbstractString} = 1),
           read(T, FitsHDU(fh[ext])))
@deprecate(readfits(fh::FITS, args...; kwds...),
           read(FitsImage, FitsIO(fh), args...; kwds...))
@deprecate(readfits(T::Type{<:FitsImage}, hdu::HDU),
           read(T, FitsHDU(hdu)))
@deprecate(readfits(T::Type{<:Array}, hdu::ImageHDU, args...),
           read(T, FitsHDU(hdu), args...))
@deprecate(readfits(::Type{FITSHeader}, hdu::HDU),
           read_header(hdu))
@deprecate(writefits(fh::FITS, args...; kwds...),
           write(FitsIO(fh), args...; kwds...))
@deprecate(write(fh::FITS, arr::AbstractArray, hdr::FitsHeader),
           write(FitsIO(fh), arr, hdr))
@deprecate(write(fh::FITS, obj::FitsImage),
           write(FitsIO(fh), obj))
@deprecate(FitsHeader(fh::FITS, args...),
           FitsHeader(FitsIO(fh), args...))
@deprecate(FitsHeader(hdu::HDU),
           FitsHeader(FitsHDU(hdu)))
@deprecate(setfitskey!(hdr::FITSHeader, args...),
           setkey!(FitsHeader(hdr), args...))
@deprecate(setfitskey!(obj::Annotated, args...),
           setkey!(obj, args...))
@deprecate(setkey!(obj::Annotated, key::AbstractString, val, com),
           setindex!(obj, (val, com), key))
@deprecate(setkey!(obj::Annotated, key::AbstractString, val),
           setindex!(obj, val, key))
@deprecate(getfitskey(obj::Annotated, key::AbstractString),
           getindex(obj, key))
@deprecate(getfitskey(T::Type, obj::Annotated, args...),
           get(T, obj, args...))
@deprecate(getfitskey(T::Type, hdr::FITSHeader, args...),
           get(T, FitsHeader(hdr), args...))
@deprecate(extname(hdu::HDU),
           extname(FitsHDU(hdu)))
@deprecate(extversion(hdu::HDU),
           extversion(FitsHDU(hdu)))
@deprecate(hduname(hdu::HDU),
           hduname(FitsHDU(hdu)))
@deprecate(hduversion(hdu::HDU),
           hduversion(FitsHDU(hdu)))
@deprecate(hduversion(hdu::FitsHDU),
           read(Int, hdu, "HDUVER", 1))
@deprecate(hduname(hdu::FitsHDU),
           read(String, hdu, "HDUNAME", nothing))
@deprecate(extversion(hdu::FitsHDU),
           read(Int, hdu, "EXTVER", 1))
@deprecate(extname(hdu::FitsHDU),
           read(String, hdu, "EXTNAME", nothing))
@deprecate(tryreadfitskey(obj, T::Type, key),
           read(T, obj, key, nothing))
@deprecate(writefits(path::AbstractString, args...; kwds...),
           write(FitsFile, path, args...; kwds...))
@deprecate(writefits!(path::AbstractString, args...; kwds...),
           write(FitsFile, path, args...; overwrite=true, kwds...))
@deprecate(readfits(T::Type{<:Annotated}, path::AbstractString, args...; kwds...),
           read(T, path, args...; kwds...))
@deprecate(readfits(path::AbstractString, args...; kwds...),
           readfits(FitsImage, path, args...; kwds...))
@deprecate(readfits(T::Type{<:Array}, path::AbstractString, args...; kwds...),
           FitsIO(path, "r") do io; read(T, io, args...; kwds...); end)
