#------------------------------------------------------------------------------
# FITS IMAGES PROPERTIES

Base.propertynames(::FitsImageHDU) = (:data_eltype, :data_size, :data_ndims,
                                      :extname, :hduname, :io, :num, :type,
                                      :xtension)

Base.getproperty(hdu::FitsImageHDU{T,N}, ::Val{:data_eltype}) where {T,N} = T
Base.getproperty(hdu::FitsImageHDU{T,N}, ::Val{:data_ndims}) where {T,N} = N
Base.getproperty(hdu::FitsImageHDU,      ::Val{:data_size}) = get_img_size(hdu)

"""
    write(io::FitsIO, FitsImageHDU, T, dims) -> hdu

creates a new primary array or image extension in FITS file `io` with a
specified data type `T` and size `dims`. If the FITS file is currently empty
then a primary array is created, otherwise a new image extension is appended to
the file.

An object to manage the new extension is returned which can be used to push
header cards and then to write the data. For example:

    hdu = write(io, FitsImageHDU, eltype(arr), size(arr))
    hdu["KEY1"] = val1
    hdu["KEY2"] = (val2, str2)
    hdu["KEY3"] = (nothing, str3)
    write(hdu, arr)

will create a new header data unit storing array `arr` with 3 header additional
header entries: one named `"KEY1"` with value `val1` and no comments, another
named `"KEY2"` with value `val2` and comment string `str2`, and yet another one
named `"KEY3"` with no value and with comment string `str3`. Note that special
names `"COMMENT"`, `"HISTORY"`, and `""` indicating commentary entries have no
associated, only a comment string, say `str` which can be specified as `str` or
as `(,str)`.

"""
function Base.write(io::FitsIO, ::Type{FitsImageHDU},
                    ::Type{T}, dims::NTuple{N,Integer}) where {T,N}
    # NOTE: All variants end up calling this type-stable version.
    check(CFITSIO.fits_create_img(io, type_to_bitpix(T), N,
                                  Ref(convert(NTuple{N,Clong}, dims)),
                                  Ref{Status}(0)))
    # The number of HDUs is only incremented after writing data.
    return FitsImageHDU{T,N}(CheckedArguments(), io, position(io))
end

# Just convert bitpix to type.
function Base.write(io::FitsIO, ::Type{FitsImageHDU}, bitpix::Integer,
                    dims::Tuple{Vararg{Integer}})
    return write(io, FitsImageHDU, type_from_bitpix(bitpix), dims)
end

# Just convert the dimensions.
function Base.write(io::FitsIO, ::Type{FitsImageHDU}, T::Union{Integer,Type},
                    dims::AbstractVector{<:Integer})
    off = firstindex(dims) - 1
    return write(io, FitsImageHDU, T, ntuple(i -> Clong(dims[i+off]), length(dims)))
end

"""
    write(io::FitsIO, hdr, arr::AbstractArray, args...) -> io

    write(io::FitsIO, arr::AbstractArray, hdr=nothing) -> io

"""
function Base.write(io::FitsIO, hdr::Union{Header,Nothing},
                    arr::AbstractArray{T,N}) where {T<:Number,N}
    # FIXME: improve type-stability
    write(push!(write(io, FitsImageHDU, T, size(arr)), hdr), arr)
    return io
end

function Base.write(io::FitsIO, arr::AbstractArray{<:Number},
                    hdr::Union{Header,Nothing} = nothing)
    return write(io, hdr, arr)
end

function Base.write(io::FitsIO, hdr::Union{Header,Nothing},
                    arr::AbstractArray{<:Number}, args...)
    write(io, hdr, arr)
    return write(io, args...)
end

function Base.write(io::FitsIO, arr::AbstractArray{<:Number},
                    hdr::Union{Header,Nothing}, args...)
    write(io, hdr, arr)
    return write(io, args...)
end

"""
    write(hdu::FitsImageHDU, arr::AbstractArray{<:Number}) -> hdu

writes all the elements of the array `arr` to the pixels of the FITS image
extension of the header data unit `hdu`. The element of `arr` are converted to
the type of the pixels. The default is to write the complete image pixels, so
the size of the destination `arr` must be identical to the dimensions of the FITS
image extension. This behavior may be changed by specifying another value than
`nothing` for the keyword `start`:

* To write a rectangular sub-image, specify keyword `start` with the
  coordinates of the first pixel to write as an `N`-tuple of integers, with `N`
  the number of dimensions of the FITS image extension. The destination `arr` may
  have fewer dimensions than `N`.

* To read consecutive pixels, specify keyword `start` with the index of the
  first pixel to read as an integer. The dimensions of the destination `arr` are
  not considered, only the number of elements of `arr` matters.

Unless writing a rectangular sub-image, keyword `null` may be used to specify
the value of undefined elements in `arr`. For integer FITS images, the FITS null
value is defined by the BLANK keyword (an error is returned if the BLANK
keyword doesn't exist). For floating point FITS images the special IEEE NaN
(Not-a-Number) value will be written into the FITS file.

"""
function Base.write(hdu::FitsImageHDU{<:Any,N},
                    arr::AbstractArray{T,L};
                    start::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
                    null::Union{Number,Nothing} = nothing) where {T<:Number,L,N}
    type = type_to_code(arr) # clash if unsupported element type
    dims = get_img_size(hdu)
    if start isa Tuple && null isa Nothing
        # Write a rectangular sub-image.
        ndims(arr) ≤ N || error("too many dimensions in rectangular sub-image to write")
        fpix = convert(NTuple{N,Clong}, start)::NTuple{N,Clong}
        lpix = last_pixel(fpix, size(arr))::NTuple{N,Clong}
        for k in 1:N
            1 ≤ fpix[k] ≤ lpix[k] ≤ dims[k] || error("out of range rectangular sub-image")
        end
        check(CFITSIO.fits_write_subset(hdu, type, Ref(fpix), Ref(lpix), dense_array(arr), Ref{Status}(0)))
    elseif start isa Union{Integer,Nothing}
        # Write a range of consecutive pixels.
        npix = length(arr)
        fpix = 1
        if start === nothing
            size(arr) == dims || error("not same dimensions")
        else
            fpix = oftype(fpix, start)
        end
        1 ≤ fpix && fpix - 1 + npix ≤ prod(dims) || error("out of range interval of pixels")
        if null === nothing
            check(CFITSIO.fits_write_img(hdu, type, fpix, npix, dense_array(arr), Ref{Status}(0)))
        else
            check(CFITSIO.fits_write_imgnull(hdu, type, fpix, npix, dense_array(arr),
                                             Ref{eltype(arr)}(null), Ref{Status}(0)))
        end
    else
        error("invalid combination of options")
    end
    return hdu
end

function Base.read(hdu::FitsImageHDU{T,N}) where {T,N}
    return read(Array{T,N}, hdu)
end

function Base.read(::Type{Array}, hdu::FitsImageHDU{T,N}) where {T,N}
    return read(Array{T,N}, hdu)
end

function Base.read(::Type{Array{T}}, hdu::FitsImageHDU{<:Any,N}) where {T<:Number,N}
    return read(Array{T,N}, hdu)
end

function Base.read(::Type{Array{T,N}}, hdu::FitsImageHDU{<:Any,N}) where {T<:Number,N}
    arr = Array{T,N}(undef, get_img_size(hdu))
    return read!(arr, hdu)
end

"""
    read!(arr::DenseArray{<:Number}, hdu::FitsImageHDU; kwds...) -> arr

overwrites all the elements of the array `arr` with pixel values read from the
FITS image extension of the header data unit `hdu` and returns `arr`. Pixel
values are converted to the element type of `arr`. The default is to read the
complete image, so the size of the destination `arr` must be identical to the
dimensions of the FITS image extension. This behavior may be changed by
specifying another value than `nothing` for the keyword `start`:

* To read a rectangular sub-image, specify keyword `start` with the coordinates
  of the first pixel to read as an `N`-tuple of integers, with `N` the number
  of dimensions of the FITS image extension. Optionally, keyword `step` may be
  specified as an `N`-tuple of integers to indicate the increment along each
  dimensions. The destination `arr` may have fewer dimensions than `N`.

* To read consecutive pixels, specify keyword `start` with the index of the
  first pixel to read as an integer. The dimensions of the destination `arr` are
  not considered, only the number of elements of `arr` matters.

Keyword `anynull` may be specified with a reference to a boolean
(`Ref{Bool}()`) to retrieve whether any of the read pixels is undefined.

Keyword `null` may be specified with a reference to a value of the same type as
the elements of the destination `arr` (`Ref{eltype(arr)}()`) to retrieve the value
of undefined pixels. Unless reading a rectangular sub-image, keyword `null` may
be set with an array of `Bool` of same size as `arr` and which will be set to
`1` for undefined pixels and to `0` elsewhere.

Output arrays `arr` and `null` must have contiguous elements, in other words,
they must be *dense arrays*.

"""
function Base.read!(arr::DenseArray{T,L},
                    hdu::FitsImageHDU{<:Any,N};
                    null::Union{DenseArray{Bool,L},Ref{T},Nothing} = nothing,
                    anynull::Union{Nothing,Ref{Bool}} = nothing,
                    start::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
                    step::Union{NTuple{N,Integer},Nothing} = nothing) where {T<:Number,L,N}
    type = type_to_code(T) # clash if unsupported element type
    dims = get_img_size(hdu)
    len = length(arr)
    anynul = Ref{Cint}(anynull isa Ref{Bool} && len > 0)
    if start isa Tuple && null isa Union{Ref{T},Nothing}
        # Read a rectangular sub-image.
        ndims(arr) ≤ N || error("too many dimensions in rectangular sub-image to read")
        if step === nothing
            ipix = NTuple(i -> one(Clong), Val(N))::NTuple{N,Clong}
        else
            ipix = convert(NTuple{N,Clong}, step)::NTuple{N,Clong}
        end
        fpix = convert(NTuple{N,Clong}, start)::NTuple{N,Clong}
        lpix = last_pixel(fpix, ipix, size(arr))::NTuple{N,Clong}
        for k in 1:N
            (ipix[k] ≥ 0 ?
                (fpix[k] ≥ 1) & (lpix[k] ≤ dims[k]) :
                (lpix[k] ≥ 1) & (fpix[k] ≤ dims[k])) || error("out of range rectangular sub-image")
        end
        if len > 0
            _null = (null === nothing ? C_NULL : null)
            check(CFITSIO.fits_read_subset(
                hdu, type, Ref(fpix), Ref(lpix), Ref(ipix), _null, arr, anynul, Ref{Status}(0)))
        end
    elseif start isa Union{Integer,Nothing} && step isa Nothing
        # Read a range of consecutive pixels.
        fpix = 1
        if start === nothing
            size(arr) == dims || error("not same dimensions")
        else
            fpix = oftype(fpix, start)
        end
        1 ≤ fpix && fpix - 1 + len ≤ prod(dims) || error("out of range interval of pixels")
        if null isa DenseArray{Bool,L}
            axes(null) == axes(arr) || error("incompatible array axes")
            if len > 0
                # NOTE: `GC.@protect null` is not needed here.
                _null = logical_pointer(null)
                check(CFITSIO.fits_read_imgnull(
                    hdu, type, fpix, len, _null, arr, anynul, Ref{Status}(0)))
                fix_logicals!(null)
            end
        elseif len > 0
            _null = (null === nothing ? C_NULL : null)
            check(CFITSIO.fits_read_img(
                hdu, type, fpix, len, _null, arr, anynul, Ref{Status}(0)))
        end
    else
        error("invalid combination of options")
    end
    if anynull isa Ref{Bool}
        anynull[] = !iszero(anynul[])
    end
    return arr
end

last_pixel(fpix::NTuple{N,T}, dims::Dims{L}) where {L,N,T<:Integer} =
    ntuple(k -> k > L ? fpix[k] : fpix[k] + T(dims[k] - 1)::L*ipix[k], Val(N))
last_pixel(fpix::NTuple{N,T}, ipix::NTuple{N,T}, dims::Dims{L}) where {L,N,T<:Integer} =
    ntuple(k -> k > L ? fpix[k] : fpix[k] + T(dims[k] - 1)::L*ipix[k], Val(N))

dense_array(arr::DenseArray) = arr
dense_array(arr::AbstractArray{T,N}) where {T,N} = convert(Array{T,N}, arr)

function get_img_type(f::Union{FitsIO,FitsImageHDU})
    bitpix = Ref{Cint}()
    check(CFITSIO.fits_get_img_type(f, bitpix, Ref{Status}(0)))
    return Int(bitpix[])
end

function get_img_equivtype(f::Union{FitsIO,FitsImageHDU})
    bitpix = Ref{Cint}()
    check(CFITSIO.fits_get_img_equivtype(f, bitpix, Ref{Status}(0)))
    return Int(bitpix[])
end

function get_img_dim(f::Union{FitsIO,FitsImageHDU})
    naxis = Ref{Cint}()
    check(CFITSIO.fits_get_img_dim(f, naxis, Ref{Status}(0)))
    return Int(naxis[])
end

function get_img_size(hdu::FitsImageHDU{T,N}) where {T,N}
    io = hdu.io # small optimization to avoid multiple seeks
    N == get_img_dim(io) || error("number of dimensions has changed, rebuild the HDU")
    dims = Ref{NTuple{N,Clong}}()
    check(CFITSIO.fits_get_img_size(hdu, N, dims, Ref{Status}(0)))
    return Dims{N}(dims[])
end
