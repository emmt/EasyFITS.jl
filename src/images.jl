#------------------------------------------------------------------------------
# FITS IMAGES PROPERTIES

Base.propertynames(::FitsImageHDU) = (:data_eltype, :data_size, :data_ndims,
                                      :extname, :hduname, :file, :number, :type,
                                      :xtension)

Base.getproperty(hdu::FitsImageHDU{T,N}, ::Val{:data_eltype}) where {T,N} = T
Base.getproperty(hdu::FitsImageHDU{T,N}, ::Val{:data_ndims}) where {T,N} = N
Base.getproperty(hdu::FitsImageHDU,      ::Val{:data_size}) = get_img_size(hdu)

function get_img_type(f::Union{FitsFile,FitsImageHDU})
    bitpix = Ref{Cint}()
    check(CFITSIO.fits_get_img_type(f, bitpix, Ref{Status}(0)))
    return Int(bitpix[])
end

function get_img_equivtype(f::Union{FitsFile,FitsImageHDU})
    bitpix = Ref{Cint}()
    check(CFITSIO.fits_get_img_equivtype(f, bitpix, Ref{Status}(0)))
    return Int(bitpix[])
end

function get_img_dim(f::Union{FitsFile,FitsImageHDU})
    naxis = Ref{Cint}()
    check(CFITSIO.fits_get_img_dim(f, naxis, Ref{Status}(0)))
    return Int(naxis[])
end

function get_img_size(hdu::FitsImageHDU{T,N}) where {T,N}
    file = get_file_at(hdu) # small optimization to avoid multiple seeks
    N == get_img_dim(file) || error("number of dimensions has changed, rebuild the HDU")
    dims = Ref{NTuple{N,Clong}}()
    check(CFITSIO.fits_get_img_size(file, N, dims, Ref{Status}(0)))
    return Dims{N}(dims[])
end

#------------------------------------------------------------------------------
# READING IMAGES

"""
    read(R::Type = Array, hdu::FitsImageHDU) -> arr::R

reads the FITS image extension in `hdu`. Optional argument `R` is to restrict
the ouput type and improve type-stability. For example:

    arr = convert(Array{Float32}, read(hdu))
    arr = read(Array{Float32}, hdu)
    arr = read(Array{Float32,2}, hdu)

yield similar results if `hdu` is a 2-D image extension but the 2nd example is
more efficient than the 1st one as no temporary array is nedded if the pixel
type is not equivalent to `Float32` and the 3rd example is completely
type-stable.

Keywords `anynull` and `null` may be specified to deal with undefined pixel
values (see documentation for `read!`).

"""
function read(hdu::FitsImageHDU{T,N}; kwds...) where {T,N}
    return read(Array{T,N}, hdu; kwds...)
end

function read(::Type{Array}, hdu::FitsImageHDU{T,N}; kwds...) where {T,N}
    return read(Array{T,N}, hdu; kwds...)
end

function read(::Type{Array{T}}, hdu::FitsImageHDU{<:Any,N}; kwds...) where {T<:Number,N}
    return read(Array{T,N}, hdu; kwds...)
end

function read(::Type{Array{T,N}}, hdu::FitsImageHDU{<:Any,N};
              null::Union{DenseArray{Bool,N},Ref{T},Nothing} = nothing,
              anynull::Union{Nothing,Ref{Bool}} = nothing) where {T<:Number,N}
    arr = Array{T,N}(undef, get_img_size(hdu))
    return read!(arr, hdu; null, anynull)
end

"""
    read(R::Type = Array, hdu::FitsImageHDU, inds...) -> arr::R

reads a rectangular sub-region of the FITS image extension in `hdu` defined by
the indices `inds...`. Keywords `anynull` and `null` may be specified to deal
with undefined pixel values (see documentation for `read!`).

The result is similar to:

    arr = read(R, hdu)[inds...]

but may be more efficient as no array other than the result is allocated and
fewer values are read.

"""
function read(hdu::FitsImageHDU{T}, inds::SubArrayIndex...; kwds...) where {T}
    return read(Array{T}, hdu, inds; kwds...)
end

function read(hdu::FitsImageHDU{T}, inds::SubArrayIndices; kwds...) where {T}
    return read(Array{T}, hdu, inds; kwds...)
end

function read(::Type{R}, hdu::FitsImageHDU, inds::SubArrayIndex...; kwds...) where {R<:Array}
    return read(R, hdu, inds; kwds...)
end

function read(::Type{Array}, hdu::FitsImageHDU{T}, inds::SubArrayIndices; kwds...) where {T}
    return read(Array{T}, hdu, inds; kwds...)
end

function read(::Type{Array{T}}, hdu::FitsImageHDU, inds::SubArrayIndices; kwds...) where {T}
    N = count(i -> !isa(i, Integer), inds) # count number of output dimensions
    return read(Array{T,N}, hdu, inds; kwds...)
end

function read(::Type{Array{T,N}}, hdu::FitsImageHDU, inds::SubArrayIndices;
              null::Union{DenseArray{Bool,N},Ref{T},Nothing} = nothing,
              anynull::Union{Nothing,Ref{Bool}} = nothing) where {T,N}
    img_dims = get_img_size(hdu)
    arr_dims, first, step, last = subarray_params(img_dims, inds)
    length(arr_dims) == N || throw(DimensionMismatch(
        "given indices yield $(length(arr_dims)) dimension(s) not N=$N"))
    return read!(new_array(T, arr_dims), hdu; null, anynull, first, step, last)
end

"""
    read!(arr::DenseArray{<:Number}, hdu::FitsImageHDU, inds...) -> arr

overwrites `arr` with a rectangular sub-region of the FITS image extension in
`hdu` defined by the indices `inds...`. The destination array must have the
same dimensions as the rectangular sub-region to read (the same rules as for
sub-indexing an array are applied to determine the dimensions of the
sub-region). Keywords `anynull` and `null` may be specified to deal with
undefined pixel values.

"""
function read!(arr::DenseArray{<:Number}, hdu::FitsImageHDU,
               inds::SubArrayIndex...; kwds...)
    return read!(arr, hdu, inds; kwds...)
end

function read!(arr::DenseArray{T,L}, hdu::FitsImageHDU{<:Any,N},
               inds::SubArrayIndices{N};
               null::Union{DenseArray{Bool,L},Ref{T},Nothing} = nothing,
               anynull::Union{Nothing,Ref{Bool}} = nothing) where {T<:Number,L,N}
    img_dims = get_img_size(hdu)
    arr_dims, first, step, last = subarray_params(img_dims, inds)
    size(arr) == arr_dims || throw(DimensionMismatch("output array has invalid dimensions"))
    return read!(arr, hdu; null, anynull, first, step, last)
end

"""
    read!(arr::DenseArray{<:Number}, hdu::FitsImageHDU) -> arr

overwrites all the elements of the array `arr` with pixel values read from the
FITS image extension of the header data unit `hdu` and returns `arr`. Pixel
values are converted to the element type of `arr`. The default is to read the
complete image. This behavior may be changed by specifying another value than
`nothing` for the keywords `first` and/or `last`:

* To read a rectangular sub-image, set keywords `first` and `last` with
  `N`-tuple of integers indicating the coordinates of the first and last pixels
  to read. Optionally, keyword `step` may be set to an `N`-tuple of integers to
  indicate the increment along each dimensions. Increments must be positive.
  Here `N` is the number of dimensions of the FITS image extension.

* To read consecutive pixels, specify at least one of the keywords `first`
  and/or `last` with an integer indicating the index of the first and/or last
  pixels to read.

When at least one of the keywords `first` and/or `last` is not `nothing`, the
dimensions of the destination `arr` are not considered. In any case, the number
of elements of `arr` must be equal to the number of pixels to read.

Keyword `anynull` may be specified with a reference to a boolean
(`Ref{Bool}()`) to retrieve whether any of the read pixels is undefined.

Keyword `null` may be specified with a reference to a value of the same type as
the elements of the destination `arr` (`Ref{eltype(arr)}()`) to retrieve the
value of undefined pixels. Unless reading a rectangular sub-image, keyword
`null` may be set with an array of `Bool` of same size as `arr` and which will
be set to `true` for undefined pixels and to `false` elsewhere.

Output arrays `arr` and `null` must have contiguous elements, in other words,
they must be *dense arrays*.

"""
function read!(arr::DenseArray{T,L},
               hdu::FitsImageHDU{<:Any,N};
               null::Union{DenseArray{Bool,L},Ref{T},Nothing} = nothing,
               anynull::Union{Nothing,Ref{Bool}} = nothing,
               first::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
               last::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
               step::Union{NTuple{N,Integer},Nothing} = nothing) where {T<:Number,L,N}
    type = pixeltype_to_code(T) # clash if unsupported pixel type
    dims = get_img_size(hdu)
    len = length(arr)
    anynul = Ref{Cint}(anynull isa Ref{Bool} && len > 0)
    if first isa Tuple && last isa Tuple && null isa Union{Ref{T},Nothing}
        # Read a rectangular sub-image. NOTE: First and last pxiels must both
        # be specified because it is not possible to guess one given the other
        # and the size of the array to read as it may have fewer dimensions
        # than the image.
        if step === nothing
            ipix = NTuple(i -> one(Clong), Val(N))::NTuple{N,Clong}
        else
            ipix = as(NTuple{N,Clong}, step)
        end
        fpix = as(NTuple{N,Clong}, first)
        lpix = as(NTuple{N,Clong}, last)
        npix = 1
        for k in 1:N
            let r = fpix[k]:ipix[k]:lpix[k]
                is_valid_subindex(dims[k], r) || bad_argument("out of range rectangular sub-image")
                npix *= Int(length(r))::Int
            end
        end
        npix == len || bad_argument("rectangular sub-image and output array have different lengths")
        if len > 0
            _null = (null === nothing ? C_NULL : null)
            check(CFITSIO.fits_read_subset(
                hdu, type, Ref(fpix), Ref(lpix), Ref(ipix), _null, arr, anynul, Ref{Status}(0)))
        end
    elseif first isa Union{Integer,Nothing} && last isa Union{Integer,Nothing} && step isa Nothing
        # Read a range of consecutive pixels. NOTE: The indices of the first
        # and last pixel to write can be both optional in this case.
        if first === nothing
            if last === nothing
                size(arr) == dims || throw(DimensionMismatch(
                    "image extension and output array have differente sizes"))
                fpix = 1
            else first === nothing
                fpix = Int(fpix, len + 1 - last)::Int
            end
        else
            last === nothing || last + 1 - first == len || throw(DimensionMismatch(
                "specified pixel range and output array have differente lengths"))
            fpix = Int(first)::Int
        end
        1 ≤ fpix && fpix - 1 + len ≤ prod(dims) || bad_argument("out of range interval of pixels")
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

#------------------------------------------------------------------------------
# WRITING FITS IMAGES

"""
    write(file::FitsFile, FitsImageHDU{T}, dims...=()) -> hdu
    write(file::FitsFile, FitsImageHDU, T::Type=UInt8, dims...=()) -> hdu
    write(file::FitsFile, FitsImageHDU, bitpix::Integer, dims=()) -> hdu

create a new primary array or image extension in FITS file `file` with a
specified pixel type `T` and size `dims...`. If the FITS file is currently
empty then a primary array is created, otherwise a new image extension is
appended to the file. Pixel type can be specified as a numeric type `T` or as
an integer BITPIX code `bitpix`.

An object to manage the new extension is returned which can be used to push
header cards and then to write the data.

For example:

    hdu = write(file, FitsImageHDU, eltype(arr), size(arr))
    hdu["KEY1"] = val1             # add a 1st header record
    hdu["KEY2"] = (val2, str2)     # add a 2nd header record
    hdu["KEY3"] = (nothing, str3)  # add a 3rd header record
    write(hdu, arr)                # write data

will create a new Header Data Unit (HDU) storing array `arr` with 3 additional
header records: one named `"KEY1"` with value `val1` and no comments, another
named `"KEY2"` with value `val2` and comment string `str2`, and yet another one
named `"KEY3"` with no value and with comment string `str3`. Note that special
names `"COMMENT"`, `"HISTORY"`, and `""` indicating commentary entries have no
associated, only a comment string, say `str` which can be specified as `str` or
as `(,str)`.

"""
function write(file::FitsFile, ::Type{FitsImageHDU{T}},
               dims::NTuple{N,Integer} = ()) where {T,N}
    # NOTE: All variants end up calling this type-stable version.
    check(CFITSIO.fits_create_img(file, type_to_bitpix(T), N,
                                  Ref(as(NTuple{N,Clong}, dims)),
                                  Ref{Status}(0)))
    # The number of HDUs as returned by fits_get_num_hdus is only incremented
    # after writing data.
    n = position(file)
    if length(file) < n
        setfield!(file, :nhdus, n)
    end
    return FitsImageHDU{T,N}(BareBuild(), file, n)
end

function write(file::FitsFile, ::Type{FitsImageHDU},
               ::Type{T} = UInt8, dims::Union{Tuple{Vararg{Integer}},
                                              AbstractVector{<:Integer}} = ()) where {T}
    return write(file, FitsImageHDU{T}, dims)
end

function write(file::FitsFile, ::Type{FitsImageHDU{T,N}},
               dims::Union{Tuple{Vararg{Integer}},
                           AbstractVector{<:Integer}} = ()) where {T,N}
    length(dims) == N || throw(DimensionMismatch("incompatible number of dimensions"))
    return write(file, FitsImageHDU{T}, dims)
end

# Just convert bitpix to type.
function write(file::FitsFile, ::Type{FitsImageHDU}, bitpix::Integer,
               dims::Union{Tuple{Vararg{Integer}},
                           AbstractVector{<:Integer}} = ())
    return write(file, FitsImageHDU, type_from_bitpix(bitpix), dims)
end

# Just pack the dimensions.
function write(file::FitsFile, ::Type{FitsImageHDU{T}},
               dims::Integer...) where {T}
    return write(file, FitsImageHDU{T}, dims)
end
function write(file::FitsFile, ::Type{FitsImageHDU}, ::Type{T},
               dims::Integer...) where {T}
    return write(file, FitsImageHDU, T, dims)
end

# Just convert the dimensions.
function write(file::FitsFile, ::Type{FitsImageHDU{T}},
               dims::AbstractVector{<:Integer}) where {T}
    off = firstindex(dims) - 1
    return write(file, FitsImageHDU{T}, ntuple(i -> Clong(dims[i+off]), length(dims)))
end
function write(file::FitsFile, ::Type{FitsImageHDU}, T::Union{Integer,Type},
               dims::AbstractVector{<:Integer})
    off = firstindex(dims) - 1
    return write(file, FitsImageHDU, T, ntuple(i -> Clong(dims[i+off]), length(dims)))
end

"""
    write(file::FitsFile, hdr, arr::AbstractArray, args...) -> file

    write(file::FitsFile, arr::AbstractArray, hdr=nothing) -> file

"""
function write(file::FitsFile, hdr::Union{Header,Nothing},
               arr::AbstractArray{T,N}) where {T<:Number,N}
    # FIXME: improve type-stability
    write(merge!(write(file, FitsImageHDU, T, size(arr)), hdr), arr)
    return file
end

function write(file::FitsFile, arr::AbstractArray{<:Number},
               hdr::Union{Header,Nothing} = nothing)
    return write(file, hdr, arr)
end

function write(file::FitsFile, hdr::Union{Header,Nothing},
               arr::AbstractArray{<:Number}, args...)
    write(file, hdr, arr)
    return write(file, args...)
end

function write(file::FitsFile, arr::AbstractArray{<:Number},
               hdr::Union{Header,Nothing}, args...)
    write(file, hdr, arr)
    return write(file, args...)
end

"""
    write(hdu::FitsImageHDU, arr::AbstractArray{<:Number}) -> hdu

writes all the elements of the array `arr` to the pixels of the FITS image
extension of the header data unit `hdu`. The element of `arr` are converted to
the type of the pixels. The default is to write the complete image pixels, so
the size of the destination `arr` must be identical to the dimensions of the FITS
image extension. This behavior may be changed by specifying another value than
`nothing` for the keywords `first` and/or `last`:

* To write a rectangular sub-image, specify keywords `first` and `last` with
  the coordinates of the first and last pixels to write as an `N`-tuple of
  integers, with `N` the number of dimensions of the FITS image extension.

* To write consecutive pixels, specify at least one of the keywords `first`
  and/or `last` with the index of the first and/or last pixels to write as an
  integer.

When at least one of the keywords `first` and/or `last` is not `nothing`, the
dimensions of `arr` are not considered. In any case, the number of elements of
`arr` must be equal to the number of pixels to write.

Unless writing a rectangular sub-image, keyword `null` may be used to specify
the value of undefined elements in `arr`. For integer FITS images, the FITS null
value is defined by the BLANK keyword (an error is returned if the BLANK
keyword doesn't exist). For floating point FITS images the special IEEE NaN
(Not-a-Number) value will be written into the FITS file.

"""
function write(hdu::FitsImageHDU{<:Any,N},
               arr::AbstractArray{T,L};
               first::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
               last::Union{Integer,NTuple{N,Integer},Nothing} = nothing,
               null::Union{Number,Nothing} = nothing) where {T<:Number,L,N}
    type = pixeltype_to_code(arr) # clash if unsupported pixel type
    dims = get_img_size(hdu)
    len = length(arr)
    if first isa Tuple && last isa Tuple && null isa Nothing
        # Write a rectangular sub-image. NOTE: First and last pxiels must both
        # be specified because it is not possible to guess one given the other
        # and the size of the array to write as it may have fewer dimensions
        # than the image.
        fpix = as(NTuple{N,Clong}, first)
        lpix = as(NTuple{N,Clong}, last)
        npix = 1
        for k in 1:N
            let r = fpix[k]:lpix[k]
                is_valid_subindex(dims[k], r) || bad_argument("out of range rectangular sub-image")
                npix *= Int(length(r))::Int
            end
        end
        npix == len || bad_argument("rectangular sub-image and input array have different lengths")
        if len > 0
            check(CFITSIO.fits_write_subset(
                hdu, type, Ref(fpix), Ref(lpix), dense_array(arr), Ref{Status}(0)))
        end
    elseif first isa Union{Integer,Nothing}
        # Write a range of consecutive pixels. NOTE: The indices of the first
        # and last pixel to write can be both optional in this case.
        if first === nothing
            if last === nothing
                size(arr) == dims || throw(DimensionMismatch(
                    "image extension and input array have differente sizes"))
                fpix = 1
            else first === nothing
                fpix = Int(fpix, len + 1 - last)::Int
            end
        else
            last === nothing || last + 1 - first == len || throw(DimensionMismatch(
                "specified pixel range and input array have differente lengths"))
            fpix = Int(first)::Int
        end
        1 ≤ fpix && fpix - 1 + len ≤ prod(dims) || bad_argument("out of range interval of pixels")
        if len > 0
            if null === nothing
                check(CFITSIO.fits_write_img(
                    hdu, type, fpix, len, dense_array(arr), Ref{Status}(0)))
            else
                check(CFITSIO.fits_write_imgnull(
                    hdu, type, fpix, len, dense_array(arr), Ref{eltype(arr)}(null),
                    Ref{Status}(0)))
            end
        end
    else
        error("invalid combination of options")
    end
    return hdu
end
