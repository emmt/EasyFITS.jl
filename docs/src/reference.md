# Reference

## Types

```@docs
FitsHDU
FitsError
FitsLogic
```

## Aliases

```@docs
EasyFITS.Header
EasyFITS.TableData
```

## Direct reading/writing of FITS files

```@docs
readfits
readfits!
writefits
writefits!
```

## FITS files

Missing `EasyFITS.write!`, `read(::FitsFile, ...)`, `read!(::DenseArray{<:Number},::FitsImageHDU,::SubArrayIndex...)`
`read!(::DenseArray{<:Number},::FitsImageHDU)`


```@docs
FitsFile
openfits
pathof(::FitsFile)
filemode(::FitsFile)
isreadonly(::FitsFile)
iswritable(::FitsFile)
seek(::FitsFile, ::Integer)
seekstart(::FitsFile)
seekend(::FitsFile)
position(::FitsFile)
flush(::FitsFile)
eachmatch(::Any, ::FitsFile)
```

## FITS Image HDUs

```@docs
FitsImageHDU
read(::FitsImageHDU)
read!(::Array, ::FitsImageHDU)
write(::FitsImageHDU, ::Array)
```

## FITS Table HDUs

```@docs
FitsTableHDU
read(::FitsTableHDU)
read(::Type{Vector}, ::FitsTableHDU)
read(::FitsTableHDU, ::String)
read!(::Dict, ::FitsTableHDU)
read!(::Array, ::FitsTableHDU, ::String)
```

## FITS header

```@docs
FitsHeader(::FitsHDU)
```

## FITS Header Data Units (HDUs)

```@docs
nameof(::FitsHDU)
EasyFITS.is_named
```

## Utilities

```@docs
EasyFITS.OutputCstring
EasyFITS.cfitsio_errmsg
```
