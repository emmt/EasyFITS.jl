# Reference

## ltoh
## ntoh
## readchomp
## read!

```@docs
EasyFITS.Header
EasyFITS.TableData
FitsHDU
FitsLogic
```

## FITS Files

```@docs
FitsFile
openfits
open(::Type{FitsFile}, ::AbstractString)
open(::Function, ::Type{FitsFile})
readfits
readfits!
read(::Type{FitsFile}, ::AbstractString)
write!
writefits
writefits!
isopen(::FitsFile)
close(::FitsFile)
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

## FITS Header Data Units (HDUs)

```@docs
nameof(::FitsHDU)
EasyFITS.is_named
```
