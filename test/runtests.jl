module EasyFITSTests

using Test
using FITSIO
using EasyFITS

function samevalues(A::AbstractArray, B::AbstractArray)
    @assert !Base.has_offset_axes(A, B)
    @assert size(A) == size(B)
    for i ∈ eachindex(A, B)
        A[i] != B[i] && return false
    end
    return true
end

function maxabsdif(A::AbstractArray, B::AbstractArray)
    @assert !Base.has_offset_axes(A, B)
    @assert size(A) == size(B)
    T = promote_type(eltype(A), eltype(B))
    res = zero(T)
    for i ∈ eachindex(A, B)
        res = max(res, abs(A[i] - B[i]))
    end
    return res
end

dat1 = rand(Float32, 3, 4, 5)
hdr1 = EasyFITS.header()
setfitskey!(hdr1, "HDUNAME", "HDU-ONE", "First custom HDU")
setfitskey!(hdr1, "GA",  42,       "Some integer keyword")
setfitskey!(hdr1, "BU",  "Shadok", "Some string keyword")
setfitskey!(hdr1, "ZO",  3.1415,   "Some real keyword")
setfitskey!(hdr1, "MEU", true,     "Some boolean keyword")

dat2 = rand(Int32, 7, 8)
hdr2 = EasyFITS.header()
setfitskey!(hdr2, "HDUNAME", "HDU-TWO", "Second custom HDU")
setfitskey!(hdr2, "GA",  -42,       "Some integer keyword")
setfitskey!(hdr2, "BU",  "Gibi",    "Some string keyword")
setfitskey!(hdr2, "ZO",  sqrt(2),   "Some real keyword")
setfitskey!(hdr2, "MEU", false,     "Some boolean keyword")

path = "test.fits"
createfits!(path) do io
    write(io, dat1; header=hdr1)
    write(io, dat2; header=hdr2)
end

@testset "Low-level" begin
    @test_throws ErrorException close(createfits(path; overwrite=false))

    openfits(path) do io
        @test EasyFITS.find(hdu -> EasyFITS.hduname(hdu) == "HDU-THREE", io) === nothing
        @test EasyFITS.findfirst(hdu -> EasyFITS.hduname(hdu) == "HDU-THREE", io) === nothing
        @test EasyFITS.findlast(hdu -> EasyFITS.hduname(hdu) == "HDU-THREE", io) === nothing
        @test EasyFITS.findnext(hdu -> EasyFITS.hduname(hdu) == "HDU-THREE", io, 2) === nothing
        @test EasyFITS.findprev(hdu -> EasyFITS.hduname(hdu) == "HDU-THREE", io, 2) === nothing
        @test EasyFITS.find(hdu -> EasyFITS.hduname(hdu) == "HDU-ONE", io) == 1
        @test EasyFITS.find(hdu -> EasyFITS.hduname(hdu) == "HDU-TWO", io) == 2
        @test EasyFITS.findfirst(hdu -> EasyFITS.hduname(hdu) == "HDU-ONE", io) == 1
        @test EasyFITS.findlast(hdu -> EasyFITS.hduname(hdu) == "HDU-ONE", io) == 1
        @test EasyFITS.findnext(hdu -> EasyFITS.hduname(hdu) == "HDU-ONE", io, 2) === nothing
        @test EasyFITS.findprev(hdu -> EasyFITS.hduname(hdu) == "HDU-ONE", io, 2) === 1
        hdu = io[1]
        @test tryreadfitskey(hdu, Int,     "GA") == 42
        @test tryreadfitskey(hdu, String,  "BU") == "Shadok"
        @test tryreadfitskey(hdu, Float64, "ZO") ≈ 3.1415
        @test tryreadfitskey(hdu, Bool,    "MEU") == true
        dat = read(hdu)
        @test eltype(dat) == eltype(dat1)
        @test size(dat) == size(dat1)
        @test samevalues(dat, dat1)
        hdu = io[2]
        @test tryreadfitskey(hdu, Int,     "GA") == -42
        @test tryreadfitskey(hdu, String,  "BU") == "Gibi"
        @test tryreadfitskey(hdu, Float64, "ZO") ≈ sqrt(2)
        @test tryreadfitskey(hdu, Bool,    "MEU") == false
        dat = read(hdu)
        @test eltype(dat) == eltype(dat2)
        @test size(dat) == size(dat2)
        @test samevalues(dat, dat2)
    end
end
@testset "High-level" begin
    A1 = readfits(path)
    @test isa(A1["GA"],  Int)     && A1["GA"] == 42
    @test isa(A1["BU"],  String)  && A1["BU"] == "Shadok"
    @test isa(A1["ZO"],  Float64) && A1["ZO"] ≈ 3.1415
    @test isa(A1["MEU"], Bool)    && A1["MEU"] == true
    @test A1.GA  == A1["GA"]
    @test A1.BU  == A1["BU"]
    @test A1.ZO  == A1["ZO"]
    @test A1.MEU == A1["MEU"]
    A1["ZO"] = 15
    @test isa(A1["ZO"],  Int) && A1["ZO"] == 15
    @test A1.ZO  == A1["ZO"]
    A1.ZO = 172
    @test isa(A1["ZO"],  Int) && A1["ZO"] == 172
    @test A1.ZO  == A1["ZO"]
    @test getfitscomment(A1, "ZO") == "Some real keyword"
    #
    setfitskey!(A1, "ZO", 21, "New comment")
    @test isa(A1["ZO"],  Int) && A1["ZO"] == 21
    @test getfitscomment(A1, "ZO") == "New comment"
    #
    A1["GA"] = 127, "Another comment"
    @test A1["GA"] == 127 && getfitscomment(A1, "GA") == "Another comment"
    A1["GA"] = 128, nothing
    @test A1["GA"] == 128 && getfitscomment(A1, "GA") == ""
    A1["GA"] = (129, )
    @test A1["GA"] == 129 && getfitscomment(A1, "GA") == ""
    #
    A1.ZO = 220.0, "Yet another comment"
    @test A1.ZO == 220 && getfitscomment(A1, "ZO") == "Yet another comment"
    A1.ZO = 221.0, nothing
    @test A1.ZO == 221 && getfitscomment(A1, "ZO") == ""
    A1.ZO = (223.0, )
    @test A1.ZO == 223 && getfitscomment(A1, "ZO") == ""
    #
    @test eltype(A1) == eltype(dat1)
    @test ndims(A1) == ndims(dat1)
    @test size(A1) == size(dat1)
    @test samevalues(A1, dat1)
    A2 = readfits(path, 2)
    @test isa(A2["GA"],  Int)     && A2["GA"] == -42
    @test isa(A2["BU"],  String)  && A2["BU"] == "Gibi"
    @test isa(A2["ZO"],  Float64) && A2["ZO"] ≈ sqrt(2)
    @test isa(A2["MEU"], Bool)    && A2["MEU"] == false
    @test eltype(A2) == eltype(dat2)
    @test ndims(A2) == ndims(dat2)
    @test size(A2) == size(dat2)
    @test samevalues(A2, dat2)
end

#rm(path)

end # module
