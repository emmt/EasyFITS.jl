module EasyFITSTests

using Test
using FITSIO
using EasyFITS
using EasyFITS: find

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

path = joinpath(@__DIR__, "test.fits")

dat1 = rand(Float32, 3, 4, 5)
hdr1 = FitsHeader()

dat2 = rand(Int32, 7, 8)
hdr2 = FitsHeader(
    "HDUNAME" => ("HDU-TWO", "Second custom HDU"),
    "GA"      => (-42,       "Some integer keyword"),
    "BU"      => ("Gibi",    "Some string keyword"),
    "ZO"      => (sqrt(2),   "Some real keyword"),
    "MEU"     => (false,     "Some boolean keyword"))

dat3 = rand(Int32, 11)
img3 = FitsImage(dat3,
                 HDUNAME = ("HDU-THREE", "Third custom HDU"),
                 GA      = (-42,         "Some integer keyword"),
                 BU      = ("Gibi",      "Some string keyword"),
                 ZO      = (sqrt(2),     "Some real keyword"),
                 MEU     = (false,       "Some boolean keyword"))

@testset "Create FITS File" begin
    @test length(hdr1) == 0
    hdr1["HDUNAME"] = ("HDU-ONE", "First custom HDU")
    hdr1["GA"]      = (42,        "Some integer keyword")
    hdr1["BU"]      = ("Shadok",  "Some string keyword")
    hdr1["ZO"]      = (3.1415,    "Some real keyword")
    hdr1["MEU"]     = (true,      "Some boolean keyword")

    #FitsIO(path, "w!") do io
    #    write(io, dat1, hdr1)
    #    write(io, dat2, hdr2)
    #    write(io, img3)
    #end
    write(FitsFile, path, (dat1, hdr1), (dat2, hdr2), img3; overwrite=true)

    # Check overwrite is forbidden.
    @test_throws ErrorException close(FitsIO(path, "w"))
end


hduname(obj::Union{FitsHDU,FitsIO}) =
    get(String, obj, "HDUNAME", nothing)

@testset "Low-level" begin
    FitsIO(path) do io
        @test find(hdu -> hduname(hdu) == "HDU-MISSING", io) === nothing
        @test findfirst(hdu -> hduname(hdu) == "HDU-MISSING", io) === nothing
        @test findlast(hdu -> hduname(hdu) == "HDU-MISSING", io) === nothing
        @test findnext(hdu -> hduname(hdu) == "HDU-MISSING", io, 2) === nothing
        @test findprev(hdu -> hduname(hdu) == "HDU-MISSING", io, 2) === nothing
        @test find(hdu -> hduname(hdu) == "HDU-ONE", io) == 1
        @test find(hdu -> hduname(hdu) == "HDU-TWO", io) == 2
        @test findfirst(hdu -> hduname(hdu) == "HDU-ONE", io) == 1
        @test findlast(hdu -> hduname(hdu) == "HDU-ONE", io) == 1
        @test findnext(hdu -> hduname(hdu) == "HDU-ONE", io, 2) === nothing
        @test findprev(hdu -> hduname(hdu) == "HDU-ONE", io, 2) === 1
        hdu = io[1]
        @test get(Int,     hdu, "GA") == 42
        @test get(String,  hdu, "BU") == "Shadok"
        @test get(Float64, hdu, "ZO") ≈ 3.1415
        @test get(Bool,    hdu, "MEU") == true
        dat = read(hdu)
        @test eltype(dat) == eltype(dat1)
        @test size(dat) == size(dat1)
        @test samevalues(dat, dat1)
        hdu = io[2]
        @test get(Int,     hdu, "GA") == -42
        @test get(String,  hdu, "BU") == "Gibi"
        @test get(Float64, hdu, "ZO") ≈ sqrt(2)
        @test get(Bool,    hdu, "MEU") == false
        dat = read(hdu)
        @test eltype(dat) == eltype(dat2)
        @test size(dat) == size(dat2)
        @test samevalues(dat, dat2)
    end
end

@testset "High-level" begin
    # Read array+header data and check indexing of array and header.
    A1 = read(FitsImage, path)
    @test isa(A1, FitsImage)
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
    @test get(FitsComment, A1, "ZO") == "Some real keyword"
    #
    A1["ZO"] = (21, "New comment")
    @test isa(A1["ZO"],  Int) && A1["ZO"] == 21
    @test get(FitsComment, A1, "ZO") == "New comment"
    #
    A1["GA"] = 127, "Another comment"
    @test A1["GA"] == 127 && get(FitsComment, A1, "GA") == "Another comment"
    A1["GA"] = 128, nothing
    @test A1["GA"] == 128 && get(FitsComment, A1, "GA") == ""
    A1["GA"] = (129, )
    @test A1["GA"] == 129 && get(FitsComment, A1, "GA") == ""
    #
    A1.ZO = 220.0, "Yet another comment"
    @test A1.ZO == 220 && get(FitsComment, A1, "ZO") == "Yet another comment"
    A1.ZO = 221.0, nothing
    @test A1.ZO == 221 && get(FitsComment, A1, "ZO") == ""
    A1.ZO = (223.0, )
    @test A1.ZO == 223 && get(FitsComment, A1, "ZO") == ""
    #
    @test eltype(A1) == eltype(dat1)
    @test ndims(A1) == ndims(dat1)
    @test size(A1) == size(dat1)
    @test samevalues(A1, dat1)
    A2 = read(FitsImage, path, 2)
    @test isa(A2, FitsImage)
    @test isa(A2["GA"],  Int)     && A2["GA"] == -42
    @test isa(A2["BU"],  String)  && A2["BU"] == "Gibi"
    @test isa(A2["ZO"],  Float64) && A2["ZO"] ≈ sqrt(2)
    @test isa(A2["MEU"], Bool)    && A2["MEU"] == false
    @test eltype(A2) == eltype(dat2)
    @test ndims(A2) == ndims(dat2)
    @test size(A2) == size(dat2)
    @test samevalues(A2, dat2)
    # Read headers.
    H1 = read(FitsHeader, path, 1)
    @test isa(H1, FitsHeader)
    @test isa(H1["GA"],  Int)     && H1["GA"] == 42
    @test isa(H1["BU"],  String)  && H1["BU"] == "Shadok"
    @test isa(H1["ZO"],  Float64) && H1["ZO"] ≈ 3.1415
    @test isa(H1["MEU"], Bool)    && H1["MEU"] == true
    @test H1.GA  == H1["GA"]
    @test H1.BU  == H1["BU"]
    @test H1.ZO  == H1["ZO"]
    @test H1.MEU == H1["MEU"]
    H1["ZO"] = 15
    @test isa(H1["ZO"],  Int) && H1["ZO"] == 15
    @test H1.ZO  == H1["ZO"]
    H1.ZO = 172
    @test isa(H1["ZO"],  Int) && H1["ZO"] == 172
    @test H1.ZO  == H1["ZO"]
    @test get(FitsComment, H1, "ZO") == "Some real keyword"
    #
    H2 = read(FitsHeader, path, 2)
    @test isa(H2, FitsHeader)
    # Read array data with contraints.
    A3 = read(FitsArray, path, 1)
    @test isa(A3, Array{eltype(dat1),ndims(dat1)})
    @test size(A3) == size(dat1)
    @test samevalues(A3, dat1)
    A4 = read(FitsArray{Int}, path, 2)
    @test isa(A4, Array{Int,ndims(dat2)})
    @test size(A4) == size(dat2)
    @test samevalues(A4, dat2)
    @test_throws MethodError read(FitsArray{Int,ndims(dat2)+1}, path, 2)
    A5 = read(FitsArray{Int,ndims(dat2)}, path, 2)
    @test isa(A5, Array{Int,ndims(dat2)})
    @test size(A5) == size(dat2)
    @test samevalues(A5, dat2)
    # Read array+header data with contraints.
    A6 = read(FitsImage, path, 1)
    @test isa(A6, FitsImage{eltype(dat1),ndims(dat1)})
    @test size(A6) == size(dat1)
    @test samevalues(A6, dat1)
    A7 = read(FitsImage{Float64}, path, 1)
    @test isa(A7, FitsImage{Float64,ndims(dat1)})
    @test size(A7) == size(dat1)
    @test samevalues(A7, dat1)
    @test_throws MethodError read(FitsImage{Int,ndims(dat2)+1}, path, 2)
    A8 = read(FitsImage{Int,ndims(dat2)}, path, 2)
    @test isa(A8, FitsImage{Int,ndims(dat2)})
    @test size(A8) == size(dat2)
    @test samevalues(A8, dat2)

end

#rm(path)

end # module
