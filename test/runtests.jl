module TestingEasyFITS

using Test, EasyFITS, Dates, DataFrames

using EasyFITS: FitsInteger, FitsFloat, FitsComplex

column_values(arg::AbstractArray) = arg
column_values(arg::Tuple{AbstractArray,AbstractString}) = arg[1]

# Yield type for which we are sure that no possible conversion is implemented.
other_type(type::FitsCardType) =
    type === FITS_LOGICAL   ? Missing :
    type === FITS_INTEGER   ? Missing :
    type === FITS_FLOAT     ? Missing :
    type === FITS_COMPLEX   ? Missing :
    type === FITS_STRING    ? Missing :
    type === FITS_COMMENT   ? Missing : Int

@testset "BITPIX" begin
    let type_to_bitpix = EasyFITS.type_to_bitpix,
        type_from_bitpix = EasyFITS.type_from_bitpix
        # Standard BITPIX types.
        @test type_to_bitpix(UInt8)   ===   8
        @test type_to_bitpix(Int16)   ===  16
        @test type_to_bitpix(Int32)   ===  32
        @test type_to_bitpix(Int64)   ===  64
        @test type_to_bitpix(Float32) === -32
        @test type_to_bitpix(Float64) === -64
        @test_throws MethodError type_to_bitpix(String)
        @test type_from_bitpix(  8) === UInt8
        @test type_from_bitpix( 16) === Int16
        @test type_from_bitpix( 32) === Int32
        @test type_from_bitpix( 64) === Int64
        @test type_from_bitpix(-32) === Float32
        @test type_from_bitpix(-64) === Float64
        @test_throws ArgumentError type_from_bitpix(typemin(Int))
        @test type_to_bitpix(Int16[]) === 16

        # Non-standard pixel types.
        @test type_from_bitpix(type_to_bitpix(Int8)) === Int8
        @test type_from_bitpix(type_to_bitpix(UInt16)) === UInt16
        @test type_from_bitpix(type_to_bitpix(UInt32)) === UInt32
        @test type_from_bitpix(type_to_bitpix(UInt64)) === UInt64
    end
end

@testset "Data types" begin
    let type_to_code = EasyFITS.type_to_code,
        type_from_code = EasyFITS.type_from_code
        @test type_from_code(type_to_code(AbstractString)) === String
        @test type_from_code(type_to_code(SubString)) === String
        @test type_from_code(type_to_code(String)) === String
        @test type_from_code(type_to_code(Bool)) === Bool
        @test type_from_code(type_to_code(Int8)) === Int8
        @test type_from_code(type_to_code(UInt8)) === UInt8
        @test type_from_code(type_to_code(Int16)) === Int16
        @test type_from_code(type_to_code(UInt16)) === UInt16
        @test type_from_code(type_to_code(Int32)) === Int32
        @test type_from_code(type_to_code(UInt32)) === UInt32
        @test type_from_code(type_to_code(Int64)) === Int64
        @test type_from_code(type_to_code(UInt64)) === UInt64
        @test type_from_code(type_to_code(Float32)) === Float32
        @test type_from_code(type_to_code(Float64)) === Float64
        @test type_from_code(type_to_code(ComplexF32)) === ComplexF32
        @test type_from_code(type_to_code(ComplexF64)) === ComplexF64
        @test_throws MethodError type_to_code(Dict)
        @test_throws ArgumentError type_from_code(typemin(Int))
    end
end

@testset "TFORM types" begin
    let Bit = EasyFITS.Bit,
        type_to_letter = EasyFITS.type_to_letter,
        type_from_letter = EasyFITS.type_from_letter
        @test type_to_letter(AbstractString) === 'A'
        @test type_to_letter(SubString)      === 'A'
        @test type_to_letter(String)         === 'A'
        @test type_to_letter(Bool)           === 'L'
        @test type_to_letter(UInt8)          === 'B'
        @test type_to_letter(Int16)          === 'I'
        @test type_to_letter(Int32)          === 'J'
        @test type_to_letter(Int64)          === 'K'
        @test type_to_letter(Float32)        === 'E'
        @test type_to_letter(Float64)        === 'D'
        @test type_to_letter(ComplexF32)     === 'C'
        @test type_to_letter(ComplexF64)     === 'M'
        @test type_to_letter(Bit)            === 'X'
        @test type_to_letter(Int8)           === 'S'
        @test type_to_letter(UInt16)         === 'U'
        @test type_to_letter(UInt32)         === 'V'
        @test type_to_letter(UInt64)         === 'W'
        @test_throws MethodError type_to_letter(Dict)
        @test type_to_letter(Bool[])         === 'L'

        @test type_from_letter('A') === String
        @test type_from_letter('L') === Bool
        @test type_from_letter('B') === UInt8
        @test type_from_letter('I') === Int16
        @test type_from_letter('J') === Int32
        @test type_from_letter('K') === Int64
        @test type_from_letter('E') === Float32
        @test type_from_letter('D') === Float64
        @test type_from_letter('C') === ComplexF32
        @test type_from_letter('M') === ComplexF64
        @test type_from_letter('X') === Bit
        @test type_from_letter('S') === Int8
        @test type_from_letter('U') === UInt16
        @test type_from_letter('V') === UInt32
        @test type_from_letter('W') === UInt64
        @test_throws ArgumentError type_from_letter('@')
        @test type_from_letter(Int('L')) === Bool

        @test type_from_letter(type_to_letter(AbstractString)) === String
        @test type_from_letter(type_to_letter(SubString)) === String
        @test type_from_letter(type_to_letter(String)) === String
        @test type_from_letter(type_to_letter(Bool)) === Bool
        @test type_from_letter(type_to_letter(Bit)) === Bit
        @test type_from_letter(type_to_letter(Int8)) === Int8
        @test type_from_letter(type_to_letter(UInt8)) === UInt8
        @test type_from_letter(type_to_letter(Int16)) === Int16
        @test type_from_letter(type_to_letter(UInt16)) === UInt16
        @test type_from_letter(type_to_letter(Int32)) === Int32
        @test type_from_letter(type_to_letter(UInt32)) === UInt32
        @test type_from_letter(type_to_letter(Int64)) === Int64
        @test type_from_letter(type_to_letter(UInt64)) === UInt64
        @test type_from_letter(type_to_letter(Float32)) === Float32
        @test type_from_letter(type_to_letter(Float64)) === Float64
        @test type_from_letter(type_to_letter(ComplexF32)) === ComplexF32
        @test type_from_letter(type_to_letter(ComplexF64)) === ComplexF64
    end
end

@testset "Utilities" begin
    @test EasyFITS.CFITSIO_VERSION isa VersionNumber
    let new_array = EasyFITS.new_array
        let A = new_array(Float32, 4, 5, 6)
            @test eltype(A) === Float32
            @test size(A) === (4, 5, 6)
        end
        let A = new_array(Int16, (4, 5, 6))
            @test eltype(A) === Int16
            @test size(A) === (4, 5, 6)
        end
        let A = new_array(UInt32, [4, 5, 6])
            @test eltype(A) === UInt32
            @test size(A) === (4, 5, 6)
        end
        let A = new_array(Float32, Val(3), 4, 5, 6)
            @test eltype(A) === Float32
            @test size(A) === (4, 5, 6)
        end
        let A = new_array(Int16, Val(3), (4, 5, 6))
            @test eltype(A) === Int16
            @test size(A) === (4, 5, 6)
        end
        let A = new_array(UInt32, Val(3), [4, 5, 6])
            @test eltype(A) === UInt32
            @test size(A) === (4, 5, 6)
        end
        @test_throws DimensionMismatch new_array(Float32, Val(2),  4, 5, 6)
        @test_throws DimensionMismatch new_array(Int16,   Val(2), (4, 5, 6))
        @test_throws DimensionMismatch new_array(UInt32,  Val(2), [4, 5, 6])
    end
    let dense_array = EasyFITS.dense_array,
        A = convert(Array, reshape(1:24, 2,3,4)),
        B = view(A, :, 1:2:3, :)
        @test dense_array(A) === A
        @test isa(dense_array(B), Array)
    end
    let string_length = EasyFITS.string_length
        @test string_length("") == 0
        @test string_length(" ") == 1
        @test string_length("  ") == 1
        @test string_length("a") == 1
        @test string_length("a ") == 1
        @test string_length("a  ") == 1
        @test string_length(" a") == 2
        @test string_length(" a ") == 2
        @test string_length(" a  ") == 2
    end
end

# Get comment as a string from header card settings.
get_comment(dat::Any) = ""
get_comment(dat::Tuple{Any}) = ""
get_comment(dat::Tuple{Any,Nothing}) = ""
get_comment(dat::Tuple{Any,AbstractString}) = dat[2]

function get_units(comment::AbstractString)
    m = match(r"^\[\s*(.*?)\s*\]", comment)
    m === nothing ? "" : m.captures[1]
end

suppress_units(comment::AbstractString) =
    replace(comment, r"^\[.*?\]\s*" => "")

cards_1 = (KEY_B1 = (true,  "This is true"),
           KEY_B2 = (false, "This is false"),
           KEY_B3 = false,
           KEY_B4 = (true, ""),
           KEY_I1 = 123, # no comment
           KEY_I2 = (0x9, nothing),
           KEY_I3 = (Int16(42), "[cm] with units"),
           KEY_I4 = (-77, "[] no units"),
           KEY_I5 = (101, "[  foo /  bar ] strip leading/trailing spaces"),
           KEY_F1 = (-1.2f-2, "Single precision"),
           KEY_F2 = (-3.7, "Double precision"),
           KEY_F3 = (11//7, "Rational number"),
           KEY_F4 = (pi, "Irrational number"),
           KEY_C1 = (complex(-2,3), "Integer complex"),
           KEY_C2 = (complex(2.1f0,-3.7f0), "Single precision complex"),
           KEY_C3 = (complex(2.1,-3.7), "Double precision complex"),
           KEY_S1 = " <- significant space here",
           KEY_S2 = ("", "Empty string"),
           KEY_S3 = (" ", "Another empty string"),
           KEY_S4 = ("'oops!'", "String with quotes"),
           KEY_S5 = " <- keep this, not these ->    ",
           KEY_S6 = (SubString("Hello world!", 7, 12), "A sub-string"),
           KEY_N2 = (undef, "undefined value"),
           KEY_N3 = (missing, "missing value"),
           COMMENT = "simple comment",
           HISTORY = "historical comment",
           VERY_LONG_KEY_NAME = (1, "Should use HIERARCH convention"),
           UNCOMMENTED_VERY_LONG_KEY_NAME = 2.0,
           )

cards_2 = [String(key) => val for (key,val) in pairs(cards_1)]
tempfile, io = mktemp(; cleanup=false)
close(io)

@testset "FITS Headers" begin
    #@test FitsCardType(Any)             === FITS_UNKNOWN
    @test FitsCardType(typeof(undef))    === FITS_UNDEFINED
    @test FitsCardType(Missing)          === FITS_UNDEFINED
    @test FitsCardType(Bool)             === FITS_LOGICAL
    @test FitsCardType(Int)              === FITS_INTEGER
    @test FitsCardType(Float32)          === FITS_FLOAT
    @test FitsCardType(Complex{Int32})   === FITS_COMPLEX
    @test FitsCardType(Complex{Float64}) === FITS_COMPLEX
    @test FitsCardType(String)           === FITS_STRING
    @test FitsCardType(Nothing)          === FITS_COMMENT
end

@testset "FITS Images" begin
    # Write a simple FITS image.
    A = convert(Array{Int16}, reshape(1:60, 3,4,5))
    openfits(tempfile, "w!") do file
        @test file isa FitsFile
        @test !isreadable(file)
        @test !isreadonly(file)
        @test iswritable(file)
        @test position(file) == 1
        @test length(file) == 0

        # Add a simple IMAGE extension.
        let hdu = write(file, FitsImageHDU, eltype(A), size(A))
            @test length(file) == 1
            @test hdu === last(file)
            @test firstindex(hdu) == 1
            @test lastindex(hdu) == length(hdu)
            @test_throws BoundsError hdu[firstindex(hdu) - 1]
            @test_throws KeyError hdu["DUMMY"]
            let card = hdu["SIMPLE"]
                @test card isa FitsCard
                @test card === hdu[firstindex(hdu)]
                @test card.key === Fits"SIMPLE"
                @test card.type === FITS_LOGICAL
                @test card.value() === card.logical
                @test card.value() === true
            end
            let card = hdu["BITPIX"]
                @test card === hdu[firstindex(hdu) + 1]
                @test card.type === FITS_INTEGER
                @test card.value() === card.integer
                @test card.value() == 16
            end
            let card = hdu["NAXIS"]
                @test card == hdu[firstindex(hdu) + 2]
                @test card.type == FITS_INTEGER
                @test card.value() === card.integer
                @test card.value() == 3
            end
            reset(hdu) # reset before testing incremental search
            for i in 1:ndims(A)+1
                local card = get(hdu, "NAXIS#", nothing)
                @test (card === nothing) == (i > ndims(A))
                card === nothing && break
                @test card.type == FITS_INTEGER
                @test card.value() === card.integer
                @test card.value() == size(A, i)
            end
            @test  haskey(hdu, "NAXIS")
            @test !haskey(hdu, "NO-SUCH-KEY")
            @test  haskey(hdu, firstindex(hdu))
            @test !haskey(hdu, lastindex(hdu) + 1)
            write(hdu, A)
        end
        @test position(file) == 1
        @test length(file) == 1
        for pass in 1:2
            # Add IMAGE extensions with no data, just header cards.
            local cards = pass == 1 ? cards_1 : cards_2
            local hdu = write(file, FitsImageHDU)
            len = length(hdu)
            merge!(hdu, cards)
            @test length(hdu) == len + length(cards)
            let buf = IOBuffer()
                # Exercise `show`.
                for i in eachindex(hdu)
                    if pass == 1
                        show(buf, hdu[i])
                    else
                        show(buf, MIME"text/plain"(), hdu[i])
                    end
                end
                # Just test something.
                @test length(String(take!(buf))) ≥ length(hdu)
            end
            if pass == 1
                # Check exponent conversion, i.e. 'd' or 'D' -> 'E' when
                # parsing a floating-point value.
                key = "TEST"
                val = 1e-200
                hdu[key] = (val, "Quite a tiny value")
                card = hdu[key]
                @test card.type == FITS_FLOAT
                # FIXME: i = findfirst(isequal(UInt8('=')), card)
                # FIXME: @test i == 9
                # FIXME: j = findnext(b -> b ∈ (UInt8('e'), UInt8('E'), UInt8('d'), UInt8('D')), card, i)
                # FIXME: bak = card.value()
                # FIXME: @test bak ≈ val
                # FIXME: for c in "dDeE"
                # FIXME:     card[j] = c
                # FIXME:     @test card.value() ≈ val
                # FIXME:     @test card.value() === bak
                # FIXME: end
                delete!(hdu, key)
            end
            for (key, dat) in (cards isa AbstractVector ? cards : pairs(cards))
                local card = get(hdu, key, nothing)
                @test card isa FitsCard
                @test_throws ErrorException convert(other_type(card.type), card.value)
                if uppercase(rstrip(String(key))) ∈ ("HISTORY","COMMENT","")
                    local com = dat isa Tuple{Nothing,AbstractString} ? dat[2] :
                        dat isa AbstractString ? dat :
                        dat isa Nothing || dat isa Tuple{Nothing,Nothing} ? "" :
                        error("bad card specification: $(repr(key)) => $(repr(dat))")
                    @test card.type == FITS_COMMENT
                    @test card.value() === nothing
                    @test card.comment == com
                else
                    local val = dat isa Tuple ? dat[1] : dat
                    local com = dat isa Tuple{Any,AbstractString} ? dat[2] : ""
                    if val isa Bool
                        @test card.type == FITS_LOGICAL
                        @test card.logical === card.value()
                        @test card.logical === val
                    elseif val isa Integer
                        @test card.type == FITS_INTEGER
                        @test card.integer === card.value()
                        @test card.integer == val
                    elseif val isa Real
                        @test card.type == FITS_FLOAT
                        @test card.float === card.value()
                        @test card.float ≈ val
                    elseif val isa Complex
                        @test card.type == FITS_COMPLEX
                        @test card.complex === card.value()
                        @test card.complex ≈ val
                    elseif val isa Complex
                        @test card.type == FITS_COMPLEX
                        @test card.complex === card.value()
                        @test card.complex ≈ val
                    elseif val isa AbstractString
                        @test card.type == FITS_STRING
                        @test card.string === card.value()
                        @test card.string == rstrip(card.string)
                        @test card.string == rstrip(val)
                    elseif val isa Union{Missing,UndefInitializer}
                        @test card.type == FITS_UNDEFINED
                        @test card.value() === missing
                    else
                        error("bad card specification: $(repr(key)) => $(repr(dat))")
                    end
                    @test card.comment == com
                    @test card.units == get_units(com)
                    @test card.unitless == suppress_units(com)
                end
            end
        end
    end
    # Read the data.
    B = readfits(tempfile)
    @test eltype(B) == eltype(A)
    @test size(B) == size(A)
    @test B == A
    openfits(tempfile) do file
        hdu = file[2]
        @test hdu["KEY_B1"].value() === true
        @test hdu["KEY_B2"].value() === false
        # Read sub-regions.
        hdu = file[1]
        @test read(hdu, :, :, :) == A
        let inds = (1, axes(A)[2], :)
            @test read(hdu, inds...) == A[inds...]
        end
        let inds = (:, axes(A)[2], 3)
            @test read(hdu, inds...) == A[inds...]
        end
        let inds = (:, 1:2:size(A,2), :)
            @test read(hdu, inds...) == A[inds...]
        end
        let inds = (2, 3, :)
            @test read(hdu, inds...) == A[inds...]
        end
    end
    # Write empty image.
    @test_throws Exception writefits!(tempfile, FitsHeader(), []) # TODO: MethodError expected
    writefits!(tempfile, FitsHeader(), Int32[])
    arr = readfits(Array, tempfile)
    @test length(arr) == 0
    @test eltype(arr) == Int32
    # Read FITS header.
    @test read(FitsHeader, tempfile) isa FitsHeader
end

@testset "FITS Tables" begin
    # Type conversion for reading.
    let eltype_to_read = EasyFITS.eltype_to_read
        @inferred UInt8   eltype_to_read(String)
        @inferred UInt8   eltype_to_read(String,UInt8)
        @inferred UInt8   eltype_to_read(String,String)
        @inferred UInt8   eltype_to_read(UInt8)
        @inferred Int16   eltype_to_read(UInt8,Int16)
        @inferred Float64 eltype_to_read(Float64)
        @inferred Float32 eltype_to_read(Float64,Float32)
        @inferred Float64 eltype_to_read(Float64,Rational)
        @test_throws ArgumentError eltype_to_read(String,Int8)
        @test_throws ArgumentError eltype_to_read(Float64,String)
        @test_throws ArgumentError eltype_to_read(Float64,Int)
        @test_throws ArgumentError eltype_to_read(ComplexF32,Int)
    end

    # Low-level API.
    openfits(tempfile, "w!") do file
        @test length(file) == 0
        hdu = write(file, FitsTableHDU, ["Col#1" => ('E', "m/s"),
                                         "Col#2" => ('D', "Hz")])
        @test length(file) == 2 # a table cannot not be the primary HDU
        @test hdu === last(file)
        @test hdu.file === file
        @test hdu.ncols == 2
        @test ncol(hdu) == 2
        @test hdu.first_column == 1
        @test hdu.columns == hdu.first_column:hdu.last_column
        @test length(hdu.columns) == hdu.ncols
        @test hdu.nrows == 0
        @test nrow(hdu) == 0
        @test hdu.first_row == 1
        @test hdu.rows == hdu.first_row:hdu.last_row
        @test length(hdu.rows) == hdu.nrows
        @test hdu.data_ndims == 2
        @test hdu.data_size == (hdu.nrows, hdu.ncols)
        @test hdu.data_axes == (hdu.first_row:hdu.last_row,
                                hdu.first_column:hdu.last_column)
        @test hdu.column_names == ["Col#1", "Col#2"]
        c1 = hdu.column_number("Col#1")
        @test hdu.column_name(c1) == "COL#1"
        @test hdu.column_name(c1; case=true) == "Col#1"
        c2 = hdu.column_number("Col#2")
        @test hdu.column_name(c2) == "COL#2"
        @test hdu.column_name(c2; case=true) == "Col#2"
        @test hdu.column_units("Col#1") == "m/s"
        @test hdu.column_units("Col#2") == "Hz"
        @test hdu[:tunit1].value() == "m/s"
        @test hdu[:tunit2].value() == "Hz"
        for col ∈ (1, "Col#1")
            local x = read(hdu, col)
            @test x isa Vector{Cfloat}
            @test size(x) == (0,)
        end
        for col ∈ (2, "col#2  ") # case and trailing spaces are ignored
            local x = read(hdu, col)
            @test x isa Vector{Cdouble}
            @test size(x) == (0,)
        end
        n1 = 4
        x1 = 2:2:2*n1
        y1 = 3:3:3*n1
        write(hdu, 1 => x1)
        @test hdu.nrows == nrow(hdu) == n1
        @test hdu.ncols == ncol(hdu) == 2
        let a = read(hdu, 1)
            @test a isa Vector{Cfloat}
            @test size(a) == (n1,)
            @test a == x1
        end
        write(hdu, "col#2" => y1)
        @test hdu.nrows == nrow(hdu) == n1
        @test hdu.ncols == ncol(hdu) == 2
        let a = read(hdu, 2)
            @test a isa Vector{Cdouble}
            @test size(a) == (n1,)
            @test a == y1
        end
        n2 = 3
        x2 = 5:5:5*n2
        y2 = 7:7:7*n2
        write(hdu, 1 => x2; first = hdu.nrows + 1)
        @test hdu.nrows == nrow(hdu) == n1 + n2
        @test hdu.ncols == ncol(hdu) == 2
        let a = read(hdu, 1)
            @test a isa Vector{Cfloat}
            @test size(a) == (n1 + n2,)
            @test a == vcat(x1, x2)
        end
        write(hdu, 2 => y2; first = n1 + 1)
        @test hdu.nrows == nrow(hdu) == n1 + n2
        @test hdu.ncols == ncol(hdu) == 2
        let a = read(hdu, 2)
            @test a isa Vector{Cdouble}
            @test size(a) == (n1 + n2,)
            @test a == vcat(y1, y2)
        end
        n3 = 2
        x3 = -1:-1:-1*n3
        y3 = 11:11:11*n3
        write(hdu, "Col#2" => y3, 1 => x3; first = hdu.nrows)
        @test hdu.nrows == nrow(hdu) == n1 + n2 + n3 - 1
        @test hdu.ncols == ncol(hdu) == 2
        let a = read(hdu, 1)
            @test a isa Vector{Cfloat}
            @test size(a) == (hdu.nrows,)
            @test a == vcat(x1, x2[1:end-1], x3)
        end
        let a = read(hdu, 2)
            @test a isa Vector{Cdouble}
            @test size(a) == (hdu.nrows,)
            @test a == vcat(y1, y2[1:end-1], y3)
        end
        # Read table as a dictionary.
        dict = read(hdu)
        @test dict isa Dict{String,<:Array}
        @test sort(collect(keys(dict))) == sort(map(uppercase, hdu.column_names))
        # Read table as a dictionary (rewriting names) and as a vector.
        dict = read(hdu; rename=lowercase)
        @test dict isa Dict{String,<:Array}
        @test sort(collect(keys(dict))) == sort(map(lowercase, hdu.column_names))
        vect = read(Vector, hdu)
        @test vect isa Vector{<:Array}
        @test length(vect) == length(dict)
        @test vect[1] == dict[lowercase(hdu.column_names[1])]
        @test vect[2] == dict[lowercase(hdu.column_names[2])]
        # Read table with units.
        dict = read(hdu; units=String, rename=identity)
        @test dict isa Dict{String,<:Tuple{<:Array,String}}
        vect = read(Vector, hdu; units=String)
        @test vect isa Vector{Tuple{<:Array,String}}
        @test length(vect) == length(dict)
        @test vect[1] == dict[hdu.column_names[1]]
        @test vect[2] == dict[hdu.column_names[2]]
        # Read some columns and some rows.
        cols = read(Vector, hdu) # to have all columns data
        dict = read(hdu, (hdu.column_names[2],), 2:3; rename=identity)
        @test dict isa Dict{String,<:Array}
        @test dict[hdu.column_names[2]] == cols[2][2:3]
        dict = read(hdu, reverse(hdu.column_names), 5:5; rename=identity)
        @test dict isa Dict{String,<:Array}
        @test dict[hdu.column_names[1]] == cols[1][5:5]
        @test dict[hdu.column_names[2]] == cols[2][5:5]
        dict = read(hdu, :, 5; rename=identity) # read a single row
        @test dict isa Dict{String,<:Array}
        @test dict[hdu.column_names[1]] == fill(cols[1][5])
        @test dict[hdu.column_names[2]] == fill(cols[2][5])
    end

    # High-level API.
    table = (Index = 1:4,
             Speed = (Float32.(-1:2:5), "cm/s"),
             Length = (Int16.(8:3:17), "mm"),
             Misc = reshape(Int8.(1:24), (2,3,4)),
             #= Name = ["Alice","Bob","Casimir","Duff"] =#)
    writefits!(tempfile, nothing, table)
    dict = readfits(tempfile, ext=2)
    @test Set(keys(dict)) == Set(map(uppercase∘String, keys(table)))
    @test dict["INDEX"]   == column_values(table[:Index])
    @test dict["SPEED"]   == column_values(table[:Speed])
    @test dict["LENGTH"]  == column_values(table[:Length])
    @test dict["MISC"]    == column_values(table[:Misc])
    rows = 2:3
    dict = readfits(tempfile, :, rows; ext=2, case=true)
    @test Set(keys(dict)) == Set(map(String, keys(table)))
    @test dict["Index"]   == column_values(table[:Index])[rows]
    @test dict["Speed"]   == column_values(table[:Speed])[rows]
    @test dict["Length"]  == column_values(table[:Length])[rows]
    @test dict["Misc"]    == column_values(table[:Misc])[:,:,rows]
    for col in keys(table)
        @test readfits(tempfile, col; ext=2) == column_values(table[col])
    end

    # Complex example.
    arr = rand(Float32, (3,4,5));
    nrows = 20;
    inds = 1:nrows;
    speed = rand(Float64, nrows);
    mass = rand(Float32, nrows);
    position = rand(Float32, 3, nrows);
    phase = (1:7) .// 3;
    amplitude = exp.(-1:-1:-7);
    name = ["Row#$i" for i in 1:length(phase)]; # a vector of strings
    x = amplitude.*cos.(phase);
    y = amplitude.*sin.(phase);
    xy = hcat(x,y)';
    label = Array{String,2}(undef, 2, length(x));
    for i in 1:length(x);
        label[1,i] = "x$i";
        label[2,i] = "y$i";
    end
    date = now();
    writefits!(
        tempfile,
        #-----------------------------------------------------------------
        # First HDU must be a FITS "image", but data may be empty.
        #
        # Header part as a vector of `key=>val` or `key=>(val,com)` pairs:
        ["DATE"    => (date, "date of creation"),
         "HISTORY" => "This file has been produced by EasyFITS",
         "USER"    => "John Doe"],
        # Data part as an array:
        arr,
        #-----------------------------------------------------------------
        # Second HDU, here a FITS "table".
        #
        # Header part of 2nd HDU as a tuple of pairs:
        ("EXTNAME" => ("MY-EXTENSION", "Name of this extension"),
         "EXTVER"  => (1, "Version of this extension")),
        # Data part is a table in the form of a named tuple:
        (Speed    = (speed, "km/s"),  # this column has units
         Indices  = inds,             # not this one
         Mass     = (mass, "kg"),
         Position = (position, "cm")),
        #-----------------------------------------------------------------
        # Third HDU, another FITS "table".
        #
        # Header part of 3rd HDU as a named tuple (note that keywords must
        # be in uppercase letters):
        (EXTNAME = ("MY-OTHER-EXTENSION", "Name of this other extension"),
         EXTVER  = (1, "Version of this other extension"),
         COMMENT = "This is an interesting comment"),
        # Data part is a table in the form of a vector of pairs (colum names
        # can be strings or symbols but not a mixture):
        [:phase => ((180/π).*phase, "deg"),
         :amplitude => (amplitude, "V"),
         :name => name,
         :xy => (xy, "V"),
         :label => label])
    h1 = read(FitsHeader, tempfile)
    @test h1 isa FitsHeader
    @test h1["DATE"].value(DateTime) === date
    @test h1["USER"].value() == "John Doe"
    @test readfits(tempfile) == arr
    h2 = read(FitsHeader, tempfile, ext=2)
    @test h2 isa FitsHeader
    @test h2["EXTNAME"].value() == "MY-EXTENSION"
    x2 = readfits(tempfile, ext=2)
    @test x2 isa Dict{String}
    @test x2["SPEED"] == speed
    @test x2["INDICES"] == inds
    @test x2["MASS"] == mass
    @test x2["POSITION"] == position
    h3 = read(FitsHeader, tempfile, ext="MY-OTHER-EXTENSION")
    @test h3 isa FitsHeader
    @test h3["EXTNAME"].value() == "MY-OTHER-EXTENSION"
    x3 = readfits(tempfile, ext="MY-OTHER-EXTENSION")
    @test x3 isa Dict{String}
    @test x3["PHASE"] == (180/π).*phase
    @test x3["AMPLITUDE"] == amplitude
    @test x3["NAME"] == name
    @test x3["XY"] == xy
    @test x3["LABEL"] == label

    # write (version with columns given as varargs)
    openfits(tempfile, "w!") do fitsfile
        @test write(fitsfile, FitsTableHDU, :col1 => Float32) isa FitsTableHDU
        @test write(fitsfile, FitsTableHDU, :col1 => Float32, :col2 => Int) isa FitsTableHDU
    end

    # String columns, check dimensions and write
    openfits(tempfile, "w!") do fitsfile
        # column of single character strings
        let data = ["a","b","c"],
            hdu = @inferred write(fitsfile, FitsTableHDU, [:col1 => (String, 1)])
            @test hdu isa FitsTableHDU
            @test write(hdu, :col1 => data) isa FitsTableHDU
            @test read(hdu, :col1) == data
            bytes = read(Array{UInt8}, hdu, :col1)
            @test bytes isa Array{UInt8}
            @test size(bytes) == (1,3)
        end

        # 10-char strings
        let data = ["abcdefghij", "abcd", "abcdefghi"],
            hdu = @inferred write(fitsfile, FitsTableHDU, :col1 => (String, 10))
            @test hdu isa FitsTableHDU
            @test write(hdu, :col1 => data) isa FitsTableHDU
            @test read(hdu, :col1) == data
        end

        # pairs of 4-char strings
        let data = ["abcd" ; "defg" ;; "hijk" ; "lmno" ;; "pq" ; "rstu"],
            hdu = @inferred write(fitsfile, FitsTableHDU, :col1 => (String, (4,2)))
            @test hdu isa FitsTableHDU
            @test write(hdu, :col1 => data) isa FitsTableHDU
            @test read(hdu, :col1) == data
        end

        # dimension is mandatory
        @test_throws ArgumentError write(fitsfile, FitsTableHDU, :col1 => String)

        # overlarge String is a failure
        let hdu = @inferred write(fitsfile, FitsTableHDU, :col1 => (String, 3))
            @test_throws ArgumentError write(hdu, :col1 => ["abcd"])
        end
    end
end

end # module

nothing
