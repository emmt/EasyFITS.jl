module TestingEasyFITS

using Test, EasyFITS, DataFrames

using EasyFITS: FitsInteger, FitsFloat, FitsComplex

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
    @test EasyFITS.library_version() isa VersionNumber
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
        @test_throws ArgumentError new_array(Float32, Val(2),  4, 5, 6)
        @test_throws ArgumentError new_array(Int16,   Val(2), (4, 5, 6))
        @test_throws ArgumentError new_array(UInt32,  Val(2), [4, 5, 6])
    end
    let dense_array = EasyFITS.dense_array,
        A = convert(Array, reshape(1:24, 2,3,4)),
        B = view(A, :, 1:2:3, :)
        @test dense_array(A) === A
        @test isa(dense_array(B), Array)
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

cards_1 = (key_b1 = (true,  "This is true"),
           key_b2 = (false, "This is false"),
           key_b3 = false,
           key_b4 = (true, ""),
           key_i1 = 123, # no comment
           key_i2 = (0x9, nothing),
           key_i3 = (Int16(42), "[cm] with units"),
           key_i4 = (-77, "[] no units"),
           key_i5 = (101, "[  foo /  bar ] strip leading/trailing spaces"),
           key_f1 = (-1.2f-2, "Single precision"),
           key_f2 = (-3.7, "Double precision"),
           key_f3 = (11//7, "Rational number"),
           key_f4 = (pi, "Irrational number"),
           key_c1 = (complex(-2,3), "Integer complex"),
           key_c2 = (complex(2.1f0,-3.7f0), "Single precision complex"),
           key_c3 = (complex(2.1,-3.7), "Double precision complex"),
           key_s1 = " <- significant space here",
           key_s2 = ("", "Empty string"),
           key_s3 = (" ", "Another empty string"),
           key_s4 = ("'oops!'", "String with quotes"),
           key_s5 = " <- keep this, not these ->    ",
           key_s6 = (SubString("Hello world!", 7, 12), "A sub-string"),
           key_n2 = (undef, "undefined value"),
           key_n3 = (missing, "missing value"),
           comment = "simple comment",
           history = "historical comment",
           very_long_key_name = (1, "Should use HIERARCH convention"),
           uncommented_very_long_key_name = 2.0,
           )

cards_2 = [uppercase(String(key)) => val for (key,val) in pairs(cards_1)]
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
    # Write empty image
    @test_throws Exception writefits!(tempfile, FitsHeader(), []) #TODO: MethodError expected
    writefits!(tempfile, FitsHeader(), Int32[])
    arr = readfits(Array, tempfile)
    @test length(arr) == 0
    @test eltype(arr) == Int32
    # Read FITS header
    @test read(FitsHeader, tempfile) isa FitsHeader
end

@testset "FITS Tables" begin
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
        @test hdu.last_column == hdu.ncols + 1 - hdu.first_column
        @test hdu.nrows == 0
        @test nrow(hdu) == 0
        @test hdu.first_row == 1
        @test hdu.last_row == hdu.nrows + 1 - hdu.first_row
        @test hdu.data_ndims == 2
        @test hdu.data_size == (hdu.nrows, hdu.ncols)
        @test hdu.data_axes == (hdu.first_row:hdu.last_row,
                                hdu.first_column:hdu.last_column)
        @test hdu.column_names == ["Col#1", "Col#2"]
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
    end
end

end # module

nothing
