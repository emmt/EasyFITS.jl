module TestingEasyFITS

using Test, EasyFITS, DataFrames

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

    # new_array
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
    m = match(r"^\[(.*?)\s*\]", comment)
    m === nothing && return ""
    return m.captures[1]
end

suppress_units(comment::AbstractString) =
    replace(comment, r"^\[.*?\]\s*" => "")

cards_1 = (key_b1 = (true,  "This is true"),
           key_b2 = (false, "This is false"),
           key_b3 = false,
           key_b4 = (true, ""),
           key_i1 = 123, # no comment
           key_i2 = (0x9, nothing),
           key_i3 = (42, "[cm] with units"),
           key_i4 = (-77, "[] no units"),
           key_i5 = (101, "[  foo / bar ] trailing spaces do not matter"),
           key_f1 = (-1.2f0, "Single precision"),
           key_f2 = (-3.7, "Double precision"),
           key_f3 = (11//7, "Rational number"),
           key_f4 = (pi, "Irrational number"),
           key_c1 = (complex(2,3), "Integer complex"),
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
           )

cards_2 = [uppercase(String(key)) => val for (key,val) in pairs(cards_1)]

@testset "FITS Images" begin
    # Write a simple FITS image.
    A = convert(Array{Int16}, reshape(1:60, 3,4,5))
    @test fits"test1.fits" === FitsFile("test1.fits")
    open(fits"test1.fits", "w!") do io
        @test io isa FitsIO
        @test !isreadable(io)
        @test !isreadonly(io)
        @test iswritable(io)
        @test position(io) == 1
        @test length(io) == 0

        # Add a simple IMAGE extension.
        let hdu = write(io, FitsImageHDU, eltype(A), size(A))
            @test length(io) == 1
            @test hdu === last(io)
            @test firstindex(hdu) == 1
            @test lastindex(hdu) == length(hdu)
            @test_throws KeyError hdu[firstindex(hdu) - 1]
            # FIXME: @test_throws KeyError hdu[lastindex(hdu) + 1]
            @test_throws KeyError hdu[lastindex(hdu) + 5000]
            @test_throws KeyError hdu["DUMMY"]
            let card = hdu["SIMPLE"]
                @test card == hdu[firstindex(hdu)]
                @test card.type == FITS_LOGICAL
                @test card.value.logical === card.value.parsed
                @test card.value.logical === true
            end
            let card = hdu["BITPIX"]
                @test card == hdu[firstindex(hdu) + 1]
                @test card.type == FITS_INTEGER
                @test card.value.integer === card.value.parsed
                @test card.value.integer == 16
            end
            let card = hdu["NAXIS"]
                @test card == hdu[firstindex(hdu) + 2]
                @test card.type == FITS_INTEGER
                @test card.value.integer === card.value.parsed
                @test card.value.integer == 3
            end
            reset(hdu) # reset before testing incremental search
            for i in 1:ndims(A)+1
                local card = get(hdu, "NAXIS#", nothing)
                @test (card === nothing) == (i > ndims(A))
                card === nothing && break
                @test card.type == FITS_INTEGER
                @test card.value.integer === card.value.parsed
                @test card.value.integer == size(A, i)
            end
            write(hdu, A)
        end
        @test position(io) == 1
        @test length(io) == 1
        # Add IMAGE extensions with no data, just headers.
        for pass in 1:2
            local cards = pass == 1 ? cards_1 : cards_2
            local hdu = write(io, FitsImageHDU)
            len = length(hdu)
            push!(hdu, cards)
            @test length(hdu) == len + length(cards)
            for (key, dat) in (cards isa AbstractVector ? cards : pairs(cards))
                local card = get(hdu, key, nothing)
                @test card isa EasyFITS.FitsCard
                if uppercase(rstrip(String(key))) ∈ ("HISTORY","COMMENT","")
                    local com = dat isa Tuple{Nothing,AbstractString} ? dat[2] :
                        dat isa AbstractString ? dat :
                        dat isa Nothing || dat isa Tuple{Nothing,Nothing} ? "" :
                        error("bad card specification: $(repr(key)) => $(repr(dat))")
                    @test card.type == FITS_COMMENT
                    @test card.value.parsed === nothing
                    @test card.comment == com
                else
                    local val = dat isa Tuple ? dat[1] : dat
                    local com = dat isa Tuple{Any,AbstractString} ? dat[2] : ""
                    if val isa Bool
                        @test card.type == FITS_LOGICAL
                        @test card.value.logical === card.value.parsed
                        @test card.value.logical === val
                    elseif val isa Integer
                        @test card.type == FITS_INTEGER
                        @test card.value.integer === card.value.parsed
                        @test card.value.integer == val
                    elseif val isa Real
                        @test card.type == FITS_FLOAT
                        @test card.value.float === card.value.parsed
                        @test card.value.float ≈ val
                    elseif val isa Complex
                        @test card.type == FITS_COMPLEX
                        @test card.value.complex === card.value.parsed
                        @test card.value.complex ≈ val
                    elseif val isa Complex
                        @test card.type == FITS_COMPLEX
                        @test card.value.complex === card.value.parsed
                        @test card.value.complex ≈ val
                    elseif val isa AbstractString
                        @test card.type == FITS_STRING
                        @test card.value.string === card.value.parsed
                        @test card.value.string == rstrip(card.value.string)
                        @test card.value.string == rstrip(val)
                    elseif val isa Union{Missing,UndefInitializer}
                        @test card.type == FITS_UNDEFINED
                        @test card.value.parsed === missing
                    else
                        error("bad card specification: $(repr(key)) => $(repr(dat))")
                    end
                    @test card.comment == com
                    @test card.comment.units == get_units(com)
                    @test card.comment.unitless == suppress_units(com)
                end
            end
        end
    end
    # Read the data.
    B = readfits("test1.fits")
    @test eltype(B) == eltype(A)
    @test size(B) == size(A)
    @test B == A
    openfits("test1.fits") do io
        hdu = io[2]
        @test hdu["KEY_B1"].value.parsed === true
        @test hdu["KEY_B2"].value.parsed === false
        # Read sub-regions.
        hdu = io[1]
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
end

@testset "FITS Tables" begin
    openfits("test2.fits", "w!") do io
        @test length(io) == 0
        hdu = write(io, FitsTableHDU, ["Col#1" => ('E', "m/s"),
                                       "Col#2" => ('D', "Hz")])
        @test length(io) == 2 # a table cannot not be the primary HDU
        @test hdu === last(io)
        @test hdu.ncols == 2
        @test hdu.nrows == 0
        @test ncol(hdu) == 2
        @test nrow(hdu) == 0
        @test hdu.column_names == ["Col#1", "Col#2"]
        @test hdu[:tunit1].value.parsed == "m/s"
        @test hdu[:tunit2].value.parsed == "Hz"
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
        write(hdu, 1 => x1)
        @test hdu.nrows == nrow(hdu) == n1
        @test hdu.ncols == ncol(hdu) == 2
        let y = read(hdu, 1)
            @test y isa Vector{Cfloat}
            @test size(y) == (n1,)
            @test y == x1
        end
        x2 = 3:3:3*n1
        write(hdu, "col#2" => x2)
        @test hdu.nrows == nrow(hdu) == n1
        @test hdu.ncols == ncol(hdu) == 2
        let y = read(hdu, 2)
            @test y isa Vector{Cdouble}
            @test size(y) == (n1,)
            @test y == x2
        end
        n2 = 3
        y1 = 5:5:5*n2
        write(hdu, 1 => y1; first = n1 + 1)
        @test hdu.nrows == nrow(hdu) == n1 + n2
        @test hdu.ncols == ncol(hdu) == 2
        let y = read(hdu, 1)
            @test y isa Vector{Cfloat}
            @test size(y) == (n1 + n2,)
            @test y == vcat(x1, y1)
        end
        y2 = 7:7:7*n2
        write(hdu, 2 => y2; first = n1 + 1)
        @test hdu.nrows == nrow(hdu) == n1 + n2
        @test hdu.ncols == ncol(hdu) == 2
        let y = read(hdu, 2)
            @test y isa Vector{Cdouble}
            @test size(y) == (n1 + n2,)
            @test y == vcat(x2, y2)
        end
    end
end

end # module

nothing
