module TestingEasyFITS

using Test, EasyFITS

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
@testset "Arrays utilities" begin
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
@testset "FITS files" begin
    @test fits"test1.fits" === FitsFile("test1.fits")
    open(fits"test1.fits", "w!") do io
        @test io isa FitsIO
        @test !isreadable(io)
        @test !isreadonly(io)
        @test iswritable(io)
        @test position(io) == 1
        @test length(io) == 0
        let A = convert(Array{Int16}, reshape(1:60, 3,4,5)),
            hdu = write(io, FitsImageHDU, eltype(A), size(A))
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
            hdu["GIZMO"] = ("O'Brian", "Gizmo name")
            let card = hdu["GIZMO"]
                @test card.type == FITS_STRING
                @test card.value.string === card.value.parsed
                @test card.value.string == "O'Brian"
                @test card.comment == "Gizmo name"
                @test card.comment.unitless == "Gizmo name"
                @test card.comment.units == ""
            end
            hdu["COUNT"] = (Int16(42), "Gizmo counts")
            let card = hdu["COUNT"]
                @test card.type == FITS_INTEGER
                @test card.value.integer === card.value.parsed
                @test card.value.integer == 42
                @test card.comment == "Gizmo counts"
                @test card.comment.unitless == "Gizmo counts"
                @test card.comment.units == ""
            end
            hdu["SPEED"] = (π, "[m/s] Gizmo speed")
            let card = hdu["SPEED"]
                @test card.type == FITS_FLOAT
                @test card.value.float === card.value.parsed
                @test card.value.float ≈ π
                @test card.comment == "[m/s] Gizmo speed"
                @test card.comment.unitless == "Gizmo speed"
                @test card.comment.units == "m/s"
            end
            hdu["LOGIC"] = false
            let card = hdu["LOGIC"]
                @test card.type == FITS_LOGICAL
                @test card.value.logical === card.value.parsed
                @test card.value.logical == false
                @test card.comment == ""
                @test card.comment.unitless == ""
                @test card.comment.units == ""
            end
            hdu["CMPLX"] = (2 + 1im, "[] Gizmo complex")
            let card = hdu["CMPLX"]
                @test card.type == FITS_COMPLEX
                @test card.value.complex === card.value.parsed
                @test card.value.complex == 2 + 1im
                @test card.comment == "[] Gizmo complex"
                @test card.comment.unitless == "Gizmo complex"
                @test card.comment.units == ""
            end
            write(hdu, A)
        end
        @test position(io) == 1
        @test length(io) == 1
    end
end

end # module

nothing
