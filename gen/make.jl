using Clang.Generators
import CFITSIO_jll: libcfitsio

cd(@__DIR__)

artifact_dir = normpath(dirname(dirname(libcfitsio)))
include_dir = joinpath(artifact_dir, "include")

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()
push!(args, "-I$include_dir")

headers = [joinpath(include_dir, header) for header in ("fitsio.h",)]
# there is also an experimental `detect_headers` function for auto-detecting top-level headers in the directory
# headers = detect_headers(clang_dir, args)

# create context
ctx = create_context(headers, args, options)

# run generator
build!(ctx)

# function to rewrite code
function rewrite(code::AbstractString)
    code = String(code)

    m = match(r"^ *const +LONGLONG *= *(\w+) *$"m, code)
    if m !== nothing
        T = m.captures[1]
        code = replace(code, r"\bLONGLONG\b" => T)
        code = replace(code, r"\bLLONG_MAX\b" => "typemax($T)")
        code = replace(code, r"\bLLONG_MIN\b" => "typemin($T)")
    end

    m = match(r"^ *const +ULONGLONG *= *(\w+) *$"m, code)
    if m !== nothing
        code = replace(code, r"\bULONGLONG\b" => m.captures[1])
    end

    dict = Dict{String,String}()
    index = firstindex(code)
    stop = lastindex(code)
    while index ≤ stop
        m = match(r"^ *const +(fits_\w+) *= *(ff\w+) *$()"m, code, index)
        m === nothing && break
        dict[m.captures[2]] = m.captures[1]
        index = last(m.offsets)
    end
    for (short_name, long_name) in dict
        code = replace(code, Regex("^ *function +$(short_name) *\\(", "m") => "function $(long_name)(")
        code = replace(code, Regex("^ *const +$(long_name) *= *$(short_name) *\$", "m") => "")
    end

    for rule in (
        # replace some system types like __pid_t by pid_t
        r"\b__(pid|off)_t\b" => s"\1_t",

        # after rewriting some symbols, remove aliases like "const $(ident) = $(ident)"
        r"^ *const +(\w+) *= *\1 *$"m => "",

        # comment unused aliases
        r"^ *(const +(ffcpimg|fits_(compress_img|decompress_img|write_nulrows)) *=.*?) *$"m => s"# \1",

        # single precision constants
        r"([0-9]+\.[0-9]+)[eE]([-+]?[0-9]+)F\b" => s"\1f\2",

        # suppress empty lines
        r"(\n?\r|\n)\1\1+"m => "\n\n")

        code = replace(code, rule)
    end
    return code
end

# load the generated code to rewrite part of it
output_file_path = options["general"]["output_file_path"]
code = read(output_file_path, String)

# save rewritten generated code
write(output_file_path, rewrite(code))
