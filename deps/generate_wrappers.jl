module GenerateWrappers

using CFITSIO_jll

# The dynamic library.
const libcfitsio = CFITSIO_jll.libcfitsio

# Where are the C header files?
const default_incdir = joinpath(CFITSIO_jll.artifact_dir, "include")

# Names to skip.
const skip_names = Set(["_FITSIO_H", "CFITS_API",
                       "mipsFortran", "_MIPS_SZLONG",
                       "OFF_T", "USE_LL_SUFFIX",
                       # `LONGLONG` is `Int64`, usually `Clonglong`
                       "LONGLONG", "LONGLONG_TYPE", "LONGLONG_MIN", "LONGLONG_MAX",
                       # `INT32BIT` is `Int32`, usually `Cint`
                       "INT32BIT",
                       # Iterators.
                       "InputCol", "InputOutputCol", "OutputCol", "TemporaryCol"])

# Table of type correspondances. It should be valid for different architectures.
const TYPES = Dict(
    # `fitsfile` is an opaque structure (only its address is needed).
    "fitsfile" => "fitsfile",
    # In CFITSIO, `INT32BIT` is a 32-bit integer.
    "INT32BIT" => "Int32",
    "unsigned INT32BIT" => "UInt32",
    # In CFITSIO, `LONGLONG` is a 64-bit integer.
    "LONGLONG" => "Int64",
    "ULONGLONG" => "UInt64",
    "unsigned LONGLONG" => "UInt64",
    # `char` is converted to `Cchar` (for strings) and `Ptr{Cchar}` may be replaced by
    # `Cstring`.
    "char" => "Cchar",
    # In CFITSIO, `signed char` and `unsigned char` are respectively used for signed and
    # unsigned 8-bit integers.
    "signed char" => "Int8",
    "unsigned char" => "UInt8",
    # Other C types have a suitable alias in Julia.
    "short" => "Cshort",
    "unsigned short" => "Cushort",
    "int" => "Cint",
    "unsigned int" => "Cuint",
    "long" => "Clong",
    "unsigned long" => "Culong",
    "double" => "Cdouble",
    "float" => "Cfloat",
    "size_t" => "Csize_t",
    "void" => "Cvoid")

# List of reserved keywords in Julia. Note that `abstract`, `mutable`, or `type` are not
# in this list as they must be combined with another keyword.
const RESERVED = Set([
    "baremodule",
    "begin",
    "break",
    "catch",
    "const",
    "continue",
    "do",
    "else",
    "elseif",
    "end",
    "export",
    "false",
    "finally",
    "for",
    "function",
    "global",
    "if",
    "import",
    "let",
    "local",
    "macro",
    "module",
    "quote",
    "return",
    "struct",
    "true",
    "try",
    "using",
    "while",
])

const REQUIRED_FUNCTIONS = [
    "fits_close_file",
    "fits_create_diskfile",
    "fits_create_file",
    "fits_create_img",
    "fits_create_tbl",
    "fits_delete_key",
    "fits_delete_record",
    "fits_flush_buffer",
    "fits_get_colname",
    "fits_get_colnum",
    "fits_get_errstatus",
    "fits_get_hdrspace",
    "fits_get_hdu_num",
    "fits_get_img_dim",
    "fits_get_img_equivtype",
    "fits_get_img_size",
    "fits_get_img_type",
    "fits_get_num_cols",
    "fits_get_num_hdus",
    "fits_get_num_rowsll",
    "fits_get_version",
    "fits_modify_comment",
    "fits_movabs_hdu",
    "fits_open_diskfile",
    "fits_open_file",
    "fits_read_card",
    "fits_read_col_<type_suffix>",
    "fits_read_colnull_<type_suffix>",
    "fits_read_img",
    "fits_read_imgnull",
    "fits_read_record",
    "fits_read_subset",
    "fits_read_tdimll",
    "fits_update_card",
    "fits_write_col_<type_suffix>",
    "fits_write_colnull_<type_suffix>",
    "fits_write_comment",
    "fits_write_date",
    "fits_write_history",
    "fits_write_img",
    "fits_write_imgnull",
    "fits_write_record",
    "fits_write_subset",
    "fits_write_tdim",
    "fits_write_tdimll"]

const COLUMN_TYPES = Dict(
    "str" => String,
    "log" => Bool,
    "byt" => UInt8,
    "sbyt" => Int8,
    "usht" => Cushort,
    "sht" => Cshort,
    "ulng" => Culong,
    "lng" => Clong,
    "ulnglng" => Culonglong,
    "lnglng" => Clonglong,
    "uint" => Cuint,
    "int" => Cint,
    "flt" => Cfloat,
    "dbl" => Cdouble,
    "cmp" => Complex{Cfloat},
    "dblcmp" => Complex{Cdouble})

function check_wrappers(filename::AbstractString = joinpath(@__DIR__, "cfitsio.jl"))
    # List of functions.
    funcs = Set{String}()
    for line in readlines(filename)
        m = match(r"^ *function +(\w+) *\(", line)
        if m !== nothing
            push!(funcs, m.captures[1])
            continue
        end
        m = match(r"^ *(\w+) *\(.*\) *=", line)
        if m !== nothing
            push!(funcs, m.captures[1])
            continue
        end
    end

    for func in REQUIRED_FUNCTIONS
        if endswith(func, "<type_suffix>")
            base = SubString(func, firstindex(func) : prevind(func, findfirst('<', func)))
            for sfx in keys(COLUMN_TYPES)
                base*sfx ∈ funcs || println(stderr, "function `$base$sfx` is missing")
            end
        else
            func ∈ funcs || println(stderr, "function `$func` is missing")
        end
    end
    nothing
end

function extract_comment(str::AbstractString)
    # Remove leading and trailing spaces.
    i_first = firstindex(str)
    i_last = lastindex(str)
    while i_first ≤ i_last && isspace(str[i_first])
        i_first = nextind(str, i_first)
    end
    while i_first < i_last && isspace(str[i_last])
        i_last = prevind(str, i_last)
    end
    if i_first > i_last
        return ""
    elseif str[i_first] == '/'
        i = nextind(str, i_first)
        if i ≤ i_last
            c = str[i]
            if c == '/' || c == '*'
                # Skip spaces after begin-of-comment marker.
                while true
                    i = nextind(str, i)
                    (i ≤ i_last && isspace(str[i])) || break
                end
                c == '/' && return String(SubString(str, i:i_last))
                # Find end-of-comment marker.
                j = prevind(str, i_last)
                if j ≥ i && str[j] == '*' && str[i_last] == '/'
                    # Skip spaces before end-of-comment marker.
                    while true
                        j = prevind(str, j)
                        (j ≥ i && isspace(str[j])) || break
                    end
                    # Return comment with ordinary spaces and with borders and leading or
                    # trailing spaces removed.
                    return String(replace(SubString(str, i:j),
                                          r"\t" => " ",
                                          r"^[ *]+"m => "",
                                          r"[ *]+$"m => ""))
                end
            end
        end
    end
    error("invalid comment: \"$(str)\"")
end

# FIXME: Return type may be more complex than just `int`.
function convert_cfunction(;
                           io::IO = stdout,
                           rtype::AbstractString,
                           func::AbstractString,
                           alias::AbstractString = func,
                           args::AbstractString,
                           char_array_is_not_string::Bool = false,
                           quiet::Bool = false)
    # Convert return type.
    rtype, success = if haskey(TYPES, rtype)
        TYPES[rtype], true
    else
        success = false
        quiet || @warn "unknown return C type `$(rtype)` for function `$func`"
        "<unknown>", false
    end

    # Remove block comments. FIXME: What about // comments?
    args = replace(args, r"/\*.*?\*/"s => " ")

    # Replace any sequence of spaces (including newlines, etc.) by a single space.
    args = replace(args, r"\s+" => " ")

    # Replace pointers to function by simple `void` pointers in argument list.
    args = replace(args, r"(^|,)[a-zA-Z_0-9* ]*\([ *]*(\w+) *\) *\([^)]*\) *(,|$)" => s"\1void*\2\3")

    # Parse argument list.
    arglist = Tuple{String,String,Bool}[]
    if args != "void"
        # Split arguments and convert them.
        words = SubString{String}[]
        for (argnum, arg) in enumerate(split(args, ","))
            # Count number of pointer levels.
            nptrs = 0
            for c in arg
                if (c == '*')|(c == '[')
                    nptrs += 1
                end
            end

            # Set stop index to strip trailing array dimensions `[...]` if any.
            i = findfirst('[', arg)
            stop = i === nothing ? lastindex(arg) : prevind(arg, i)

            # Parse argument definition: type and name.
            is_const = false
            empty!(words)
            i = firstindex(arg)
            while i ≤ stop
                m = match(r"\w+", arg, i)
                m === nothing && break
                i = m.offset + ncodeunits(m.match)
                if m.match == "const"
                    is_const = true
                else
                    push!(words, m.match)
                end
            end
            argname = if length(words) > 1
                name = pop!(words)
                name ∈ RESERVED ? convert(eltype(words), "_$(name)") : name
            else
                convert(eltype(words), "_arg_$(argnum)")
            end
            if length(words) < 1
                quiet || @warn "missing type argument n° $(argnum) of function `$func`"
                success = false
                argtype = "<unknown>"
            else
                argtype = String(join(words, " "))
                if !haskey(TYPES, argtype)
                    quiet || @warn "unknown argument n° $(argnum) type `$(argtype)` for function `$func`"
                    success = false
                else
                    argtype = TYPES[argtype]
                    if nptrs > 0 && !char_array_is_not_string && argtype == "Cchar"
                        argtype = oftype(argtype, "Cstring")
                        nptrs -= 1
                    end
                    while nptrs > 0
                        argtype = "Ptr{"*argtype*"}"
                        nptrs -= 1
                    end
                end
            end
            push!(arglist, (argname, argtype, is_const))
        end
    end
    if success
        print(io, alias, "(")
        cnt = 0
        for (argname, argtype, is_const) in arglist
            (cnt += 1) > 1 && print(io, ", ")
            print(io, argname)
        end
        print(io, ") =\n")
        print(io, "    @ccall libcfitsio.", func, "(")
        cnt = 0
        for (argname, argtype, is_const) in arglist
            (cnt += 1) > 1 && print(io, ", ")
            print(io, argname, "::")
            is_const && print(io, "#= const =#")
            print(io, argtype)
        end
        print(io, ")::", rtype, "\n")
    end
    return success
end

function generate_wrappers(filename::AbstractString = joinpath(@__DIR__, "cfitsio.jl");
                           kwds...)
    open(filename, "w") do io
         generate_wrappers(io; kwds...)
    end
    nothing
end

function generate_wrappers(io::IO; incdir::AbstractString = default_incdir, kwds...)
    # Parse <longnam.h> and to memorize `#define long_name short_name` aliases.
    aliases = Dict{String,String}()
    open(joinpath(incdir, "longnam.h"), "r") do file
        while !eof(file)
            # Read next line.
            str = readline(file)

            # Simple alias.
            m = match(r"^\s*#\s*define\s+(fits_[a-z_0-9]*)\s+([a-z][a-z_0-9]*)\s*$", str)
            if m !== nothing
                aliases[m.captures[2]] = m.captures[1]
            end

            # Function-like alias.
            m = match(r"^\s*#\s*define\s+(fits_[a-z_0-9]*) *\( *(.*?) *\) +(\w+) *\( *(.*?) *\)", str)
            if m !== nothing
                aliases[m.captures[3]*"("*m.captures[4]*")"] =
                    m.captures[1]*"("*m.captures[2]*")"
            end
        end
    end

    # Preamble.
    print(io,
          "module CFITSIO\n",
          "\nusing CFITSIO_jll\n",
          "\n# Dynamic library.\n",
          "const libcfitsio = CFITSIO_jll.libcfitsio\n",
          "\n# Opaque structure.\n",
          "abstract type fitsfile end\n")

    # Parse <fitsio.h> using regular expressions.
    buf = read(joinpath(incdir, "fitsio.h"), String)

    # Replace CR-LF and CR by LF, tab by spaces, and strip trailing spaces.
    buf = replace(buf, r"\r\n?" => "\n", "\t" => " ", r" +$"m => "")

    # Convert macro definitions in <fitsio.h> into constants.
    print(io, "\n# Constants.\n")
    i = firstindex(buf)
    while true
        m = match(r"^ *# *define +(\w+) +([^\n]+?) *(|//[^\n]*|/\*.*?\*/)$"ms, buf, i)
        m === nothing && break
        name = m.captures[1]
        if !(name ∈ skip_names)
            print(io, "const ", name, " = ")
            print(io, replace(m.captures[2],
                              # Remove parentheses around expressions.
                              r"^\s*\(\s*(.+?)\s*\)\s*$" => s"\1",
                              # Single precision floating-point constants.
                              r"^\s*([-+]?[0-9]\.[0-9]+)[eE]([-+]?[0-9]+)[fF]\s*$" => s"\1f\2",
                              # Double precision floating-point constants.
                              r"^\s*([-+]?[0-9]\.[0-9]+)[eE]([-+]?[0-9]+)\s*$" => s"\1e\2",
                              # Versions.
                              r"^\s*([0-9]+\.[0-9]+\.[0-9])\s*$" => s"v\"\1\""))
            comm = extract_comment(m.captures[3])
            isempty(comm) || print(io, " # ", replace(comm, "\n" => "\n# "))
            print(io, "\n")
        end
        i = m.offset + sizeof(m.match)
    end

    print(io, "\n# Methods.\n")

    # Function-like macros in <longnam.h>.
    for key in keys(aliases)
        findfirst('(', key) === nothing && continue
        print(io, aliases[key], " =\n    ", key, "\n\n")
    end

    # Convert functions definitions in <fitsio.h> into `@ccall`.
    i = firstindex(buf)
    while true
        m = match(r"^ *(\w+)\s+CFITS_API\s+(\w+)\s*\((.*?)\)\s*;\s*"ms, buf, i)
        m === nothing && break
        #println("func def: ", m.captures[2], " / ", get(aliases, m.captures[2], "<none>"))
        rtype = m.captures[1]
        func = m.captures[2]
        args = m.captures[3]
        if convert_cfunction(; io, rtype, func, args, alias = get(aliases, func, func), kwds...)
            print(io, "\n")
        end
        i = m.offset + sizeof(m.match)
    end

    print(io, "end # module CFITSIO\n")
    nothing
end

end # module

if isinteractive()
    println("You may call:\n")
    println("    GenerateWrappers.generate_wrappers(; quiet=true/false)\n")
    println("to (re-)generate file `cfitsio.jl`")
end

nothing
