using Libdl
using Clang.Generators
using Clang.LibClang.Clang_jll
using CFITSIO_jll

import Clang

# NOTE There is a bug in Clang < 0.18.0 that may prevent building EasyFITS.
# See https://github.com/JuliaInterop/Clang.jl/issues/450 for more informations.
#
# The following is a hack to patch the code.
#
# Other possibility to replace the rather crude test (but requires Compat):
#
#     import Compat
#     if Compat.pkgversion(Clang) < v"0.18.0"
#         ...
#     end
#
if isdefined(Clang, :JLLEnvs) && isdefined(Clang.JLLEnvs, :get_system_shard_key)
    @warn "A patched version of `Clang.JLLEnvs.get_system_shard_key` will be used for building `EasyFITS` (see https://github.com/JuliaInterop/Clang.jl/issues/450)"
    function Clang.JLLEnvs.get_system_shard_key(triple::String)
        @assert triple ∈ Clang.JLLEnvs.JLL_ENV_ARTIFACT_TRIPLES "Wrong JLL target triple: $triple. Please choose a triple listed in `JLL_ENV_ARTIFACT_TRIPLES`."
        platform_keys = filter(collect(keys(Clang.JLLEnvs.JLL_ENV_SHARDS))) do key
            startswith(key, "$(Clang.JLLEnvs.JLL_ENV_SYSTEM_SHARD_NAME)-$triple.") &&
            endswith(key, "$(Clang.JLLEnvs.JLL_ENV_HOST_TRIPLE).unpacked")
        end
        return first(platform_keys)
    end
end

# We wrap everything into a function to avoid having undefined variables...
function build_deps()

    prefix = CFITSIO_jll.artifact_dir
    incdir = joinpath(prefix, "include")

    # Header files.
    headers = map(x -> joinpath(incdir, x), [
        #"drvrsmem.h",
        "fitsio2.h",
        "fitsio.h",
        #"longnam.h",
    ])

    # List of (unique and in order) include directories.
    include_dirs = String[]
    for dir in Iterators.map(dirname, headers)
        dir in include_dirs || push!(include_dirs, dir)
    end

    # The rest is pretty standard.
    cd(@__DIR__)
    options = load_options(joinpath(@__DIR__, "generator.toml"))
    args = get_default_args()
    for dir in include_dirs
        push!(args, "-I$dir")
    end
    ctx = create_context(headers, args, options)
    build!(ctx)

    # Rewrite destination file.
    dest_file = options["general"]["output_file_path"]
    code = readlines(dest_file)
    for repl in [
        r"^\s*const\s+UINT64_MAX\s*=.*$" => "const UINT64_MAX = typemax(UInt64)",
        r"^\s*const\s+UINT32_MAX\s*=.*$" => "const UINT32_MAX = typemax(UInt32)",
        r"^\s*const\s+LONGLONG_MAX\s*=.*$" => "const LONGLONG_MAX = typemax(Clonglong)",
        r"^\s*const\s+LONGLONG_MIN\s*=.*$" => "const LONGLONG_MIN = typemin(Clonglong)",
        r"^\s*const\s+FLOATNULLVALUE\s*=\s*(.*)e(.*)F(.*)" => s"const FLOATNULLVALUE = \1f\2\3",
        ]
        for i in eachindex(code)
            code[i] = replace(code[i], repl)
        end
    end
    open(dest_file, "w") do io
        foreach(line -> println(io, line), code)
    end
end

# Run the build script.
build_deps()
