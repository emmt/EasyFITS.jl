srcdir = .
oldcode = $(srcdir)/../deps.bak/deps.jl

default: cfitsio.jl

cfitsio.jl: build.jl generate_wrappers.jl
	julia --project -e 'using Pkg; Pkg.instantiate();'
	julia --project build.jl

ccalls-old.jl: $(oldcode)
	grep '^ *@ccall' "$<" \
	    | sed -e 's/^  *//;s/ULONGLONG/UInt64/g;s/LONGLONG/Int64/g;s/Cuchar/UInt8/g;' \
	    | sort > "$@"
ccalls-new.jl: $(srcdir)/cfitsio.jl
	grep '^ *@ccall' "$<" \
	    | sed -e 's/^  *//;s/#= *const *=#//g;' \
	    | sort > "$@"
