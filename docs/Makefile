SOURCES = \
    make.jl \
    src/files.md \
    src/hdus.md \
    src/headers.md \
    src/images.md \
    src/tables.md \
    src/install.md \
    src/grammar.md \
    src/library.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

force:
	julia --color=yes make.jl

build/index.html: $(SOURCES)
	make force

.PHOONY: default clean force
