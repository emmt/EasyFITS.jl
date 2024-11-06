SOURCES = \
    make.jl \
    src/index.md \
    src/structure.md \
    src/reading.md \
    src/writing.md \
    src/files.md \
    src/images.md \
    src/tables.md \
    src/library.md \
    src/links.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

force:
	julia --color=yes make.jl

build/index.html: $(SOURCES)
	make force

.PHOONY: default clean force
