using Documenter

push!(LOAD_PATH, "../src/")

using EasyFITS
import EasyFITS:
    Rows, Columns, ColumnIdent, ColumnName, ColumnData, ColumnDataPair,
    TableData, ImageData, Header, SubArrayIndex

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "Easy access to FITS files in Julia",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = ["index.md", "files.md", "hdus.md", "headers.md", "images.md",
             "tables.md", "install.md", "grammar.md", "library.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/EasyFITS.jl.git",
    )
end
