using Documenter

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using EasyFITS
import EasyFITS:
    Rows, Columns, ColumnIdent, ColumnName, ColumnData,
    TableData, ImageData, Header, SubArrayIndex
using FITSHeaders

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "Easy access to FITS files in Julia",
    format = Documenter.HTML(
        edit_link = "main",
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = [
        "index.md",
        "structure.md",
        "reading.md",
        "writing.md",
        "files.md",
        "images.md",
        "tables.md",
        "reference.md",
        "links.md",
    ]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/EasyFITS.jl.git",
    )
end
