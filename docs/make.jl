using Pkg
using Documenter, XAM

makedocs(
    format = Documenter.HTML(
        edit_link = "develop"
    ),
    modules = [XAM, XAM.SAM, XAM.BAM],
    sitename = "XAM.jl",
    pages = [
        "Home" => "index.md",
        "SAM and BAM" => "man/hts-files.md",
        "API Reference" =>  "man/api.md"
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)
deploydocs(
    repo = "github.com/BioJulia/XAM.jl.git",
    devbranch = "develop",
    push_preview = true
)
