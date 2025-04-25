using TOML
using Documenter, XAM

authors = let
    project_path = joinpath(Base.pkgdir(XAM), "Project.toml")
    toml = TOML.parsefile(project_path)
    authors_with_email = join(toml["authors"], ", ")
    authors = replace(authors_with_email, r" <.*?>" => "")
    authors * ", The BioJulia Organisation, and other contributors."
end

makedocs(
    checkdocs = :all,
    # linkcheck = true,
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
    authors = authors
)

deploydocs(
    repo = "github.com/BioJulia/XAM.jl.git",
    devbranch = "develop",
    push_preview = true
)
