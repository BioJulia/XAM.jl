using Pkg

Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, XAM

makedocs(
    format = Documenter.HTML(
        edit_branch = "develop"
    ),
    modules = [XAM, XAM.SAM, XAM.BAM],
    sitename = "XAM.jl",
    pages = [
        "Home" => "index.md",
        "SAM and BAM" => "hts-files.md",
        "API Reference" => [
            "Public" => "api/public.md"
        ],
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)
deploydocs(
    repo = "github.com/BioJulia/XAM.jl.git",
    devbranch = "develop"
)
