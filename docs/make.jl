using Documenter
using PruntTrajectories

makedocs(
    sitename = "PruntTrajectories.jl",
    modules = [PruntTrajectories],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://benchung.github.io/PruntTrajectories.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Example" => "example.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/BenChung/PruntTrajectories.jl.git",
    devbranch = "main",
    push_preview = true,
)
