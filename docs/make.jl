using Documenter
using SplitVOF

makedocs(
    modules = [SplitVOF],
    authors = "PenguinxCutCell and contributors",
    repo = "https://github.com/PenguinxCutCell/SplitVOF.jl/blob/{commit}{path}#{line}",
    sitename = "SplitVOF.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/SplitVOF.jl",
        repolink = "https://github.com/PenguinxCutCell/SplitVOF.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Grid And Params" => "grid-and-params.md",
        "Initialization" => "initialization.md",
        "Reconstruction And Advection" => "reconstruction-and-advection.md",
        "Validation And Metrics" => "validation-and-metrics.md",
        "Reference" => "95-reference.md",
    ],
    pagesonly = true,
    warnonly = true,
)

deploydocs(
    repo = "github.com/PenguinxCutCell/SplitVOF.jl",
    push_preview = true,
)
