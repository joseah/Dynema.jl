using Documenter
using Dynema

DocMeta.setdocmeta!(Dynema, :DocTestSetup, :(using Dynema); recursive = true)

makedocs(
    sitename = "Dynema.jl",
    authors = "Jose Alquicira-Hernandez",
    modules = [Dynema],
    workdir = @__DIR__,  # so `include("src/assets/....jl")` in @example/@setup blocks
                         # resolves the same way regardless of which page runs it
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://joseah.github.io/Dynema.jl",
        edit_link = "main",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,        # logo already includes the "Dynema" wordmark
        ansicolor = false,               # strip ANSI color instead of inlining <span> styles
        size_threshold = 600 * 2^10,      # 600 KiB hard limit (default 200 KiB)
        size_threshold_warn = 300 * 2^10, # 300 KiB warn limit (default 100 KiB)
        size_threshold_ignore = ["tutorials/simulation.md"],  # fallback: skip the check for this page entirely
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "CLI: Command-line overview" => "tutorials/command_line.md",
            "CLI: Main effect (context-independent)" => "tutorials/cli_main_effect.md",
            "CLI: Interaction effect (context-dependent)" => "tutorials/cli_interaction_effect.md",
            "CLI: Total effect (main + interaction)" => "tutorials/cli_total_effect.md",
            "Julia API: Download demo data" => "tutorials/simulation.md",
            "Julia API: Main effect (context-independent)" => "tutorials/main_effect.md",
            "Julia API: Interaction effect (context-dependent)" => "tutorials/interaction.md",
            "Julia API: Total effect (main + interaction)" => "tutorials/total_effect.md",
        ],
        "API Reference" => "functions.md",
        "Internals" => "internals.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/joseah/Dynema.jl.git",
    devbranch = "main",
    push_preview = true,
)


