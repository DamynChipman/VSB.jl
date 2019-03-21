using Documenter, VSB

makedocs(
    format = Documenter.HTML(),
    sitename="VSB.jl",
    pages = ["index.md", "guide.md", "docstrings.md"],
    modules = [VSB]
)

# deploydocs(
#     repo = "github.com/camperD/VSB.jl.git",
# )
