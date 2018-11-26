using Documenter
push!(LOAD_PATH,"../src/")
using Check3DPackings, Chen1995

makedocs(
  modules = [Check3DPackings, Chen1995],
  format = :html,
  sitename = "HBDModules.jl",
  doctest = true,
  prettyurls = false
)

