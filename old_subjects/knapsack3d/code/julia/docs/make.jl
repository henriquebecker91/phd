using Documenter
push!(LOAD_PATH,"../src/")
using Check3DPackings, Chen1995, Kurpel2018

makedocs(
  modules = [Check3DPackings, Chen1995, Kurpel2018],
  format = :html,
  sitename = "HBDModules.jl",
  doctest = true,
  html_prettyurls = false
)

