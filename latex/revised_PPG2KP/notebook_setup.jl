import Pkg
Pkg.activate(".")
Pkg.instantiate()
using WebIO
WebIO.install_jupyter_nbextension()
using DataFrames
using DataFramesMeta
# using Gadfly # Removed for now because downgrades DataFrames version
using IJulia
using Weave
using CSV
using Revise
using TableView
using Printf
using PrettyTables
using Formatting
