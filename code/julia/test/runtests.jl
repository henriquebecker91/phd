#!/usr/bin/env julia

using Test
push!(LOAD_PATH, "../src")

include("CornerPointsTests.jl")
#include("ModelsTests.jl")
#include("Check3DPackingsTests.jl")
#include("ExtremumPlanesTests.jl")

