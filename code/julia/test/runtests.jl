#!/usr/bin/env julia

using Test
push!(LOAD_PATH, "../src")

include("ModelsTests.jl")
include("Check3DPackingsTests.jl")

