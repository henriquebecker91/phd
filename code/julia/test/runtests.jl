#!/usr/bin/env julia

using Test
push!(LOAD_PATH, "../src")

include("Chen1995Tests.jl")
include("Check3DPackingsTests.jl")

