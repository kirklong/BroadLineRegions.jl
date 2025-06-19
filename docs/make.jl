#!/usr/bin/env julia
push!(LOAD_PATH, "../src/")
using BLR
using Documenter 

makedocs(
    sitename = "BLR.jl",
    authors = "Kirk Long",
    modules = [BLR],
    pages = [
        "Home" => "index.md",
        "Usage and Examples" => "usage_examples.md",
        "API" => "api.md",]
)

deploydocs(
    repo = "github.com/kirklong/BLR.jl",
)