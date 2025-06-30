#!/usr/bin/env julia
push!(LOAD_PATH, "../src/")
using BroadLineRegions
using Documenter 

makedocs(
    sitename = "BroadLineRegions.jl",
    authors = "Kirk Long",
    modules = [BroadLineRegions],
    pages = [
        "Home" => "index.md",
        "Usage and Examples" => "usage_examples.md",
        "API" => "api.md",],
    format = Documenter.HTML(
        sidebar_sitename = false
    )
)

deploydocs(
    repo = "github.com/kirklong/BroadLineRegions.jl",
)