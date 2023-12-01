push!(LOAD_PATH,"../src/")

using Documenter
using TOFTracer2

DocMeta.setdocmeta!(TOFTracer2, :DocTestSetup, :(using TOFTracer2); recursive=true)


makedocs(;
	root    = "/home/wiebke/Documents/code_software/testingPackage_JL/TOFTracer2",
    	source  = "src",
    	build   = "build",
    	clean   = true,
    	doctest = true,
	modules = Module[TOFTracer2]
	)
