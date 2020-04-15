using Distributed
using Dates
using DelimitedFiles
using PyPlot

if (nprocs() < 2 && isdefined(Main, :usePrecaching) && usePrecaching)
    println("Adding process for file precaching!")
    addprocs(1)
end

scriptspath =  pwd()
if ! isdefined(Main, :FunctionsLoaded)
    #load includes
    include("includes/BaselineFunctions.jl")
    include("includes/InterpolationFunctions.jl")
    include("includes/DatasetPreloader.jl")
    include("includes/MasslistFunctions.jl")
    include("includes/PeakshapeFunctions.jl")
    include("includes/APiTofFunctions.jl")
    include("includes/MultipeakFunctionsAPi.jl")
    include("includes/ResultFileFunctions.jl")
    #misc module leftovers
    include("includes/TofTracer.jl")

    #Load Personal MassLibrary Definitions
    include("$(scriptspath)/manualMassLibrary.jl")

    #Load workflows
    include("$(scriptspath)/workflows/combinedMassScaleAndExtractSumSpecAPi.jl")
    include("$(scriptspath)/workflows/peakShape.jl")
    include("$(scriptspath)/workflows/deconvolutionMatrixAPi.jl")
    FunctionsLoaded = true
end
