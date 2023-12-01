module TOFTracer2

using Distributed
using Dates
using DelimitedFiles
using PyPlot

export correctMassScaleAndExtractSumSpec, baselineAndPeakshape, deconvolute, ResultFileFunctions, MasslistFunctions, PlotFunctions, InterpolationFunctions, massLibrary

if (nprocs() < 2 && isdefined(Main, :usePrecaching) && usePrecaching)
    println("Adding process for file precaching!")
    addprocs(1)
end

scriptspath = @__DIR__
print("scriptspath: ",scriptspath, "\n")

if ! isdefined(Main, :FunctionsLoaded)
    #load includes
    include("includes/BaselineFunctions.jl")
    include("includes/InterpolationFunctions.jl")
    include("includes/DatasetPreloader.jl")
    include("includes/MasslistFunctions.jl")
    include("includes/PeakshapeFunctions.jl")
    include("includes/TofFunctions.jl")
    include("includes/MultipeakFunctions.jl")
    include("includes/ResultFileFunctions.jl")
    #misc module leftovers
    include("includes/TofTracer.jl")

    #Load Personal MassLibrary Definitions
    include("includes/manualMassLibrary.jl")

    #Load workflows
    include("$(scriptspath)/workflows/combinedMassScaleAndExtractSumSpec.jl")
    include("$(scriptspath)/workflows/peakShape.jl")
    include("$(scriptspath)/workflows/deconvolutionMatrix.jl")
    FunctionsLoaded = true
end

end
