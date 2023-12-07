module TOFTracer2

using Distributed
using Dates
using DelimitedFiles
using PyPlot

export correctMassScaleAndExtractSumSpec, correctMassScaleAndExtractSumSpecAPi, baselineAndPeakshape, deconvolute, ResultFileFunctions, MasslistFunctions, PlotFunctions, InterpolationFunctions, massLibrary

if (nprocs() < 2 && isdefined(Main, :usePrecaching) && usePrecaching)
    println("Adding process for file precaching!")
    addprocs(1)
end

scriptspath = @__DIR__
print("scriptspath: ",scriptspath, "\n")

if ! isdefined(Main, :FunctionsLoaded)
    #load includes
    include(joinpath("includes","BaselineFunctions.jl"))
    include(joinpath("includes","InterpolationFunctions.jl"))
    include(joinpath("includes","DatasetPreloader.jl"))
    include(joinpath("includes","MasslistFunctions.jl"))
    include(joinpath("includes","PeakshapeFunctions.jl"))
    include(joinpath("includes","TofFunctions.jl"))
    include(joinpath("includes","MultipeakFunctions.jl"))
    include(joinpath("includes","ResultFileFunctions.jl"))
    include(joinpath("includes","PlotFunctions.jl"))
    include(joinpath("includes","ExportFunctions.jl"))
    include(joinpath("includes","CalibrationFunctions.jl"))

    #Load Personal MassLibrary Definitions
    include(joinpath("includes","manualMassLibrary.jl"))

    #Load workflows
    include(joinpath("workflows","combinedMassScaleAndExtractSumSpec.jl"))
    include(joinpath("workflows","peakShape.jl"))
    include(joinpath("workflows","deconvolutionMatrix.jl"))
    
    #Load stuff for APITOF?
    #include(joinpath("includes","APiTofFunctions.jl"))
    #include(joinpath("includes","MultipeakFunctionsAPi.jl"))
    #include(joinpath("workflows","combinedMassScaleAndExtractSumSpecAPi.jl"))
    
    FunctionsLoaded = true
end

end
