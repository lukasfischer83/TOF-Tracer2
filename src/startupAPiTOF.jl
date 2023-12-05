module TOFTracer2_APITOF

using Distributed
using Dates
using DelimitedFiles
using PyPlot

if (nprocs() < 2 && isdefined(Main, :usePrecaching) && usePrecaching)
    println("Adding process for file precaching!")
    addprocs(1)
end

scriptspath =  @__DIR__
print("scriptspath: ",scriptspath, "\n")

if ! isdefined(Main, :FunctionsLoaded)
    #load includes
    include(joinpath("includes","BaselineFunctions.jl"))
    include(joinpath("includes","InterpolationFunctions.jl"))
    include(joinpath("includes","DatasetPreloader.jl"))
    include(joinpath("includes","MasslistFunctions.jl"))
    include(joinpath("includes","PeakshapeFunctions.jl"))
    include(joinpath("includes","ResultFileFunctions.jl"))
    #misc module leftovers
    include(joinpath("includes","TofTracer.jl"))

    #Load Personal MassLibrary Definitions
    include("includes","manualMassLibrary.jl"))

    #Load stuff for APITOF
    include(joinpath("includes","APiTofFunctions.jl"))
    include(joinpath("includes","MultipeakFunctionsAPi.jl"))
    include(joinpath("workflows","combinedMassScaleAndExtractSumSpecAPi.jl"))
    FunctionsLoaded = true
end

end
