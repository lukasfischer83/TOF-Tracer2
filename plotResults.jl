push!(LOAD_PATH, pwd())
include("startup.jl")

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions

file = "ExampleFiles/TOFDATA/results/_result.hdf5"
plotHighTimeRes = true # Plot every datapoint or only file averages
plotFittedInsteadOfSummed = true # Use multi peak fitted data instead of raw
smoothing = 1 # Average n samples, 1 for raw
plotsymbol = ".-"
isobarToPlot = 0
#timedelay = Dates.Hour(0) # CLOUD12, ...
timedelay = Dates.Hour(1) # CLOUDX, CLOUD11

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,17,14)
backgroundEnd = DateTime(2016,10,02,17,20)


include("manualMassLibrary.jl")

massesToPlot = [
# Examples for selecting what to plot:
MasslistFunctions.massFromComposition(C=10,H=16,O=2)
APINENE[1]
]

measResult = ResultFileFunctions.loadResults(file, massesToLoad=massesToPlot, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed, massMatchTolerance=0.01)
if isobarToPlot != 0
  isobarResult = ResultFileFunctions.loadResults(file, massesToLoad=[isobarToPlot+0.3], massMatchTolerance=0.5, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed)
  measResult=joinResultsMasses(measResult, isobarResult)
end

measResult.Times = measResult.Times .- timedelay

fig=figure()
ax = subplot(111)


if (backgroundSubstractionMode == 0)
  background=0
elseif (backgroundSubstractionMode == 1)
  background = minimum(InterpolationFunctions.averageSamples(measResult.Traces,smoothing),dims=1)
elseif backgroundSubstractionMode == 2
  import Statistics
  background = Statistics.mean(measResult.Traces[(measResult.Times.>backgroundStart) .& (measResult.Times.<backgroundEnd),:],dims=1)
end

bgCorrectedTraces = measResult.Traces .- background

semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)
#semilogy(Dates.unix2datetime(InterpolationFunctions.averageSamples(Dates.datetime2unix(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)

startTimeString = Dates.format(measResult.Times[1],"yyyy/mm/dd")
endTimeString = Dates.format(measResult.Times[end],"yyyy/mm/dd")
title("$startTimeString - $endTimeString")
xlabel("Time [UTC]")
ylabel("Signal [CPS]")

legStrings = []
for i = 1:length(measResult.MasslistMasses)
  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i]))")
end

box = ax.get_position()
cols = 1

majorformatter = matplotlib.dates.DateFormatter("%m/%d %H:%M")
ax.xaxis.set_major_formatter(majorformatter)
legend(legStrings)
grid()

tight_layout()
