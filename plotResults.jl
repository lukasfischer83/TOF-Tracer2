push!(LOAD_PATH, pwd())
include("startup.jl")

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
import .InterpolationFunctions
import  .MasslistFunctions
using .ResultFileFunctions

file = "ExampleFiles/TOFDATA/results/_result.hdf5"
plotHighTimeRes = true
plotFittedInsteadOfSummed = true # Use multi peak fitted data instead of raw
smoothing = 1 #Average
reconstruct = false
plotsymbol = ".-"
isobarToPlot = 0
#timedelay = Dates.Hour(0) # CLOUD12, ...
timedelay = Dates.Hour(1) # CLOUDX, CLOUD11


if reconstruct == true
  import MultipeakFunctions
  # Used for Peak Reconstruction
  println("Reconstructing Spectrum from fitted single peaks...")
  massAxis = HDF5.h5read(file, "MassAxis")
  masslistElements = HDF5.h5read(file,"ElementNames")
  peakShapesCenterMass = HDF5.h5read(file, "MassDepPeakshapeCenterMasses")
  peakShapesY = HDF5.h5read(file, "MassDepPeakshape")
  counts = mean(traces,1)
  avgSpec = HDF5.h5read(file,"AvgSpectrum") - HDF5.h5read(file,"AvgBaseline")
  reconstructedSpectrum = reconstructSpectrum(massAxis, masses, masslistElements, compositions, counts, peakShapesCenterMass, peakShapesY)
  figure()
  semilogy(massAxis, avgSpec, label="Original")
  semilogy(massAxis, reconstructedSpectrum/10, label="Reconstruction")
  legend()
end

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2017,10,18,20)
backgroundEnd = DateTime(2017,10,18,20,45)


include("manualMassLibrary.jl")

massesToPlot = [
#NAPHTHA[1]
#massFromComposition(C=10,H=6,O=2)
#massFromComposition(C=10,H=8,O=2)
#massFromComposition(C=10,H=8,O=3)
#massFromComposition(C=10,H=6,O=4)

#massFromComposition(C=8,H=6,O=2)

#masses[(masses.>185) & (masses.<186)]
#massFromComposition( C=8,H=12,O=6)
#massFromComposition(C=19,H=28,O=11)
#massFromComposition(C=9,H=14,O=4)
#H3O[1]
#H3OH2O[1]
#H3OH2OH2O[1]

APINENE[1]
#massFromComposition(C=8, H=10)
#massFromComposition(C=6, H=12, O=1)
#massFromComposition(C=4, H=8, O=1)

#DMS[1]
#DMSO[1]
#DMSO2[1]
#MSIA[1]






#NAPHTHA[1]
#NH3[1]
#DMA[1]
#TMA[1]
#ACETONITRILE[1]
#ISOPRENE[1]
#massFromComposition(C=2, H=5, N=1, O=1)
#massFromComposition(C=2, H=7, N=1, O=2)
#massFromComposition(C=4, H=6, N=2)
#massFromComposition(C=6, H=11, N=1)
#PINONALDEHYDE[1]
#PYRIDINE[1]
#massFromComposition(C=10,H=16,O=1)
#massFromComposition(C=20,H=32,O=7)
#massFromComposition(C=14,H=30,O=2)
#massFromComposition(C=16,H=32,O=2)
#massFromComposition(C=12,H=24,O=2)
#massFromComposition(C=22,H=42,O=4)
#massFromComposition(C=4,H=6,O=2)
#massFromComposition( C=9,H=14,O=3)
#massFromComposition( C=10,H=18,O=2)
#padhypan
#PINONALDEHYDEPAN[1]
#OrgNitratNO[1]
#NORPINONALDEHYDE[1]
#NORPINONALDEHYDEPAN[1]
#PINONICACID[1]
#PINICACID[1]
#ACETIC[1]
#ACETICFRAG[1]
#massFromComposition(C=16,H=40,O=5)
#massFromComposition(C=10,H=14,O=5)
#massFromComposition(C=19,H=28,O=9)
#massFromComposition(C=20,H=30,O=8)

#massFromComposition(C=10,H=16,O=3)
#massFromComposition(C=10,H=19,O=3,N=1)

############# GROWTH STUFFS ############
#massFromComposition(C=10,H=16,O=3)
#massFromComposition(C=10,H=14,O=3)
#massFromComposition(C=9,H=14,O=3)

#massFromComposition(C=10,H=16,O=4)
#massFromComposition(C=9,H=14,O=4)
########################################


#massFromComposition(C=10,H=14,O=3)
#massFromComposition(C=10,H=17,O=3,N=1)

#massFromComposition(C=10,H=14,O=6)

#Pinonaldehyde-Surface-Reaction-Fragments
#createCompound(C=4, H=6, O=2)[1]
#createCompound(C=6, H=10, O=2)[1]
#Pinonaldehyde-Surface-Reaction-Fragments + Amines
#createCompound(C=4, H=7, N=1, O=1)[1]
#masses[263]
#masses[270]
#createCompound(C=5, H=9, N=1, O=1)[1]
#createCompound(C=6, H=11, N=1, O=1)[1]

#BCARY[1]
#HEXANONE[1]
#Nitrates
#C15H23NO4[1]
#C13H19NO6[1]
#C15H25NO4[1]
#C15H21NO5[1]
#C15H23NO5[1]
#C15H23NO6[1]
#C15H17NO7[1]

#C15H22O2[1]
#C15H24O2[1]
#C15H24O3[1]
#C15H26O4[1]
#Slow / Sticky
#C14H20O3[1]
#C13H20O4[1]
#C14H22O4[1]
#C15H22O4[1]
#C15H24O4[1]
#C14H25NO4[1]
#C15H29NO6[1]
#dimer1[1]
#dimer2[1]
#dimer3[1]
#dimer4[1]
#dimer5[1]

# Hyde
#massFromComposition(C=5,H=9,O=1, N=1)
#massFromComposition(C=6,H=10,O=1)
#TMB[1]
#massFromComposition(C=3,H=6)
]

#massesToPlot = vcat(massesToPlot,masses[(masses.>isobarToPlot) & (masses.<isobarToPlot+1)])
#massesToPlot = vcat(massesToPlot,masses[compositions[1,:].>=18])
#massesToPlot = vcat(massesToPlot,masses[sortedIndices[1:3]])
#massesToPlot = vcat(massesToPlot,masses[sortedIndices[(allMeans[sortedIndices].>10) & (volatility.>0.6) & (volatility.<1)]])

measResult = ResultFileFunctions.loadResults(file, massesToLoad=massesToPlot, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed, massMatchTolerance=0.01)
if isobarToPlot != 0
  isobarResult = ResultFileFunctions.loadResults(file, massesToLoad=[isobarToPlot+0.3], massMatchTolerance=0.5, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed)
  measResult=joinResultsMasses(measResult, isobarResult)
end

for i=1:length(measResult.Times)
  measResult.Times[i] = measResult.Times[i] - timedelay
end

fig=figure()
ax = subplot(111)


if (backgroundSubstractionMode == 0)
  background=0
elseif (backgroundSubstractionMode == 1)
  background = minimum(InterpolationFunctions.averageSamples(traces,smoothing),1)
elseif backgroundSubstractionMode == 2
  background = mean(measResult.Traces[(measResult.Times.>backgroundStart) & (measResult.Times.<backgroundEnd),:],1)
end

bgCorrectedTraces = measResult.Traces.-background

semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)
#semilogy(Dates.unix2datetime(InterpolationFunctions.averageSamples(Dates.datetime2unix(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)

startTimeString = Dates.format(measResult.Times[1],"yyyy/mm/dd")
endTimeString = Dates.format(measResult.Times[end],"yyyy/mm/dd")
title("$startTimeString - $endTimeString")
xlabel("Time [UTC]")
ylabel("Signal estimate [ppb]")

legStrings = []
for i = 1:length(measResult.MasslistMasses)
  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i]))")
end

box = ax[:get_position]()
#ax[:set_position]([box[:x0], box[:y0], box[:width] * 0.8, box[:height]])
cols = 1

majorformatter = matplotlib[:dates][:DateFormatter]("%m/%d %H:%M")
ax[:xaxis][:set_major_formatter](majorformatter)
#legend(legStrings, loc="center left", bbox_to_anchor=(1, 0.5), ncol=cols)
legend(legStrings)
#legend(legStrings, ncol=cols)
#ax[:yaxis][:grid]("on", which="minor")
grid()

#plotStages(measResult.Times[1], measResult.Times[end])
#plotStagesNames(measResult.Times[1], measResult.Times[end])


#=
include("$(pwd())/hydeFunctions.jl")

subplot(412, sharex=ax)
raint, raind = getAVAATrace("tconc", "HYY_DMPS")
s = (raint .> measResult.Times[1]) & (raint .< measResult.Times[end])
raint = raint[s]
raind = raind[s]
plot(raint, raind, "o-")

subplot(413, sharex=ax)
par = getAVAATrace("PAR", "HYY_META")
s = (par[1] .> measResult.Times[1]) & (par[1] .< measResult.Times[end])
par = (par[1][s] , par[2][s])

plot(par[1],par[2])

subplot(414)
winddir =  getAVAATrace("WDU336", "HYY_META")
s=(winddir[1] .> measResult.Times[1]) & (winddir[1] .< measResult.Times[end])
winddir = (winddir[1][s], winddir[2][s])

windspeed =  getAVAATrace("WSU336", "HYY_META")
s=(windspeed[1] .> measResult.Times[1]) & (windspeed[1] .< measResult.Times[end])
windspeed = (windspeed[1][s], windspeed[2][s])

windspeed=(windspeed[1][1:length(winddir[1])], windspeed[2][1:length(winddir[1])])
barbs(Dates.datetime2unix(winddir[1]), zeros(length(winddir[1])), 2*windspeed[2].*sin(3.14*winddir[2]/360), 2*windspeed[2].*cos(3.14*winddir[2]/360))
=#


#title("Biogenic VOCs")
tight_layout()
