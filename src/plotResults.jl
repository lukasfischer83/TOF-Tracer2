
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2

file = joinpath(pwd(),"ExampleFiles","TOFDATA","results","_result.hdf5")

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

massesToPlot = [
# Examples for selecting what to plot:
MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(C=10,H=16,O=2)
massLibrary.APINENE[1]
]


tracesFig = PlotFunctions.plotTracesFromHDF5(file;
				    	massesToPlot=massesToPlot,
				   	plotHighTimeRes = false,
				    	smoothing = 1,
				    	backgroundSubstractionMode = backgroundSubstractionMode,
				    	bg = (backgroundStart,backgroundEnd),
				    	timedelay = Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
				    	isobarToPlot = 0,
				    	plotsymbol = ".-"
				    	)


