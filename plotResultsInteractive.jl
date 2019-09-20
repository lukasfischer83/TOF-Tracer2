push!(LOAD_PATH, pwd())
include("startup.jl")

#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
import PyPlot
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions


file = "/dataFast/HYYDE2016Data/running2/calibs/results/_result.hdf5"
availableTraces = ResultFileFunctions.loadResults(file,masslistOnly=true)
linewidth = 2

fig=PyPlot.figure()
concAx = PyPlot.subplot(1,1,1)
data =   ResultFileFunctions.getTraces(file,massIndices=[1])
concAx.semilogy(availableTraces.Times, data, linewidth=linewidth, label="m/z $(round(availableTraces.MasslistMasses[1],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,1]))","b")
concAx.set_ylabel("CPS")
concAx.set_xlabel("Time")
concAx.grid()
concAx.legend(loc=1)
lastPlottedIndex = 1
activeIndices = []

function changeLastPlotTo(i)
    global lastPlottedIndex, data
    if length(concAx.lines) > 0
        concAx.lines[end].remove()
    end
    data = ResultFileFunctions.getTraces(file,massIndices=[i])
    concAx.semilogy(availableTraces.Times, data, linewidth=linewidth, label="m/z $(round(availableTraces.MasslistMasses[i],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,i]))","b")
    concAx.legend(loc=1)
    global lastPlottedIndex = i
    println("Setting $lastPlottedIndex as last plot")
end

function addNewPlot()
    global lastPlottedIndex, data
    println("Adding $lastPlottedIndex to plots")
    if length(concAx.lines) > 0
        concAx.lines[end].remove()
    end
    data = ResultFileFunctions.getTraces(file,massIndices=[lastPlottedIndex])
    concAx.semilogy(availableTraces.Times, data, linewidth=linewidth, label="m/z $(round(availableTraces.MasslistMasses[lastPlottedIndex],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,lastPlottedIndex]))")
    concAx.semilogy(availableTraces.Times, data, linewidth=linewidth, label="m/z $(round(availableTraces.MasslistMasses[lastPlottedIndex],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,lastPlottedIndex]))","b")
    concAx.legend(loc=1)
end

function keyReleaseHandler(event)
    global activeIndices
    k=event.key
    if k == " "
        println("Adding Plot permanently")
        addNewPlot()
        push!(activeIndices,lastPlottedIndex)
    end
end
eventStepCumulator=0
function scrollHandler(event)
    global lastPlottedIndex, eventStepCumulator
    if event.step != 0
        println("Scrolled $(event[:step])")
    end
    eventStepCumulator += event.step # mouse pads scroll fractional
    if (eventStepCumulator >= 1) || (eventStepCumulator <=-0.6)
        plotIndex = clamp(lastPlottedIndex + eventStepCumulator, 1, length(availableTraces.MasslistMasses))
        changeLastPlotTo(Int64(floor(plotIndex)))
        eventStepCumulator = 0
    end
end

fig.canvas.mpl_connect("scroll_event", scrollHandler)
fig.canvas.mpl_connect("key_release_event", keyReleaseHandler)
