


module PlotFunctions
	import  ..MasslistFunctions
	using PyPlot

	export massDefectPlot, bananaPlot, addClickToggle

	function massDefectPlot(masses, compositions, concentrations, colors, plotTitle, colorCodeTitle; dotSize = 10, maxMass = 450, maxDefect = 0.25, minConc = 0.02, sumformulas = false)
	  figure()

	  h2o = MasslistFunctions.createCompound(H=2,O=1, Hplus=0)
	  o = MasslistFunctions.createCompound(O=1, Hplus=0)

	  for i=1:length(masses)
	    m = masses[i]
	    adduct = h2o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="lightblue", zorder=1)
	    end
	    adduct = o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="red", zorder=1)
	    end
	    if sumformulas text(masses[i],masses[i]-round(masses[i]),sumFormulaStringFromCompositionArray(compositions[:,i]), color="grey", clip_on=true, verticalalignment="center", size=10, zorder=100) end
	  end

	  scatter(masses, masses-round.(masses),dotSize*log.(concentrations./minConc), colors, zorder=10, linewidths=0.5)
	  xlim(0,maxMass)
	  ylim(0,maxDefect)
	  cb=colorbar()
	  cb["ax"]["set_ylabel"](colorCodeTitle)
	  xlabel("Mass [amu]")
	  ylabel("Kendrick Mass Defect")
	  title(plotTitle)
	  grid("on")
	end

	function bananaPlot(xbins,ybins,meshdataXY;subplotAx=0)
	    println("Not implemented yet!!!")
	    #if subplotAx == 0
	    #    figure()
	    #    subplotAx = subplot(1,1,1)
	    #end
	    pcolormesh([DateTime(2018,05,05),DateTime(2018,05,06),DateTime(2018,05,07)],[1,7,9],[1 2 3; 4 5 6; 7 8 9][1:end-1,1:end-1])
	    #yscale("log")
	end

	function addClickToggle(plotAxis)
	    f = plotAxis.get_figure()
	    plotLines = plotAxis.get_lines()
	    legLines = plotAxis.get_legend().get_lines()
	    for ll in legLines
	        ll.set_picker(5)
	    end
	    linesDict = Dict(zip(legLines,plotLines))
	    function onLegendPick(event)
	        legLine = event.artist
	        plotLine = linesDict[legLine]
	        isVisible = plotLine.get_visible()
	        plotLine.set_visible(! isVisible)
	        if ! isVisible
	            legLine.set_alpha(1.0)
	        else
	            legLine.set_alpha(0.1)
	        end
	        f.canvas.draw()
	    end
	    f.canvas.mpl_connect("pick_event",onLegendPick)
	end
	function getVisibleTraces(axis)
		return [ axis.get_lines()[i].get_visible() for i in 1:length(axis.get_lines())]
	end

end


module CalibrationFunctions
	import ..InterpolationFunctions
	import PyPlot

	export generateCalibFactorTrace, generateBgTraces, interpolateBgTraces
	function getMeanOfQuantile(samples,quant)
	    q = quantile(samples,quant)
	    return mean(samples[samples.<q])
	end
	function generateBgTraces(times, traces; slices=10, quant=0.05)
	    dt = false
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
		dt = true
	    end
	    l=size(traces,1)
	    w=size(traces,2)
	    m=slices
	    bgtraces = Array{Float64}(size(traces,1), size(traces,2))
	    bgtimes = [mean(times[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),1]) for n = 1:m]
	    for i=1:w
		bgtraces[:,i] = InterpolationFunctions.interpolate(times, bgtimes, [getMeanOfQuantile(traces[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),i],0.08) for n = 1:m])
	    end
	    if dt
		bgtimes = Dates.unix2datetime.(bgtimes)
	    end
	    return bgtraces
	end

	function interpolateBgTraces(times, bgtimes, bgvalues)
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
	    end
	    if typeof(bgtimes[1]) == DateTime
		bgtimes = Dates.datetime2unix.(bgtimes)
	    end

	    nMasses = size(bgvalues,2)
	    bgtraces = Array{Float64}(size(times,1), nMasses)

	    for i=1:nMasses
		bgtraces[:,i] = InterpolationFunctions.interpolate(times, bgtimes, bgvalues[:,i])
	    end

	    return bgtraces
	end

	function generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)
	    if typeof(traceTimes[1]) == DateTime
		println("Convertung traceTimes to unixtime")
		traceTimes = Dates.datetime2unix(traceTimes)
	    end
	    if typeof(calibTimes[1]) == DateTime
		println("Convertung calibTimes to unixtime")
		calibTimes = Dates.datetime2unix(calibTimes)
	    end
	    if typeof(transitionTimes[1]) == DateTime
		println("Convertung transitionTimes to unixtime")
		transitionTimes = Dates.datetime2unix(transitionTimes)
	    end
	  # Stuff might be unsorted
	  traceTimesSorted = sort(traceTimes)
	  sortedCalibIndices = sortperm(calibTimes)
	  calibTimesSorted = calibTimes[sortedCalibIndices]
	  calibValuesSorted = calibValues[sortedCalibIndices]
	  transitionTimesSorted = sort(transitionTimes)

	  piecewiseTimes = Array{Float64}(0)
	  piecewiseValues = Array{Float64}(0)
	  currentValue = 0
	  currentCalibIndex = 1
	  currentTransitionIndex = 1
	  currentPiecewiseIndex = 1

	  # Drop transitions before start of trace
	  while (transitionTimesSorted[currentTransitionIndex] < traceTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
	    currentTransitionIndex += 1
	  end
	  # Drop transitions before first calibration
	  while (transitionTimesSorted[currentTransitionIndex] < calibTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
	    currentTransitionIndex += 1
	  end
	  println("First transition within data range: $currentTransitionIndex - $(transitionTimesSorted[currentTransitionIndex])")

	  # If no calib was done before start of trace, insert copy of first calib
	  println("No Calib before start of data range found, copying first successive calib")
	  push!(piecewiseTimes, minimum([traceTimes[1] calibTimes[1]]))
	  push!(piecewiseValues, calibValuesSorted[1])
	  currentPiecewiseIndex += 1

	  # Iterate over calibs and transitions
	  while (currentCalibIndex <= length(calibTimesSorted))
	    if currentTransitionIndex < length(transitionTimesSorted)
	      println("CurrCalib: $currentCalibIndex ($(calibTimesSorted[currentCalibIndex])), CurrTransition: $currentTransitionIndex ($(transitionTimesSorted[currentTransitionIndex]))")
	    end
	    # define piecewise values depending on what comes next, calib or transition
	    if ((currentTransitionIndex > length(transitionTimesSorted)) || (calibTimesSorted[currentCalibIndex] < transitionTimesSorted[currentTransitionIndex])) && (currentCalibIndex <= length(calibTimesSorted))
	      # add one point for a calib
	      println("Adding calib point  $(calibValuesSorted[currentCalibIndex]) at $(calibTimesSorted[currentCalibIndex])")
	      push!(piecewiseTimes, calibTimesSorted[currentCalibIndex])
	      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
	      currentPiecewiseIndex += 1
	      currentCalibIndex += 1
	    else
	      # skip all following transitions that come before a calibration, since there are no calibs in between
	      while (currentTransitionIndex < length(transitionTimesSorted)) && (transitionTimesSorted[currentTransitionIndex+1] < calibTimesSorted[currentCalibIndex])
		println("Skipping transition $currentTransitionIndex, there are no calibs in between!")
		currentTransitionIndex +=1
	      end
	      # add one point with last value before transition
	      println("Adding pre-transition point $(piecewiseValues[end]) at $(transitionTimesSorted[currentTransitionIndex])")
	      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] - 1e-99)
	      push!(piecewiseValues, piecewiseValues[end])
	      currentPiecewiseIndex += 1

	      # add one point with next value after transition
	      println("Adding post-transition point $(calibValuesSorted[currentCalibIndex]) at $(transitionTimesSorted[currentTransitionIndex])")
	      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] + 1e-99)
	      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
	      currentPiecewiseIndex += 1

	      currentTransitionIndex += 1
	    end
	  end
	  PyPlot.figure()
	  PyPlot.plot(piecewiseTimes, piecewiseValues, "x-")

	  return InterpolationFunctions.interpolate(traceTimes, piecewiseTimes, piecewiseValues)
	end
end



module ExportFunctions
	using ..InterpolationFunctions
	using Dates
	import ..MasslistFunctions

	export exportTracesCSV, exportTracesCSVLossCorr, toMatlabTime, fromMatlabTime

	function exportTracesCSV(saveFolderPath, elementNames, compositions, times, traces; average=0)
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  f = open("$saveFolderPath/ptr3compositions.txt", "w")
	  writedlm(f, hcat(["Mass" "SumFormula"],reshape(elementNames,(1,length(elementNames)))))
	  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions'))
	  close(f)
	  f = open("$saveFolderPath/ptr3traces.csv", "w")
	  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    writedlm(f, hcat(times ,traces))
	  else
	    writedlm(f, hcat(averageSamples(times,average) ,averageSamples(traces,average)))
	  end
	  close(f)
	end

	function exportTracesCSVLossCorr(saveFolderPath, elementNames, compositions, times, traces, lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes; average=0)
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  f = open("$saveFolderPath/ptr3compositions.txt", "w")
	  writedlm(f, hcat(["Mass"	"SumFormula"],reshape(elementNames,(1,length(elementNames))),["LossFactor"	"LossFactorError"	"CorrFactor"	"CorrFactorErr"	"CorrNotes"]))
	  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions', lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes))
	  close(f)
	  f = open("$saveFolderPath/ptr3tracesInletLossCorr.csv", "w")
	  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    writedlm(f, hcat(times , (corrfactor.*traces' )' ))
	  else
	    writedlm(f, hcat(averageSamples(times,average) ,(corrfactor.*(averageSamples(traces,average))' )' ))
	  end
	  close(f)
	end

	################ EXAMPLE from Matlab: 737551.673515479 should be 06-May-2019 16:09:51 #############
	function toMatlabTime(t::Dates.DateTime)
	    timespan = ((t+Dates.Day(1)) - Dates.DateTime(0,1,1,0,0,0))
	    millis::Float64 = Dates.value(timespan)
	    return millis/24/3600/1000
	end

	function fromMatlabTime(timestamp::Number)
	    days=Int(floor(timestamp))
	    millisecondsRemainder = Int(round((timestamp-days)*24*3600*1000))
	    return Dates.DateTime(0,1,1,0,0,0)+Dates.Day(days)+Dates.Millisecond(millisecondsRemainder)-Dates.Day(1)
	end

end
