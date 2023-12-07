
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
