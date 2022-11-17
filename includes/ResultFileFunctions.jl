module ResultFileFunctions
	import HDF5
	using Dates
	import Statistics
	import  ..MasslistFunctions
	import ..InterpolationFunctions

	export MeasurementResult, joinResultsTime, joinResultsMasses, getTraces, getTimetraces, transposeStickCps, getNbrTraces, getTraceSamples, getNbrTraceSamples, findChangingMasses, findVaryingMasses, saturationFromComposition, getIndicesInTimeframe

	#type MeasurementResult
	mutable struct MeasurementResult
	    Times
	    MasslistMasses
	    MasslistElements
	    MasslistElementsMasses
	    MasslistCompositions
	    Traces
	end

	function joinResultsTime(firstResult::MeasurementResult, secondResult::MeasurementResult)
	    if (length(firstResult.Times) > 0) & (length(secondResult.Times) > 0)
		if (firstResult.MasslistMasses == secondResult.MasslistMasses)
		    firstResult.Times = vcat(firstResult.Times, secondResult.Times)
		    firstResult.MasslistCompositions = hcat(firstResult.MasslistCompositions, secondResult.MasslistCompositions)
		    firstResult.Traces = hcat(firstResult.Traces, secondResult.Traces)
		else
		println("Masses did not match, could not merge results!")
		end
	    end
	return firstResult
	end

	function joinResultsMasses(firstResult::MeasurementResult, secondResult::MeasurementResult)
	    if (length(firstResult.MasslistMasses) > 0) & (length(secondResult.MasslistMasses) > 0)
		if (firstResult.Times == secondResult.Times)
		    firstResult.MasslistMasses = vcat(firstResult.MasslistMasses, secondResult.MasslistMasses)
		    firstResult.MasslistCompositions = hcat(firstResult.MasslistCompositions, secondResult.MasslistCompositions)
		    firstResult.Traces = hcat(firstResult.Traces, secondResult.Traces)
		else
		println("Times did not match, could not merge results!")
		end
	    end
	return firstResult
	end

	function loadResults(filename; useAveragesOnly = false, raw = false, startTime::DateTime = DateTime(0), endTime::DateTime = DateTime(3000), massesToLoad=Array{Float64,1}(), massMatchTolerance = 0.00001, masslistOnly = false)
	    println("Loading Times:")
	  if useAveragesOnly
	    timesUnix = HDF5.h5read(filename, "AvgStickCpsTimes")
	  else
	      try
		  timesUnix = HDF5.h5read(filename, "Times")
	      catch
		  println("Could not load high res time, maybe this is an average-only result file?")
		  return []
	      end
	  end

	println("Looking for errors in time axis")

	  ############## Check for Errors in time axis ###################
	#  if ! issorted(timesUnix)
	#    println("Time in result file is not in accending order, selection might not be complete!")
	#  end
	  if (isnan(timesUnix[1])) | (timesUnix[1] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[1] > Dates.datetime2unix(DateTime(2030,1,1)))
	      timesUnix[1] = timesUnix[1] - (timesUnix[3]-timesUnix[2])
	  end
	  if (isnan(timesUnix[length(timesUnix)]))| (timesUnix[length(timesUnix)] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[length(timesUnix)] > Dates.datetime2unix(DateTime(2030,1,1)))
	      timesUnix[length(timesUnix)] = timesUnix[length(timesUnix)] + (timesUnix[length(timesUnix)-1]-timesUnix[length(timesUnix)-2])
	  end

	  for i=2:length(timesUnix)-1
	      if (isnan(timesUnix[i])) | (timesUnix[i] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[i] > Dates.datetime2unix(DateTime(2030,1,1)))
		  timesUnix[i] = 0.5*(timesUnix[i-1]+timesUnix[i+1])
	      end
	  end

	  println("Selecting Times and Masses:")


	  masslistCompositions = HDF5.h5read(filename, "ElementalCompositions")
	  masslistMasses = HDF5.h5read(filename, "MassList")
	  masslistElements = HDF5.h5read(filename, "ElementNames")
	  masslistElementsMasses = HDF5.h5read(filename, "ElementMasses")

	  if masslistOnly
		  return MeasurementResult(Dates.unix2datetime.(timesUnix), masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions, Array{Float64}(undef,0,0))
	  end

	  if length(massesToLoad) > 0
	      selectionMassesIndices = Array{Int,1}()
	      for i=1:length(masslistMasses)
		for j=1:length(massesToLoad)
		  if (isapprox(masslistMasses[i], massesToLoad[j], atol=massMatchTolerance))
		    push!(selectionMassesIndices, i)
		  end
		end
	      end
	  else
	      selectionMassesIndices=1:length(masslistMasses)
	  end

	  selectionTimeIndexStart = searchsortedfirst(timesUnix, Dates.datetime2unix(startTime))
	  println("Start index: $selectionTimeIndexStart")
	  selectionTimeIndexEnd = searchsortedlast(timesUnix, Dates.datetime2unix(endTime))
	  println("End index: $selectionTimeIndexEnd")

	  if (selectionTimeIndexStart > 0) & (selectionTimeIndexEnd > selectionTimeIndexStart)
	    selectedTimesUnix = timesUnix[selectionTimeIndexStart:selectionTimeIndexEnd]
	  else
	    selectedTimesUnix = timesUnix
	    selectionTimeIndexStart = 1
	    selectionTimeIndexEnd = length(timesUnix)
	  end

	  ############# Load Traces here ##################

	  selectedTimes = Dates.unix2datetime.(selectedTimesUnix)

	  selectedMasslistMasses = masslistMasses[selectionMassesIndices]
	  selectedMassesCompositions = masslistCompositions[:,selectionMassesIndices]

	  if length(selectionMassesIndices) == 0
	    return MeasurementResult([], selectedMasslistMasses, masslistElements, masslistElementsMasses, Array{Float64}(undef,0,0), Array{Float64}(undef,0,0))
	  end

	  println("Loading Traces:")

	  traces = getTraces(filename, timeIndexStart=selectionTimeIndexStart, timeIndexEnd=selectionTimeIndexEnd, massIndices=selectionMassesIndices, raw=raw, useAveragesOnly=useAveragesOnly)

	  println("Loaded $(size(traces)) traces")
	  return MeasurementResult(selectedTimes, selectedMasslistMasses, masslistElements, masslistElementsMasses, selectedMassesCompositions, traces)
	end
###########################MeasurementResult end #########################

	function getIndicesInTimeframe(filename, startTime::DateTime, endTime::DateTime)
	  times = HDF5.h5read(filename, "Times")
	  return (1:length(times))[(times.>startTime) & (times.<endTime)]
	end

	function getTraces(filename; timeIndexStart=1, timeIndexEnd=0, massIndices=nothing, raw=false,  useAveragesOnly = false)
	  dsTexists = false

	  fh = HDF5.h5open(filename,"r")
	  if useAveragesOnly
	    if raw
	      ds = fh["AvgStickCps"]
	    else
	      ds = fh["CorrAvgStickCps"]
	    end
	  else
	    if raw
	      ds = fh["StickCps"]
	    else
	      ds = fh["CorrStickCps"]
	      if haskey(fh, "CorrStickCpsT")
		  	dsTexists = true
		  	println("A transposed dataset for faster loading is available.")
		  	dsT = fh["CorrStickCpsT"]
	      else
		  	println("Transposed dataset for faster loading is NOT available. You can create one with ResultFileFunctions.transposeStickCps(filename).")
	      end
	    end
	  end
	  result = 0

	  if timeIndexEnd == 0
	    timeIndexEnd = size(ds)[1]
	    println("No end time given, using all $(size(ds)[1]) elements")
	  end

	  if massIndices == nothing
	    massIndices = ( 1:size(ds)[2])
	    println("No mass selection given, selecting all $(size(ds)[2]) masses.")
	  end

	  if typeof(massIndices) == BitArray{1}
	    massIndices = ( 1:length(massIndices) )[massIndices]
	  end

	  if ((timeIndexStart>0) & (timeIndexEnd<=size(ds)[1])) & (minimum(massIndices)>0) & (maximum(massIndices)<=size(ds)[2])
	    # result = similar(ds[1,1],(timeIndexEnd-timeIndexStart+1),length(massIndices))
	    result = similar(ds,((timeIndexEnd-timeIndexStart+1),length(massIndices)))
	    
		println("Creating result Matrix with size $(size(result))")
	    if length(massIndices) < size(ds)[2]
		println("Subset of masses requested")
		if dsTexists
		    println("Loading from transposed dataset.")
		    for j=1:length(massIndices)
		        result[:,j] = dsT[massIndices[j],timeIndexStart:timeIndexEnd]
		    end
		else
		    for j=1:length(massIndices)
		        result[:,j] = ds[timeIndexStart:timeIndexEnd,massIndices[j]]
		        #println("Loaded $((timeIndexEnd-timeIndexStart)*j) samples")
		    end
		end
	    else
		println("Many masses requested, slicing time first for better performance")
		for j=timeIndexStart:timeIndexEnd
		    tmp = ds[j,:]
		    result[j-timeIndexStart+1,:] = tmp[massIndices]
		    #println("Loaded $(length(massIndices)*j) samples")
		end
	    end
	  end

	  close(fh)
	  return result
	end


	function getTimetraces(filename, indices; raw=false)
	  fh = HDF5.h5open(filename,"r")
	  if raw
	    ds = fh["StickCps"]
	  else
	    ds = fh["CorrStickCps"]
	  end

	  result = 0

	  if typeof(indices) == BitArray{1}
	    indices = ( 1:length(indices) )[indices]
	  end

	  if (typeof(indices) == Array{Int,1}) | (typeof(indices) == Array{Integer})
	    if (minimum(indices)>0) & (maximum(indices)<size(ds)[2])
	      result = similar(ds,(size(ds)[1],length(indices)))
	      for i=1:length(indices)
		result[:,i] = ds[:,i]
	      end
	    end
	  end

	  if typeof(indices) == UnitRange
	    if (minimum(indices)>0) & (maximum(indices)<size(ds)[2])
	      result = ds[:,indices]
	    end
	  end

	  close(fh)
	  return result
	end

	function transposeStickCps(filename)
	  fh = HDF5.h5open(filename,"r+")
	  if ! haskey(fh,"CorrStickCps")
	      println("CorrStickCps not present, skipping transpose.")
	      HDF5.close(fh)
	      return
	  end
	  if haskey(fh,"CorrStickCpsT")
	      HDF5.delete_attribute(fh,"CorrStickCpsT")
	  end
	  dsOld = fh["CorrStickCps"]
	  nbrSpectra, nbrMasses = size(dsOld)
	  dsNew = HDF5.create_dataset(fh, "CorrStickCpsT", HDF5.datatype(Float32), HDF5.dataspace(nbrMasses, nbrSpectra), chunk=(1,nbrSpectra), compress=3)
	  spectraAtOnce = 10000
	  for i=1:Int(floor(nbrSpectra/spectraAtOnce))
	      startIdx = (i-1)*spectraAtOnce + 1
	      endIdx = i*spectraAtOnce
	      print("Transposing Spectra $startIdx to $endIdx for faster loading...")
	      a=transpose(deepcopy(dsOld[startIdx:endIdx,:]))
	      for j=1:nbrMasses
		  dsNew[j,startIdx:endIdx] = a[j,:]
	      end
		  a=nothing
		  GC.gc()
	      println("DONE")
	  end
	  startIdx = (Int(floor(nbrSpectra/spectraAtOnce)))*spectraAtOnce + 1
	  endIdx = nbrSpectra
	  print("Transposing Spectra $startIdx to $endIdx for faster loading...")
	  a=transpose(deepcopy(dsOld[startIdx:endIdx,:]))
	  for j=1:nbrMasses
	      dsNew[j,startIdx:endIdx] = a[j,:]
	  end
	  a=nothing
	  println("DONE")
	  close(fh)
	  GC.gc()
	end


	function getTraceSamples(filename, indices; raw=false)
	  fh = HDF5.h5open(filename,"r")
	  if raw
	    ds = fh["StickCps"]
	  else
	    ds = fh["CorrStickCps"]
	  end
	  result = 0

	  if typeof(indices) == BitArray{1}
	    indices = ( 1:length(indices) )[indices]
	  end

	  if (typeof(indices) == Array{Int,1}) | (typeof(indices) == Array{Integer})
	    if (minimum(indices)>0) & (maximum(indices)<=size(ds)[1])
	      result = similar(ds,(length(indices),size(ds)[2]))
	      for i=1:length(indices)
		result[i,:] = ds[i,:]
	      end
	    end
	  end

	  if typeof(indices) == UnitRange{Int}
	    if (minimum(indices)>0) & (maximum(indices)<=size(ds)[1])
	      result = ds[indices,:]
	    end
	  end

	  close(fh)
	  if result == 0
	    println("Failed getting Traces Subset!")
	  end
	  return result
	end

	function getNbrTraces(filename)
	  fh = HDF5.h5open(filename,"r")
	  ds = fh["StickCps"]
	  nbr = size(ds)[2]
	  close(fh)
	  return nbr
	end

	function getNbrTraceSamples(filename)
	  fh = HDF5.h5open(filename,"r")
	  ds = fh["StickCps"]
	  nbr = size(ds)[1]
	  close(fh)
	  return nbr
	end

	function findChangingMasses(masses, compositions, traces, times, bgTimesSelection, signalTimesSelection; minOxygen = 0, sigmaThreshold=3, sorting = "mean", noNitrogen = true, onlySaneMasses = true, filterCrosstalkMasses=true)
	  bgTraces = traces[bgTimesSelection,:]
	  sigTraces = traces[signalTimesSelection,:]
	  meansBG = mean(bgTraces,1)[1,:]
	  meansBG[meansBG.<=0] = 1e-99
	  means = mean(sigTraces,1)[1,:] - meansBG
	  println("Means received: $(size(means))")
	  if size(bgTraces,1) < 3
	    println("Not enough BG points for standard deviation! plotHighTimeRes = true ?")
	  end

	  stderror = std(bgTraces,1)[1,:] ./ sqrt(size(bgTraces,2))

	  selMasses =  (means .> sigmaThreshold*stderror)
	  if filterCrosstalkMasses
	      # Filter out crosstalk masses
	      s = MasslistFunctions.filterMassListByContribution2(masses, means, 5000, 0.05)
	      selMasses = selMasses & s
	  end
	  println("\nRemoving $(length(masses[!s])) masses")

	  if noNitrogen == true
	    selMasses = selMasses & (compositions[5,:] .== 0)
	  end

	  selMasses = selMasses & (compositions[6,:] .>= minOxygen)
	  #selMasses = selMasses| ((masses.>137.1) & (masses.<137.15))

	  ## Filter only sane masses
	  if onlySaneMasses
	    selMasses = selMasses & (compositions[3,:] .> 1*compositions[1,:]) # H:C > 1.3
	    selMasses = selMasses & (compositions[3,:] .< 2.2*compositions[1,:]) # H:C < 2.2
	    selMasses = selMasses & (compositions[6,:] .< 1.5*compositions[1,:]) # O:C < 1.5
	  end
	  println("Selected $(sum(ones(length(masses))[selMasses])) masses!")

	  selIndices = linearindices(masses)[selMasses]
	  if sorting == "mean"
	      println("Sorting masses by mean value!")
	      sorter = sortperm(means[selMasses], rev=true)
	  elseif sorting == "mass"
	      println("Sorting masses!")
	    sorter = sortperm(masses[selMasses], rev=false)
	  else
	      println("Not sorting.")
	    sorter = linearindices(masses)
	  end

	  return selIndices[sorter], means[selIndices[sorter]], stderror[selIndices[sorter]]
	end

	function findVaryingMasses(masses, compositions, traces; minOxygen = 0, sigmaThreshold=3, sorting = "mean", noNitrogen = true, onlySaneMasses = true, filterCrosstalkMasses=true)
	  means = Statistics.mean(traces,dims=1)[1,:]
	  means[means .<= 0] .= 1e-99

	  println("Means received: $(size(means))")
	  if size(traces,1) < 3
	    println("Not enough points for variance! plotHighTimeRes = true ?")
	  end

	  tracesPrime = traces - InterpolationFunctions.smooth(traces, 0.01)
	  stderror = Statistics.std(tracesPrime,dims=1)[1,:] ./ sqrt(size(traces,2))
	  variances = Statistics.var(traces,dims=1)[1,:]
	  lastCount = length(masses)

	  selMasses =  (means .> sigmaThreshold*stderror)
	  println("Removed $(lastCount - sum(selMasses)) masses because mean<sigmaThreshold*stderr")
	  lastCount = sum(selMasses)
	  selMasses = selMasses .& (variances .> means .* sigmaThreshold)

	  println("Removed $(lastCount - sum(selMasses)) masses because variance<sigmaThreshold*mean")
	  lastCount = sum(selMasses)

	  if filterCrosstalkMasses
	      # Filter out crosstalk masses
	      s = MasslistFunctions.filterMassListByContribution2(masses, means, 5000, 0.05)
	      selMasses = selMasses .& s
	  end

	  println("\nRemoved $(lastCount - sum(selMasses)) masses because of crosstalk with neighbors")
	  lastCount = sum(selMasses)


	  println("\nRemoving $(length(masses[.!s])) masses")

	  if noNitrogen == true
	    selMasses = selMasses .& (compositions[5,:] .== 0)
	  end

	  selMasses = selMasses .& (compositions[6,:] .>= minOxygen)
	  #selMasses = selMasses .| ((masses .> 137.1) .& (masses .< 137.15))

	  ## Filter only sane masses
	  if onlySaneMasses
	    selMasses = selMasses .& (compositions[3,:] .> 1*compositions[1,:]) # H:C > 1.3
	    selMasses = selMasses .& (compositions[3,:] .< 2.2*compositions[1,:]) # H:C < 2.2
	    selMasses = selMasses .& (compositions[6,:] .< 1.5*compositions[1,:]) # O:C < 1.5
	  end
	  println("Selected $(sum(ones(length(masses))[selMasses])) masses!")
	  selIndices = (1:length(masses))[selMasses]
	  if sorting == "mean"
	      println("Sorting masses by mean value!")
	      sorter = sortperm(means[selMasses], rev=true)
	  elseif sorting == "mass"
	      println("Sorting masses!")
	    sorter = sortperm(masses[selMasses], rev=false)
	  else
	      println("Not sorting.")
	    sorter = linearindices(masses)
	  end

	  return selIndices[sorter], means[selIndices[sorter]], stderror[selIndices[sorter]]
	end

	function saturationFromComposition(compositions)
	  #return exp(10)*50000*exp.(-compositions[6,:]*5).*exp.(-compositions[1,:])
	  n_C = float(compositions[1,:])
	  n_O = float(compositions[6,:])
	  b_add = n_C
	  b_add[n_C .<= 15] = 0.92  # monomers
	  b_add[n_C .> 15] = 1.15 # dimers

	  return 10 .^((25-n_C)*0.475-n_O.*(2.3-b_add)-2 .*(n_C.*n_O./(n_C+n_O))*(-0.3))
	end
end
