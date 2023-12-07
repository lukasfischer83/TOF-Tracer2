"The module PeakshapeFunctions is a module needed for the processing workflow and collects all functions that are related to determining the mass-dependent peakshape from the raw data."
module PeakshapeFunctions
	import PyPlot
	import SharedArrays
	import ..InterpolationFunctions
	import Statistics

	export findPeakIndices, calculatePeakshapes, getLocalPeakshape

	"""
		findPeakIndices(massAxis, avgSpectrum, baseline, baselineNoise; kwargs...)
		
	searches an average spectrum for peaks significantly higher than the baseline noise and returns their interpolated indices in the average spectrum
	"""
	  function findPeakIndices(massAxis, avgSpectrum, baseline, baselineNoise; noiseThreshold = 80, oddEven="both", signalLimit = 1)
	    totalMax = maximum(avgSpectrum)
	    peakIndices = Array{Float64,1}()
	    for j=searchsortedfirst(massAxis,10):length(massAxis)-1
	      if ( avgSpectrum[j] > baseline[j] + baselineNoise[j] * noiseThreshold * 400/(massAxis[j]^0.6)) && (avgSpectrum[j] > avgSpectrum[j-1]) && (avgSpectrum[j] > avgSpectrum[j+1]) && (avgSpectrum[j] < totalMax*signalLimit)
		  if (oddEven == "both") | ((oddEven == "even") & iseven(Int(round(massAxis[j],digits=0)))) | ((oddEven == "odd") & !iseven(Int(round(massAxis[j],digits=0))))
		      ############ put interpolated exact peak position in peak list ###
		      push!(peakIndices,InterpolationFunctions.interpolatedMax(j,avgSpectrum))
		  end
	      end
	    end
	    return peakIndices
	  end

	"""
		calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices; kwargs...)
		
	calculates the per-massregion peakshapes by normalizing and shifting all peaks within a massregion onto each other. 
	The resulting peakshape is the inner envelope (the lower quantile) of the peakshape data points.
	
	returns an array of the mass regions' center masses and an array of the corresponding peakshapes
	"""
	  function calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices; nbrMassRegions = 10, peakWindowWidth = peakWindowWidth, quantileValue = 0.05, regionStretch=1)
	    peakMasses = InterpolationFunctions.interpolate(peakIndices, massAxis)
	    peakValues = InterpolationFunctions.interpolate(peakIndices,baselineCorrectedAvgSpec)

	    peakShapesY = Array{Float64}(undef,2*peakWindowWidth+1, nbrMassRegions)
	    peakShapesCenterMass = Array{Float64}(undef,nbrMassRegions)

	    PyPlot.figure()
	    for massRegion = 1:nbrMassRegions
	      ax=PyPlot.subplot(Int(ceil(nbrMassRegions/4)),4,massRegion) #change Int() Leander

	      #peakshapeRangeStart = (massRegion-1) *peakMasses[end] / nbrMassRegions
	      #peakshapeRangeEnd = (massRegion) *peakMasses[end] / nbrMassRegions
	      peakshapeRangeStart = (massRegion-1)^regionStretch *peakMasses[end] / nbrMassRegions^regionStretch
	      peakshapeRangeEnd = (massRegion)^regionStretch *peakMasses[end] / nbrMassRegions^regionStretch

	      peakShapesCenterMass[massRegion] = (peakshapeRangeEnd + peakshapeRangeStart)/2
	      sel = (peakMasses.>peakshapeRangeStart) .& (peakMasses.<peakshapeRangeEnd) .& (peakIndices.> peakWindowWidth) .& (peakIndices .< length(baselineCorrectedAvgSpec)-peakWindowWidth)

	      peakIndicesInRegion = peakIndices[sel]
	      peakMassesInRegion = peakMasses[sel]
	      peakValuesInRegion = peakValues[sel]

	      peakWood = SharedArrays.SharedArray{Float64}(2*peakWindowWidth+1,length(peakMassesInRegion))
	      fill(peakWood, 0)
	      for i=1:length(peakMassesInRegion)
		  ############# interpolate exact peak position ###################
		lowIdx = Int64(floor(peakIndicesInRegion[i]))
		highIdx = lowIdx+1
		fract = peakIndicesInRegion[i]-lowIdx

		peakWood[:,i] = (
		(1-fract)*view(baselineCorrectedAvgSpec, lowIdx-peakWindowWidth : lowIdx+peakWindowWidth) +
		fract*view(baselineCorrectedAvgSpec, highIdx-peakWindowWidth : highIdx+peakWindowWidth))/peakValuesInRegion[i]
	      end

	      peakshapeY = Array{Float64}(undef,2*peakWindowWidth+1)
	      for i=1:2*peakWindowWidth+1
		#peakshapeY[i] = minimum(peakWood[i,:])
		values = peakWood[i,:]
		#upperThreshold = quantile(values,quantileValue*peakshapeRangeEnd/100)
		###### quantile based #############
		#upperThreshold = quantile(values,quantileValue*peakshapeRangeEnd/100)
		#peakshapeY[i] = mean(values[values.<upperThreshold])
		#lowest n perent based ############
		values = sort(values[values.>0])
		if length(values) > 0
		  peakshapeY[i] = Statistics.mean(values[1 : Int64(ceil(length(values)*quantileValue))])
		else
		  peakshapeY[i] = 0
		  if i == peakWindowWidth + 1
		    peakshapeY[i] = 1
		  end
		end
	      end
	      peakshapeY[peakshapeY.<0] .= 0
	      peakshapeY[peakWindowWidth+1] = 1
	      peakShapesY[:,massRegion] = peakshapeY / sum(peakshapeY)
	      if (size(peakWood,2) > 0)
		  PyPlot.semilogy(peakWood)
	      end
	      PyPlot.semilogy(peakshapeY, "o-", label="$peakshapeRangeStart - $peakshapeRangeEnd")
	      ax["set_title"]("$(round(peakshapeRangeStart,digits=1)) - $(round(peakshapeRangeEnd,digits=1))")
	    end
	    PyPlot.suptitle("Mass dependent Peakshapes")
	    return peakShapesCenterMass, peakShapesY
	  end

	"""
		getLocalPeakshape(mass, peakShapesCenterMass, peakShapesY)
		
	calculates the local peakshape by linearly interpolating between the nearest peakshapes-per-massregion, 
	with the contribution depending on the distance of the local mass to each mass-region's centermass
	
	returns an array containing the local peakshape
	"""
	  function getLocalPeakshape(mass, peakShapesCenterMass, peakShapesY)
	    peakShapeCenterMassIndexHigh = searchsortedfirst(peakShapesCenterMass, mass)
	    if (peakShapeCenterMassIndexHigh == 1)
	      localPeakshape = peakShapesY[:,1]
	    elseif (peakShapeCenterMassIndexHigh > length(peakShapesCenterMass))
	      localPeakshape = peakShapesY[:,end]
	    else
	      peakShapeCenterMassIndexLow = peakShapeCenterMassIndexHigh - 1
	      peakShapesCenterMassDistance = peakShapesCenterMass[peakShapeCenterMassIndexHigh] - peakShapesCenterMass[peakShapeCenterMassIndexLow]
	      contribLow =  (peakShapesCenterMass[peakShapeCenterMassIndexHigh] - mass) / peakShapesCenterMassDistance
	      contribHigh = 1 - contribLow
	      localPeakshape = peakShapesY[:,peakShapeCenterMassIndexLow]*contribLow + peakShapesY[:,peakShapeCenterMassIndexHigh]*contribHigh
	    end
	    return localPeakshape
	  end

end
