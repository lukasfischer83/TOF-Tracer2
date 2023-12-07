
module MultipeakFunctionsAPi
	import ..PeakshapeFunctions
	import ..MasslistFunctions
	import ..APiTOFFunctions
	#import ..TOFFunctions
	import ..InterpolationFunctions

	export calculateCrossTalkMatrix, reconstructSpectrum, sumBins


	function createPeakPattern(mass, composition, massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
	    safetyMarginMasses = 2
	    peakShapeSize = size(peakShapesY,1)
	    # Get isotopes
	    if sum(composition) > 0
	      isotopeMasses, isotopeMasslistElements, isotopeCompositions, isotopeAbundances = MasslistFunctions.isotopesFromCompositionArray(composition)
	    else
	      #println("Unidentified Peak at mass $(mass)")
	      isotopeMasses = [mass]
	      isotopeAbundances = [1]
	    end
	    isotopeMassesIdx = APiTOFFunctions.mass2timebin(isotopeMasses, massScaleMode, massScaleParameters)
	    # Create Peak Pattern
	    # how many points minimum are needed?
	    nLocalPeakPatternPoints = Int(ceil(APiTOFFunctions.mass2timebin(maximum(isotopeMasses) + safetyMarginMasses, massScaleMode, massScaleParameters) - APiTOFFunctions.mass2timebin(minimum(isotopeMasses) - safetyMarginMasses, massScaleMode, massScaleParameters)))
	    localPeakPattern = zeros(nLocalPeakPatternPoints)

	    # where will be zero index? isotope could have lower mass!
	    startIndex = APiTOFFunctions.mass2timebin(minimum(isotopeMasses) - safetyMarginMasses, massScaleMode, massScaleParameters)
	    for isotope = 1:length(isotopeMasses)
	      localPeakshape = PeakshapeFunctions.getLocalPeakshape(isotopeMasses[isotope], peakShapesCenterMass, peakShapesY)
	      localPeakshape = localPeakshape * isotopeAbundances[isotope]
	      indexRelativeToStartIdx = isotopeMassesIdx[isotope] - startIndex
	      indexRelativeToPeakshapeStart = indexRelativeToStartIdx - (peakShapeSize-1)/2
	      InterpolationFunctions.addArraysShiftedInterpolated(localPeakPattern, localPeakshape, indexRelativeToPeakshapeStart-1)
	    end
	    return startIndex, localPeakPattern
	end


	function calculateCrossTalkMatrix(masses, centerindices, lowindices, highindices, massScaleMode, massScaleParameters, compositions, peakShapesCenterMass, peakShapesY)
	    mtrx = Array{Float64}(undef,length(masses), length(masses))
	    fill(mtrx,0)
	    for i = 1:length(masses)
		startIndex, localPeakPattern = createPeakPattern(masses[i], compositions[:,i], massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
		for j = 1:length(masses)
		    if ((masses[j] > masses[i]-1) && (masses[j] < masses[i]+3.5))
		        overlap = InterpolationFunctions.interpolatedSum(lowindices[j]-startIndex-1, highindices[j]-startIndex-1, localPeakPattern)
		        mtrx[j,i] = overlap
		    else
		        mtrx[j,i] = 0
		    end
		end
	    end
	    return mtrx
	end

	function sumBins(massAxis, binWidth, spectrum, masses)
	  stickRaw = Array{Float64}(undef,length(masses))
	  fill(stickRaw,0)
	  for i=1:length(masses)
	    midx = searchsortedfirst(massAxis, masses[i])-1
	    stickRaw[i] = sum(spectrum[midx-binWidth:midx+binWidth])
	  end
	  return stickRaw
	end

	function reconstructSpectrum(massAxis, massScaleMode, massScaleParameters, masses, compositions, counts, peakShapesCenterMass, peakShapesY)

	  reconstructedSpectrum = zeros(length(massAxis))
	  fill(reconstructedSpectrum, 0)
	  for i=1:length(masses)
	      startIndex, localPeakPattern = createPeakPattern(masses[i], compositions[:,i], massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
	      InterpolationFunctions.addArraysShiftedInterpolated(reconstructedSpectrum, counts[i]*localPeakPattern, startIndex)
	  end
	  return reconstructedSpectrum
	end


end
