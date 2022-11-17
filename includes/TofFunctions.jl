

module TOFFunctions
	import HDF5
	import DSP
	using Distributed
	using Dates
	import ..DatasetPreloader
	import ..InterpolationFunctions
	import PyPlot

	SPEC_CACHE_SIZE_LIMIT = 5e8


	export clearCache, mass2timebin, timebin2mass, getMassCalibParametersFromFile, getSubSpectraCount, getSubSpectrumFromFile, getSpecMultiplicator, getSubSpectrumTimeFromFile, getAvgSpectrumFromFile, getTimeFromFile, validateHDF5Files, setMassScaleReferenceSpectrum, recalibrateMassScale

	debuglevel = 3

	mutable struct h5cache
	    filename
		filehandle
	    content
	end


	timeCache = h5cache("",0,0)
	spectraCache = h5cache("",0,0)

	function clearCache()
		timeCache = h5cache("",0,0)
		spectraCache = h5cache("",0,0)
	end

	

	function mass2timebin(mass::Number,mode,parameters)
	  if mode == 0
	    return parameters[1]*sqrt(mass) + parameters[2]
	  end
	  if mode == 1
	    return parameters[1]/sqrt(mass) + parameters[2]
	  end
	  if mode == 2
	    return parameters[1]*mass^parameters[3] + parameters[2]
	  end
	end
	function mass2timebin(mass::AbstractArray,mode,parameters)
	  ret = Array{Float64}(undef,length(mass))
	  Threads.@threads for i = 1:length(mass)
	    ret[i] = TOFFunctions.mass2timebin(mass[i],mode,parameters)
	  end
	  return ret
	end

	function timebin2mass(time::Number,mode,parameters)
	  if mode == 0
	    return ((time-parameters[2])/parameters[1])^2
	  end
	  if mode == 1
	    return (parameters[1]/(time-parameters[2]))^2
	  end
	  if mode == 2
	    return ((time-parameters[2])/parameters[1])^(1/parameters[3])
	  end
	end

	function timebin2mass(time::AbstractArray,mode,parameters)
	  ret = Array{Float64}(undef,length(time))
	  Threads.@threads for i = 1:length(time)
	    ret[i] = timebin2mass(time[i],mode,parameters)
	  end
	  return ret
	end

	function getMassCalibParametersFromFile(filename)
	  attributesFullSpectra = HDF5.h5readattr(filename, "/FullSpectra")

	  mcm = attributesFullSpectra["MassCalibMode"]
	  massCalibMode = mcm[1]
	  massCalibParameters = []

	  if (massCalibMode == 0)
	    p1 = attributesFullSpectra["MassCalibration p1"]
	    p2 = attributesFullSpectra["MassCalibration p2"]
	    massCalibParameters = [p1[1] p2[1]]
	  end

	  if (massCalibMode == 2)
	    p1 = attributesFullSpectra["MassCalibration p1"]
	    p2 = attributesFullSpectra["MassCalibration p2"]
	    p3 = attributesFullSpectra["MassCalibration p3"]
	    massCalibParameters = [p1[1] p2[1] p3[1]]
	  end
	  #HDF5.close(fh)
	  return massCalibMode, massCalibParameters
	end

	function getAvgSpectrumFromFile(filename)
	  #println("Getting avg spectrum from $filename")
	  attributesRoot = HDF5.h5readattr(filename, "/")
	  attributesFullSpectra = HDF5.h5readattr(filename, "/FullSpectra")
	  
	  H5NbrWrites = Float32.(attributesRoot["NbrWrites"])
	  H5NbrBufs = Float32.(attributesRoot["NbrBufs"])
	  H5NbrWaveForms = Float32.(attributesRoot["NbrWaveforms"])
	  H5TofPeriod = Float32.(HDF5.h5readattr(filename, "/TimingData")["TofPeriod"])
	  H5NbrSegments = Float32.(attributesRoot["NbrSegments"])
	  H5NbrBlocks = Float32.(attributesRoot["NbrBlocks"])

	  H5SampleInterval = Float32.(attributesFullSpectra["SampleInterval"]) .* 1e9
	  H5SingleIonSignal = Float32.(attributesFullSpectra["Single Ion Signal"])

	  H5inttime = H5NbrWaveForms .* H5TofPeriod .* H5NbrSegments .* H5NbrBlocks .* H5NbrBufs .* H5NbrWrites .* 1e-9
	  avgSpectrum = HDF5.h5read(filename, "FullSpectra/SumSpectrum") .* H5SampleInterval ./ H5SingleIonSignal ./ H5inttime

	  println("Read from File:\n   Writes: $H5NbrWrites, Bufs: $H5NbrBufs, WaveForms: $H5NbrWaveForms, TofPeriod: $H5TofPeriod, Segments: $H5NbrSegments, Blocks: $H5NbrBlocks, SampleInterval: $H5SampleInterval, SIS: $H5SingleIonSignal")
	  println("Calculated:\n   Avg Spec Multiplier = $(H5SampleInterval./H5SingleIonSignal./H5inttime)\nAvg Spec Inttime = $H5inttime")
	  return avgSpectrum
	end

	function getSubSpectraCount(filename)
	  #println("Getting sub spectrum count from $filename")
	  fh = HDF5.h5open(filename,"r")
	  ds = fh["/FullSpectra/TofData"]
	  l=size(ds)[3]
	  m=size(ds)[4]
	  HDF5.close(fh)
	  return l*m
	end

	function getSubSpectrumFromFile(filename, index; openWholeFile = true, preloadFile = "")
		if openWholeFile # We have RRRRAAAAMMMMM !
		  if filename == ""
			return 0
		  end
		  if filename == DatasetPreloader.getPreloadFileName() # We discovered a preloaded file ready to fetch
		    if spectraCache.filename != filename
		      spectraCache.content = DatasetPreloader.load(filename)
			  println("File has dims: $(size(spectraCache.content))")
		      spectraCache.filename = filename
			  GC.gc()
		    end
		  end
		  if spectraCache.filename != filename # we still have no spectrum, seems we have to load directly
		    spectraCache.content = DatasetPreloader.load(filename)
			println("File has dims: $(size(spectraCache.content))")
		    spectraCache.filename = filename
			GC.gc()
		  end
		  if preloadFile != "" && preloadFile != DatasetPreloader.getPreloadFileName() # nothing seems to be preloading, starting preload
			  if Distributed.nprocs() > 1
				  DatasetPreloader.preload(preloadFile)
			  end
			GC.gc()
		  end
		  if spectraCache.content == 0
			  return 0
		  end
		  (l,m) = Base._ind2sub((size(spectraCache.content)[3],size(spectraCache.content)[4]), index)
		  subSpectrum = deepcopy(spectraCache.content[:,1,l,m])
		  return subSpectrum
	  else # We don't have RAM, try to keep memory footprint low.
		  if spectraCache.filename != filename
			  #close old file if open
			  if spectraCache.filehandle != 0
				  try
					  println("closing file $filename")
					  HDF5.close(spectraCache.filehandle)
					  spectraCache.filehandle = 0
				  catch
					  println("Cannot close filehandle of $(spectraCache.filename)")
				  end
			  end
			  #open new file
			  try
				  println("opening file $filename for partial reading")
				  fh = HDF5.h5open(filename)
				  spectraCache.filehandle = fh
				  spectraCache.filename = filename
				  spectraCache.content = fh["/FullSpectra/TofData"]
				  println("File has dims: $(size(spectraCache.content))")
			  catch
				  println("Cannot open file $filename")
			  end
		  end
		  (l,m) = Base._ind2sub((size(spectraCache.content)[3],size(spectraCache.content)[4]), index)
		  subSpectrum = deepcopy(spectraCache.content[:,1,l,m])
		  GC.gc()
		  return subSpectrum
	  end
	end

	function getSpecMultiplicator(filename)
	  attributesRoot = HDF5.h5readattr(filename, "/")
	  attributesFullSpectra = HDF5.h5readattr(filename, "/FullSpectra")

	  #H5NbrWrites::Float32 = attributesRoot["NbrWrites"]
	  #H5NbrBufs::Float32 = attributesRoot["NbrBufs"]
	  H5NbrWaveForms::Float32 = attributesRoot["NbrWaveforms"][1]
	  H5TofPeriod::Float32 = HDF5.h5readattr(filename, "/TimingData")["TofPeriod"][1]
	  H5NbrSegments::Float32 = attributesRoot["NbrSegments"][1]
	  H5NbrBlocks::Float32 = attributesRoot["NbrBlocks"][1]

	  H5SampleInterval::Float32 = attributesFullSpectra["SampleInterval"][1] .* 1.0e9
	  H5SingleIonSignal::Float32 = attributesFullSpectra["Single Ion Signal"][1]

	  #H5inttime::Float32 = H5NbrWaveForms.*H5TofPeriod.*H5NbrSegments.*H5NbrBlocks .* 1.0e-9
	  H5inttime::Float32 = H5NbrWaveForms.*H5TofPeriod.* 1.0e-9 # laut Tanner keine Segments und Blocks
	  #println("H5inttime: $H5inttime")
	  #println("H5SampleInterval: $H5SampleInterval")
	  #println("H5SingleIonSignal: $H5SingleIonSignal")
	  return (H5SampleInterval./H5SingleIonSignal./H5inttime)[1] #orig
	  #return (1.0/H5inttime)/H5SingleIonSignal
	end

	function getSubSpectrumTimeFromFile(filename, index)
	    if filename != timeCache.filename
	      fh = HDF5.h5open(filename,"r")
	      allSpectraTimes = HDF5.h5read(filename, "/TimingData/BufTimes")
	      ds = fh["/TimingData/BufTimes"]
	      timeCache.filename = filename
	      timeCache.content = ds[:,:]
	    end
	    (l,m) = Base._ind2sub((size(timeCache.content)[1],size(timeCache.content)[2]), index)
	    #println("Getting sub spectrum Time ($l,$m) from $filename")
	    time = timeCache.content[l,m]
	    return time
	end


	function getTimeFromFile(filename)

	    time_windowsTimestamp = HDF5.h5read(filename, "/AcquisitionLog/Log")[1][1]
	    tUnix = time_windowsTimestamp/(10.0*1000.0*1000.0)-11644473600.0
	    #println("tUnix has size $(size(tUnix))")
	    t = Dates.unix2datetime(tUnix)#[1]
	    return t
	end

	function validateHDF5Files(filepath, files)
	  if debuglevel > 0   println("$(length(files)) files found, checking if valid.") end
	  validFiles = []
	  badFiles = []
	  startTimes = []
	  nFiles = size(files,1)

	  for j=1:nFiles
	     try
		totalPath = joinpath(filepath, files[j])
		if debuglevel > 1 print("Checking $totalPath ... ") end
		if (!HDF5.ishdf5(totalPath))
		  if debuglevel > 0   println("Bad File: $totalPath") end
		  push!(badFiles,files[j])
		else
		  fh = HDF5.h5open(totalPath,"r")
		  ds = fh["TimingData/BufTimes"]
		end

		if (length(ds) > 0)
		  if (ds[size(ds,1), size(ds,2)][1] > 1e-99) # Last timestamp seems to be very small on corrupted files
		    if debuglevel > 1 println("OK") end
			ds = fh["/FullSpectra/TofData"]
			if ds[1,1,1,size(ds,4)][1] >= 0
			    push!(validFiles,files[j])
			    push!(startTimes, getTimeFromFile(totalPath))
			else
				push!(badFiles,files[j])
			end
		  end
		else
		  if debuglevel > 0   println("Bad File: $totalPath") end
		  push!(badFiles,files[j])
		end
		HDF5.close(fh)
	    catch e
		println("File seems corrupt: $(files[j]): $e")
	    end
	  end
	  if debuglevel > 0   println("$(length(files)-length(validFiles)) files removed.") end
	  println("Finished validation of HDF5 Files.")
	  return validFiles, sortperm(startTimes)
	end



	############# Mass recalibration stuff ##############################

	m_regionsToMatch = []
	m_regionMaxMatchCoeffs = []
	m_calibRegions = []
	m_searchWidth = 0
	m_referenceMassScaleMode = 0
	m_referenceMassScaleParameters = []
	m_referenceSpectrum = []

	#plot stuff
	crIndStart = 0
	crIndEnd = 0
	crMassIndicesOriginal = []
	crOriginalMasses = []
	m_plotControlMass = false

	function setMassScaleReferenceSpectrum(referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters; plotControlMass=false, testRangeStart=0, testRangeEnd=0)
	  println("Setting mass scale reference spectrum.")
	  global m_regionsToMatch = []
	  global m_regionMaxMatchCoeffs = []
	  global m_referenceSpectrum = referenceSpectrum
	  global m_calibRegions = calibRegions
	  println("calibRegions: $m_calibRegions")

	  global m_searchWidth = searchWidth
	  global m_referenceMassScaleMode = referenceMassScaleMode
	  global m_referenceMassScaleParameters = referenceMassScaleParameters

	  global m_plotControlMass = plotControlMass

	  for region in m_calibRegions
	    referenceIndexStart::Int64 = round(TOFFunctions.mass2timebin(region -  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
	    referenceIndexEnd::Int64 = round(TOFFunctions.mass2timebin(region +  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
	    regionToMatch = m_referenceSpectrum[referenceIndexStart:referenceIndexEnd]

	    push!(m_regionsToMatch,  regionToMatch)

	    # calc self correlation coeff max for later normalization of correlation quality
	    correlation = (DSP.xcorr(regionToMatch,regionToMatch))
	    maximumCorrelation = findmax(correlation)
	    intensity = maximumCorrelation[1]
	    push!(m_regionMaxMatchCoeffs,intensity)

	  end
	  if m_plotControlMass
	    global crIndStart = round(TOFFunctions.mass2timebin(testRangeStart,referenceMassScaleMode,referenceMassScaleParameters))
	    global crIndEnd = round(TOFFunctions.mass2timebin(testRangeEnd,referenceMassScaleMode,referenceMassScaleParameters))
	    global crMassIndicesOriginal = collect(crIndStart:1:crIndEnd)
	    global crOriginalMasses = zeros(size(crMassIndicesOriginal))
	    for i = 1 : length(crOriginalMasses)
	      crOriginalMasses[i] = timebin2mass(crMassIndicesOriginal[i],referenceMassScaleMode,referenceMassScaleParameters)
	    end
	  end
	end

	function recalibrateMassScale(spectrum, referenceSpectrum, calibRegions, searchWidth, massCalibMode, massCalibParameters)
	success = true

	A = Array{Float64}(undef,length(m_calibRegions),2)
	A[:,2] .= 1
	B = Array{Float64}(undef,length(m_calibRegions),1)

	timebinshifts = zeros(length(m_calibRegions))
	intensities = zeros(length(m_calibRegions))

	#@sync @parallel
	for regionindex=1:length(m_calibRegions)
	    region = m_calibRegions[regionindex]


	    indexStart::Int64 = round(TOFFunctions.mass2timebin(region - m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
	    indexEnd::Int64 = round(TOFFunctions.mass2timebin(region + m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
	    regionToSearch = spectrum[indexStart:indexEnd]
	    correlation = (DSP.xcorr(convert(Array{Float64,1}, regionToSearch),convert(Array{Float64,1}, m_regionsToMatch[regionindex])))
	    maximumCorrelation = findmax(correlation)
	    if debuglevel > 3   println("maximumCorrelation at $(maximumCorrelation)") end
	    shift = 0;
	    intensity = 0;
	    if (maximumCorrelation[2]>1 && maximumCorrelation[2]<length(correlation))
	      inta=correlation[maximumCorrelation[2]-1]
	      intb=maximumCorrelation[1]
	      intc=correlation[maximumCorrelation[2]+1]
	      minAC = min(inta, intc)
	      shiftDeltaInterpolated = (intc - inta)/(inta + intb + intc - 3*minAC) # Tricky: -0.5 if a=b, +0.5 if b=c, +0.0 if a=c, interpol in between
	      #plot(correlation[maximumCorrelation[2]-10:maximumCorrelation[2]+10])
	      if debuglevel > 3   println("Interpolated Shift: $shiftDeltaInterpolated") end
	      shift = shiftDeltaInterpolated + maximumCorrelation[2] - (length(regionToSearch) + length(m_regionsToMatch[regionindex]))/2
	      intensity = maximumCorrelation[1]/m_regionMaxMatchCoeffs[regionindex]
	      timebinshifts[regionindex] = shift
	      intensities[regionindex] = intensity
	    else
	      timebinshifts[regionindex] = 0
	      intensities[regionindex] = 0

	    end
	    if debuglevel > 2   println("Mass $region found shifted by $shift timebins with correlation coeff $(intensity)") end
	    if (intensity < 0.05)
	      success = false
	      println("Could not correctly match m$region")
	    end
	    A[regionindex,1] = sqrt(region)
	    B[regionindex] = TOFFunctions.mass2timebin(region, m_referenceMassScaleMode,m_referenceMassScaleParameters)  + shift  + 0.5
	  end

	  if success
	    if debuglevel > 3   println("A: $A") end
	    if debuglevel > 3   println("B: $B") end
	    newParams = \(A,B)
	    if debuglevel > 3   println("New parameters: $(newParams)") end
	  else
	    newParams = m_referenceMassScaleParameters
	    if debuglevel > 3   println("Using WRONG calib parameters.") end
	  end

	  if (m_plotControlMass == true)
	    indexesExact = TOFFunctions.mass2timebin(crOriginalMasses, m_referenceMassScaleMode, newParams)
	    crNewInterpolatedValues = InterpolationFunctions.interpolate(indexesExact,spectrum)
	    PyPlot.plot(crOriginalMasses, crNewInterpolatedValues/maximum(crNewInterpolatedValues),".-")
	  end
	return newParams, success, timebinshifts, intensities
	end

end
