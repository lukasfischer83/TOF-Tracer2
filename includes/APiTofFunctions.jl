

module APiTOFFunctions
	import HDF5
	import DSP
	using Distributed
	using Dates
	using Statistics
	import ..DatasetPreloader
	import ..InterpolationFunctions
	import PyPlot
	import LsqFit

	SPEC_CACHE_SIZE_LIMIT = 5e8


	export mass2timebin, timebin2mass, getMassCalibParametersFromFile, getSubSpectraCount, getSubSpectrumFromFile, getSpecMultiplicator, getSubSpectrumTimeFromFile, getAvgSpectrumFromFile, getTimeFromFile, validateHDF5Files, setMassScaleReferenceSpectrum, recalibrateMassScale

	debuglevel = 3

	mutable struct h5cache
	    filename
		filehandle
	    content
	end


	timeCache = h5cache("",0,0)
	spectraCache = h5cache("",0,0)


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
	    ret[i] = APiTOFFunctions.mass2timebin(mass[i],mode,parameters)
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
		fh = HDF5.h5open(filename,"r")
	    if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
		  try obj0 = fh["fitParams"]
			mcm = 2
		  catch
			mcm = 0
		  end
	      massCalibMode = mcm[1]
	      massCalibParameters = []
		  p1 = []
		  p2 = []
		  p3 = []
		  try obj2 = fh["fitParams"]
			  println("##############  Use fitParams data  ###############")
			  pV = HDF5.h5read(filename,"/fitParams")#fh["/fitParams"]#HDF5.h5read(referenceFile,"/PROCESSED/MassCalibration/ParametersData")
			  p1 = pV[1]#Statistics.median(pV[1,:]) # TODO: exchange median as soon as a better way of handling is found
			  p2 = pV[2]#Statistics.median(pV[2,:])
			  p3 = pV[3]
			  println("+++++++++ using fitParams: $pV +++++++++")
		  catch
			  println("##############  Use Viewer data  ###############")
			  try obj = fh["PROCESSED"]
				  println("##############  Use PROCESSED data  ###############")
 	 	          pV = fh["/PROCESSED/MassCalibration/ParametersData"]
		          p1 = Statistics.median(pV[1,:]) # TODO: exchange median as soon as a better way of handling is found
		          p2 = Statistics.median(pV[2,:])
				  p3 = true
			  catch
				  println("########## No PROCESSED data available. Use raw data instead. ##########")
		          pV = fh["/CALdata/Spectrum"]
		          p1 = pV[1,1]
		          p2 = pV[2,1]
			  	  p3 = true
			  end
	      end
		  if p3 == true
			  massCalibParameters = [p1[1] p2[1]]
		  else
			  massCalibParameters = [p1[1] p2[1] p3[1]]
		  end
	    else
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
		end
	  HDF5.close(fh)
	  return massCalibMode, massCalibParameters
	end

	function getAvgSpectrumFromFile(filename)
	  fh = HDF5.h5open(filename,"r")
	  if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
 		avgSpectrum = HDF5.h5read(filename, "SPECdata/AverageSpec")
  		return avgSpectrum
	  else
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

	  	#println("Read from File:\n   Writes: $H5NbrWrites, Bufs: $H5NbrBufs, WaveForms: $H5NbrWaveForms, TofPeriod: $H5TofPeriod, Segments: $H5NbrSegments, Blocks: $H5NbrBlocks, SampleInterval: $H5SampleInterval, SIS: $H5SingleIonSignal")
		println("Calculated:\n   Avg Spec Multiplier = $(H5SampleInterval./H5SingleIonSignal./H5inttime), Avg Spec Inttime = $H5inttime, Sum of Avg Spec = $(sum(avgSpectrum))")
	  	return avgSpectrum
	  end
	end

	function getSubSpectraCount(filename)
	  #println("Getting sub spectrum count from $filename")
	  fh = HDF5.h5open(filename,"r")
	  if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
	      ds = fh["/SPECdata/Intensities"]
	      l=1#size(ds)[1]
	      m=size(ds)[2]
	      HDF5.close(fh)
	      return l*m
	  else
	  	fh = HDF5.h5open(filename,"r")
	  	ds = fh["/FullSpectra/TofData"]
	  	l=size(ds)[3]
	  	m=size(ds)[4]
	  	HDF5.close(fh)
	  	return l*m
	  end
	end

	function getSubSpectrumFromFile(filename, index; openWholeFile = true, preloadFile = "")
		if openWholeFile # We have RRRRAAAAMMMMM !
			println("Filename: $(filename)")
		  if filename == ""
			return 0
		  end
		  if filename == DatasetPreloader.getPreloadFileName() # We discovered a preloaded file ready to fetch
		    if spectraCache.filename != filename
		      spectraCache.content = DatasetPreloader.load(filename)
			  #println("File has dims 1: $(size(spectraCache.content))")
		      spectraCache.filename = filename
			  GC.gc()
		    end
		  end
		  if spectraCache.filename != filename # we still have no spectrum, seems we have to load directly
		    spectraCache.content = DatasetPreloader.load(filename)
			#println("File has dims 2: $(size(spectraCache.content))")
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
		  #(l,m) = Base._ind2sub((size(spectraCache.content)[3],size(spectraCache.content)[4]), index)
		  subSpectrum = deepcopy(spectraCache.content[:,index])
#		  subSpectrum = deepcopy(spectraCache.content[:,1,l,m])
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
				  spectraCache.content = fh["/SPECdata/Intensities"] #fh["/FullSpectra/TofData"]
				  #println("File has dims: $(size(spectraCache.content))")
			  catch
				  println("Cannot open file $filename")
			  end
		  end
#		  (l,m) = Base._ind2sub((size(spectraCache.content)[3],size(spectraCache.content)[4]), index)
#		  subSpectrum = deepcopy(spectraCache.content[:,1,l,m])
		  subSpectrum = deepcopy(spectraCache.content[:,index])
		  GC.gc()
		  return subSpectrum
	  end
	end

	function getSpecMultiplicator(filename)
		fh = HDF5.h5open(filename,"r")
  		if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
    	return 1
  		else
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
	 	HDF5.close(fh)
	end

	function getSubSpectrumTimeFromFile(filename, index)
		fh = HDF5.h5open(filename,"r")
        if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
          SpectraTimes = HDF5.h5read(filename, "/SPECdata/Times")
          ds = fh["/SPECdata/Times"]
          timeCache.filename = filename
          timeCache.content = ds[1,index]#allSpectraTimes #ds[:,:]
          attributesRoot = HDF5.h5readattr(filename, "/")
          time = timeCache.content[1]*attributesRoot["Single Spec Duration (ms)"]./1000
         else
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
		end
		HDF5.close(fh)
	    return time
	end


	function getTimeFromFile(filename)
		fh = HDF5.h5open(filename,"r")
		if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
  		attributesRoot = HDF5.h5readattr(filename, "/")
  		if HDF5.exists(HDF5.attrs(fh), "FileCreatedTime_UTC")
	  	time_LabVIEWTimestamp = attributesRoot["FileCreatedTime_UTC"]
	  	println("### Using UTC time stamp ###")
  		else
	  	time_LabVIEWTimestamp = attributesRoot["FileCreatedTime"]
	  	println("### Using local time stamp ###")
  		end
  		tUnix = time_LabVIEWTimestamp[1]-((1970-1904)*365+(1970-1904)/4+0.5)*24*3600 # don't forget the leap years
  		t = Dates.unix2datetime(tUnix)
  		HDF5.close(fh)
		else
	    time_windowsTimestamp = HDF5.h5read(filename, "/AcquisitionLog/Log")[1].data[1]
	    tUnix = time_windowsTimestamp/(10.0*1000.0*1000.0)-11644473600.0
	    t = Dates.unix2datetime(tUnix)#[1]
		end
	    return t
	end

	function validateHDF5Files(filepath, files)
	  if debuglevel > 0   println("$(length(files)) files found, checking if valid.") end
	  validFiles = []
	  badFiles = []
	  startTimes = []
	  nFiles = size(files,1)
	  #println("Checking filepath $filepath ... ")
	  #println("Checking files $files ... ")
	  for j=1:nFiles
	     try
			 totalPath = joinpath(filepath, files[j])
			 if debuglevel > 1 print("Checking $totalPath ... ") end
			 if (!HDF5.ishdf5(totalPath))
		  		if debuglevel > 0   println("Bad File: $totalPath") end
		  		push!(badFiles,files[j])
			 else
		  		fh = HDF5.h5open(totalPath,"r")
		  		if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
			  		println("Ok, ioniTOF file format detected.")
					SpectraTimes = HDF5.h5read(totalPath, "/SPECdata/Times")
					if size(SpectraTimes[:,1])[1]==4
				  		n=3
			  		else
				  		n=2
			  		end
					ds = SpectraTimes[n,:]
					println("size of ds: $(size(ds))")
		  		else
		  			ds = fh["TimingData/BufTimes"]
	  	 		end
			end

			if (length(ds) > 0)
				if HDF5.exists(HDF5.attrs(fh), "InstrumentType")
					if (ds[end,1] > 1e-99) #(ds[end,][end,] > 1e-99) #(ds[end,end][end,end] > 1e-99) # Last timestamp seems to be very small on corrupted files TODO
						if debuglevel > 1 println("OK ioniAPiTOF file") end
						push!(validFiles,files[j])
						push!(startTimes,getTimeFromFile(totalPath))
					else
						push!(badFiles,files[j])
					end
				else
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
		  		end
			else
		  		if debuglevel > 0   println("Bad File: $totalPath") end
		  		push!(badFiles,files[j])
			end
	    catch e
			println("File seems corrupt: $(files[j]), ERROR: $e")
	    end
	  end
	  if debuglevel > 0   println("$(length(files)-length(validFiles)) files removed.") end
	  println("Finished validation of HDF5 Files.")
	  #println("Number of valid files: $(length(validFiles)).")
	  #println("Number of bad files: $(length(badFiles)).")
	  if length(validFiles)==0 	  println("No valid files...") end
	  return validFiles, sortperm(startTimes)
	end



	############# Mass recalibration stuff ##############################

	m_regionsToMatch = []
	m_regionMaxMatchCoeffs = []
	m_calibRegions = []
	m_searchWidth = 0
	m_referenceMassScaleMode = []#0
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
	  println("m_referenceMassScaleMode: $m_referenceMassScaleMode")

	  global m_plotControlMass = plotControlMass

	  for region in m_calibRegions
	    referenceIndexStart::Int64 = round(APiTOFFunctions.mass2timebin(region -  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
	    referenceIndexEnd::Int64 = round(APiTOFFunctions.mass2timebin(region +  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
	    regionToMatch = m_referenceSpectrum[referenceIndexStart:referenceIndexEnd]

	    push!(m_regionsToMatch,  regionToMatch)

	    # calc self correlation coeff max for later normalization of correlation quality
	    correlation = (DSP.xcorr(regionToMatch,regionToMatch))
	    maximumCorrelation = findmax(correlation)
	    intensity = maximumCorrelation[1]
	    push!(m_regionMaxMatchCoeffs,intensity)

	  end
	  if m_plotControlMass
	    global crIndStart = round(APiTOFFunctions.mass2timebin(testRangeStart,referenceMassScaleMode,referenceMassScaleParameters))
	    global crIndEnd = round(APiTOFFunctions.mass2timebin(testRangeEnd,referenceMassScaleMode,referenceMassScaleParameters))
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


	    indexStart::Int64 = round(APiTOFFunctions.mass2timebin(region - m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
	    indexEnd::Int64 = round(APiTOFFunctions.mass2timebin(region + m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
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
	    B[regionindex] = APiTOFFunctions.mass2timebin(region, m_referenceMassScaleMode,m_referenceMassScaleParameters)  + shift  + 0.5
	  end

	  if success
	    if debuglevel > 3   println("A: $A") end
	    if debuglevel > 3   println("B: $B") end
		if m_referenceMassScaleMode == 0
		  println("A: $A, B: $B")
	      newParams = \(A,B)
		end
		if m_referenceMassScaleMode == 2
		  @. model_tb2m(x, p) = ((x-p[2])/p[1])^(1/p[3])
		  p0_tb2m = [15000.0, -60000.0, 0.5]#[minimum(maxPeakVal)/1.4,maximum(maxPeakVal)]
		  println(" Size of B: $(size(B[:])), size calibRegions: $(size(m_calibRegions[:]))")
		  fit_res = LsqFit.curve_fit(model_tb2m, B[:], m_calibRegions[:], p0_tb2m )
		  newParams = fit_res.param
		end
	    if debuglevel > 3   println("New parameters: $(newParams)") end
	  else
	    newParams = m_referenceMassScaleParameters
	    if debuglevel > 3   println("Using WRONG calib parameters.") end
	  end

	  if (m_plotControlMass == true)
		println("New parameters: $(newParams), refMassCalMode: $m_referenceMassScaleMode")
		m_MassScaleMode = m_referenceMassScaleMode # 0 # instead of using m_referenceMassScaleMode
	    indexesExact = APiTOFFunctions.mass2timebin(crOriginalMasses, m_MassScaleMode, newParams) # (crOriginalMasses, m_referenceMassScaleMode, newParams)
	    crNewInterpolatedValues = InterpolationFunctions.interpolate(indexesExact,spectrum)
	    PyPlot.plot(crOriginalMasses, crNewInterpolatedValues/maximum(crNewInterpolatedValues),".-")
	  end
	return newParams, success, timebinshifts, intensities
	end

end
