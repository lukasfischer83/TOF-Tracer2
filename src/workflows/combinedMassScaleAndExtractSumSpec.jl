
using Distributed
import HDF5
import Statistics
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
import PyPlot
#import .TOFFunctions
#import .InterpolationFunctions

function correctMassScaleAndExtractSumSpec(
  filepath,
  masslistMasses,
  masslistElements,
  masslistElementsMasses,
  masslistCompositions,
  referenceFile,
  calibRegions; # Regions which are pattern matched for mass shifts
  filefilterRegexp=r"\.h5$",
  outputfilename="results/_result.hdf5",
  firstNFiles=0, # only analyze the first N files, set to 0 for all files
  lastNFiles=0, # only analyze the last N files, set to 0 for all files
  debuglevel=3,
  searchWidth = 0.7, # Width of Mass Scale Search Regions in AMU
  dynamicMassScaleCorrection = true,
  recalibInterval = 60,
  peakshapeRegions = 10,
  createTotalAvg = true, # Needed for manualPeakFitter, not needed for timetraces
  onlyUseAverages = false, # Fast mode using only total file averages instead of individual spectra in file
  storeSubSpectra = true,
  plotControlMass = true, # plot the control mass peak for each file after correction
  testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
  testRangeEnd = 137.5,
  massBorderCalculation = 1, # How to calculate borders? 0 = Cernter -0.1 to Center + 0.4,  1 = based on resolution and peak distance, 2 = constant bin width
  binWidth = 6,
  resolution = 7500,
  filePrecaching = true, # Precache whole next file. Only use when you have enough RAM and no SWAP has to be used.
  openWholeFile = true
  )


  if ! dynamicMassScaleCorrection #TODO: dummy fix, should be treated in workflow
      recalibInterval = 1e99
  end

  if (isfile(joinpath(filepath, outputfilename)))
    mv(joinpath(filepath, outputfilename), joinpath(filepath, "$outputfilename.bak"), force=true)
    if debuglevel > 0   println("Found existing result file, moving to '$outputfilename.bak'") end
  end

  #files = filter(r"2016-07-04.*\.h5", readdir(filepath)) #one day only
  allFiles = readdir(filepath)
  files = filter(s->occursin(filefilterRegexp, s), allFiles)


  nFiles = size(files,1)

  if firstNFiles != 0
    if length(files) > firstNFiles
      files = files[1:firstNFiles]
    end
  end

  if lastNFiles != 0
    if length(files) > lastNFiles
      files = files[length(files)-lastNFiles:end]
    end
  end

  nFiles = size(files,1)

  validFiles, timeSortIndices = TOFFunctions.validateHDF5Files(filepath, files)
  if !issorted(timeSortIndices)
      println("Reordering Files according to start acquisition time!")
  else
      println("File order seems fine, continuing...")
  end
  files = validFiles[timeSortIndices]
  nFiles = size(files,1)


  if (referenceFile == "")
    referenceFile = joinpath(filepath, files[1])
    if (debuglevel > 1)
      println("Reference File Autodetected (using first file in folder): $referenceFile")
    end
  end


  referenceSpectrum = TOFFunctions.getAvgSpectrumFromFile(referenceFile)
  referenceMassScaleMode, referenceMassScaleParameters = TOFFunctions.getMassCalibParametersFromFile(referenceFile)
  TOFFunctions.setMassScaleReferenceSpectrum(referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters, plotControlMass=plotControlMass, testRangeStart=testRangeStart, testRangeEnd=testRangeEnd)
  referenceMassAxis = []
  referenceMassAxis = HDF5.h5read(referenceFile, "FullSpectra/MassAxis")

  # Tofwerk does not update the mass axis after manual recalibration in TofDAQ Viewer.
  # ==> recalculate based on calib parameters.
  println("Recalculating mass axis")
  for i=1:length(referenceMassAxis)
      referenceMassAxis[i] = TOFFunctions.timebin2mass(i, referenceMassScaleMode, referenceMassScaleParameters)
  end

  totalAvgSpectrum = Array{Float64}(undef,length(referenceSpectrum))
  totalAvgSubSpectrum = Array{Float64}(undef,length(referenceSpectrum))
  totalMinSpectrum = Array{Float64}(undef,length(referenceSpectrum))
  totalMaxSpectrum = Array{Float64}(undef,length(referenceSpectrum))

  #@sync @parallel
  for i = 1 : length(referenceSpectrum)
    totalAvgSpectrum[i] = 0
    totalAvgSubSpectrum[i] = 0
    totalMinSpectrum[i] = 1e99
    totalMaxSpectrum[i] = 0
  end

  nMasses=length(masslistMasses)
  println("Calculating Stick CPS for $nMasses masses.")
  timebinshifts = zeros(Float64,(length(calibRegions),nFiles))
  intensities = zeros(Float64,(length(calibRegions),nFiles))
  monitorTimetrace = zeros(nFiles)

  #time = SharedArray{DateTime}(undef,nFiles)
  #stickcps=SharedArray{Float64}(undef,nMasses,nFiles)
  time = Array{DateTime}(undef,nFiles)
  stickcps=Array{Float64}(undef,nMasses,nFiles)

  # Calculate Integration Borders
  if (massBorderCalculation == 0)
  mlow = masslistMasses - 0.1
  mhigh = masslistMasses + 0.4
  elseif  (massBorderCalculation == 1)
  #mlow = masslistMasses .* resolution./(resolution+0.5)
  #mhigh = masslistMasses .* resolution./(resolution-0.5)
  massborders = MasslistFunctions.IntegrationBorders(masslistMasses; resolution=resolution)
  end


  if (plotControlMass == true)
    PyPlot.figure()
  end

  if (createTotalAvg == true)
      interpolatedSpectrum = Array{Float64}(undef,length(referenceSpectrum))
  end

  ############## Check output path and remove existing files #####################
  if (!isdir("$filepath/results"))
    mkdir("$filepath/results")
  end
  outfilepath = joinpath(filepath, outputfilename)
  if isfile(outfilepath)
    rm(outfilepath)
  end

  fid = HDF5.h5open(outfilepath, "w")

  ################################################################################

  ############## open outputfile and create extendable SumSpecs Array ############
  if (createTotalAvg == true)

    dsAvgSumWidth = length(referenceSpectrum)
    dspaceAvgSumSpecs = HDF5.dataspace((dsAvgSumWidth,1)::Dims, max_dims=(dsAvgSumWidth,typemax(Int64)))
    dtypeAvgSumSpecs = HDF5.datatype(Float32)
    dsetAvgSumSpecs = HDF5.create_dataset(fid, "SumSpecs", dtypeAvgSumSpecs, dspaceAvgSumSpecs; chunk=(dsAvgSumWidth,1), compress=3)
  end

  if (!onlyUseAverages)
    specCount = 0
    dsStickCpsWidth = nMasses
    dspaceStickCps = HDF5.dataspace((1,dsStickCpsWidth)::Dims, max_dims=(typemax(Int64),dsStickCpsWidth))
    dtypeStickCps = HDF5.datatype(Float32)
    dsetStickCps = HDF5.create_dataset(fid, "StickCps", dtypeStickCps, dspaceStickCps; chunk=(1,dsStickCpsWidth), compress=3)
    dsetStickCpsErr = HDF5.create_dataset(fid, "StickCpsErr", dtypeStickCps, dspaceStickCps, chunk=(1,dsStickCpsWidth), compress=3)

    dspaceTimes = HDF5.dataspace((1,)::Dims, max_dims=(typemax(Int64),))
    dtypeTimes = HDF5.datatype(Float64)
    dsetTimes = HDF5.create_dataset(fid, "Times", dtypeTimes, dspaceTimes; chunk=(1,))

  end
  ################################################################################
  badFiles = Array{String,1}()

  for j=1:nFiles
    totalPath = joinpath(filepath, files[j])
    if j < nFiles && filePrecaching
        totalPrecachePath = joinpath(filepath, files[j+1])
    else
        totalPrecachePath = ""
    end
    if debuglevel > 0   println("Processing File $j/$nFiles :  $totalPath") end

    fileIsBad = false;
    massAxis = []
    massAxis = HDF5.h5read(totalPath, "FullSpectra/MassAxis")

    time[j] = TOFFunctions.getTimeFromFile(totalPath)

    avgSpectrum = TOFFunctions.getAvgSpectrumFromFile(totalPath)
    newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(avgSpectrum, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)

    #if debuglevel > 0   println("Processing File $totalPath") end
    if (createTotalAvg == true && !fileIsBad)
      print("Spectrum interpolation for total average... ")
      indexesExact = TOFFunctions.mass2timebin(referenceMassAxis, referenceMassScaleMode, newParams)
      interpolatedSpectrum = InterpolationFunctions.interpolate(indexesExact,avgSpectrum)
      totalAvgSpectrum += interpolatedSpectrum
      # calculate min and max spectra
      #@sync @parallel
      for bin=1:length(interpolatedSpectrum)
        if interpolatedSpectrum[bin] > totalMaxSpectrum[bin]
          totalMaxSpectrum[bin] = interpolatedSpectrum[bin]
        end
        if interpolatedSpectrum[bin] < totalMinSpectrum[bin]
          totalMinSpectrum[bin] = interpolatedSpectrum[bin]
        end
      end
      # Write to hdf5, line by line, so there is no limit to number of files that can fit in RAM
      if (storeSubSpectra == true)
        currDims = size(dsetAvgSumSpecs)[2]
        dsetAvgSumSpecs[:,currDims] = interpolatedSpectrum
        if (j != nFiles)
          HDF5.set_extent_dims(dsetAvgSumSpecs, (dsAvgSumWidth,currDims+1)::Dims)
        end
      end
      println("DONE")
    end
    if debuglevel > 1   println() end
    print("Average Stick integration... ")
    for i=(1:nMasses)
      #if debuglevel > 0   println("Processing region mass($(mlow[i]):$(mhigh[i])) --> timebin($(round(TOFFunctions.mass2timebin(mlow[i],massCalibMode,newParams))):$(round(TOFFunctions.mass2timebin(mhigh[i],massCalibMode,newParams))))") end
      if (massBorderCalculation == 2)
        centerIndex = Int64(floor(TOFFunctions.mass2timebin(masslistMasses[i], referenceMassScaleMode, newParams)))
        if ((centerIndex-binWidth) > 0) && ((centerIndex+binWidth)<length(avgSpectrum))
          raw = sum(view(avgSpectrum, (centerIndex-binWidth : centerIndex + binWidth)))
        else
          raw = 0
        end
      else
        subIdxStartExact=TOFFunctions.mass2timebin(massborders.lowMass[i],referenceMassScaleMode,newParams)
        subIdxEndExact = TOFFunctions.mass2timebin(massborders.highMass[i],referenceMassScaleMode,newParams)
        raw = InterpolationFunctions.interpolatedSum(subIdxStartExact,subIdxEndExact,avgSpectrum)
      end
      stickcps[i,j] = raw
    end
    println("DONE")
    if (!onlyUseAverages)
      fileInternalLocalAvg = 0 #clear big data and references for garbage collection
      GC.gc()
      spectrumMultFactor = TOFFunctions.getSpecMultiplicator(totalPath)
      println("   spectrumMultFactor = $spectrumMultFactor")
      #subSpecStickCps = SharedArray{Float32}(undef,nMasses)
      subSpecStickCps = Array{Float32}(undef,nMasses)
      totalSubSpectra = TOFFunctions.getSubSpectraCount(totalPath)

      ################## prepare arrays for averaging #########################
      fileInternalLocalAvg = TOFFunctions.getSubSpectrumFromFile(totalPath,1, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
      if fileInternalLocalAvg == 0 # Skip file if the spectra could not be loaded
        println("Error: Could not load Subspectra from file $totalPath, skipping whole file!")
        continue
      end
      #println("referenceSpectrum is a $(summary(referenceSpectrum))")
      #println("SubSpec is a $(summary(fileInternalLocalAvg))")
      #println("spectrumMultFactor is a $(summary(spectrumMultFactor))")
      print("Sub spectrum Stick integration... ")

      fill!(fileInternalLocalAvg,0)
      fileInternalLocalAvgCount = 0
      print("Precalculating average for mass scale correction... ")
      ################## Get First Set of Spectra for first mass scale calib
      if totalSubSpectra >= recalibInterval
        # Prepare first calib beforehead
        fileInternalLocalAvg = TOFFunctions.getSubSpectrumFromFile(totalPath,1, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
        for avgIdx=2:minimum([recalibInterval totalSubSpectra])
          fileInternalLocalAvg += TOFFunctions.getSubSpectrumFromFile(totalPath,avgIdx, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
        end
        newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(fileInternalLocalAvg, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)
        ######
        fill!(fileInternalLocalAvg,0)
        fileInternalLocalAvgCount = 0
      end
      #######################################################################

      ################### Loop over Subspectra Start ########################
      ################### Averaging and Mass Scale Calib ####################
      print("Recalibrating at: ")
      specCount += totalSubSpectra
      for subSpecIdx=1:totalSubSpectra
        subSpectrum = TOFFunctions.getSubSpectrumFromFile(totalPath,subSpecIdx, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
        fileInternalLocalAvg += subSpectrum
        fileInternalLocalAvgCount += 1

        if (fileInternalLocalAvgCount > recalibInterval) || fileInternalLocalAvgCount == totalSubSpectra
          print("$subSpecIdx...")
          newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(fileInternalLocalAvg, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)

          indexesExact = TOFFunctions.mass2timebin(referenceMassAxis, referenceMassScaleMode, newParams)
          interpolatedSpectrum = InterpolationFunctions.interpolate(indexesExact,fileInternalLocalAvg)

          totalAvgSubSpectrum += interpolatedSpectrum*spectrumMultFactor

          fill!(fileInternalLocalAvg,0)
          fileInternalLocalAvgCount = 0
        end
        ################## Peak Integration ################################
        #@sync @parallel
        mbIndicesLowMass = TOFFunctions.mass2timebin(massborders.lowMass,referenceMassScaleMode,newParams)
        mbIndicesHighMass = TOFFunctions.mass2timebin(massborders.highMass,referenceMassScaleMode,newParams)
        Threads.@threads for i=(1:nMasses)
          if (mod(i,100) == 0)
            #println("Processing mass $(masslistMasses[i])")
          end
          if (massBorderCalculation == 2)
            centerIndex = Int64(floor(TOFFunctions.mass2timebin(masslistMasses[i], referenceMassScaleMode, newParams)))
            if ((centerIndex-binWidth) > 0) && ((centerIndex+binWidth)<length(avgSpectrum))
              raw = sum(view(subSpectrum, (centerIndex-binWidth : centerIndex + binWidth)))
            else
              raw = 0
            end
            subSpecStickCps[i]=raw
          else
            #if debuglevel > 0   println("Processing region mass($(mlow[i]):$(mhigh[i])) --> timebin($(round(TOFFunctions.mass2timebin(mlow[i],massCalibMode,newParams))):$(round(TOFFunctions.mass2timebin(mhigh[i],massCalibMode,newParams))))") end
            subSpecStickCps[i]= InterpolationFunctions.interpolatedSum(mbIndicesLowMass[i],mbIndicesHighMass[i],subSpectrum)
          end
        end
        rawTime = TOFFunctions.getSubSpectrumTimeFromFile(totalPath,subSpecIdx)
        currDimsTime = size(dsetTimes)[1]
        absTimeOfCurrSample = Dates.datetime2unix(time[j]) + rawTime
        dsetTimes[currDimsTime] = absTimeOfCurrSample
        dsetStickCps[currDimsTime,:] = subSpecStickCps*spectrumMultFactor
        if (currDimsTime == 1)
          deltaT = 5 #rawTime # Not correct???
        else
          deltaT = absTimeOfCurrSample - dsetTimes[currDimsTime-1][1]
        end
        #println("DeltaT: $deltaT")
        dsetStickCpsErr[currDimsTime,:] = sqrt.(abs.(subSpecStickCps./deltaT))
        if !((j==nFiles) && (subSpecIdx== TOFFunctions.getSubSpectraCount(totalPath)))
          HDF5.set_extent_dims(dsetTimes, (currDimsTime+1,)::Dims)
          HDF5.set_extent_dims(dsetStickCps, (currDimsTime+1,dsStickCpsWidth)::Dims)
          HDF5.set_extent_dims(dsetStickCpsErr, (currDimsTime+1,dsStickCpsWidth)::Dims)
        end
      end
      println("DONE")
    end
    println("###################################################################\n\n")
  end
  close(fid)

  stickcpsFlat = zeros(nFiles,nMasses)
  stickcpsErrFlat = zeros(nFiles,nMasses)

  timeFlat = zeros(nFiles)
  timeFlat = Dates.datetime2unix.(time)
  # We need more than one timestamp to calculate a difference
  if length(timeFlat) > 1
      deltaT = Statistics.median(diff(timeFlat))
  else
      deltaT = 1
  end

  totalAvgSpectrumFlat = zeros(length(totalAvgSpectrum))
  totalAvgSubSpectrumFlat = zeros(length(totalAvgSubSpectrum))


  for i=1:nMasses, j=1:nFiles
    stickcpsFlat[j,i] = stickcps[i,j]
    stickcpsErrFlat[j,i] = sqrt(abs(stickcps[i,j])/deltaT)
  end

  for i=1:length(totalAvgSpectrum)
    totalAvgSpectrumFlat[i] = totalAvgSpectrum[i]/nFiles
  end
  if !onlyUseAverages
    for i=1:length(totalAvgSpectrum)
      totalAvgSubSpectrumFlat[i] = totalAvgSubSpectrum[i]/specCount
    end
  end
  HDF5.h5write(outfilepath, "AvgStickCps", stickcpsFlat)
  HDF5.h5write(outfilepath, "AvgStickCpsErr", stickcpsErrFlat)
  HDF5.h5write(outfilepath, "MassList", masslistMasses)
  println("Wrote $(length(masslistMasses)) masses to file.")
  HDF5.h5write(outfilepath, "MassListIntegrationBordersLow", massborders.lowMass)
  HDF5.h5write(outfilepath, "MassListIntegrationBordersHigh", massborders.highMass)

  HDF5.h5write(outfilepath, "MassListIdx", TOFFunctions.mass2timebin(masslistMasses, referenceMassScaleMode, referenceMassScaleParameters))
  HDF5.h5write(outfilepath, "MassListIntegrationBordersIdxLow", TOFFunctions.mass2timebin(massborders.lowMass, referenceMassScaleMode, referenceMassScaleParameters))
  HDF5.h5write(outfilepath, "MassListIntegrationBordersIdxHigh", TOFFunctions.mass2timebin(massborders.highMass, referenceMassScaleMode, referenceMassScaleParameters))
  HDF5.h5write(outfilepath, "MassCalibMode", referenceMassScaleMode)
  HDF5.h5write(outfilepath, "MassCalibParameters", referenceMassScaleParameters)

  HDF5.h5write(outfilepath, "AvgStickCpsTimes", timeFlat)
  HDF5.h5write(outfilepath, "MassAxis", referenceMassAxis)
  if onlyUseAverages
    HDF5.h5write(outfilepath, "AvgSpectrum", totalAvgSpectrumFlat)
  else
    HDF5.h5write(outfilepath, "AvgSpectrum", totalAvgSubSpectrumFlat)
  end

  if (createTotalAvg == true)
    HDF5.h5write(outfilepath, "SumSpecMax", totalMaxSpectrum)
    HDF5.h5write(outfilepath, "SumSpecMin", totalMinSpectrum)
  end
  HDF5.h5write(outfilepath, "ElementNames", masslistElements)
  HDF5.h5write(outfilepath, "ElementMasses", masslistElementsMasses)
  compFlat = Array{Int64}(undef,length(masslistCompositions[1]),length(masslistCompositions))
  [compFlat[:,i]=masslistCompositions[i] for i=1:length(masslistCompositions)]
  HDF5.h5write(outfilepath, "ElementalCompositions", compFlat)

  ###### free some mem #######
  TOFFunctions.clearCache()

end
