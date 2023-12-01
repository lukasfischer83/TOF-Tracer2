import HDF5
import PyPlot
import SparseArrays
import .MultipeakFunctionsAPi
import .MasslistFunctions
import .ResultFileFunctions
import .InterpolationFunctions

#include("masslistFunctions.jl")

function deconvoluteAPi(
  filepath;
  outputfilename="results/_result.hdf5",
  binWidth = 6,
  calcTransposed = false
  )
  #binWidth += 1 # Maybe crosstalk is calculated from "<" while summing is done from "<=" ?? Gives better crostalk removal.
  file = joinpath(filepath, outputfilename)

  peakShapesY = HDF5.h5read(file, "MassDepPeakshape")
  peakShapesY = peakShapesY ./ sum(peakShapesY,dims=1) # Normalize!

  peakShapesCenterMass = HDF5.h5read(file, "MassDepPeakshapeCenterMasses")
  totalAvgSpectrum = HDF5.h5read(file, "AvgSpectrum") - HDF5.h5read(file,"AvgBaseline")
  massAxis = HDF5.h5read(file, "MassAxis")
  masslistElements = HDF5.h5read(file,"ElementNames")
  compositionsOrig = HDF5.h5read(file, "ElementalCompositions")
  massesOrig = HDF5.h5read(file, "MassList")
  massBordersLowOrig = HDF5.h5read(file, "MassListIntegrationBordersLow")
  massBordersHighOrig = HDF5.h5read(file, "MassListIntegrationBordersHigh")
  massCenterIdxOrig = HDF5.h5read(file, "MassListIdx")
  massLowIdxOrig = HDF5.h5read(file, "MassListIntegrationBordersIdxLow")
  massHighIdxOrig = HDF5.h5read(file, "MassListIntegrationBordersIdxHigh")
  massScaleMode = HDF5.h5read(file, "MassCalibMode")
  massScaleParameters = HDF5.h5read(file, "MassCalibParameters")
  if unique(massesOrig) != massesOrig
    println("Multiple entries of the same mass --> will produce singular matrix!!!")
  end

  PyPlot.figure()
  ax = PyPlot.subplot(111)

    selector = (massesOrig .> 0)
    masses = massesOrig[selector]
    massBordersLow = massBordersLowOrig[selector]
    massBordersHigh = massBordersHighOrig[selector]
    massCenterIdx = massCenterIdxOrig[selector]
    massLowIdx = massLowIdxOrig[selector]
    massHighIdx = massHighIdxOrig[selector]

    compositions = compositionsOrig[:,selector]




    print("Populating matrix for inversion of linear system...")
    mtrx = MultipeakFunctionsAPi.calculateCrossTalkMatrix(masses, massCenterIdx, massLowIdx, massHighIdx, massScaleMode, massScaleParameters, compositions, peakShapesCenterMass, peakShapesY)
    stickRaw = [InterpolationFunctions.interpolatedSum(massLowIdx[i], massHighIdx[i], totalAvgSpectrum) for i=1:length(masses)]
    println(" DONE")



    print("Inverting Matrix...")
    #deconvolutionMatrix = SparseArrays.sparse(inv(mtrx)) # sparse on one processor is as fast as dense on mp, --> choosing less power consumption.
    deconvolutionMatrix = (inv(mtrx)) # there seems to be no multithread support for sparse matrix multiplication, still sparse is same speed
    println(" DONE")

    print("Applying deconvolution kernel...")
    counts = deconvolutionMatrix * Float64.(stickRaw)
    println(" DONE")

  print("Reconstructing Spectrum for visual check...")

  PyPlot.plot(massAxis,totalAvgSpectrum, "-o", label="Original", color="r")
  #reconstructedSpectrum = reconstructSpectrum(massAxis, masses[(masses.>158) & (masses.<162)], masslistElements, compositions[:,(masses.>158) & (masses.<162)], counts[(masses.>158) & (masses.<162)], peakShapesCenterMass, peakShapesY)
  reconstructedSpectrum = MultipeakFunctionsAPi.reconstructSpectrum(massAxis, massScaleMode, massScaleParameters, masses, compositions, counts, peakShapesCenterMass, peakShapesY)

  PyPlot.plot(massAxis, reconstructedSpectrum, label="Fit", color="b")
  PyPlot.plot(massAxis, totalAvgSpectrum-reconstructedSpectrum, label="Residual", color="g")
  assyErrorX = [(masses-massBordersLow)'; (massBordersHigh-masses)']
  y = InterpolationFunctions.interpolate(masses, massAxis, totalAvgSpectrum)
  PyPlot.errorbar(APiTOFFunctions.timebin2mass(massCenterIdx, massScaleMode, massScaleParameters),y,xerr=assyErrorX, fmt="o")
  PyPlot.errorbar(masses,y,xerr=assyErrorX, fmt="x")
  #=
  fittedPeaks = Array{Float64}(length(masses),2001)
  for i=1:length(masses)
    approxMassIndex = searchsortedfirst(massAxis,masses[i])
    fittedPeaks[i,:] = reconstructSpectrum(massAxis[approxMassIndex-300 : approxMassIndex+1700], masses[i], masslistElements, compositions[:,i], counts[i], peakShapesCenterMass, peakShapesY)
    plot(massAxis[approxMassIndex-300 : approxMassIndex+1700],fittedPeaks[i,:],"--", color="green")
  end
  =#

  PyPlot.legend()
  ax.set_ylim([minimum(totalAvgSpectrum),maximum(totalAvgSpectrum)])
  println(" DONE")



  ############ WRITE OUTPUT TO FILE ##############################################
  haveStickCps = false
  # Correct timetraces

  fh = HDF5.h5open(file,"r+")

  if haskey(fh,"CorrStickCps")
    HDF5.delete_attribute(fh,"CorrStickCps")
  end

  if haskey(fh,"CorrStickCpsErrors") 
    HDF5.delete_attribute(fh,"CorrStickCpsErrors")
  end

  if haskey(fh,"StickCps")
    haveStickCps = true

    # Create empty Dataspace
    nbrSpectra = ResultFileFunctions.getNbrTraceSamples(file)
    dset = HDF5.create_dataset(fh, "CorrStickCps", HDF5.datatype(Float32), HDF5.dataspace(nbrSpectra, length(masses)); chunk=(1,length(masses)), compress=3)

    toProcessLow = 0
    toProcessHigh = 0
    while toProcessHigh < nbrSpectra
      toProcessLow = toProcessHigh + 1
      toProcessHigh = toProcessLow + 9999
      if toProcessHigh > nbrSpectra
        toProcessHigh = nbrSpectra
      end
      print("Correcting spectrum $toProcessLow to $toProcessHigh of $nbrSpectra: Loading...")
      #samplesSubRange = convert(SharedArray, ResultFileFunctions.getTraceSamples(file,toProcessLow:toProcessHigh, raw=true)[:,selector])
      samplesSubRange = ResultFileFunctions.getTraceSamples(file,toProcessLow:toProcessHigh, raw=true)[:,selector]
      print("Deconvoluting...")
      for i=1:(toProcessHigh - toProcessLow + 1)
        #traces[i,:] = deconvolutionMatrix * traces[i,:]
        #tracesErrors[i,:] = abs(deconvolutionMatrix) * sqrt(abs(traces[i,:])/5
        dset[toProcessLow - 1 + i,:] = deconvolutionMatrix *samplesSubRange[i,:]
      end
      println("DONE")
    end
    #HDF5.h5write(file, "CorrStickCps", convert(Array,traces))
    #HDF5.h5write(file, "CorrStickCpsErrors", tracesErrors)
    HDF5.write_attribute(file, "CorrStickCps", dset)
  end
  HDF5.close(fh)

  if calcTransposed
      ResultFileFunctions.transposeStickCps(file)
  end

  fh = HDF5.h5open(file,"r+")
  if haskey(fh, "AvgStickCps")

    traces = HDF5.h5read(file, "AvgStickCps")[:,selector]
    tracesErrors = similar(traces)
    for i=1:size(traces,1)
      #traces[i,:] = deconvolutionMatrix \ traces[i,:]
      traces[i,:] = deconvolutionMatrix * traces[i,:]
      #tracesErrors[i,:] = abs(deconvolutionMatrix) * sqrt(abs(traces[i,:])/3600)
  end

  if haskey(fh,"CorrAvgStickCps")
    HDF5.delete_attribute(fh,"CorrAvgStickCps")
  end
  if haskey(fh,"CorrAvgStickCpsErrors")
    HDF5.delete_attribute(fh,"CorrAvgStickCpsErrors")
  end

  HDF5.h5write(file, "CorrAvgStickCps", traces)
  HDF5.h5write(file, "CorrAvgStickCpsErrors", tracesErrors)
  end
  
  HDF5.close(fh)

  
  return deconvolutionMatrix
end
