push!(LOAD_PATH, pwd())

import PyPlot
import HDF5
import .InterpolationFunctions
import .BaselineFunctions
import .PeakshapeFunctions

function baselineAndPeakshape(
  filepath;
  outputfilename="results/_result.hdf5",
  peakshapeRegions=8,
  peakshapeRegionStretch=0.5,
  peakshapeQuantileValue = 0.1,
  peakfindingNoiseThresholdValue = 25,
  peakWindowWidth = 200,
  peakfindingSignalLimit = 0.2
  )

  file = joinpath(filepath, outputfilename)
  massAxis = HDF5.h5read(file, "MassAxis")
  avgSpectrum = HDF5.h5read(file, "AvgSpectrum")

  baselinePoints, baselineValues, baselineNoise = BaselineFunctions.calculateBaseline(massAxis,avgSpectrum,baselinePointWidth = 0.8, threshold=0.2)
  baselineNoiseInterpolated = InterpolationFunctions.interpolate(massAxis, baselinePoints, baselineNoise)
  baselineInterpolated = InterpolationFunctions.interpolate(massAxis, baselinePoints, baselineValues)
  baselineCorrectedAvgSpec = avgSpectrum[:,1] - baselineInterpolated;


  peakIndices = PeakshapeFunctions.findPeakIndices(massAxis, avgSpectrum, baselineInterpolated, baselineNoiseInterpolated, noiseThreshold = peakfindingNoiseThresholdValue, signalLimit = peakfindingSignalLimit)
  peakShapesCenterMass, peakShapesY = PeakshapeFunctions.calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices, nbrMassRegions = peakshapeRegions,regionStretch=peakshapeRegionStretch, peakWindowWidth = peakWindowWidth, quantileValue = peakshapeQuantileValue)

  PyPlot.figure()
  PyPlot.title("Baseline Correction")
  PyPlot.semilogy(baselinePoints,baselineValues,".-")
  PyPlot.semilogy(baselinePoints,baselineValues + baselineNoise,".")
  #semilogy(peakMasses, peakValues, "x")
  PyPlot.semilogy(massAxis,avgSpectrum)


  ############ delete h5 data that will be overwritten ###########
  fh = HDF5.h5open(file,"r+")
  for attr in ["AvgBaseline", "MassDepPeakshape", "MassDepPeakshapeCenterMasses", "BaseLines"]
    if haskey(fh, attr)
      HDF5.delete_attribute(fh,attr)
    end
  end
  HDF5.close(fh)


  HDF5.h5write(file, "AvgBaseline", baselineInterpolated)
  HDF5.h5write(file, "MassDepPeakshape", peakShapesY)
  HDF5.h5write(file, "MassDepPeakshapeCenterMasses", peakShapesCenterMass)

  fh = HDF5.h5open(file,"r")
  ds = fh["SumSpecs"]
  spectra = size(ds)[2]
  HDF5.close(fh)
  if (spectra > 0)
    SumSpecs = HDF5.h5read(file, "SumSpecs")
    fill!(SumSpecs,0)
    HDF5.h5write(file, "BaseLines", SumSpecs)

  else
    HDF5.h5write(file, "SumSpecs", avgSpectrum)
    fill!(avgSpectrum,0)
    HDF5.h5write(file, "BaseLines", avgSpectrum)
  end

  println("DONE")
end
