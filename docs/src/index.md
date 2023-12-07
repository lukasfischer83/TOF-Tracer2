# TOFTracer2 Documentation

The TOFTracer2 project is an open-source tool for time-of-flieght mass spectrometry data postprocessing.

It provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation.

```@contents
```

## Processing Workflow


### PeakshapeFunctions
```@docs
TOFTracer2.PeakshapeFunctions
```

```@docs
TOFTracer2.PeakshapeFunctions.findPeakIndices(massAxis, avgSpectrum, baseline, baselineNoise; kwargs...)
TOFTracer2.PeakshapeFunctions.calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices; kwargs...)
TOFTracer2.PeakshapeFunctions.getLocalPeakshape(mass, peakShapesCenterMass, peakShapesY)
```

## PlotFunctions
```@docs
PlotFunctions
```
### massDefectPlot()
```@docs
TOFTracer2.PlotFunctions.massDefectPlot(masses, compositions, concentrations, colors, plotTitle, colorCodeTitle; kwargs...)
```
### plotTracesFromHDF5()
```@docs
TOFTracer2.PlotFunctions.plotTracesFromHDF5(file, massesToPlot; kwargs...)
```

## Index

```@index
```
