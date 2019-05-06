# Julia lang data processing scripts for Tofwerk based PTR TOF Data
This software provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation. Update 20/06/2018: New feature: Process data of IONICON ioniAPi-TOF

## Prerequisites
Julia >0.5
Packages: HDF5, PyPlot

## Getting Started:
Try the example in "processingProjects/processingProject-example.jl", plot the result with plotResults.jl
