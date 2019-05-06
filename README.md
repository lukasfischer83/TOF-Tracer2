# Julia lang data processing scripts for Tofwerk based PTR TOF Data
This software provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation.

## Prerequisites
Julia > 1.1
Packages: HDF5, PyPlot, DSP

## Getting Started:
start julia from base directory ("TOF-Tracer2/")
add packages ( press "]" to enter Pkg mode, run "add PyPlot HDF5 DSP", leave Pkg mode with backspace)
Try the example with include("processingProjects/processingProject-example.jl"), plot the result with plotResults.jl
