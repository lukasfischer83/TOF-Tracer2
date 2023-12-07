# Julia lang data processing scripts for Tofwerk based PTR TOF Data



This software provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation.

## Prerequisites

Python3
Julia > 1.6 (ubuntu: add it to your PATH)
Packages: HDF5, PyPlot, DSP, LsqFit (for ApiTOF data analysis)

currently tested with:

Julia Version 1.6.5 and Version 1.8.1
HDF5 v0.16.12 - v0.16.15
DSP v0.7.7
PyPlot v2.11.0
LsqFit v0.13.0

## if you want an editor:
Visual Studio Code (https://code.visualstudio.com/docs/setup/linux )
open VS Code, add julia extension
(see here: https://www.julia-vscode.org/docs/stable/gettingstarted/ )

## Getting Started:

start julia from base directory ("TOFTracer2/")
add packages ( press "]" to enter Pkg mode, run "add PyPlot HDF5 DSP", leave Pkg mode with backspace)
Try the example with include("processingProjects/processingProject-example.jl") -> change all paths to the paths on your PC,


plot the result with plotResults.jl

# Testing the package with
activate the package with 

	] activate .

then 
run tests in package mode ```(TOFTracer2) pkg>``` with 

	test [--coverage]
	


