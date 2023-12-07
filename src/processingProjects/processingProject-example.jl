using TOFTracer2

fp = joinpath(pwd(),"ExampleFiles","TOFDATA") # All files in this path will be processed
filefilterRegexp = r"\.h5$"
rf = joinpath(fp,"2016-10-02-19h15m05.h5")  # The mass scale from this file defines the mass scale of all
masslist = MasslistFunctions.loadMasslist(joinpath(pwd(),"ExampleFiles","MASSLISTS","exampleMassList.csv"))
cr = [59 391]
goodSignal2Noise = true # PTR3: = true, STOF: = false


# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist


s = (masslistMasses.>0) .& ( masslistMasses.<600)
masslistMasses = masslistMasses[s]
masslistCompositions = masslistCompositions[s,:]

####################### END OF SETTINGS ###############################################################


####################### Processing sequence ###########################################################

correctMassScaleAndExtractSumSpec(
    fp,
    masslistMasses,
    masslistElements,
    masslistElementsMasses,
    masslistCompositions,
    rf,
    cr,
    filefilterRegexp=filefilterRegexp,
    onlyUseAverages = true,
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    filePrecaching = false,
    openWholeFile = true,
    testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 137.5,
    recalibInterval = 60, # longer intervals, e.g. = 300 better for STOF (otherwise mass calib peaks not visible)
    resolution = 7500  # approx. 2000 for STOF
    )

if goodSignal2Noise # PTR3 case
	baselineAndPeakshape(
		fp, 
		peakshapeRegions=8,
		peakshapeRegionStretch=0.5,
		peakshapeQuantileValue = 0.1,
		peakfindingNoiseThresholdValue = 25,
		peakfindingSignalLimit = 0.2
		)
else # STOF case
	baselineAndPeakshape(
		fp, 
		peakshapeRegions=4,
		peakshapeRegionStretch=1,
		peakshapeQuantileValue = 0.2, 		
		peakfindingNoiseThresholdValue = 2,	
		peakfindingSignalLimit = 0.01		
		)
end

mtrx = deconvolute(
    fp,
    calcTransposed = true,
    APITOF = false
    )
