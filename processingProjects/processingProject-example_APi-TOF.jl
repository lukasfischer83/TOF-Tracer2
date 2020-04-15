#cd("/home/user/github/TOF-Tracer2/")

include("$(pwd())/startupAPiTOF.jl")

fp = "./ExampleFiles/APiTOFDATA/" # All files in this path will be processed
filefilterRegexp = r"\.h5$"
rf = "./ExampleFiles/APiTOFDATA/Data_09_46_20_04mm.h5"  # The mass scale from this file defines the mass scale of all

masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/exampleMasslist.csv")
cr = [37.028406 55.038971 282.118343] # pos
#cr = [62 125] # neg

# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>17) .& ( masslistMasses.<900)
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
    testRangeStart = 54.5, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 55.5,
    recalibInterval = 900,
    resolution = 5500,
    firstNFiles= 1,
    lastNFiles = 0,
    UMRpeaks = true, # create Unit Mass Resolution peaks
    UMRmasses = length(s),
    UMRmassLow = 0.2,
    UMRmassHigh = 0.4,
    umrBinWidth =100
    #binWidth = 10
    )
#
baselineAndPeakshape(
    fp,
    peakshapeRegions=4,
    peakshapeRegionStretch=1.0,
    peakshapeQuantileValue = 0.7,
    peakfindingNoiseThresholdValue = 5,
    peakWindowWidth = 100,
    peakfindingSignalLimit=0.5
    )
    #
#
mtrx = deconvolute(
    fp,
    calcTransposed = false
    )
#
