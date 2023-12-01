using TOFTracer2

using Test
using HDF5
# using Suppressor
# using Documenter

@testset "allTests" begin

	@testset "IntegrationBorders" begin
	    @test_throws MethodError TOFTracer2.MasslistFunctions.IntegrationBorders()
	    cmasses = [59.1,63.05]
	    I = TOFTracer2.MasslistFunctions.IntegrationBorders(cmasses)
	    @test round.(I.lowMass,digits=5) == [59.0803,63.0303]
	    @test round.(I.highMass,digits=5) == [59.11970,63.07102]
	    @test I.centerMass == cmasses
	    
	    I = TOFTracer2.MasslistFunctions.IntegrationBorders(cmasses;resolution=10000)
	    @test round.(I.lowMass,digits=5) == [59.09409,63.04409]
	    @test round.(I.highMass,digits=5) == [ 59.10591,63.0563]
	    @test I.centerMass == cmasses
	end

	@testset "processingProject" begin
	    fpfiles = joinpath("..","ExampleFiles")
	    fp = joinpath(fpfiles,"TOFDATA")
	    print("path to files: ", fp, "\n")
	    @test isdir(fp)
	    rf = joinpath(fp,"2016-10-02-19h15m05.h5")
	    @test isfile(rf)
	    masslistfp = joinpath(fpfiles,"MASSLISTS","exampleMassList.csv")
	    @test isfile(masslistfp)
	    cr = [59 391]
	    useAveragesOnly = true
	    masslist = TOFTracer2.MasslistFunctions.loadMasslist(masslistfp)
	    (masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
	    s = (masslistMasses.>0) .& ( masslistMasses.<600)
	    masslistMasses = masslistMasses[s]
	    masslistCompositions = masslistCompositions[s,:]
	    correctMassScaleAndExtractSumSpec(
		fp,
		masslistMasses,
		masslistElements,
		masslistElementsMasses,
		masslistCompositions,
		rf,
		cr,
		filefilterRegexp=r"\.h5$",
		onlyUseAverages = useAveragesOnly,
		plotControlMass = false,
		firstNFiles=0,
		lastNFiles = 0,
		filePrecaching = false,
		openWholeFile = true,
		testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
		testRangeEnd = 137.5,
		)
	    baselineAndPeakshape(
		fp,
		peakshapeRegions=8,
		peakshapeQuantileValue = 0.1,
		peakfindingNoiseThresholdValue = 25
		)
	    mtrx = deconvolute(
		fp,
		calcTransposed = true
		)
	    @test size(masslist[1])[1] == size(mtrx)[1] == size(mtrx)[2]
	    @test length(mtrx.rowval) <= size(mtrx)[1]*size(mtrx)[2]/2
	    outfile = joinpath(fp,"results","_result.hdf5")
	    @test HDF5.ishdf5(outfile)
	    measResult = ResultFileFunctions.loadResults(outfile,useAveragesOnly=true)
	end

end
