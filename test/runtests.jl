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

	include("./test_processingWorkflow.jl")

end
