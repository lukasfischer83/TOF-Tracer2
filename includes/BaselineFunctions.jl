module BaselineFunctions

	import Statistics

	export calculateBaseline

	function calculateBaseline(massAxis, avgSpectrum; baselinePointWidth = 0.3, threshold=0.5)

	  baselinePoints = ceil(massAxis[1])+0.6:baselinePointWidth:floor(massAxis[end]-1)-0.4
	  baselineValues = similar(baselinePoints)
	  baselineNoise = similar(baselinePoints)

	  for i=1:length(baselinePoints)
	    startIdx = searchsortedfirst(massAxis,baselinePoints[i]-baselinePointWidth)
	    endIdx = searchsortedfirst(massAxis,baselinePoints[i]+baselinePointWidth)
	    subSet = view(avgSpectrum, startIdx:endIdx)
	    threshold = Statistics.quantile(subSet,0.2)
	    baselineSamples = Array{Float64,1}()
	    #fill!(baselineSamples,0.0)
	    for j=1:length(subSet)
	      if subSet[j] <= threshold
		push!(baselineSamples,subSet[j])
	      end
	    end
	    baselineValues[i] = Statistics.mean(baselineSamples)
	    if length(baselineSamples) > 3
	      baselineNoise[i] = Statistics.std(baselineSamples)
	    else
	      baselineNoise[i] = baselineValues[i]
	    end
	    if (baselineValues[i] > 10000 || baselineValues == 0)
	      println("Strange Baseline at $(baselinePoints[i]):\n$baselineSamples")
	    end

	  end
	  return baselinePoints,baselineValues, baselineNoise
	end
end


