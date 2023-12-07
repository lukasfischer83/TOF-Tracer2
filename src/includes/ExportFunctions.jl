module ExportFunctions
	using ..InterpolationFunctions
	using Dates
	import ..MasslistFunctions

	export exportTracesCSV, exportTracesCSVLossCorr, toMatlabTime, fromMatlabTime

	function exportTracesCSV(saveFolderPath, elementNames, compositions, times, traces; average=0)
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  f = open(joinpath(saveFolderPath,"ptr3compositions.txt"), "w")
	  writedlm(f, hcat(["Mass" "SumFormula"],reshape(elementNames,(1,length(elementNames)))))
	  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions'))
	  close(f)
	  f = open(joinpath(saveFolderPath,"ptr3traces.csv"), "w")
	  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    writedlm(f, hcat(times ,traces))
	  else
	    writedlm(f, hcat(averageSamples(times,average) ,averageSamples(traces,average)))
	  end
	  close(f)
	end

	function exportTracesCSVLossCorr(saveFolderPath, elementNames, compositions, times, traces, lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes; average=0)
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  f = open(joinpath(saveFolderPath,"ptr3compositions.txt"), "w")
	  writedlm(f, hcat(["Mass"	"SumFormula"],reshape(elementNames,(1,length(elementNames))),["LossFactor"	"LossFactorError"	"CorrFactor"	"CorrFactorErr"	"CorrNotes"]))
	  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions', lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes))
	  close(f)
	  f = open(joinpath(saveFolderPath,"ptr3tracesInletLossCorr.csv"), "w")
	  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    writedlm(f, hcat(times , (corrfactor.*traces' )' ))
	  else
	    writedlm(f, hcat(averageSamples(times,average) ,(corrfactor.*(averageSamples(traces,average))' )' ))
	  end
	  close(f)
	end

	################ EXAMPLE from Matlab: 737551.673515479 should be 06-May-2019 16:09:51 #############
	function toMatlabTime(t::Dates.DateTime)
	    timespan = ((t+Dates.Day(1)) - Dates.DateTime(0,1,1,0,0,0))
	    millis::Float64 = Dates.value(timespan)
	    return millis/24/3600/1000
	end

	function fromMatlabTime(timestamp::Number)
	    days=Int(floor(timestamp))
	    millisecondsRemainder = Int(round((timestamp-days)*24*3600*1000))
	    return Dates.DateTime(0,1,1,0,0,0)+Dates.Day(days)+Dates.Millisecond(millisecondsRemainder)-Dates.Day(1)
	end

end
