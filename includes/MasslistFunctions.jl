module MasslistFunctions
	using HDF5
	using DelimitedFiles

	export IntegrationBorders, masslistPos, createMassList, loadMasslist, createCompound, massFromComposition, massFromCompositionArray, massFromCompositionArrayList, isotopesFromComposition, isotopesFromCompositionArray, sumFormulaStringFromCompositionArray, sumFormulaStringListFromCompositionArrayList, filterMassListByContribution1, filterMassListByContribution2, findClosestMassIndex, inCompositions, findInCompositions

	struct IntegrationBorders
	    centerMass::Array{Float64}
	    lowMass::Array{Float64}
	    highMass::Array{Float64}

	    function IntegrationBorders(masslist::Array{Float64,1}, resolution)
		borders = new(similar(masslist), similar(masslist),similar(masslist))
		if length(masslist) < 2
		    # return just one
		    borders.lowMass[1] = masslist[1]-masslist[1]/resolution
		    borders.highMass[1] = masslist[1]+masslist[1]/resolution
		    borders.centerMass[1] = masslist[1]
		    return borders
		end
		#set center masses
		#1st and last don't have two neighbors
		borders.lowMass[1] = masslist[1]-(masslist[1]/resolution)
		borders.highMass[end] = masslist[end]+(masslist[end]/resolution)
		borders.centerMass[1] = masslist[1]
		borders.centerMass[end] = masslist[end]

		for i=1:length(masslist)-1 # working on gaps between masses, not on masses
		    borders.centerMass[i] = masslist[i]
		    if (2*masslist[i]/resolution) < (masslist[i+1]-masslist[i]) # distance to next mass is bigger than 2*resolution --> use resolution border
		        borders.highMass[i] = masslist[i]+(masslist[i]/resolution)
		        borders.lowMass[i+1] = masslist[i+1]-(masslist[i]/resolution)
		    else # use half the distance as border
		        center = (masslist[i] + masslist[i+1])/2
		        borders.highMass[i] = center
		        borders.lowMass[i+1] = center
		    end
		end
		return borders
	    end
	    function IntegrationBorders(masslist::Array{Float64,1})
		return IntegrationBorders(masslist, 3000.0)
	    end
	end

	massC=12
	massC13=13.00335
	massH=1.00783
	massHplus=1.007276
	massN=14.00307
	massO=15.99492
	massO18=17.99916
	massS=31.97207
	nC13=0
	nO18=0

	abundanceC13 = 0.010816
	abundanceO18 = 0.002005

	function masslistPos(x)
	  if (x>0)
	    return x
	  else
	    return 0
	  end
	end

	function inCompositions(composition, compositionList)
	  for i=1:size(compositionList,2)
	    if composition == compositionList[:,i]
	      return true
	    end
	  end
	  return false
	end

	masslistElements = ["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"]
	masslistElementMasses = [massC, massC13, massH, massHplus, massN, massO, massO18, massS]
	function createMassList(; C=0:0, O=0:0, N=0:0, S=0:0, nHplus=1, allowRadicals=false)
	  masses = [37.02]
	  masslistCompositions = [[0 0 4 1 0 2 0 0]]
	  for nC in C
	    for nN in N
	      for nS in S
		for nO in O[O.<=nC+1]
		  if allowRadicals
		    possibleH = masslistPos(round(nC/2)*2-4):1:2*nC+2+nN*3
		  else
		    possibleH = masslistPos(round(nC/2)*2-4):2:2*nC+2+nN*3
		  end
		  if (isodd(nN))
		    possibleH = possibleH + 1
		  end
		  for nH in possibleH
		    append!(masses,nS*massS+nC*massC+nN*massN+nO*massO+nH*massH+nHplus*massHplus)
		    push!(masslistCompositions, [nC nC13 nH nHplus nN nO nO18 nS])
		  end
		end
	      end
	    end
	  end
	  sortIndices = sortperm(masses)
	  return masses[sortIndices],masslistElements,masslistElementMasses, masslistCompositions[sortIndices]
	end

	function loadMasslist(filename)
	    fileIsValidResultfile = (endswith(filename, ".hdf5") || endswith(filename, ".h5"))
	    if fileIsValidResultfile
		println("Trying to load masslist from result hdf5 file ($filename)!")
		masses =  HDF5.h5read(filename, "MassList")
		masslistElementsLoaded = HDF5.h5read(filename,"ElementNames")
		compRaw = HDF5.h5read(filename, "ElementalCompositions")
		masslistCompositions = []
		for i=1:size(compRaw,2)
		    push!(masslistCompositions,compRaw[:,i])
		end
		return masses,masslistElementsLoaded, masslistElementMasses, masslistCompositions
	    else
	      println("Trying to load masslist from csv file ($filename)!")
	      list = readdlm(filename, '\t', skipstart=7)
	      masses = Array{Float64,1}()
	      masslistCompositions = []
	      print("Unidentified Peaks found in masslist $filename: ")
	      for i=1:size(list,1)
		if sum(list[i,1:8]) > 0
		  push!(masses,massFromCompositionArray(list[i,1:8]))
		  push!(masslistCompositions, list[i,1:8])
		else
		    if list[i,9] > 0
		    print("$(list[i,9]) ")
		    push!(masses,list[i,9])
		    push!(masslistCompositions, list[i,1:8])
		    end
		end
	      end
	      sortIndices = sortperm(masses)
	      return masses[sortIndices],masslistElements, masslistElementMasses, masslistCompositions[sortIndices]
	  end
	end

	function createCompound(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  mass = C*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + C13*massC13 + O18*massO18
	  composition = [C, C13, H, Hplus, N, O, O18, S]
	  return mass, masslistElements, masslistElementMasses, composition
	end

	function massFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  mass = C*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + C13*massC13 + O18*massO18
	  #composition = [C C13 H Hplus N O O18 S]
	  return mass
	end

	function isotopesFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  masses = Array{Float64,1}()
	  compositions = []
	  abundances = Array{Float64,1}()
	  for carbIsotopeNr = 0:minimum([C,2])
	    mass = (C-carbIsotopeNr)*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + (C13+carbIsotopeNr)*massC13 + O18*massO18
	    abundance = (abundanceC13^carbIsotopeNr)*binomial(Int64(C),Int64(carbIsotopeNr))
	    composition = [C-carbIsotopeNr C13+carbIsotopeNr H Hplus N O O18 S]
	    push!(masses,mass)
	    push!(compositions, composition)
	    push!(abundances, abundance)
	  end
	  if (length(masses) == 1 && masses[1] == 0)
	    println("Strange composition found, maybe unidentified compound: $(compositions[1])")
	  end
	  return masses, masslistElements, compositions, abundances
	end

	function isotopesFromCompositionArray(composition)
	  #[C C13 H Hplus N O O18 S]
	  return isotopesFromComposition(C=composition[1],
	  C13=composition[2],
	  H=composition[3],
	  Hplus=composition[4],
	  N=composition[5],
	  O=composition[6],
	  O18=composition[7],
	  S=composition[8])
	end

	function massFromCompositionArray(composition)
	  #mass = composition[1]*massC + composition[2]*massC13 +composition[3]*massH + composition[4]*massHplus
	  #+ composition[5]*massN + composition[6]*massO + composition[7]*massO18 + composition[8]*massS
	  mass = sum(masslistElementMasses .* composition)
	  #composition = [C C13 H Hplus N O O18 S]
	  return mass
	end

	function massFromCompositionArrayList(compositions)
	  ret = Array{Float32,1}()
	  for i=1:size(compositions,2)
	    push!(ret,massFromCompositionArray(compositions[:,i]))
	  end
	  return ret
	end


	function sumFormulaStringFromCompositionArray(composition)
	  if (composition[1]>1)
	    CString = "C$(composition[1])"
	  elseif composition[1] == 1
	    CString = "C"
	  else
	    CString = ""
	  end

	  if (composition[3]>1)
	    HString = "H$(composition[3])"
	  elseif composition[3] == 1
	    HString = "H"
	  else
	    HString = ""
	  end

	  if (composition[6]>1)
	    OString = "O$(composition[6])"
	  elseif composition[6] == 1
	    OString = "O"
	  else
	    OString = ""
	  end

	  if (composition[5]>1)
	    NString = "N$(composition[5])"
	  elseif composition[5] == 1
	    NString = "N"
	  else
	    NString = ""
	  end

	  if (composition[8]>1)
	    SString = "S$(composition[8])"
	  elseif composition[8] == 1
	    SString = "S"
	  else
	    SString = ""
	  end

	  return "$CString$HString$OString$NString$SString"
	end

	function sumFormulaStringListFromCompositionArrayList(compositions; showMass = false)
	  ret = Array{String,1}()
	  for i=1:size(compositions,2)
	    if showMass
	      push!(ret,"$(round(massFromCompositionArray(compositions[:,i]),4)) - $(sumFormulaStringFromCompositionArray(compositions[:,i]))")
	    else
	      push!(ret,sumFormulaStringFromCompositionArray(compositions[:,i]))
	    end
	  end
	  return ret
	end

	function filterMassListByContribution1(masses, countrates, relThreshold)
	  selected = Array{Bool}(undef,length(masses))
	  for i=1:ceil(maximum(masses))
	    sel = (masses.>i-0.3) & (masses.<i+0.7)
	    quant = maximum(countrates[sel])*threshold
	    selected[sel] = countrates[sel].>quant
	  end
	  return selected
	end

	function filterMassListByContribution2(masses, countrates, resolution, relThreshold)
	  selected = Array{Bool}(undef,length(masses))
	  fill!(selected,true)

	  for i=2:length(masses)-1
	    if abs(masses[i]-masses[i+1]) < masses[i]/4*resolution # if neighbor is within resolution
	      th =  relThreshold * masses[i]/(resolution*abs(masses[i]-masses[i+1]))
	      print("\nMass $(masses[i]) close to neighbor, th=$th  ")
	      if countrates[i] < (countrates[i+1] * th) # deselect if contribution is lower than relThreshold*neighbor
		print("REMOVED")
		selected[i] = false
	      end
	    end
	    if abs(masses[i]-masses[i-1]) < masses[i]/4*resolution # if neighbor is within resolution
	      th = relThreshold * masses[i]/(resolution*abs(masses[i]-masses[i-1]))
	      print("\nMass $(masses[i]) close to neighbor, th=$th  ")
	      if countrates[i] < (countrates[i-1] * th) # deselect if contribution is lower than relThreshold*neighbor
		print("REMOVED")
		selected[i] = false
	      end
	    end
	    if countrates[i] < 0
	      selected[i] = false
	    end
	  end
	  return selected
	end


	function OScFromComposition(composition)
	  return (2*composition[6]-composition[3])/composition[1] # OSc = 2*O/C - (H/C)
	end
	function OScFromCompositionArray(compositions)
	  ret = Array{Float64,1}()
	  for i=1:size(compositions,2)
	      push!(ret,OScFromComposition(compositions[:,i]))
	  end
	end

	function findClosestMassIndex(mass, masses)
	    closestIndex::Int64 = 1
	    for i=1:length(masses)
		if abs(masses[i]-mass)<abs(masses[closestIndex]-mass)
		    closestIndex = i
		end
	    end
	    return closestIndex
	end

	function findInCompositions(compositions, composition)
	    for i=1:size(compositions,2)
		if(compositions[:,i] == composition)
		    return i
		end
	    end
	    return 0
	end

end

