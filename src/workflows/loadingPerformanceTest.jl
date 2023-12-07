usePrecaching = true

using Distributed
import Statistics
import .TOFFunctions

function averageAllSpectra(
  filepath;
  filefilterRegexp=r"\.h5$",
  filePrecaching = true, # Precache whole next file. Only use when you have enough RAM and no SWAP has to be used.
  openWholeFile = true
  )

    allFiles = readdir(filepath)
    files = filter(s->occursin(filefilterRegexp, s), allFiles)
    nFiles = size(files,1)

    validFiles, timeSortIndices = TOFFunctions.validateHDF5Files(filepath, files)

    if !issorted(timeSortIndices)
    println("Reordering Files according to start acquisition time!")
    else
    println("File order seems fine, continuing...")
    end
    files = validFiles[timeSortIndices]
    nFiles = size(files,1)

    averagedSpectra = Array{Any,1}()

    @time for j=1:5
        println("working on $filepath")
        totalPath = joinpath(filepath, files[j])
        if j < nFiles && filePrecaching
            totalPrecachePath = joinpath(filepath, files[j+1])
        else
            totalPrecachePath = ""
        end

        totalSubSpectra = TOFFunctions.getSubSpectraCount(totalPath)
        sumspec = TOFFunctions.getSubSpectrumFromFile(totalPath,1, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
        for i=2:totalSubSpectra
            sumspec = sumspec .+ TOFFunctions.getSubSpectrumFromFile(totalPath,i, preloadFile=totalPrecachePath, openWholeFile = openWholeFile)
        end
        push!(averagedSpectra, sumspec./totalSubSpectra)
    end
    return averagedSpectra
end
