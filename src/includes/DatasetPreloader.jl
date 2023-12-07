
module DatasetPreloader
	using HDF5
	using Distributed
	export preload, load, getPreloadFileName

	preloadFileName = ""
	preloadFuture = 0

	function readFromHDF(filename)
		GC.gc()
		try
			return h5read(filename,"/FullSpectra/TofData")
		catch
			return 0
		end
	end

	getPreloadFileName() = preloadFileName

	function preload(filename)
		global preloadFileName = filename
		println("preloading file $filename")
		global preloadFuture = remotecall(readFromHDF,2,filename)
	end

	function load(filename)
		global preloadFileName
		global preloadFuture
		if (preloadFileName != "" && preloadFileName == filename && preloadFuture != 0)
			println("fetching $filename from future")
			ds = fetch(preloadFuture)
			preloadFuture = 0
			preloadFilename = ""
			return ds
		else
			println("reading $filename directly")
			preloadFuture = 0
			preloadFilename = ""
			return readFromHDF(filename)
		end
	end
end
