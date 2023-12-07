### libs
using PyCall
using HDF5
using PyPlot, Colors
using Statistics
###
# path
file = joinpath("path","to","file.h5")
filesavepath = joinpath("path","to","resultfile.h5")
# loading raw data
avgSpectrum = HDF5.h5read(file, "SPECdata/AverageSpec")
timebins = Array{Float64}(undef,length(avgSpectrum),1)
timebins = range(1,stop=length(avgSpectrum),length=length(avgSpectrum))
rawMassCalParams = HDF5.h5read(file, "CALdata/Spectrum")[1:2,2]
rawMassAxis = Array{Float64}(undef,length(avgSpectrum),1)
rawMassAxis = ((timebins .- rawMassCalParams[2])./rawMassCalParams[1]).^2
### plotting raw data
fig = figure(figsize=(12,8));
ax1 = fig.add_subplot(2,1,1);
semilogy(rawMassAxis,avgSpectrum)
xlabel("m/z [Th]")
ylabel("intensity [a.u.]")
title("To get the timebin interval of the desired peaks for mass axis calibration:\n click left and right of each peak in the lower panel.")
ax2 = fig.add_subplot(2,1,2);
semilogy(timebins,avgSpectrum)
xlabel("time bin [#]")
ylabel("intensity [a.u.]")
#savefig("$filesavepath plot_file.png", papertype="A4",orientation="landscape")
### processing data
using LsqFit
@. model(x, p) = p[1]*exp((x-p[2])^2/p[3]) # Gaussian distribution
#@. model(x, p) = p[1]/(p[1]^2 + (x-p[2])^2 ) # Lorentzian distribution, p[1]>0, -inf<p[2]<inf
massCalCompounds = [37.028406, 55.038971, 282.118343]

# select TIMEBIN intervals for mass axis calibration
# -> click left and right of each peak to get the interval
nSteps = length(massCalCompounds) # number of peaks for mass axis calibration
gdata = Array{Float32}(undef,(2*nSteps),1)
for i =1:2*nSteps
    gin = ginput()
    gdata[i] = gin[1][1]
    ax2.plot([gdata[i], gdata[i]], [minimum(avgSpectrum), maximum(avgSpectrum)],"r-.")
end
# find interval indices
maxPeakVal = Any[]
IdxT = Array{Float64}(undef,length(avgSpectrum),1)
for l=1:1
    for i = 1:2:2*nSteps
        specInt = Any[]
        idxSpecInt = Any[]
        for m = 1:length(timebins)
            if gdata[i] < timebins[m] && timebins[m] < gdata[i+1]
                push!(specInt,avgSpectrum[m])
                push!(idxSpecInt,timebins[m])
            end
        end
        p0 = [maximum(specInt), idxSpecInt[findall(x->x==maximum(specInt),specInt)[1]], -maximum(specInt)*1000] # p0 for Guassian distribution
        #p0 = [maximum(specInt)*7, mean(idxSpecInt)] # p0 for Lorentzian distribution
        #lb = [maximum(specInt)-0.1*maximum(specInt), mean(idxSpecInt)-0.01*mean(idxSpecInt),-maximum(specInt)*10000]
        #ub = [maximum(specInt)+0.1*maximum(specInt), mean(idxSpecInt)+0.01*mean(idxSpecInt),Inf]
        fit = LsqFit.curve_fit(model, idxSpecInt, specInt, p0)#, lower=lb, upper=ub )
        push!(maxPeakVal, fit.param[2])
        figure(figsize=(10,7));
        plot(idxSpecInt,specInt)
        plot(idxSpecInt,model(idxSpecInt,fit.param) )
    end
    #return specInt, idxSpecInt, fit, maxPeakVal
end
# 3-params
@. model_tb2m(x, p) = ((x-p[2])/p[1])^(1/p[3])
p0_tb2m = [15000.0, -60000.0, 0.5] # [minimum(maxPeakVal)/1.4,maximum(maxPeakVal)]
fit_tb2m = LsqFit.curve_fit(model_tb2m, maxPeakVal, massCalCompounds, p0_tb2m ) # , p0_bounds, lower=lb, upper=ub)
tb2m = ((timebins.-fit_tb2m.param[2]) ./fit_tb2m.param[1]).^(1/fit_tb2m.param[3])
### plot result
figure(figsize=(12,8));
semilogy(rawMassAxis,avgSpectrum, label="old mass calibration")
semilogy(tb2m,avgSpectrum, label="new mass calibration")
xlabel("m/z [Th]")
ylabel("intensity [a.u.]")
legend()

### save data to APi HDF5 file
processedfilepath = file
fid = HDF5.h5open(processedfilepath, "cw")
HDF5.h5write(processedfilepath, "fitParams", fit_tb2m.param)
HDF5.close(fid)
### test if writing was succesful
fitiPara = HDF5.h5read(processedfilepath, "/fitParams")
