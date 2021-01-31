#--------------------------------------------------------------------------
# Script to:
# Subsample draws from a posterior draw
#--------------------------------------------------------------------------

using Atom
# using Distributions
# using Clustering
# using Optim
# using LsqFit
# using Distances
# import StatsBase.sample
# import StatsBase.pweights
using CSV
using DataFrames
#--------------------------------------------------------------------------
# Load the draws
#   each row is a draw for the parameter vector
#   each column is a vector of draws for a specific parameter
#--------------------------------------------------------------------------

nSpacing = 100;
nFilter = "COPF"; # Choose COPFexactZLB or KF or BSPF
nPrior  = "1";
nMhrun  = "1";
nDataSet = "Sim3" # US Or Sim1


nFileStrMH = nFilter * "_Prior" * nPrior * "_" * nDataSet * "_Mhrun" * nMhrun

readDir           = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/"
readNameDraws     = nFileStrMH * "_loop_Draws.csv"

mDraws            = CSV.read(readDir * readNameDraws, DataFrame;header=0)
mDraws = convert(Array,mDraws)

nSizeDraw = size(mDraws)[1]

nSizeSubsample = Int(floor(nSizeDraw/nSpacing))

mSubSample = zeros(nSizeSubsample,size(mDraws)[2])
# Create Subsample
for iSubDraw = 1:nSizeSubsample
    mSubSample[iSubDraw,:] = mDraws[iSubDraw*nSpacing,:]
end

#--------------------------------------------------------------------------
# Save the subsample
#--------------------------------------------------------------------------
saveDir      = readDir
saveName     = nFileStrMH * "_SubsampleSpacing" * "$(nSpacing)" * ".csv"
CSV.write(saveDir * saveName,  DataFrame(mSubSample), writeheader=false)
