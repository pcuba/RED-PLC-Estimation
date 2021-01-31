#--------------------------------------------------------------------------
# Script to:
# Generate draws from prior using RWMH
#--------------------------------------------------------------------------

using Atom
using Distributions
import StatsBase.sample
import StatsBase.pweights
using CSV
using DataFrames

clearconsole()
#--------------------------------------------------------------------------
# Load the draws
#   each row is a draw for the parameter vector
#   each column is a vector of draws for a specific parameter
#--------------------------------------------------------------------------

nFilter = "KF"; # Choose COPFexactZLB or COPF or KF or BSPF
nPrior  = "1";
nMhrun  = "1";
nDataSet = "Sim3" # US Or Sim1

nFileStrMH = nFilter * "_Prior" * nPrior * "_" * nDataSet * "_Mhrun" * nMhrun

readDir           = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/"
readNameDraws     = nFileStrMH * "_Draws.csv"
readNameLopPost   = nFileStrMH * "_LogPosterior.csv"

mDraws            = CSV.read(readDir * readNameDraws, DataFrame;header=0)
mDraws = convert(Array,mDraws)

vLogPost          = CSV.read(readDir * readNameLopPost, DataFrame;header=0)
vLogPost = convert(Array,vLogPost)

# Get Max Log Posterior index and maximum posterior density estimator
nMaxLogPOstId = findmax(vLogPost)[2][1]

vMPDEstimator = mDraws[nMaxLogPOstId,:]

# Get CovVar Matrix from draws
mCovMat = cov(mDraws)
saveDir      = readDir
saveName     = nFileStrMH * "_CovVarMatrix.csv"
CSV.write(saveDir * saveName,  DataFrame(mCovMat), header=false)

#--------------------------------------------------------------------------
# Compute means and variances
#--------------------------------------------------------------------------
vMeanTheta   = mean(mDraws; dims = 1)
vStdDevTheta = std(mDraws; dims = 1)

#--------------------------------------------------------------------------
# HPD intervals: take as input posterior draws
#--------------------------------------------------------------------------
include("$(pwd())/FunctionsEstimation/" * "fHPDinterval.jl")
nHPDprob = 0.90;
mHPDinterval = fHPDinterval(mDraws,nHPDprob)

#--------------------------------------------------------------------------
# Save the results
#--------------------------------------------------------------------------
mOutput = [vMeanTheta; vStdDevTheta; mHPDinterval; vMPDEstimator']
saveDir      = readDir
saveName     = nFileStrMH * "_Summary.csv"
CSV.write(saveDir * saveName,  DataFrame(mOutput, :auto), header=false)
