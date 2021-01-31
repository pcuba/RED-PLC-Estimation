#--------------------------------------------------------------------------
# Script to:
# Generate draws from prior by direct sampling
#--------------------------------------------------------------------------

# Current version: 01/05/2021

using Atom
using Distributions
using Random
using CSV
using DataFrames

clearconsole()

#--------------------------------------------------------------------------
# Select prior specification
#--------------------------------------------------------------------------
nNamePriorSpec = "Prior1"
include("PriorSpec/" * nNamePriorSpec * ".jl")
# prior specification is in mPriorSpec

#--------------------------------------------------------------------------
# Include relevant procedures
#--------------------------------------------------------------------------
include("FunctionsEstimation/fPriorDraws.jl")

#--------------------------------------------------------------------------
# Direct sampling from the prior
#--------------------------------------------------------------------------
nNumDraw = 10000;  # Select number of draws

nRNG = Random.seed!(123)
mPriorDraws = fPriorDraws(nNumDraw,mPriorSpec)

saveDir      = "$(pwd())/PriorDraws/"
saveName     = "DirectDraws" * nNamePriorSpec * ".csv"
CSV.write(saveDir * saveName,  DataFrame(mPriorDraws, :auto), header=false)

mDraws        = convert(Array,mPriorDraws)

# Get CovVar Matrix from draws
mCovMat      = cov(mDraws)
saveName     = "CovVarMatrixDirectDraws" * nNamePriorSpec * ".csv"
CSV.write(saveDir * saveName,  DataFrame(mCovMat, :auto), header=false)

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
mOutput = [vMeanTheta; vStdDevTheta; mHPDinterval]
saveName     = "SummaryDirectDraws" * nNamePriorSpec * ".csv"
CSV.write(saveDir * saveName,  DataFrame(mOutput, :auto), header=false)
