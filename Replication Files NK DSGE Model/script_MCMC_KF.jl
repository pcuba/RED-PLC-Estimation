
using Atom
using Distributions
using StatsBase
using DataFrames
using CSV
using Random
using LinearAlgebra


clearconsole()

#--------------------------------------------------------------------------
# Include procedures for model solution
#--------------------------------------------------------------------------
readDirSol = "$(pwd())/FunctionsSolution/"

include(readDirSol * "fLinearSolve.jl")
include(readDirSol * "fSGU_gx_hx.jl")
include(readDirSol * "fSGUGetOrder1.jl")
include(readDirSol * "fSolveOptions.jl")

#--------------------------------------------------------------------------
# Include procedures for model estimation
#--------------------------------------------------------------------------
readDir = "$(pwd())/FunctionsEstimation/"
include(readDir * "fCreateParamDict.jl")
include(readDir * "fAddSteadyStateParamDict.jl")
include(readDir * "fCanonicalFormLinear.jl")
include(readDir * "fStateSpaceFormLinear.jl")
# include(readDir * "fStateSpaceSimulationLinear.jl")
include(readDir * "fPriorDensity.jl")
include(readDir * "fPriorDraws.jl")
include(readDir * "run_MCMC_KF.jl")
include(readDir * "fKF.jl")


ME_scale  = 0.0;
nSimRWMH  = 110000;
mhrun     = 1;
nBurnRWMH = 10000;      # number of simulations for the Metropolis-Hasting algorithm
silent    = 1;
nOutput   = 10000;  # print output every nOutput (>1) steps (use = 0 to print each step)
c         = 0.3;     # scaling of proposal covariance matrix for Sim1 data
                       # use 0.002 if constructed from Prior
                       # use 0.3 if constructed from Posterior
nSistResFilter = 0; # = 1 systematic resampling in filter, multinomial otherwise


# Choose Prior

nNamePriorSpec = "Prior1"

# Options for the Proposal Covariance Matrix

nLoadCovVar = 1;        # If 0 uses the file "SummaryDirectDrawsPrior1.csv" under PriorDraws directory,
                        # If 1 uses the saved draws from a previous run's posterior draws

# Choose Dataset

nDataSet = "Sim3" # US or SimX


#--------------------------------------------------------------------------
# Load data set
# date obs ygrowth (%) infl (annual %) R (annual %)
# 100*logcy_raw 100*logcy 100*logcy4
#--------------------------------------------------------------------------

if nDataSet == "US"
   mData    = CSV.read("Data/ACS-RED-USDATA-G.csv",DataFrame)
   nStartDate = 1
   nEndDate   = 140 # 140 is full sample, 96 is preZLB
   mData = convert(Array, mData);
   R_quantile = quantile(mData[nStartDate:nEndDate,5]./400 .+ 1,0.90) # R data quantile

elseif nDataSet == "Sim1"

   mData    = CSV.read("Data/COPFexactZLB_DGP1_Data.csv",DataFrame)
   nStartDate = 401
   nEndDate   = 540
   mData = convert(Array, mData);
   R_quantile = 1.01

elseif nDataSet == "Sim3"

   mData    = CSV.read("Data/COPFexactZLB_DGP3_Data.csv",DataFrame)
   nStartDate = 5901
   nEndDate   = 6040
   mData = convert(Array, mData);
   R_quantile = 1.01

end

nFileStr   = "KF_" * nNamePriorSpec * "_" * nDataSet
nFileStrMH = "KF_" * nNamePriorSpec * "_" * nDataSet * "_Mhrun" *"$(mhrun)"
savedir      = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/"
try mkdir(savedir) catch; end  # Create directory, if it exists skip creation
#--------------------------------------------------------------------------
# Construct proposal covariance matrix
#--------------------------------------------------------------------------

if nLoadCovVar == 0
    readDir       = "$(pwd())/PriorDraws/"
    readName      = "SummaryDirectDrawsPrior1.csv"
    mSummaryDraws = CSV.read(readDir * readName, DataFrame ; datarow=1)
    mSummaryDraws = convert(Array,mSummaryDraws)
    vStdDevTheta  = mSummaryDraws[2,:]
    mSigmaProp    = Diagonal(vStdDevTheta.^2)
else
    readDir       = "$(pwd())/PosteriorDraws/" * nFileStr * "_Mhrun0/"
    readName      = nFileStr * "_Mhrun0_CovVarMatrix.csv"
    mSigmaProp    = CSV.read(readDir * readName, DataFrame; datarow=1)
    mSigmaProp    = convert(Array,mSigmaProp)
end

#--------------------------------------------------------------------------
# Initial Value for MCMC
#--------------------------------------------------------------------------

readDir      = "$(pwd())/PosteriorDraws/"
readName     = nFileStrMH * "_InitialValue.csv"
vParameters  = CSV.read(readDir * readName, DataFrame; datarow=1 , limit=1)
vParameters  = convert(Array,vParameters)'


run_MCMC_KF(ME_scale, nSimRWMH, mhrun, nBurnRWMH,
           silent, nOutput, c, nNamePriorSpec, mData, nStartDate, nEndDate,
           nFileStr, nFileStrMH, vParameters, mSigmaProp)
