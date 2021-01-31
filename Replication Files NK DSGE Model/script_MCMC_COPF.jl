
# Current version: 01/07/2021

using Atom
using Distributions
using StatsBase
using DataFrames
using CSV
using Random
using LinearAlgebra

#clearconsole()

#--------------------------------------------------------------------------
# Set Options
#--------------------------------------------------------------------------

# Estimation Options

mhrun     = 1;
nM        = 400; # number of particles
ME_scale  = 0.15;
nSimRWMH  = 22000;
nBurnRWMH = 2000;  # number of simulations for the Metropolis-Hasting algorithm
silent    = 1;
nOutput   = 100;  # print output every nOutput (>1) steps (use = 0 to print each step)
c         = 0.1; # scaling of proposal covariance matrix
                 # mhrun1 = 0.1

rnd_off   = 0;    # control random seed
nSistResFilter = 0; # = 1 systematic resampling in filter, multinomial otherwise
pZLB = 0

# Choose Prior

nNamePriorSpec = "Prior1"

# Choose Dataset

nDataSet = "Sim3" # US, US2000 or Sim1

# Choose proposal covariance matrix (if = 0 uses SummaryDirectDrawsPrior1.csv in the prior draws folder
#                                    if = 1 uses the results from a previous run specified below)

nLoadCovVar = 1;

#--------------------------------------------------------------------------
# Include procedures for model solution
#--------------------------------------------------------------------------
readDirSol = "$(pwd())/FunctionsSolution/"

# Solution Options

include(readDirSol * "fSolveOptions.jl")

# Generic SGU Functions

include(readDirSol * "fSGU_gx_hx.jl")

# Generic Smolyak Functions

include(readDirSol * "fSmolyak_Elem_isotrop.jl")
include(readDirSol * "fSmolyak_Grid.jl")
include(readDirSol * "fGenSmolyakGrid.jl")
include(readDirSol * "fScaleSmolyakGrid.jl")

# Generic PLC Functions

include(readDirSol * "fPLCInitialize.jl")

# ZLB Model Specific SGU and PLC Functions

include(readDirSol * "fSGUGetOrder1.jl")
include(readDirSol * "fPLCCoefs.jl")
include(readDirSol * "fPLCSolve.jl")
include(readDirSol * "fPLCRestrictions.jl")
include(readDirSol * "fPLCSystemJacob.jl")
include(readDirSol * "fPLCGetJacobian.jl")
include(readDirSol * "fPLCPreJacobian.jl")

# Solver

include(readDirSol * "fMyLevenbergMarquardt.jl")

# Integration

include(readDirSol * "fMonomials_2_N4.jl")

# Misc

include(readDirSol * "vprint.jl")

#--------------------------------------------------------------------------
# Include procedures for model estimation
#--------------------------------------------------------------------------
readDir = "FunctionsEstimation/"
include(readDir * "fCreateParamDict.jl")
include(readDir * "fAddSteadyStateParamDict.jl")
include(readDir * "fCanonicalForm.jl")
include(readDir * "fStateSpaceForm.jl")
include(readDir * "fStateSpaceSimulationaug.jl")
include(readDir * "fCOPF.jl")
include(readDir * "fCondDist_v5.jl")
include(readDir * "fPriorDraws.jl")
include(readDir * "fPriorDensity.jl")
include(readDir * "run_MCMC_COPF.jl")
include(readDir * "fSystematicResampling.jl")



#--------------------------------------------------------------------------
# Load data set
# date obs ygrowth (%) infl (annual %) R (annual %)
# 100*logcy_raw 100*logcy 100*logcy4
#--------------------------------------------------------------------------

if nDataSet == "US"
   mData    = CSV.read("Data/ACS-RED-USDATA-G.csv",DataFrame)
   nStartDate = 1
   nEndDate   = 140 # 140 is full sample, 96 is pre-2008
   mData = convert(Array, mData);
   R_quantile = quantile(mData[nStartDate:nEndDate,5]./400 .+ 1,0.90) # R data 90th percentile for grid

elseif nDataSet == "US2000"
      mData    = CSV.read("Data/ACS-RED-USDATA-G.csv",DataFrame)
      nStartDate = 65
      nEndDate   = 140 # 140 is full sample, 96 is pre-2008
      mData = convert(Array, mData);
      R_quantile = quantile(mData[nStartDate:nEndDate,5]./400 .+ 1,0.90) # R data 90th percentile for grid

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

#--------------------------------------------------------------------------
# Load Solution Options
#--------------------------------------------------------------------------

options_int_use,options_float_use,sPLC = @eval fSolveOptions()

# Create Generic Smolyak Grid (will be stretched according to exact parameters every time)

Smolyak_fixed_grid = fGenSmolyakGrid(options_int_use["Smolyak_d"],options_int_use["Smolyak_mu"],options_float_use["Smolyak_sd"])

# Add R_quantile to the options to be carried

options_float_use["R_quantile"] = R_quantile;

# Create File Names

nFileStr   = "KF_" * nNamePriorSpec * "_" * nDataSet
nFileStrMH = "COPF_" * nNamePriorSpec * "_" * nDataSet * "_Mhrun" *"$(mhrun)"
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
    readDir       = "$(pwd())/PosteriorDraws/" * nFileStr * "_Mhrun1/"
    readName      = nFileStr * "_Mhrun1_CovVarMatrix.csv"
    mSigmaProp    = CSV.read(readDir * readName, DataFrame ; datarow=1)
    mSigmaProp    = convert(Array,mSigmaProp)
end

#--------------------------------------------------------------------------
# Initial Value for MCMC
#--------------------------------------------------------------------------

readDir      = "$(pwd())/PosteriorDraws/"
readName     = nFileStrMH * "_InitialValue.csv"
vParameters  = CSV.read(readDir * readName, DataFrame , datarow=1 , limit=1)
vParameters  = convert(Array,vParameters)'

run_MCMC_COPF(ME_scale, nM, nSimRWMH, mhrun, nBurnRWMH,
           silent, nOutput, c, nNamePriorSpec, mData, nStartDate, nEndDate,
           nFileStr, nFileStrMH, vParameters, mSigmaProp, rnd_off,
           Smolyak_fixed_grid,options_int_use,options_float_use,sPLC)
