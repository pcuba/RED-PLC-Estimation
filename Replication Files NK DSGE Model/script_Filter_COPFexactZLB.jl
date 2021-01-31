
# Current version: 01/07/2021

using Atom
using Distributions
using StatsBase
using DataFrames
using CSV
using Random
using LinearAlgebra

#--------------------------------------------------------------------------
# Set Options
#--------------------------------------------------------------------------

nM        = 250; # number of particles
ME_scale  = 0.001;
mhrun     = 1;
silent    = 1;
nSistResFilter = 0; # = 1 systematic resampling in filter, multinomial otherwise
pZLB      = 0;
rnd_off   = 0;

# What parameters to use? (MAP or Mean)

ind_param = "MAP";

# Choose Prior

nNamePriorSpec = "Prior1"

# Choose Dataset

nDataSet = "US" # US, US2000 or Sim1

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
include(readDir * "fCOPFaugExactZLB.jl")
include(readDir * "fCondDistExactZLB_v5.jl")
include(readDir * "fPriorDraws.jl")
include(readDir * "fPriorDensity.jl")
include(readDir * "run_filter_COPFexactZLB.jl")
include(readDir * "fSystematicResampling.jl")


#--------------------------------------------------------------------------
# Load data set
# date obs ygrowth (%) infl (annual %) R (annual %)
# 100*logcy_raw 100*logcy 100*logcy4
#--------------------------------------------------------------------------
nDataSet = "US" # US or Sim1

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

else # if nDataSet == "Sim1"
   mData    = CSV.read("Data/COPFexactZLB_DGP1_Data.csv",DataFrame)
   nStartDate = 401
   nEndDate   = 540
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

nFileStrMH   = "COPFexactZLB_" * nNamePriorSpec * "_" * nDataSet * "_Mhrun" * "$(mhrun)"
nFileStrFilt = "COPFexactZLB_Filtered_" * ind_param * "_" * nNamePriorSpec * "_" * nDataSet * "_Mhrun" * "$(mhrun)"
savedir = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/Filtered/"
try mkdir(savedir) catch; end  # Create directory, if it exists skip creation

run_filter_COPFexactZLB(ME_scale, nM, mhrun, silent, mData, nFileStrMH, nFileStrFilt,
                        rnd_off,Smolyak_fixed_grid,options_int_use,options_float_use,sPLC,ind_param)
