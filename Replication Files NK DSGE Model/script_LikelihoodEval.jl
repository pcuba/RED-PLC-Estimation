
# Current version: 01/08/2021

using Atom
using Distributions
using StatsBase
using DataFrames
using CSV
using Random
using LinearAlgebra

# clearconsole()

#--------------------------------------------------------------------------
# Set Options
#--------------------------------------------------------------------------

nSimEval   = 100;
nParamDraw = 1;
ME_scale   = 0.1;

nFilter1   = "BSPF" # KF, COPFexactZLB or BSPF
nM_1       = 1000; # number of particles

nFilter2   = "COPF" # KF, COPFexactZLB or BSPF
nM_2       = 400; # number of particles

nFilter3   = "BSPF" # KF, COPFexactZLB or BSPF
nM_3       = 10000; # number of particles


# Create the vectors for the loop

# vFilter    = [ nFilter1 nFilter2  ]
# vM         = [ nM_1 nM_2]

vFilter    = [ nFilter3 ]
vM         = [  nM_3 ]


nSistResFilter = 0; # = 1 systematic resampling in filter, multinomial otherwise
pZLB = 0;
rnd_off = 0;

# Choose Dataset

nDataSet = "Sim3" # US, US2000 or Sim1


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
include(readDir * "fCanonicalFormLinear.jl")
include(readDir * "fStateSpaceForm.jl")
include(readDir * "fStateSpaceFormLinear.jl")
include(readDir * "fStateSpaceSimulationaug.jl")
include(readDir * "fBSPF.jl")
include(readDir * "fCOPFaugExactZLB.jl")
include(readDir * "fCOPF.jl")
include(readDir * "fKF.jl")
include(readDir * "fCondDistExactZLB_v5.jl")
include(readDir * "fCondDist_v5.jl")
include(readDir * "fPriorDraws.jl")
include(readDir * "fPriorDensity.jl")
include(readDir * "run_LikelihoodEval.jl")
include(readDir * "fSystematicResampling.jl")


#--------------------------------------------------------------------------
# Load parameter draws to evaluate the filters
#--------------------------------------------------------------------------
#nNamePriorSpec = "Prior1"

mDraw = CSV.read("LikelihoodEvaluations/COPF_Prior1_Sim3_Mhrun1_SubsampleSpacing100.csv", DataFrame; header=0)
mDraw = convert(Array,mDraw)
mDraw = mDraw[1:nParamDraw,:]

#--------------------------------------------------------------------------
# Load data set
# date obs ygrowth (%) infl (annual %) R (annual %)
# 100*logcy_raw 100*logcy 100*logcy4
#--------------------------------------------------------------------------

if nDataSet == "Sim1"
   mData    = CSV.read("Data/COPFexactZLB_DGP1_Data.csv",DataFrame)
   mState   = CSV.read("Data/COPFexactZLB_DGP1_States.csv",DataFrame)
   nStartDate = 401
   nEndDate   = 540 # 540
   mData      = convert(Array, mData);
   mState     = convert(Array, mState);
   vInitS     = mState[(nStartDate-1), 3:end]
   R_quantile = 1.01
   vTrueParam = convert(Array,CSV.read("Data/COPFexactZLB_DGP1_Param.csv", DataFrame; header=0))
   mDraw[1,:] = vTrueParam
   mDraw[:,21:end] = ones(nParamDraw,1) .* vInitS'

elseif nDataSet == "Sim3"
      mData    = CSV.read("Data/COPFexactZLB_DGP3_Data.csv",DataFrame)
      mState   = CSV.read("Data/COPFexactZLB_DGP3_States.csv",DataFrame)
      nStartDate = 5901
      nEndDate   = 6040
      mData      = convert(Array, mData);
      mState     = convert(Array, mState);
      vInitS     = mState[(nStartDate-1), 3:end]
      R_quantile = 1.01
      vTrueParam = convert(Array,CSV.read("Data/COPFexactZLB_DGP3_Param.csv", DataFrame; header=0))
      mDraw[1,:] = vTrueParam
      mDraw[:,21:end] = ones(nParamDraw,1) .* vInitS'

elseif  nDataSet == "US2007Q4"
   mData    = CSV.read("Data/ACS-RED-USDATA-G.csv",DataFrame)
   nStartDate = 1
   nEndDate   = 96 # 140 is full sample, 96 is preZLB
   mData = convert(Array, mData);
   R_quantile = quantile(mData[nStartDate:nEndDate,5]./400 .+ 1,0.90) # R data quantile

else # if  nDataSet == "US"
  mData    = CSV.read("Data/ACS-RED-USDATA-G.csv",DataFrame)
  nStartDate = 1
  nEndDate   = 140 # 140 or 95
  mData = convert(Array, mData);
  R_quantile = quantile(mData[nStartDate:nEndDate,5]./400 .+ 1,0.90) # R data quantile

end


#--------------------------------------------------------------------------
# Load Solution Options
#--------------------------------------------------------------------------

options_int_use,options_float_use,sPLC = @eval fSolveOptions()

# Create Generic Smolyak Grid (will be stretched according to exact parameters every time)

Smolyak_fixed_grid = fGenSmolyakGrid(options_int_use["Smolyak_d"],options_int_use["Smolyak_mu"],options_float_use["Smolyak_sd"])

# Add R_quantile to the options to be carried

options_float_use["R_quantile"] = R_quantile;

run_LikelihoodEval(nSimEval, ME_scale, vFilter, vM, nDataSet,
                mData, nStartDate, nEndDate, mDraw, nParamDraw,
                nSistResFilter,
                pZLB, rnd_off,
                Smolyak_fixed_grid,options_int_use,options_float_use,sPLC)
