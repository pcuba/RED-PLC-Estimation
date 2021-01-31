# Current version: 01/07/2021


using Atom
using Distributions
using StatsBase
using DataFrames
using CSV
using Random
using LinearAlgebra


#--------------------------------------------------------------------------
# Inputs
#--------------------------------------------------------------------------
nLoadParam = "COPFexactZLB"
nParamDGP = "DGP3" # Name of file with in \PriorSpec with parameter values
nSim     = 11000; # Number of simulated observations
nBurnSim = 100; # Number of initial simulations to discard



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
include(readDir * "fStateSpaceForm.jl")
include(readDir * "fAddSteadyStateParamDict.jl")
include(readDir * "fCanonicalForm.jl")
include(readDir * "fStateSpaceSimulationaug.jl")


#--------------------------------------------------------------------------
# Load Solution Options
#--------------------------------------------------------------------------

options_int_use,options_float_use,sPLC = @eval fSolveOptions()

# Create Generic Smolyak Grid (will be stretched according to exact parameters every time)

Smolyak_fixed_grid = fGenSmolyakGrid(options_int_use["Smolyak_d"],options_int_use["Smolyak_mu"],options_float_use["Smolyak_sd"])

# Add R_quantile to the options to be carried

R_quantile = 1.01; # Arbitrary

options_float_use["R_quantile"] = R_quantile;


#--------------------------------------------------------------------------
# Script Starts
#--------------------------------------------------------------------------
nPath = pwd()
rnd_off = 12
Random.seed!(rnd_off+1)
savedir      = "$(pwd())/Data/"

#--------------------------------------------------------------------------
# Load parameters
#--------------------------------------------------------------------------
readDir     = "$(pwd())/Data/"
vParameters = CSV.read(readDir * nLoadParam * "_" * nParamDGP * "_Param.csv" , DataFrame,
              datarow=1 , limit=1)
vParameters = convert(Array,vParameters)'

#--------------------------------------------------------------------------
# Solve model and create state-space representation
#--------------------------------------------------------------------------
sParam            = fCreateParamDict(vParameters)
flag_out, sParam_StateSpace = fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

#--------------------------------------------------------------------------
# Simulate Model
#--------------------------------------------------------------------------
vInitS     = sParam_StateSpace["vInitS"];
pSigmaTE   = sParam_StateSpace["pSigmaTE"];
pA         = sParam_StateSpace["pA"];
pA0        = sParam_StateSpace["pA0"];
D_EpsToEta = sParam_StateSpace["D_EpsToEta"]

pSigmaTE  = sParam_StateSpace["pSigmaTE"];
pSigmaEta = sParam_StateSpace["pSigmaTE"];

eta1BarCoeff = sParam_StateSpace["eta1BarCoeff"];
pPhi1_0      = sParam_StateSpace["Phi01"];
Phi01        = pPhi1_0
pPhi1_1      = sParam_StateSpace["Phi11"];
Phi11        = pPhi1_1
pPhi1_Eta    = sParam_StateSpace["PhiEta1"];
PhiEta1      = pPhi1_Eta
pPhi2_0      = sParam_StateSpace["Phi02"];
Phi02        = pPhi2_0
pPhi2_1      = sParam_StateSpace["Phi12"];
Phi12        = pPhi2_1
pPhi2_Eta    = sParam_StateSpace["PhiEta2"];
PhiEta2      = pPhi2_Eta

# Draw all the shocks
mShockEps = rand(MvNormal(zeros(size(pSigmaTE,1)),pSigmaTE), nSim)';
mShockEta = (mShockEps*D_EpsToEta');

# Create matrices
nNumStates = size(pA,2)
nNumObs    = size(pA,1)
mState     = zeros(nSim,nNumStates)
mSimObs    = zeros(nSim,nNumObs)

for iTime = 1:nSim
    if iTime == 1
        (mState[iTime,:], mSimObs[iTime,:]) = fStateSpaceSimulationaug(eta1BarCoeff,Phi01,Phi11,PhiEta1,Phi02,Phi12,PhiEta2,pA,pA0,vInitS,mShockEta[iTime,:])
    else
        (mState[iTime,:], mSimObs[iTime,:]) = fStateSpaceSimulationaug(eta1BarCoeff,Phi01,Phi11,PhiEta1,Phi02,Phi12,PhiEta2,pA,pA0,mState[iTime-1,:],mShockEta[iTime,:])
    end
end

# Discard first nBurnSim
mState  = mState[nBurnSim + 1:end, :]
mSimObs = mSimObs[nBurnSim + 1:end, :]

# Order observation like input data
mSimObs = [ collect(1:(nSim-nBurnSim)) collect(1:(nSim-nBurnSim)) mSimObs[:,1] mSimObs[:,3] mSimObs[:,4] mSimObs[:,2] mSimObs[:,2] mSimObs[:,2]]
vName   = ["date"; "date"; "ygrowth"; "infl"; "r"; "logcy4"; "logcy4"; "logcy4"]

# Save Output
savename     = nLoadParam * "_" * nParamDGP * "_Data.csv"
CSV.write(savedir * savename,  DataFrame(mSimObs, :auto), header=vName)

mState = [ collect(1:(nSim-nBurnSim)) collect(1:(nSim-nBurnSim)) mState]
vName   = ["date"; "date"; "er"; "g"; "z"; "d"; "y"; "pi"; "c"; "R"; "ylag"]

# Save Output
savename     = nLoadParam * "_" * nParamDGP * "_States.csv"
CSV.write(savedir * savename,  DataFrame(mState, :auto), header=vName)
