# -----------------------------------------------------------------------
# Simulation in log deviations from the steady state
# vStatePrev = ()
# -----------------------------------------------------------------------
function fStateSpaceSimulationaug(eta1BarCoeff::Array{Float64,2},Phi01::Array{Float64,2},Phi11::Array{Float64,2},PhiEta1::Array{Float64,2},Phi02::Array{Float64,2},Phi12::Array{Float64,2},PhiEta2::Array{Float64,2},pA,pA0,vStatePrev::Array{Float64,1}, vEta::Array{Float64,1})
#function fStateSpaceSimulation(sParam_StateSpace, vStatePrev::Adjoint{Float64,Array{Float64,2}}, vEta)
#function fStateSpaceSimulation(sParam_StateSpace, vStatePrev::Array{Float64,1}, vEta::Array{Float64,1})
#function fStateSpaceSimulation(eta1BarCoeff,Phi01,Phi11,PhiEta1,Phi02,Phi12,PhiEta2,pA,pA0,vStatePrev::Array{Float64,1}, vEta::Array{Float64,1})
#function fStateSpaceSimulation(sParam_StateSpace, vStatePrev::Array{Float64,2}, vEta)

#Changes

#=
eta1Bar type definition, remove concatenate ~35
coordinate of eta1Bar[1,1] to eta1Bar[1] in if ~37
define types of inputs ~6
Receive inputs unpackaged
=#

# Unpack parameters
#eta1BarCoeff = sParam_StateSpace["eta1BarCoeff"];
#Phi01        = sParam_StateSpace["Phi01"];
#Phi11        = sParam_StateSpace["Phi11"];
#PhiEta1      = sParam_StateSpace["PhiEta1"];
#Phi02        = sParam_StateSpace["Phi02"];
#Phi12        = sParam_StateSpace["Phi12"];
#PhiEta2      = sParam_StateSpace["PhiEta2"];
#pA           = sParam_StateSpace["pA"];
#pA0          = sParam_StateSpace["pA0"];

# Order of variable: er, g, z, d, y, pi, c, R
# Order of shocks:   er, eg, ez, ed

# Compute the regime indicator

#eta1Bar = dot(eta1BarCoeff,[1.0;vStatePrev])
#eta1Bar = eta1BarCoeff*[1.0;vStatePrev]
eta1Bar::Float64 = eta1BarCoeff[1]+dot(eta1BarCoeff[2:end],vStatePrev[1:8])

#println(size(Phi01))
#println(size(Phi11))
#println(size(vStatePrev))
#println(size(PhiEta1))
#println(size(vEta))

#vState=Array{Float64}(undef,8,1)

temp = vStatePrev[5];

if vEta[1] < eta1Bar[1]

    vState_t = Phi01 + Phi11*vStatePrev[1:8] + PhiEta1*vEta
    # vState = [ vState_t ; vStatePrev[5:5] ] # y_t-1 = y_t-1
    vState = vcat(vState_t , vStatePrev[5:5] ) # y_t-1 = y_t-1

   pRegime = 1

else

   vState_t = Phi02 + Phi12*vStatePrev[1:8] + PhiEta2*vEta
   # vState = [ vState_t ; vStatePrev[5:5] ] # y_t-1 = y_t-1
   vState = vcat(vState_t , vStatePrev[5:5] ) # y_t-1 = y_t-1

   pRegime = 2


end


vObs = pA0 + pA*vState

return vState, vObs, pRegime;

end
