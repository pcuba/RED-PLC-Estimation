# -----------------------------------------------------------------------
# particle filtering for NK model with ZLB (work 3):
#  M.E.  y_t = pB*s_t + merror_t,  merror_t ~ iid MN((0, 0), pSigmaME)
#  T.E.  s_t = f2(s_t-1) + e_t,  e_t ~ iid N((0, 0), pSigmaTE)
# -----------------------------------------------------------------------

#Changes

# Initilize arrays for iSim Loop ~66-68
# Devectorize perror and density calculations ~65-72
# Unpackage StateSpace before passing to fStateSpaceSimulation


function fCOPFaugExactZLB(mY, sParam_StateSpace, mSigmaMEtilde,ME_scale, nM, pInvSigmaMEtilde)

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

# housekeeping
ne        = size(pSigmaTE,1);
ns        = size(vInitS,1);
ny        = size(mY,2);
T         = size(mY,1);
yt        = mY;

# matrix to store restults
s_up_stat = zeros(T,ns,3);
s_lag_stat = zeros(T,ns,3);
shock_stat = zeros(T,ne,3);
regime_stat = zeros(T);

lik       = zeros(T,1);
Neff      = zeros(T,1);
mS_fore   = zeros(ns,nM);
mY_fore   = zeros(ny,nM);
mS_up     = zeros(ns,nM);

vRegime   = zeros(nM);

# initialization
#mS_up     = zeros(ns,nM)
mS_up     = repeat(vInitS,1,nM)
#println(mS_up[:,1]);
############################

weights = ones(nM,1);

#Initilize Arrays
perror    = Array{Float64}(undef,ny,nM)
logdensity= Array{Float64}(undef,nM)
adjDensity= Array{Float64}(undef,nM)

nFlagError = 0
# Rest of Steps
for tt = 1:T

    mS_lag = copy(mS_up)
    yy = yt[tt,:];

    # Propagate and and update
    (mShock, iDlog, mS_fore, mY_fore, vRegime) = fCondDistExactZLB_v5(mS_fore,mY_fore,vRegime,tt,pSigmaTE,pSigmaEta,pA,pA0,eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,mSigmaMEtilde,ME_scale,mS_up,yy,nM, pInvSigmaMEtilde);

    if isnan(maximum(iDlog))==1 || maximum(iDlog) == -Inf
        nFlagError = 1
        iDlog .= 1
    end

    # option 3
    iDlog          = iDlog'
    cMaxLogdensity = maximum(iDlog)
    adjDensity     = exp.(iDlog .- cMaxLogdensity);
    lik[tt,1]      = cMaxLogdensity + log(mean(adjDensity.*weights));

    # option 2
    #adjDensity  = iD1' + iD2';
    #lik[tt,1]   = log(mean(adjDensity.*weights));

    # option 1
    #cMaxLogdensity = maximum(logdensity)
    #adjDensity     = exp.(logdensity .- cMaxLogdensity);
    #lik[tt,1]        = cMaxLogdensity + log(mean(adjDensity.*weights));

    # Normalize weights
    weights = (adjDensity.*weights)./(mean(adjDensity.*weights));
    weights = pweights(weights);

    # Effective sample size
    Neff[tt,1] = (nM^2)/sum(weights.^2);

    # if resample == 1 && Neff(tt,1) <= N/2
    resample = 1; # always resample
    if resample == 1
        if nSistResFilter == 1
            # Sistematic Reampling function
            vParticleSelect = fSystematicResampling(weights'); # systematic resampling
            vParticleSelect = Int.(vParticleSelect)
            #vParticleSelect = convert(Vector{Int64},vParticleSelect)'
            mS_up   = dropdims(mS_fore[:,vParticleSelect']; dims = 2);
            mS_lag  = dropdims(mS_lag[:,vParticleSelect']; dims = 2);
            mShock  = dropdims(mShock[:,vParticleSelect']; dims = 2);
            vRegime = vRegime[vParticleSelect];
            weights = ones(nM,1);
        else
            # Different algorithms:
            # http://juliastats.github.io/StatsBase.jl/stable/sampling.html#Algorithms-1
            # direct_sample!   alias_sample!
            #mS_up = squeeze(mS_fore[:,sample(1:nM, weights, nM, replace = true)'],2);
            vParticleSelect = zeros(nM)
            vParticleSelect = StatsBase.sample(1:nM, weights, nM, replace = true)
            #vParticleSelect = StatsBase.direct_sample!(1:nM, weights, vParticleSelect)
            #vParticleSelect = StatsBase.alias_sample!(1:nM, weights, vParticleSelect)
            vParticleSelect = convert(Vector{Int64},vParticleSelect)'
            mS_up   = dropdims(mS_fore[:,vParticleSelect]; dims = 2);
            mS_lag  = dropdims(mS_lag[:,vParticleSelect]; dims = 2);
            mShock  = dropdims(mShock[:,vParticleSelect]; dims = 2);
            vRegime = vRegime[vParticleSelect];
            weights = ones(nM,1) ;
        end
    else
        mS_up = mS_fore;
    end

    # reduce dimension after one iteration.

    for is = 1:ns
        s_up_stat[tt,is,1] = quantile(sort(mS_up[is,:]),0.05);
        s_up_stat[tt,is,2] = mean(mS_up[is,:]);
        s_up_stat[tt,is,3] = quantile(sort(mS_up[is,:]),0.95);

        s_lag_stat[tt,is,1] = quantile(sort(mS_lag[is,:]),0.05);
        s_lag_stat[tt,is,2] = mean(mS_lag[is,:]);
        s_lag_stat[tt,is,3] = quantile(sort(mS_lag[is,:]),0.95);

    end

    for is = 1:ne

        shock_stat[tt,is,1] = quantile(sort(mShock[is,:]),0.05);
        shock_stat[tt,is,2] = mean(mShock[is,:]);
        shock_stat[tt,is,3] = quantile(sort(mShock[is,:]),0.95);

    end

    regime_stat[tt] = mean(vRegime)

end

if nFlagError == 1
    lik .= -1E20
end

return lik, s_up_stat, s_lag_stat, shock_stat, regime_stat, Neff;
end

# llik, s_up_stat, s_lag_stat, shock_stat, regime_stat, Neff;stat, shock_stat, regime_stat, Neff;
