# -----------------------------------------------------------------------
# particle filtering for NK model with ZLB (work 3):
#  M.E.  y_t = pB*s_t + merror_t,  merror_t ~ iid MN((0, 0), pSigmaME)
#  T.E.  s_t = f2(s_t-1) + e_t,  e_t ~ iid N((0, 0), pSigmaTE)
# -----------------------------------------------------------------------

#Changes

# Initilize arrays for iSim Loop ~50-51
# Devectorize perror and density calculations ~65-72
# Unpackage StateSpace before passing to fStateSpaceSimulation


function fBSPF(mY, sParam_StateSpace, mSigmaME, ME_scale, nM)

vInitS     = sParam_StateSpace["vInitS"];
pSigmaTE   = sParam_StateSpace["pSigmaTE"];
pA         = sParam_StateSpace["pA"];
pA0        = sParam_StateSpace["pA0"];
D_EpsToEta = sParam_StateSpace["D_EpsToEta"]

eta1BarCoeff = sParam_StateSpace["eta1BarCoeff"];
Phi01        = sParam_StateSpace["Phi01"];
Phi11        = sParam_StateSpace["Phi11"];
PhiEta1      = sParam_StateSpace["PhiEta1"];
Phi02        = sParam_StateSpace["Phi02"];
Phi12        = sParam_StateSpace["Phi12"];
PhiEta2      = sParam_StateSpace["PhiEta2"];


# housekeeping
ne        = size(pSigmaTE,1);
ns        = size(vInitS,1);
ny        = size(mY,2);
T         = size(mY,1);
yt        = mY;

# matrix to store restults
s_up_stat   = zeros(T,ns,3);
s_lag_stat  = zeros(T,ns,3);
shock_stat  = zeros(T,ne,3);
regime_stat = zeros(T);
lik         = zeros(T,1);
Neff        = zeros(T,1);
mS_fore     = zeros(ns,nM);
mY_fore     = zeros(ny,nM);

vRegime   = zeros(nM);

# initialization
# vInitS is the true state in period 0
mS_up     = repeat(vInitS,1,nM)

weights = ones(nM,1);

#Initilize Arrays
perror    = Array{Float64}(undef,ny,nM)
logdensity= Array{Float64}(undef,nM)
adjDensity= Array{Float64}(undef,nM)

# Rest of Steps
for tt = 1:T

    mS_lag = copy(mS_up)
    yy = yt[tt,:];

    # Propagate each particle
    vRegime   = zeros(nM);
    mShockEps = rand(MvNormal(zeros(size(pSigmaTE,1)),pSigmaTE), nM)';
    mShock = (mShockEps*D_EpsToEta')';
    for iSim = 1:nM
        (mS_fore[:,iSim], mY_fore[:,iSim], vRegime[iSim]) = fStateSpaceSimulationaug(eta1BarCoeff,Phi01,Phi11,PhiEta1,Phi02,Phi12,PhiEta2,pA,pA0, mS_up[:,iSim], mShock[:,iSim]);
        perror[:,iSim]   = yy - mY_fore[:,iSim];
        logdensity[iSim] = logpdf(MvNormal(zeros(size(yy,1)), mSigmaME), perror[:,iSim])
    end

    cMaxLogdensity = maximum(logdensity)
    adjDensity     = exp.(logdensity .- cMaxLogdensity)

    if all(logdensity .== -Inf)
        adjDensity .= 1
        println("All log densities are -Inf in period: ",tt)
    end

    # Store results
    lik[tt,1]        = cMaxLogdensity + log(mean(adjDensity.*weights));

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
        mS_up = copy(mS_fore);
    end

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

return lik, s_up_stat, s_lag_stat, shock_stat, regime_stat, Neff;
end
