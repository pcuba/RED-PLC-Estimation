# -----------------------------------------------------------------------
# Kalman filter
#  M.E.  y_t = pA0 + pA*s_t + merror_t,  merror_t ~ iid MN(0, pSigmaME)
#  T.E.  s_t = phi0 + phi1*s_t-1 + phieta*e_t,  e_t ~ iid N((0, 0), pSigmaTE)
# -----------------------------------------------------------------------

function fKF(mY, sParam_StateSpace, mSigmaME)

vInitS     = sParam_StateSpace["vInitS"];
pSigmaTE   = sParam_StateSpace["pSigmaTE"];
pA         = sParam_StateSpace["pA"];
pA0        = sParam_StateSpace["pA0"];

Phi00      = sParam_StateSpace["Phi00"];
Phi10      = sParam_StateSpace["Phi10"];
PhiEta0    = sParam_StateSpace["PhiEta0"];



# housekeeping
ne        = size(pSigmaTE,1);
ns        = size(vInitS,1);
ny        = size(mY,2);
T         = size(mY,1);
yt        = mY;

# matrix to store restults
lik       = zeros(T,1);
Neff      = zeros(T,1);
mS_fore   = zeros(T,ns);
mY_fore   = zeros(T,ny);
mP_up     = zeros(T,ns,ns);
mP_fore   = zeros(T,ns,ns);
mS_up     = zeros(T,ns);
mError    = zeros(T,ny);

# initialization
mP_init = zeros(ns,ns)
#########################

#Initilize Arrays
perror = Array{Float64}(undef,ny)
density = Array{Float64}(undef,1)

# Rest of Steps
for tt = 1:T

    yy = yt[tt,:];

    # Prediction (forecast)
    if tt == 1
        mS_fore[tt,:] = (Phi00 + Phi10*vInitS)'
        mP_fore[tt,:,:] = Phi10*mP_init*Phi10' + PhiEta0*pSigmaTE*PhiEta0'

        # Update
        mF = pA * mP_fore[tt,:,:] * pA' + mSigmaME
        mError[tt,:]  = yy - pA0 - pA*mS_fore[tt,:]

        mS_up[tt,:] = (mS_fore[tt,:] + mP_fore[tt,:,:]*pA'* (mF^(-1)) * mError[tt,:] )'
        mP_up[tt,:,:] = mP_fore[tt,:,:] - mP_fore[tt,:,:]*pA' *(mF^(-1))  * pA * mP_fore[tt,:,:]
    else
        mS_fore[tt,:] = (Phi00 + Phi10*mS_up[tt-1,:] )'
        mP_fore[tt,:,:] = Phi10*mP_up[tt-1,:,:]*Phi10' + PhiEta0*pSigmaTE*PhiEta0'

        # Update
        mF = pA * mP_fore[tt,:,:] * pA' + mSigmaME
        mError[tt,:]  = yy - pA0 - pA*mS_fore[tt,:]

        mS_up[tt,:] = (mS_fore[tt,:] + mP_fore[tt,:,:]*pA'* (mF^(-1)) * mError[tt,:] )'
        mP_up[tt,:,:] = mP_fore[tt,:,:] - mP_fore[tt,:,:]*pA' *(mF^(-1))  * pA * mP_fore[tt,:,:]

    end

    # Log likelihood
    lik[tt,1] = -0.5*ny*log(2*pi) - 0.5*log(det(mF)) - 0.5*mError[tt,:]'*(mF^(-1))*mError[tt,:];
end

return lik, mS_up;
end
