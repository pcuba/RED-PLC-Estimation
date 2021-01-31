function fLinearSolve(params_in::Dict{String}{Float64})

# Updated 01/05/2021

# Inputs : params_in : model parameters
#outputs

# Compute derivatives of f (equilibrium conditions)
nfx, nfxp, nfy, nfyp, vcv =  fSGUGetOrder1(params_in);

#Compute first-order approximation
eflag_SGU, gx,hx = fSGU_gx_hx(nfy,nfx,nfyp,nfxp);

if eflag_SGU == 1

    theta_out = [[gx[1,1];  0; gx[1,2:end]];        # pihat
                 [gx[2,1];  0; gx[2,2:end]]        # yhat
              ];

    return eflag_SGU, theta_out

else

    return eflag_SGU, nothing

end

end
