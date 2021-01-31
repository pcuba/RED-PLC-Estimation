function fPLCGetJacobian(SMAT::Array{Array{Float64,2},1}, ALP::Array{Array{Float64,2},1}, delta::Array{Float64,1}, dALP22::Array{Array{Float64,2},1}, dBET22::Array{Array{Float64,2},1}, par::Dict{String,Float64}, sPLC, state_i::Array{Float64,1}, state_i_prime::Array{Float64,1}, dXdtheta::Array{Float64,2}, dYdtheta::Array{Float64,2}, index_states_sgu_1::Array{Int64,1}, index_states_sgu_2::Array{Int64,1}, p1::Float64, y1::Float64, p2::Float64, y2::Float64, dYprime_dXprime_1::Array{Float64,2}, dYprime_dXprime_2::Array{Float64,2}, SMAT_1::Array{Float64,2}, SMAT_1_1,SMAT_2::Array{Float64,2}, SMAT_2_1)

#  Get parameters

 eta   = par["eta"] ;
 nu    = par["nu"] ;
 chi_h = par["chi_h"] ;
 cy    = par["cy"] ;
 gst   = par["gstar"] ;
 tau   = par["tau"] ;
 kappa = par["kappa"] ;
 psi1  = par["psi1"] ;
 psi2  = par["psi2"] ;
 rhor  = par["rho_r"] ;
 rhoz  = par["rho_z"];
 rhog  = par["rho_g"];
 rhod  = par["rho_d"];
 gam   = par["gamma"] ;
 pist  = par["pi_ss"] ;
 pibar = par["pibar"] ;
 bet   = par["beta"] ;
 r     = par["r"] ;
 phi   = par["phi"] ;
 b     = par["b"] ;
 inveta= 1.0/eta;
 c_ss  = par["c_ss"] ;
 y_ss  = par["y_ss"] ;

index_1= sPLC["index_1"];
index_2= sPLC["index_2"];

index_11= sPLC["index_11"];
index_12= sPLC["index_12"];
index_21= sPLC["index_21"];
index_22= sPLC["index_22"];

# %--------------------------------------------------------------------------
# % Get controls and states:
# %state_i = [Rhat0,1,yhat0,z,g,d,er];
# %--------------------------------------------------------------------------
Rlag1::Float64 = state_i[1];
ylag1::Float64 = state_i[3];
z1::Float64    = state_i[4];
g1::Float64    = state_i[5];
d1::Float64    = state_i[6];
m1::Float64    = state_i[7];

z2::Float64 = state_i_prime[4];
g2::Float64 = state_i_prime[5];
d2::Float64 = state_i_prime[6];
m2::Float64 = state_i_prime[7];

#Preallocate memory

nfxp::Array{Float64,2}=zeros(Float64,2,6)
nfy::Array{Float64,2}=zeros(Float64,2,2)
nfyp::Array{Float64,2}=zeros(Float64,2,2)


dYprime=Array{Float64,2}(undef, sPLC["my"], sPLC["ny"])

if state_i[1]>dot(delta,state_i[2:end])
    bind = 0;
else
    bind = 1;
end

# %--------------------------------------------------------------------------
# % Derivative of time-t+1 decision rules with respect to coefficients
# %--------------------------------------------------------------------------


if state_i_prime[1]>dot(delta,state_i_prime[2:end])

    @inbounds dYprime[:,1] = state_i_prime'*SMAT_1_1;
    @inbounds dYprime[:,2] = state_i_prime'*SMAT_2_1;

    dYprime_dXprime = dYprime_dXprime_1

else

    @inbounds dYprime[:,1] = state_i_prime[1]*SMAT[1][index_21,:]' + state_i_prime[2:end]'*dALP22[1];
    @inbounds dYprime[:,2] = state_i_prime[1]*SMAT[2][index_21,:]' + state_i_prime[2:end]'*dALP22[2];

    dYprime_dXprime = dYprime_dXprime_2
end

# %--------------------------------------------------------------------------
# % COMPUTE THE REGIME DEPENDENT JACOBIAN
# %--------------------------------------------------------------------------

if bind == 0  #%--Regime 1--%

    sub1::Float64=exp(d2 - d1 + m1 - p2 - z2 + Rlag1*rhor - tau*(g1 - y1) + tau*(g2 - y2) - p1*psi1*(rhor - 1) - psi2*(rhor - 1)*(y1 - ylag1 + z1))

    nfxp[1,3] = sub1;
    nfxp[1,4] = -tau*sub1;
    nfxp[2,4] = bet*nu*phi*pist*tau*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2));
    nfxp[1,5] = -sub1;
    nfxp[2,5] = bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2));

    nfy[1,1] = psi1*sub1*(rhor - 1);
    nfy[2,1] = nu*phi*pist^2*exp(2*p1) + phi*pist*exp(p1)*(pibar - pist*exp(p1)) - nu*phi*pist*exp(p1)*(pibar - pist*exp(p1))
    nfy[1,2] = -sub1*(tau - psi2*(rhor - 1));
    nfy[2,2] = bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2))*(tau - 1) - c_ss^tau*chi_h*y_ss^inveta*exp(inveta*y1 - tau*(g1 - y1))*(inveta + tau);

    nfyp[1,1] = exp(d2 - d1 + m1 - p2 - z2 + Rlag1*rhor - tau*(g1 - y1) + tau*(g2 - y2) - p1*psi1*(rhor - 1) - psi2*(rhor - 1)*(y1 - ylag1 + z1));
    nfyp[2,1] = bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2)) - bet*nu*phi*pist^2*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*exp(p2);
    nfyp[1,2] = tau*exp(d2 - d1 + m1 - p2 - z2 + Rlag1*rhor - tau*(g1 - y1) + tau*(g2 - y2) - p1*psi1*(rhor - 1) - psi2*(rhor - 1)*(y1 - ylag1 + z1));
    nfyp[2,2] = -bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2))*(tau - 1);

else         #%--Regime 0--%

    sub1=exp(d2 - d1 - p2 - z2 + log(bet/(gam*pist)) - tau*(g1 - y1) + tau*(g2 - y2))

    nfxp[1,3] = sub1;
    nfxp[1,4] = -tau*sub1;
    nfxp[2,4] = bet*nu*phi*pist*tau*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2));
    nfxp[1,5] = -sub1;
    nfxp[2,5] = bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2));

    nfy[2,1] =  nu*phi*pist^2*exp(2*p1) + phi*pist*exp(p1)*(pibar - pist*exp(p1)) - nu*phi*pist*exp(p1)*(pibar - pist*exp(p1));                
    nfy[1,2] =  -tau*sub1;
    nfy[2,2] =  bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2))*(tau - 1) - c_ss^tau*chi_h*y_ss^inveta*exp(inveta*y1 - tau*(g1 - y1))*(inveta + tau);

    nfyp[1,1] = sub1;
    nfyp[2,1] = bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2)) - bet*nu*phi*pist^2*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*exp(p2);
    nfyp[1,2] = tau*sub1;
    nfyp[2,2] = -bet*nu*phi*pist*exp(d2 - d1 + p2 - y1 + y2 - tau*(g1 - y1) + tau*(g2 - y2))*(pibar - pist*exp(p2))*(tau - 1);

end



# %--------------------------------------------------------------------------
# % COLLECT DERIVATIVE WITH RESPECT TO COEFFICIENTS
# %--------------------------------------------------------------------------

dYprimedtheta = dYprime'

# %--------------------------------------------------------------------------
# % CONSTRUCT JACOBIAN
# %--------------------------------------------------------------------------

jacout::Array{Float64,2} = nfxp*dXdtheta .+ nfy*dYdtheta .+ nfyp*(dYprime_dXprime*dXdtheta .+ dYprimedtheta);

jacout_vec = vec(jacout);

return jacout_vec
end
