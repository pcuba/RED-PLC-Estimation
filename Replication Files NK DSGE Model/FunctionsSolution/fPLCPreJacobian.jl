function fPLCPreJacobian(ALP,delta,dALP22,dBET22,par,sPLC,state_i,index_states_sgu_1,index_states_sgu_2,SMAT_1_1,SMAT_2_1)

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
 rhom  = 0.0;
 sigz  = par["sig_z"] ;
 sigg  = par["sig_g"] ;
 sigd  = par["sig_d"] ;
 sigm  = par["sig_r"] ;
 gam   = par["gamma"] ;
 pist  = par["pi_ss"] ;
 pibar = par["pibar"] ;
 bet   = par["beta"] ;
 r     = par["r"] ;
 phi   = par["phi"] ;
 b     = par["b"] ;
 inveta= 1.0/eta;
 R_ss  = par["R_ss"] ;
 c_ss  = par["c_ss"] ;
 y_ss  = par["y_ss"] ;
 pist  = par["pi_ss"] ;

index_1= sPLC["index_1"];
index_2= sPLC["index_2"];

index_11= sPLC["index_11"];
index_12= sPLC["index_12"];
index_21= sPLC["index_21"];
index_22= sPLC["index_22"];


SMAT     = sPLC["SMAT"];
GAM      = sPLC["GAM"];
PHIMAT   = sPLC["PHIMAT"];
OMEGAMAT = sPLC["OMEGAMAT"];

GAM2     = GAM[2];

# %--------------------------------------------------------------------------
# % Get controls and states:
# %state_i = [Rhat0,1,yhat0,z,g,d,er];
# %--------------------------------------------------------------------------
Rlag1 = state_i[1];
ylag1 = state_i[3];
z1    = state_i[4];
g1    = state_i[5];
d1    = state_i[6];
m1    = state_i[7];

# %--------------------------------------------------------------------------
# % Derivative of time-t decision rules with respect to coefficients
# %--------------------------------------------------------------------------
dY = zeros(sPLC["ny"],sPLC["my"]);
dX = zeros(sPLC["nx"],sPLC["my"]);

if state_i[1]>dot(delta,state_i[2:end])
    # % Approx. Policies at current point
    p1    = dot(state_i,ALP[1][index_1]);
    y1    = dot(state_i,ALP[2][index_1]);

    #  dYtheta[i}
    dY[1,:] = state_i'*SMAT_1_1;
    dY[2,:] = state_i'*SMAT_2_1;

    #  dXdtheta{i}
    dX[1,:] = state_i'*OMEGAMAT[1][index_1,:];
    dX[2,:] = state_i'*OMEGAMAT[2][index_1,:];
    dX[3,:] = state_i'*OMEGAMAT[3][index_1,:];
    dX[4,:] = state_i'*OMEGAMAT[4][index_1,:];
    dX[5,:] = state_i'*OMEGAMAT[5][index_1,:];
    dX[6,:] = state_i'*OMEGAMAT[6][index_1,:];


    bind = 0;

else
    # % Approx. Policies at current point
    p1    = dot(state_i,ALP[1][index_2]);
    y1    = dot(state_i,ALP[2][index_2]);

    # % dYdtheta{i}
    dY[1,:] = state_i[1]*SMAT[1][index_21,:]' + state_i[2:end]'*dALP22[1];
    dY[2,:] = state_i[1]*SMAT[2][index_21,:]' + state_i[2:end]'*dALP22[2];

    # % dXdtheta{i}
    dX[1,:] = state_i[1]*OMEGAMAT[1][index_21,:]' + state_i[2:end]'*dBET22[1];
    dX[2,:] = state_i[1]*OMEGAMAT[2][index_21,:]' + state_i[2:end]'*dBET22[2];
    dX[3,:] = state_i[1]*OMEGAMAT[3][index_21,:]' + state_i[2:end]'*dBET22[3];
    dX[4,:] = state_i[1]*OMEGAMAT[4][index_21,:]' + state_i[2:end]'*dBET22[4];
    dX[5,:] = state_i[1]*OMEGAMAT[5][index_21,:]' + state_i[2:end]'*dBET22[5];
    dX[6,:] = state_i[1]*OMEGAMAT[6][index_21,:]' + state_i[2:end]'*dBET22[6];

    bind = 1;

end

# %--------------------------------------------------------------------------
# % COLLECT DERIVATIVE WITH RESPECT TO COEFFICIENTS
# %--------------------------------------------------------------------------
dXdtheta = dX;
dYdtheta = dY;

dYprime_dXprime_1 = [ALP[1][index_states_sgu_1]';ALP[2][index_states_sgu_1]'];
dYprime_dXprime_2 = [ALP[1][index_states_sgu_2]';ALP[2][index_states_sgu_2]'];


return dXdtheta,dYdtheta,dYprime_dXprime_1,dYprime_dXprime_2
end
