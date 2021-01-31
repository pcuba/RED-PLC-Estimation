# -----------------------------------------------------------------------
# Computes Linear solution and transforms it into canonical form
# Regime 1 indicator:   eta1 < eta1BarCoeff'*[1,vStatePrev]
# only 1 regime
# vState = Phi00 + Phi10*vStatePrev + PhiEta0*vEta
# output: theta_out, Phi00, Phi10, PhiEta0
# -----------------------------------------------------------------------
function fCanonicalFormLinear(sParam_steady)

# solve model
#############
flag_out, theta_out = fLinearSolve(sParam_steady)

if flag_out == 1

    theta_pi_1 = theta_out[1:7,1]
    theta_yy_1 = theta_out[8:14,1]

    theta_cc_1 = theta_yy_1;
    theta_cc_1[5] = theta_cc_1[5] - 1;

    # Current state order
    # 1        2  3       4 5 6 7
    # Rhat(-1) 1 yhat(-1) z g d er
    # old state order
    # 1 er g z d yhat(-1) Rhat(-1)

    # Map to old state order
    old_order = [2;7;5;4;6;3;1];

    theta_pi_1 = theta_pi_1[old_order];
    theta_yy_1 = theta_yy_1[old_order];
    theta_cc_1 = theta_cc_1[old_order];

    # Unpack parameters
    ###################
    rho_g = sParam_steady["rho_g"];
    rho_z = sParam_steady["rho_z"];
    rho_d = sParam_steady["rho_d"];
    sig_r = sParam_steady["sig_r"];
    sig_g = sParam_steady["sig_g"];
    sig_z = sParam_steady["sig_z"];
    sig_d = sParam_steady["sig_d"];

    rho_r = sParam_steady["rho_r"];
    psi1  = sParam_steady["psi1"];
    psi2  = sParam_steady["psi2"];
    R_ss  = sParam_steady["R_ss"];

    # NO NEED IN LINEAR MODEL, USE DIRECTLY EPSILONS
    # convert epsilons into etas
    ############################
    #d_vec_norm  = sqrt(d_vec*d_vec')[1]
    #d_vec       = d_vec/d_vec_norm # normalize the length of this vector
    #d_null      = nullspace(d_vec)
    #D_EpsToEta  = [d_vec; d_null']
    #D_EtaToEps  = inv(D_EpsToEta)

    # threshold for regime indicator
    ################################
    # Order of variables: intercept, er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed

    # Matrices for the evolution of exogenous STATES
    ################################################
    # Order of variables: er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed

    Phi0ex  = zeros(4,1)

    Phi1ex  = [0.0 0.0   0.0 0.0 0.0 0.0 0.0 0.0 ;
               0.0 rho_g 0.0 0.0 0.0 0.0 0.0 0.0 ;
               0.0 0.0 rho_z 0.0 0.0 0.0 0.0 0.0 ;
               0.0 0.0 0.0 rho_d 0.0 0.0 0.0 0.0]

    PhiEpsex  = [sig_r 0.0 0.0 0.0;
                 0.0 sig_g 0.0 0.0;
                 0.0 0.0 sig_z 0.0;
                 0.0 0.0 0.0 sig_d]

    # Matrices for Regime 1
    #######################

    # generate y(t), pi(t), c(t) equations
    # Order of variables: er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed
    Phi0ypic = [ theta_yy_1[1]; theta_pi_1[1] ; theta_cc_1[1] ];

    Phi1ypic = ( [ theta_yy_1[2:end]' ; theta_pi_1[2:end]' ; theta_cc_1[2:end]']
               *[Diagonal([0,rho_g,rho_z,rho_d,1,0]) zeros(6,1) [0;0;0;0;0;1] ] );

    PhiEpsypic = ([ theta_yy_1[2:end]' ; theta_pi_1[2:end]' ; theta_cc_1[2:end]']
                  *[Diagonal([sig_r,sig_g,sig_z,sig_d]);zeros(2,4) ]);

    # generate R(t) equation
    # Order of variables: er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed

    Phi0r    = (1-rho_r)*(psi1*Phi0ypic[2,:] + psi2*Phi0ypic[1,:]);
    Phi1r    = (1-rho_r)*(psi1*Phi1ypic[2,:] + psi2*Phi1ypic[1,:]);
    PhiEpsr  = (1-rho_r)*(psi1*PhiEpsypic[2,:] + psi2*PhiEpsypic[1,:]);

    Phi1r[3] = Phi1r[3] + (1-rho_r)*psi2*rho_z;
    Phi1r[5] = Phi1r[5] - (1-rho_r)*psi2;
    Phi1r[8] = Phi1r[8] + rho_r;

    PhiEpsr[1] = PhiEpsr[1] + sig_r;
    PhiEpsr[3] = PhiEpsr[3] + (1-rho_r)*psi2*sig_z;

    # Combine the three sets of equations
    Phi00   = [ Phi0ex; Phi0ypic[1:3]; Phi0r  ]
    Phi10   = [ Phi1ex; Phi1ypic[1:3,:]; Phi1r' ]
    PhiEta0 = [ PhiEpsex; PhiEpsypic[1:3,:]; PhiEpsr' ]

    # Augment the matrices so s_t-1 includes y_t-1
    Phi00 = [Phi00; 0]
    Phi10 = [Phi10 zeros(8)]
    Phi10 = [Phi10 ; [0 0 0 0 1 0 0 0 0]]
    PhiEta0 = [PhiEta0 ; [0 0 0 0]]

    return flag_out, theta_out, Phi00, Phi10, PhiEta0;

else

    return flag_out, nothing, nothing, nothing, nothing;

end

end
