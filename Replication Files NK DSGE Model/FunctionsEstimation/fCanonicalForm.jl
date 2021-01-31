# -----------------------------------------------------------------------
# Computes PLC solution and transforms it into canonical form
# Regime 1 indicator:   eta1 < eta1BarCoeff'*[1,vStatePrev]
# Conversion of shocks: vEta = D_EpsToEta * vEps
#
# Regime j, j=1,2
# vState = Phi0j + Phi1j*vStatePrev + PhiEtaj*vEta
# output D_EpsToEta, eta1BarCoeff,
#        Phi01, Phi11, PhiEta1
#        Phi02, Phi12, PhiEta2
# -----------------------------------------------------------------------

# Updated: 01/03/21

function fCanonicalForm(sParam_steady,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

# solve model
#############
flag_out,theta_out,theta_in,res_out,sPLC_out = fPLCSolve(sParam_steady,sPLC,options_int_use,options_float_use,Smolyak_fixed_grid)

if flag_out > 0

    ## THIS NEEDS TO BE CLEANED UP EVENTUALLY

    # get coefficients of PLC decision rules
    ########################################
    #theta_pi_1, theta_yy_1, theta_cc_1, theta_pi_2, theta_yy_2, theta_cc_2, d0,d1,d2,d3,d4,d5 = get_coefs_plc(theta_out, sParam_steady);

    ALP, BET, delta, dALP22,dBET22 = fPLCCoefs(theta_out,sPLC_out);

    # MAP  PLC COEFFICIENTS
    index_1= sPLC_out["index_1"];
    index_2= sPLC_out["index_2"];

    theta_pi_1 = ALP[1][index_1];
    theta_pi_2 = ALP[1][index_2];

    theta_yy_1 = ALP[2][index_1];
    theta_yy_2 = ALP[2][index_2];

    theta_cc_1 = ALP[2][index_1];
    theta_cc_2 = ALP[2][index_2];

    theta_cc_1[5] = theta_cc_1[5] - 1;
    theta_cc_2[5] = theta_cc_2[5] - 1;

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

    theta_pi_2 = theta_pi_2[old_order];
    theta_yy_2 = theta_yy_2[old_order];
    theta_cc_2 = theta_cc_2[old_order];

    # delta function
    # 1 yhat(-1) z g d er
    d0 = delta[1]       # constant
    d1 = delta[6]       # er
    d2 = delta[4]       # g
    d3 = delta[3]       # z
    d4 = delta[5]       # d
    d5 = delta[2]       # yhat(-1)

    ################################


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

    # convert epsilons into etas
    ############################
    d_vec       = Array([d1;d2;d3;d4]')*Diagonal([sig_r,sig_g,sig_z,sig_d])
    d_vec_norm  = sqrt(d_vec*d_vec')[1]
    d_vec       = d_vec/d_vec_norm # normalize the length of this vector
    d_null      = nullspace(d_vec)
    D_EpsToEta  = [d_vec; d_null']
    D_EtaToEps  = inv(D_EpsToEta)

    # threshold for regime indicator
    ################################
    # Order of variables: intercept, er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed
    eta1BarCoeff = Array([-d0; 0; -d2*rho_g; -d3*rho_z; -d4*rho_d; -d5; 0; 0 ; 1]')/d_vec_norm

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

    # Rhats = (1-rho_r)*psi1*(t) + (1-rho_r)*psi2*y(t)
    Phi0r    = (1-rho_r)*(psi1*Phi0ypic[2,:] + psi2*Phi0ypic[1,:]);
    Phi1r    = (1-rho_r)*(psi1*Phi1ypic[2,:] + psi2*Phi1ypic[1,:]);
    PhiEpsr  = (1-rho_r)*(psi1*PhiEpsypic[2,:] + psi2*PhiEpsypic[1,:]);

    # Rhats = Rhats - (1-rho_r)*psi2*( y(t-1)+rhoz*z(t-1) ) + rho_r*Rhat_lag
    #         + (1-rho_r)*psi2*epsz(t) + sig_r*ers;
    Phi1r[3] = Phi1r[3] + (1-rho_r)*psi2*rho_z; # + (1-rho_r)*psi2*rhoz*z(t-1)
    Phi1r[5] = Phi1r[5] - (1-rho_r)*psi2; # - (1-rho_r)*psi2*y(t-1)
    Phi1r[8] = Phi1r[8] + rho_r; #rho_r*Rhat_lag

    PhiEpsr[1] = PhiEpsr[1] + sig_r; # + sig_r*epsr
    PhiEpsr[3] = PhiEpsr[3] + (1-rho_r)*psi2*sig_z; # + (1-rho_r)*psi2*sig_z*epsz

    # Combine the three sets of equations
    Phi01   = [ Phi0ex; Phi0ypic[1:3]; Phi0r ]
    Phi11   = [ Phi1ex; Phi1ypic[1:3,:]; Phi1r' ]
    PhiEta1 = [ PhiEpsex; PhiEpsypic[1:3,:]; PhiEpsr' ]*D_EtaToEps

    # Regime 2
    ##########

    # generate y(t), pi(t), c(t) equations
    # Order of variables: er, g, z, d, y, pi, c, R
    # Order of shocks: er, eg, ez, ed
    Phi0ypic = [ theta_yy_2[1]; theta_pi_2[1] ; theta_cc_2[1] ]
    Phi1ypic = ([ theta_yy_2[2:end]' ; theta_pi_2[2:end]' ; theta_cc_2[2:end]']
               *[Diagonal([0,rho_g,rho_z,rho_d,1,0]) zeros(6,1) [0;0;0;0;0;1] ] );

    PhiEpsypic = ([ theta_yy_2[2:end]' ; theta_pi_2[2:end]'; theta_cc_2[2:end]']
                 *[Diagonal([sig_r,sig_g,sig_z,sig_d]);zeros(2,4) ]);

    Phi0r    = log(1/R_ss);
    Phi1r    = zeros(8,1);
    PhiEpsr  = zeros(4,1);

    Phi02   = [ Phi0ex; Phi0ypic[1:3]; Phi0r ]
    Phi12   = [ Phi1ex; Phi1ypic[1:3,:]; Phi1r' ]
    PhiEta2 = [ PhiEpsex; PhiEpsypic[1:3,:]; PhiEpsr']*D_EtaToEps

    # c_t = y_t - g_t
    Phi01[7,:]   = Phi01[5,:] - Phi01[2,:]
    Phi11[7,:]   = Phi11[5,:] - Phi11[2,:]
    PhiEta1[7,:] = PhiEta1[5,:] - PhiEta1[2,:]
    Phi02[7,:]   = Phi02[5,:] - Phi02[2,:]
    Phi12[7,:]   = Phi12[5,:] - Phi12[2,:]
    PhiEta2[7,:] = PhiEta2[5,:] - PhiEta2[2,:]

    return flag_out, theta_out,  D_EpsToEta, eta1BarCoeff, Phi01, Phi11, PhiEta1, Phi02, Phi12, PhiEta2;

else

    return flag_out, nothing,  nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing;


end


end
