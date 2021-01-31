# This function creates the state space representation of the model

# This version: 01/03/21

function fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

#Changes
# Change vInits to Array instead of adjoint Array

# Compute steady states
#######################
sParam_SS  = fAddSteadyStateParamDict(sParam)

# Solve for canonical form
##########################
# Order of variables: er, g, z, d, y, pi, c, R
flag_out, theta_out , D_EpsToEta, eta1BarCoeff, Phi01, Phi11, PhiEta1, Phi02, Phi12, PhiEta2 = fCanonicalForm(sParam_SS,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

if flag_out > 0

    # Create dictionary
    ###################
    sParam_StateSpace = Dict()

    sParam_StateSpace["theta_out"]    = theta_out;
    sParam_StateSpace["D_EpsToEta"]   = D_EpsToEta;
    sParam_StateSpace["eta1BarCoeff"] = eta1BarCoeff;
    sParam_StateSpace["Phi01"]        = Phi01;

    sParam_StateSpace["Phi11"]        = Phi11;
    sParam_StateSpace["PhiEta1"]      = PhiEta1;
    sParam_StateSpace["Phi02"]        = Phi02;

    sParam_StateSpace["Phi12"]        = Phi12;
    sParam_StateSpace["PhiEta2"]      = PhiEta2;

    #sParam_StateSpace["vInitS"]     = [0., 0., 0., 0., 0., 0., 0., 0., 0.];
    # Initial value Steady State: # Order of variables: er, g, z, d, y, pi, c, R, ylag0

    sParam_StateSpace["vInitS"]     = [sParam["er0"],
                                           sParam["g0"],
                                           sParam["z0"],
                                           sParam["d0"],
                                           sParam["c0"] + sParam["g0"],
                                           sParam["pi0"],
                                           sParam["c0"],
                                           sParam["R0star"] - log(sParam_SS["R_ss"]),
                                           sParam["c0"] + sParam["g0"] ]

    sParam_StateSpace["pSigmaTE"]   = [1.0 0. 0. 0.;
                                       0. 1.0 0. 0.;
                                       0. 0. 1.0 0.;
                                       0. 0. 0. 1.0];

    # observables: output, consumption, inflation interest
    # states: er, g, z, d, y, pi, c, R
    # Measurement equation: loadings
    sParam_StateSpace["pA"]  = [0. 0. 100. 0. 100. 0. 0. 0. -100.;
                                0. -100. 0. 0. 0. 0. 0. 0. 0.;
                                0. 0. 0. 0. 0. 400. 0. 0. 0.;
                                0. 0. 0. 0. 0. 0. 0. 400. 0.];
    # Measurement equation: intercepts
    sParam_StateSpace["pA0"] = [100*log(sParam_SS["gamma"]);
                                -100*log(sParam_SS["g_ss"]);
                                400*log(sParam_SS["pi_ss"]);
                                400*log(sParam_SS["R_ss"])];

    return flag_out, sParam_StateSpace;

else

    return flag_out, nothing;

end

end
