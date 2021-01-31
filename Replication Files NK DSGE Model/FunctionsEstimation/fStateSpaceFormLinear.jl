# This function creates the state space representation of the model
function fStateSpaceFormLinear(sParam)

# Changes
# Change vInits to Array instead of adjoint Array

# Compute steady states
#######################
sParam_SS  = fAddSteadyStateParamDict(sParam)

# Solve for canonical form
##########################
# Order of variables: er, g, z, d, y, pi, c, R, ylag
flag_out, theta_out, Phi00, Phi10, PhiEta0 = fCanonicalFormLinear(sParam_SS)

if flag_out == 1

    # Create dictionary
    ###################
    sParam_StateSpace = Dict()

    sParam_StateSpace["theta_out"]    = theta_out;
    sParam_StateSpace["Phi00"]        = Phi00;

    sParam_StateSpace["Phi10"]        = Phi10;
    sParam_StateSpace["PhiEta0"]      = PhiEta0;

    sParam_StateSpace["vInitS"]     = [sParam["er0"],
                                       sParam["g0"],
                                       sParam["z0"],
                                       sParam["d0"],
                                       sParam["c0"] + sParam["g0"],
                                       sParam["pi0"],
                                       sParam["c0"],
                                       sParam["R0star"] - log(sParam_SS["R_ss"]),
                                       sParam["c0"] + sParam["g0"] ]

    sig_r = sParam_SS["sig_r"];
    sig_g = sParam_SS["sig_g"];
    sig_z = sParam_SS["sig_z"];
    sig_d = sParam_SS["sig_d"];

    sParam_StateSpace["pSigmaTE"]   = [1.0 0. 0. 0.;
                                       0. 1.0 0. 0.;
                                       0. 0. 1.0 0.;
                                       0. 0. 0. 1.0];

    # observables: output, consumption, inflation interest
    # states: er, g, z, d, y, pi, c, R
    # Measurement equation: loadings
    #sParam_StateSpace["pA"]  = [0. 0. 100. 0. 100. 0. 0. 0. -100. 0.;
    #                            0. 0. 0. 0. -100. 0. 100. 0. 0. 0.;
    #                            0. 0. 0. 0. 0. 400. 0. 0. 0. 0.;
    #                            0. 0. 0. 0. 0. 0. 0. 400. 0. 0.];
    # Measurement equation: intercepts
    #sParam_StateSpace["pA0"] = [100*log(sParam_SS["gamma"]);
    #                            100*log(sParam_SS["cy"]);
    #                            400*log(sParam_SS["pi_ss"]);
    #                            400*log(sParam_SS["R_ss"])];

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
