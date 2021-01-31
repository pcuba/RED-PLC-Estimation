# Function to add steady states to parameter dictionary
function fAddSteadyStateParamDict(ParaDict_in)

ParaDict_out = ParaDict_in;

#--------------------------------------------------------------------------
#                       EXTRACT PARAMETERS
#--------------------------------------------------------------------------
tau       = ParaDict_in["tau"];
kappa     = ParaDict_in["kappa"];
eta       = ParaDict_in["eta"];
nu        = ParaDict_in["nu"];
chi_h     = ParaDict_in["chi_h"];
gstar     = ParaDict_in["gstar"];
rAnet     = ParaDict_in["rAnet"];
gamQnet   = ParaDict_in["gamQnet"];
piAnet    = ParaDict_in["piAnet"];
pibarAnet = ParaDict_out["pibarAnet"];

#--------------------------------------------------------------------------
#                       TRANSFORMED PARAMETERS
#--------------------------------------------------------------------------
pistar= exp(piAnet/400);                      # Annualized inflation rate
pibar = pistar;                               # Reference inflation rate for cost adjustment
beta  = exp(-rAnet/400);                      # Disccount factor
gamma = exp(gamQnet/100);                     # gross growth rate
r     = gamma/beta;                           # Real interest rate
phi   = tau*(1-nu)/((nu * pistar^2 * kappa)); # Adjustment cost parameter
b     = 1/(2*nu);                             # Constant

ParaDict_out["pistar"] = pistar
ParaDict_out["pibar"]  = pibar
ParaDict_out["beta"]   = beta
ParaDict_out["gamma"]  = gamma
ParaDict_out["r"]      = r
ParaDict_out["phi"]    = phi
ParaDict_out["b"]      = b

#--------------------------------------------------------------------------
#                   NON-STOCHASTIC STEADY STATE
#--------------------------------------------------------------------------
pi_ss  = pistar;
c_ss   = ( (1 - nu + phi*nu*(1-beta)*pi_ss*(pi_ss-pibar)-0.5*phi*(pi_ss-pibar)^2) / (chi_h*((1/(gstar))-0.5*phi*(pi_ss-pibar)^2)^(-1/eta)) )^(1/(tau + 1/eta));
y_ss   = c_ss/(1/(gstar) - 0.5 * phi * (pi_ss -pibar)^2);
R_ss   = pi_ss * gamma/beta;
ee1_ss = (c_ss^-tau)/(gamma*pi_ss);

ParaDict_out["pi_ss"]  = pi_ss
ParaDict_out["c_ss"]   = c_ss
ParaDict_out["y_ss"]   = y_ss
ParaDict_out["R_ss"]   = R_ss
ParaDict_out["ee1_ss"] = ee1_ss
ParaDict_out["cy"]     = c_ss/y_ss
ParaDict_out["gy"]     = gstar/y_ss
ParaDict_out["g_ss"]   = gstar

return ParaDict_out

end
