function fPLCRestrictions(par,sPLC)

# Map parameters
rhor = par["rho_r"]
rhoz = par["rho_z"]
rhog = par["rho_g"]
rhod = par["rho_d"]
rhom = 0;
psi1 = par["psi1"]
psi2 = par["psi2"]
r    = par["r"]
pist = par["pi_ss"]

#  Map PLC Matrices
SMAT     = sPLC["SMAT"];
GAM      = sPLC["GAM"];
PHIMAT   = sPLC["PHIMAT"];
OMEGAMAT = sPLC["OMEGAMAT"];

# ------------------------------------------------
# PLC LINEAR RESTRICTIONS (MODEL SPECIFIC)
# ------------------------------------------------
# GAM1             Rhat(-1)

#****** HOW TO MAP GAM[1] = [0]
# RED Paper specification
GAM[1][1]  = rhor;

# GAM2       1            yhat(-1)          zhat            ghat   dhat    mphat
GAM[2]    = [log(r*pist)  -(1-rhor)*psi2    (1-rhor)*psi2   0      0        1];

# GAMY      yhat             pihat
GAM[3]    = [(1-rhor)*psi1   (1-rhor)*psi2];

# GAMX      Rhat(-1)  Yhat(-1)  zhat  ghat dhat mphat
GAM[4]    = [0        0         0     0    0    0   ];


# ------------------------------------------------
# PLC TRANSITION STATE VARIABLES (MODEL SPECIFIC)
# ------------------------------------------------

#                 Rhat(-1) 1  yhat(-1)       zhat          ghat dhat mphat
PHIMAT[1] = [rhor   0  -(1-rhor)*psi2 (1-rhor)*psi2 0    0    1      0]';    # PHI(R)
PHIMAT[2] = [ 0     0  0               0            0    0    0      0]';    # PHI(y)
PHIMAT[3] = [ 0     0  0               rhoz         0    0    0      0]';    # PHI(z)
PHIMAT[4] = [ 0     0  0               0            rhog 0    0      0]';    # PHI(g)
PHIMAT[5] = [ 0     0  0               0            0    rhod 0      0]';    # PHI(d)
PHIMAT[6] = [ 0     0  0               0            0    0    rhom   0]';    # PHI(mp)


#OMEGA(R)
#OMEGAMAT[1][1:sPLC["n"][1],1:sPLC["n"][1]]                  = (1-rhor)*psi1*I(sPLC["n"][1]);   # pi(t)
#OMEGAMAT[1][1:sPLC["n"][1],sPLC["n"][1]+1:2*sPLC["n"][1]] = (1-rhor)*psi2*I(sPLC["n"][1]);   # y(t)
OMEGAMAT[1][1:sPLC["n"][1],1:sPLC["n"][1]]                  = (1-rhor)*psi1*Matrix(I,sPLC["n"][1],sPLC["n"][1]);   # KHF Use Julia Notation
OMEGAMAT[1][1:sPLC["n"][1],sPLC["n"][1]+1:2*sPLC["n"][1]] = (1-rhor)*psi2*Matrix(I,sPLC["n"][1],sPLC["n"][1]);   # KHF Use Julia Notation



#OMEGA(y)
#OMEGAMAT[2][1:sPLC["n"][1],sPLC["n"][1]+1:2*sPLC["n"][1]] = I(sPLC["n"][1]);                # OMEGA(y)
OMEGAMAT[2][1:sPLC["n"][1],sPLC["n"][1]+1:2*sPLC["n"][1]] = Matrix(I,sPLC["n"][1],sPLC["n"][1]); #KHF Use Julia Notation               # OMEGA(y)
OMEGAMAT[2][sPLC["n"][1]+1,2*sPLC["n"][1]+2]                = 1;                                # OMEGA(y)


# Update matrices
#sPLC = Dict("GAM" => GAM, "SMAT" => SMAT , "PHIMAT" => PHIMAT, "OMEGAMAT" => OMEGAMAT);
sPLC["GAM"] = GAM;
sPLC["SMAT"] = SMAT;
sPLC["PHIMAT"] = PHIMAT;
sPLC["OMEGAMAT"] = OMEGAMAT;

return sPLC

end
