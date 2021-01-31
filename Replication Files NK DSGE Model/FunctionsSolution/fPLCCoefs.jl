function fPLCCoefs(theta_in,sPLC)

# Map PLC objects
# ny = sPLC["ny"]
# nx = sPLC["nx"]
# n  = sPLC["n"]

index_11 = sPLC["index_11"];
index_12 = sPLC["index_12"];
index_21 = sPLC["index_21"];
index_22 = sPLC["index_22"];

SMAT    = sPLC["SMAT"]
GAM      = sPLC["GAM"]
OMEGAMAT = sPLC["OMEGAMAT"]
PHIMAT   = sPLC["PHIMAT"]

# Assing GAM objecs
GAM1 = GAM[1];
GAM2 = GAM[2];
GAMY = GAM[3];
GAMX = GAM[4];

# %---------------------------------------------
# % Fill out ALP_1, ALP_21, BET_1, BET_21
# %---------------------------------------------
ALP = Vector{Array{Float64,2}}(undef, sPLC["ny"][1])
for iy=1:sPLC["ny"][1]
    ALP[iy] = zeros(2*sPLC["n"][1],1);   # Would be better to initializer with missing
    ALP[iy][1:sPLC["n"][1]+1,1] = SMAT[iy]*theta_in;
end
#
BET = Vector{Array{Float64,2}}(undef, sPLC["nx"][1])
for ix=1:sPLC["nx"][1]
     BET[ix] = zeros(2*sPLC["n"][1],1);  # Would be better to initializer with missing
     BET[ix][1:sPLC["n"][1]+1,1] = PHIMAT[ix] + OMEGAMAT[ix]*theta_in;
end

# --------------
# COMPUTE DELTA
# --------------
sum_gamy12 = zeros(1,sPLC["n"][1]-1);
sum_gamy11 = 0;

sum_gamx12 = zeros(1,sPLC["n"][1]-1);
sum_gamx11 = 0;

sum_gamysmat12 = zeros(sPLC["n"][1]-1,length(theta_in));
sum_gamysmat11 = zeros(length(theta_in),1);

for iy=1:sPLC["ny"][1]
    sum_gamy12 = sum_gamy12 + GAMY[iy]*ALP[iy][index_12]';
    sum_gamy11 = sum_gamy11 + GAMY[iy]*ALP[iy][index_11][1];

    sum_gamysmat11 = sum_gamysmat11 + GAMY[iy]*SMAT[iy][index_11,:];
    sum_gamysmat12 = sum_gamysmat12 + GAMY[iy]*SMAT[iy][index_12,:];
end

for ix=1:sPLC["nx"][1]
    sum_gamx12 = sum_gamx12 + GAMX[ix]*BET[ix][index_12]';
    sum_gamx11 = sum_gamx11 + GAMX[ix]*BET[ix][index_11][1];
end

delta = vec(-(GAM2 + sum_gamy12 + sum_gamx12)/(GAM1[1] + sum_gamy11 +sum_gamx11));

# Numerator and denominator of delta function
a_theta   = (GAM1[1] + sum_gamy11);
B_theta   = -(GAM2 + sum_gamy12);

# derivative of numerator and denominator
dB_dtheta = -sum_gamysmat12;
da_dtheta = sum_gamysmat11;

# This is d(B(theta)/a(theta)) / d theta
ddeldtheta = (a_theta*dB_dtheta - transpose(da_dtheta*B_theta) ) / (a_theta^2);

# -----------------------------------------
# IMPOSE CONTINUITY TO GET ALP_22, BET_22
# -----------------------------------------
dALP22 = Vector{Array{Float64,2}}(undef, sPLC["ny"][1])
for iy=1:sPLC["ny"][1]
    ALP[iy][sPLC["n"][1]+2:end,1] = transpose([ALP[iy][index_11][1] - ALP[iy][index_21][1]]*delta') + ALP[iy][index_12];

    dALP22[iy] = transpose((SMAT[iy][index_11,:] - SMAT[iy][index_21,:])*delta') +
        (ALP[iy][index_11,:] - ALP[iy][index_21,:])[1]*ddeldtheta + SMAT[iy][index_12,:];
end

dBET22 = Vector{Array{Float64,2}}(undef, sPLC["nx"][1])
for ix=1:sPLC["nx"][1]
    BET[ix][sPLC["n"][1]+2:end,1] = transpose([BET[ix][index_11][1] - BET[ix][index_21][1]]*delta') + BET[ix][index_12];

    dBET22[ix] = transpose((OMEGAMAT[ix][index_11,:] - OMEGAMAT[ix][index_21,:])*delta') +
        (BET[ix][index_11,:] - BET[ix][index_21,:])[1]*ddeldtheta + OMEGAMAT[ix][index_12,:];
end


return ALP, BET, delta, dALP22,dBET22


end
