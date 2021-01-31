function [ALP, BET, delta,dALP22,dBET22]=get_coefs_plc_deriv(theta_in,par)

GAM1 = par.GAM1;
GAM2 = par.GAM2;
GAMY = par.GAMY;
GAMX = par.GAMX;
SMAT = par.SMAT;
PHIMAT   = par.PHIMAT;
OMEGAMAT = par.OMEGAMAT;

index_11 = par.index_11;
index_12 = par.index_12;
index_21 = par.index_21;
index_22 = par.index_22;

%---------------------------------------------
% Fill out ALP_1, ALP_21, BET_1, BET_21
%---------------------------------------------
for iy=1:par.ny
    ALP{iy} = NaN(2*par.n,1);
    ALP{iy}(1:par.n+1,1) = SMAT{iy}*theta_in;
end

for ix=1:par.nx
    BET{ix} = NaN(2*par.n,1);
    BET{ix}(1:par.n+1,1) = PHIMAT{ix} + OMEGAMAT{ix}*theta_in;
end

%--------------
% COMPUTE DELTA
%--------------
sum_gamy12 = zeros(1,par.n-1);
sum_gamy11 = 0;

sum_gamx12 = zeros(1,par.n-1);
sum_gamx11 = 0;

for iy=1:par.ny
    sum_gamy12 = sum_gamy12 + GAMY(iy)*ALP{iy}(index_12)';
    sum_gamy11 = sum_gamy11 + GAMY(iy)*ALP{iy}(index_11);
end

for ix=1:par.nx
    sum_gamx12 = sum_gamx12 + GAMX(ix)*BET{ix}(index_12)';
    sum_gamx11 = sum_gamx11 + GAMX(ix)*BET{ix}(index_11);
end


delta = -(GAM2 + sum_gamy12 + sum_gamx12)/(GAM1 + sum_gamy11 +sum_gamx11);


% ddeldtheta is model specific need a general expression:
ddeldtheta = GAM2'*OMEGAMAT{1}(index_11,:)/(BET{1}(index_11))^2 - ...
    (BET{1}(index_11)*OMEGAMAT{1}(index_12,:) - BET{1}(index_12)*OMEGAMAT{1}(index_11,:))/(BET{1}(index_11))^2;

%-----------------------------------------
% IMPOSE CONTINUITY TO GET ALP_22, BET_22
%-----------------------------------------

for iy=1:par.ny
    ALP{iy}(par.n+2:end,1) = (ALP{iy}(index_11) - ALP{iy}(index_21))*delta' + ALP{iy}(index_12);
    
    dALP22{iy} = delta'*(SMAT{iy}(index_11,:) - SMAT{iy}(index_21,:)) + ...
        (ALP{iy}(index_11,:) - ALP{iy}(index_21,:))*ddeldtheta + SMAT{iy}(index_12,:);
    
end

for ix=1:par.nx
    BET{ix}(par.n+2:end,1) = (BET{ix}(index_11) - BET{ix}(index_21))*delta' + BET{ix}(index_12);
    
    dBET22{ix} = delta'*(OMEGAMAT{ix}(index_11,:) - OMEGAMAT{ix}(index_21,:)) + ...
        (BET{ix}(index_11,:) - BET{ix}(index_21,:))*ddeldtheta + OMEGAMAT{ix}(index_12,:);
    
end


end