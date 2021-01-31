function [ALP, BET, delta]=get_coefs_plc_general(theta_in,sModel)

GAM1 = sModel.GAM{1};
GAM2 = sModel.GAM{2};
GAMY = sModel.GAM{3};
GAMX = sModel.GAM{4};
SMAT = sModel.SMAT;
PHIMAT   = sModel.PHIMAT;
OMEGAMAT = sModel.OMEGAMAT;

index_11 = sModel.index_11;
index_12 = sModel.index_12;
index_21 = sModel.index_21;
index_22 = sModel.index_22;

%---------------------------------------------
% Fill out ALP_1, ALP_21, BET_1, BET_21
%---------------------------------------------
for iy=1:sModel.ny
    ALP{iy} = NaN(2*sModel.n,1);
    ALP{iy}(1:sModel.n+1,1) = SMAT{iy}*theta_in;
end

for ix=1:sModel.nx
    BET{ix} = NaN(2*sModel.n,1);
   BET{ix}(1:sModel.n+1,1) = PHIMAT{ix} + OMEGAMAT{ix}*theta_in;
end

%--------------
% COMPUTE DELTA
%--------------
sum_gamy12 = zeros(1,sModel.n-1);
sum_gamy11 = 0;

sum_gamx12 = zeros(1,sModel.n-1);
sum_gamx11 = 0;

for iy=1:sModel.ny
    sum_gamy12 = sum_gamy12 + GAMY(iy)*ALP{iy}(index_12)';
    sum_gamy11 = sum_gamy11 + GAMY(iy)*ALP{iy}(index_11);
end

for ix=1:sModel.nx
    sum_gamx12 = sum_gamx12 + GAMX(ix)*BET{ix}(index_12)';
    sum_gamx11 = sum_gamx11 + GAMX(ix)*BET{ix}(index_11);
end


delta = -(GAM2 + sum_gamy12 + sum_gamx12)/(GAM1 + sum_gamy11 +sum_gamx11);

%-----------------------------------------
% IMPOSE CONTINUITY TO GET ALP_22, BET_22
%-----------------------------------------
for iy=1:sModel.ny
    ALP{iy}(sModel.n+2:end,1) = (ALP{iy}(index_11) - ALP{iy}(index_21))*delta' + ALP{iy}(index_12);
end

for ix=1:sModel.nx
    BET{ix}(sModel.n+2:end,1) = (BET{ix}(index_11) - BET{ix}(index_21))*delta' + BET{ix}(index_12);
end


end
