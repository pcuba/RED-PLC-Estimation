%--------------------------------------------------------------------------
% Script to generate the Posterior Distribution of Estimated Parameters
%
% This script reproduces Table 2: Posterior Distribution (PLC / COPF)
% of "Piecewise-Linear Approximations and Filtering for DSGE Models 
% with Occasionally Binding Constraints"  Aruoba, Cuba-Borda, Higa-Flores, 
% Schorfheide and Villalvazo (2020).
%
% Written by: Pablo Cuba-Borda
% Created : January 15, 2020
% This version : January 13, 2021
%--------------------------------------------------------------------------

clear; clc;


%=== DIRECTORIES WITH POSTERIOR DRAWS ===%

% COPF
modeldir   = 'PosteriorDraws/COPFexactZLB_Prior1_US_Mhrun1/';
modelname  = 'COPFexactZLB_Prior1_US_Mhrun1_';


modeldirKF   = 'PosteriorDraws/KF_Prior1_US_Mhrun1/';
modelnameKF  = 'KF_Prior1_US_Mhrun1_';


Posterior   = csvread([modeldir modelname 'Summary.csv']);
PosteriorKF = csvread([modeldirKF modelnameKF 'Summary.csv']);

% Rows of Posterior are:
% [vMeanTheta; vStdDevTheta; mHPDinterval; vMPDEstimator']

% Select for Table
% ParaDict_out["tau"]       = vPara[1];
% ParaDict_out["kappa"]     = vPara[2];
% ParaDict_out["rho_r"]     = vPara[5];
% ParaDict_out["rho_g"]     = vPara[6];
% ParaDict_out["rho_d"]     = vPara[7];
% ParaDict_out["rho_z"]     = vPara[8];
% ParaDict_out["sig_r"]     = vPara[9];
% ParaDict_out["sig_g"]     = vPara[10];
% ParaDict_out["sig_d"]     = vPara[11];
% ParaDict_out["sig_z"]     = vPara[12];
% ParaDict_out["gstar"]     = vPara[16];
% ParaDict_out["rAnet"]     = vPara[17];
% ParaDict_out["gamQnet"]   = vPara[18];
% ParaDict_out["piAnet"]    = vPara[19];

colindex = [1 2 5 6 7 8 9 10 11 12 16 17 18 19];

Posterior(:,9) = Posterior(:,9)*100;
Posterior(:,10) = Posterior(:,10)*100;
Posterior(:,11) = Posterior(:,11)*100;
Posterior(:,12) = Posterior(:,12)*100;

PosteriorKF(:,9) = PosteriorKF(:,9)*100;
PosteriorKF(:,10) = PosteriorKF(:,10)*100;
PosteriorKF(:,11) = PosteriorKF(:,11)*100;
PosteriorKF(:,12) = PosteriorKF(:,12)*100;


colname  = {'$\tau$','$\kappa$','$\rho_R$','$\rho_g$','$\rho_d$','$\rho_z$',...
    '$100\sigma_R$','$100\sigma_g$','$100\sigma_d$','$100\sigma_z$', '$g^*$  ','rAnet  ','gamQnet  ','piAnet'};

clc;

fprintf('\n =============== COPF RESULTS ====================');
fprintf('\n Parameter \t   Mean     MAP     HPD Low  HPD High');
for jj=1:length(colindex)
fprintf('\n %s \t & %2.2f  &  %2.2f  &   %2.2f  &  %2.2f \\\\',colname{jj},Posterior([1 5 3 4],colindex(jj)));
end

fprintf('\n \n');
fprintf('\n =============== KF RESULTS ====================');
fprintf('\n Parameter \t   Mean     MAP     HPD Low  HPD High');
for jj=1:length(colindex)
fprintf('\n %s \t & %2.2f  &  %2.2f  &   %2.2f  &  %2.2f \\\\',colname{jj},PosteriorKF([1 5 3 4],colindex(jj)));
end


