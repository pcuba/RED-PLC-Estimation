function sModel = fLoadCanonicalForm(modeldir,modelname)

% LOAD MODEL PARAMETERS
sParam = fCreateParams([modeldir modelname]);


% LOAD SOLUTION MATRICES
D_EpsToEta   = csvread([modeldir modelname 'Param_D_EpsToEta.csv']);
eta1BarCoeff = csvread([modeldir modelname 'Param_eta1BarCoeff.csv']);
Phi01        = csvread([modeldir modelname 'Param_Phi01.csv']);
Phi11        = csvread([modeldir modelname 'Param_Phi11.csv']);
PhiEta1      = csvread([modeldir modelname 'Param_PhiEta1.csv']);
Phi02        = csvread([modeldir modelname 'Param_Phi02.csv']);
Phi12        = csvread([modeldir modelname 'Param_Phi12.csv']);
PhiEta2      = csvread([modeldir modelname 'Param_PhiEta2.csv']);


% VAR-COV MATRIX OF SHOCKS
pSigmaTE   = [1.0 0. 0. 0.;
              0. 1.0 0. 0.;
              0. 0. 1.0 0.;
              0. 0. 0. 1.0];

% # Measurement equation: loadings          
% # observables: output, consumption, inflation interest
% # states: er, g, z, d, y, pi, c, R, ylag
%     sParam_StateSpace["pA"]  = [0. 0. 100. 0. 100. 0. 0. 0. -100.;
%                                 0. -100. 0. 0. 0. 0. 0. 0. 0.;
%                                 0. 0. 0. 0. 0. 400. 0. 0. 0.;
%                                 0. 0. 0. 0. 0. 0. 0. 400. 0.];
%     # Measurement equation: intercepts
%     sParam_StateSpace["pA0"] = [100*log(sParam_SS["gamma"]);
%                                 -100*log(sParam_SS["g_ss"]);
%                                 400*log(sParam_SS["pi_ss"]);
%                                 400*log(sParam_SS["R_ss"])];


pA  =  [0. 0. 100. 0. 100. 0. 0. 0. -100.;
        0. -100. 0. 0. 0. 0. 0. 0. 0.;
        0. 0. 0. 0. 0. 400. 0. 0. 0.;
        0. 0. 0. 0. 0. 0. 0. 400. 0.];
   
% # Measurement equation: intercepts
pA0 = [100*log(sParam.gamma);
       -100*log(sParam.gstar);
       400*log(sParam.pi_ss);
       400*log(sParam.R_ss)];

 sModel.D_EpsToEta   = D_EpsToEta;
 sModel.eta1BarCoeff = eta1BarCoeff;
 sModel.Phi01        = Phi01;
 sModel.Phi11        = Phi11;
 sModel.PhiEta1      = PhiEta1;
 sModel.Phi02        = Phi02;
 sModel.Phi12        = Phi12;
 sModel.PhiEta2      = PhiEta2;
 sModel.pSigmaTE     = pSigmaTE;
 sModel.pA           = pA;
 sModel.pA0          = pA0;