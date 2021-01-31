

% The vector of states is:
% s_t = [er_t, g_t, z_t, d_t, y_t, pi_t, c_t, R_t, y_lag]
% The order of the shocks is:
% vEps: er, eg, ez, ed

function [vState,vObs,vInd,sSIM] = fSimulatePLC_2(sModel,sParam,vStatePrev, vEps)

% Map canonical from matrices
D_EpsToEta   = sModel.D_EpsToEta;
eta1BarCoeff = sModel.eta1BarCoeff;
Phi01        = sModel.Phi01;
Phi11        = sModel.Phi11;
PhiEta1      = sModel.PhiEta1;
Phi02        = sModel.Phi02;
Phi12        = sModel.Phi12;
PhiEta2      = sModel.PhiEta2;
pSigmaTE     = sModel.pSigmaTE;
pA           = sModel.pA;
pA0          = sModel.pA0;

y_lag = vStatePrev(9,1);

% Number of periods for the simulation
[~, nT] = size(vEps);
[nS,~]  = size(vStatePrev(1:8));
[nY,~]  = size(pA);

% Initialize matrices
vState = zeros(nS,nT);
vObs   = zeros(nY,nT);
vInd   = zeros(1,nT);
As     = zeros(1,nT);
CCs    = zeros(1,nT);
YYs    = zeros(1,nT);
GGs    = zeros(1,nT);
delys = zeros(1,nT);

for tt=1:nT

    
    eta1Bar = eta1BarCoeff(1) + eta1BarCoeff(2:end)*vStatePrev(1:8);

    vEta = D_EpsToEta*vEps(:,tt);

% Simulate states in period tt
% if vEta(1) < eta1Bar(1)
%         
%     % Non-binding
%    vState(:,tt) = Phi01 + Phi11*vStatePrev(1:8) + PhiEta1*vEta;
%    
%    ind = 1;
% else

    % Binding regime
   vState(:,tt) = Phi12*vStatePrev(1:8) + PhiEta2*vEta;

   ind = 2;
% end


vStateAugmented = [vState(:,tt);vStatePrev(5,1)];

% Simulate observables in period tt
vObs(:,tt) = pA0 + pA*vStateAugmented;

% Regime indicator in period tt
vInd(:,tt) =  ind;


% Construct level variables %
g_t = vState(2,tt);
z_t = vState(3,tt);
y_t = sParam.y_ss*exp(vState(5,tt));
c_t = sParam.c_ss*exp(vState(7,tt));

if tt==1
    A_lag    = 1;
    As(tt)   = sParam.gamma*A_lag*exp(z_t);
    delys(tt) = 100*log(sParam.gamma) + 100*(vState(5,tt) - y_lag + vState(3,tt));
else
    As(tt)   = sParam.gamma*As(tt-1)*exp(z_t);    
    delys(tt) = 100*log(sParam.gamma) + 100*(vState(5,tt) - vState(5,tt-1) + vState(3,tt));
end

if As(tt) < 1E11
    
    CCs(tt) = c_t*As(tt);
    YYs(tt) = y_t*As(tt);
    GGs(tt) = (1 - ( 1/ (sParam.gstar*exp(g_t)) )) * YYs(tt);
    
else
    
    CCs(tt)= NaN;
    YYs(tt)= NaN;
    GGs(tt)= NaN;
end



% Update States
vStatePrev = [vState(:,tt);vStatePrev(5)];

end


% Collect simulations in levels
% The vector of states is:
% s_t = [er_t, g_t, z_t, d_t, y_t, pi_t, c_t, R_t, y_lag]

%sSIM.er   = vState(1,:);
sSIM.g    = sParam.gstar*exp(vState(2,:));
sSIM.z    = exp(vState(3,:));
sSIM.d    = exp(vState(4,:));
sSIM.y    = (vState(5,:));
sSIM.pi   = (vState(6,:));
sSIM.c    = (vState(7,:));
sSIM.R    = (vState(8,:));
sSIM.BIGC = CCs;
sSIM.BIGY = YYs;
sSIM.BIGG = GGs;
sSIM.dely = delys;

% vEps: er, eg, ez, ed
sSIM.er   = vEps(1,:);
sSIM.eg   = vEps(2,:);
sSIM.ez   = vEps(3,:);
sSIM.ed   = vEps(4,:);



