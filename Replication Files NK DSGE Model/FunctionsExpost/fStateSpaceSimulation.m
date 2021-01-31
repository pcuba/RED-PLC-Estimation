
function [vState,vObs,vInd] = fStateSpaceSimulation(eta1BarCoeff,Phi01,Phi11,PhiEta1,Phi02,Phi12,PhiEta2,pA,pA0, vStatePrev, vEta)


eta1Bar = eta1BarCoeff(1) + eta1BarCoeff(2:end)*vStatePrev;


% Simulate states
if vEta(1) < eta1Bar(1)
    
    
    % Non-binding
   vState = Phi01 + Phi11*vStatePrev + PhiEta1*vEta;
   
   ind = 1;
else

    % Binding regime
   vState = Phi02 + Phi12*vStatePrev + PhiEta2*vEta;

   ind = 2;
end

% Simulate observables
vObs = pA0 + pA*vState;

% Regime indicator
vInd =  ind;





