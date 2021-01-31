function out = resid_forward_guidance(er_in,sModel,sParam,shocks,init,R_nointer,intervention_bound,objective)

% Written by S. Boragan Aruoba
% Chevy Chase, MD
% Created : August 6, 2012
% This version : October 16, 2012

nshocks   = length(er_in);
shocks_in = shocks;
shocks_in(1,1:nshocks) = er_in';

% temp = simulate_sunspot(theta_use,shocks_in,init,length(er_in),0,sunspot);
[~,~,~,temp] = fSimulatePLC(sModel,sParam,init,shocks_in);

if objective == 1       % Used for hitting the ZLB
    out = temp.R - 1;
elseif objective == 2   % Used for reducing R by the intervention bound
    out = 40000 * ( R_nointer - temp.R)  - intervention_bound; 
elseif objective == 3   % Used for minimizing the intervention R
    out = temp.R;
end

