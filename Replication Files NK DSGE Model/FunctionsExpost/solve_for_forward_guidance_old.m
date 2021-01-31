function [shocks_use_inter,cases, monotone] = solve_for_forward_guidance_old(sModel,sParam,shocks_use,init_use,length_intervention_fg,...
    length_intervention_fiscal,intervention_bound)

% Function that solves for the er that sets the interest rate equal to zero
% for a pre-determined number of periods.

% Written by S. Boragan Aruoba
% Chevy Chase, MD
% Created : August 6, 2012
% This version : October 16, 2012

global options_fzero options_fmin size_spending_intervention

init = init_use;


solution = nan .* ones(length_intervention_fg,1);
cases    = nan .* ones(length_intervention_fg,1);
monotone = nan .* ones(length_intervention_fg,1);

go_on = 1;

t = 1;



while go_on == 1 && t <= length_intervention_fg   
    
    shocks = shocks_use(:,t);
    
    if t <= length_intervention_fiscal
        
        shocks(2,1) = shocks_use(2,t) + size_spending_intervention;
       
    end     
    
    % The no-intervention case (not to be confused with the no-intervention
    % case in the main script which doesn't have the g intervention),
    % taking as given the intervention up to this
    % period as well as the current-period fiscal intervention
    
    shocks_nointer      = shocks;
    shocks_nointer(1,1) = 0;  
    
    %temp = simulate_sunspot(theta_use,shocks_nofg,init,1,0,sunspot);
    [~,~,~,temp] = fSimulatePLC(sModel,sParam,init,shocks_nointer);

    
    R_nointer = temp.R;clear temp;  
    
    % The maximum intervention case, taking as given the intervention up to
    % this period, the maximum possible er intervention (on the negative
    % side) before the decision rule starts being non-monotonic, where the
    % lower-bound for the solution is given by e_min_use
    
    e_min_use = -2;
    
    shocks_maxinter = shocks;
    
%     [out_maxinter,~,eflag] = fminbnd(@(in) resid_forward_guidance(in,theta_use,...
%         shocks,init,sunspot,R_nofg,intervention_bound,3),e_min_use,0,options_fmin);

    [out_maxinter,~,eflag] = fminbnd(@(er_in) resid_forward_guidance(er_in,sModel,sParam,...
        shocks_nointer,init,R_nointer,intervention_bound,3),e_min_use,0,options_fmin);
    
    shocks_maxinter(1,1) = out_maxinter;
    %temp = simulate_sunspot(theta_use,shocks_maxinter,init,1,0,sunspot);
    [~,~,~,temp] = fSimulatePLC(sModel,sParam,init,shocks_maxinter);

    
    R_maxinter = temp.R;clear temp;  

        % Check monotonicity
    
    if out_maxinter > e_min_use + 0.0001
        monotone(t,1) = 0;
    else
        monotone(t,1) = 1;
    end
    
    if 40000 * (R_nointer - R_maxinter) < intervention_bound   % Here the maximum intervention is not larger than the intervention bound
        
        if R_maxinter == 1   % If maximum intervention leads to ZLB, solve for the smallest er that delivers that 
            
%             [out,~,eflag] = fzero(@(in) resid_forward_guidance(in,theta_use,shocks,init,sunspot,R_nofg,intervention_bound,1),[out_maxinter,0],options_fzero);

            [out,~,eflag] = fzero(@(er_in) resid_forward_guidance(er_in,sModel,sParam,...
    shocks_nointer,init,R_nointer,intervention_bound,1),[out_maxinter,0],options_fzero);
 
            cases(t,1) = 1;
            
        else    % If R_maxinter > 1, then the solution is simply the R_maxinter
            
            cases(t,1) = 3;
            
            out = out_maxinter;
            
        end
        
    else  % If the maximum intervention is too large, solve for the maximum intervention within the bound
        
%         [out,~,eflag] = fzero(@(in) resid_forward_guidance(in,theta_use,shocks,init,sunspot,R_nofg,intervention_bound,2),[out_maxinter,0],options_fzero);
                    [out,~,eflag] = fzero(@(er_in) resid_forward_guidance(er_in,sModel,sParam,...
    shocks_nointer,init,R_nointer,intervention_bound,2),[out_maxinter,0],options_fzero);

        cases(t,1) = 2;
        
    end
    
    if eflag < 1
        
        warning('Problem -- this should not happen!')
        go_on = 0;
        
    else
        
        shocks_in = shocks;
        shocks_in(1,1) = out;
        
        solution(1,t) = shocks_in(1,1);
        
        %temp = simulate_sunspot(theta_use,shocks_in,init,1,0,sunspot);
        [init_out,~,~,~] = fSimulatePLC(sModel,sParam,init,shocks_in);

        init = [init_out;init(5)];
        
        t = t + 1;
    end
    
    
        
        
    
end

shocks_use_inter = shocks_use;

shocks_use_inter(2,1:length_intervention_fiscal) = shocks_use(2,1:length_intervention_fiscal) + size_spending_intervention;
shocks_use_inter(1,1:length_intervention_fg) = solution(1,:); 





