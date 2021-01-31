function [filtered_errs resids Emat stateval ] = myfilter(constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,solver)

global M00_

global filtered_errs_switch filtered_errs_init this_period sample_length obs_temp model_temp


%-------------------------------------
% Filter shocks
%-------------------------------------


sample_length = size(obs,1);
nerrs = size(err_list,1);
init_val = zeros(M00_.endo_nbr,1);
err_vals = zeros(nerrs,1);

resids = zeros(sample_length,nerrs);
stateval = zeros(sample_length,M00_.endo_nbr);


my_list = 1:nerrs;

maxiters = 20;
tolresidr = 1e-5;
tolsolve = 1e-10; 


options_fsolve = optimset('Display','None','MaxFunEvals',1e10,'MaxIter',1e5,'TolFun',1e-4,...
    'Algorithm','trust-region-dogleg');

for this_period=1:sample_length;
    
    
    current_obs = obs(this_period,:);
    init_val_old = init_val;
    %disp(this_period);
    
        
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        if solver ==0
        
         [ err_vals_out em ] = csolve(@(err_vals) match_function(...
            err_vals,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference),...
            err0',[],tolsolve,maxiters);

        elseif solver==1
            
        [ err_vals_out em ] = csolve_grad('match_function',...
            err0',tolsolve,maxiters,...
            err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        elseif solver==2

        [ err_vals_out_zlb ] = fsolve(@(err_vals) match_function(...
         err_vals,err_list,obs_list,current_obs,init_val,...
         constraint1_difference,constraint2_difference,...
         constraint_relax1_difference,constraint_relax2_difference),err0',options_fsolve);

        elseif solver==3

         err_vals_out = fzero(@(err_vals) match_function(...
         err_vals,err_list,obs_list,current_obs,init_val,...
         constraint1_difference,constraint2_difference,...
         constraint_relax1_difference,constraint_relax2_difference), err0');
        
        end
        
        filtered_errs(this_period,:)=err_vals_out';
        
        
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period),~,~,~,stateval(this_period,:)] = match_function(...
            err_vals_out,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
%         disp([this_period em resids(this_period,:)])
        
        if abs(resids(this_period))>0.001
            disp('I am stopping because match_function could not find the shocks that')
            disp('solve for the model`s observed variables')
            disp('')
            disp('Call test_error.m to do more debugging')
            keyboard
        end







    
    if max(abs(resids(this_period,:)))>0.05
        init_val_old=0*init_val_old;
        error('huge resids, give up')
    end
    

    
    
end




end


