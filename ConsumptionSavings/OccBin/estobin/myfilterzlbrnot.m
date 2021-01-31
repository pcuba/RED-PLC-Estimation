function [filtered_errs resids Emat requalzero ] = myfilterzlbrnot(constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,rzlb,GG,PP,QQ)

global M00_

global filtered_errs_switch filtered_errs_init this_period sample_length obs_temp model_temp

warning off
verbose = 0;

%-------------------------------------
% Filter shocks
%-------------------------------------

options_fsolve = optimset('Display','None','MaxFunEvals',1e10,'MaxIter',1e5,'TolFun',1e-4,...
    'Algorithm','trust-region-dogleg');

sample_length = size(obs,1);
nerrs = size(err_list,1);
init_val = zeros(M00_.endo_nbr,1);
err_vals = zeros(nerrs,1);
err_vals_zlb = zeros(nerrs-1,1);

resids = zeros(sample_length,nerrs);

pos_r = find(strcmp('r',cellstr(obs_list)));
pos_eps_r = find(strcmp('eps_r',cellstr(err_list)));

my_list = 1:nerrs;

err_list_zlb = setdiff(err_list,'eps_r','rows');
obs_list_zlb = setdiff(obs_list,'r','rows');
obs_list_rnot = union(obs_list_zlb,'rnot','rows');


[~, i_obs_list, ~]=intersect(M00_.endo_names,obs_list,'rows');
[~, i_obs_list_rnot, ~]=intersect(M00_.endo_names,obs_list_rnot,'rows');

[~, i1_obs_list, ~]=intersect(M00_.endo_names,obs_list,'rows');
[~, i1_obs_list_zlb, ~]=intersect(M00_.endo_names,obs_list_zlb,'rows');
[~, i1_obs_list_rnot, ~]=intersect(M00_.endo_names,obs_list_rnot,'rows');
[~, i2_err_list, ~]=intersect(M00_.exo_names,err_list,'rows');
[~, i2_err_list_zlb, ~]=intersect(M00_.exo_names,err_list_zlb,'rows');

obs_list_withrnot = union(obs_list,'rnot','rows');
obs_list_nor = setdiff(obs_list_withrnot,{'rnot';'r'},'rows');
obs_list_r = intersect(obs_list_withrnot,{'rnot';'data_r'},'rows');
[~, ii_nor, ~]=intersect(M00_.endo_names,obs_list_nor,'rows');
[~, ii_r, ~]=intersect(M00_.endo_names,obs_list_r,'rows');



maxiters = 8;
requalzero=zeros(1,sample_length);
tolresidr = 1e-6;
tolsolve = 1e-6;
tolzlb = 1e-6;


for this_period=1:sample_length;
    
    if verbose
        disp(' ')
        disp(' ')
        disp(' ')
        disp(['Period number ' num2str(this_period)])
    end
    
    current_obs = obs(this_period,:);
    init_val_old = init_val;
    
    % STEP 1: WE ARE AT THE ZLB AND DROP R AS OBSERVABLE
    % TRYING TO SOLVE FOR TWO SHOCKS ONLY
    if current_obs(pos_r)<rzlb+tolzlb
        
        if verbose
            disp(' ')
            disp('At ZLB, try csolve')
        end
        
        requalzero(this_period)=1;
        
        current_obs_zlb = current_obs(find(my_list~=pos_r));
        err0 = filtered_errs_init(this_period,1:numel(err_vals_zlb));
        
        
        [ err_vals_out_zlb em ]= csolve_grad('match_function_zlb',...
            err0',tolsolve,maxiters,...
            err_list_zlb,obs_list_zlb,current_obs_zlb,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_zlb,i2_err_list_zlb,i1_obs_list_zlb);
        
        
        err_vals_out = zeros(nerrs,1);
        err_vals_out(find(my_list~=pos_eps_r)) = err_vals_out_zlb ;
        filtered_errs(this_period,:)=err_vals_out';
        
        [resids(this_period,:), ~, init_val, Emat(:,:,this_period) ] = match_function_zlb(...
            err_vals_out,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_zlb,i2_err_list_zlb,i_obs_list);
        
        
    end
    
    
    % STEP 2: TRIED STEP 1 @ZLB. IF RESIDUALS ARE STILL LARGE, 
    
    if max(abs(resids(this_period,:)))>tolresidr && current_obs(pos_r)<rzlb+tolzlb
        
        if verbose==1
            disp(' ')
            disp('After trying csolve at ZLB, got big residuals')
            disp([ this_period NaN 100*resids(this_period,:)])
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals_zlb));
        
        [ err_vals_out_zlb ] = fsolve(@(err_vals) match_function_zlb(...
            err_vals,err_list_zlb,obs_list_zlb,current_obs_zlb,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_zlb,i2_err_list_zlb,i1_obs_list_zlb), err0',options_fsolve);
        
        err_vals_out = zeros(nerrs,1);
        err_vals_out(find(my_list~=pos_eps_r)) = err_vals_out_zlb ;
        filtered_errs(this_period,:)=err_vals_out';
        
        [resids(this_period,:), ~, init_val, Emat(:,:,this_period) ] = ...
            match_function_zlb(...
            err_vals_out,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list,i2_err_list,i_obs_list);
        
        if verbose==1
            disp(' ')
            disp('Call fsolve... done')
            disp([ this_period NaN 100*resids(this_period,:)])
            disp('Called fsolve at ZLB, compare residuals before and after ')
        end
        
    end
    
    
    
    if current_obs(pos_r)>=rzlb+tolzlb
        
        if verbose
            disp(' ')
            disp('We are not at ZLB, use fast linear filter')
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        % Exploit
        % X(t) = P*X(t-1)+Q*e(t)
        % Y(t) = G*X(t)
        % e(t) = (GG*QQ)^(-1)*(Y(t)-GG*PP*X(t-1))
        
        err_vals_out = inv(GG*QQ)*(current_obs'-GG*PP*init_val_old);
        
        %                 [ err_vals_out em ] = csolve_grad('match_function_zlb',...
        %                     err0',tolsolve,maxiters,...
        %                     err_list,obs_list_rnot,current_obs,init_val_old,...
        %                     constraint1_difference,constraint2_difference,...
        %                     constraint_relax1_difference,constraint_relax2_difference,...
        %                     i1_obs_list,i2_err_list,i_obs_list_rnot);
         
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = ...
            match_function_zlb(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list,i2_err_list,i_obs_list_rnot);
        
    end
    
    
    
    if (resids(this_period,pos_r))>tolresidr
        
        if verbose==1
            disp(' ')
            disp('I tried fsolve, but residuals are positive for interest rate only')
            disp([ this_period NaN resids(this_period,:)])
            disp('Calling csolve again, this time allowing for monetary shocks too?...')
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        [ err_vals_out em ] = csolve_grad('match_function_zlb',...
            err0',tolsolve,maxiters,...
            err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_rnot,i2_err_list,i_obs_list_rnot);
        
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = match_function_zlb(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_rnot,i2_err_list,i_obs_list_rnot);
        
        if verbose==1
            disp(' ')
            disp([ this_period NaN resids(this_period,:)])
            disp('Added monetary shocks because notional rate was above zero')
            disp(' ')
        end
        
    end
    
    
    
    if max(abs(resids(this_period,:)))>tolresidr && current_obs(pos_r)>=rzlb+tolzlb
    
        if verbose==1
            disp(' ')
            disp('Residuals still nonzero and interest rate above zero')
            disp([ this_period NaN resids(this_period,:)])
            disp(' ')
        end
            
        
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        [ err_vals_out ] = fsolve(@(err_vals) match_function_zlb(...
            err_vals,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_rnot,i2_err_list,i_obs_list_rnot),...
            err0',options_fsolve);
        
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = match_function_zlb(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference,...
            i1_obs_list_rnot,i2_err_list,i_obs_list_rnot);
        
        if verbose==1
            disp('Calling fsolve...')
            disp([ this_period NaN 100*resids(this_period,:)])
            disp('I just called fsolve, compare residuals before and after ')
            disp(' ')
        end
        
    end
    
    
    
    
    
    
    
    if max(abs(resids(this_period,:)))>0.05
        init_val_old=0*init_val_old;
        error('huge resids, give up')
    end
    
    
    if verbose
        elapsed_time(this_period,:) = toc;
    end
    
    
    
end




end