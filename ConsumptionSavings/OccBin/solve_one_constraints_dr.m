function [zdata zdataconcatenated ys_ init_out error_flag Ecurrent ] = solve_one_constraints_dr(...
    constraint1_difference,constraint_relax1_difference,...
    scalefactormod,irfshock,nperiods,curb_retrench,maxiter,init_orig)

  
  
global M00_ oo00_ 

global cof cof10  ...
       Jbarmat Jbarmat10  ...
       Dbarmat10  ...
       decrulea decruleb


 eval_endo_names

 eval_param




nvars = M00_.endo_nbr;
ys_ = oo00_.dr.ys;


endog_ = M00_.endo_names;
exog_ =  M00_.exo_names;


nshocks = size(scalefactormod,1);
if nargin<8 % default value if not passed in
    init_orig = zeros(nvars,1);
end
init = init_orig;
zdataconcatenated = zeros(nperiods,nvars);


regime1(1) =1;
regimestart1 =1;
violvecbool = zeros(nperiods+1,1);  
% This sets the first guess for when
% the constraints are going to hold.
% The variable is a boolean with two
% columns. The first column refers to
% constraint1; the second to
% constraint2.
% Each row is a period in time.
% If the boolean is true it indicates
% the relevant constraint is expected
% to evaluate to true.
% The default initial guess is
% consistent with the base model always
% holding -- equivalent to the linear
% solution.

wishlist = endog_;
nwishes = size(wishlist,1);


eval_endo_names
eval_param



for ishock = 1:nshocks
    
    
    changes=1;
    iter = 0;
    
    while (changes && iter<maxiter)
        
        iter = iter +1;
        % analyse violvec and isolate contiguous periods in the other
        % regime.
        [regime1, regimestart1]=map_regime(violvecbool);
        
        
        [zdatalinear_, Ecurrent ]=mkdatap_anticipated_1constraints_fast(nperiods,decrulea,decruleb,...
            cof,Jbarmat,...
            cof10,Jbarmat10,Dbarmat10,...
            regime1,regimestart1,...
            violvecbool,endog_,exog_,...
            irfshock,scalefactormod(ishock,:),init);
        
        
        eval_difference_script
        
        
        
        newviolvecbool = eval(constraint1_difference);
        relaxconstraint = eval(constraint_relax1_difference);
        
        
        
        
        
        
        % check if changes
        if (max(newviolvecbool(:)-violvecbool(:)>0)) || sum(relaxconstraint(find(violvecbool==1))>0)
            changes = 1;
        else
            changes = 0;
        end
        
        if curb_retrench   % apply Gauss-Sidel idea of slowing down the change in the guess
            % for the constraint -- only relax one
            % period at a time starting from the last
            % one when each of the constraints is true.
            retrench = 0*violvecbool(:);
            if ~isempty(find(relaxconstraint & violvecbool(:,1)))
                retrenchpos = max(find(relaxconstraint & violvecbool(:,1)));
                retrench(retrenchpos) = 1;
            end
            
            violvecbool = (violvecbool(:) | newviolvecbool(:))-retrench(:);
        else
            violvecbool = (violvecbool(:) | newviolvecbool(:))-(relaxconstraint(:) & violvecbool(:));
        end
        
        violvecbool = reshape(violvecbool,nperiods+1,1);
        
        
        
    end
    if changes ==1
        display('Did not converge')
        error_flag = 1;
    else 
        error_flag = 0;
    end
    
    init = zdatalinear_(1,:);
    zdataconcatenated(ishock,:)=init;
    init= init';
    
    % update the guess for constraint violations for next period
    % update is consistent with expecting no additional shocks next period
    violvecbool=[violvecbool(2:end,:);0];
    
    % overwrite guess for next period -- reset guess to expect no
    % constraint wil hold
    violvecbool=0*violvecbool;

end

init_out = init;
zdataconcatenated(ishock+1:end,:)=zdatalinear_(2:nperiods-ishock+1,:);

zdata = mkdata_fullendog(nperiods,decrulea,decruleb,endog_,exog_,irfshock,scalefactormod,init_orig);

