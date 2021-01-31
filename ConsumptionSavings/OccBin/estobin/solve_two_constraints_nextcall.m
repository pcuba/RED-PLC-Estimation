function [zdata zdataconcatenated ys_ init_out error_flag Ecurrent newviolvecbool relaxconstraint iter] = solve_two_constraints_nextcall(...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    scalefactormod,irfshock,nperiods,curb_retrench,maxiter,init_orig)

  
  
global M00_ oo00_ 

global cof cof10 cof01 cof11 ...
       Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 ...
       Dbarmat10 Dbarmat01 Dbarmat11 ...
       decrulea decruleb


 eval_endo_names

 eval_param





nvars = M00_.endo_nbr;
ys_ = oo00_.dr.ys;


endog_ = M00_.endo_names;
exog_ =  M00_.exo_names;


nshocks = size(scalefactormod,1);
if nargin<10 % default value if not passed in
    init_orig = zeros(nvars,1);
end
init = init_orig;
zdataconcatenated = zeros(nperiods,nvars);


regime1(1) = 1;
regime2(1) = 1;
regimestart1 = 1;
regimestart2 = 1;
violvecbool = zeros(nperiods+1,2);  % This sets the first guess for when
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
    
    while (changes & iter<maxiter)
        
        iter = iter +1;
        
        % analyse violvec and isolate contiguous periods in the other
        % regime.
        [regime1 regimestart1]=map_regime(violvecbool(:,1));
        [regime2 regimestart2]=map_regime(violvecbool(:,2));
        
        
        [zdatalinear_ Ecurrent ]=mkdatap_anticipated_2constraints_fast(nperiods,decrulea,decruleb,...
            cof,Jbarmat,...
            cof10,Jbarmat10,Dbarmat10,...
            cof01,Jbarmat01,Dbarmat01,...
            cof11,Jbarmat11,Dbarmat11,...
            regime1,regimestart1,...
            regime2,regimestart2,...
            violvecbool,endog_,exog_,...
            irfshock,scalefactormod(ishock,:),init);
        
        
        eval_difference_script
        
        
        
        newviolvecbool1 = eval(constraint1_difference);
        relaxconstraint1 = eval(constraint_relax1_difference);
        
        newviolvecbool2 = eval(constraint2_difference);
        relaxconstraint2 = eval(constraint_relax2_difference);
                
        
        newviolvecbool = [newviolvecbool1;newviolvecbool2];
        relaxconstraint = [relaxconstraint1;relaxconstraint2];
        
        
        
        % check if changes
        if (max(newviolvecbool(:)-violvecbool(:)>0)) | sum(relaxconstraint(find(violvecbool==1))>0)
            changes = 1;
        else
            changes = 0;
        end
        
%         iter
%          keyboard
        
        if curb_retrench   % apply Gauss-Sidel idea of slowing down the change in the guess
            % for the constraint -- only relax one
            % period at a time starting from the last
            % one when each of the constraints is true.
            retrench = 0*violvecbool(:);
            if ~isempty(find(relaxconstraint1 & violvecbool(:,1)))
                retrenchpos = max(find(relaxconstraint1 & violvecbool(:,1)));
                retrench(retrenchpos) = 1;
            end
            if ~isempty(find(relaxconstraint2 & violvecbool(:,2)))
                retrenchpos = max(find(relaxconstraint2 & violvecbool(:,2)));
                retrench(retrenchpos+nperiods+1) = 1;
            end
            violvecbool = (violvecbool(:) | newviolvecbool(:))-retrench(:);
        else
            violvecbool = (violvecbool(:) | newviolvecbool(:))-(relaxconstraint(:) & violvecbool(:));
        end
        
        violvecbool = reshape(violvecbool,nperiods+1,2);
        
        
        
    end
    if changes ==1
%         display('Did not converge')
        error_flag = 1;
    else 
        error_flag = 0;
    end
    
    init = zdatalinear_(1,:);
    zdataconcatenated(ishock,:)=init;
    init= init';
    
    % update the guess for constraint violations for next period
    % update is consistent with expecting no additional shocks next period
    violvecbool=[violvecbool(2:end,:);zeros(1,2)];
    
    % overwrite guess for next period -- reset guess to expect no
    % constraint wil hold
    violvecbool=0*violvecbool;
    
end

init_out = init;
zdataconcatenated(ishock+1:end,:)=zdatalinear_(2:nperiods-ishock+1,:);

zdata = mkdata_fullendog(nperiods,decrulea,decruleb,endog_,exog_,irfshock,scalefactormod,init_orig);



