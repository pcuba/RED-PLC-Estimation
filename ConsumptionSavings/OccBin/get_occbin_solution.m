function [nvars,ys_,endog_,exog_,params, decrulea,decruleb,cof,...
    Jbarmat,cof10,Jbarmat10,Dbarmat10] = get_occbin_solution(modnam_00_,modnam_10_,solver,paramvec_)


global M_ oo_ oo00_  M00_ M10_  options_

eval_param;
eval(modnam_00_)
oo00_ = oo_;
M00_ = M_;

% 
% % Pass parameters
% M_.params(1) = paramvec_(1);
% M_.params(2) = paramvec_(2);
% M_.params(3) = paramvec_(3);
% M_.params(4) = paramvec_(4);
% M_.params(5) = paramvec_(5);
% M_.params(6) = paramvec_(6);
% 
% % Resolve Reference Model at New Parameter Vector
% check_flag = 0;
% try
% [~,~,M00_,~,oo00_] = resol(check_flag,M_,options_,oo_);
% catch
%     keyboard;
% end

% FOR COMPATIBILITY WITH DYNARE 4.5 or higher
oo00_.dr.nstatic = M00_.nstatic;
oo00_.dr.nfwrd = M00_.nfwrd;


zdatass = oo00_.dr.ys;

[hm1,h,hl1,Jbarmat] = get_deriv(M00_,zdatass);
cof = [hm1,h,hl1];


M10_.params = M00_.params;
[hm1,h,hl1,Jbarmat10,resid] = get_deriv(M10_,zdatass);
cof10 = [hm1,h,hl1];
Dbarmat10 = resid;


[decrulea,decruleb]=get_pq(oo00_.dr);


nvars          = M00_.endo_nbr;
ys_            = oo00_.dr.ys;
endog_         = M00_.endo_names;
exog_          = M00_.exo_names;


% Evaluate model parameters
% eval_param
params(1)=M00_.params(1);   % RHO
params(2)=M00_.params(2);   % BETA
params(3)=M00_.params(3);   % M
params(4)=M00_.params(4);   % R
params(5)=M00_.params(5);   % STD_U
params(6)=M00_.params(6);   % GAMMAC



