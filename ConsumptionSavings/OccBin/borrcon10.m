%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'borrcon10';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'eps_u';
M_.exo_names_tex = 'eps\_u';
M_.exo_names_long = 'eps_u';
M_.endo_names = 'b';
M_.endo_names_tex = 'b';
M_.endo_names_long = 'b';
M_.endo_names = char(M_.endo_names, 'bnot');
M_.endo_names_tex = char(M_.endo_names_tex, 'bnot');
M_.endo_names_long = char(M_.endo_names_long, 'bnot');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'ec');
M_.endo_names_tex = char(M_.endo_names_tex, 'ec');
M_.endo_names_long = char(M_.endo_names_long, 'ec');
M_.endo_names = char(M_.endo_names, 'lb');
M_.endo_names_tex = char(M_.endo_names_tex, 'lb');
M_.endo_names_long = char(M_.endo_names_long, 'lb');
M_.endo_names = char(M_.endo_names, 'maxlev');
M_.endo_names_tex = char(M_.endo_names_tex, 'maxlev');
M_.endo_names_long = char(M_.endo_names_long, 'maxlev');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_partitions = struct();
M_.param_names = 'RHO';
M_.param_names_tex = 'RHO';
M_.param_names_long = 'RHO';
M_.param_names = char(M_.param_names, 'BETA');
M_.param_names_tex = char(M_.param_names_tex, 'BETA');
M_.param_names_long = char(M_.param_names_long, 'BETA');
M_.param_names = char(M_.param_names, 'M');
M_.param_names_tex = char(M_.param_names_tex, 'M');
M_.param_names_long = char(M_.param_names_long, 'M');
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names_long = char(M_.param_names_long, 'R');
M_.param_names = char(M_.param_names, 'STD_U');
M_.param_names_tex = char(M_.param_names_tex, 'STD\_U');
M_.param_names_long = char(M_.param_names_long, 'STD_U');
M_.param_names = char(M_.param_names, 'GAMMAC');
M_.param_names_tex = char(M_.param_names_tex, 'GAMMAC');
M_.param_names_long = char(M_.param_names_long, 'GAMMAC');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 7;
M_.param_nbr = 6;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('borrcon10_static');
erase_compiled_function('borrcon10_dynamic');
M_.orig_eq_nbr = 7;
M_.eq_nbr = 7;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 3 0;
 0 4 0;
 0 5 10;
 0 6 0;
 0 7 0;
 0 8 0;
 2 9 0;]';
M_.nstatic = 4;
M_.nfwrd   = 1;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.NNZDerivatives = [18; 4; -1];
save('borrcon10_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('borrcon10_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('borrcon10_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('borrcon10_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('borrcon10_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('borrcon10_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('borrcon10_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
