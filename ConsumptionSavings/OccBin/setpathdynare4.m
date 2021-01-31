% Adapt to local installation
% We recommend to use Dynare 4.3.1
% addpath C:\dynare\4.3.1\matlab\
%
% For Dynare >4.5 you need to modify estobin script solve_one_constraint
% To make sure you pass the number of static and forward looking variables
% For example:
% oo00_.dr.nstatic = M00_.nstatic;
% oo00_.dr.nfwrd = M00_.nfwrd;
%======================================================================

dir1='/Applications/Dynare/4.5.7/matlab';
dir2='/Users/pcb/Documents/MATLAB/occbin_20130531/toolkit_files';

path(dir1,path);
path(dir2,path);


dynare_config

