% This is a script that prior to estimation does the following
% 1) Select variables used in estimation
% 2) creates eval_zdata_script.m file 
% 3) creates eval_endo_names.m file containing steady-state values of variables
% 4) creates eval_param.m containing parameter values
% 5) processes constraint and creates eval_difference_script


%------------------------------------------------------------------
% 1) Select variables used in estimation
%------------------------------------------------------------------

eval(['dynare ',modnam_00,' nolog noclearall'])

wishlist_ = M_.endo_names;
nwishes_ = M_.endo_nbr;
[~, i1, ~]=intersect(wishlist_,obs_list,'rows');



%------------------------------------------------------------------
% 2) Create eval_zdata_script.m
%------------------------------------------------------------------

disp('create zdata script')
fid = fopen('eval_zdata_script.m','wt');
for i_indx_=i1'
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_l=zdatal(:,' num2str(i_indx_) ');\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_p=zdatap(:,' num2str(i_indx_) ');\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_ss=zdatass(' num2str(i_indx_) ');\n']);
end
fclose(fid);


%------------------------------------------------------------------
% 3) Create eval_endo_names.m
%------------------------------------------------------------------

fid = fopen('eval_endo_names.m','wt');
for i_indx_=1:M_.endo_nbr
  fprintf(fid,[deblank(M_.endo_names(i_indx_,:)),'_ss=oo00_.dr.ys(' num2str(i_indx_) ');\n']);
end
fclose(fid);


%------------------------------------------------------------------
% 4) Create eval_param.m
%------------------------------------------------------------------

fid = fopen('eval_param.m','wt');
for i_indx_=1:size(M_.param_names,1);
  fprintf(fid,[deblank(M_.param_names(i_indx_,:)),'=M00_.params(' num2str(i_indx_) ');\n']);
end
fclose(fid);




%---------------------------------------
% 5) Create all the matrices used by OccBin
% processes the constraint so as to uppend a suffix to each
%---------------------------------------
solve_two_constraints_firstcall(modnam_00,modnam_10,modnam_01,modnam_11);

[constraint1_difference iendo1]= process_constraint_with_tokens(constraint1,'_difference',M00_.endo_names,0);
[constraint_relax1_difference iendo2]= process_constraint_with_tokens(constraint_relax1,'_difference',M00_.endo_names,0);
[constraint2_difference iendo3]= process_constraint_with_tokens(constraint2,'_difference',M00_.endo_names,0);
[constraint_relax2_difference iendo4]= process_constraint_with_tokens(constraint_relax2,'_difference',M00_.endo_names,0);

iendo_constraint = union(union(union(iendo1,iendo2),iendo3),iendo4);
iendo_constraint = iendo_constraint(:);

fid = fopen('eval_difference_script.m','wt');
for i_indx_=iendo_constraint'
	fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_difference=zdatalinear_(:,' num2str(i_indx_) ');\n']);
end
fclose(fid);
