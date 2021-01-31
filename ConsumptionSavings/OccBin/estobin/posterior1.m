function [posterior filtered_errs like prior resids stateval ] =posterior1(current_params,params_labels,lowerbound,upperbound,...
    modnam_00_,modnam_10_,modnam_01_,modnam_11_,...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    err_list,obs_list,obs,ntrain,codes, p6, p7, IPRIOR,solver)



global M_ oo_ oo00_  M00_ M10_  M01_  M11_
global cof cof10 cof01 cof11 Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 Dbarmat10 Dbarmat01 Dbarmat11 decrulea decruleb
global filtered_errs_init 
global datavec irep xstory fstory

if size(current_params,1)<size(current_params,2)
    current_params=current_params';
    params=current_params;
else
    params = min(max(lowerbound,current_params),upperbound);
end


% save to disk -- will be read by the parameter file invoked
% when the model is solved again.
for i=1:numel(params_labels)
    x=cell2mat(params_labels(i));
    evalc([ x '= params(' num2str(i) ')']) ;
end
eval([ 'save PARAM_EXTRA_BABY ' cell2mat(params_labels') ]);


eval(modnam_00_)
oo00_ = oo_;
M00_ = M_;



zdatass = oo00_.dr.ys;

[hm1,h,hl1,Jbarmat] = get_deriv(M00_,zdatass);
cof = [hm1,h,hl1];


M10_.params = M00_.params;
[hm1,h,hl1,Jbarmat10,resid] = get_deriv(M10_,zdatass);
cof10 = [hm1,h,hl1];
Dbarmat10 = resid;

if isempty(modnam_01_)==0
    
    M01_.params = M00_.params;
    [hm1,h,hl1,Jbarmat01,resid] = get_deriv(M01_,zdatass);
    cof01 = [hm1,h,hl1];
    Dbarmat01 = resid;
    
    M11_.params = M00_.params;
    [hm1,h,hl1,Jbarmat11,resid] = get_deriv(M11_,zdatass);
    cof11 = [hm1,h,hl1];
    Dbarmat11 = resid;
    
else
    
    cof01=[];
    cof11=[];
    
    Dbarmat01 = [];
    Dbarmat11 = [];
    
    Jbarmat01=[];
    Jbarmat11=[];
    
end

[decrulea,decruleb]=get_pq(oo00_.dr);


%---------------------------------------------
% Calculate likelihood
%---------------------------------------------


sample_length = size(obs,1);
nerrs = size(err_list,1);

filtered_errs_init = zeros(sample_length,nerrs);

[filtered_errs resids Emat stateval] = myfilter(constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,solver);
nobs=size(filtered_errs,1);


[~, ~, ishocksfe ]=intersect(err_list,err_list,'rows');


% If shocks are among the params_labels, estimate them
if numel(findstr('STD',cell2mat(params_labels(:)')))==nerrs
    
    
    disp('Will estimate shocks without concentrating')
    for i=1:size(err_list,1)
        eval( [ 'COVMAT1(i,i) = STD_' upper(err_list(i,5)) '^2;'] )
    end
    
else
    
    % Otherwise concentrate likelihood
    disp('Will concentrate likelihood')
    for i=1:size(err_list,1)
        COVMAT1(i,i) = std(filtered_errs(:,i))^2 ;
        evalc([ 'STD_' upper(err_list(i,5)) '= COVMAT1(i,i)^0.5']) ;
    end
    COVMAT0=COVMAT1(ishocksfe_zlb,ishocksfe_zlb');
    
end



%-------------------------------------
% Calculate the selection matrix
%-------------------------------------

selector_matrix1=zeros(size(obs_list,1),M_.endo_nbr);
for iobs=1:size(obs_list,1)
    [~, ~, iobscols1]=intersect(obs_list,M_.endo_names,'rows');
    selector_matrix1(iobs,iobscols1(iobs))=1;
end

[~, ~, ishocks ]=intersect(err_list,M_.exo_names,'rows');




likeall=0;

for t = 1:nobs
    
    Gmat1  = selector_matrix1*Emat(:,ishocks,t);
    log_det_jacobian(t) = log(det(COVMAT1)) + 2*log(abs(det(Gmat1)));
    trace_term(t) = filtered_errs(t,ishocksfe)*inv(COVMAT1)*filtered_errs(t,ishocksfe)';
    
    likei(t,1) = log_det_jacobian(t)/2 + trace_term(t)/2;
    likeall = likeall + likei(t);
    
end

like = sum(likei(ntrain+1:end));


if max(abs(params-current_params))>1e-8
    disp('Penalize params outside bound')
    like = like + 1e6*max(abs(params-current_params)) ;
end

maxresid = max(abs(resids(:)));
if maxresid>1e-3
    disp('Penalize failure of residuals to be zero')
    like = like + sum(resids(:).^2)*1e7;
end


[prior] = -priordens(params, codes, p6, p7, lowerbound, upperbound,1);
if prior == Inf
    % If parameters outside prior bound, minus prior density is very large
    prior= 1e8;
end



if isinf(like)==1
    like=1e8;
end


%---------------------------------------------
% Calculate posterior
%---------------------------------------------

% remember that the likelihood has already been multiplied by -1
% hence, posterior is -1 times the log of the prior
posterior = like+IPRIOR*prior;



%-------------------------------------
% Display info on screen
%-------------------------------------

disp(['Current minus posterior is ' num2str(-posterior) ' , iteration # ' num2str(irep) ])
datavec(irep,:) = [ posterior ];
fstory(irep,:) = [ posterior like prior ];
xstory(irep,:) = [ current_params ];

% PCB Commented This: 
% % save datavec datavec fstory xstory params* IPRIOR
% % 
% % if datavec(irep,1)==min(datavec(:,1))
% %     disp('minimum found, save into mle_estimates_temp_test')
% %     params1=params;
% %     save mle_estimates_temp_test like posterior prior datavec filtered_errs err_list obs_list obs params*
% % end

irep = irep+1;

disp(' ')