%------------------------------------------
% Script that feeds filtered shocks back into model
% -- add observation block in model, check that final_samples matches obs
%------------------------------------------

init0=zeros(M00_.endo_nbr,1);


sample_length = size(obs,1);
load PARAM_EXTRA_CALIBRATED
for i=1:numel(params_labels)
  evalc([ cell2mat(params_labels(i)) '= params1(' num2str(i) ')']) ;
end
eval([ 'save PARAM_EXTRA_BABY ' cell2mat(params_labels') ]);




filtered_errs2=filtered_errs;
[amax imax]=min(obs(:,strmatch('ctot',obs_list)));

extra_t=0;


[zdatal zdatap zdatass oo00_ M00_ ] = solve_two_constraints_fast2_temp1(...
  modnam_00,modnam_10,modnam_01,modnam_11,...
  constraint1, constraint2,...
  constraint_relax1, constraint_relax2,...
  filtered_errs2,err_list,sample_length+extra_t,curb_retrench,maxiter,init0);
tt_obs2=[tt_obs
  tt_obs(end)+(1:extra_t)'/4];

for i=1:M00_.endo_nbr
  eval([deblank(M00_.endo_names(i,:)),'_l=zdatal(1:sample_length+extra_t,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_p=zdatap(1:sample_length+extra_t,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_ss=zdatass(i);']);
end


figure
final_sample = [];
%This is plotting the observable variable (in this case - consumption)
for index=1:size(obs_list,1)
  subplot(2,3,index)
  eval([ 'final_sample(:,index) = ' deblank(obs_list(index,:)) '_p;'])
  plot(tt_obs,obs(:,index),'bo-'); hold on
  plot(tt_obs2,eval([deblank(obs_list(index,:)) '_p']),'rs-','Linewidth',1)
  title(obs_list(index,:))
end
legend('Actual','Filtered')

%We also want to plot b, lb
load fakedata2
subplot(2,3,size(obs_list,1)+1)
plot(tt_obs,b_pp,'bo-'); hold on
plot(tt_obs2,b_p,'rs-','Linewidth',1)
title('b')

subplot(2,3,size(obs_list,1)+2)
plot(tt_obs,lb_pp,'bo-'); hold on
plot(tt_obs2,lb_p,'rs-','Linewidth',1)
title('lb')

for index=1:size(err_list,1)
  subplot(2,3,3+index)
%   load fakedata2 shocks_pp
  plot(tt_obs,shocks_pp(tt_obs,index),'bo-'); hold on
  plot(tt_obs2,filtered_errs(:,index),'rs-','Linewidth',1)
  title(err_list(index,:))
end


