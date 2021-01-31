% check_metropolis.m
% This is a function that displays output of metropolis-hastings for specific folders
% e.g. 
% check_metropolis({'C:\Dropbox\E\occbin_Estimation\borrcon',...
% 'C:\Dropbox\E\occbin_Estimation\borrcon_alt1',...
% 'C:\Dropbox\E\occbin_Estimation\borrcon_alt2'},'borrcon00',15)

function mmm = check_metropolis(foldername,modfilename,cutoff1)



nx=size(foldername,2); 

for ix=1:nx
    folder=foldername{ix}; 
    


eval([ 'load ' folder '\mle_estimates_fminsearch params1 params_labels IPRIOR'])

eval([ 'load ' folder '\metropolis_chain theta_history fval_history indx'])


% Find names of _results file
zzz=dir([folder '\*_results.mat']);
xx=zzz.name;
eval([ 'load ' folder '\' xx ' M_ ']);

eval([ 'load ' folder '\params_matrix params_matrix '])


figure(100*ix)
for ii=1:numel(params_labels)
    subplot(3,3,ii)
    [f,xi] = ksdensity(theta_history(ii,2:indx),'kernel','triangle','npoints',200); 
	xi(f<0.0001)=[];
    f(f<0.0001)=[];
    plot(xi,f); hold on
    
    params_labels = (params_matrix(:,1));
    
    plab = cell2mat(params_labels(ii,:));
    if plab(1:2)=='ST'; axis tight; end
    
    params0 =   cell2mat(params_matrix(:,2));
    params_lo = cell2mat(params_matrix(:,3));
    params_hi = cell2mat(params_matrix(:,4));
    params_mean = cell2mat(params_matrix(:,7));
    params_std = cell2mat(params_matrix(:,8));
    dist_names = params_matrix(ii,6);
    codes = dist_names2codes(params_matrix(ii,6));
    [p6 p7] = get_dist_inputs(codes,params_mean(ii),params_std(ii));
    hold on
    if IPRIOR==1
    myplot_priors2(codes,params_lo(ii),params_hi(ii),p6,p7,params_labels(ii,:))
    end
    hold on
    plot_lines_green(params1(ii));
    v = axis;
    title(cell2mat(params_labels(ii,:)))
    xlim([0.97*v(1) 1.03*v(2)])
end

if IPRIOR==1
legend('Posterior','Prior','mode')
else
legend('Posterior (Likelihood)','mode')
end
subplot(3,3,ii+1)
plot(-fval_history(1:indx))
title([ 'The posterior draws ' num2str(max(-fval_history),'%0.4f')])

figure(100*ix+1)
for ii=1:numel(params_labels)
    subplot(3,3,ii)
    plot(theta_history(ii,2:indx)); hold on; axis tight
    plab = cell2mat(params_labels(ii,:));
    if plab(1:2)=='ST'; axis tight; end    
    title(cell2mat(params_labels(ii,:)))
end



[ a b ]=max(-fval_history);

mode_metropolis = theta_history(:,b);
isprior=1;
if IPRIOR==0; isprior=NaN; end

disp('-----------')
disp(folder)
disp(['Number of runs is ' num2str(indx)])
disp(' ')
disp('PARAMETER   MODE MODE-METRO    PRIOR-DIS   PRIOR-MEAN PRIOR-SD POST-MEAN  10%      50%     90%')
for ii=1:numel(params_labels)
    trspaces=blanks(10-size(char(params_labels(ii,:)),2));
    trspaces2=blanks(13-size(char(params_matrix(ii,6)),2));
    disp([ char(params_labels(ii,:)) trspaces ' ' ...
        num2str(params1(ii),'%0.4f') '   '  ...
        num2str(mode_metropolis(ii),'%0.4f') '      '  ...
        char(params_matrix(ii,6)) trspaces2 '  ' ...
        num2str(isprior*cell2mat(params_matrix(ii,7)),'%0.4f')   '    ' ...
        num2str(isprior*cell2mat(params_matrix(ii,8)),'%0.4f')   '      ' ...
        num2str(mean(theta_history(ii,1:indx)),'%0.4f')   '   ' ...
        num2str(prctile(theta_history(ii,1:indx),10),'%0.4f')   '   ' ...
        num2str(prctile(theta_history(ii,1:indx),50),'%0.4f')  '   ' ...
                num2str(prctile(theta_history(ii,1:indx),90),'%0.4f')  ])
end
disp(' ')
disp('Calibrated parameters')
[ calibrated_names ical]=setdiff(M_.param_names,char(params_labels),'rows');
disp('-----------------')
for i=1:numel(ical)
    trspaces3=blanks(11-size(char(M_.param_names(ical(i),:)),2)) ;
    disp([  M_.param_names(ical(i),:)  '=' trspaces3  num2str(M_.params(ical(i),:),'%0.6f') ';'   ])
    i=i+1;
end

disp(' ')
disp(' ')
disp('Laplace approximation')
% disp('gtt: The variance covariance matrix of the estimates')
gtt=cov(theta_history(:,1:indx)');
% disp('gtt1: The correlation matrix of the estimates')
gtt1=corr(theta_history(:,1:indx)');
% See page 121 of http://mx1.texlips.net/download/yada.pdf

ndraws=numel(fval_history);

i=1;
for trunc=0.1:0.1:0.9
md(i,ix) = mhm_marginal_density(theta_history(:,(0.2*indx+1:indx)),fval_history(0.2*indx+1:indx),trunc);
i=i+1;
end
disp(['Marginal data density is ' num2str(median(md(:,ix)))])

end

