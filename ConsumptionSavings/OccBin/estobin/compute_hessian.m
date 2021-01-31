function [ hessian_reg stdh_reg hessian_fmin stdh_fmin ] = ...
    compute_hessian(xstory,fstory,cutoff);

likstory=-fstory(:,1);
    
[yy]=likstory(likstory>max(likstory)-cutoff);
[xx]=xstory(find(likstory>max(likstory)-cutoff),:);
[cc]=yy.^0;

n = size(xstory,2);

[B]=regress(yy,[cc demean(xx) demean(xx).^2/2 ]);
diag_hessian=-B(end-n+1:end);
hessian_reg=diag(diag_hessian);

stdh_reg = 1./abs(sqrt(diag_hessian));

hessian_fmin = inv(cov(xx));

stdh_fmin = diag(inv(hessian_fmin)).^0.5;

