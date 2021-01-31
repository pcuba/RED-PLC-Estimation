function [resid] = rprob(c,posterior_func,P,mode,target_average,init_draws,varargin)



% draw a candidate from a jumping distribution; use N(theta0,c^2 sigmaTilde)
% to draw from jumping distribution, use Cholesky approach:
% y = mu + P' e     where e ~ N(0,I) and P'P = c^2 sigmaTilde 
randn('seed',1)
rand('seed',1)
theta = mode;
nparams = size(theta,1);
fval= feval(posterior_func,theta,varargin{:});

accept = 0;
for indx =2:init_draws
    theta_new = theta + c^2*P'*randn(nparams,1);
    fval_new= feval(posterior_func,theta_new,varargin{:});
     
    % check if new draw is accepted or rejected
    if fval_new<=1e10;
        r= min(exp(-fval_new+fval),1);     % probability of acceptance
        r= max(r,0);
        if r<0
            fprintf('r<0')
            break
        end
        toss = rand;
        if (toss<r) 
            fval = fval_new;
            theta = theta_new;
            accept = accept+1;
        end
        
    end
end
accept_average = accept/init_draws;

resid = accept_average-target_average;
