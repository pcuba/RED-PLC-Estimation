
function [theta_history, fval_history, accept_average]= ...
    mhalgo(posterior_func,mode,Hessian,target_average,init_draws,ndraws,load_switch,jump,varargin)

save_every_n_draws = 100;
filename = 'metropolis_chain';

P = chol(inv(Hessian));

% find c for jumping distribution -- open up or close the distribution to
% hit target acceptance rate.
options = optimset('display','iter','TolFun',1e-2,'TolX',1e-2,'MaxIter',1000,'MaxFunEvals',10000);
cinit= 1;
% c=fzero('rprob',cinit,options,posterior_func,P,mode,target_average,init_draws,varargin{:})
c=jump;

nparams = size(mode,1);

if load_switch
    load(filename,'fval_history','theta_history','c','accept_ratio','accept','indx', 'toss', 'toss_normal')
    s = RandStream('mt19937ar','Seed',indx);
    RandStream.setGlobalStream(s);
    theta_history = [theta_history(:,1:indx) zeros(nparams,ndraws)];
    fval_history = [fval_history(:,1:indx) zeros(1,ndraws)];
    accept_ratio = [accept_ratio(1:indx), zeros(1,ndraws)];
    fval = fval_history(:,indx);
    theta = theta_history(:,indx);
    
    toss= [toss(1:indx) rand(1,ndraws)];
    toss_normal=[toss_normal(:,1:indx) randn(nparams,ndraws)];
    indx_init = indx;
    
    
else
    if exist('metropolis','file')
        if isunix
            !rm metropolis.mat
        else
            !del metropolis.mat
        end
    end
    theta = mode;
    
    fval= feval(posterior_func,theta,varargin{:});
    
    theta_history = zeros(nparams,ndraws);
    fval_history = zeros(1,ndraws);
    accept_ratio = zeros(1,ndraws);
    
    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s);
    toss=rand(1,ndraws);
    toss_normal=randn(nparams,ndraws);
    indx_init = 0;
    accept = 0;
    
end


for indx = indx_init+1:indx_init+ndraws
    theta_new = theta + c^2*P'*toss_normal(:,indx);
    fval_new= feval(posterior_func,theta_new,varargin{:});
    
    % check if new draw is accepted or rejected
    if fval_new<=1e10;
        r = exp(-fval_new+fval)     % probability of acceptance
        if (toss(indx)<r)
            disp('accept')
            fval = fval_new;
            theta = theta_new;
            accept = accept+1;
            success(indx)=1;
        else
            disp('reject')
            success(indx)=0;
        end
        accept_ratio(indx)=accept/indx;
        disp([ 'Current ar = ' num2str(accept_ratio(indx),2) ])
        theta_history(:,indx)=theta;
        fval_history(indx)=fval;
    end
    
    if mod(indx,save_every_n_draws)==0
        save(filename,'fval_history','theta_history','c','accept_ratio','accept','indx', 'toss', 'toss_normal')
    end
    
end
save(filename,'fval_history','theta_history','c','accept_ratio','accept','indx', 'toss', 'toss_normal')
accept_average = accept/ndraws;




