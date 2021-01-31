function [zdata]=mkdata_fullendog(nperiods,decrulea,decruleb,endog_,exog_,irfshock,scalefactormod,init)




% given decision rule 
neqs = size(endog_,1);

if  nargin<8
   init = zeros(neqs,1);
end

if  nargin<7
    scalefactormod=1;
end

if nargin<6
    error('Not enough inputs')
end

history = zeros(neqs,nperiods+1);
exog_ = cellstr(exog_);
    nshocks = size(irfshock,1);
    for i = 1:nshocks
        shockpos = find(strcmp(irfshock(i,:),exog_));
        if ~isempty(shockpos)
            irfshockpos(i) = shockpos;
        else
            error(['Shock ',irfshock(i,:),' is not in the model']);
        end
    end


% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history = zeros(neqs,nperiods);
history(:,1)= init;

lengthshock = size(scalefactormod,1);

errvec = zeros(size(exog_,1),1);

for i = 2:nperiods+1
    if i<=(lengthshock+1)
        for j = 1:nshocks
            errvec(irfshockpos(j)) = scalefactormod(i-1,j);
        end
        history(:,i) = decrulea * history(:,i-1)+decruleb*errvec;
    else
    % update endogenous variables
    history(:,i) = decrulea * history(:,i-1);
    end
end


history=history';

zdata = history(2:end,:);