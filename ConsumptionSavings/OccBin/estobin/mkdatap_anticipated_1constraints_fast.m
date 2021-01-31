function [zdata Ecurrent ]=mkdatap_anticipated_1constraints_fast(nperiods,decrulea,decruleb,...
    cof,Jbarmat,...
    cof10,Jbarmat10,Dbarmat10,...
    regime1,regimestart1,...
    violvecbool,endog_,exog_,...
    irfshock,scalefactormod,init)


nvars = size(endog_,1);


if nargin<16
    init=zeros(nvars,1);
end

if nargin<15;
    scalefactormod=1;
end


nshocks = size(irfshock,1);

exog_cell = cellstr(exog_);

for i = 1:nshocks
    % Matteo: shock names are different length, you need to you strcmp,
    % not find(strcmp)
    shockpos = (strmatch(irfshock(i,:),exog_cell));
    
    %shockpos = (strmatch(irfshock(i,:),exog_cell)); %erase potentialle again
    
    if ~isempty(shockpos)
        irfshockpos(i) = shockpos;
    else
        error(['Shock ',irfshock(i,:),' is not in the model']);
    end
end



Cbarmat = cof(:,1:nvars);
Bbarmat = cof(:,nvars+1:2*nvars);
Abarmat = cof(:,2*nvars+1:3*nvars);


% cofstar contains the system for the model when the constraint binds


Cbarmat10 = cof10(:,1:nvars);
Bbarmat10 = cof10(:,nvars+1:2*nvars);
Abarmat10 = cof10(:,2*nvars+1:3*nvars);

% get the time-dependent decision rules
nregimes1 = length(regime1);

Tmax = regimestart1(nregimes1)-1;  %-1 Tmax is the position of the last period
% when the constraint binds
 
if Tmax > 0
    P = zeros(nvars,nvars,Tmax);
    D = zeros(nvars,Tmax);
    
%     invmat = inv((Astarbarmat*decrulea+Bstarbarmat));
%     P(:,:,Tmax) = -invmat*Cstarbarmat;
%     D(:,Tmax) = -invmat*Dstarbarmat;
    %XXX fix next three lines
    invmat = inv((Abarmat10*decrulea+Bbarmat10));
    P(:,:,Tmax) = -invmat*Cbarmat10;
    D(:,Tmax) = -invmat*Dbarmat10;  
    
    
    
    
    for i = Tmax-1:-1:1        
        
        if violvecbool(i,1)
            invmat = inv(Bbarmat10+Abarmat10*P(:,:,i+1));
            P(:,:,i)=-invmat*Cbarmat10;
            D(:,i) = -invmat*(Abarmat10*D(:,i+1)+Dbarmat10); 
        else
            invmat = inv(Bbarmat+Abarmat*P(:,:,i+1));
            P(:,:,i)=-invmat*Cbarmat;
            D(:,i) = -invmat*(Abarmat*D(:,i+1));
        end
        
    end

    
% Double check the appropriate invmat in each case
% right now -- inherited from previous loop
if Tmax > 1
    if  violvecbool(1,1) 
        E = -invmat*Jbarmat10;
    else
        E = -invmat*Jbarmat;
    end
    
else  % Tmax is equal to 1    
%     invmat = inv((Astarbarmat*decrulea+Bstarbarmat));
%     E = -invmat*Jstarbarmat;
    
        invmat = inv((Abarmat10*decrulea+Bbarmat10));
        E = -invmat*Jbarmat10;
            
    
    
end

Ecurrent = E;

else
  
  Ecurrent = decruleb;

end

% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history = zeros(nvars,nperiods+1);
history(:,1) = init;
errvec = zeros(size(exog_,1),1);

for i = 1:nshocks
    errvec(irfshockpos(i)) = scalefactormod(i);
end

% deal with shocks
irfpos =1;
if irfpos <=Tmax
    history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
        D(:,irfpos) + E*errvec;
else
    history(:,irfpos+1) = decrulea*history(:,irfpos)+decruleb*errvec;
end

% all other periods
for irfpos=2:nperiods+1
    if irfpos <=Tmax
        history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
            D(:,irfpos);
    else
        history(:,irfpos+1) = decrulea*history(:,irfpos);
    end
end


history=history';
zdata = history(2:end,:);