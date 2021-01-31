% Function the simulates consumption-savings model using using the PLC
% solution algorithm 
%
% Inputs:
%           theta_in : vector of coefficients from model solution.
%           sModel   : structure with PLC matrices
%           par      : structure with model parameters
%           sGrid    : structure with integration nodes and weights
%           shocks_in: structure with realizations of the exogenous shocks
%           init     : structure with initial values of state vector
%           T        : number of simulations to be produced
%           Tdrop    : burn-in period for simulation
%           ind_extra: binary indictor to compute expectations
%
%--------------------------------------------------------------------------

function outsim = simulate_borrcon_plc(theta_in,sModel,par,sGrid,shocks_in,init,T,Tdrop,ind_extra,stringnore)

%--------------------------------------------------------------------------
% MAP PARAMETERS
%--------------------------------------------------------------------------
R      = par.R;
BET    = par.BET;
GAMMAC = par.GAMMAC;
STD_U  = par.STD_U;
RHO    = par.RHO;
%--------------------------------------------------------------------------
% Map Coefficients
%--------------------------------------------------------------------------

if strcmp('ignore',stringnore)
    ignore=1;
else
    ignore=0;
end

% MAP INITIAL GUESS TO COEFFICIENTS FOR PLC

theta_use = [theta_in; 0];

[ALP,~,delta] = get_coefs_plc_general(theta_use,sModel);

%--------------------------------------------------------------------------
% Map Integration Nodes and Weights
%--------------------------------------------------------------------------

if ind_extra
    % INTEGRATION NODES
    EA_prime = sGrid.epsi_nodes(:,1);
    WEIGHT   = sGrid.weight_nodes;
    nNodes   = length(EA_prime);
    integrand_lam = NaN(1,nNodes);
end

%--------------------------------------------------------------------------
%   Preallocate vectors
%--------------------------------------------------------------------------
TT = T + Tdrop;                                   % Total # of observations

Bs    = zeros(TT,1);
Cs    = zeros(TT,1);
LAMs  = zeros(TT,1);
Zs    = zeros(TT,1);
Ys    = zeros(TT,1);
binds = zeros(TT,1);

%--------------------------------------------------------------------------
% Map shocks
%--------------------------------------------------------------------------
eas = shocks_in.ea;

%--------------------------------------------------------------------------
%   Main Simulation Code
%--------------------------------------------------------------------------

for tt = 1 : TT
    
    if tt == 1
        
        Blag     = init.Blag;
        Zlag     = init.Zlag;
        
    else
        
        Blag     = Bs(tt-1);
        Zlag     = Zs(tt-1);
        
    end
    
    % Here Z = ln(Y)
    Zs(tt)   = RHO*Zlag + STD_U*eas(tt);
    
    % Get Income in level
    Ys(tt)   = exp(Zs(tt));
    
    % State partition
    x1 = Blag;
    X2 = [1  Ys(tt)]';
    
    % Basis function
    state_i = [x1; X2];
    
    % ignore switch: solve reference regime
    if ignore
        % Decision rules
        Bs(tt)  = ALP{1}(sModel.index_2)'*state_i;          % Marginal value of Net Worth
        
        % Update binding indicator
        binds(tt) = 1;        
        
    else % solve PLC
        
        if x1<delta*X2    % Reference regime: constraint slack
            
            % Decision rules
            Bs(tt)  = ALP{1}(sModel.index_1)'*state_i;
            
            LAMs(tt) = 0;
            % Update binding indicator
            binds(tt) = 0;
            
        else
            
            % Decision rules
            Bs(tt)  = ALP{1}(sModel.index_2)'*state_i;          
            
            % Update binding indicator
            binds(tt) = 1;
        end
    end
    
    % CONSUMPTION
    Cs(tt) = Ys(tt) - R*Blag + Bs(tt);
    
    
    % Compute objects that require evaluating expectations
    if ind_extra
        
        B   = Bs(tt);
        
        for count = 1:nNodes
            
            Z_prime = RHO*Zs(tt) + EA_prime(count);
            
            % Get Income in level
            Y_prime = exp(Z_prime);
            
            % State partition
            x1_prime = B;
            X2_prime = [1  Y_prime]';
            
            % Basis function
            state_prime = [x1_prime; X2_prime];
            
            % Ignore switch: solve reference regime
            if ignore
                % Decision rules
                B_prime  = ALP{1}(sModel.index_2)'*state_prime;

            else
                if x1_prime < delta*X2_prime
                    
                    % Decision rules
                    B_prime  = ALP{1}(sModel.index_1)'*state_prime;
                    
                else
                    
                    % Decision rules
                    B_prime  = ALP{1}(sModel.index_2)'*state_prime;
                    
                end
            end
            
            
            % TRANSFORM TO LEVELS
            C_prime = Y_prime - R*B + B_prime;
            
            
            % Integrand for residual equations
            integrand_lam(:,count)    = (C_prime)^-GAMMAC;
            
        end
        
        % Integrand for residual equations
        
        sol_int = integrand_lam * WEIGHT;
        
        % Compute multiplier when constraint binds
        if binds(tt)
            LAMs(tt) = Cs(tt)^-GAMMAC - BET*R*sol_int;
        end
        
    end
    
    
end

%--------------------------------------------------------------------------
% Trim simulation and output results
%--------------------------------------------------------------------------
% Level variables
outsim.C      = Cs(Tdrop+1:T+Tdrop);
outsim.Y      = Ys(Tdrop+1:T+Tdrop);
outsim.Z      = Zs(Tdrop+1:T+Tdrop);
outsim.B      = Bs(Tdrop+1:T+Tdrop);
outsim.LAM    = LAMs(Tdrop+1:T+Tdrop);
outsim.bind   = binds(Tdrop+1:T+Tdrop);


% Shock innovations
outsim.ea    = eas(Tdrop+1:T+Tdrop);


% Get lagged simulations of endogenous states %
if Tdrop==0
    
    outsim.Blag(1,1)  = init.Blag;
    outsim.Zlag(1,1)  = init.Zlag;

    outsim.Blag(2:T,1) = Bs(1:end-1);
    outsim.Zlag(2:T,1) = Zs(1:end-1);

else
    
    outsim.Blag     = Bs(Tdrop:end-1);
    outsim.Zlag     = Zs(Tdrop:end-1);
    
end

end

