% Residual function for the simple GK model

function G = system_borrcon_plc(theta_in,sModel,sGrid,par)

%---------------------------------------
% Map parameters
%---------------------------------------
R      = par.R;
BET    = par.BET;
M      = par.M;
GAMMAC = par.GAMMAC;
RHO    = par.RHO;

% Map grid
Y_in = sGrid.Y;
B_in = sGrid.Blag;

% INTEGRATION NODES
EA_prime = sGrid.epsi_nodes(:,1);
WEIGHT   = sGrid.weight_nodes;
nNodes   = length(EA_prime);

% Number of grid points
nGrid    = length(Y_in);
%nb = length(B_in);

%---------------------------------------
% Housekeeping
%---------------------------------------
integrand_1 = NaN(1,nNodes);
res_1 = NaN(nGrid,1);


%---------------------------------------
% Map coefficients
%---------------------------------------
theta_use = [theta_in; 0];

[ALP,~,delta] = get_coefs_plc_general(theta_use,sModel);


%-----------------------------------------------------
% This loop computes the residual for each pair (k,a)
%-----------------------------------------------------

% Counter to loop over collocation points
row = 1;

for igrid = 1:nGrid
    
        
        % States
        Blag = B_in(igrid);
        Y    = Y_in(igrid);
        
        
        % State partition
        x1 = Blag;
        X2 = [1  Y]';
        
        % Basis function
        state_i = [x1; X2];
        
        
        
        if x1<delta*X2     % Reference regime: slack
            
            % Decision rules
            B  = ALP{1}(sModel.index_1)'*state_i;          
            
            % Constraint is slack
            active = 0;
            
        else
            
            
            % Decision rules
            B  = ALP{1}(sModel.index_2)'*state_i;          % Marginal value of Net Worth
            
            % Constraint is binding
            active  = 1;
            
        end
        
        % CONSUMPTION
        c_use = Y - R*Blag + B;
        
        
        %============================
        % Compute Expectations
        %============================
        for iq = 1:nNodes
            
            % Z_prime = log(Y_prime)
            Z_prime = RHO*log(Y) + EA_prime(iq);
            
            % Get Income in level
            Y_prime = exp(Z_prime);

            % State partition
            x1_prime = B;
            X2_prime = [1  Y_prime]';
            
            % Basis function
            state_i_prime = [x1_prime; X2_prime];
            
            
            if x1_prime <delta*X2_prime
                
                % Decision rules
                Bprime  = ALP{1}(sModel.index_1)'*state_i_prime;          
                
            else
                
                % Decision rules
                Bprime  = ALP{1}(sModel.index_2)'*state_i_prime;          
            end
            
            % TRANSFORM TO LEVELS
            c_prime = Y_prime + Bprime - R*B;
            
            % Integrand for residual equations
            integrand_1(1,iq) = (c_prime^-GAMMAC);
            
        end
        
        % Solve integral for residuals
        sol_integ_1  = (integrand_1)*WEIGHT;
                
                
        % Compute Residual for Euler Equation when constraint is slack
        if active == 0
            res_1(row,:) = c_use^-GAMMAC - BET*R*sol_integ_1;
            lambda = 0;
        else
            lambda = c_use^-GAMMAC - BET*R*sol_integ_1;

            res_1(row,:) = lambda*(B - M*Y_in(igrid));
        end
        
        row = row + 1;
        
    
end

G     = res_1;
