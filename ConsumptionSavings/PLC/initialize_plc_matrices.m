function [GAM, SMAT, PHIMAT, OMEGAMAT, THETA,sModel] = initialize_plc_matrices(par,sModel)

% HOUSEKEEPING AND MATRIX INITIALIZATION
sModel.n        = sModel.nx + 1;             % Number of elements in state vector (including constant)
sModel.my       = (sModel.n+1)*sModel.ny;    % Number of free coefficients.
sModel.index_11 = 1;                         % Index corresponding to alpha_11
sModel.index_12 = 2:sModel.n;                % Index corresponding to alpha_12
sModel.index_21 = sModel.n+1;                % Index corresponding to beta_21
sModel.index_22 = sModel.n+2:2*sModel.n;     % Index corresponding to beta_22

sModel.index_1  = [sModel.index_11 sModel.index_12];
sModel.index_2  = [sModel.index_21 sModel.index_22];

% INITIALIZE PLC MATRICES
GAM{1} = zeros(1,1);               % This is gamma_1
GAM{2} = zeros(sModel.n-1,1);      % This is gamma_2
GAM{3} = zeros(sModel.ny,1);       % This is gamma_y
GAM{4} = zeros(sModel.nx,1);       % This is gamma_x


for iy=1:sModel.ny
    SMAT{iy} = zeros(sModel.n+1,sModel.my);
    SMAT{iy}(1:sModel.n,1+sModel.n*(iy-1):sModel.n*iy) = eye(sModel.n);
    SMAT{iy}(sModel.n+1,sModel.ny*sModel.n+iy) = 1;
end

for ix=1:sModel.nx
    PHIMAT{ix} = zeros(sModel.n+1,1)';
end

for ix=1:sModel.nx
    OMEGAMAT{ix} = zeros(sModel.n+1,sModel.my);
end

THETA     = zeros(sModel.my,1);


% COLLECT MODEL MATRICES
sModel.GAM      = GAM;
sModel.SMAT     = SMAT;
sModel.PHIMAT   = PHIMAT;
sModel.OMEGAMAT = OMEGAMAT;

end

