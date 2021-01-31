%--------------------------------------------------------------------------
% This code reproduces results in Appendix F. Consumption-Savings Model with 
% Borrowing Constraint% Consumption - Savings Model with Occasionally Binding Constraints" 
% described in the  paper: "Piecewise-Linear Approximations and Filtering 
% for DSGE Models with Occasionally Binding Constraints," by Aruoba,
% Cuba-Borda, Higa-Flores, Schorfheide and Villavazo (2020).
% 
% Written by: Pablo Cuba-Borda
% Federal Reserve Board. Washington, D.C. 
%
% Created: 07/01/2020.
% Last modified: 09/11/2020.
%--------------------------------------------------------------------------

clear all; warning off; close all; 
%-------------------------------------------
% Housekeeping
%-------------------------------------------

addpath External/
addpath OccBin/
addpath OccBin/estobin
addpath FiPIt/
addpath PLC/

% restoredefaultpath
setpathdynare4
warning off

global oo00_  M00_ M10_  

%-------------------------------------------
% MODEL PARAMETERS
%-------------------------------------------
R    = 1.05;
BETA = 0.945;
RHO  = 0.9;
STD_U= 0.010;
M = 1;
GAMMAC = 1; 

% STEADY STATE
Bss = M;
Css = 1 + M - R*M;
LBss=(1-BETA*R)/Css^GAMMAC;
Yss =1;


save PARAM_EXTRA_BABY R BETA RHO STD_U GAMMAC
paramvec_ = [R BETA RHO STD_U M GAMMAC];

% Grid for Decision Rules and Integreation Nodes
%-----------------------------------------------
nb=200;         % Grid points for debt
nz=11;          % Income grid in terms of z=log(Y)

% Rouwenhorst discretization 
[logY, P] = rouwenhorst(RHO,STD_U,nz);
Ygrid = exp(logY);
Zgrid = logY;

% Borrowing grid
bmin  = 0.75*M;
bmax  = M*Ygrid(end) ;
Bgrid = linspace(bmin,bmax,nb);

%==========================================================================
% OCCBIN SOLUTION
%==========================================================================
cd OccBin

% 0: uses csolve; 1: csolve with gradient; 2: fsolve
solver=0;

set(0,'DefaultLineLineWidth',2)
randn('seed',1);
format compact
obs_list = char('c');
err_list = char('eps_u');

modnam     = 'borrcon00';   % Constrained model
modnamstar = 'borrcon10';   % Unconstrained model

%--------------------------------------------------------------------------
% DEFINE CONSTRAINTS
%
% Note that lb_ss is the steady-state value of the multiplier lb (see borrcon_steadystate.m)
% lb_ss is calculated automatically as additional auxiliary parameter
% around line 49 of solve_one_constraint
%--------------------------------------------------------------------------
constraint       = 'lb<-lb_ss';
constraint_relax ='b>bnot';

% SOLUTION OPTIONS
%--------------------------------------------------------------------------
nperiods_solve = 50;        % Number of forward simulations for shooting solution

fprintf('\n ============================');
fprintf('\n Solving with OccBin algorithm')
fprintf('\n ============================');

tic
constraint1 = 'lb<-lb_ss';
constraint_relax1 = 'b>bnot';
curb_retrench =0;
maxiter = 20;

% Create occbin objects
solve_one_constraints_firstcall(modnam,modnamstar);

% Store steady state
zdatass = oo00_.dr.ys;

% Create constraints
[constraint1_difference, iendo1]= process_constraint_with_tokens(constraint1,'_difference',M00_.endo_names,0);
[constraint_relax1_difference, iendo2]= process_constraint_with_tokens(constraint_relax1,'_difference',M00_.endo_names,0);

% Compute OccBin Solution Matrices
[nvars,ys_,endog_,exog_,params, decrulea,decruleb,cof,...
    Jbarmat,cof10,Jbarmat10,Dbarmat10] = get_occbin_solution(modnam,modnamstar,solver,paramvec_);
%
% Compute OccBin Related Objects
[~, i1, ~]=intersect(M00_.endo_names,obs_list,'rows');
[~, i2, ~]=intersect(M00_.exo_names,err_list,'rows');
    
% Current observation (keep for compatibility)
current_obs = [];

% Initial Vector: Steady State
init_val_old  = zeros(length(zdatass),1);

% COMPUTE OCCBIN DECISION RULE
%---------------------------------------------------
for iz=1:length(Zgrid)
    for ib=1:length(Bgrid)
        
        err_vals        = Zgrid(iz);
        init_val_old(1) = (Bgrid(ib)-zdatass(1));
        
        [zdata, zdataconcatenated, ys_, init_out, error_flag, Ecurrent ] = solve_one_constraints_dr(...
            constraint1_difference,constraint_relax1_difference,...
            err_vals',err_list,nperiods_solve,curb_retrench,maxiter,init_val_old);
        
        DR.B(iz,ib)  = zdataconcatenated(1,1) + zdatass(1);
        DR.C(iz,ib)  = zdataconcatenated(1,3) + zdatass(3);
        DR.LB(iz,ib) = (zdataconcatenated(1,5) + zdatass(5));
        DR.EC(iz,ib) = (zdataconcatenated(1,4) + zdatass(4));
        DR.Z(iz,ib)  = zdataconcatenated(1,7) + zdatass(7);
        
        DR_LINEAR.B(iz,ib)  = zdata(1,1) + zdatass(1);
        DR_LINEAR.C(iz,ib)  = zdata(1,3) + zdatass(3);
        DR_LINEAR.EC(iz,ib) = zdata(1,4) + zdatass(4);
        DR_LINEAR.LB(iz,ib) =(zdata(1,5) + zdatass(5));
        DR_LINEAR.Z(iz,ib)  = zdata(1,7) + zdatass(7);
    end
end
time1 = toc;
fprintf('\n *** OccBin solution time: %4.4f seconds *** \n',time1);

cd ..

%==============================================
% SOLVE MODEL TIME ITERATION: FiPIT
%==============================================

% Initial guess
[Ym, Bm ] = ndgrid(Ygrid,Bgrid);
CdecOld = max(1e-200, -R*Bm + (1+M)*Ym);
BdecOld = CdecOld + R*Bm - Ym;

% Initialize objects
Bdec_use = BdecOld;
dist=100;

% Start solution
fprintf('\n ============================');
fprintf('\n Solving with FiPIt algorithm')
fprintf('\n iter, distance')
fprintf('\n ============================');

tic % start timing FiPIt
count = 1;
while dist>1E-10

    Bdec_up = solve_FiPit(Bgrid,P,Ygrid,Bdec_use,[STD_U;RHO;GAMMAC; R; M;BETA]);


    Bdec_new = Bdec_up*0.25 + Bdec_use*(1-0.25);

    
    dist = norm(abs(Bdec_up - Bdec_use));

    if mod(count,10)==0;fprintf('\n %i, dist = %4.11f',count,dist);end
    % update for next iteration
    Bdec_use = Bdec_new;
    
    count = count + 1;
end
% Store decision rules
KbTI = Bdec_use;
KcTI = KbTI - R*Bm + Ym;

% Compute decision rule for multiplier
Bdec = KbTI;
Cdec = Bdec - R*Bm + Ym;

for iz=1:nz
    for ib=1:nb
        
        b_use = Bgrid(ib);
        z_use = Ygrid(iz);
        
        b_prime = Bdec(iz,ib);
        c_use   = z_use - R*b_use + Bdec(iz,ib);
        
        if (b_prime < bmin)
            b_prime = bmin;
            c_use = b_prime - R*b_use + z_use;
        end
        
        
        % Linear interpolation
        % allterp211 linear inter/extrapolation (2 states, 1 policies, 1 stoch comp)
        % Inputs:
        %   x*      :   Grid
        %   x*i     :   Point to evaluate
        %   pf*     :   Policy function
        % Outputs:
        %   o*      :   Interpolated/extrapolated values of dimension x*ipts
        
        % Bdec is in z x B space. 
        for iq=1:nz
            
            int_Ct_prime(iq) = allterp211(Bgrid, Ygrid, Bdec(iq,ib),Ygrid(iq), Cdec');
                        
            int_Ct_prime(iq) = max(int_Ct_prime(iq),1E-20)^-GAMMAC;
            
        end
                
        sol_Ct2 = P(iz,:)*(int_Ct_prime');
                
        Ldec(iz,ib) = max(0,c_use^-GAMMAC - (BETA*R*sol_Ct2));              
        
    end
end
time_FiPIt=toc;
fprintf('\n ');
fprintf('\n *** FiPIt solution time: %4.4f seconds *** \n',time_FiPIt);

%==========================================================================
% SOLVE USING PLC ALGORITHM
%==========================================================================

% SOLUTION OPTIONS
opt.solve_symbolic = 0;     % [1] Resolve symbolic objects [0] Use stored objects. Set to [1] when model changes.
opt.ind_analytic   = 0;     % [1] Analytic derivatives     [0] Numerical Derivatives
opt.smolyakPC      = 0;     % [1] Use smolyak on PC        [0] Smolyak on rectangular grid
opt.smolyak_mu     = 2;     % Smolyak approximation order  
opt.printfigs      = 0;     % [1] Print to pdf             [0] Print figures to pdf

% SIMULATION OPTIONS
opt.Tsim           = 100000;       % Number of periods for simulating ergodic distribution.
opt.Tdrop          = 150;          % Number of initial simulations that will be dropped.

% Parameter Structure
par.R      = R;
par.BET    = BETA;
par.M      = M;
par.STD_U  = STD_U;
par.RHO    = RHO;
par.GAMMAC = GAMMAC;

% USER INPUT: PLC-SETUP (MODEL SPECIFIC)
%--------------------------------------------------------------------------
sModel.nf = 2;                     % Number of equilibrium conditions (include law of motions for exogenous states)
sModel.nx = 2;                     % Number of state variables
sModel.ny = 1;                     % Number of control variables
sModel.nexo = 1;                   % Number of exogenous variables
sModel.nQ   = 5;                   % Number of GH quadrature nodes for integration
sModel.intrule = 0;                % 0: GH, 1: Monomial_1, 2: Monomial_2 

sModel.y_names = {'B'};            % List of control variables
sModel.x_names = {'Blag','Y'};     % List of state variables 

% INITIALIZE PLC MATRICES (REQUIERES USER INPUT)
[GAM, SMAT, PHIMAT, OMEGAMAT, THETA, sModel] = initialize_plc_matrices(par,sModel);

% USER INPUT: PLC LINEAR RESTRICTIONS (MODEL SPECIFIC)
%--------------------------------------------------------------------------
% GAM1             % Blag
sModel.GAM{1}    = [0];

% GAM2              1  Y
sModel.GAM{2}    = [0 par.M];

% GAMY               Y
sModel.GAM{3}    = [-1];

% GAMX
sModel.GAM{4}    = [0 0];

% PLC TRANSITION STATE VARIABLES (MODEL SPECIFIC)
%------------------------------------------------
sModel.PHIMAT{1} = [ 0  0  0 0 ]';                % Bprime
sModel.PHIMAT{2} = [ 0  0  par.RHO 0]';           % z

sModel.OMEGAMAT{1}(1:sModel.n,1:sModel.n)     = eye(sModel.n);      %OMEGA(B)
sModel.OMEGAMAT{1}(sModel.n+1,sModel.n+1:end) = 0;                  %OMEGA(B)

sModel.OMEGAMAT{2}(1:sModel.n,1:sModel.n)     = eye(sModel.n);      %OMEGA(Z)
sModel.OMEGAMAT{2}(sModel.n+1,sModel.n+1:end) = 0;                  %OMEGA(Z)

% Because alp_21{mu} = 0 to satisfy continuity, we impose that here:
sModel.SMAT{1}(sModel.n+1,sModel.ny*sModel.n+1:end) = 0;


% SOLVER OPTIONS
script_solver_options

% QUADRATURE NODES AND WEIGHTS
vcv = diag([par.STD_U].^2);

if sModel.intrule == 0
    [sGrid.n_nodes, sGrid.epsi_nodes, sGrid.weight_nodes] = GH_Quadrature(sModel.nQ,sModel.nexo,vcv);
elseif sModel.intrule == 1
    [sGrid.n_nodes, sGrid.epsi_nodes, sGrid.weight_nodes] = Monomials_1(sModel.nexo,vcv);
elseif sModel.intrule == 2
    [sGrid.n_nodes, sGrid.epsi_nodes, sGrid.weight_nodes] = Monomials_2(sModel.nexo,vcv);
end

%--------------------------------------------------------------------------
% ARRANGE INITIAL GUESS: STATE VECTOR FOLLOWS THIS CONVENTION  [x1 X2]
% WITH X2 INCLUDING A CONSTANT.
%
%    x1 = B1;
%    X2 = [1  A1]';
%
% LIST OF CONTROL VARIABLES APPROXIMATED USING PLC RULES
% {'B'};
%
% ORDER OF VARIABLES
%  1   2   3   4  
% [B1 A1 BLAG1 LB1];
%--------------------------------------------------------------------------

% Coefficients unconstrained regime [B1, 1, A1]
% Obtained from running linear approximation
theta0 = [0.8609 0.3356 -0.1847]';

% SIMULATE ERGODIC SET
%--------------------------------------------------------------------------
% SHOCKS TO SIMULATE ERGODIC SET
rng(123);
shocks_use.ea = randn(opt.Tsim+opt.Tdrop,1);

%INITIAL VALUES FOR SIMULATION
init_use.Blag    = Bss;     % This is Blag
init_use.Zlag    = 0;       % This is log(Y)

% Simulate using PLC with initial guess
simplc0 = simulate_borrcon_plc(theta0,sModel,par,sGrid,shocks_use,init_use,opt.Tsim,opt.Tdrop,1,'');

% COMPUTE DENSITIES OF STATE VARIABLES
%--------------------------------------------------------------------------
tic % start timing PLC solution
Blagsim  = simplc0.Blag;
Ysim     = simplc0.Y;

[sGrid.fb, sGrid.xb] = ksdensity(Blagsim);
sGrid.distb          = prctile(Blagsim,[1 99]);

[sGrid.fy, sGrid.xy] = ksdensity(Ysim);
sGrid.disty         = prctile(Ysim,[1 99]);

% SMOLYAK GRID
%--------------------------------------------------------------------------

% BOUNDS FOR SMOLYAK RECTANGULAR GRID
sGrid.ymax = sGrid.disty(1); sGrid.ymin = sGrid.disty(2);
sGrid.bmax = sGrid.distb(1); sGrid.bmin = sGrid.distb(2);

sGrid.smolyak_use= Smolyak_Grid(sModel.nx,opt.smolyak_mu,Smolyak_Elem_Isotrop(sModel.nx,opt.smolyak_mu));

sGrid.nM = length(sGrid.smolyak_use);


% MAP TO GRID
sGrid.Blag  = (sGrid.bmax - sGrid.bmin)*(sGrid.smolyak_use(:,1) + 1)/2 + sGrid.bmin;
sGrid.Y     = (sGrid.ymax - sGrid.ymin)*(sGrid.smolyak_use(:,2) + 1)/2 + sGrid.ymin;

GRIDUSE = [sGrid.Blag sGrid.Y];

% PLC SOLUTION ALGORITHM
%------------------------------------------
[theta_out,resnorm,res_use,exitflag,output,lambda,jacobian] = lsqnonlin(@(theta_in) system_borrcon_plc(theta_in,sModel,sGrid,par), theta0,[],[],opt.nojac);

% PLC DECISION RULES
%------------------------------------------

theta_use = [theta_out; 0];

[ALP_dr, ~, delta_dr]=get_coefs_plc_general(theta_use,sModel);

Kb_plc = NaN(nz,nb);

for iz=1:nz
    for ib=1:nb
    
        % States
        Blag = Bgrid(ib);
        Y    = Ygrid(iz);
        
        % State partition
        x1 = Blag;
        X2 = [1  Y]';
        
        % Basis function
        state_i = [x1; X2];
        
        % Compute decision rule 
        if x1<delta_dr*X2                 
            % Decision rules
            Kb_plc(iz,ib)  = ALP_dr{1}(sModel.index_1)'*state_i;                          
        else                        
            % Decision rules
            Kb_plc(iz,ib)  = ALP_dr{1}(sModel.index_2)'*state_i;                          
        end
                
    end
end

% Compute consumption decision rule from resource constraint
Kc_plc = Kb_plc - R*Bm + Ym;

% Store decision rules
DRPLC.B = Kb_plc; 
DRPLC.C = Kc_plc;
DRPLC.Ygrid = Ygrid';
DRPLC.Bgrid = Bgrid';

time_plc = toc;

fprintf('\n ************************************** \n');
fprintf('\n Total elapsed time seconds = %4.4f \n', time_plc);
fprintf('\n ************************************** \n');


%==============================================
% SIMULATE FiPIt, OCCBIN AND PLC
%==============================================
rng(123);
nperiods = 10500;
sequence = zscore(randn(nperiods,1))*STD_U;

% Simulating models
fprintf('\n')
fprintf('\n ============================');
fprintf('\n Simulating models using model specific states')
fprintf('\n ============================');

% SIMULATE FiPIt 
%---------------------------------------------------
fprintf('\n *** Simulating FiPIt ***')
Blag = 1;
Zlag = 0;

for tt=1:nperiods
    
    % Simulate FiPIT
    Z_use = RHO*Zlag + sequence(tt,1);
    
    Y_use = exp(Z_use);
    
    SIMTI.C(tt,1) = allterp211(Bgrid, Ygrid, Blag,Y_use, KcTI');
    
    SIMTI.B(tt,1) = allterp211(Bgrid, Ygrid, Blag,Y_use, KbTI');
        
    SIMTI.LAM(tt,1) = allterp211(Bgrid, Ygrid, Blag,Y_use, Ldec');
      
    SIMTI.Z(tt,1)  = Z_use;

    SIMTI.Y(tt,1)  = exp(Z_use);

    % Update States
    Blag = SIMTI.B(tt,1);
    Zlag = SIMTI.Z(tt,1);
    
    %
    if mod(tt,1000)==0; fprintf('\n sim = %i out of %i',tt,nperiods); end

end

% SIMULATE OCCBIN 
%---------------------------------------------------
fprintf('\n *** Simulating OccBin ***')
cd OccBin
% Draw shocks for simulation
err_vals        = sequence;
init_val_old(1) = 0;

[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  sequence,err_list,nperiods_solve);

% COLLECT OCCBIN SIMULATION
SIMOCC.B  = zdatapiecewise(:,1) + zdatass(1);
SIMOCC.C  = zdatapiecewise(:,3) + zdatass(3);
SIMOCC.LAM= zdatapiecewise(:,5) + zdatass(5);
SIMOCC.EC = zdatapiecewise(:,4) + zdatass(4);
SIMOCC.Y  = zdatapiecewise(:,7) + zdatass(7);

cd ..

% SIMULATE PLC 
%---------------------------------------------------
fprintf('\n *** Simulating PLC ***')
Blag = 1;
Zlag = 0;
for tt=1:nperiods
    
    % Simulate PLC
    Z_use = RHO*Zlag + sequence(tt,1);
    
    Y_use = exp(Z_use);
    
    % Simulate PLC
    shocks_use.ea    = sequence(tt,1)/STD_U;
    init_use.Blag    = Blag;             % This is Blag
    init_use.Zlag    = Zlag;            % This is log(Y)

    outsim = simulate_borrcon_plc(theta_out,sModel,par,sGrid,shocks_use,init_use,1,0,1,'');

    SIMPLC.C(tt,1)      = outsim.C;
    SIMPLC.Y(tt,1)      = outsim.Y;
    SIMPLC.Z(tt,1)      = outsim.Z;
    SIMPLC.B(tt,1)      = outsim.B;
    SIMPLC.LAM(tt,1)    = outsim.LAM;
    SIMPLC.bind(tt,1)   = outsim.bind;
        
    % Update States
    Blag = SIMPLC.B(tt,1);
    Zlag = SIMPLC.Z(tt,1);
    
    %
    if mod(tt,1000)==0; fprintf('\n sim = %i out of %i',tt,nperiods); end

end


%==============================================
% RESULTS AND FIGURES
%==============================================

% COMPUTE TABLE OF BINDING PERIODS ACROSS SIMULATION 
%----------------------------------------------------

% COLLECT MULTIPLIERS
LAM.OCC=SIMOCC.LAM(501:10500);
LAM.TI=max(0,SIMTI.LAM(501:10500));
LAM.PLC=max(0,SIMPLC.LAM(501:10500));

% Number of simulations
Tobs = length(LAM.TI);

% Find Binding FiPIt Periods
index_bind = find(LAM.TI>0);
QQ1 = quantile(LAM.TI(index_bind),0.25);
QQ2 = quantile(LAM.TI(index_bind),0.5);
QQ3 = quantile(LAM.TI(index_bind),0.75);

% FIND PERIODS IN WHICH CONSTRAINT BINDS PER QUARTILE
TI_B_QQ1 = find(LAM.TI>0 & LAM.TI<=QQ1);
TI_B_QQ2 = find(LAM.TI>QQ1 & LAM.TI<=QQ2);
TI_B_QQ3 = find(LAM.TI>QQ2 & LAM.TI<=QQ3);
TI_B_QQ4 = find(LAM.TI>QQ3);

% FIND PERIODS IN WHICH PLC AND GLOBAL BIND PER QUARTILE
PLC_B_QQ1 = sum(LAM.PLC(TI_B_QQ1)>0)/length(TI_B_QQ1);
PLC_B_QQ2 = sum(LAM.PLC(TI_B_QQ2)>0)/length(TI_B_QQ2);
PLC_B_QQ3 = sum(LAM.PLC(TI_B_QQ3)>0)/length(TI_B_QQ3);
PLC_B_QQ4 = sum(LAM.PLC(TI_B_QQ4)>0)/length(TI_B_QQ4);

% FIND PERIODS IN WHICH OCC AND GLOBAL BIND PER QUARTILE
OCC_B_QQ1 = sum(LAM.OCC(TI_B_QQ1)>0)/length(TI_B_QQ1);
OCC_B_QQ2 = sum(LAM.OCC(TI_B_QQ2)>0)/length(TI_B_QQ2);
OCC_B_QQ3 = sum(LAM.OCC(TI_B_QQ3)>0)/length(TI_B_QQ3);
OCC_B_QQ4 = sum(LAM.OCC(TI_B_QQ4)>0)/length(TI_B_QQ4);

% STORE RESULTS IN MATRIX
MAT2(1,:) = [PLC_B_QQ1 PLC_B_QQ2 PLC_B_QQ3 PLC_B_QQ4];
MAT2(2,:) = [OCC_B_QQ1 OCC_B_QQ2 OCC_B_QQ3 OCC_B_QQ4];
TableBinding = array2table(round(MAT2*100,1),'VariableNames',{'Q1','Q2','Q3','Q4'},'RowNames',{'PLC','OCCBIN'});
fprintf('\n')
fprintf('\n =============================================================');
fprintf('\n Fraction of binding simulations per quartile relative to FiPIt');
fprintf('\n (Solution specific state variables)');
fprintf('\n ============================================================= \n');
disp(TableBinding) 

% FIGURE OPTIONS
%--------------------------------------

% Compute kernel density of state variables B and Y 
[~,Bti] = ksdensity(SIMTI.B(501:10:end));
[~,Yti] = ksdensity(SIMPLC.Y(501:10:end));

% Fitted distributions
pdB = fitdist(SIMTI.B(501:10:end),'Kernel','Kernel','epanechnikov');
pdY = fitdist(SIMPLC.Y(501:10:end),'Kernel','Kernel','epanechnikov');

% Compute 95% interval of state variables
blim_use = icdf(pdB,[0.025, 0.975]);
ylim_use = icdf(pdY,[0.025, 0.975]);

% Plotting options
CT  =  cbrewer('div', 'RdYlGn', 8);
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')  
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultTextFontSize',20);
set(0,'DefaultLineLineWidth', 3);

% FIGURE 1. DECISION RULES
%--------------------------------------

fig = figure('PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.1 0.1 10.5 8]); 
subplot(221)
shadedplot(ylim_use,[0 0],[1.5 1.5],[0.9 0.9 0.9]); hold on; grid off;
plot(Ygrid,DR.C(:,156),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(Ygrid,KcTI(:,156),'-','color',CT(8,:))
plot(Ygrid,DRPLC.C(:,156),'-k')
plots=get(gca, 'Children');
l1 = legend(plots([1,2, 3]), {'PLC','FiPIt', 'OccBin'},'Orientation','Horizontal'); legend boxoff
xlabel('Income $(Y)$'); title('$C$');%ylabel('$C$')
axis tight; ylim([0.82, 1.04]); 

subplot(222)
shadedplot(ylim_use,[0 0],[1.5 1.5],[0.9 0.9 0.9]); hold on; grid off;
plot(Ygrid,DR.B(:,156),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(Ygrid,KbTI(:,156),'-','color',CT(8,:))
plot(Ygrid,DRPLC.B(:,156),'-k')
xlabel('Income $(Y)$'); title('$B^\prime$');%ylabel('$B^\prime$');
axis tight;  ylim([0.92, 1.03]); 

subplot(223)
shadedplot(blim_use,[0 0],[1.5 1.5],[0.9 0.9 0.9]); hold on; grid off;
plot(Bgrid,DR.C(6,:),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(Bgrid,KcTI(6,:),'-','color',CT(8,:))
plot(Bgrid,DRPLC.C(6,:),'-k')
xlabel('Debt $(B)$'); title('$C$');%ylabel('$C$')
axis tight; ylim([0.82, 1.04]); xlim([0.85 1.07]);

subplot(224)
shadedplot(blim_use,[0 0],[1.5 1.5],[0.9 0.9 0.9]); hold on; grid off;
plot(Bgrid,DR.B(6,:),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(Bgrid,KbTI(6,:),'-','color',CT(8,:))
plot(Bgrid,DRPLC.B(6,:),'-k')
xlabel('Debt $(B)$'); title('$B^\prime$'); %ylabel('$B^\prime$')
axis tight;  ylim([0.88 1.01]);  xlim([0.85 1.07]); 

set(l1,'Position',[0.001 0.001 1.6 0.015]);

print('-dpdf','-r300','Figures/fig_dr_occbin_plc_fipit.pdf');

% FIGURE 1. SIMULATED PATHS
%--------------------------------------

t0 = 301;
t1 = 50;
fig1 = figure('PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.1 0.1 10.5 8]); 
subplot(221)
plot(SIMOCC.C(t0:t0+t1),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(SIMTI.C(t0:t0+t1),'-','color',CT(8,:));
plot(SIMPLC.C(t0:t0+t1),'k-'); hold on; 
plots=get(gca, 'Children');
l1 = legend(plots([1,2, 3]), {'PLC','FiPIt', 'OccBin'},'Orientation','Horizontal'); legend boxoff
xlabel('time'); title('Consumption $(C)$'); 
axis tight; %ylim([0.89 1.03]);

subplot(222)
plot(SIMOCC.B(t0:t0+t1),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(SIMTI.B(t0:t0+t1),'-','color',CT(8,:)); hold on; 
plot(SIMPLC.B(t0:t0+t1),'k-'); hold on; 
xlabel('time'); title('Borrowing $(B^\prime)$');
axis tight; %ylim([0.95 1.05]);

subplot(223)
plot(SIMOCC.LAM(t0:t0+t1),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(SIMTI.LAM(t0:t0+t1),'color',CT(8,:),'LineWidth',2); hold on; 
plot(SIMPLC.LAM(t0:t0+t1),'k-'); hold on; 
xlabel('time'); title('Multiplier $(\lambda)$');
axis tight; 

subplot(224)
plot(SIMOCC.Y(t0:t0+t1),'-','color',CT(1,:),'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on; 
plot(SIMTI.Y(t0:t0+t1),'-','color',CT(8,:),'LineWidth',2); hold on; 
plot(SIMPLC.Y(t0:t0+t1),'k-'); hold on; 
xlabel('time'); title('Income $(Y)$');
axis tight; %ylim([0.95 1.05]);
set(l1,'Position',[0.01 0.01 1.7 0.015]);
print('-dpdf','-r300','Figures/fig_simulated_paths.pdf');


% FIGURE 1. ERGODIC DISTRIBUTION OF ASSETS
%-----------------------------------------
[focc,xocc] = ksdensity(-SIMOCC.B(501:10:end));
[fplc,xplc] = ksdensity(-SIMPLC.B(501:10:end));
[fti,xti] = ksdensity(-SIMTI.B(501:10:end));
fig2 = figure('PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.1 0.1 10.5 8]); 
plot(xocc,focc,'-','color',CT(1,:),'LineWidth',3,'MarkerSize',5,'MarkerFaceColor',CT(1,:)); hold on 
plot(xti,fti,'-','color',CT(8,:),'LineWidth',3); hold on 
plot(xplc,fplc,'k-','LineWidth',3); hold on 
plots=get(gca, 'Children');
l1 = legend(plots([1,2, 3]), {'PLC','FiPIt', 'OccBin'},'Orientation','Horizontal'); legend boxoff
xlabel('Assets'); title('Ergodic Distribution of Assets'); 
axis tight; %ylim([0.89 1.03]);
print('-dpdf','-r300','Figures/fig_asset_distribution.pdf');

