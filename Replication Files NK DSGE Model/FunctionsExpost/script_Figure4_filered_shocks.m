% Script to implement the Forward Guidance exercise in NK with ZLB constraint.
%
% This script reproduces Figure 4: Ex-Post Policy Analysis 
% of "Piecewise-Linear Approximations and Filtering for DSGE Models 
% with Occasionally Binding Constraints"  Aruoba, Cuba-Borda, Higa-Flores, 
% Schorfheide and Villalvazo (2020).
%
% Written by: S. Boragan Aruoba and Pablo Cuba-Borda
% Created : August 6, 2012
% This version : July 29, 2020
%
%--------------------------------------------------------------------------
clear all; clc;  close all;

global Tsim Tirf Tdrop
global ind_scale ind_impose
global ind_period_one
global plotpath  start_plot Tsim_plot
global size_shock pos   ncoefs
global tag_name
global ind_percentile
global options_fzero
global e_min_use e_max_use
global size_spending_intervention


%=== OPTIONS FOR MODEL SELECTION ===%

ind_pick_star = 1;
fname = [];
%--------------------------------------------------------------------------
%                        OPTIONS FOR MODEL SELECTION: SOLUTION AND
%                        POSTERIOR DRAWS FROM SMC SAMPLER:
% Model vintage in paper is: _Prior1_US_RClin1_GridSmolyak_Mhrun1

% Draws and model solution are stored in the following subfolders:
% Subfolders: /Filtered             : Benchamark draws + benchmark solution
%--------------------------------------------------------------------------
modeldir   = './COPFexactZLB_Prior1_US_RClin1_GridSmolyak_Mhrun1/Filtered/';
resultsdir = './COPFexactZLB_Prior1_US_RClin1_GridSmolyak_Mhrun1/Results/';
modelname  = 'COPFexactZLB_Filtered_Prior1_US_RClin1_GridSmolyak_Mhrun1_';

%=== OPTIONS FOR RESULTS ===%

period_use = 2009.25;              % Set period for experiment on '2009Q2'
period_plot_start = 2007.75;       % Set period to start plot
ind_add_shocks    = 0;             % If 1 then the g shock is added and monetary policy is used to set R = 1 or reduce it as much as possible;
                                   % If 0 then the g shock is substracted and monetary policy shock is set to zero.


Tsim  = 8;                         % Number of periods to keep track after initial period

length_intervention_fg = 8;        % The number of period where monetary policy forces R = 1.
length_intervention_fiscal = 1;    % (keep at 1 for now) The number of periods where government spending increases.

size_spending_intervention = 0.027*100; %The increase in government spending to implement during the intervention.  (a one-time increase in period 1)

intervention_bound = 100;          % The maximum size of the intervention in annualized basis points

max_er_shock = 2;                  % Maximum size of the monetary policy shock (in s.d. units)

ind_no_er_shock = 1;               % If set to 1 there are no monetary policy shocks except for the intervention.

low_prctile  = 20;                 % The percentiles to tabulate
high_prctile = 80;

plotpath = strcat(pwd,'/');        % Save .pdf and .txt results in current directory

Tdrop = 0;
ind_impose = 1;                    % Keep at 1.
ind_dyn_shocks = 0;                % Keep at 0
seedname  = 'mrg32k3a';            % String to generate RandomNumber Streams

%=== OPTIONS FOR THE EXERCISE ===%

Tirf  = Tsim;                      % Length of Impulse Responses
Nrep  = 1000;                      % Number of repetitions to construct IRFs.
ind_irf_mean = 0;                  % [1] Mean IRF [0] Median IRF
ind_irf_band = 1;                  % [1] Produces 10-90 Bands
% ind_percentile = 0;              % [1] Uses point-wise percentiles [0] Picks the percentiles based on impact response


e_min_use = - max_er_shock;
e_max_use =   max_er_shock;



%--------------------------------------------------------------------------
%                       EXPERIMENT DATES
%--------------------------------------------------------------------------

options_fzero = optimset('display','off','MaxFunEvals',100,'MaxIter',20,'TolFun',1e-12,'TolX',1e-8,'FunValCheck','On');
options_fmin = optimset('display','off','MaxFunEvals',1000,'MaxIter',100,'TolFun',1e-12,'TolX',1e-8,'FunValCheck','On');


% ############### LOAD COEFFICIENTS AND PARAMETERS ########################


%% ############### LOAD COEFFICIENTS AND PARAMETERS ########################

%--------------------------------------------------------------------------
%                   LOAD PARAMETERS AND MODEL SOLUTION
%--------------------------------------------------------------------------


sModel = fLoadCanonicalForm(modeldir,modelname);
 
sParam = fCreateParams([modeldir modelname]);

fnames = fieldnames(sParam);
for ii=1:length(fnames)
    eval([ char(fnames(ii)) ' = sParam.' char(fnames(ii)) ';'])       
end


ind_dyn_shocks_sim = -1;            % Keep for compatibility, keep at -1.
Tsim_sim  = 10000;                  % Keep for compatibility
Tdrop_sim = 150;                    % Keep for compatibility

% GET NUMBER OF OPEN FIGURES
nopenfigs = get(0,'Children');      % Count open figures, if any.



%--------------------------------------------------------------------------
%                   LOAD AND ARRANGE FILTERED STATES
% The COPF_FilteredStates.csv contains the following:
% E[s(t)|Y_{1:t}], hat-eta_t, hat-eps_t, regime indicator (n = 1 vs b = 2)
% 
% where the vector of states is: 
% s_t = [er_t, g_t, z_t, d_t, y_t, pi_t, c_t, R_t]
% and the order of the shocks is:
% er, eg, ez, ed
%--------------------------------------------------------------------------

RAW   = csvread([modeldir modelname 'FilteredStates_lag.csv']);
%%% LOAD t-1 FILTERED STATES %%%
% s_t-1 = [er_t-1, g_t-1, z_t-1, d_t-1, y_t-1, pi_t-1, c_t-1, R_t-1|Y_{1:t}]
LOAD_STATES_LAGGED.er      = RAW(:,1);
LOAD_STATES_LAGGED.g       = RAW(:,2);
LOAD_STATES_LAGGED.z       = RAW(:,3);
LOAD_STATES_LAGGED.d       = RAW(:,4);
LOAD_STATES_LAGGED.y       = RAW(:,5);
LOAD_STATES_LAGGED.pi      = RAW(:,6);
LOAD_STATES_LAGGED.c       = RAW(:,7);
LOAD_STATES_LAGGED.R       = RAW(:,8);
clearvars RAW;

% LOAD FILTERED SHOCKS
RAW   = csvread([modeldir modelname 'FilteredEps.csv']);
LOAD_SHOCKS.heps_er = RAW(:,1);
LOAD_SHOCKS.heps_eg = RAW(:,2);
LOAD_SHOCKS.heps_ez = RAW(:,3);
LOAD_SHOCKS.heps_ed = RAW(:,4);
clearvars RAW;

% LOAD TIME-t STATES
RAW   = csvread([modeldir modelname 'FilteredStates.csv']);
LOAD_STATES.R     = RAW(:,8);
LOAD_STATES.ylag  = RAW(:,9);
clearvars RAW;


% LOAD DATA
DATA   = csvread([modeldir modelname 'Data.csv']);
Year   = DATA(:,5);

% MAP DATA
DATA_GDPDEF.dely = DATA(:,1);
DATA_GDPDEF.pi   = DATA(:,3);
DATA_GDPDEF.R    = DATA(:,4);
G_SHOCKS         = LOAD_STATES_LAGGED.g + log(sParam.gstar);

% PERIOD FOR EXPERIMENTS

pick_initial_period = find(Year==period_use);
pick_0_period = find(Year==period_plot_start);


% INITIALIZE STATES:  E[s(t-1)|Y_{1:t}] + y_t-1
init_0 = [LOAD_STATES_LAGGED.er(pick_0_period)
    LOAD_STATES_LAGGED.g(pick_0_period)
    LOAD_STATES_LAGGED.z(pick_0_period)
    LOAD_STATES_LAGGED.d(pick_0_period)
    LOAD_STATES_LAGGED.y(pick_0_period)
    LOAD_STATES_LAGGED.pi(pick_0_period)
    LOAD_STATES_LAGGED.c(pick_0_period)
    LOAD_STATES_LAGGED.R(pick_0_period)
    LOAD_STATES.ylag(pick_0_period)];


init_use = [LOAD_STATES_LAGGED.er(pick_initial_period)
    LOAD_STATES_LAGGED.g(pick_initial_period)
    LOAD_STATES_LAGGED.z(pick_initial_period)
    LOAD_STATES_LAGGED.d(pick_initial_period)
    LOAD_STATES_LAGGED.y(pick_initial_period)
    LOAD_STATES_LAGGED.pi(pick_initial_period)
    LOAD_STATES_LAGGED.c(pick_initial_period)
    LOAD_STATES_LAGGED.R(pick_initial_period)
    LOAD_STATES.ylag(pick_initial_period)];


% GET PATH OF REALIZED SHOCKS %
shocks_realized(1,:) = LOAD_SHOCKS.heps_er(pick_0_period:pick_initial_period + Tsim-1);
shocks_realized(2,:) = LOAD_SHOCKS.heps_eg(pick_0_period:pick_initial_period + Tsim-1);
shocks_realized(3,:) = LOAD_SHOCKS.heps_ez(pick_0_period:pick_initial_period + Tsim-1);
shocks_realized(4,:) = LOAD_SHOCKS.heps_ed(pick_0_period:pick_initial_period + Tsim-1);

zeta_star = 1 - (1 / gstar);


%################# PRODUCE RESULTS ########################################

rep_problem = 1;

if ind_add_shocks == 1
    
    % Get new vector of realized shocks %
    
    shocks_use = shocks_realized;

    % GET SHOCKS FOR FG (FISCAL + MONETARY)
    [shocks_use_inter,cases,monotone]  = solve_for_forward_guidance_old(sModel,sParam,...
            shocks_use,init_use,length_intervention_fg,...
            length_intervention_fiscal,intervention_bound);

    % GET ONLY FISCAL
    
    shocks_use_justfiscal = shocks_use;
    
    shocks_use_justfiscal(2,1:length_intervention_fiscal) = shocks_use(2,1:length_intervention_fiscal) + size_spending_intervention;
    
    % SIMULATE PATHS %
    
    [~,~,~,SIM]            = fSimulatePLC(sModel,sParam,init_use,shocks_realized);
    [~,~,~,SIM_INTER]      = fSimulatePLC(sModel,sParam,init_use,shocks_use_inter);
    [~,~,~,SIM_JUSTFISCAL] = fSimulatePLC(sModel,sParam,init_use,shocks_use_justfiscal);


elseif ind_add_shocks == 0
    
    % In this case, the shocks that are observed are assumed to contain the
    % g increase and the reduction in eps_r. 
    
    % GET SHOCKS FOR FG (FISCAL + MONETARY)
    
    shocks_use_inter = shocks_realized;
    
    % Number of periods before ARRA shock;
    tWindow = pick_initial_period - pick_0_period;
    
    % GET ONLY FISCAL    
    shocks_use_justfiscal = shocks_realized;
    shocks_use_justfiscal(1,tWindow+1:tWindow+length_intervention_fg) = zeros(length_intervention_fg,1);

    % NO INTERVENTION    
    shocks_use = shocks_realized; 
    shocks_use(1,tWindow+1:tWindow+length_intervention_fg)     = zeros(length_intervention_fg,1);
    shocks_use(2,tWindow+1:tWindow+length_intervention_fiscal) = shocks_realized(2,tWindow+1:tWindow+length_intervention_fiscal) - size_spending_intervention;

    % NO ed, ez, eg, SHOCKS AFTER 2009Q2
    shocks_use_zero = shocks_realized;

    shocks_use_zero(1,tWindow+1:tWindow+length_intervention_fg) = 0;
    shocks_use_zero(2,tWindow+1:tWindow+length_intervention_fg) = 0;
    shocks_use_zero(3,tWindow+1:tWindow+length_intervention_fg) = 0;
    shocks_use_zero(4,tWindow+1:tWindow+length_intervention_fg) = 0;
    
    % SIMULATE PATHS %    
    [~,~,~,SIM]            = fSimulatePLC(sModel,sParam,init_0,shocks_use);
    [~,~,~,SIM_INTER]      = fSimulatePLC(sModel,sParam,init_0,shocks_use_inter);
    [~,~,~,SIM_JUSTFISCAL] = fSimulatePLC(sModel,sParam,init_0,shocks_use_justfiscal);
    [~,~,~,SIM_ZERO]       = fSimulatePLC(sModel,sParam,init_0,shocks_use_zero);


end
    
if sum(isnan(SIM_INTER.R)) == 0    
    
    if length_intervention_fg > 0 && ind_add_shocks == 1
        
        SIM_INTER.cases = cases;
        
        SIM_INTER.monotone = monotone;
        
    end
    
 % PREPARE MULTIPLIERS FG%
    SIM_INTER.MU_NUM   = (SIM_INTER.BIGY - SIM.BIGY);
    SIM_INTER.MU_DENOM = (SIM_INTER.BIGG - SIM.BIGG);
    SIM_INTER.MU       =  SIM_INTER.MU_NUM ./ SIM_INTER.MU_DENOM;
    
    % CUMULATE EFFECTS %
    SIM_INTER.MU_NUM_CUM   = cumsum(SIM_INTER.MU_NUM);
    SIM_INTER.MU_DENOM_CUM = cumsum(SIM_INTER.MU_DENOM);
    SIM_INTER.MU_CUM       = SIM_INTER.MU_NUM_CUM ./  SIM_INTER.MU_DENOM_CUM;
        
    % PREPARE MULTIPLIERS JUST FISCAL %
    SIM_JUSTFISCAL.MU_NUM   = (SIM_JUSTFISCAL.BIGY - SIM.BIGY);
    SIM_JUSTFISCAL.MU_DENOM = (SIM_JUSTFISCAL.BIGG - SIM.BIGG);
    SIM_JUSTFISCAL.MU       =  SIM_JUSTFISCAL.MU_NUM ./ SIM_JUSTFISCAL.MU_DENOM;
    
    % CUMULATE EFFECTS %
    SIM_JUSTFISCAL.MU_NUM_CUM   = cumsum(SIM_JUSTFISCAL.MU_NUM);
    SIM_JUSTFISCAL.MU_DENOM_CUM = cumsum(SIM_JUSTFISCAL.MU_DENOM);
    SIM_JUSTFISCAL.MU_CUM       = SIM_JUSTFISCAL.MU_NUM_CUM ./  SIM_JUSTFISCAL.MU_DENOM_CUM;
    
else
    
    SIM_PROBLEM.y(:,rep_problem)  = SIM.y;
    SIM_PROBLEM.R(:,rep_problem)  = SIM.R;
    SIM_PROBLEM.pi(:,rep_problem) = SIM.pi;
    SIM_PROBLEM.c(:,rep_problem)  = SIM.c;
    SIM_PROBLEM.z(:,rep_problem)  = SIM.z;
    SIM_PROBLEM.g(:,rep_problem)  = SIM.g;
    SIM_PROBLEM.er(:,rep_problem) = SIM.er;
    
    rep_problem  = rep_problem + 1;
    
end
    
%------------------------------
% PRODUCE SERIES FOR PLOTING 
%------------------------------

% Level of Intervention
    
plot_ep_star_R_INTER    = 400 * (SIM_INTER.R -1);
plot_ep_star_pi_INTER   = 400 * (SIM_INTER.pi -1);
plot_ep_star_Y_INTER    = SIM_INTER.BIGY;
plot_ep_star_dely_INTER = SIM_INTER.dely;
plot_ep_star_er_INTER   = SIM_INTER.er;
plot_ep_star_eg_INTER   = SIM_INTER.eg;
plot_ep_star_g_INTER    = SIM_INTER.g;
plot_ep_star_z_INTER    = SIM_INTER.z;
plot_ep_star_d_INTER    = SIM_INTER.d;

% Intervention Relative to No Intervention

plot_ep_star_Y_change_relnopol  = 100 * (SIM_INTER.BIGY ./ SIM.BIGY - 1);
plot_ep_star_pi_change_relnopol = 400 * (SIM_INTER.pi - SIM.pi);
plot_ep_star_R_change_relnopol  = 400 * (SIM_INTER.R - SIM.R);

% Intervention Relative to Just Fiscal Intervention

plot_ep_star_Y_change_relfiscal  = 100 * (SIM_INTER.BIGY ./ SIM_JUSTFISCAL.BIGY - 1);
plot_ep_star_pi_change_relfiscal = 400 * (SIM_INTER.pi - SIM_JUSTFISCAL.pi);
plot_ep_star_R_change_relfiscal  = 400 * (SIM_INTER.R - SIM_JUSTFISCAL.R);

% Fiscal Intervention Relative to Benchmark

plot_ep_star_Y_change_fiscalrelnopol  = 100 * (SIM_JUSTFISCAL.BIGY ./ SIM.BIGY - 1);
plot_ep_star_pi_change_fiscalrelnopol = 400 * (SIM_JUSTFISCAL.pi - SIM.pi);
plot_ep_star_R_change_fiscalrelnopol  = 400 * (SIM_JUSTFISCAL.R - SIM.R);

% Level of Just Fiscal Intervention

plot_ep_star_er_JUSTFISCAL = SIM_JUSTFISCAL.er;
plot_ep_star_eg_JUSTFISCAL = SIM_JUSTFISCAL.eg;
plot_ep_star_R_JUSTFISCAL  = 400 * (SIM_JUSTFISCAL.R -1 );
plot_ep_star_pi_JUSTFISCAL = 400 * (SIM_JUSTFISCAL.pi -1);
plot_ep_star_Y_JUSTFISCAL = SIM_JUSTFISCAL.BIGY;
plot_ep_star_dely_JUSTFISCAL = SIM_JUSTFISCAL.dely;

% Level of No Policy Intervention

plot_ep_star_er_NOINTER = SIM.er;
plot_ep_star_eg_NOINTER = SIM.eg;
plot_ep_star_ez_NOINTER = SIM.ez;
plot_ep_star_R_NOINTER = 400 * (SIM.R - 1);
plot_ep_star_pi_NOINTER  = 400 * (SIM.pi -1);
plot_ep_star_Y_NOINTER = SIM.BIGY;
plot_ep_star_dely_NOINTER = SIM.dely;
plot_ep_star_g_NOINTER = SIM.g;
plot_ep_star_z_NOINTER = SIM.z;
plot_ep_star_d_NOINTER = SIM.d;

% Level of No Policy Intervention

plot_ep_star_g_ZERO = SIM_ZERO.g;
plot_ep_star_z_ZERO = SIM_ZERO.z;
plot_ep_star_d_ZERO = SIM_ZERO.d;



plot_ep_star_init = init_use;

    
%% PLOT FIGURES %

%-----------------------------------------------
% Figure 4. Ex post Policy Analysis
%-----------------------------------------------

% Set figure x-axis
xt = (1:tWindow+Tirf);
Tperiods = tWindow+Tirf;
mydates = {'07Q4','','','','08Q4','','','','09Q4','','','','10Q4',''};


% unconditional standard deviations
sig_dhat = sig_d/sqrt(1-rho_d^2);
sig_zhat = sig_z/sqrt(1-rho_z^2);
sig_ghat = sig_g/sqrt(1-rho_g^2);

figure(1); clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.2 0.2 10 8.5]);
subplot(3,2,2)
plot(xt,plot_ep_star_er_INTER' ,'k','linewidth',2);hold on;
plot(xt,plot_ep_star_er_NOINTER' ,'b','linewidth',2);
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Monetary Policy Shock Innovations ($\epsilon_{R,t}$)','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -2 2]);set(gca,'FontSize',12); box off; grid on; 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
h_legend = legend('Intervention','No Intervention','location','ne');
set(h_legend,'FontSize',12); legend box off
vline(7,'r-');

subplot(3,2,1)
plot(xt,plot_ep_star_eg_INTER','k','linewidth',2);hold on;
plot(xt,plot_ep_star_eg_NOINTER','b','linewidth',2);hold on;
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Gov. Spending Shock Innovations ($\epsilon_{g,t}$)','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -2 4]); set(gca,'FontSize',12); box off; grid on;
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
h_legend = legend('Intervention','No Intervention','location','ne');
set(h_legend,'FontSize',12); legend box off
vline(7,'r-');

subplot(3,3,4)
plot(xt,log(plot_ep_star_g_INTER'/sParam.gstar)/sig_ghat,'k-','linewidth',2);hold on;
plot(xt,log(plot_ep_star_g_ZERO'/sParam.gstar)/sig_ghat,'r--','linewidth',2);hold on;
plot(xt,log(plot_ep_star_g_NOINTER'/sParam.gstar)/sig_ghat,'b','linewidth',2);hold on;
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Gov. Spending Shock','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -1 2]); set(gca,'FontSize',12); box off; grid on;
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

subplot(3,3,5)
plot(xt,log(plot_ep_star_z_INTER')/sig_zhat,'k-','linewidth',2);hold on;
plot(xt,log(plot_ep_star_z_ZERO')/sig_zhat,'r--','linewidth',2);hold on;
plot(xt,log(plot_ep_star_z_NOINTER')/sig_zhat,'b','linewidth',2);hold on;
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Technology Growth Shock','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -4.5 2]); set(gca,'FontSize',12);  box off; grid on; 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

subplot(3,3,6)
plot(xt,log(plot_ep_star_d_INTER')/sig_dhat,'k-','linewidth',2);hold on;
plot(xt,log(plot_ep_star_d_ZERO')/sig_dhat,'r--','linewidth',2);hold on;
plot(xt,log(plot_ep_star_d_NOINTER')/sig_dhat,'b','linewidth',2);hold on;
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Discount Factor Shock','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -3 0]); set(gca,'FontSize',12);  box off; grid on; 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

subplot(3,3,7)
plot(xt,cumsum(plot_ep_star_dely_INTER)','k','linewidth',2);hold on;
plot(xt,cumsum(plot_ep_star_dely_NOINTER)','b','linewidth',2);
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Output (log-level)','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -6 1]); set(gca,'FontSize',12); box off; grid on;
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

subplot(3,3,8)
plot(xt,plot_ep_star_pi_INTER','k','linewidth',2);hold on;
plot(xt,plot_ep_star_pi_NOINTER','b','linewidth',2);
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Inflation (\% annualized)','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -1 4]);set(gca,'FontSize',12); box off; grid on;
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

subplot(3,3,9)
plot(xt,plot_ep_star_R_INTER','k','linewidth',2);hold on;
plot(xt,plot_ep_star_R_NOINTER','b','linewidth',2);
plot(xt,zeros(Tperiods,1)','k--','linewidth',1);
set(gca,'XTick',(1:14)); set(gca,'XTickLabel',mydates);
title('Interest Rate (\% annualized)','fontsize',16,'fontweight','bold','Interpreter','Latex')
axis([1 Tperiods -0.5 5]); set(gca,'FontSize',12); box off; grid on; 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
vline(7,'r-');

%---------------------------------------------------------------
% Figure A-1. Calibration of ARRA Fiscal Policy Intervenation 
%---------------------------------------------------------------

init_arra = init_use*0;
shocks_arra = zeros(4,12);
shocks_arra(2,1:length_intervention_fiscal) = size_spending_intervention;

zet_ss = 1-1/sParam.gstar;

zet_arra = [36 /3488 
    18/ 3533
8/ 3568 
24/ 3603 
26/ 3644 
21/ 3684
19/ 3704
20/ 3751
17/ 3791 
12/ 3830
9/ 3870 
8/ 3899];

zethat_arra = (zet_arra/zet_ss);

ghat_arra = (sParam.gstar-1)*zethat_arra;

[~,~,~,SIM_ARRA]            = fSimulatePLC(sModel,sParam,init_arra,shocks_arra);

figure(2); clf;
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
plot([1:12],ghat_arra(1:12),'b-.','linewidth',3); hold on;
plot([1:12],log(SIM_ARRA.g/sParam.gstar),'r-','linewidth',3);hold on;
set(gca,'XTick',(1:12)); set(gca,'XTickLabel',{'09Q2','','09Q4','','10Q2','','10Q4','','11Q2','','11Q4',''});
set(gca,'YTick',[0.002:0.003:0.016]); set(gca,'YTickLabel',{'','0.0050','0.0080','0.011','0.014'});
axis tight;  ylim([0.002 0.014]);
title('ARRA Comparison ($\hat{g}_t$)','fontsize',16,'fontweight','bold','Interpreter','Latex')
h_legend=legend('Received','Simulated','location','north','Orientation', 'Horizontal'); legend boxoff
set(h_legend,'FontSize',14); legend boxoff;
set(gca,'FontSize',14)


print(figure(1),'-dpdf','-r400',[resultsdir 'Figure_Shocks_and_ExPost_Intervention_v3.pdf'],'-fillpage');
print(figure(2),'-dpdf','-r400',[resultsdir 'Figure_ARRA_ghat.pdf']);
