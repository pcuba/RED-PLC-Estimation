%# Function to create a dictionary with parameter values
function ParaDict_out = fCreateParams(vPara)

vPara = csvread([vPara 'Param.csv']);

ParaDict_out.tau     = vPara(1);
ParaDict_out.kappa   = vPara(2);
ParaDict_out.psi1    = vPara(3);
ParaDict_out.psi2    = vPara(4);
ParaDict_out.rho_r   = vPara(5);
ParaDict_out.rho_g   = vPara(6);
ParaDict_out.rho_d   = vPara(7);
ParaDict_out.rho_z   = vPara(8);
ParaDict_out.sig_r   = vPara(9);
ParaDict_out.sig_g   = vPara(10);
ParaDict_out.sig_d   = vPara(11);
ParaDict_out.sig_z   = vPara(12);
ParaDict_out.eta     = vPara(13);
ParaDict_out.nu      = vPara(14);
ParaDict_out.chi_h   = vPara(15);
ParaDict_out.gstar   = vPara(16);
ParaDict_out.rAnet   = vPara(17);
ParaDict_out.gamQnet = vPara(18);
ParaDict_out.piAnet  = vPara(19);
ParaDict_out.pibarAnet = vPara(20);


% #--------------------------------------------------------------------------
% #                       TRANSFORMED PARAMETERS
% #--------------------------------------------------------------------------
piAnet  = ParaDict_out.piAnet;
rAnet   = ParaDict_out.rAnet;
gamQnet = ParaDict_out.gamQnet;
tau     = ParaDict_out.tau;
nu      = ParaDict_out.nu;
kappa   = ParaDict_out.kappa;
chi_h   = ParaDict_out.chi_h;
gstar   = ParaDict_out.gstar;
eta     = ParaDict_out.eta;


pistar= exp(piAnet/400);                      % Annualized inflation rate
pibar = pistar;                               % Reference inflation rate for cost adjustment
beta  = exp(-rAnet/400);                      % Disccount factor
gamma = exp(gamQnet/100);                     % gross growth rate
r     = gamma/beta;                           % Real interest rate
phi   = tau*(1-nu)/((nu * pistar^2 * kappa)); % Adjustment cost parameter
b     = 1/(2*nu);                             % Constant

ParaDict_out.pistar = pistar;
ParaDict_out.pibar  = pibar;
ParaDict_out.beta   = beta;
ParaDict_out.gamma  = gamma;
ParaDict_out.r      = r;
ParaDict_out.phi    = phi;
ParaDict_out.b      = b;

% #--------------------------------------------------------------------------
% #                   NON-STOCHASTIC STEADY STATE
% #--------------------------------------------------------------------------
pi_ss  = pistar;
c_ss   = ( (1 - nu + phi*nu*(1-beta)*pi_ss*(pi_ss-pibar)-0.5*phi*(pi_ss-pibar)^2) / (chi_h*((1/(gstar))-0.5*phi*(pi_ss-pibar)^2)^(-1/eta)) )^(1/(tau + 1/eta));
y_ss   = c_ss/(1/(gstar) - 0.5 * phi * (pi_ss -pibar)^2);
R_ss   = pi_ss * gamma/beta;
ee1_ss = (c_ss^-tau)/(gamma*pi_ss);

ParaDict_out.pi_ss = pi_ss;
ParaDict_out.c_ss  = c_ss;
ParaDict_out.y_ss  = y_ss;
ParaDict_out.R_ss  = R_ss;
ParaDict_out.ee1_ss= ee1_ss;
ParaDict_out.cy    = c_ss/y_ss;
