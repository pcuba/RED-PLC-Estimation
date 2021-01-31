% This function checks that path of shocks does not generate a boundary
% problem with the exogenous state variables.

% Used with the forward guidance code 

function V_out = fn_get_shocks_fg(Nrep,Tshocks,init)
global seedname
global size_shock 

[s_er s_ez s_eg] = RandStream.create(seedname,'NumStreams',3);

% Draw 10000 shocks for simulations and 100 rep %

count = 1;

while count<=Nrep
    
    init_use.zlag = init.zlag;
    init_use.glag = init.glag;
    
    er_use = sig_r*randn(s_er, Tshocks, 1);
    ez_use = sig_z*randn(s_ez, Tshocks, 1);
    eg_use = sig_g*randn(s_eg, Tshocks, 1);
    
    er_keep = er_use;
    ez_keep = ez_use;
    eg_keep = eg_use;
    
    if abs(size_shock) > 1
        er_use(1) = er_use(1) + size_shock*sig_r;
        ez_use(1) = ez_use(1) + size_shock*sig_z;
        eg_use(1) = eg_use(1) + size_shock*sig_g;
    end
    
    zs = zeros(Tshocks,1);
    gs = zeros(Tshocks,1);
    ers = zeros(Tshocks,1);
    
    for i=1:Tshocks
        
        if i == 1
            z_lag = init_use.zlag; %0
            g_lag = init_use.glag; %log(gstar);
        else
            z_lag = zs(i-1);
            g_lag = gs(i-1);
        end
        
        
        zs(i) = rho_z*z_lag + ez_use(i);
        gs(i) = (1-rho_g)*log(gstar) + rho_g*g_lag + eg_use(i);
        ers(i)   = er_use(i);
    end
    
    % Now check that simulation of g is inside R_min R_max
    emax_temp = max(ers);
    emin_temp = min(ers);

    gmax_temp = max(gs);
    gmin_temp = min(gs);

    zmax_temp = max(zs);
    zmin_temp = min(zs);

    
    if gmax_temp<=g_max && gmin_temp > g_min ...
       && zmax_temp<=z_max && zmin_temp > z_min ...
       && emax_temp<=e_max && emin_temp > e_min

        V_out.ER(:,count) = er_keep; 
        V_out.EZ(:,count) = ez_keep; 
        V_out.EG(:,count) = eg_keep; 
        
        count=count+1;    
    else
        fprintf('\n Discarding simulation at count = %i... \n',count);
    end
    
end