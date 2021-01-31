% This function checks that path of shocks does not generate a boundary
% problem with the exogenous state variables.
% The order of the shocks is:
% er, eg, ez, ed
function V_out = fn_get_shocks(Nrep,Tshocks,ind_no_er_shock)
global seedname

[s_er, s_eg, s_ez, s_ed] = RandStream.create(seedname,'NumStreams',4);

% Draw 10000 shocks for simulations and 100 rep %

count = 1;

while count<=Nrep
    
    er_use = randn(s_er, Tshocks, 1);
    ez_use = randn(s_ez, Tshocks, 1);
    eg_use = randn(s_eg, Tshocks, 1);
    ed_use = randn(s_ed, Tshocks, 1);
        

    V_out.ER(:,count) = er_use;
    V_out.EG(:,count) = eg_use;
    V_out.EZ(:,count) = ez_use;
    V_out.ED(:,count) = ed_use;
    
    if ind_no_er_shock == 1
    
    V_out.ER(:,count) = 0*er_use;
    
    end
    
    count=count+1;
    
end