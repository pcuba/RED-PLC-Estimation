# Calls Smolyak_Elem_Isotrop.jl and Smolyak_Grid.jl to produce a Smolyak GRID
# Scales the dimensions given an input vector.

function fScaleSmolyakGrid(fixed_grid,params_use,options_int_in,options_float_in)

    d      = options_int_in["Smolyak_d"]
    mu     = options_int_in["Smolyak_mu"]
    sd     = options_float_in["Smolyak_sd"]
    R_quantile = options_float_in["R_quantile"]
    d_min_norm = options_float_in["Smolyak_sd_min_d"]
    d_max_norm = options_float_in["Smolyak_sd_max_d"]

    Smol_grid = zeros(size(fixed_grid))

    R_ss=params_use["R_ss"]
    R_min=log(1/R_ss)
    R_max=log(R_quantile/R_ss)

    rho_z  = params_use["rho_z"];
    rho_g  = params_use["rho_g"];
    rho_r  = params_use["rho_r"];
    rho_d  = params_use["rho_d"];

    sig_z  = params_use["sig_z"];
    sig_g  = params_use["sig_g"];
    sig_r  = params_use["sig_r"];
    sig_d  = params_use["sig_d"];

    z_min=-sd*sig_z/sqrt(1-rho_z^2)
    z_max=sd*sig_z/sqrt(1-rho_z^2)

    d_min=d_min_norm*sig_d/sqrt(1-rho_d^2)
    d_max=d_max_norm*sig_d/sqrt(1-rho_d^2)

    g_min=-sd*sig_g/sqrt(1-rho_g^2)
    g_max=sd*sig_g/sqrt(1-rho_g^2)

    er_min=-sd*sig_r
    er_max=sd*sig_r

    y_max=z_max;
    y_min=z_min;

    max_min=[R_max R_min ;y_max y_min; er_max er_min; g_max g_min; z_max  z_min; d_max  d_min]

    #Asymetrically adjust R dimension
    for j=1:size(Smol_grid,1)
    if fixed_grid[j,1]<=0
    Smol_grid[j,1]=(fixed_grid[j,1]+1)*(0-R_min)+R_min
    else
    Smol_grid[j,1]=(fixed_grid[j,1]+0)*(R_max-0)+0
    end
    end

    for i=2:Int(d-1)
        max_t=max_min[i,1]
        min_t=max_min[i,2]
        Smol_grid[:,i]=(fixed_grid[:,i].+1)/2 .*(max_t-min_t) .+min_t
    end

    for j=1:size(Smol_grid,1)
    if fixed_grid[j,6]<=0
    Smol_grid[j,6]=(fixed_grid[j,6]+1)*(0-d_min)+d_min
    else
    Smol_grid[j,6]=(fixed_grid[j,6]+0)*(d_max-0)+0
    end
    end

return Smol_grid

end
