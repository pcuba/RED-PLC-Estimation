# Function to create a dictionary with parameter values
function fCreateParamDict(vPara)

ParaDict_out=Dict{String,Float64}()

ParaDict_out["tau"]       = vPara[1];
ParaDict_out["kappa"]     = vPara[2];
ParaDict_out["psi1"]      = vPara[3];
ParaDict_out["psi2"]      = vPara[4];
ParaDict_out["rho_r"]     = vPara[5];
ParaDict_out["rho_g"]     = vPara[6];
ParaDict_out["rho_d"]     = vPara[7];
ParaDict_out["rho_z"]     = vPara[8];
ParaDict_out["sig_r"]     = vPara[9];
ParaDict_out["sig_g"]     = vPara[10];
ParaDict_out["sig_d"]     = vPara[11];
ParaDict_out["sig_z"]     = vPara[12];
ParaDict_out["eta"]       = vPara[13];
ParaDict_out["nu"]        = vPara[14];
ParaDict_out["chi_h"]     = vPara[15];
ParaDict_out["gstar"]     = vPara[16];
ParaDict_out["rAnet"]     = vPara[17];
ParaDict_out["gamQnet"]   = vPara[18];
ParaDict_out["piAnet"]    = vPara[19];
ParaDict_out["pibarAnet"] = vPara[20];
ParaDict_out["er0"]       = vPara[21];
ParaDict_out["g0"]        = vPara[22];
ParaDict_out["z0"]        = vPara[23];
ParaDict_out["d0"]        = vPara[24];
ParaDict_out["y0"]        = vPara[25];
ParaDict_out["pi0"]       = vPara[26];
ParaDict_out["c0"]        = vPara[27];
ParaDict_out["R0star"]    = vPara[28];
ParaDict_out["ylag0"]     = vPara[29];

return ParaDict_out

end
