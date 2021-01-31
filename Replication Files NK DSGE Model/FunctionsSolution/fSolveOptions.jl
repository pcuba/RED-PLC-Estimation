# Function to generate the solution options
# "int" for integers, "float" for real numbers
# First Version: 2018-10-25
# Current Version: 01/04/2021
# Version Notes:
#=

=#
function fSolveOptions()

options_int_out=Dict{String,Int64}()
options_float_out=Dict{String,Float64}()

# Solution Algorithm Options

options_float_out["tolx"]       = 1e-8  ; # Step tolerance
options_float_out["tolg"]       = 1e-8  ; # Gradient tolerance
options_float_out["tolf"]       = 1e-8  ; # Gradient tolerance

# PLC Options

options_int_out["maxIt"]        = 3000  ; # Maximum solve iterations

# Smolyak Grid Options

options_int_out["Smolyak_d"]          = 6;      # Dimension
options_int_out["Smolyak_mu"]         = 2;      # Parameter for grid
options_float_out["Smolyak_sd"]       = 1.282;  # The standard deviation multiple to use for stochastic variables (yields 90th percentile)

options_float_out["Smolyak_sd_min_d"]   = -4.0      # For the d grid, lower bound multiple for stand deviation
options_float_out["Smolyak_sd_max_d"]   = 1.282     # For the d grid, upper bound multiple for stand deviation

# Declare model variables
sModel= Dict{String,Int64}()
sModel["nf"]   = 10;             # Number of equilibrium conditions (include law of motions for exogenous states)
sModel["nx"]   = 6;              # Number of state variables (without constant)
sModel["ny"]   = 2;              # Number of control variables
sModel["nexo"] = 4;              # Number of exogenous variables

# Initialize PLC matrices
sPLC = fPLCInitialize(sModel)

return options_int_out,options_float_out,sPLC

end
