# -----------------------------------------------------
# Inputs : params_in : model parameters
#          options_int_in : solutions options (integers)
#          options_float_in : solutions options (floats)
#
#outputs
#=
-theta_out: final solution theta
-time_out: Time of the PLC Solution
-flag_mat: PLC Solution flags, 1=clustering error, 2=Non-cconvergence in solution
=#
#
# Updated : 01/03/21
# -----------------------------------------------------

function fPLCSolve(params_in::Dict{String}{Float64},
    sPLC_in::Dict{String,Any},
    options_int_in::Dict{String}{Int64},
    options_float_in::Dict{String}{Float64},
    Smolyak_fixed_grid::Array{Float64,2}=zeros(1,1))

# -----------------------------------------------------
# Impose PLC restrictions model and parameter specific
# -----------------------------------------------------
sPLC_use = fPLCRestrictions(params_in,sPLC_in)

# ----------------------------------
# Initial Guess
# Solve log-linearized model
# ----------------------------------

# Compute derivatives of f (equilibrium conditions)
nfx, nfxp, nfy, nfyp, vcv =  fSGUGetOrder1(params_in);

#Compute first-order approximation
eflag_SGU, gx,hx = fSGU_gx_hx(nfy,nfx,nfyp,nfxp);

if eflag_SGU == 1

    # COLLECT COEFFICIENTS FROM LINEAR SOLUTION
    alpha0 = Vector{Any}(undef, sPLC["ny"][1])
    alpha0[1] = [gx[1,1];  0; gx[1,2:end]];        # pihat (1)
    alpha0[2] = [gx[2,1];  0; gx[2,2:end]];        # yhat  (2)

    # ARRANGE INITIAL GUESS (ONLY FREE COEFFICIENTS)
    theta_0 = [alpha0[1];        # pihat (nb)
              alpha0[2];        # yhat  (nb)
              alpha0[1][1];     # pihat (b)
              alpha0[2][1];     # yhat  (b)
              ];
    # --------------------------------------------------------------------------
    #                         Call Model Solutions
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # OPTIONS FOR SOLVER
    # --------------------------------------------------------------------------
    M = size(Smolyak_fixed_grid,1)

    tolx_in  = options_float_in["tolx"]
    tolf_in  = options_float_in["tolf"]
    tolg_in  = options_float_in["tolg"]
    maxIt_in = options_int_in["maxIt"]

    # --------------------------------------------------
    # Integration Rule (Monomial) [Order: EZ EG ED ER]
    # --------------------------------------------------
    n_nodes,nodes, weight = fMonomials_2_N4(vcv);

    EZ_nodes=nodes[:,1]
    EG_nodes=nodes[:,2]
    ED_nodes=nodes[:,3]
    ER_nodes=nodes[:,4]
    GH_weight=vec(weight)

    # Create Grid (Asymmetric Smolyak -- based on parameter values, stretch the fixed_grid)

    GRID = fScaleSmolyakGrid(Smolyak_fixed_grid,params_in,options_int_in,options_float_in)

    RHATLAGGRID = GRID[:,1];
    YHATLAGGRID = GRID[:,2];
    EGRID       = GRID[:,3];
    GGRID       = GRID[:,4];
    ZGRID       = GRID[:,5];
    DGRID       = GRID[:,6];

    # --------------------------------------------------------------------------
    # Minimization
    # --------------------------------------------------------------------------

    e_grid=vec(EGRID)
    g_grid=vec(GGRID)
    d_grid=vec(DGRID)
    z_grid=vec(ZGRID)
    Rhat_grid=vec(RHATLAGGRID)
    yhat_grid=vec(YHATLAGGRID)

    # ------------------------------------------------
    # ANONYMOUS FUNCTION THAT COMPUTES RESIDUALS AND JACOBIAN
    # ------------------------------------------------

    system(p::Array{Float64,1}) = fPLCSystemJacob(p,params_in,sPLC_use,
    e_grid,g_grid,z_grid,d_grid,Rhat_grid,yhat_grid,
    EG_nodes,ER_nodes,EZ_nodes,ED_nodes,GH_weight)

    flag = 0;

    # ---------------------------------
    # MINIMIZE OBJECTIVE FUNCTION
    # ---------------------------------

    theta_out,res_out,iter_c,conv,out_x,out_res,out_conv = fMyLevenbergMarquardt(system,theta_0,tolG=tolg_in,tolX=tolx_in,tolF=tolf_in,maxIter=maxIt_in,show_trace=true,lambda=1)

     if iter_c == maxIt_in
         flag = 0
     elseif ~conv
         flag = -1;
     else
         flag = 1;
     end

     return flag,theta_out,theta_0,res_out,sPLC_use

 else
     flag = -10;

     println("First-order solution failed, skipping PLC solution")

     return flag, nothing, nothing, nothing, nothing

 end

end
