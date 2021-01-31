function fPLCSystemJacob(theta_in::Array{Float64,1}, par::Dict{String}{Float64},
sPLC::Dict{String,Any},e_grid::Array{Float64,1},g_grid::Array{Float64,1},z_grid::Array{Float64,1},
d_grid::Array{Float64,1},Rhat_grid::Array{Float64,1},yhat_grid::Array{Float64,1},
EG_nodes::Array{Float64,1},ER_nodes::Array{Float64,1},EZ_nodes::Array{Float64,1},
ED_nodes::Array{Float64,1},weight::Array{Float64,1})

# ---------------------------------------
# Map parameters
# ---------------------------------------
tau   = par["tau"];
nu    = par["nu"];
chi_h = par["chi_h"];
phi   = par["phi"];
pibar = par["pibar"];
pist  = par["pi_ss"];
gst   = par["gstar"];
rho_g = par["rho_g"];
rho_z = par["rho_z"];
rho_d = par["rho_d"];
bet   = par["beta"];
eta   = par["eta"];
inveta= 1.0/eta;
y_ss  = par["y_ss"];
c_ss  = par["c_ss"];
cy    = par["cy"];

# ---------------------------------------
#  Housekeeping
# ---------------------------------------

# Allocate memory for objects

X2 = ones(Float64,6)
X2_prime = ones(Float64,6)
state_i_prime = ones(Float64,7)
state_i = ones(Float64,7)

nGrid = length(e_grid);

res1 = zeros(nGrid);
res2 = zeros(nGrid);

# Allocate memory for integrands

# Integreation nodes
nNodes = length(EZ_nodes);

integrand_1 = zeros(nNodes,1);
integrand_2 = zeros(nNodes,1);

# -----------------------------------
#  Jacobian
# -----------------------------------

JACOB = zeros(16*2,nGrid);
integrand_jacob = zeros(16*2,nNodes);

# ---------------------------------------
# Map coefficients and PLC objecs
# ---------------------------------------

ALP, BET, delta,dALP22,dBET22 = fPLCCoefs(theta_in,sPLC)

SMAT=sPLC["SMAT"];
SMAT_1=SMAT[1];
SMAT_2=SMAT[2];

index_1= sPLC["index_1"];
index_2= sPLC["index_2"];

SMAT_1_1=SMAT[1][index_1,:];
SMAT_2_1=SMAT[2][index_1,:];

index_states_sgu_1::Array{Int64,1}    = vcat(index_1[1], index_1[3:end]) ;
index_states_sgu_2::Array{Int64,1}    = vcat(index_2[1], index_2[3:end]) ;

# MAIN LOOP
for i=1:nGrid

    # States
    er::Float64   = e_grid[i];
    g::Float64    = g_grid[i];
    z::Float64    = z_grid[i];
    d::Float64    = d_grid[i];
    Rhat0::Float64= Rhat_grid[i];
    yhat0::Float64= yhat_grid[i];

    # State partition
    x1 = Rhat0;
    X2[2] = yhat0;
    X2[3] = z;
    X2[4] = g;
    X2[5] = d;
    X2[6] = er;

    # Basis function
    state_i[1] = Rhat0;
    state_i[3] = yhat0;
    state_i[4] = z;
    state_i[5] = g;
    state_i[6] = d;
    state_i[7] = er;

    # Check ZLB
    if x1>dot(delta,X2)

        # Approx. Policies at current point
        pihat = dot(state_i,ALP[1][index_1]);
        yhat  = dot(state_i,ALP[2][index_1]);
        Rhat  = dot(state_i,BET[1][index_1]);
        bind  = 0;

    else

        # Approx. Policies at current point
        pihat = dot(state_i,ALP[1][index_2]);
        yhat  = dot(state_i,ALP[2][index_2]);
        Rhat  = dot(state_i,BET[1][index_2]);
        bind  = 1;

    end

    chat  = yhat - g;

    # Compute endogenous objects
    ff = ( 1 - nu -chi_h*(c_ss^tau)*(y_ss^inveta)*exp(tau*chat+inveta*yhat) + nu*phi*pist*(pist*exp(pihat)-pibar)*exp(pihat)) - 0.5*phi*(pist*exp(pihat) - pibar)^2;   # OK

    # Compute elements for derivative
    dXdTheta,dYdTheta,dYprime_dXprime_1,dYprime_dXprime_2 = fPLCPreJacobian(ALP,delta,dALP22,dBET22,par,sPLC,state_i,index_states_sgu_1,index_states_sgu_2,SMAT_1_1,SMAT_2_1)

    # Integrals
    for q=1:nNodes

        # Get all primes at current node
        g_prime::Float64  =  rho_g*g + EG_nodes[q];
        z_prime::Float64  =  rho_z*z + EZ_nodes[q];
        d_prime::Float64  =  rho_d*d + ED_nodes[q];
        er_prime::Float64 =  ER_nodes[q];

        # State partition
        x1_prime = Rhat;
        X2_prime[2] = yhat;
        X2_prime[3] = z_prime;
        X2_prime[4] = g_prime;
        X2_prime[5] = d_prime;
        X2_prime[6] = er_prime;

        # Basis function
        state_i_prime[1] = Rhat;
        state_i_prime[3] = yhat;
        state_i_prime[4] = z_prime;
        state_i_prime[5] = g_prime;
        state_i_prime[6] = d_prime;
        state_i_prime[7] = er_prime;

        # Check ZLB at t+1
        if x1_prime > dot(delta,X2_prime)

            # Approximate t+1 policies
            pihat_prime = dot(state_i_prime,ALP[1][index_1]);
            yhat_prime  = dot(state_i_prime,ALP[2][index_1]);

        else

            # Approximate t+1 policies
            pihat_prime = dot(state_i_prime,ALP[1][index_2]);
            yhat_prime  = dot(state_i_prime,ALP[2][index_2]);

        end

        chat_prime  = yhat_prime - g_prime;

        # Compute integrands at current node
        @inbounds integrand_1[q] =  exp(Rhat - pihat_prime - z_prime + d_prime - d + tau*chat - tau*chat_prime);  # OK
        @inbounds integrand_2[q] =  exp(-tau*chat_prime + tau*chat + d_prime - d + yhat_prime - yhat + pihat_prime)*(pist*exp(pihat_prime)-pibar);  # OK 

        @inbounds integrand_jacob[:,q] = fPLCGetJacobian(SMAT,ALP,delta,dALP22,dBET22,par,sPLC,state_i,state_i_prime,dXdTheta,dYdTheta,index_states_sgu_1,index_states_sgu_2,pihat,yhat,pihat_prime,yhat_prime,dYprime_dXprime_1,dYprime_dXprime_2,SMAT_1,SMAT_1_1,SMAT_2,SMAT_2_1)

    end

    # Computing Expectations
    sol_int1 = dot(integrand_1,weight);
    sol_int2 = dot(integrand_2,weight);
    sol_int_jac = integrand_jacob*weight;

    # Residual Function
    res1[i] = 1  - sol_int1;  # OK
    res2[i] = ff - nu*phi*bet*pist*sol_int2;  # OK

    # Jacobian
    JACOB[:,i] = sol_int_jac

end

# Residuals to be minimized

G = [res1;res2];

# Reshape Jacobian

JACOUT = reshape(JACOB',nGrid*2,16)

return G, JACOUT

end
