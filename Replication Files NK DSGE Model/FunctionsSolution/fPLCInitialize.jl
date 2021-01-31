# ------------------------------------
# Function to initialize PLC matrices
# ------------------------------------

function fPLCInitialize(sModel)


### Change this script to spit out index_11.... to the SPLC structure
# HOUSEKEEPING AND MATRIX INITIALIZATION
nx       = sModel["nx"][1];
ny       = sModel["ny"][1];
n        = nx + 1;                   # Number of elements in state vector (including constant)
my       = (n+1)*ny;               # Number of free coefficients.
index_11 = 1;                                       # Index corresponding to alpha_11
index_12 = [2:n;]                                   # Index corresponding to alpha_12
index_21 = n+1;                                     # Index corresponding to beta_21
index_22 = [n+2:2*n;];                              # Index corresponding to beta_22
index_1  = [index_11; index_12];
index_2  = [index_21; index_22];

# INITIALIZE PLC MATRICES
# GAM{1}=gamma_1, GAM{2}=gamma_2, GAM{3}=gamma_y, GAM{4}=gamma_x
GAM = [zeros(1,1), zeros(n-1,1),zeros(ny,1),zeros(nx,1)]

#**** NOT SURE HOW TO TRANSLATE THIS

SMAT = Vector{Array{Float64,2}}(undef, ny)
for iy=1:ny
    SMAT[iy] = zeros(n+1,my);
    #SMAT[iy][1:n,1+n*(iy-1):n*iy] = I(n);
    SMAT[iy][1:n,1+n*(iy-1):n*iy] = Matrix(I,n,n); # Use Julia notation for Identity
    SMAT[iy][n+1,ny*n+iy] = 1;
end

PHIMAT = Vector{Array{Float64,2}}(undef, nx)
 for ix=1:nx
     PHIMAT[ix] = zeros(n+1,1);
 end

OMEGAMAT = Vector{Array{Float64,2}}(undef, nx)
 for ix=1:nx
     OMEGAMAT[ix] = zeros(n+1,my);
 end

 THETA     = zeros(my,1);

# Create dictionary
sPLC = Dict("GAM" => GAM, "SMAT" => SMAT , "PHIMAT" => PHIMAT, "OMEGAMAT" => OMEGAMAT,
"index_11" => index_11,"index_12"=> index_12,"index_21"=> index_21,
"index_22"=> index_22,"index_1"=> index_1,"index_2"=> index_2,
"n"=>n,"nx"=>nx,"ny"=>ny,"my"=>my);

return sPLC

end
