# Piecewise-linear Approximations and Filtering for DSGE Models with Occasionally-Binding Constraints

Replication Files for: 

"Piecewise-linear Approximations and Filtering for DSGE Models with Occasionally-Binding Constraints".
Review of Economic Dynamics. (Forthcoming)

by

[Boragan Aruoba](http://aruoba.econ.umd.edu) (Maryland), 
[Pablo Cuba-Borda](http://www.pcubaborda.net) (Board of Governors), 
[Kenji Higa-Flores]() (Maryland), 
[Sergio Villalvazo](https://www.sergiovillalvazo.com) (Penn), and 
[Frank Schorfheide](http://sites.sas.upenn.edu/schorf/) (Penn)

## List of Scripts
These scripts, run in this order, generate the results presented in the paper. All results are produced using Julia v1.5.3. 
Plots are produced using the Matlab R2020a scripts in this directory.

1. __Generate data__
* script_SimulateModel.jl

2. __Figure 1__
* script_LikelihoodEval.jl: nSimEval=100, nParamDraw=1, run for each ME scale separately
* Graph_LoglhApprox_Densities.m: run for each ME scale separately 

3. __Figure 2__
script_LikelihoodEval.jl: nSimEval=100, nParamDraw=100, run for each ME scale separately
Graph_LoglhApprox_Stddev.m: run for each ME separately

4. __Figure 3__
* script_PriorDrawsDirect.jl to generate draws from the prior for proposal covariance matrix
* script_MCMC_KF.jl to generate draws from the posterior of the linearized DSGE using the Kalman Filter (KF) to evaluate the likelihood function; mhrun=0 with proposal = prior covariance matrix, initial parameters from the true DGP and c=0.002. 110,000 draws, 10,000 dropped. 
* script_SummarizeDraws to compute a posterior covariance matrix based on the KF run; need to specify input file with nFilter, nPrior, nMHrun, nDataSet
* script_MCMC_KF.jl to generate draws from the posterior of the linearized DSGE using the Kalman Filter (KF) to evaluate the likelihood function; mhrun=1 with proposal = posterior covariance matrix from previous run , initial parameters from the posterior mean of the previous run and c=0.3. 110,000 draws, 10,000 dropped. 
* script_SummarizeDraws to compute a posterior covariance matrix based on the KF run; need to specify input file with nFilter, nPrior, nMHrun, nDataSet
* script_MCMC_BSPF.jl to generate draws from the posterior of the nonlinear DSGE using the BSPF to evaluate the likelihood function; mhrun = 2 with proposal = posterior VCV and initial parameters from posterior mean from KF mhrun=1.
* script_MCMC_COPF.jl to generate draws from the posterior of the nonlinear DSGE using the BSPF to evaluate the likelihood function; mhrun = 2 with proposal = posterior VCV and initial parameters from posterior mean from KF mhrun=1.
* Top row: Graph_PostDraws_Densities_Compare.m
* Middle row: Graph_PostDraws_ACFs_Compare.m
* Bottom row: Graph_PostDraws_ACF_Scatters.m

5. __Table 2__
* script_PriorDrawsDirect.jl to generate draws from the prior for proposal covariance matrix [Already one for Figure 3]
* script_MCMC_KF.jl to generate draws from the posterior of the linearized DSGE using the Kalman Filter (KF) to evaluate the likelihood function; mhrun=0 with proposal = prior covariance matrix and c=0.002. (110,000 draws, 10,000 dropped)
* script_SummarizeDraws to compute a posterior covariance matrix based on the KF run; need to specify input file with nFilter, nPrior, nMHrun, nDataSet
* script_MCMC_KF.jl to generate draws from the posterior of the linearized DSGE using the Kalman Filter (KF) to evaluate the likelihood function; mhrun=1 with proposal = posterior covariance matrix from previous run, initial parameters are the posterior mean from previous run and c=0.3 (110,000 draws, 10,000 dropped) 
* script_SummarizeDraws to compute a posterior covariance matrix based on the KF run; need to specify input file with nFilter, nPrior, nMHrun, nDataSet
* script_MCMC_COPFexactZLB to generate draws from the posterior of the nonlinear DSGE model using the COPF to evaluate the likelihood function; M=180, c=0.2, ME=0.001, 55,000 draws, discard 5000.
* script_SummarizeDraws to compute a posterior covariance matrix based on the COPF run; need to specify input file with nFilter, nPrior, nMHrun, nDataSet

6. __Figure 4__
*  script_Filter_COPFexactZLB.jl to generate filtered states conditional on the MAP or posterior mean estimates.
*  script_ExPost_Analysis.m to produce intervention vs no-intervention paths; need to specify directory containing filtered states conditional on posterior mean or MAP, point to output directory from script_Filter_COPFexactZLB.jl

7. __Table 3__
* script_ExPost_FiscalMultiplier.m compute Ex-Post fiscal multiplier and interest path under alternative; need to specify directory containing filtered states conditional on posterior mean or MAP, point to output directory from script_Filter_COPFexactZLB.jl
 
