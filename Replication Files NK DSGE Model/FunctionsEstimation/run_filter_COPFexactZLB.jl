function run_filter_COPFexactZLB(ME_scale, nM, mhrun, silent, mData, nFileStrMH, nFileStrFilt,
                        rnd_off,Smolyak_fixed_grid,options_int_use,options_float_use,sPLC,ind_param)

    time_init_proc=time_ns()
    nPath = pwd()
    savedir = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/Filtered/"
    loaddir = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/"

    #--------------------------------------------------------------------------
    # Load prior specification
    #--------------------------------------------------------------------------
    include(nPath*"/PriorSpec/" * nNamePriorSpec * ".jl")

    #--------------------------------------------------------------------------
    # Configure data matrix
    #--------------------------------------------------------------------------
    vDate = mData[nStartDate:nEndDate,2]
    mY    = mData[nStartDate:nEndDate,[3,8,4,5]]
    # ygrowth (percent), 100*logcy, infl (annual %), R (annual %)
    println(" Start Date: $(vDate[1])")
    println(" End Date: $(vDate[end])")

    mSigmaMEtilde = [var(mY[:,1]) 0 0 0;
                     0 var(mY[:,2]) 0 0;
                     0 0 var(mY[:,3]) 0;
                     0 0 0 var(mY[:,4]) ];

    mSigmaME     = ME_scale*mSigmaMEtilde;

    pInvSigmaMEtilde = mSigmaMEtilde^(-1);


    #--------------------------------------------------------------------------
    # Load parameters
    #--------------------------------------------------------------------------

    mSummaryDraws = CSV.read(loaddir * nFileStrMH * "_Summary.csv", DataFrame; datarow=1)
    mSummaryDraws = convert(Array,mSummaryDraws)

    if ind_param == "MAP"

        vParameters   = mSummaryDraws[5,:]      # Use the MAP

    elseif ind_param == "Mean"

        vParameters   = mSummaryDraws[1,:]    # Use the posterior means

    end

    #--------------------------------------------------------------------------
    # Prior / Model Solution / State-space Representation / Likelihood
    # at initial parameter value
    #--------------------------------------------------------------------------

    time_init_loop = time_ns()

    # evaluate prior at initial value
    Random.seed!(rnd_off+1)
    nLogPriorDens = fPriorDensity(vParameters,mPriorSpec)

    # solve model and create state-space representation
    sParam            = fCreateParamDict(vParameters)
    flag_out, sParam_StateSpace = fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

    if flag_out < 1
        error("Model did not solve")
    end

    # initalize output matrix for filtered StateSpace
    Phi11                = sParam_StateSpace["Phi11"];
    ns                   = size(Phi11,1)+1 # for y(t-1)
    mFilteredStates      = Array{Float64,2}(undef, size(mY,1), ns)
    mFilteredStates_lag  = Array{Float64,2}(undef, size(mY,1), ns)
    mFilteredShocks      = Array{Float64,2}(undef, size(mY,1), 4)
    mReconstructedStates = Array{Float64,2}(undef, size(mY,1), ns-1) # drop y(t-1)
    vReconstructedRegime = Array{Float64,1}(undef, size(mY,1))

    # evaluate likelihood function
    Random.seed!(rnd_off+10+2)
    (lik, s_up_stat, s_lag_stat, Eta_stat, vRegime, Neff) = fCOPFaugExactZLB(mY, sParam_StateSpace, mSigmaMEtilde,ME_scale, nM, pInvSigmaMEtilde);
    nLogLiki = sum(lik)

    # s_lag_stat_i, shock_stat_i, regime_stat_i,


    # original: compute objective function at initial value
    nLogPost            = nLogLiki + nLogPriorDens

    # # experiment: compute likelihood at initial value
    # nLogPost            = nLogLiki
    # ############################

    mFilteredStates     = s_up_stat[:,:,2]
    mFilteredStates_lag = s_lag_stat[:,:,2]
    mFilteredEtas       = Eta_stat[:,:,2]

    # Unpack the Matrices
    vInitS     = sParam_StateSpace["vInitS"];
    pSigmaTE   = sParam_StateSpace["pSigmaTE"];
    pA         = sParam_StateSpace["pA"];
    pA0        = sParam_StateSpace["pA0"];
    D_EpsToEta = sParam_StateSpace["D_EpsToEta"]

    pSigmaTE  = sParam_StateSpace["pSigmaTE"];
    pSigmaEta = sParam_StateSpace["pSigmaTE"];

    eta1BarCoeff = sParam_StateSpace["eta1BarCoeff"];
    pPhi1_0      = sParam_StateSpace["Phi01"];
    Phi01        = pPhi1_0
    pPhi1_1      = sParam_StateSpace["Phi11"];
    Phi11        = pPhi1_1

    print("Eigenvalues Phi11")
    print(eigvals(Phi11))

    pPhi1_Eta    = sParam_StateSpace["PhiEta1"];
    PhiEta1      = pPhi1_Eta
    pPhi2_0      = sParam_StateSpace["Phi02"];
    Phi02        = pPhi2_0
    pPhi2_1      = sParam_StateSpace["Phi12"];
    Phi12        = pPhi2_1

    print("Eigenvalues Phi12")
    print(eigvals(Phi12))
    pPhi2_Eta    = sParam_StateSpace["PhiEta2"];
    PhiEta2      = pPhi2_Eta

    # convert etas into eps
    mFilteredEps = (inv(D_EpsToEta) * mFilteredEtas')'

    # check whether we can reproduce the states by forward simulation with epsilon
    mReconstructedStates[1,:] = mFilteredStates[1,1:end-1]
    vReconstructedRegime[1]   = vRegime[1]

    for iTime = 2:size(mY,1)
        if mFilteredEtas[iTime,1] <= (eta1BarCoeff*[1;mFilteredStates_lag[iTime,1:end-1] ])[1,1]
           vReconstructedRegime[iTime] = 1
        else
           vReconstructedRegime[iTime] = 2
        end

        if vRegime[iTime] == 1
           mReconstructedStates[iTime,:] = (pPhi1_0 + pPhi1_1*mFilteredStates_lag[iTime,1:end-1]
           + pPhi1_Eta*D_EpsToEta*mFilteredEps[iTime,:])
        else
           mReconstructedStates[iTime,:] = (pPhi2_0 + pPhi2_1*mFilteredStates_lag[iTime,1:end-1]
           + pPhi2_Eta*D_EpsToEta*mFilteredEps[iTime,:])
        end

    end

    mReconstructedStates = (mReconstructedStates - mFilteredStates[:,1:end-1])./abs.(mFilteredStates[:,1:end-1])

    # Create data comparable series
    mObsOut = (pA0 .+  pA*mFilteredStates')'

    # Create output level
    vOutput = mObsOut[:,1]
    for iTime = 1:size(vDate)[1]
        vOutput[iTime] = pA0[1]*iTime + 100*sum(mFilteredStates[1:iTime,3]) + 100*sum(mFilteredStates[1:iTime,5]) - 100*sum(mFilteredStates[1:iTime,9])
    end

    mObsOut = [mObsOut vOutput]

    vIntRateLag = [NaN ; (mY[1:end-1,4] .-pA0[4])./400]

    mXStates = [vIntRateLag ones(size(vIntRateLag)[1]) mFilteredStates[:,9] mFilteredStates[:,3] mFilteredStates[:,4] mFilteredStates[:,2] mFilteredStates[:,1]  ]

    savename     = nFileStrFilt * "_XStates" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mXStates), header=false)

    savename     = nFileStrFilt * "_ObsOut" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mObsOut), header=false)

    savename     = nFileStrFilt * "_Param" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mParameters=vParameters), header=false)

    savename     = nFileStrFilt * "_Data" * ".csv"
    CSV.write(savedir * savename,  DataFrame([mY vDate], :auto), header=false)

    savename     = nFileStrFilt * "_FilteredStates" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mFilteredStates, :auto), header=false)

    savename     = nFileStrFilt * "_FilteredEtas" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mFilteredEtas, :auto), header=false)

    savename     = nFileStrFilt * "_FilteredEps" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mFilteredEps, :auto), header=false)

    savename     = nFileStrFilt * "_FilteredStates_lag" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mFilteredStates_lag, :auto), header=false)

    savename     = nFileStrFilt * "_ReconstructedStates" * ".csv"
    CSV.write(savedir * savename,  DataFrame(mReconstructedStates, :auto), header=false)

    savename     = nFileStrFilt * "_Regime" * ".csv"
    CSV.write(savedir * savename,  DataFrame([vRegime vReconstructedRegime], :auto), header=false)

    savename     = nFileStrFilt * "_LogPosterior" * ".csv"
    CSV.write(savedir * savename,  DataFrame(LogPost = nLogPost), header=false)

    # Phi Matrices
    savename     = nFileStrFilt * "_Param_Phi01" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["Phi01"], :auto), header=false)
    savename     = nFileStrFilt * "_Param_Phi11" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["Phi11"], :auto), header=false)
    savename     = nFileStrFilt * "_Param_PhiEta1" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["PhiEta1"], :auto), header=false)

    savename     = nFileStrFilt * "_Param_Phi02" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["Phi02"], :auto), header=false)
    savename     = nFileStrFilt * "_Param_Phi12" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["Phi12"], :auto), header=false)
    savename     = nFileStrFilt * "_Param_PhiEta2" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["PhiEta2"], :auto), header=false)

    # Zeta Function
    savename     = nFileStrFilt * "_Param_eta1BarCoeff" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["eta1BarCoeff"], :auto), header=false)

    # Zeta Function
    savename     = nFileStrFilt * "_Param_D_EpsToEta" * ".csv"
    CSV.write(savedir * savename,  DataFrame(sParam_StateSpace["D_EpsToEta"], :auto), header=false)




end
