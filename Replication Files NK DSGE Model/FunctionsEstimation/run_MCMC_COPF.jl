function run_MCMC_COPF(ME_scale, nM, nSimRWMH, mhrun, nBurnRWMH,
           silent, nOutput, c, nNamePriorSpec, mData, nStartDate, nEndDate,
           nFileStr, nFileStrMH, vParameters, mSigmaProp, rnd_off,
           Smolyak_fixed_grid,options_int_use,options_float_use,sPLC)

    time_init_proc=time_ns()
    nPath = pwd()
    savedir = "$(pwd())/PosteriorDraws/" * nFileStrMH * "/"

    #--------------------------------------------------------------------------
    # Load prior specification
    #--------------------------------------------------------------------------
    include(nPath * "/PriorSpec/" * nNamePriorSpec * ".jl")

    #--------------------------------------------------------------------------
    # Configure data set and measurement errors
    #--------------------------------------------------------------------------

    vDate = mData[nStartDate:nEndDate,2]
    mY    = mData[nStartDate:nEndDate,[3,8,4,5]]


    # ygrowth (percent), log (gy), infl (annual %), R (annual %)
    println(" Start Date: $(vDate[1])")
    println(" End Date: $(vDate[end])")

    mSigmaMEtilde = [var(mY[:,1]) 0 0 0;
                     0 var(mY[:,2]) 0 0;
                     0 0 var(mY[:,3]) 0;
                     0 0 0 var(mY[:,4]) ];

    pInvSigmaMEtilde = mSigmaMEtilde^(-1);

    mSigmaME     = ME_scale*mSigmaMEtilde;
    #--------------------------------------------------------------------------
    # Configure proposal covariance matrix
    #--------------------------------------------------------------------------

    # set values in covariance matrix to zero that correspond
    # to fixed parameters
    vMask      = mPriorSpec[:,4]
    vFix       = mPriorSpec[:,5]
    vMaskInv   = ones(size(vMask)) - vMask
    mSigmaProp = mSigmaProp .* (vMaskInv*vMaskInv')

    #--------------------------------------------------------------------------
    # Initialize RWMH
    #--------------------------------------------------------------------------

    # Initial simulation of RWMH
    # initialize output matrices
    nStructParam = size(vParameters,1)
    mStructParam           = Array{Float64,2}(undef,nSimRWMH, nStructParam)
    vLogPosterior          = Array{Float64,2}(undef,nSimRWMH,2)
    vLogPosteriorCandidate = Array{Float64,2}(undef,nSimRWMH,2)
    mStructParamCandidate  = Array{Float64,2}(undef,nSimRWMH, nStructParam)

    #--------------------------------------------------------------------------
    # Prior / Model Solution / State-space Representation / Likelihood
    # at initial parameter value
    #--------------------------------------------------------------------------

    iSim = 1
    time_init_loop = time_ns()

    # evaluate prior at initial value
    Random.seed!(rnd_off+1)
    nLogPriorDens = fPriorDensity(vParameters,mPriorSpec)

    # solve model and create state-space representation
    sParam = fCreateParamDict(vParameters)
    flag_out, sParam_StateSpace = fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

    if flag_out < 1
        error("Failed to solve the model at the initial parameters")
    end

    # compute max eigenvalues
    # Phi11 = sParam_StateSpace["Phi11"];
    # nMaxEv11 = maximum(abs.(eigvals(Phi11)))
    # Phi12 = sParam_StateSpace["Phi12"];
    # nMaxEv12 = maximum(abs.(eigvals(Phi12)))

    # evaluate likelihood function
    Random.seed!(rnd_off+iSim*10+2)
    (lik_i, s_up_stat_i, s_lag_stat_i, shock_stat_i, regime_stat_i, Neff) = fCOPF(mY, sParam_StateSpace, mSigmaMEtilde, ME_scale, nM, pInvSigmaMEtilde);
    nLogLiki = sum(lik_i)


    # compute objective function at initial value
    nObj          = nLogLiki + nLogPriorDens

    time_loop = signed(time_ns()-time_init_loop)/1000000000
    println("")
    println(" Post at Initial Para:  $nObj")
    println(" Elapsed time: $time_loop")


    #--------------------------------------------------------------------------
    # Initialization of Loop
    #--------------------------------------------------------------------------

    # initializations for loop
    nCount      = 0.0;
    nAccept     = 0.0;
    nAcceptRate = 0.0;
    nAlpha      = 0.0;
    mStructParam[1,:]          = vParameters'
    mStructParamCandidate[1,:] = vParameters'
    vLogPosterior[1,1]         = nObj;
    vLogPosterior[1,2]         = 1;
    # vLogPosterior[1,3]         = nMaxEv11;
    # vLogPosterior[1,4]         = nMaxEv12;
    vLogPosteriorCandidate[1,1]  = nObj;
    vLogPosteriorCandidate[1,2]= 1;
    # vLogPosteriorCandidate[1,3]= nMaxEv11;
    # vLogPosteriorCandidate[1,4]= nMaxEv12;

    # sqrt of proposal covariance Matrix
    mPropVarSVD   = svd(Matrix(c*mSigmaProp))
    mPropVar_U    = mPropVarSVD.U
    mPropVar_S    = mPropVarSVD.S
    mPropVar_sqrt = mPropVar_U * ( Matrix(1.0I, nStructParam,nStructParam).* (mPropVar_S.^(1/2)) )

    #--------------------------------------------------------------------------
    # RWMH Loop
    #--------------------------------------------------------------------------

    counter = 1
    if silent == 0
       println("")
       println(" Bayesian Computations: RWMH... ")
       println("")
    end
    time_init_loop=time_ns()

    for iSim = 2:nSimRWMH

        # Draw a THETA candidate
        Random.seed!(rnd_off+iSim*10+0)
        nCandidateShock = mPropVar_sqrt*Random.rand(MvNormal(zeros(nStructParam),Matrix(1.0I, nStructParam,nStructParam)))
        vParametersCandidate = mStructParam[iSim-1,:] + nCandidateShock
        # make sure parameter vector is consistent w fixed values
        vParametersCandidate = vMaskInv.*vParametersCandidate + vMask.*vFix

        nParametersCandidateValid = ( (vParametersCandidate[2]<=1)*(vParametersCandidate[2]>=0)*(vParametersCandidate[5]<=1)*(vParametersCandidate[6]<=1)*
            (vParametersCandidate[7]<=1)*(vParametersCandidate[8]<=1)*(vParametersCandidate[9]>0)*
            (vParametersCandidate[10]>0)*(vParametersCandidate[11]>0)*(vParametersCandidate[12]>0) )

        if  nParametersCandidateValid == 1

            # Solve model
            # Random.seed!(rnd_off+iSim*10+1)
            Random.seed!(rnd_off+1)
            sParam = fCreateParamDict(vParametersCandidate)

            flag_out, sParam_StateSpace = fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

            if flag_out < 1
                nParametersCandidateValid = 0
            end
        end

        if  nParametersCandidateValid == 1

            # Evaluate prior at proposed draw
            nLogPriorDensCandidate = fPriorDensity(vParametersCandidate,mPriorSpec)

            # compute max eigenvalues
            # Phi11 = sParam_StateSpace["Phi11"];
            # nMaxEv11 = maximum(abs.(eigvals(Phi11)))
            # Phi12 = sParam_StateSpace["Phi12"];
            # nMaxEv12 = maximum(abs.(eigvals(Phi12)))

            # Evaluate likelihood function
            #println("*** Evaluate Likelihood Fcn ")
            Random.seed!(rnd_off+iSim*10+2)
            (lik_i, s_up_stat_i, s_lag_stat_i, shock_stat_i, regime_stat_i, Neff) = fCOPF(mY, sParam_StateSpace, mSigmaMEtilde, ME_scale, nM,pInvSigmaMEtilde);
            nLogLikiCandidate = sum(lik_i)

            # objective function
            nObjCandidate = nLogPriorDensCandidate + nLogLikiCandidate
            nAlpha = min(1,exp(nObjCandidate-nObj))

        else # candidate of THETA outside of bounds
            nLogPriorDensCandidate = -1E20
            nLogLikiCandidate      = -1E20
            nObjCandidate          = nLogPriorDensCandidate + nLogLikiCandidate
            nAlpha                 = 0
        end

        # MH step
        vLogPosteriorCandidate[iSim,1] = nObjCandidate;
        vLogPosteriorCandidate[iSim,2]= 1;
        # vLogPosteriorCandidate[iSim,3]= nMaxEv11;
        # vLogPosteriorCandidate[iSim,4]= nMaxEv12;

        Random.seed!(rnd_off+iSim*10+3)
        iRand = Random.rand(1);
        if iRand[1,1] <= nAlpha[1]
              mStructParam[iSim,:] = vParametersCandidate';
              nAccept              = nAccept + 1;
              nObj                 = nObjCandidate;
              vLogPosterior[iSim,1] = nObjCandidate;
              vLogPosterior[iSim,2]= 1;
              # vLogPosterior[iSim,3]= nMaxEv11;
              # vLogPosterior[iSim,4]= nMaxEv12;
              nCount               = nCount + 1;
        else
              mStructParam[iSim,:] = mStructParam[iSim-1,:];
              vLogPosterior[iSim]  = nObj;
              vLogPosterior[iSim,2]= 0;
              # vLogPosterior[iSim,3:4]= vLogPosterior[iSim-1,3:4];
        end

        nAcceptRate = nAccept/iSim;
        mStructParamCandidate[iSim,:] = vParametersCandidate';

        counter = counter + 1
        if counter == nOutput || nOutput == 0
           time_loop = signed(time_ns()-time_init_loop)/1000000000
           println("")
           println(" Draw number:  $iSim")
           println(" Remaining draws:  $(nSimRWMH-iSim)")
           println(" Post at Current Para:  $nObj")
           println(" Post at Proposed Para: $nObjCandidate")
           println(" Acceptance rate: $nAcceptRate")
           println(" Elapsed time: $time_loop")
           counter = 0
           time_init_loop = time_ns()

           savename     = nFileStrMH * "_loop_LogPosteriorCandidate.csv"
           CSV.write(savedir * savename,  DataFrame(LogPost = vLogPosteriorCandidate[iSim-nOutput+1:iSim]), header=false, append=true)

           savename     = nFileStrMH * "_loop_CandidateDraws.csv"
           CSV.write(savedir * savename,  DataFrame(mStructParamCandidate[iSim-nOutput+1:iSim,:], :auto), header=false, append=true)

           savename     = nFileStrMH * "_loop_Draws.csv"
           CSV.write(savedir * savename,  DataFrame(mStructParam[iSim-nOutput+1:iSim,:], :auto), header=false, append=true)

           savename     = nFileStrMH * "_loop_LogPosterior.csv"
           CSV.write(savedir * savename,  DataFrame(LogPost = vLogPosterior[iSim-nOutput+1:iSim]), header=false, append=true)

        end

    end

    #--------------------------------------------------------------------------
    # End of RWMH Loop
    #--------------------------------------------------------------------------

    mStructParam          = mStructParam[nBurnRWMH+1:end,:];
    mStructParamCandidate = mStructParamCandidate[nBurnRWMH+1:end,:];
    vLogPosterior         = vLogPosterior[nBurnRWMH+1:end,:];
    vLogPosteriorCandidate= vLogPosteriorCandidate[nBurnRWMH+1:end,:];

    # output

    savename     = nFileStrMH * "_Draws.csv"
    CSV.write(savedir * savename,  DataFrame(mStructParam, :auto), header=false)

    savename     = nFileStrMH * "_CandidateDraws.csv"
    CSV.write(savedir * savename,  DataFrame(mStructParamCandidate, :auto), header=false)

    savename     = nFileStrMH * "_LogPosterior.csv"
    CSV.write(savedir * savename,  DataFrame(vLogPosterior), header=false)

    savename     = nFileStrMH * "_LogPosteriorCandidate.csv"
    CSV.write(savedir * savename,  DataFrame(vLogPosteriorCandidate),header=false)

    time_total=signed(time_ns()-time_init_proc)/1000000000

    println("")
    println(" Total time: $time_total")

end
