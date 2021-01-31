function run_LikelihoodEval(nSimEval, ME_scale, vFilter, vM, nDataSet,
                mData, nStartDate, nEndDate, mDraw, nParamDraw,
                nSistResFilter,
                pZLB, rnd_off,
                Smolyak_fixed_grid,options_int_use,options_float_use,sPLC)


    time_init_proc=time_ns()
    nPath = pwd()
    if ME_scale < 0.01
        nMEpercent = 0
    else
        nMEpercent = convert(Int,ME_scale*100)
    end
    nMEpercent = string(nMEpercent)
    if length(nMEpercent) == 1
        nMEpercent = "0"*nMEpercent
    end
    savedir = "$(pwd())/LikelihoodEvaluations/" * "ME" * "$(nMEpercent)" * "_" * nDataSet * "/"

    # Create directory where everything will be saved, if it exists skip creation
    try mkdir(savedir) catch; end


    #--------------------------------------------------------------------------
    # Load prior specification
    #--------------------------------------------------------------------------
    # include(nPath*"/PriorSpec/" * nNamePriorSpec * ".jl")

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

    pInvSigmaMEtilde = mSigmaMEtilde^(-1);

    mSigmaME     = ME_scale*mSigmaMEtilde;

    #--------------------------------------------------------------------------
    # Loop filters and draws
    #--------------------------------------------------------------------------

    nData = size(mY,1)
    iFilterNum = 0

    for iFilter in vFilter

        #--------------------------------------------------------------------------
        # Create output matrices
        #--------------------------------------------------------------------------
        mTimePrior    = Array{Float64,2}(undef,nSimEval+1,nParamDraw)*0
        mLogPosterior = Array{Float64,2}(undef,nSimEval,nParamDraw)*0
        mLogStats     = Array{Float64,2}(undef,2,nParamDraw)*0

        iFilterNum = iFilterNum + 1
        nM = vM[iFilterNum]

        for iDraw = 1:nParamDraw

            time_init_loop=time_ns()
            mLogLikeHInc  = Array{Float64,2}(undef,nSimEval,nData+1)*0

            #--------------------------------------------------------------------------
            # Load parameters
            #--------------------------------------------------------------------------

            vParameters = mDraw[iDraw,:]

            nParametersValid = ( (vParameters[2]<=1)*(vParameters[2]>=0)*(vParameters[5]<=1)*(vParameters[6]<=1)*
                    (vParameters[7]<=1)*(vParameters[8]<=1)*(vParameters[9]>0)*
                    (vParameters[10]>0)*(vParameters[11]>0)*(vParameters[12]>0) )

            #--------------------------------------------------------------------------
            # Model Solution / State-space Representation / Likelihood
            #--------------------------------------------------------------------------

            if  nParametersValid == 1

                # Solve model
                # Random.seed!(rnd_off+iSim*10+1)
                Random.seed!(rnd_off+1)
                sParam = fCreateParamDict(vParameters)

                flag_out, sParam_StateSpace = fStateSpaceForm(sParam,options_int_use,options_float_use,Smolyak_fixed_grid,sPLC)

                if flag_out < 1
                    nParametersValid = 0
                end
            end

            if nParametersValid == 1

                # evaluate likelihood function multiple times

                for iSim = 1:nSimEval
                    iTimeEval=time_ns()
                    Random.seed!(rnd_off+iSim*10+2)
                    if iFilter == "BSPF"
                       (lik, s_up_stat, s_lag_stat, Eta_stat, vRegime, Neff) = fBSPF(mY, sParam_StateSpace, mSigmaME, ME_scale, nM);
                    elseif iFilter == "COPFexactZLB"
                       (lik, s_up_stat, s_lag_stat, Eta_stat, vRegime, Neff) = fCOPFaugExactZLB(mY, sParam_StateSpace, mSigmaMEtilde, ME_scale, nM, pInvSigmaMEtilde);
                    elseif iFilter == "COPF"
                       (lik, s_up_stat, s_lag_stat, Eta_stat, vRegime, Neff) = fCOPF(mY, sParam_StateSpace, mSigmaMEtilde, ME_scale, nM, pInvSigmaMEtilde);
                    elseif iFilter == "KF"
                       (lik, s_up_stat) = fKF(mY, sParam_StateSpace, mSigmaME);
                    end
                    # use for conditional likelihood lik = lik[2:end,:]
                    mLogPosterior[iSim,iDraw] = sum(lik)
                    mLogLikeHInc[iSim,1:end-1] = lik'
                    mLogLikeHInc[iSim,end] = sum(lik)
                    mTimePrior[iSim,iDraw] = signed(time_ns()-iTimeEval)/1000000000
                end

                nMean = mean(mLogPosterior[:,iDraw]) + log(mean(exp.(mLogPosterior[:,iDraw] .- mean(mLogPosterior[:,iDraw]))))
                if nMean == Inf || nMean == -Inf
                    mTimePrior[end,iDraw] = NaN
                    mLogStats[1,iDraw] = NaN
                    mLogStats[2,iDraw] = NaN
                else
                    mTimePrior[end,iDraw] = mean(mTimePrior[2:end-1,iDraw])
                    mLogStats[1,iDraw] = std(mLogPosterior[:,iDraw], corrected = false)
                    mLogStats[2,iDraw] = mean(mLogPosterior[:,iDraw]) +
                    log(mean(exp.(mLogPosterior[:,iDraw] .- mean(mLogPosterior[:,iDraw]))))
                end

                time_loop = signed(time_ns()-time_init_loop)/1000000000
                println("")
                println("Parameter number:  $iDraw")
                println("Filter: "*iFilter )
                println("Elapsed time: $time_loop")

            end # if parameter is valid

            savename     = nDataSet * "_" *iFilter * "$nM" * "_IncrLogLH_Draw" * "$iDraw" * ".csv"
            CSV.write(savedir * savename,  DataFrame(mLogLikeHInc), header=false)

        end # End parameter loop

        savename     = iFilter * "$nM" * "_ParamDraws"  * ".csv"
        CSV.write(savedir * savename,  DataFrame(mDraw[1:nParamDraw,:], :auto), header=false)

        savename     = iFilter * "$nM" * "_LogLH"  * ".csv"
        CSV.write(savedir * savename,  DataFrame(mLogPosterior), header=false)

        savename     = iFilter * "$nM" * "_Time"  * ".csv"
        CSV.write(savedir * savename,  DataFrame(mTimePrior), header=false)

        savename     = iFilter * "$nM" * "_LogStats"  * ".csv"
        CSV.write(savedir * savename,  DataFrame(mLogStats, :auto), header=false)

    end # filter loop

end # function
