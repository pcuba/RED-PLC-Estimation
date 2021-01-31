function fCondDistExactZLB_v5(mS_fore,mY_fore,vRegime,tt,pSigmaTE,pSigmaEta,pA,pA0,eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,mSigmaMEtilde,pKappa,mS_up,yy,nM, pInvSigmaMEtilde)


    nSmallME     =1E-8

    if pKappa < nSmallME
       print("ME is too small!!!")
    end

    # Measurement errors
    pSigmaME     = pKappa*mSigmaMEtilde;
    pA_Augmented = copy(pA)
    pA_NonAug    = copy(pA[:,1:8])

    nME               = size(pSigmaEta,1);
    pDetSigmaEta      = det(pSigmaEta)

    ny         = size(yy,1);
    ns         = size(mS_up,1);
    ne         = size(pSigmaEta,1);

    nZeta      = zeros(1,nM);
    iDlog      = zeros(1,nM);
    iLambda    = zeros(1,nM);
    mSimEta    = zeros(ne,nM);

    # Check ZLB and get log densities

    if yy[4] > pZLB # ZLB is not binding

        pOmegaUpBar1tilde = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi1_Eta)^(-1);
        pOmegaUpBar1      = pKappa*pOmegaUpBar1tilde;
        pOmegaBar1_2cond1 = Array(cholesky(Hermitian(pOmegaUpBar1[2:ne,2:ne]
                            - pOmegaUpBar1tilde[2:ne,1]*pOmegaUpBar1tilde[1,1]^(-1)*(pOmegaUpBar1[1,2:ne])' )));
        pInvSigmaPred_1   = (pKappa*mSigmaMEtilde + pA_NonAug*pPhi1_Eta*pPhi1_Eta'pA_NonAug')^(-1)
        pDetSigmaPred_1   = det(pKappa*mSigmaMEtilde + pA_NonAug*pPhi1_Eta*pPhi1_Eta'pA_NonAug')

        pEtaBar1_temp = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi1_Eta)^(-1)*pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde

        for iSim = 1:nM # iterate over particles

            iLambda[iSim] = 1 # set regime indicator
            # eta1Bar = eta1BarCoeff*[1.0;mS_up[1:8,iSim]]
            eta1Bar = eta1BarCoeff*vcat(ones(1),mS_up[1:8,iSim]);
            nZeta[1,iSim] = eta1Bar[1,1]
            # pEtaBar1 = pEtaBar1_temp*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )
            pEtaBar1 = pEtaBar1_temp*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )
            # compute log density
            # iTempExp11 =
            #     ( -(1/2)*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )'*
            #     pInvSigmaPred_1*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] ) )[1,1]

            iTempExp11 =
                ( -(1/2)*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )'*
                pInvSigmaPred_1*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) ) )[1,1]

            iTempExp11 = -(ny/2)*log(2*pi) - (1/2)*log(pDetSigmaPred_1) + iTempExp11

            # compute truncation constant
            iTempCDF1 = (nZeta[1,iSim] - pEtaBar1[1] )/((pOmegaUpBar1[1,1])^(1/2))
            iTempNormallogCDF1 = logcdf(Normal(0,1), iTempCDF1  )
            iTempNormalCDF1    = exp(iTempNormallogCDF1)

            # combine the terms
            iDlog[1,iSim] = iTempExp11 + iTempNormallogCDF1

            # Generate eta draws
            mSimEta[1,iSim]    = rand(truncated(Normal(pEtaBar1[1], (pOmegaUpBar1[1,1])^(1/2)), -Inf, nZeta[1,iSim] ) )
            pEtaBar1_2cond1    = pEtaBar1[2:ne] + pOmegaUpBar1[2:ne,1]*pOmegaUpBar1[1,1]^(-1)*(mSimEta[1,iSim] - pEtaBar1[1]);
            mSimEta[2:ne,iSim] = rand(MvNormal(pEtaBar1_2cond1, pOmegaBar1_2cond1));
            (mS_fore[:,iSim], mY_fore[:,iSim], vRegime[iSim]) = fStateSpaceSimulationaug(eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,pA,pA0, mS_up[:,iSim], mSimEta[:,iSim]);

        end

    else # ZLB is binding

        pInvSigmaPred_2exactZLB = (pKappa*mSigmaMEtilde[1:end-1,1:end-1] + pA_NonAug[1:end-1,:]*pPhi2_Eta*pPhi2_Eta'pA_NonAug[1:end-1,:]')^(-1)
        pDetSigmaPred_2exactZLB = det(pKappa*mSigmaMEtilde[1:end-1,1:end-1] + pA_NonAug[1:end-1,:]*pPhi2_Eta*pPhi2_Eta'pA_NonAug[1:end-1,:]')

        pOmegaUpBar2tilde = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi2_Eta)^(-1);
        pOmegaUpBar2      = pKappa*pOmegaUpBar2tilde;
        pOmegaBar2_2cond1 = Array(cholesky(Hermitian(pOmegaUpBar2[2:ne,2:ne]
                            - pOmegaUpBar2tilde[2:ne,1]*pOmegaUpBar2tilde[1,1]^(-1)*(pOmegaUpBar2[1,2:ne])' )));

        pEtaBar2_temp = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi2_Eta)^(-1)*pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde

        for iSim = 1:nM # iterate over particles

            iLambda[iSim] = 0
            # eta1Bar = eta1BarCoeff*[1.0;mS_up[1:8,iSim]]
            eta1Bar = eta1BarCoeff*vcat(ones(1),mS_up[1:8,iSim]);
            nZeta[1,iSim] = eta1Bar[1,1]

            # pEtaBar2 = pEtaBar2_temp*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*[pPhi2_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )
            pEtaBar2 = pEtaBar2_temp*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )
            # Compute log density of the binding with exact ZLB
            iTempExp21exactZLB = (-1/(2)*(yy[1:end-1] - pA0[1:end-1] - pA_NonAug[1:end-1,:]*pPhi2_0 - pA_Augmented[1:end-1,:]*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5,iSim] ) )'*
                pInvSigmaPred_2exactZLB*(yy[1:end-1] - pA0[1:end-1] - pA_NonAug[1:end-1,:]*pPhi2_0 - pA_Augmented[1:end-1,:]*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5,iSim] ) ) )[1,1]
                
            iTempExp21exactZLB = -((ny-1)/2)*log(2*pi) - (1/2)*log(pDetSigmaPred_2exactZLB) + iTempExp21exactZLB

            # Compute weights on density components
            iTempCDF2 = (nZeta[1,iSim] - pEtaBar2[1] )/ ((pOmegaUpBar2[1,1])^(1/2))
            iTempNormallogCDF2 = logcdf(Normal(0,1), -iTempCDF2  )

            # combine terms
            iDlog[1,iSim] = iTempExp21exactZLB + iTempNormallogCDF2

            # Generate eta draws
            mSimEta[1,iSim]    = rand(truncated(Normal(pEtaBar2[1], (pOmegaUpBar2[1,1])^(1/2)), nZeta[1,iSim] , Inf ))
            pEtaBar2_2cond1    = pEtaBar2[2:ne] + pOmegaUpBar2[2:ne,1]*pOmegaUpBar2[1,1]^(-1)*(mSimEta[1,iSim] - pEtaBar2[1]);
            mSimEta[2:ne,iSim] = rand(MvNormal(pEtaBar2_2cond1, pOmegaBar2_2cond1));

            (mS_fore[:,iSim], mY_fore[:,iSim], vRegime[iSim]) = fStateSpaceSimulationaug(eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,pA,pA0, mS_up[:,iSim], mSimEta[:,iSim]);

        end # particle iteration

    end # ZLB binding vs. nonbinding

    return mSimEta, iDlog, mS_fore, mY_fore, vRegime;
end # End of function
