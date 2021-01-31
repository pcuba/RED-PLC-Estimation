function fCondDist_v5(mS_fore,mY_fore,vRegime,tt,pSigmaTE,pSigmaEta,pA,pA0,eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,pSigmaMEtilde,pKappa,mS_up,yy,nM,pInvSigmaMEtilde)

nSmallME     = 10^-8
# Measurement errors

pSigmaME     = pKappa*pSigmaMEtilde;

pA_Augmented = copy(pA)
pA_NonAug = copy(pA[:,1:8])

nME               = size(pSigmaEta,1);
ne                = size(pSigmaTE,1);

pOmegaUpBar1tilde = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi1_Eta)^(-1);
pOmegaUpBar2tilde = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi2_Eta)^(-1);
pOmegaUpBar1      = pKappa*pOmegaUpBar1tilde;
pOmegaUpBar2      = pKappa*pOmegaUpBar2tilde;

pOmegaBar1_2cond1 = Array(cholesky(Hermitian(pOmegaUpBar1[2:ne,2:ne]
                    - pOmegaUpBar1tilde[2:ne,1]*pOmegaUpBar1tilde[1,1]^(-1)*(pOmegaUpBar1[1,2:ne])' )));
pOmegaBar2_2cond1 = Array(cholesky(Hermitian(pOmegaUpBar2[2:ne,2:ne]
                    - pOmegaUpBar2tilde[2:ne,1]*pOmegaUpBar2tilde[1,1]^(-1)*(pOmegaUpBar2[1,2:ne])' )));

pDetSigmaEta      = det(pSigmaEta)
pDetOmegaUpBar1   = det(pOmegaUpBar1)
pDetOmegaUpBar2   = det(pOmegaUpBar2)


pInvSigmaPred_1   = (pKappa*pSigmaMEtilde + pA_NonAug*pPhi1_Eta*pPhi1_Eta'pA_NonAug')^(-1)
pInvSigmaPred_2   = (pKappa*pSigmaMEtilde + pA_NonAug*pPhi2_Eta*pPhi2_Eta'pA_NonAug')^(-1)
pDetSigmaPred_1   = det(pKappa*pSigmaMEtilde + pA_NonAug*pPhi1_Eta*pPhi1_Eta'pA_NonAug')
pDetSigmaPred_2   = det(pKappa*pSigmaMEtilde + pA_NonAug*pPhi2_Eta*pPhi2_Eta'pA_NonAug')

ny         = size(yy,1);
ns         = size(mS_up,1);
ne         = size(pSigmaEta,1);
mS_fore    = zeros(ns,nM);

nZeta      = zeros(1,nM);
iDlog      = zeros(1,nM);
iLambda    = zeros(1,nM);
mSimEta    = zeros(ne,nM);

pEtaBar1_temp = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi1_Eta)^(-1)*pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde;
pEtaBar2_temp = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde*pA_NonAug*pPhi2_Eta)^(-1)*pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde;
for iSim = 1:nM

    # eta1Bar = eta1BarCoeff*[1.0;mS_up[1:8,iSim]]
    eta1Bar = eta1BarCoeff*vcat(ones(1),mS_up[1:8,iSim]);
    nZeta[1,iSim] = eta1Bar[1,1]

    # pEtaBar1 = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi1_Eta'*pA_NonAug'*(pSigmaMEtilde)^(-1)*pA_NonAug*pPhi1_Eta)^(-1)*pPhi1_Eta'*pA_NonAug'*pInvSigmaMEtilde*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )
    # pEtaBar2 = (pKappa*Matrix{Float64}(I, nME, nME) + pPhi2_Eta'*pA_NonAug'*(pSigmaMEtilde)^(-1)*pA_NonAug*pPhi2_Eta)^(-1)*pPhi2_Eta'*pA_NonAug'*pInvSigmaMEtilde*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*[pPhi2_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )


    # pEtaBar1 = pEtaBar1_temp*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )
    # pEtaBar2 = pEtaBar2_temp*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*[pPhi2_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )

    pEtaBar1 = pEtaBar1_temp*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )
    pEtaBar2 = pEtaBar2_temp*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )

    # Compute weights on density components

    if pKappa < nSmallME

        if nZeta[1,iSim] - pEtaBar1[1] >= 0
            iTempNormalCDF1 = 1
        else
            iTempNormalCDF1 = 0
        end

        if nZeta[1,iSim] - pEtaBar2[1] >= 0
            iTempNormalCDF2 = 0
        else
            iTempNormalCDF2 = 1
        end

    else
        iTempCDF1 = (nZeta[1,iSim] - pEtaBar1[1] )/((pOmegaUpBar1[1,1])^(1/2))
        #iTempNormalCDF1 = cdf(Normal(0,1), iTempCDF1  )
        iTempNormallogCDF1 = logcdf(Normal(0,1), iTempCDF1  )
        iTempNormalCDF1    = exp(iTempNormallogCDF1)

        iTempCDF2 = (nZeta[1,iSim] - pEtaBar2[1] )/ ((pOmegaUpBar2[1,1])^(1/2))
        #iTempNormalCDF2 = 1 - cdf(Normal(0,1), iTempCDF2  )
        iTempNormallogCDF2 = logcdf(Normal(0,1), -iTempCDF2  )
        iTempNormalCDF2    = exp(iTempNormallogCDF2)
    end

    # compute log densities

    # iTempExp11 =
    #     ( -(1/2)*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )'*
    #     pInvSigmaPred_1*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*[pPhi1_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] ) )[1,1]

    iTempExp11 =
        ( -(1/2)*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )'*
        pInvSigmaPred_1*(yy - pA0 - pA_NonAug*pPhi1_0 - pA_Augmented*vcat(pPhi1_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) ) )[1,1]

    iTempExp11 = -(ny/2)*log(2*pi) - (1/2)*log(pDetSigmaPred_1) + iTempExp11

    # iTempExp21 = (-1/(2)*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*[pPhi2_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] )'*
    #     pInvSigmaPred_2*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*[pPhi2_1*mS_up[1:8,iSim] ; mS_up[5,iSim] ] ) )[1,1]

    iTempExp21 = (-1/(2)*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) )'*
        pInvSigmaPred_2*(yy - pA0 - pA_NonAug*pPhi2_0 - pA_Augmented*vcat(pPhi2_1*mS_up[1:8,iSim] , mS_up[5:5,iSim] ) ) )[1,1]

    iTempExp21 = -(ny/2)*log(2*pi) - (1/2)*log(pDetSigmaPred_2) + iTempExp21

    # compute log weighted density differential and log particle weight

    if iTempNormalCDF2 == 0 && iTempNormalCDF1 > 0
        iDlog[1,iSim] = iTempExp11 + iTempNormallogCDF1
        iLambda[iSim] = 1
        #if tt == 53 && iSim == 1
        #    println(1)
        #    println(iDlog[1,1])
        #end
    elseif iTempNormalCDF1 == 0 && iTempNormalCDF2 > 0
        iDlog[1,iSim] = iTempExp21 + iTempNormallogCDF2
        iLambda[iSim] = 0
        #if tt == 53 && iSim == 1
        #    println(2)
        #    println(iDlog[1,1])
        #end
    elseif iTempNormalCDF1 == 0 && iTempNormalCDF2 == 0

        iDifLogD2D1 =  0
        iDlog[1,iSim] = NaN
        iLambda[iSim] = 0

        #if tt == 53 && iSim == 1
        #    println(3)
        #    println(iTempNormalCDF1)
        #    println(iTempNormalCDF2)
        #    println(iDifLogD2D1)
        #end

    elseif ( iTempExp11 + log( iTempNormalCDF1 ) ) >= ( iTempExp21 + log( iTempNormalCDF2 ) )
        iDifLogD2D1   =  ( iTempExp21 + iTempNormallogCDF2 -
                         iTempExp11 - iTempNormallogCDF1 )
        iDlog[1,iSim] =  ( iTempExp11 + iTempNormallogCDF1 +
                         log( 1 + exp(iDifLogD2D1) ) )
        iLambda[iSim] = 1 / (1 + exp( iDifLogD2D1 ) )
        #if tt == 53 && iSim == 1
        #    println(4)
        #    println(iTempNormalCDF1)
        #    println(iTempNormalCDF2)
        #    println(iDifLogD2D1)
        #end
    else
        iDifLogD1D2 =  ( iTempExp11 + iTempNormallogCDF1 -
                       iTempExp21 -  iTempNormallogCDF2 )
        iDlog[1,iSim] =  ( iTempExp21 + iTempNormallogCDF2 +
                         log( 1 + exp(iDifLogD1D2) ) )
        iLambda[iSim] = 1 / (1 + exp( -iDifLogD1D2 ) )
        #if tt == 53 && iSim == 1
        #    println(5)
        #    println(iDlog[1,1])
        #end
    end


    # Generate eta draws

    iRand = rand(1);
    if iRand[1,1] <= iLambda[1,iSim]

         mSimEta[1,iSim]    = rand(truncated(Normal(pEtaBar1[1], (pOmegaUpBar1[1,1])^(1/2)), -Inf, nZeta[1,iSim] ) )
         pEtaBar1_2cond1    = pEtaBar1[2:ne] + pOmegaUpBar1[2:ne,1]*pOmegaUpBar1[1,1]^(-1)*(mSimEta[1,iSim] - pEtaBar1[1]);
         mSimEta[2:ne,iSim] = rand(MvNormal(pEtaBar1_2cond1, pOmegaBar1_2cond1));

    else

         mSimEta[1,iSim]    = rand(truncated(Normal(pEtaBar2[1], (pOmegaUpBar2[1,1])^(1/2)), nZeta[1,iSim] , Inf ))
         pEtaBar2_2cond1    = pEtaBar2[2:ne] + pOmegaUpBar2[2:ne,1]*pOmegaUpBar2[1,1]^(-1)*(mSimEta[1,iSim] - pEtaBar2[1]);
         mSimEta[2:ne,iSim] = rand(MvNormal(pEtaBar2_2cond1, pOmegaBar2_2cond1));

    end

    (mS_fore[:,iSim], mY_fore[:,iSim], vRegime[iSim]) = fStateSpaceSimulationaug(eta1BarCoeff,pPhi1_0,pPhi1_1,pPhi1_Eta,pPhi2_0,pPhi2_1,pPhi2_Eta,pA,pA0, mS_up[:,iSim], mSimEta[:,iSim]);

end
return mSimEta, iDlog, mS_fore, mY_fore, vRegime;
end # End of function
