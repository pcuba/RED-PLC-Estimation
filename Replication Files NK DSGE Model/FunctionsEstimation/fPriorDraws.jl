
# Function to draw from any prior specification vector of parameters Theta
function fPriorDraws(nNumDraw,mPriorSpec)
    mSpec = mPriorSpec
    nNumTheta = size(mSpec,1)

    mPriorDraws = zeros(nNumDraw,nNumTheta)

    vShape = mSpec[:,1]
    vArg1  = mSpec[:,2]
    vArg2  = mSpec[:,3]
    vMask  = mSpec[:,4]
    vFix   = mSpec[:,5]
    vMaskInv = ones(size(vMask)) - vMask
    vShape = vShape.*vMaskInv

    for iTheta = 1 : nNumTheta

        if vMask[iTheta] == 1
            mPriorDraws[:,iTheta] .= vFix[iTheta]
        else

            if  mSpec[iTheta, 1] == 1.0 # BETA
                a = (1-vArg1[iTheta])*vArg1[iTheta]^2/vArg2[iTheta]^2 - vArg1[iTheta]
                b = a*(1/vArg1[iTheta] - 1)
                mPriorDraws[:,iTheta] = rand(Beta(a,b),nNumDraw)
            elseif mSpec[iTheta, 1] == 2.0 # GAMMA
                b = vArg2[iTheta]^2/vArg1[iTheta]
                a = vArg1[iTheta]/b
                mPriorDraws[:,iTheta] = rand(Gamma(a,b),nNumDraw)
            elseif mSpec[iTheta, 1] == 3 # NORMAL (Second argument is StdD)
                mPriorDraws[:,iTheta] = rand(Normal(vArg1[iTheta],vArg2[iTheta]),nNumDraw)
            elseif mSpec[iTheta, 1] == 4 # INVGAMMA
                a = vArg2[iTheta]/2
                b = vArg1[iTheta]^2*(a)
                mPriorDraws[:,iTheta] = sqrt.(rand(InverseGamma(a,b),nNumDraw))
            else mSpec[iTheta, 1] == 5 # UNIFORM
                mPriorDraws[:,iTheta] = rand(Uniform(vArg1[iTheta],vArg2[iTheta]),nNumDraw)
            end
        end
    end

    return mPriorDraws
end
