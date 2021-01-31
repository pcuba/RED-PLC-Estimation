
# Function to get log density of priors
function fPriorDensity(vTheta,mPriorSpec)
    mSpec     = mPriorSpec
    nNumTheta = size(mSpec,1)
    vDensity  = zeros(nNumTheta,1)

    vShape = mSpec[:,1]
    vArg1  = mSpec[:,2]
    vArg2  = mSpec[:,3]
    vMask  = mSpec[:,4]
    vFix   = mSpec[:,5]
    vMaskInv = ones(size(vMask)) - vMask
    vShape = vShape.*vMaskInv

    for iTheta = 1 : nNumTheta

        if vMask[iTheta] == 1
            vDensity[iTheta,1] = 1
        else

            if     mSpec[iTheta, 1] == 1.0 # BETA
                a = (1-vArg1[iTheta])*vArg1[iTheta]^2/vArg2[iTheta]^2 - vArg1[iTheta]
                b = a*(1/vArg1[iTheta] - 1)
                vDensity[iTheta,1] = pdf(Beta(a,b),vTheta[iTheta,1])
            elseif mSpec[iTheta, 1] == 2.0 # GAMMA
                b = vArg2[iTheta]^2/vArg1[iTheta]
                a = vArg1[iTheta]/b
                vDensity[iTheta,1] = pdf(Gamma(a,b),vTheta[iTheta,1])
            elseif mSpec[iTheta, 1] == 3 # NORMAL (Second argument is StdD)
                vDensity[iTheta,1] = pdf(Normal(vArg1[iTheta],vArg2[iTheta]),vTheta[iTheta,1])
            elseif mSpec[iTheta, 1] == 4 # INVGAMMA
                a = vArg2[iTheta]/2
                b = vArg1[iTheta]^2*(a)
                vDensity[iTheta,1] = 2*vTheta[iTheta,1]*pdf(InverseGamma(a,b),vTheta[iTheta,1]^2)
            else   mSpec[iTheta, 1] == 5 # UNIFORM
                vDensity[iTheta,1] = pdf(Uniform(vArg1[iTheta],vArg2[iTheta]),vTheta[iTheta,1])
            end
        end
    end

    nDensity = 1.0
    for iTheta = 1 : nNumTheta
        nDensity = nDensity * vDensity[iTheta,1]
    end

    nLogDens = log(nDensity)

    return nLogDens
end
