
# Function to get the HPD interval
function fHPDinterval(mDraws,nHPDprob)

    nParamDim = size(mDraws,2)
    nDrawDim = size(mDraws,1)

    mHPDband = zeros(2,nParamDim)
    nWidth = convert(Int64,floor((nHPDprob*nDrawDim)))

    for iParam = 1:nParamDim

        vDrawParSort = sort(mDraws[:,iParam], rev = true)

        nMinWidth = vDrawParSort[1] - vDrawParSort[nWidth]
        nBUP = 1
        nIter = 2
        while nIter <= nDrawDim - nWidth + 1
            nNewWidth = vDrawParSort[nIter] - vDrawParSort[nIter + nWidth - 1]
            if nNewWidth < nMinWidth
                nBUP = nIter
                nMinWidth = nNewWidth
            end

            nIter = nIter + 1

        end

        mHPDband[2,iParam] = vDrawParSort[nBUP]
        mHPDband[1,iParam] = vDrawParSort[nBUP + nWidth - 1]


    end

    return mHPDband
end
