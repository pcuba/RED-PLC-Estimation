
function fSystematicResampling(w)

    np = length(w);

    w = w';
    cw = cumsum(w);
    uu = zeros(np,1);
    csi=rand(1)[1];

    for j=1:np

        uu[j] = j-1 + csi;

    end

        indx = zeros(np, 1);

    for i = 1:np

        u = uu[i] / np;

        j=1;
        while j <= np
            if (u < cw[j])
                break
            end

           j = j + 1;

        end

        indx[i] = j;

    end


return indx
end
