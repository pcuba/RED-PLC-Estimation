function [input6 input7] = get_dist_inputs(codes,means,stds)

nparams= length(codes);
input6 = zeros(nparams,1);
input7 = input6;


for i = 1:nparams
    if codes(i) == 1 % BETA
        mu = means(i);
        sigma2 = stds(i)^2;
        xsol = fsolve(@(x) [x(1)/(x(1)+x(2))-mu; x(1)*x(2)/((x(1)+x(2))^2*(x(1)+x(2)+1))-sigma2],[0.5 0.1],optimset('Display','off'));
        input6(i) = xsol(1);
        input7(i) = xsol(2);
    elseif codes(i) == 2 % GAMMA
        [input6(i) input7(i)]=gamma_specification(means(i),stds(i));
    elseif codes(i) == 3 % NORMAL
        input6(i) = means(i);
        input7(i) = stds(i);
    elseif codes(i) == 4 % INVERSE GAMMA
        [input6(i) input7(i)]=inverse_gamma_specification(means(i), stds(i),1,0);
    elseif codes(i) == 5 % UNIFORM
        input6(i) = nan; % dist is only based on bounds
        input7(i) = nan; % dist is only based on bounds 
    end
end
