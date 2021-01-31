function md = mhm_marginal_density(thetaHistory,fvalHistory,truncation_point)
                              
fvalHistory=fvalHistory;                        

chain_length=size(thetaHistory,2);
nparams = size(thetaHistory,1);

critval = chi2inv(truncation_point,nparams);

mu = mean(thetaHistory,2);

varcov = cov(thetaHistory');
detVarcov = det(varcov);

invVarcov = inv(varcov);


md = 0;
for i = 1:chain_length

this_theta = thetaHistory(:,i);

quadratic_form = (this_theta-mu)'*invVarcov*(this_theta-mu);

if quadratic_form < critval
    
    numerator = -log(truncation_point)...
        +(nparams/2)*log(2*pi)...
        -0.5*log(detVarcov)...
        -0.5*quadratic_form ;
    md =md+exp(numerator+fvalHistory(i)-fvalHistory(1));

end


end

md=-log(md/chain_length)-fvalHistory(1);


