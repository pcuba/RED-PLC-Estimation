function [kappa theta] = gamma_specification(mean,std);

theta = std^2/mean;
kappa = mean/theta;
