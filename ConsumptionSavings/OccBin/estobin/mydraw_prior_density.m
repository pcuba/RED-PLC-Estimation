function [x,f,abscissa,dens,binf,bsup] = mydraw_prior_density(indx,pshape,p3,p4,p6,p7)
% Computes values of prior densities at many points (before plotting)
%
% INPUTS
%    indx          [integer]    Parameter number.
%    bayestopt_    [structure]  Describes the prior beliefs.
%    
% OUTPUTS
%    x             [double]     Row vector, subset of 'abscissa' such as the density is less than 10
%    f             [double]     Row vector, subset of 'dens' such as the density is less than 10
%    abscissa      [double]     Row vector, abscissa 
%    dens          [double]     Row vector, density
%    binf:         [double]     Scalar, first element of x
%    bsup:         [double]     Scalar, last element of x


% Copyright (C) 2004-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


truncprior = 1e-3;
steps = 200;

switch pshape(indx)
  case 1% Beta prior
    density = @(x,a,b,aa,bb) betapdf((x-aa)/(bb-aa), a, b)/(bb-aa);
    infbound = betainv(truncprior,p6(indx),p7(indx))*(p4(indx)-p3(indx))+p3(indx);
    supbound = betainv(1-truncprior,p6(indx),p7(indx))*(p4(indx)-p3(indx))+p3(indx);
    abscissa = linspace(infbound,supbound,steps);
    dens = density(abscissa,p6(indx),p7(indx),p3(indx),p4(indx));
  case 2% Generalized Gamma prior
    density = @(x,a,b,c) gampdf(x-c,a,b);
    try
        infbound = gaminv(truncprior,p6(indx),p7(indx))+p3(indx);
        supbound = gaminv(1-truncprior,p6(indx),p7(indx))+p3(indx);
    catch
        % Workaround for ticket #161
        if exist('OCTAVE_VERSION')
            error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt_.name{indx} ', or use another shape'])
        else
            rethrow(lasterror)
        end
    end
    abscissa = linspace(infbound,supbound,steps);
    dens = density(abscissa,p6(indx),p7(indx),p3(indx));
  case 3% Gaussian prior
    infbound = norminv(truncprior,p6(indx),p7(indx)); 
    supbound = norminv(1-truncprior,p6(indx),p7(indx));
    abscissa = linspace(infbound,supbound,steps);
    dens = normpdf(abscissa,p6(indx),p7(indx));  
  case 4% Inverse-gamma of type 1 prior
    
        infbound = 1/sqrt(gaminv(1-10*truncprior, p7(indx)/2, 2/p6(indx)))+p3(indx);
        supbound = 1/sqrt(gaminv(10*truncprior, p7(indx)/2, 2/p6(indx)))+p3(indx);
    abscissa = linspace(infbound,supbound,steps);
    dens = exp(lpdfig1(abscissa-p3(indx),p6(indx),p7(indx)));  
  case 5% Uniform prior
    infbound = p6(indx);
    supbound = p7(indx);
    abscissa = linspace(infbound,supbound,steps);
    dens = ones(1, steps) / (supbound-infbound);
  case 6% Inverse-gamma of type 2 prior
    
        infbound = 1/(gaminv(1-10*truncprior, p7(indx)/2, 2/p6(indx)))+p3(indx);
        supbound = 1/(gaminv(10*truncprior, p7(indx)/2, 2/p6(indx)))+p3(indx);
    
    abscissa = linspace(infbound,supbound,steps);
    dens = exp(lpdfig2(abscissa-p3(indx),p6(indx),p7(indx)));
  otherwise
    error(sprintf('draw_prior_density: unknown distribution shape (index %d, type %d)', indx, pshape(indx)));
end 

if pshape(indx) ~= 5 
    [junk,k1] = max(dens);
    if k1 == 1 || k1 == length(dens)
        k = find(dens > 10);
        dens(k) = NaN;
    end
end
binf = abscissa(1);
bsup = abscissa(end);
x = abscissa;
f = dens;

