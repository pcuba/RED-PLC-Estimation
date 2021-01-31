function myplot_priors2(pshape,p3,p4,p6,p7,param_labels)

% function plot_priors
% plots prior density
%
% INPUTS
%    o bayestopt_  [structure]ff
%    o M_          [structure]
%    o options_    [structure]
%
% OUTPUTS
%    None
%
% SPECIAL REQUIREMENTS
%    None

% Copyright (C) 2004-2012 Dynare Team
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



figurename = 'Priors';
npar = length(pshape);
[nbplt,nr,nc,lr,lc,nstar] = pltorg2(npar);

     
[x,f] = mydraw_prior_density(1,pshape,p3,p4,p6,p7);

x(f<0.001)=[]
f(f<0.001)=[]
hh = plot(x,f,'k--','linewidth',2);
set(hh,'color',[0.7 0.7 0.7]);
box on
axit tight
drawnow


end
    
