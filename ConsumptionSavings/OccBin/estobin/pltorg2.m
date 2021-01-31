function [nbplt,nr,nc,lr,lc,nstar] = pltorg(number)

% Copyright (C) 2004-2008 Dynare Team
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

nrstar = 5;
ncstar = 4;
nstar  = nrstar*ncstar;
nbplt  = 0;
nr     = 0;
nc     = 0;
lr     = 0;
lc     = 0;
if number == 1
    nbplt = 1;
    nr    = 1;
    nc    = 1;
elseif number == 2
    nbplt = 1;
    nr    = 2;
    nc    = 1;
elseif number == 3
    nbplt = 1;
    nr    = 3;
    nc    = 1;
elseif number == 4
    nbplt = 1;
    nr    = 2;
    nc    = 2;
elseif number == 5
    nbplt = 1;
    nr    = 3;
    nc    = 2;
elseif number == 6
    nbplt = 1;
    nr    = 3;
    nc    = 2;    
elseif number == 7
    nbplt = 1;
    nr    = 3;
    nc    = 3;    
elseif number == 8
    nbplt = 1;
    nr    = 3;
    nc    = 3;
elseif number == 9
    nbplt = 1;
    nr    = 3;
    nc    = 3;
elseif number == 10
    nbplt = 1;
    nr    = 4;
    nc    = 3;
elseif number == 11
    nbplt = 1;
    nr    = 4;
    nc    = 3;
elseif number == 12
    nbplt = 1;
    nr    = 4;
    nc    = 3;
else
    if number/nstar == round(number/nstar)
        nbplt = number/nstar;
        nr    = nrstar;
        nc    = ncstar;
        lr    = nr;
        lc    = nc; 
    else
        nbplt = ceil(number/nstar);
        nr    = nrstar;
        nc    = ncstar;
        reste = number-(nbplt-1)*nstar;
        if reste == 1
            lr    = 1;
            lc    = 1;
        elseif reste == 2
            lr    = 2;
            lc    = 1;
        elseif reste == 3
            lr    = 3;
            lc    = 1;
        elseif reste == 4
            lr    = 2;
            lc    = 2;
        elseif reste == 5
            lr    = 3;
            lc    = 2;
        elseif reste == 6
            lr    = 3;
            lc    = 2;    
        elseif reste == 7
            lr    = 3;
            lc    = 3;    
        elseif reste == 8
            lr    = 3;
            lc    = 3;
        end
    end
end