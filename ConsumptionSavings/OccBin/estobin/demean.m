function c = demean(x)
% Removes the mean of each column of a matrix.
 
%@info:
%! @deftypefn {Function File} {@var{c} =} demean (@var{x})
%! @anchor{demean}
%! This function removes the mean of each column of a matrix.
%!
%! @strong{Inputs}
%! @table @var
%! @item x
%! Matlab matrix (T-by-N).
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item c
%! Matlab matrix (T-by-N). The demeaned x matrix.
%! @end table
%! 
%! @strong{This function is called by:} 
%! @ref{compute_cova}, @ref{compute_acov}, @ref{compute_std}.
%! 
%! @strong{This function calls:} 
%! @ref{ndim},
%!    
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
% stephane DOT adjemian AT ens DOT fr
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
   
if ndim(x)==1
    c = x-nanmean(x);
elseif ndim(x)==2
    c = bsxfun(@minus,x,nanmean(x));
else
    error('descriptive_statistics::demean:: This function is not implemented for arrays with dimension greater than two!')
end