function im = OSEM_im(im, varargin)
%OSEM_IM Computes the Ordered Subsets Expectation Maximization (OSEM)
%estimates
%
% Examples:
%   im = OSEM_im(im, A, epps, uu, Summ)
%   im = OSEM_im(im, RHS, Summ, epps)
% INPUTS:
%   im = The current estimate
%   A = The (transpose of the) (sparse) system matrix at current subset
%   epps = Small constant to prevent division by zero
%   uu = Measurements at current subset
%   Summ = sum(A,2)
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero. 
%   is_transposed = true if A matrix is the transpose of it, false if not
%   RHS = The right hand side of OSEM (RHS = A'*(uu./(A*im + SinDelayed))) 
%
% OUTPUTS:
%   im = The updated estimate
% 
% See also ROSEM_subiter, BSREM_subiter, DRAMA_subiter, RBI_subiter,
% COSEM_im, ACOSEM_im, ECOSEM_im, MBSREM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 7
    if varargin{6}
        im = (im./varargin{4}).*(varargin{1}*(varargin{3}./(varargin{1}'*im + varargin{2} + varargin{5})) + varargin{2});
    else
        im = (im./varargin{4}).*(varargin{1}'*(varargin{3}./(varargin{1}*im + varargin{2} + varargin{5})) + varargin{2});
    end
elseif nargin == 4
    im = im./varargin{2} .* (varargin{1} + varargin{3});
else
    error('Invalid number of input arguments')
end