function im = MLEM_im(im, Summ, epps, varargin)
%MLEM_im Computes the Maximum Likelihood Expectation Maximization (MLEM)
%estimates
%
% Example:
%   im = MLEM_im(im, A, epps, uu, Summ)
% INPUTS:
%   im = The current estimate
%   Summ = sum(A,2) (normalization)
%   epps = Small constant to prevent division by zero
%   A = The (transpose of the) (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   is_transposed = true if A matrix is the transpose of it, false if not
%
% OUTPUTS:
%   im = The updated estimate
% 
% See also OSEM_im, ROSEM_subiter, BSREM_subiter, DRAMA_subiter,
% RBI_subiter, COSEM_im, ACOSEM_im, ECOSEM_im, MBSREM

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
if nargin == 6
    if varargin{3}
        im = (im./Summ).*(varargin{1}*(varargin{2}./(varargin{1}'*im + epps)));
    else
        im = (im./Summ).*(varargin{1}'*(varargin{2}./(varargin{1}*im + epps)));
    end
elseif nargin == 4
    im = (im./Summ).*(varargin{1} + epps);
else
    error('Invalid number of input arguments')
end