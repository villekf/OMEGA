function im = RBI_subiter(im, A, uu, epps, Summ, D, SinDelayed, is_transposed, varargin)
%RBI_SUBITER Computes the Rescaled Block Iterative (RBI) estimate for the
%current subset/sub-iteration
%
% Examples:
%   im = RBI_subiter(im, A, uu, epps, Summ, D, SinDelayed)
%   im = RBI_subiter(im, A, uu, epps, Summ, D, SinDelayed, beta, dU)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   epps = Small constant to prevent division by zero
%   Summ = sum(A,2)
%   D = Sum of the complete data system matrix divided by the number of
%   subsets (D = sum(B,2)/subsets, where B = [A_1,A_2,...,A_subsets]
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero. 
%   is_transposed = true if A matrix is the transpose of it, false if not
%   beta = Regularization parameter (0, if no regularization)
%   dU = Gradient of the prior
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_iter, BSREM_subiter, MBSREM, COSEM_OSL, ROSEM_subiter,
% OSL_OSEM

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

%%%% This version is from: Block-iterative techniques for fast 4D
%%%% reconstruction using a priori motion models in gated cardiac SPECT,
%%%% David S Lalush and Benjamin M W Tsui, 1998 Phys. Med. Biol. 43 875
if nargin > 8
    Summ = max((Summ + varargin{1} .* varargin{2}) ./ (D + varargin{1} .* varargin{2}));
else
    Summ = max(Summ ./ D);
end
if is_transposed
    if nargin > 8
        im = im + (1/Summ).*(im./(D + varargin{1} .* varargin{2})).*(A*(uu./(A'*im + epps + SinDelayed) - 1) - varargin{1} .* varargin{2});
    else
        im = im + (1/Summ).*(im./(D)).*(A*(uu./(A'*im + epps + SinDelayed) - 1));
    end
else
    if nargin > 8
        im = im + (1/Summ).*(im./(D + varargin{1} .* varargin{2})).*(A'*(uu./(A*im + epps + SinDelayed) - 1) - varargin{1} .* varargin{2});
    else
        im = im + (1/Summ).*(im./(D)).*(A'*(uu./(A*im + epps + SinDelayed) - 1));
    end
end

%%%% This is the original version
% Summ = max(Summ + varargin{1} .* varargin{2});
% if is_transposed
%     if nargin > 8
%         im = im + (im./Summ).*(A*(uu./(A'*im + epps + SinDelayed) - 1) - varargin{1} .* varargin{2});
%     else
%         im = im + (im./Summ).*(A*(uu./(A'*im + epps + SinDelayed) - 1));
%     end
% else
%     if nargin > 8
%         im = im + (im./Summ).*(A'*(uu./(A*im + epps + SinDelayed) - 1) - varargin{1} .* varargin{2});
%     else
%         im = im + (im./Summ).*(A'*(uu./(A*im + epps + SinDelayed) - 1));
%     end
% end