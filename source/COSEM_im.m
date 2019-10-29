function [im, C_co] = COSEM_im(im, A, epps, uu, C_co, D, osa_iter, SinDelayed, is_transposed)
%COSEM_IM Computes the Complete-data OSEM (COSEM) estimates
%
% Example:
%   [im, C_co] = COSEM_im(im, A, epps, uu, C_co, D, osa_iter, SinDelayed, is_transposed)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   epps = Small constant to prevent division by zero
%   uu = Measurements at current subset
%   C_co = Complete-data matrix
%   D = Sum of the complete data system matrix (D = sum(B,2), where B =
%   [A_1,A_2,...,A_subsets]
%   osa_iter = Current subset (sub-iteration)
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero. 
%   is_transposed = true if A matrix is the transpose of it, false if not
%
% OUTPUTS:
%   im = The updated estimate
%   C_co = Updated complete-data matrix
%
%   See also ECOSEM_im, ACOSEM_im, COSEM_OSL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
% C_co(:,osa_iter) = full(sum(bsxfun(@times,uu./(A'*im + epps),bsxfun(@times,A,im))');
% apu = C_co(:,osa_iter);
if is_transposed
    C_co(:,osa_iter) = full(sum((spdiags(im,0,size(A,1),size(A,1))*A) * spdiags(uu./(A'*im + epps + SinDelayed),0,size(A,2),size(A,2)),2));
else
    C_co(:,osa_iter) = full(sum((spdiags(im,0,size(A,2),size(A,2))*A') * spdiags(uu./(A*im + epps + SinDelayed),0,size(A,1),size(A,1)),2));
end
% C_co(:,osa_iter) = full(sum((spdiags(im,0,size(A,1),size(A,1))*A),2));
% B = B + C_co(:,osa_iter) - apu;
% im = (B)./D;
% im = (apu)*sum(uu)/sum(A'*apu+epps);
im = (sum(C_co,2)./D);