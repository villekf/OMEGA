function [im, C_aco] = ACOSEM_im(im, A, epps, uu, C_aco, D, h, osa_iter, SinDelayed, is_transposed)
%ACOSEM_IM Computes the Accelerated COSEM (ACOSEM) estimates
%
% Example:
%   [im, C_aco] = ACOSEM_im(im, A, epps, uu, C_aco, D, h, osa_iter, SinDelayed, is_transposed)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   epps = Small constant to prevent division by zero
%   uu = Measurements at current subset
%   C_aco = Complete-data matrix for ACOSEM
%   D = Sum of the complete data system matrix (D = sum(B,2), where B =
%   [A_1,A_2,...,A_subsets]
%   h = Acceleration parameter
%   osa_iter = Current subset (sub-iteration)
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero. 
%   is_transposed = true if A matrix is the transpose of it, false if not
%
% OUTPUTS:
%   im = The updated estimate
%   C_aco = Updated complete-data matrix
%
%   See also COSEM_im, ECOSEM_im, COSEM_OSL

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
%C_aco(:,osa_iter) = full(sum(spdiags(uu./(A'*pz_acosem_apu + epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',pz_acosem_apu'.^(1/h))))');
if is_transposed
    C_aco(:,osa_iter) = full(sum((spdiags(im.^(1/h),0,size(A,1),size(A,1))*A) * spdiags(uu./(A'*im + epps + SinDelayed),0,size(A,2),size(A,2)),2));
else
    C_aco(:,osa_iter) = full(sum((spdiags(im.^(1/h),0,size(A,2),size(A,2))*A') * spdiags(uu./(A*im + epps + SinDelayed),0,size(A,1),size(A,1)),2));
end
im = (sum(C_aco,2)./D).^h;
% im = sum(C_aco,2);
% disp(num2str(sum(A'*im)))
if is_transposed
    im = (im)*(sum(uu)/(sum(A'*im)));
else
    im = (im)*(sum(uu)/(sum(A*im)));
end
    im(im < 0) = epps;