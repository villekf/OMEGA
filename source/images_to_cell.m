function [pz, varargout] = images_to_cell(im_vectors, llo, pz, options)
%IMAGES_TO_CELL Save the images at current time-step to a cell vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019 Ville-Veikko Wettenhovi
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
gg = 1;
if options.mlem
    pz{gg, llo} = im_vectors.MLEM;
end
gg = gg + 1;
if options.osem
    pz{gg, llo} = im_vectors.OSEM;
end
gg = gg + 1;
if options.mramla
    pz{gg, llo} = im_vectors.MRAMLA;
end
gg = gg + 1;
if options.ramla
    pz{gg, llo} = im_vectors.RAMLA;
end
gg = gg + 1;
if options.rosem
    pz{gg, llo} = im_vectors.ROSEM;
end
gg = gg + 1;
if options.rbi
    pz{gg, llo} = im_vectors.RBI;
end
gg = gg + 1;
if options.drama
    pz{gg, llo} = im_vectors.DRAMA;
end
gg = gg + 1;
if options.cosem
    pz{gg, llo} = im_vectors.COSEM;
end
gg = gg + 1;
if options.ecosem
    pz{gg, llo} = im_vectors.ECOSEM;
end
gg = gg + 1;
if options.acosem
    pz{gg, llo} = im_vectors.ACOSEM;
end

gg = gg + 1;
if options.MRP && options.OSL_OSEM
    pz{gg, llo} = im_vectors.MRP_OSL;
end
gg = gg + 1;
if options.MRP && options.OSL_MLEM
    pz{gg, llo} = im_vectors.MRP_MLEM;
end
gg = gg + 1;
if options.MRP && options.BSREM
    pz{gg, llo} = im_vectors.MRP_BSREM;
end
gg = gg + 1;
if options.MRP && options.MBSREM
    pz{gg, llo} = im_vectors.MRP_MBSREM;
end
gg = gg + 1;
if options.MRP && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.MRP_ROSEM;
end
gg = gg + 1;
if options.MRP && options.RBI_OSL
    pz{gg, llo} = im_vectors.MRP_RBI;
end
gg = gg + 1;
if options.MRP && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.MRP_COSEM;
end

gg = gg + 1;
if options.quad && options.OSL_OSEM
    pz{gg, llo} = im_vectors.Quad_OSL;
end
gg = gg + 1;
if options.quad && options.OSL_MLEM
    pz{gg, llo} = im_vectors.Quad_MLEM;
end
gg = gg + 1;
if options.quad && options.BSREM
    pz{gg, llo} = im_vectors.Quad_BSREM;
end
gg = gg + 1;
if options.quad && options.MBSREM
    pz{gg, llo} = im_vectors.Quad_MBSREM;
end
gg = gg + 1;
if options.quad && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.Quad_ROSEM;
end
gg = gg + 1;
if options.quad && options.RBI_OSL
    pz{gg, llo} = im_vectors.Quad_RBI;
end
gg = gg + 1;
if options.quad && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.Quad_COSEM;
end

gg = gg + 1;
if options.Huber && options.OSL_OSEM
    pz{gg, llo} = im_vectors.Huber_OSL;
end
gg = gg + 1;
if options.Huber && options.OSL_MLEM
    pz{gg, llo} = im_vectors.Huber_MLEM;
end
gg = gg + 1;
if options.Huber && options.BSREM
    pz{gg, llo} = im_vectors.Huber_BSREM;
end
gg = gg + 1;
if options.Huber && options.MBSREM
    pz{gg, llo} = im_vectors.Huber_MBSREM;
end
gg = gg + 1;
if options.Huber && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.Huber_ROSEM;
end
gg = gg + 1;
if options.Huber && options.RBI_OSL
    pz{gg, llo} = im_vectors.Huber_RBI;
end
gg = gg + 1;
if options.Huber && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.Huber_COSEM;
end

gg = gg + 1;
if options.L && options.OSL_OSEM
    pz{gg, llo} = im_vectors.L_OSL;
end
gg = gg + 1;
if options.L && options.OSL_MLEM
    pz{gg, llo} = im_vectors.L_MLEM;
end
gg = gg + 1;
if options.L && options.BSREM
    pz{gg, llo} = im_vectors.L_BSREM;
end
gg = gg + 1;
if options.L && options.MBSREM
    pz{gg, llo} = im_vectors.L_MBSREM;
end
gg = gg + 1;
if options.L && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.L_ROSEM;
end
gg = gg + 1;
if options.L && options.RBI_OSL
    pz{gg, llo} = im_vectors.L_RBI;
end
gg = gg + 1;
if options.L && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.L_COSEM;
end

gg = gg + 1;
if options.FMH && options.OSL_OSEM
    pz{gg, llo} = im_vectors.FMH_OSL;
end
gg = gg + 1;
if options.FMH && options.OSL_MLEM
    pz{gg, llo} = im_vectors.FMH_MLEM;
end
gg = gg + 1;
if options.FMH && options.BSREM
    pz{gg, llo} = im_vectors.FMH_BSREM;
end
gg = gg + 1;
if options.FMH && options.MBSREM
    pz{gg, llo} = im_vectors.FMH_MBSREM;
end
gg = gg + 1;
if options.FMH && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.FMH_ROSEM;
end
gg = gg + 1;
if options.FMH && options.RBI_OSL
    pz{gg, llo} = im_vectors.FMH_RBI;
end
gg = gg + 1;
if options.FMH && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.FMH_COSEM;
end

gg = gg + 1;
if options.weighted_mean && options.OSL_OSEM
    pz{gg, llo} = im_vectors.Weighted_OSL;
end
gg = gg + 1;
if options.weighted_mean && options.OSL_MLEM
    pz{gg, llo} = im_vectors.Weighted_MLEM;
end
gg = gg + 1;
if options.weighted_mean && options.BSREM
    pz{gg, llo} = im_vectors.Weighted_BSREM;
end
gg = gg + 1;
if options.weighted_mean && options.MBSREM
    pz{gg, llo} = im_vectors.Weighted_MBSREM;
end
gg = gg + 1;
if options.weighted_mean && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.Weighted_ROSEM;
end
gg = gg + 1;
if options.weighted_mean && options.RBI_OSL
    pz{gg, llo} = im_vectors.Weighted_RBI;
end
gg = gg + 1;
if options.weighted_mean && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.Weighted_COSEM;
end

gg = gg + 1;
if options.TV && options.OSL_OSEM
    pz{gg, llo} = im_vectors.TV_OSL;
end
gg = gg + 1;
if options.TV && options.OSL_MLEM
    pz{gg, llo} = im_vectors.TV_MLEM;
end
gg = gg + 1;
if options.TV && options.BSREM
    pz{gg, llo} = im_vectors.TV_BSREM;
end
gg = gg + 1;
if options.TV && options.MBSREM
    pz{gg, llo} = im_vectors.TV_MBSREM;
end
gg = gg + 1;
if options.TV && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.TV_ROSEM;
end
gg = gg + 1;
if options.TV && options.RBI_OSL
    pz{gg, llo} = im_vectors.TV_RBI;
end
gg = gg + 1;
if options.TV && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.TV_COSEM;
end

gg = gg + 1;
if options.AD && options.OSL_OSEM
    pz{gg, llo} = im_vectors.AD_OSL;
end
gg = gg + 1;
if options.AD && options.OSL_MLEM
    pz{gg, llo} = im_vectors.AD_MLEM;
end
gg = gg + 1;
if options.AD && options.BSREM
    pz{gg, llo} = im_vectors.AD_BSREM;
end
gg = gg + 1;
if options.AD && options.MBSREM
    pz{gg, llo} = im_vectors.AD_MBSREM;
end
gg = gg + 1;
if options.AD && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.AD_ROSEM;
end
gg = gg + 1;
if options.AD && options.RBI_OSL
    pz{gg, llo} = im_vectors.AD_RBI;
end
gg = gg + 1;
if options.AD && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.AD_COSEM;
end

gg = gg + 1;
if options.APLS && options.OSL_OSEM
    pz{gg, llo} = im_vectors.APLS_OSL;
end
gg = gg + 1;
if options.APLS && options.OSL_MLEM
    pz{gg, llo} = im_vectors.APLS_MLEM;
end
gg = gg + 1;
if options.APLS && options.BSREM
    pz{gg, llo} = im_vectors.APLS_BSREM;
end
gg = gg + 1;
if options.APLS && options.MBSREM
    pz{gg, llo} = im_vectors.APLS_MBSREM;
end
gg = gg + 1;
if options.APLS && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.APLS_ROSEM;
end
gg = gg + 1;
if options.APLS && options.RBI_OSL
    pz{gg, llo} = im_vectors.APLS_RBI;
end
gg = gg + 1;
if options.APLS && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.APLS_COSEM;
end

gg = gg + 1;
if options.TGV && options.OSL_OSEM
    pz{gg, llo} = im_vectors.TGV_OSL;
end
gg = gg + 1;
if options.TGV && options.OSL_MLEM
    pz{gg, llo} = im_vectors.TGV_MLEM;
end
gg = gg + 1;
if options.TGV && options.BSREM
    pz{gg, llo} = im_vectors.TGV_BSREM;
end
gg = gg + 1;
if options.TGV && options.MBSREM
    pz{gg, llo} = im_vectors.TGV_MBSREM;
end
gg = gg + 1;
if options.TGV && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.TGV_ROSEM;
end
gg = gg + 1;
if options.TGV && options.RBI_OSL
    pz{gg, llo} = im_vectors.TGV_RBI;
end
gg = gg + 1;
if options.TGV && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.TGV_COSEM;
end

gg = gg + 1;
if options.NLM && options.OSL_OSEM
    pz{gg, llo} = im_vectors.NLM_OSL;
end
gg = gg + 1;
if options.NLM && options.OSL_MLEM
    pz{gg, llo} = im_vectors.NLM_MLEM;
end
gg = gg + 1;
if options.NLM && options.BSREM
    pz{gg, llo} = im_vectors.NLM_BSREM;
end
gg = gg + 1;
if options.NLM && options.MBSREM
    pz{gg, llo} = im_vectors.NLM_MBSREM;
end
gg = gg + 1;
if options.NLM && options.ROSEM_MAP
    pz{gg, llo} = im_vectors.NLM_ROSEM;
end
gg = gg + 1;
if options.NLM && options.RBI_OSL
    pz{gg, llo} = im_vectors.NLM_RBI;
end
gg = gg + 1;
if options.NLM && any(options.COSEM_OSL)
    pz{gg, llo} = im_vectors.NLM_COSEM;
end
if nargout > 1
    varargout{1} = gg;
end