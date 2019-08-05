function im_vectors = form_image_vectors(options, N)
%FORM_IMAGE_VECTORS This function simply forms the applicable image vectors
%for each reconstruction method
%
%
%
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

if options.reconstruction_method ~= 2
    
    
    MLEM_bool = options.OSL_MLEM || options.mlem;
    OS_bool = options.osem || options.rosem || options.ramla || options.OSL_OSEM || options.BSREM || options.ROSEM_MAP;
    
    if options.reconstruction_method == 4 && OS_bool
        im_vectors.OSEM_apu = ones(N, 1);
    end
    if options.reconstruction_method == 4 && MLEM_bool
        im_vectors.MLEM_apu = ones(N, 1);
    end
    if options.osem
        im_vectors.OSEM = ones(N,options.Niter + 1);
        im_vectors.OSEM(:,1) = options.x0(:);
        im_vectors.OSEM_apu = im_vectors.OSEM(:,1);
    end
    if options.mlem
        im_vectors.MLEM = ones(N,options.Niter + 1);
        im_vectors.MLEM(:,1) = options.x0(:);
        im_vectors.MLEM_apu = im_vectors.MLEM(:,1);
    end
    if options.mramla
        im_vectors.MRAMLA = ones(N,options.Niter + 1);
        im_vectors.MRAMLA(:,1) = options.x0(:);
        im_vectors.MRAMLA_apu = im_vectors.MRAMLA(:,1);
    end
    if options.ramla
        im_vectors.RAMLA = ones(N,options.Niter + 1);
        im_vectors.RAMLA(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.RAMLA_apu = im_vectors.RAMLA(:,1);
        end
    end
    if options.rosem
        im_vectors.ROSEM = ones(N,options.Niter + 1);
        im_vectors.ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.ROSEM_apu = im_vectors.ROSEM(:,1);
        end
    end
    if options.rbi
        im_vectors.RBI = ones(N,options.Niter + 1);
        im_vectors.RBI(:,1) = options.x0(:);
        im_vectors.RBI_apu = im_vectors.RBI(:,1);
    end
    if options.drama
        im_vectors.DRAMA = ones(N,options.Niter + 1);
        im_vectors.DRAMA(:,1) = options.x0(:);
        im_vectors.DRAMA_apu = im_vectors.DRAMA(:,1);
    end
    if options.cosem
        im_vectors.COSEM = ones(N,options.Niter + 1);
        im_vectors.COSEM(:,1) = options.x0(:);
        im_vectors.COSEM_apu = im_vectors.COSEM(:,1);
    end
    if options.ecosem
        im_vectors.ECOSEM = ones(N,options.Niter + 1);
        im_vectors.ECOSEM(:,1) = options.x0(:);
        im_vectors.ECOSEM_apu = im_vectors.ECOSEM(:,1);
        if ~options.osem
            im_vectors.OSEM_apu = options.x0(:);
        end
        if ~options.cosem
            im_vectors.COSEM_apu = options.x0(:);
        end
    end
    if options.acosem
        im_vectors.ACOSEM = ones(N,options.Niter + 1);
        im_vectors.ACOSEM(:,1) = options.x0(:);
        im_vectors.ACOSEM_apu = im_vectors.ACOSEM(:,1);
    end
    
    if options.MRP && options.OSL_OSEM
        im_vectors.MRP_OSL = ones(N,options.Niter + 1);
        im_vectors.MRP_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.MRP_OSL_apu = im_vectors.MRP_OSL(:,1);
        end
    end
    if options.MRP && options.OSL_MLEM
        im_vectors.MRP_MLEM = ones(N,options.Niter + 1);
        im_vectors.MRP_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.MRP_MLEM_apu = im_vectors.MRP_MLEM(:,1);
        end
    end
    if options.MRP && options.MBSREM
        im_vectors.MRP_MBSREM = ones(N,options.Niter + 1);
        im_vectors.MRP_MBSREM(:,1) = options.x0(:);
        im_vectors.MRP_MBSREM_apu = im_vectors.MRP_MBSREM(:,1);
    end
    if options.MRP && options.BSREM
        im_vectors.MRP_BSREM = ones(N,options.Niter + 1);
        im_vectors.MRP_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.MRP_BSREM_apu = im_vectors.MRP_BSREM(:,1);
        end
    end
    if options.MRP && options.ROSEM_MAP
        im_vectors.MRP_ROSEM = ones(N,options.Niter + 1);
        im_vectors.MRP_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.MRP_ROSEM_apu = im_vectors.MRP_ROSEM(:,1);
        end
    end
    if options.MRP && options.RBI_MAP
        im_vectors.MRP_RBI = ones(N,options.Niter + 1);
        im_vectors.MRP_RBI(:,1) = options.x0(:);
        im_vectors.MRP_RBI_apu = im_vectors.MRP_RBI(:,1);
    end
    if options.MRP && any(options.COSEM_MAP)
        im_vectors.MRP_COSEM = ones(N,options.Niter + 1);
        im_vectors.MRP_COSEM(:,1) = options.x0(:);
        im_vectors.MRP_COSEM_apu = im_vectors.MRP_COSEM(:,1);
    end
    
    if options.quad && options.OSL_OSEM
        im_vectors.Quad_OSL = ones(N,options.Niter + 1);
        im_vectors.Quad_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Quad_OSL_apu = im_vectors.Quad_OSL(:,1);
        end
    end
    if options.quad && options.OSL_MLEM
        im_vectors.Quad_MLEM = ones(N,options.Niter + 1);
        im_vectors.Quad_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Quad_MLEM_apu = im_vectors.Quad_MLEM(:,1);
        end
    end
    if options.quad && options.MBSREM
        im_vectors.Quad_MBSREM = ones(N,options.Niter + 1);
        im_vectors.Quad_MBSREM(:,1) = options.x0(:);
        im_vectors.Quad_MBSREM_apu = im_vectors.Quad_MBSREM(:,1);
    end
    if options.quad && options.BSREM
        im_vectors.Quad_BSREM = ones(N,options.Niter + 1);
        im_vectors.Quad_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Quad_BSREM_apu = im_vectors.Quad_BSREM(:,1);
        end
    end
    if options.quad && options.ROSEM_MAP
        im_vectors.Quad_ROSEM = ones(N,options.Niter + 1);
        im_vectors.Quad_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Quad_ROSEM_apu = im_vectors.Quad_ROSEM(:,1);
        end
    end
    if options.quad && options.RBI_MAP
        im_vectors.Quad_RBI = ones(N,options.Niter + 1);
        im_vectors.Quad_RBI(:,1) = options.x0(:);
        im_vectors.Quad_RBI_apu = im_vectors.Quad_RBI(:,1);
    end
    if options.quad && any(options.COSEM_MAP)
        im_vectors.Quad_COSEM = ones(N,options.Niter + 1);
        im_vectors.Quad_COSEM(:,1) = options.x0(:);
        im_vectors.Quad_COSEM_apu = im_vectors.Quad_COSEM(:,1);
    end
    
    if options.L && options.OSL_OSEM
        im_vectors.L_OSL = ones(N,options.Niter + 1);
        im_vectors.L_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.L_OSL_apu = im_vectors.L_OSL(:,1);
        end
    end
    if options.L && options.OSL_MLEM
        im_vectors.L_MLEM = ones(N,options.Niter + 1);
        im_vectors.L_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.L_MLEM_apu = im_vectors.L_MLEM(:,1);
        end
    end
    if options.L && options.MBSREM
        im_vectors.L_MBSREM = ones(N,options.Niter + 1);
        im_vectors.L_MBSREM(:,1) = options.x0(:);
        im_vectors.L_MBSREM_apu = im_vectors.L_MBSREM(:,1);
    end
    if options.L && options.BSREM
        im_vectors.L_BSREM = ones(N,options.Niter + 1);
        im_vectors.L_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.L_BSREM_apu = im_vectors.L_BSREM(:,1);
        end
    end
    if options.L && options.ROSEM_MAP
        im_vectors.L_ROSEM = ones(N,options.Niter + 1);
        im_vectors.L_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.L_ROSEM_apu = im_vectors.L_ROSEM(:,1);
        end
    end
    if options.L && options.RBI_MAP
        im_vectors.L_RBI = ones(N,options.Niter + 1);
        im_vectors.L_RBI(:,1) = options.x0(:);
        im_vectors.L_RBI_apu = im_vectors.L_RBI(:,1);
    end
    if options.L && any(options.COSEM_MAP)
        im_vectors.L_COSEM = ones(N,options.Niter + 1);
        im_vectors.L_COSEM(:,1) = options.x0(:);
        im_vectors.L_COSEM_apu = im_vectors.L_COSEM(:,1);
    end
    
    if options.FMH && options.OSL_OSEM
        im_vectors.FMH_OSL = ones(N,options.Niter + 1);
        im_vectors.FMH_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.FMH_OSL_apu = im_vectors.FMH_OSL(:,1);
        end
    end
    if options.FMH && options.OSL_MLEM
        im_vectors.FMH_MLEM = ones(N,options.Niter + 1);
        im_vectors.FMH_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.FMH_MLEM_apu = im_vectors.FMH_MLEM(:,1);
        end
    end
    if options.FMH && options.MBSREM
        im_vectors.FMH_MBSREM = ones(N,options.Niter + 1);
        im_vectors.FMH_MBSREM(:,1) = options.x0(:);
        im_vectors.FMH_MBSREM_apu = im_vectors.FMH_MBSREM(:,1);
    end
    if options.FMH && options.BSREM
        im_vectors.FMH_BSREM = ones(N,options.Niter + 1);
        im_vectors.FMH_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.FMH_BSREM_apu = im_vectors.FMH_BSREM(:,1);
        end
    end
    if options.FMH && options.ROSEM_MAP
        im_vectors.FMH_ROSEM = ones(N,options.Niter + 1);
        im_vectors.FMH_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.FMH_ROSEM_apu = im_vectors.FMH_ROSEM(:,1);
        end
    end
    if options.FMH && options.RBI_MAP
        im_vectors.FMH_RBI = ones(N,options.Niter + 1);
        im_vectors.FMH_RBI(:,1) = options.x0(:);
        im_vectors.FMH_RBI_apu = im_vectors.FMH_RBI(:,1);
    end
    if options.FMH && any(options.COSEM_MAP)
        im_vectors.FMH_COSEM = ones(N,options.Niter + 1);
        im_vectors.FMH_COSEM(:,1) = options.x0(:);
        im_vectors.FMH_COSEM_apu = im_vectors.FMH_COSEM(:,1);
    end
    
    if options.weighted_mean && options.OSL_OSEM
        im_vectors.Weighted_OSL = ones(N,options.Niter + 1);
        im_vectors.Weighted_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Weighted_OSL_apu = im_vectors.Weighted_OSL(:,1);
        end
    end
    if options.weighted_mean && options.OSL_MLEM
        im_vectors.Weighted_MLEM = ones(N,options.Niter + 1);
        im_vectors.Weighted_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Weighted_MLEM_apu = im_vectors.Weighted_MLEM(:,1);
        end
    end
    if options.weighted_mean && options.MBSREM
        im_vectors.Weighted_MBSREM = ones(N,options.Niter + 1);
        im_vectors.Weighted_MBSREM(:,1) = options.x0(:);
        im_vectors.Weighted_MBSREM_apu = im_vectors.Weighted_MBSREM(:,1);
    end
    if options.weighted_mean && options.BSREM
        im_vectors.Weighted_BSREM = ones(N,options.Niter + 1);
        im_vectors.Weighted_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Weighted_BSREM_apu = im_vectors.Weighted_BSREM(:,1);
        end
    end
    if options.weighted_mean && options.ROSEM_MAP
        im_vectors.Weighted_ROSEM = ones(N,options.Niter + 1);
        im_vectors.Weighted_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.Weighted_ROSEM_apu = im_vectors.Weighted_ROSEM(:,1);
        end
    end
    if options.weighted_mean && options.RBI_MAP
        im_vectors.Weighted_RBI = ones(N,options.Niter + 1);
        im_vectors.Weighted_RBI(:,1) = options.x0(:);
        im_vectors.Weighted_RBI_apu = im_vectors.Weighted_RBI(:,1);
    end
    if options.weighted_mean && any(options.COSEM_MAP)
        im_vectors.Weighted_COSEM = ones(N,options.Niter + 1);
        im_vectors.Weighted_COSEM(:,1) = options.x0(:);
        im_vectors.Weighted_COSEM_apu = im_vectors.Weighted_COSEM(:,1);
    end
    
    if options.TV && options.OSL_OSEM
        im_vectors.TV_OSL = ones(N,options.Niter + 1);
        im_vectors.TV_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TV_OSL_apu = im_vectors.TV_OSL(:,1);
        end
    end
    if options.TV && options.OSL_MLEM
        im_vectors.TV_MLEM = ones(N,options.Niter + 1);
        im_vectors.TV_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TV_MLEM_apu = im_vectors.TV_MLEM(:,1);
        end
    end
    if options.TV && options.MBSREM
        im_vectors.TV_MBSREM = ones(N,options.Niter + 1);
        im_vectors.TV_MBSREM(:,1) = options.x0(:);
        im_vectors.TV_MBSREM_apu = im_vectors.TV_MBSREM(:,1);
    end
    if options.TV && options.BSREM
        im_vectors.TV_BSREM = ones(N,options.Niter + 1);
        im_vectors.TV_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TV_BSREM_apu = im_vectors.TV_BSREM(:,1);
        end
    end
    if options.TV && options.ROSEM_MAP
        im_vectors.TV_ROSEM = ones(N,options.Niter + 1);
        im_vectors.TV_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TV_ROSEM_apu = im_vectors.TV_ROSEM(:,1);
        end
    end
    if options.TV && options.RBI_MAP
        im_vectors.TV_RBI = ones(N,options.Niter + 1);
        im_vectors.TV_RBI(:,1) = options.x0(:);
        im_vectors.TV_RBI_apu = im_vectors.TV_RBI(:,1);
    end
    if options.TV && any(options.COSEM_MAP)
        im_vectors.TV_COSEM = ones(N,options.Niter + 1);
        im_vectors.TV_COSEM(:,1) = options.x0(:);
        im_vectors.TV_COSEM_apu = im_vectors.TV_COSEM(:,1);
    end
    
    if options.AD && options.OSL_OSEM
        im_vectors.AD_OSL = ones(N,options.Niter + 1);
        im_vectors.AD_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.AD_OSL_apu = im_vectors.AD_OSL(:,1);
        end
    end
    if options.AD && options.OSL_MLEM
        im_vectors.AD_MLEM = ones(N,options.Niter + 1);
        im_vectors.AD_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.AD_MLEM_apu = im_vectors.AD_MLEM(:,1);
        end
    end
    if options.AD && options.MBSREM
        im_vectors.AD_MBSREM = ones(N,options.Niter + 1);
        im_vectors.AD_MBSREM(:,1) = options.x0(:);
        im_vectors.AD_MBSREM_apu = im_vectors.AD_MBSREM(:,1);
    end
    if options.AD && options.BSREM
        im_vectors.AD_BSREM = ones(N,options.Niter + 1);
        im_vectors.AD_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.AD_BSREM_apu = im_vectors.AD_BSREM(:,1);
        end
    end
    if options.AD && options.ROSEM_MAP
        im_vectors.AD_ROSEM = ones(N,options.Niter + 1);
        im_vectors.AD_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.AD_ROSEM_apu = im_vectors.AD_ROSEM(:,1);
        end
    end
    if options.AD && options.RBI_MAP
        im_vectors.AD_RBI = ones(N,options.Niter + 1);
        im_vectors.AD_RBI(:,1) = options.x0(:);
        im_vectors.AD_RBI_apu = im_vectors.AD_RBI(:,1);
    end
    if options.AD && any(options.COSEM_MAP)
        im_vectors.AD_COSEM = ones(N,options.Niter + 1);
        im_vectors.AD_COSEM(:,1) = options.x0(:);
        im_vectors.AD_COSEM_apu = im_vectors.AD_COSEM(:,1);
    end
    
    if options.APLS && options.OSL_OSEM
        im_vectors.APLS_OSL = ones(N,options.Niter + 1);
        im_vectors.APLS_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.APLS_OSL_apu = im_vectors.APLS_OSL(:,1);
        end
    end
    if options.APLS && options.OSL_MLEM
        im_vectors.APLS_MLEM = ones(N,options.Niter + 1);
        im_vectors.APLS_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.APLS_MLEM_apu = im_vectors.APLS_MLEM(:,1);
        end
    end
    if options.APLS && options.MBSREM
        im_vectors.APLS_MBSREM = ones(N,options.Niter + 1);
        im_vectors.APLS_MBSREM(:,1) = options.x0(:);
        im_vectors.APLS_MBSREM_apu = im_vectors.APLS_MBSREM(:,1);
    end
    if options.APLS && options.BSREM
        im_vectors.APLS_BSREM = ones(N,options.Niter + 1);
        im_vectors.APLS_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.APLS_BSREM_apu = im_vectors.APLS_BSREM(:,1);
        end
    end
    if options.APLS && options.ROSEM_MAP
        im_vectors.APLS_ROSEM = ones(N,options.Niter + 1);
        im_vectors.APLS_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.APLS_ROSEM_apu = im_vectors.APLS_ROSEM(:,1);
        end
    end
    if options.APLS && options.RBI_MAP
        im_vectors.APLS_RBI = ones(N,options.Niter + 1);
        im_vectors.APLS_RBI(:,1) = options.x0(:);
        im_vectors.APLS_RBI_apu = im_vectors.APLS_RBI(:,1);
    end
    if options.APLS && any(options.COSEM_MAP)
        im_vectors.APLS_COSEM = ones(N,options.Niter + 1);
        im_vectors.APLS_COSEM(:,1) = options.x0(:);
        im_vectors.APLS_COSEM_apu = im_vectors.APLS_COSEM(:,1);
    end
    
    if options.TGV && options.OSL_OSEM
        im_vectors.TGV_OSL = ones(N,options.Niter + 1);
        im_vectors.TGV_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TGV_OSL_apu = im_vectors.TGV_OSL(:,1);
        end
    end
    if options.TGV && options.OSL_MLEM
        im_vectors.TGV_MLEM = ones(N,options.Niter + 1);
        im_vectors.TGV_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TGV_MLEM_apu = im_vectors.TGV_MLEM(:,1);
        end
    end
    if options.TGV && options.MBSREM
        im_vectors.TGV_MBSREM = ones(N,options.Niter + 1);
        im_vectors.TGV_MBSREM(:,1) = options.x0(:);
        im_vectors.TGV_MBSREM_apu = im_vectors.TGV_MBSREM(:,1);
    end
    if options.TGV && options.BSREM
        im_vectors.TGV_BSREM = ones(N,options.Niter + 1);
        im_vectors.TGV_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TGV_BSREM_apu = im_vectors.TGV_BSREM(:,1);
        end
    end
    if options.TGV && options.ROSEM_MAP
        im_vectors.TGV_ROSEM = ones(N,options.Niter + 1);
        im_vectors.TGV_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.TGV_ROSEM_apu = im_vectors.TGV_ROSEM(:,1);
        end
    end
    if options.TGV && options.RBI_MAP
        im_vectors.TGV_RBI = ones(N,options.Niter + 1);
        im_vectors.TGV_RBI(:,1) = options.x0(:);
        im_vectors.TGV_RBI_apu = im_vectors.TGV_RBI(:,1);
    end
    if options.TGV && any(options.COSEM_MAP)
        im_vectors.TGV_COSEM = ones(N,options.Niter + 1);
        im_vectors.TGV_COSEM(:,1) = options.x0(:);
        im_vectors.TGV_COSEM_apu = im_vectors.TGV_COSEM(:,1);
    end
    
    if options.NLM && options.OSL_OSEM
        im_vectors.NLM_OSL = ones(N,options.Niter + 1);
        im_vectors.NLM_OSL(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.NLM_OSL_apu = im_vectors.NLM_OSL(:,1);
        end
    end
    if options.NLM && options.OSL_MLEM
        im_vectors.NLM_MLEM = ones(N,options.Niter + 1);
        im_vectors.NLM_MLEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.NLM_MLEM_apu = im_vectors.NLM_MLEM(:,1);
        end
    end
    if options.NLM && options.MBSREM
        im_vectors.NLM_MBSREM = ones(N,options.Niter + 1);
        im_vectors.NLM_MBSREM(:,1) = options.x0(:);
        im_vectors.NLM_MBSREM_apu = im_vectors.NLM_MBSREM(:,1);
    end
    if options.NLM && options.BSREM
        im_vectors.NLM_BSREM = ones(N,options.Niter + 1);
        im_vectors.NLM_BSREM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.NLM_BSREM_apu = im_vectors.NLM_BSREM(:,1);
        end
    end
    if options.NLM && options.ROSEM_MAP
        im_vectors.NLM_ROSEM = ones(N,options.Niter + 1);
        im_vectors.NLM_ROSEM(:,1) = options.x0(:);
        if ~options.reconstruction_method == 4
            im_vectors.NLM_ROSEM_apu = im_vectors.NLM_ROSEM(:,1);
        end
    end
    if options.NLM && options.RBI_MAP
        im_vectors.NLM_RBI = ones(N,options.Niter + 1);
        im_vectors.NLM_RBI(:,1) = options.x0(:);
        im_vectors.NLM_RBI_apu = im_vectors.NLM_RBI(:,1);
    end
    if options.NLM && any(options.COSEM_MAP)
        im_vectors.NLM_COSEM = ones(N,options.Niter + 1);
        im_vectors.NLM_COSEM(:,1) = options.x0(:);
        im_vectors.NLM_COSEM_apu = im_vectors.NLM_COSEM(:,1);
    end
    
else
    if options.osem
        im_vectors.OSEM = ones(N,options.Niter + 1,'single');
        im_vectors.OSEM(:,1) = options.x0(:);
        im_vectors.OSEM_apu = im_vectors.OSEM(:,1);
    end
    if options.mlem
        im_vectors.MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.MLEM(:,1) = options.x0(:);
        im_vectors.MLEM_apu = im_vectors.MLEM(:,1);
    end
    if options.mramla
        im_vectors.MRAMLA = ones(N,options.Niter + 1,'single');
        im_vectors.MRAMLA(:,1) = options.x0(:);
        im_vectors.MRAMLA_apu = im_vectors.MRAMLA(:,1);
    end
    if options.ramla
        im_vectors.RAMLA = ones(N,options.Niter + 1,'single');
        im_vectors.RAMLA(:,1) = options.x0(:);
        im_vectors.RAMLA_apu = im_vectors.RAMLA(:,1);
    end
    if options.rosem
        im_vectors.ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.ROSEM(:,1) = options.x0(:);
        im_vectors.ROSEM_apu = im_vectors.ROSEM(:,1);
    end
    if options.rbi
        im_vectors.RBI = ones(N,options.Niter + 1,'single');
        im_vectors.RBI(:,1) = options.x0(:);
        im_vectors.RBI_apu = im_vectors.RBI(:,1);
    end
    if options.drama
        im_vectors.DRAMA = ones(N,options.Niter + 1,'single');
        im_vectors.DRAMA(:,1) = options.x0(:);
        im_vectors.DRAMA_apu = im_vectors.DRAMA(:,1);
    end
    if options.cosem
        im_vectors.COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.COSEM(:,1) = options.x0(:);
        im_vectors.COSEM_apu = im_vectors.COSEM(:,1);
    end
    if options.ecosem
        im_vectors.ECOSEM = ones(N,options.Niter + 1,'single');
        im_vectors.ECOSEM(:,1) = options.x0(:);
        im_vectors.ECOSEM_apu = im_vectors.ECOSEM(:,1);
        if ~options.osem
            im_vectors.OSEM_apu = options.x0(:);
        end
        if ~options.cosem
            im_vectors.COSEM_apu = options.x0(:);
        end
    end
    if options.acosem
        im_vectors.ACOSEM = ones(N,options.Niter + 1,'single');
        im_vectors.ACOSEM(:,1) = options.x0(:);
        im_vectors.ACOSEM_apu = im_vectors.ACOSEM(:,1);
    end
    
    if options.MRP && options.OSL_OSEM
        im_vectors.MRP_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_OSL(:,1) = options.x0(:);
        im_vectors.MRP_OSL_apu = im_vectors.MRP_OSL(:,1);
    end
    if options.MRP && options.OSL_MLEM
        im_vectors.MRP_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_MLEM(:,1) = options.x0(:);
        im_vectors.MRP_MLEM_apu = im_vectors.MRP_MLEM(:,1);
    end
    if options.MRP && options.MBSREM
        im_vectors.MRP_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_MBSREM(:,1) = options.x0(:);
        im_vectors.MRP_MBSREM_apu = im_vectors.MRP_MBSREM(:,1);
    end
    if options.MRP && options.BSREM
        im_vectors.MRP_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_BSREM(:,1) = options.x0(:);
        im_vectors.MRP_BSREM_apu = im_vectors.MRP_BSREM(:,1);
    end
    if options.MRP && options.ROSEM_MAP
        im_vectors.MRP_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_ROSEM(:,1) = options.x0(:);
        im_vectors.MRP_ROSEM_apu = im_vectors.MRP_ROSEM(:,1);
    end
    if options.MRP && options.RBI_MAP
        im_vectors.MRP_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_RBI(:,1) = options.x0(:);
        im_vectors.MRP_RBI_apu = im_vectors.MRP_RBI(:,1);
    end
    if options.MRP && any(options.COSEM_MAP)
        im_vectors.MRP_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.MRP_COSEM(:,1) = options.x0(:);
        im_vectors.MRP_COSEM_apu = im_vectors.MRP_COSEM(:,1);
    end
    
    if options.quad && options.OSL_OSEM
        im_vectors.Quad_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_OSL(:,1) = options.x0(:);
        im_vectors.Quad_OSL_apu = im_vectors.Quad_OSL(:,1);
    end
    if options.quad && options.OSL_MLEM
        im_vectors.Quad_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_MLEM(:,1) = options.x0(:);
        im_vectors.Quad_MLEM_apu = im_vectors.Quad_MLEM(:,1);
    end
    if options.quad && options.MBSREM
        im_vectors.Quad_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_MBSREM(:,1) = options.x0(:);
        im_vectors.Quad_MBSREM_apu = im_vectors.Quad_MBSREM(:,1);
    end
    if options.quad && options.BSREM
        im_vectors.Quad_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_BSREM(:,1) = options.x0(:);
        im_vectors.Quad_BSREM_apu = im_vectors.Quad_BSREM(:,1);
    end
    if options.quad && options.ROSEM_MAP
        im_vectors.Quad_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_ROSEM(:,1) = options.x0(:);
        im_vectors.Quad_ROSEM_apu = im_vectors.Quad_ROSEM(:,1);
    end
    if options.quad && options.RBI_MAP
        im_vectors.Quad_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_RBI(:,1) = options.x0(:);
        im_vectors.Quad_RBI_apu = im_vectors.Quad_RBI(:,1);
    end
    if options.quad && any(options.COSEM_MAP)
        im_vectors.Quad_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.Quad_COSEM(:,1) = options.x0(:);
        im_vectors.Quad_COSEM_apu = im_vectors.Quad_COSEM(:,1);
    end
    
    if options.L && options.OSL_OSEM
        im_vectors.L_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.L_OSL(:,1) = options.x0(:);
        im_vectors.L_OSL_apu = im_vectors.L_OSL(:,1);
    end
    if options.L && options.OSL_MLEM
        im_vectors.L_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.L_MLEM(:,1) = options.x0(:);
        im_vectors.L_MLEM_apu = im_vectors.L_MLEM(:,1);
    end
    if options.L && options.MBSREM
        im_vectors.L_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.L_MBSREM(:,1) = options.x0(:);
        im_vectors.L_MBSREM_apu = im_vectors.L_MBSREM(:,1);
    end
    if options.L && options.BSREM
        im_vectors.L_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.L_BSREM(:,1) = options.x0(:);
        im_vectors.L_BSREM_apu = im_vectors.L_BSREM(:,1);
    end
    if options.L && options.ROSEM_MAP
        im_vectors.L_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.L_ROSEM(:,1) = options.x0(:);
        im_vectors.L_ROSEM_apu = im_vectors.L_ROSEM(:,1);
    end
    if options.L && options.RBI_MAP
        im_vectors.L_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.L_RBI(:,1) = options.x0(:);
        im_vectors.L_RBI_apu = im_vectors.L_RBI(:,1);
    end
    if options.L && any(options.COSEM_MAP)
        im_vectors.L_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.L_COSEM(:,1) = options.x0(:);
        im_vectors.L_COSEM_apu = im_vectors.L_COSEM(:,1);
    end
    
    if options.FMH && options.OSL_OSEM
        im_vectors.FMH_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_OSL(:,1) = options.x0(:);
        im_vectors.FMH_OSL_apu = im_vectors.FMH_OSL(:,1);
    end
    if options.FMH && options.OSL_MLEM
        im_vectors.FMH_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_MLEM(:,1) = options.x0(:);
        im_vectors.FMH_MLEM_apu = im_vectors.FMH_MLEM(:,1);
    end
    if options.FMH && options.MBSREM
        im_vectors.FMH_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_MBSREM(:,1) = options.x0(:);
        im_vectors.FMH_MBSREM_apu = im_vectors.FMH_MBSREM(:,1);
    end
    if options.FMH && options.BSREM
        im_vectors.FMH_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_BSREM(:,1) = options.x0(:);
        im_vectors.FMH_BSREM_apu = im_vectors.FMH_BSREM(:,1);
    end
    if options.FMH && options.ROSEM_MAP
        im_vectors.FMH_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_ROSEM(:,1) = options.x0(:);
        im_vectors.FMH_ROSEM_apu = im_vectors.FMH_ROSEM(:,1);
    end
    if options.FMH && options.RBI_MAP
        im_vectors.FMH_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_RBI(:,1) = options.x0(:);
        im_vectors.FMH_RBI_apu = im_vectors.FMH_RBI(:,1);
    end
    if options.FMH && any(options.COSEM_MAP)
        im_vectors.FMH_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.FMH_COSEM(:,1) = options.x0(:);
        im_vectors.FMH_COSEM_apu = im_vectors.FMH_COSEM(:,1);
    end
    
    if options.weighted_mean && options.OSL_OSEM
        im_vectors.Weighted_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_OSL(:,1) = options.x0(:);
        im_vectors.Weighted_OSL_apu = im_vectors.Weighted_OSL(:,1);
    end
    if options.weighted_mean && options.OSL_MLEM
        im_vectors.Weighted_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_MLEM(:,1) = options.x0(:);
        im_vectors.Weighted_MLEM_apu = im_vectors.Weighted_MLEM(:,1);
    end
    if options.weighted_mean && options.MBSREM
        im_vectors.Weighted_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_MBSREM(:,1) = options.x0(:);
        im_vectors.Weighted_MBSREM_apu = im_vectors.Weighted_MBSREM(:,1);
    end
    if options.weighted_mean && options.BSREM
        im_vectors.Weighted_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_BSREM(:,1) = options.x0(:);
        im_vectors.Weighted_BSREM_apu = im_vectors.Weighted_BSREM(:,1);
    end
    if options.weighted_mean && options.ROSEM_MAP
        im_vectors.Weighted_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_ROSEM(:,1) = options.x0(:);
        im_vectors.Weighted_ROSEM_apu = im_vectors.Weighted_ROSEM(:,1);
    end
    if options.weighted_mean && options.RBI_MAP
        im_vectors.Weighted_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_RBI(:,1) = options.x0(:);
        im_vectors.Weighted_RBI_apu = im_vectors.Weighted_RBI(:,1);
    end
    if options.weighted_mean && any(options.COSEM_MAP)
        im_vectors.Weighted_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.Weighted_COSEM(:,1) = options.x0(:);
        im_vectors.Weighted_COSEM_apu = im_vectors.Weighted_COSEM(:,1);
    end
    
    if options.TV && options.OSL_OSEM
        im_vectors.TV_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.TV_OSL(:,1) = options.x0(:);
        im_vectors.TV_OSL_apu = im_vectors.TV_OSL(:,1);
    end
    if options.TV && options.OSL_MLEM
        im_vectors.TV_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.TV_MLEM(:,1) = options.x0(:);
        im_vectors.TV_MLEM_apu = im_vectors.TV_MLEM(:,1);
    end
    if options.TV && options.MBSREM
        im_vectors.TV_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.TV_MBSREM(:,1) = options.x0(:);
        im_vectors.TV_MBSREM_apu = im_vectors.TV_MBSREM(:,1);
    end
    if options.TV && options.BSREM
        im_vectors.TV_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.TV_BSREM(:,1) = options.x0(:);
        im_vectors.TV_BSREM_apu = im_vectors.TV_BSREM(:,1);
    end
    if options.TV && options.ROSEM_MAP
        im_vectors.TV_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.TV_ROSEM(:,1) = options.x0(:);
        im_vectors.TV_ROSEM_apu = im_vectors.TV_ROSEM(:,1);
    end
    if options.TV && options.RBI_MAP
        im_vectors.TV_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.TV_RBI(:,1) = options.x0(:);
        im_vectors.TV_RBI_apu = im_vectors.TV_RBI(:,1);
    end
    if options.TV && any(options.COSEM_MAP)
        im_vectors.TV_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.TV_COSEM(:,1) = options.x0(:);
        im_vectors.TV_COSEM_apu = im_vectors.TV_COSEM(:,1);
    end
    
    if options.AD && options.OSL_OSEM
        im_vectors.AD_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.AD_OSL(:,1) = options.x0(:);
        im_vectors.AD_OSL_apu = im_vectors.AD_OSL(:,1);
    end
    if options.AD && options.OSL_MLEM
        im_vectors.AD_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.AD_MLEM(:,1) = options.x0(:);
        im_vectors.AD_MLEM_apu = im_vectors.AD_MLEM(:,1);
    end
    if options.AD && options.MBSREM
        im_vectors.AD_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.AD_MBSREM(:,1) = options.x0(:);
        im_vectors.AD_MBSREM_apu = im_vectors.AD_MBSREM(:,1);
    end
    if options.AD && options.BSREM
        im_vectors.AD_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.AD_BSREM(:,1) = options.x0(:);
        im_vectors.AD_BSREM_apu = im_vectors.AD_BSREM(:,1);
    end
    if options.AD && options.ROSEM_MAP
        im_vectors.AD_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.AD_ROSEM(:,1) = options.x0(:);
        im_vectors.AD_ROSEM_apu = im_vectors.AD_ROSEM(:,1);
    end
    if options.AD && options.RBI_MAP
        im_vectors.AD_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.AD_RBI(:,1) = options.x0(:);
        im_vectors.AD_RBI_apu = im_vectors.AD_RBI(:,1);
    end
    if options.AD && any(options.COSEM_MAP)
        im_vectors.AD_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.AD_COSEM(:,1) = options.x0(:);
        im_vectors.AD_COSEM_apu = im_vectors.AD_COSEM(:,1);
    end
    
    if options.APLS && options.OSL_OSEM
        im_vectors.APLS_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_OSL(:,1) = options.x0(:);
        im_vectors.APLS_OSL_apu = im_vectors.APLS_OSL(:,1);
    end
    if options.APLS && options.OSL_MLEM
        im_vectors.APLS_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_MLEM(:,1) = options.x0(:);
        im_vectors.APLS_MLEM_apu = im_vectors.APLS_MLEM(:,1);
    end
    if options.APLS && options.MBSREM
        im_vectors.APLS_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_MBSREM(:,1) = options.x0(:);
        im_vectors.APLS_MBSREM_apu = im_vectors.APLS_MBSREM(:,1);
    end
    if options.APLS && options.BSREM
        im_vectors.APLS_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_BSREM(:,1) = options.x0(:);
        im_vectors.APLS_BSREM_apu = im_vectors.APLS_BSREM(:,1);
    end
    if options.APLS && options.ROSEM_MAP
        im_vectors.APLS_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_ROSEM(:,1) = options.x0(:);
        im_vectors.APLS_ROSEM_apu = im_vectors.APLS_ROSEM(:,1);
    end
    if options.APLS && options.RBI_MAP
        im_vectors.APLS_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_RBI(:,1) = options.x0(:);
        im_vectors.APLS_RBI_apu = im_vectors.APLS_RBI(:,1);
    end
    if options.APLS && any(options.COSEM_MAP)
        im_vectors.APLS_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.APLS_COSEM(:,1) = options.x0(:);
        im_vectors.APLS_COSEM_apu = im_vectors.APLS_COSEM(:,1);
    end
    
    if options.TGV && options.OSL_OSEM
        im_vectors.TGV_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_OSL(:,1) = options.x0(:);
        im_vectors.TGV_OSL_apu = im_vectors.TGV_OSL(:,1);
    end
    if options.TGV && options.OSL_MLEM
        im_vectors.TGV_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_MLEM(:,1) = options.x0(:);
        im_vectors.TGV_MLEM_apu = im_vectors.TGV_MLEM(:,1);
    end
    if options.TGV && options.MBSREM
        im_vectors.TGV_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_MBSREM(:,1) = options.x0(:);
        im_vectors.TGV_MBSREM_apu = im_vectors.TGV_MBSREM(:,1);
    end
    if options.TGV && options.BSREM
        im_vectors.TGV_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_BSREM(:,1) = options.x0(:);
        im_vectors.TGV_BSREM_apu = im_vectors.TGV_BSREM(:,1);
    end
    if options.TGV && options.ROSEM_MAP
        im_vectors.TGV_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_ROSEM(:,1) = options.x0(:);
        im_vectors.TGV_ROSEM_apu = im_vectors.TGV_ROSEM(:,1);
    end
    if options.TGV && options.RBI_MAP
        im_vectors.TGV_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_RBI(:,1) = options.x0(:);
        im_vectors.TGV_RBI_apu = im_vectors.TGV_RBI(:,1);
    end
    if options.TGV && any(options.COSEM_MAP)
        im_vectors.TGV_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.TGV_COSEM(:,1) = options.x0(:);
        im_vectors.TGV_COSEM_apu = im_vectors.TGV_COSEM(:,1);
    end
    
    if options.NLM && options.OSL_OSEM
        im_vectors.NLM_OSL = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_OSL(:,1) = options.x0(:);
        im_vectors.NLM_OSL_apu = im_vectors.NLM_OSL(:,1);
    end
    if options.NLM && options.OSL_MLEM
        im_vectors.NLM_MLEM = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_MLEM(:,1) = options.x0(:);
        im_vectors.NLM_MLEM_apu = im_vectors.NLM_MLEM(:,1);
    end
    if options.NLM && options.MBSREM
        im_vectors.NLM_MBSREM = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_MBSREM(:,1) = options.x0(:);
        im_vectors.NLM_MBSREM_apu = im_vectors.NLM_MBSREM(:,1);
    end
    if options.NLM && options.BSREM
        im_vectors.NLM_BSREM = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_BSREM(:,1) = options.x0(:);
        im_vectors.NLM_BSREM_apu = im_vectors.NLM_BSREM(:,1);
    end
    if options.NLM && options.ROSEM_MAP
        im_vectors.NLM_ROSEM = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_ROSEM(:,1) = options.x0(:);
        im_vectors.NLM_ROSEM_apu = im_vectors.NLM_ROSEM(:,1);
    end
    if options.NLM && options.RBI_MAP
        im_vectors.NLM_RBI = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_RBI(:,1) = options.x0(:);
        im_vectors.NLM_RBI_apu = im_vectors.NLM_RBI(:,1);
    end
    if options.NLM && any(options.COSEM_MAP)
        im_vectors.NLM_COSEM = ones(N,options.Niter + 1,'single');
        im_vectors.NLM_COSEM(:,1) = options.x0(:);
        im_vectors.NLM_COSEM_apu = im_vectors.NLM_COSEM(:,1);
    end
end