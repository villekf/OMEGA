function im_vectors = reshape_vectors(im_vectors, options)
%RESHAPE_VECTORS This function simply reshapes the image vectors into
%actual images (matrices)

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
if options.save_iter
    Niter = options.Niter;
else
    Niter = 0;
end
if options.MLEM
    im_vectors.MLEM = reshape(im_vectors.MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.OSEM
    im_vectors.OSEM = reshape(im_vectors.OSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.MRAMLA
    im_vectors.MRAMLA = reshape(im_vectors.MRAMLA,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.RAMLA
    im_vectors.RAMLA = reshape(im_vectors.RAMLA,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.ROSEM
    im_vectors.ROSEM = reshape(im_vectors.ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.RBI
    im_vectors.RBI = reshape(im_vectors.RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.DRAMA
    im_vectors.DRAMA = reshape(im_vectors.DRAMA,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.COSEM
    im_vectors.COSEM = reshape(im_vectors.COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.ECOSEM
    im_vectors.ECOSEM = reshape(im_vectors.ECOSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.ACOSEM
    im_vectors.ACOSEM = reshape(im_vectors.ACOSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.MRP && options.OSL_OSEM
    im_vectors.MRP_OSL = reshape(im_vectors.MRP_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && options.OSL_MLEM
    im_vectors.MRP_OSL_MLEM = reshape(im_vectors.MRP_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && options.MBSREM
    im_vectors.MRP_MBSREM = reshape(im_vectors.MRP_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && options.BSREM
    im_vectors.MRP_BSREM = reshape(im_vectors.MRP_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && options.ROSEM_MAP
    im_vectors.MRP_ROSEM = reshape(im_vectors.MRP_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && options.OSL_RBI
    im_vectors.MRP_RBI = reshape(im_vectors.MRP_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.MRP && any(options.OSL_COSEM)
    im_vectors.MRP_COSEM = reshape(im_vectors.MRP_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.quad && options.OSL_OSEM
    im_vectors.Quad_OSL = reshape(im_vectors.Quad_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && options.OSL_MLEM
    im_vectors.Quad_OSL_MLEM = reshape(im_vectors.Quad_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && options.MBSREM
    im_vectors.Quad_MBSREM = reshape(im_vectors.Quad_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && options.BSREM
    im_vectors.Quad_BSREM = reshape(im_vectors.Quad_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && options.ROSEM_MAP
    im_vectors.Quad_ROSEM = reshape(im_vectors.Quad_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && options.OSL_RBI
    im_vectors.Quad_RBI = reshape(im_vectors.Quad_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.quad && any(options.OSL_COSEM)
    im_vectors.Quad_COSEM = reshape(im_vectors.Quad_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.Huber && options.OSL_OSEM
    im_vectors.Huber_OSL = reshape(im_vectors.Huber_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && options.OSL_MLEM
    im_vectors.Huber_OSL_MLEM = reshape(im_vectors.Huber_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && options.MBSREM
    im_vectors.Huber_MBSREM = reshape(im_vectors.Huber_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && options.BSREM
    im_vectors.Huber_BSREM = reshape(im_vectors.Huber_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && options.ROSEM_MAP
    im_vectors.Huber_ROSEM = reshape(im_vectors.Huber_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && options.OSL_RBI
    im_vectors.Huber_RBI = reshape(im_vectors.Huber_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.Huber && any(options.OSL_COSEM)
    im_vectors.Huber_COSEM = reshape(im_vectors.Huber_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.L && options.OSL_OSEM
    im_vectors.L_OSL = reshape(im_vectors.L_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && options.OSL_MLEM
    im_vectors.L_OSL_MLEM = reshape(im_vectors.L_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && options.MBSREM
    im_vectors.L_MBSREM = reshape(im_vectors.L_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && options.BSREM
    im_vectors.L_BSREM = reshape(im_vectors.L_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && options.ROSEM_MAP
    im_vectors.L_ROSEM = reshape(im_vectors.L_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && options.OSL_RBI
    im_vectors.L_RBI = reshape(im_vectors.L_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.L && any(options.OSL_COSEM)
    im_vectors.L_COSEM = reshape(im_vectors.L_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.FMH && options.OSL_OSEM
    im_vectors.FMH_OSL = reshape(im_vectors.FMH_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && options.OSL_MLEM
    im_vectors.FMH_OSL_MLEM = reshape(im_vectors.FMH_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && options.MBSREM
    im_vectors.FMH_MBSREM = reshape(im_vectors.FMH_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && options.BSREM
    im_vectors.FMH_BSREM = reshape(im_vectors.FMH_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && options.ROSEM_MAP
    im_vectors.FMH_ROSEM = reshape(im_vectors.FMH_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && options.OSL_RBI
    im_vectors.FMH_RBI = reshape(im_vectors.FMH_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.FMH && any(options.OSL_COSEM)
    im_vectors.FMH_COSEM = reshape(im_vectors.FMH_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.weighted_mean && options.OSL_OSEM
    im_vectors.Weighted_OSL = reshape(im_vectors.Weighted_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && options.OSL_MLEM
    im_vectors.Weighted_OSL_MLEM = reshape(im_vectors.Weighted_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && options.MBSREM
    im_vectors.Weighted_MBSREM = reshape(im_vectors.Weighted_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && options.BSREM
    im_vectors.Weighted_BSREM = reshape(im_vectors.Weighted_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && options.ROSEM_MAP
    im_vectors.Weighted_ROSEM = reshape(im_vectors.Weighted_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && options.OSL_RBI
    im_vectors.Weighted_RBI = reshape(im_vectors.Weighted_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.weighted_mean && any(options.OSL_COSEM)
    im_vectors.Weighted_COSEM = reshape(im_vectors.Weighted_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.TV && options.OSL_OSEM
    im_vectors.TV_OSL = reshape(im_vectors.TV_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && options.OSL_MLEM
    im_vectors.TV_OSL_MLEM = reshape(im_vectors.TV_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && options.MBSREM
    im_vectors.TV_MBSREM = reshape(im_vectors.TV_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && options.BSREM
    im_vectors.TV_BSREM = reshape(im_vectors.TV_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && options.ROSEM_MAP
    im_vectors.TV_ROSEM = reshape(im_vectors.TV_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && options.OSL_RBI
    im_vectors.TV_RBI = reshape(im_vectors.TV_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TV && any(options.OSL_COSEM)
    im_vectors.TV_COSEM = reshape(im_vectors.TV_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.AD && options.OSL_OSEM
    im_vectors.AD_OSL = reshape(im_vectors.AD_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && options.OSL_MLEM
    im_vectors.AD_OSL_MLEM = reshape(im_vectors.AD_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && options.MBSREM
    im_vectors.AD_MBSREM = reshape(im_vectors.AD_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && options.BSREM
    im_vectors.AD_BSREM = reshape(im_vectors.AD_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && options.ROSEM_MAP
    im_vectors.AD_ROSEM = reshape(im_vectors.AD_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && options.OSL_RBI
    im_vectors.AD_RBI = reshape(im_vectors.AD_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.AD && any(options.OSL_COSEM)
    im_vectors.AD_COSEM = reshape(im_vectors.AD_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.APLS && options.OSL_OSEM
    im_vectors.APLS_OSL = reshape(im_vectors.APLS_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && options.OSL_MLEM
    im_vectors.APLS_OSL_MLEM = reshape(im_vectors.APLS_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && options.MBSREM
    im_vectors.APLS_MBSREM = reshape(im_vectors.APLS_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && options.BSREM
    im_vectors.APLS_BSREM = reshape(im_vectors.APLS_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && options.ROSEM_MAP
    im_vectors.APLS_ROSEM = reshape(im_vectors.APLS_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && options.OSL_RBI
    im_vectors.APLS_RBI = reshape(im_vectors.APLS_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.APLS && any(options.OSL_COSEM)
    im_vectors.APLS_COSEM = reshape(im_vectors.APLS_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.TGV && options.OSL_OSEM
    im_vectors.TGV_OSL = reshape(im_vectors.TGV_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && options.OSL_MLEM
    im_vectors.TGV_OSL_MLEM = reshape(im_vectors.TGV_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && options.MBSREM
    im_vectors.TGV_MBSREM = reshape(im_vectors.TGV_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && options.BSREM
    im_vectors.TGV_BSREM = reshape(im_vectors.TGV_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && options.ROSEM_MAP
    im_vectors.TGV_ROSEM = reshape(im_vectors.TGV_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && options.OSL_RBI
    im_vectors.TGV_RBI = reshape(im_vectors.TGV_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.TGV && any(options.OSL_COSEM)
    im_vectors.TGV_COSEM = reshape(im_vectors.TGV_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end

if options.NLM && options.OSL_OSEM
    im_vectors.NLM_OSL = reshape(im_vectors.NLM_OSL,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && options.OSL_MLEM
    im_vectors.NLM_OSL_MLEM = reshape(im_vectors.NLM_MLEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && options.MBSREM
    im_vectors.NLM_MBSREM = reshape(im_vectors.NLM_MBSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && options.BSREM
    im_vectors.NLM_BSREM = reshape(im_vectors.NLM_BSREM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && options.ROSEM_MAP
    im_vectors.NLM_ROSEM = reshape(im_vectors.NLM_ROSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && options.OSL_RBI
    im_vectors.NLM_RBI = reshape(im_vectors.NLM_RBI,options.Nx,options.Ny,options.Nz,Niter+1);
end
if options.NLM && any(options.OSL_COSEM)
    im_vectors.NLM_COSEM = reshape(im_vectors.NLM_COSEM,options.Nx,options.Ny,options.Nz,Niter+1);
end
