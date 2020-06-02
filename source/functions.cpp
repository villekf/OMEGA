/**************************************************************************
* ArrayFire functions. Used by implementation 2.
*
* Copyright(C) 2020 Ville - Veikko Wettenhovi
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include "functions.hpp"

// Loads the input data and forms device data variables
void form_data_variables(AF_im_vectors & vec, Beta & beta, Weighting & w_vec, const mxArray *options, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t Niter, const af::array &x0, const uint32_t im_dim, const size_t koko_l, const RecMethods &MethodList, TVdata &data, 
	const uint32_t subsets, const uint32_t osa_iter0, const bool use_psf)
{
	// Load the necessary variables if the corresponding reconstruction method is used and set the initial value
	if (MethodList.MLEM) {
		vec.MLEM = af::constant(0.f, im_dim, Niter + 1);
		vec.MLEM(af::span, 0) = x0;
		
	}
	
	if (MethodList.OSEM) {
		vec.OSEM = af::constant(0.f, im_dim, Niter + 1);
		vec.OSEM(af::span, 0) = x0;
		
	}
	
	if (MethodList.MRAMLA) {
		vec.MRAMLA = af::constant(0.f, im_dim, Niter + 1);
		vec.MRAMLA(af::span, 0) = x0;
		
	}
	
	if (MethodList.RAMLA) {
		vec.RAMLA = af::constant(0.f, im_dim, Niter + 1);
		vec.RAMLA(af::span, 0) = x0;
		
	}
	
	if (MethodList.ROSEM) {
		vec.ROSEM = af::constant(0.f, im_dim, Niter + 1);
		vec.ROSEM(af::span, 0) = x0;
		
	}
	
	if (MethodList.RBI) {
		vec.RBI = af::constant(0.f, im_dim, Niter + 1);
		vec.RBI(af::span, 0) = x0;
		
	}
	
	if (MethodList.DRAMA) {
		vec.DRAMA = af::constant(0.f, im_dim, Niter + 1);
		vec.DRAMA(af::span, 0) = x0;
		
		// Relaxation parameter
		w_vec.lambda_DRAMA = (float*)mxGetData(mxGetField(options, 0, "lam_drama"));
	}
	
	if (MethodList.COSEM) {
		vec.COSEM = af::constant(0.f, im_dim, Niter + 1);
		vec.COSEM(af::span, 0) = x0;
		
		// Complete data
		vec.C_co = af::constant(0.f, im_dim, subsets);
	}
	if (MethodList.ECOSEM) {
		vec.ECOSEM = af::constant(0.f, im_dim, Niter + 1);
		vec.ECOSEM(af::span, 0) = x0;
		
		if (!MethodList.COSEM) {
			
			// Complete data
			vec.C_co = af::constant(0.f, im_dim, subsets);
		}
	}
	if (MethodList.ACOSEM) {
		vec.ACOSEM = af::constant(0.f, im_dim, Niter + 1);
		vec.ACOSEM(af::span, 0) = x0;
		
		// Complete data
		vec.C_aco = af::constant(0.f, im_dim, subsets);
	}

	if (MethodList.OSLCOSEM > 0)
		vec.C_osl = af::constant(0.f, im_dim, subsets);

	// Load the regularization parameter as well if the prior is used
	if (MethodList.MRP) {
		if (MethodList.OSLOSEM) {
			vec.MRP_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_OSEM(af::span, 0) = x0;
			
			beta.MRP_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.MRP_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_MLEM(af::span, 0) = x0;
			
			beta.MRP_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.MRP_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_MBSREM(af::span, 0) = x0;
			
			beta.MRP_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.MRP_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_BSREM(af::span, 0) = x0;
			
			beta.MRP_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.MRP_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_ROSEM(af::span, 0) = x0;
			
			beta.MRP_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.MRP_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_RBI(af::span, 0) = x0;
			
			beta.MRP_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.MRP_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.MRP_COSEM(af::span, 0) = x0;
			
			beta.MRP_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_mrp_cosem"));
		}
	}
	
	if (MethodList.Quad) {
		if (MethodList.OSLOSEM) {
			vec.Quad_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_OSEM(af::span, 0) = x0;
			
			beta.Quad_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.Quad_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_MLEM(af::span, 0) = x0;
			
			beta.Quad_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.Quad_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_MBSREM(af::span, 0) = x0;
			
			beta.Quad_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.Quad_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_BSREM(af::span, 0) = x0;
			
			beta.Quad_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.Quad_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_ROSEM(af::span, 0) = x0;
			
			beta.Quad_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.Quad_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_RBI(af::span, 0) = x0;
			
			beta.Quad_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Quad_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Quad_COSEM(af::span, 0) = x0;
			
			beta.Quad_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_quad_cosem"));
		}
	}

	if (MethodList.Huber) {
		if (MethodList.OSLOSEM) {
			vec.Huber_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_OSEM(af::span, 0) = x0;

			beta.Huber_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.Huber_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_MLEM(af::span, 0) = x0;

			beta.Huber_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.Huber_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_MBSREM(af::span, 0) = x0;

			beta.Huber_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.Huber_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_BSREM(af::span, 0) = x0;

			beta.Huber_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.Huber_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_ROSEM(af::span, 0) = x0;

			beta.Huber_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.Huber_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_RBI(af::span, 0) = x0;

			beta.Huber_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Huber_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Huber_COSEM(af::span, 0) = x0;

			beta.Huber_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_huber_cosem"));
		}
	}
	
	if (MethodList.L) {
		if (MethodList.OSLOSEM) {
			vec.L_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_OSEM(af::span, 0) = x0;
			
			beta.L_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.L_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_MLEM(af::span, 0) = x0;
			
			beta.L_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.L_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_MBSREM(af::span, 0) = x0;
			
			beta.L_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.L_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_BSREM(af::span, 0) = x0;
			
			beta.L_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.L_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_ROSEM(af::span, 0) = x0;
			
			beta.L_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.L_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.L_RBI(af::span, 0) = x0;
			
			beta.L_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_L_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.L_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.L_COSEM(af::span, 0) = x0;
			
			beta.L_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_L_cosem"));
		}
	}
	
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM) {
			vec.FMH_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_OSEM(af::span, 0) = x0;
			
			beta.FMH_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.FMH_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_MLEM(af::span, 0) = x0;
			
			beta.FMH_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.FMH_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_MBSREM(af::span, 0) = x0;
			
			beta.FMH_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.FMH_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_BSREM(af::span, 0) = x0;
			
			beta.FMH_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.FMH_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_ROSEM(af::span, 0) = x0;
			
			beta.FMH_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.FMH_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_RBI(af::span, 0) = x0;
			
			beta.FMH_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.FMH_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.FMH_COSEM(af::span, 0) = x0;
			
			beta.FMH_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_fmh_cosem"));
		}
	}
	
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM) {
			vec.Weighted_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_OSEM(af::span, 0) = x0;
			
			beta.Weighted_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.Weighted_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_MLEM(af::span, 0) = x0;
			
			beta.Weighted_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.Weighted_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_MBSREM(af::span, 0) = x0;
			
			beta.Weighted_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.Weighted_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_BSREM(af::span, 0) = x0;
			
			beta.Weighted_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.Weighted_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_ROSEM(af::span, 0) = x0;
			
			beta.Weighted_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.Weighted_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_RBI(af::span, 0) = x0;
			
			beta.Weighted_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Weighted_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.Weighted_COSEM(af::span, 0) = x0;
			
			beta.Weighted_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_weighted_cosem"));
		}
	}
	
	if (MethodList.TV) {
		if (MethodList.OSLOSEM) {
			vec.TV_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_OSEM(af::span, 0) = x0;
			
			beta.TV_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.TV_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_MLEM(af::span, 0) = x0;
			
			beta.TV_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.TV_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_MBSREM(af::span, 0) = x0;
			
			beta.TV_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.TV_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_BSREM(af::span, 0) = x0;
			
			beta.TV_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.TV_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_ROSEM(af::span, 0) = x0;
			
			beta.TV_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.TV_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_RBI(af::span, 0) = x0;
			
			beta.TV_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.TV_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TV_COSEM(af::span, 0) = x0;
			
			beta.TV_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TV_cosem"));
		}
	}
	
	if (MethodList.AD) {
		if (MethodList.OSLOSEM) {
			vec.AD_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_OSEM(af::span, 0) = x0;
			
			beta.AD_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.AD_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_MLEM(af::span, 0) = x0;
			
			beta.AD_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.AD_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_MBSREM(af::span, 0) = x0;
			
			beta.AD_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_bsrem"));
		}
		if (MethodList.BSREM) {
			vec.AD_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_BSREM(af::span, 0) = x0;
			
			beta.AD_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_mbsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.AD_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_ROSEM(af::span, 0) = x0;
			
			beta.AD_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.AD_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_RBI(af::span, 0) = x0;
			
			beta.AD_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.AD_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.AD_COSEM(af::span, 0) = x0;
			
			beta.AD_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_ad_cosem"));
		}
	}
	
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM) {
			vec.APLS_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_OSEM(af::span, 0) = x0;
			
			beta.APLS_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.APLS_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_MLEM(af::span, 0) = x0;
			
			beta.APLS_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.APLS_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_MBSREM(af::span, 0) = x0;
			
			beta.APLS_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_bsrem"));
		}
		if (MethodList.BSREM) {
			vec.APLS_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_BSREM(af::span, 0) = x0;
			
			beta.APLS_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_mbsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.APLS_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_ROSEM(af::span, 0) = x0;
			
			beta.APLS_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.APLS_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_RBI(af::span, 0) = x0;
			
			beta.APLS_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.APLS_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.APLS_COSEM(af::span, 0) = x0;
			
			beta.APLS_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_APLS_cosem"));
		}
	}

	if (MethodList.TGV) {
		if (MethodList.OSLOSEM) {
			vec.TGV_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_OSEM(af::span, 0) = x0;
			
			beta.TGV_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.TGV_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_MLEM(af::span, 0) = x0;
			
			beta.TGV_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.TGV_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_MBSREM(af::span, 0) = x0;
			
			beta.TGV_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.TGV_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_BSREM(af::span, 0) = x0;
			
			beta.TGV_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.TGV_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_ROSEM(af::span, 0) = x0;
			
			beta.TGV_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.TGV_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_RBI(af::span, 0) = x0;
			
			beta.TGV_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.TGV_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.TGV_COSEM(af::span, 0) = x0;
			
			beta.TGV_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_TGV_cosem"));
		}
	}

	if (MethodList.NLM) {
		if (MethodList.OSLOSEM) {
			vec.NLM_OSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_OSEM(af::span, 0) = x0;

			beta.NLM_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.NLM_MLEM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_MLEM(af::span, 0) = x0;

			beta.NLM_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.NLM_MBSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_MBSREM(af::span, 0) = x0;

			beta.NLM_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.NLM_BSREM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_BSREM(af::span, 0) = x0;

			beta.NLM_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.NLM_ROSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_ROSEM(af::span, 0) = x0;

			beta.NLM_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.NLM_RBI = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_RBI(af::span, 0) = x0;

			beta.NLM_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.NLM_COSEM = af::constant(0.f, im_dim, Niter + 1);
			vec.NLM_COSEM(af::span, 0) = x0;

			beta.NLM_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_NLM_cosem"));
		}
	}

	// Load TV related input data
	if (MethodList.TV && MethodList.MAP) {
		// Is anatomical reference image used
		data.TV_use_anatomical = (bool)mxGetScalar(mxGetField(options, 0, "TV_use_anatomical"));
		// Tau-value
		data.tau = (float)mxGetScalar(mxGetField(options, 0, "tau"));
		// "Smoothing" parameter, prevents zero values in the square root
		data.TVsmoothing = (float)mxGetScalar(mxGetField(options, 0, "TVsmoothing"));
		// The type of TV prior used
		data.TVtype = (int32_t)mxGetScalar(mxGetField(options, 0, "TVtype"));
		// If anatomical prior is used, load the necessary coefficients
		if (data.TV_use_anatomical) {
			mxArray* TVdata_init = mxGetField(options, 0, "TVdata");
			if (data.TVtype == 1) {
				data.s1 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s1")), afHost);
				data.s2 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s2")), afHost);
				data.s3 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s3")), afHost);
				data.s4 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s4")), afHost);
				data.s5 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s5")), afHost);
				data.s6 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s6")), afHost);
				data.s7 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s7")), afHost);
				data.s8 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s8")), afHost);
				data.s9 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s9")), afHost);
			}
			else {
				data.reference_image = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "reference_image")), afHost);
			}
			data.T = (float)mxGetScalar(mxGetField(options, 0, "T"));
			data.C = (float)mxGetScalar(mxGetField(options, 0, "C"));
		}
		// Additional weights for the TV type 3
		if (data.TVtype == 3 && !MethodList.Quad) {
			w_vec.Ndx = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndx"));
			w_vec.Ndy = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndy"));
			w_vec.Ndz = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndz"));
			w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
			w_vec.weights_TV = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_quad")), afHost);
		}
		if (data.TVtype == 3)
			data.C = (float)mxGetScalar(mxGetField(options, 0, "C"));
	}
	// General variables for neighborhood-based methods
	if ((MethodList.L || MethodList.FMH || MethodList.WeightedMean || MethodList.Quad || MethodList.Huber || MethodList.MRP || MethodList.NLM || (data.TVtype == 3 && MethodList.TV)) && MethodList.MAP) {
		// Neighborhood size
		w_vec.Ndx = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndx"));
		w_vec.Ndy = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndy"));
		w_vec.Ndz = (uint32_t)mxGetScalar(mxGetField(options, 0, "Ndz"));
		// Is normalization used in MRP, FMH, L, weighted mean or AD
		w_vec.med_no_norm = (bool)mxGetScalar(mxGetField(options, 0, "med_no_norm"));
		w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
	}
	if ((MethodList.L || MethodList.FMH || MethodList.MRP || (data.TVtype == 3 && MethodList.TV)) && MethodList.MAP) {
		// Index values for the neighborhood
		w_vec.tr_offsets = af::array(im_dim, w_vec.dimmu, (uint32_t*)mxGetData(mxGetField(options, 0, "tr_offsets")), afHost);
		if (MethodList.FMH || MethodList.Quad || MethodList.Huber)
			w_vec.inffi = (uint32_t)mxGetScalar(mxGetField(options, 0, "inffi"));
	}
	// Weights for the various priors
	if (MethodList.Quad && MethodList.MAP) {
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_quad")), afHost);	
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_quad * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_quad = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);

	}
	if (MethodList.Huber && MethodList.MAP) {
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_huber")), afHost);
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_huber * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_huber = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
		w_vec.huber_delta = (float)mxGetScalar(mxGetField(options, 0, "huber_delta"));
	}
	if (MethodList.L && MethodList.MAP)
		w_vec.a_L = af::array(w_vec.dimmu, (float*)mxGetData(mxGetField(options, 0, "a_L")), afHost);
	if (MethodList.FMH && MethodList.MAP) {
		if (Nz == 1 || w_vec.Ndz == 0)
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, (float*)mxGetData(mxGetField(options, 0, "fmh_weights")), afHost);
		else
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, (float*)mxGetData(mxGetField(options, 0, "fmh_weights")), afHost);
		w_vec.alku_fmh = static_cast<int>((uint32_t)mxGetScalar(mxGetField(options, 0, "inffi")));
	}
	if (MethodList.WeightedMean && MethodList.MAP) {
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, (float*)mxGetData(mxGetField(options, 0, "weighted_weights")), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
		// Type of mean used (arithmetic, harmonic or geometric)
		w_vec.mean_type = (int)mxGetScalar(mxGetField(options, 0, "mean_type"));
		// Sum of the weights
		w_vec.w_sum = (float)mxGetScalar(mxGetField(options, 0, "w_sum"));
	}
	if (MethodList.AD && MethodList.MAP) {
		// Time-step value
		w_vec.TimeStepAD = (float)mxGetScalar(mxGetField(options, 0, "TimeStepAD"));
		// Conductance (edge value)
		w_vec.KAD = (float)mxGetScalar(mxGetField(options, 0, "KAD"));
		// Number of AD iterations
		w_vec.NiterAD = (uint32_t)mxGetScalar(mxGetField(options, 0, "NiterAD"));
		// Flux type
		int Flux = (int)mxGetScalar(mxGetField(options, 0, "FluxType"));
		// Diffusion type
		int Diffusion = (int)mxGetScalar(mxGetField(options, 0, "DiffusionType"));
		if (Flux == 1)
			w_vec.FluxType = AF_FLUX_EXPONENTIAL;
		else if (Flux == 2)
			w_vec.FluxType = AF_FLUX_QUADRATIC;
		if (Diffusion == 1)
			w_vec.DiffusionType = AF_DIFFUSION_GRAD;
		else if (Diffusion == 2)
			w_vec.DiffusionType = AF_DIFFUSION_MCDE;
		if (!MethodList.L && !MethodList.FMH && !MethodList.WeightedMean && !MethodList.Quad && !MethodList.MRP)
			w_vec.med_no_norm = (bool)mxGetScalar(mxGetField(options, 0, "med_no_norm"));
	}
	if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.COSEM
		|| MethodList.ECOSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0)
		w_vec.MBSREM_prepass = (bool)mxGetScalar(mxGetField(options, 0, "MBSREM_prepass"));
	if (MethodList.MRAMLA || MethodList.MBSREM) {
		// Relaxation parameter
		w_vec.lambda_MBSREM = (float*)mxGetData(mxGetField(options, 0, "lam_mbsrem"));
		// Upper bound
		w_vec.U = (float)mxGetScalar(mxGetField(options, 0, "U"));
	}
	if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.COSEM
		|| MethodList.ECOSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0) {
		// Sum of the rows (measurements) of the system matrix
		//w_vec.D = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "D")), afHost);
		w_vec.D = af::constant(0.f, im_dim, 1);
		// For manual determination of the upper bound
		if (MethodList.MRAMLA || MethodList.MBSREM)
			w_vec.Amin = af::constant(0.f, koko_l, 1);
			//w_vec.Amin = af::array(koko_l, (float*)mxGetData(mxGetField(options, 0, "Amin")), afHost);
	}
	if (MethodList.APLS && MethodList.MAP) {
		// Eta value
		data.eta = (float)mxGetScalar(mxGetField(options, 0, "eta"));
		// Smoothing value
		data.APLSsmoothing = (float)mxGetScalar(mxGetField(options, 0, "APLSsmoothing"));
		// Anatomical reference image
		data.APLSReference = af::array(Nx, Ny, Nz, (float*)mxGetData(mxGetField(options, 0, "APLS_ref_image")), afHost);
		//data.APLSReference = af::moddims(data.APLSReference, Nx, Ny, Nz);
	}
	// Relaxation parameters
	if (MethodList.RAMLA || MethodList.BSREM)
		w_vec.lambda_BSREM = (float*)mxGetData(mxGetField(options, 0, "lam"));
	if (MethodList.ROSEM || MethodList.ROSEMMAP)
		w_vec.lambda_ROSEM = (float*)mxGetData(mxGetField(options, 0, "lam_rosem"));
	// Power factor for ACOSEM
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1)
		w_vec.h_ACOSEM = (float)mxGetScalar(mxGetField(options, 0, "h"));
	if (MethodList.TGV && MethodList.MAP) {
		data.TGVAlpha = (float)mxGetScalar(mxGetField(options, 0, "alphaTGV"));
		data.TGVBeta = (float)mxGetScalar(mxGetField(options, 0, "betaTGV"));
		data.NiterTGV = (uint32_t)mxGetScalar(mxGetField(options, 0, "NiterTGV"));
	}
	if (MethodList.NLM && MethodList.MAP) {
		w_vec.NLM_anatomical = (bool)mxGetScalar(mxGetField(options, 0, "NLM_use_anatomical"));
		//w_vec.NLM_gauss = (float)mxGetScalar(mxGetField(options, 0, "NLM_gauss"));
		w_vec.NLTV = (bool)mxGetScalar(mxGetField(options, 0, "NLTV"));
		w_vec.NLM_MRP = (bool)mxGetScalar(mxGetField(options, 0, "NLM_MRP"));
		if (w_vec.NLM_anatomical)
			w_vec.NLM_ref = af::array(Nx, Ny, Nz, (float*)mxGetData(mxGetField(options, 0, "NLM_ref")), afHost);
		w_vec.h2 = (float)mxGetScalar(mxGetField(options, 0, "sigma"));
		w_vec.h2 = w_vec.h2 * w_vec.h2;
		w_vec.Nlx = (uint32_t)mxGetScalar(mxGetField(options, 0, "Nlx"));
		w_vec.Nly = (uint32_t)mxGetScalar(mxGetField(options, 0, "Nly"));
		w_vec.Nlz = (uint32_t)mxGetScalar(mxGetField(options, 0, "Nlz"));
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), (float*)mxGetData(mxGetField(options, 0, "gaussianNLM")), afHost);
	}
	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			vec.custom_OSEM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "custom_osl_apu")), afHost);
			w_vec.dU_OSEM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_OSEM")), afHost);

			beta.custom_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_osem"));
		}
		if (MethodList.OSLMLEM) {
			vec.custom_MLEM = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_mlem_apu")), afHost);
			w_vec.dU_MLEM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_MLEM")), afHost);

			beta.custom_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_mlem"));
		}
		if (MethodList.MBSREM) {
			vec.custom_MBSREM = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_mbsrem_apu")), afHost);
			w_vec.dU_MBSREM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_MBSREM")), afHost);

			beta.custom_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_mbsrem"));
		}
		if (MethodList.BSREM) {
			vec.custom_BSREM = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_bsrem_apu")), afHost);
			w_vec.dU_BSREM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_BSREM")), afHost);

			beta.custom_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_bsrem"));
		}
		if (MethodList.ROSEMMAP) {
			vec.custom_ROSEM = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_rosem_apu")), afHost);
			w_vec.dU_ROSEM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_ROSEM")), afHost);

			beta.custom_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_rosem"));
		}
		if (MethodList.RBIOSL) {
			vec.custom_RBI = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_rbi_apu")), afHost);
			w_vec.dU_RBI = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_RBI")), afHost);

			beta.custom_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_rbi"));
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.custom_COSEM = af::array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "custom_cosem_apu")), afHost);
			w_vec.dU_COSEM = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "grad_COSEM")), afHost);

			beta.custom_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_cosem"));
			vec.C_osl = af::array(im_dim, subsets, (float*)mxGetData(mxGetField(options, 0, "C_osl")), afHost);
		}
		if ((MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI))
			w_vec.D = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "D")), afHost);
	}
	if (use_psf) {
		w_vec.g_dim_x = (uint32_t)mxGetScalar(mxGetField(options, 0, "g_dim_x"));
		w_vec.g_dim_y = (uint32_t)mxGetScalar(mxGetField(options, 0, "g_dim_y"));
		w_vec.g_dim_z = (uint32_t)mxGetScalar(mxGetField(options, 0, "g_dim_z"));
		w_vec.deconvolution = (bool)mxGetScalar(mxGetField(options, 0, "deblurring"));
	}
}

// Obtain the reconstruction methods used
void get_rec_methods(const mxArray * options, RecMethods &MethodList)
{
	MethodList.MLEM = (bool)mxGetScalar(mxGetField(options, 0, "mlem"));
	MethodList.OSEM = (bool)mxGetScalar(mxGetField(options, 0, "osem"));
	MethodList.RAMLA = (bool)mxGetScalar(mxGetField(options, 0, "ramla"));
	MethodList.MRAMLA = (bool)mxGetScalar(mxGetField(options, 0, "mramla"));
	MethodList.ROSEM = (bool)mxGetScalar(mxGetField(options, 0, "rosem"));
	MethodList.RBI = (bool)mxGetScalar(mxGetField(options, 0, "rbi"));
	MethodList.DRAMA = (bool)mxGetScalar(mxGetField(options, 0, "drama"));
	MethodList.COSEM = (bool)mxGetScalar(mxGetField(options, 0, "cosem"));
	MethodList.ECOSEM = (bool)mxGetScalar(mxGetField(options, 0, "ecosem"));
	MethodList.ACOSEM = (bool)mxGetScalar(mxGetField(options, 0, "acosem"));

	MethodList.MRP = (bool)mxGetScalar(mxGetField(options, 0, "MRP"));
	MethodList.Quad = (bool)mxGetScalar(mxGetField(options, 0, "quad"));
	MethodList.Huber = (bool)mxGetScalar(mxGetField(options, 0, "Huber"));
	MethodList.L = (bool)mxGetScalar(mxGetField(options, 0, "L"));
	MethodList.FMH = (bool)mxGetScalar(mxGetField(options, 0, "FMH"));
	MethodList.WeightedMean = (bool)mxGetScalar(mxGetField(options, 0, "weighted_mean"));
	MethodList.TV = (bool)mxGetScalar(mxGetField(options, 0, "TV"));
	MethodList.AD = (bool)mxGetScalar(mxGetField(options, 0, "AD"));
	MethodList.APLS = (bool)mxGetScalar(mxGetField(options, 0, "APLS"));
	MethodList.TGV = (bool)mxGetScalar(mxGetField(options, 0, "TGV"));
	MethodList.NLM = (bool)mxGetScalar(mxGetField(options, 0, "NLM"));

	MethodList.OSLMLEM = (bool)mxGetScalar(mxGetField(options, 0, "OSL_MLEM"));
	MethodList.OSLOSEM = (bool)mxGetScalar(mxGetField(options, 0, "OSL_OSEM"));
	MethodList.BSREM = (bool)mxGetScalar(mxGetField(options, 0, "BSREM"));
	MethodList.MBSREM = (bool)mxGetScalar(mxGetField(options, 0, "MBSREM"));
	MethodList.ROSEMMAP = (bool)mxGetScalar(mxGetField(options, 0, "ROSEM_MAP"));
	MethodList.RBIOSL = (bool)mxGetScalar(mxGetField(options, 0, "RBI_OSL"));
	MethodList.OSLCOSEM = (uint32_t)mxGetScalar(mxGetField(options, 0, "COSEM_OSL"));

	MethodList.MAP = (bool)mxGetScalar(mxGetField(options, 0, "MAP"));

	MethodList.CUSTOM = (bool)mxGetScalar(mxGetField(options, 0, "custom"));
}

// Create the MATLAB output
// Creates the mxArrays that are later placed in the cell array
// Creates the array pointers to the mxArrays
void create_matlab_output(matlabArrays &ArrayList, const mwSize *dimmi, const RecMethods &MethodList, const uint32_t dim_n)
{
	//mwSize dimu[2] = { 128 * 128 * 63, 9 };
	// Output arrays for the MATLAB cell array
	if (MethodList.MLEM)
		ArrayList.mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSEM)
		ArrayList.osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRAMLA)
		ArrayList.ramlaM = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.ramlaM = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.RAMLA)
		ArrayList.ramla = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.ramla = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ROSEM)
		ArrayList.rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.RBI)
		ArrayList.rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.DRAMA)
		ArrayList.drama = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.drama = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.COSEM)
		ArrayList.cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ECOSEM)
		ArrayList.ecosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.ecosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ACOSEM)
		ArrayList.acosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.acosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.MRP && MethodList.OSLOSEM)
		ArrayList.mrp_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.OSLMLEM)
		ArrayList.mrp_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.BSREM)
		ArrayList.mrp_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.MBSREM)
		ArrayList.mrp_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.ROSEMMAP)
		ArrayList.mrp_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.RBIOSL)
		ArrayList.mrp_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
		ArrayList.mrp_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.mrp_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.Quad && MethodList.OSLOSEM)
		ArrayList.quad_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.OSLMLEM)
		ArrayList.quad_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.BSREM)
		ArrayList.quad_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.MBSREM)
		ArrayList.quad_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.ROSEMMAP)
		ArrayList.quad_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.RBIOSL)
		ArrayList.quad_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
		ArrayList.quad_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.quad_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.Huber && MethodList.OSLOSEM)
		ArrayList.Huber_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.OSLMLEM)
		ArrayList.Huber_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.BSREM)
		ArrayList.Huber_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.MBSREM)
		ArrayList.Huber_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.ROSEMMAP)
		ArrayList.Huber_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.RBIOSL)
		ArrayList.Huber_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
		ArrayList.Huber_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.Huber_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.L && MethodList.OSLOSEM)
		ArrayList.L_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.OSLMLEM)
		ArrayList.L_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.BSREM)
		ArrayList.L_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.MBSREM)
		ArrayList.L_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.ROSEMMAP)
		ArrayList.L_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.RBIOSL)
		ArrayList.L_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.L && MethodList.OSLCOSEM > 0)
		ArrayList.L_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.L_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.FMH && MethodList.OSLOSEM)
		ArrayList.fmh_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.OSLMLEM)
		ArrayList.fmh_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.BSREM)
		ArrayList.fmh_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.MBSREM)
		ArrayList.fmh_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.ROSEMMAP)
		ArrayList.fmh_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.RBIOSL)
		ArrayList.fmh_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
		ArrayList.fmh_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.fmh_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.WeightedMean && MethodList.OSLOSEM)
		ArrayList.weighted_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.OSLMLEM)
		ArrayList.weighted_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.BSREM)
		ArrayList.weighted_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.MBSREM)
		ArrayList.weighted_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
		ArrayList.weighted_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.RBIOSL)
		ArrayList.weighted_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
		ArrayList.weighted_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.weighted_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.TV && MethodList.OSLOSEM)
		ArrayList.TV_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.OSLMLEM)
		ArrayList.TV_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.BSREM)
		ArrayList.TV_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.MBSREM)
		ArrayList.TV_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.ROSEMMAP)
		ArrayList.TV_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.RBIOSL)
		ArrayList.TV_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.TV && MethodList.OSLCOSEM > 0)
		ArrayList.TV_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TV_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.AD && MethodList.OSLOSEM)
		ArrayList.AD_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.OSLMLEM)
		ArrayList.AD_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.BSREM)
		ArrayList.AD_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.MBSREM)
		ArrayList.AD_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.ROSEMMAP)
		ArrayList.AD_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.RBIOSL)
		ArrayList.AD_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.AD && MethodList.OSLCOSEM > 0)
		ArrayList.AD_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.AD_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.APLS && MethodList.OSLOSEM)
		ArrayList.APLS_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.OSLMLEM)
		ArrayList.APLS_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.BSREM)
		ArrayList.APLS_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.MBSREM)
		ArrayList.APLS_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.ROSEMMAP)
		ArrayList.APLS_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.RBIOSL)
		ArrayList.APLS_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
		ArrayList.APLS_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.APLS_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.OSLOSEM)
		ArrayList.TGV_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLMLEM)
		ArrayList.TGV_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.BSREM)
		ArrayList.TGV_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MBSREM)
		ArrayList.TGV_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ROSEMMAP)
		ArrayList.TGV_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.RBIOSL)
		ArrayList.TGV_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLCOSEM > 0)
		ArrayList.TGV_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.TGV_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.OSLOSEM)
		ArrayList.NLM_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLMLEM)
		ArrayList.NLM_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.BSREM)
		ArrayList.NLM_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MBSREM)
		ArrayList.NLM_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ROSEMMAP)
		ArrayList.NLM_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.RBIOSL)
		ArrayList.NLM_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLCOSEM > 0)
		ArrayList.NLM_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
	else
		ArrayList.NLM_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			ArrayList.custom_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		}
		else
			ArrayList.custom_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.OSLMLEM)
			ArrayList.custom_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.BSREM)
			ArrayList.custom_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.MBSREM)
			ArrayList.custom_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.ROSEMMAP)
			ArrayList.custom_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.RBIOSL)
			ArrayList.custom_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		if (MethodList.OSLCOSEM > 0)
			ArrayList.custom_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
		else
			ArrayList.custom_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	}


	if (MethodList.MLEM)
		ArrayList.ele_ml = (float*)mxGetData(ArrayList.mlem);
	if (MethodList.OSEM)
		ArrayList.ele_os = (float*)mxGetData(ArrayList.osem);
	if (MethodList.MRAMLA)
		ArrayList.ele_ramlaM = (float*)mxGetData(ArrayList.ramlaM);
	if (MethodList.RAMLA)
		ArrayList.ele_ramla = (float*)mxGetData(ArrayList.ramla);
	if (MethodList.ROSEM)
		ArrayList.ele_rosem = (float*)mxGetData(ArrayList.rosem);
	if (MethodList.RBI)
		ArrayList.ele_rbi = (float*)mxGetData(ArrayList.rbi);
	if (MethodList.DRAMA)
		ArrayList.ele_drama = (float*)mxGetData(ArrayList.drama);
	if (MethodList.COSEM)
		ArrayList.ele_cosem = (float*)mxGetData(ArrayList.cosem);
	if (MethodList.ECOSEM)
		ArrayList.ele_ecosem = (float*)mxGetData(ArrayList.ecosem);
	if (MethodList.ACOSEM)
		ArrayList.ele_acosem = (float*)mxGetData(ArrayList.acosem);

	if (MethodList.OSLOSEM && MethodList.MRP)
		ArrayList.ele_mrp_osem = (float*)mxGetData(ArrayList.mrp_osem);
	if (MethodList.OSLMLEM && MethodList.MRP)
		ArrayList.ele_mrp_mlem = (float*)mxGetData(ArrayList.mrp_mlem);
	if (MethodList.MRP && MethodList.BSREM)
		ArrayList.ele_mrp_bsrem = (float*)mxGetData(ArrayList.mrp_bsrem);
	if (MethodList.MRP && MethodList.MBSREM)
		ArrayList.ele_mrp_mbsrem = (float*)mxGetData(ArrayList.mrp_mbsrem);
	if (MethodList.MRP && MethodList.ROSEMMAP)
		ArrayList.ele_mrp_rosem = (float*)mxGetData(ArrayList.mrp_rosem);
	if (MethodList.MRP && MethodList.RBIOSL)
		ArrayList.ele_mrp_rbi = (float*)mxGetData(ArrayList.mrp_rbi);
	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
		ArrayList.ele_mrp_cosem = (float*)mxGetData(ArrayList.mrp_cosem);

	if (MethodList.Quad && MethodList.OSLOSEM)
		ArrayList.ele_quad_osem = (float*)mxGetData(ArrayList.quad_osem);
	if (MethodList.OSLMLEM && MethodList.Quad)
		ArrayList.ele_quad_mlem = (float*)mxGetData(ArrayList.quad_mlem);
	if (MethodList.Quad && MethodList.BSREM)
		ArrayList.ele_quad_bsrem = (float*)mxGetData(ArrayList.quad_bsrem);
	if (MethodList.Quad && MethodList.MBSREM)
		ArrayList.ele_quad_mbsrem = (float*)mxGetData(ArrayList.quad_mbsrem);
	if (MethodList.Quad && MethodList.ROSEMMAP)
		ArrayList.ele_quad_rosem = (float*)mxGetData(ArrayList.quad_rosem);
	if (MethodList.Quad && MethodList.RBIOSL)
		ArrayList.ele_quad_rbi = (float*)mxGetData(ArrayList.quad_rbi);
	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
		ArrayList.ele_quad_cosem = (float*)mxGetData(ArrayList.quad_cosem);

	if (MethodList.Huber && MethodList.OSLOSEM)
		ArrayList.ele_Huber_osem = (float*)mxGetData(ArrayList.Huber_osem);
	if (MethodList.OSLMLEM && MethodList.Huber)
		ArrayList.ele_Huber_mlem = (float*)mxGetData(ArrayList.Huber_mlem);
	if (MethodList.Huber && MethodList.BSREM)
		ArrayList.ele_Huber_bsrem = (float*)mxGetData(ArrayList.Huber_bsrem);
	if (MethodList.Huber && MethodList.MBSREM)
		ArrayList.ele_Huber_mbsrem = (float*)mxGetData(ArrayList.Huber_mbsrem);
	if (MethodList.Huber && MethodList.ROSEMMAP)
		ArrayList.ele_Huber_rosem = (float*)mxGetData(ArrayList.Huber_rosem);
	if (MethodList.Huber && MethodList.RBIOSL)
		ArrayList.ele_Huber_rbi = (float*)mxGetData(ArrayList.Huber_rbi);
	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
		ArrayList.ele_Huber_cosem = (float*)mxGetData(ArrayList.Huber_cosem);

	if (MethodList.L && MethodList.OSLOSEM)
		ArrayList.ele_L_osem = (float*)mxGetData(ArrayList.L_osem);
	if (MethodList.OSLMLEM && MethodList.L)
		ArrayList.ele_L_mlem = (float*)mxGetData(ArrayList.L_mlem);
	if (MethodList.L && MethodList.BSREM)
		ArrayList.ele_L_bsrem = (float*)mxGetData(ArrayList.L_bsrem);
	if (MethodList.L && MethodList.MBSREM)
		ArrayList.ele_L_mbsrem = (float*)mxGetData(ArrayList.L_mbsrem);
	if (MethodList.L && MethodList.ROSEMMAP)
		ArrayList.ele_L_rosem = (float*)mxGetData(ArrayList.L_rosem);
	if (MethodList.L && MethodList.RBIOSL)
		ArrayList.ele_L_rbi = (float*)mxGetData(ArrayList.L_rbi);
	if (MethodList.L && MethodList.OSLCOSEM > 0)
		ArrayList.ele_L_cosem = (float*)mxGetData(ArrayList.L_cosem);

	if (MethodList.FMH && MethodList.OSLOSEM)
		ArrayList.ele_fmh_osem = (float*)mxGetData(ArrayList.fmh_osem);
	if (MethodList.OSLMLEM && MethodList.FMH)
		ArrayList.ele_fmh_mlem = (float*)mxGetData(ArrayList.fmh_mlem);
	if (MethodList.FMH && MethodList.BSREM)
		ArrayList.ele_fmh_bsrem = (float*)mxGetData(ArrayList.fmh_bsrem);
	if (MethodList.FMH && MethodList.MBSREM)
		ArrayList.ele_fmh_mbsrem = (float*)mxGetData(ArrayList.fmh_mbsrem);
	if (MethodList.FMH && MethodList.ROSEMMAP)
		ArrayList.ele_fmh_rosem = (float*)mxGetData(ArrayList.fmh_rosem);
	if (MethodList.FMH && MethodList.RBIOSL)
		ArrayList.ele_fmh_rbi = (float*)mxGetData(ArrayList.fmh_rbi);
	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
		ArrayList.ele_fmh_cosem = (float*)mxGetData(ArrayList.fmh_cosem);

	if (MethodList.WeightedMean && MethodList.OSLOSEM)
		ArrayList.ele_weighted_osem = (float*)mxGetData(ArrayList.weighted_osem);
	if (MethodList.OSLMLEM && MethodList.WeightedMean)
		ArrayList.ele_weighted_mlem = (float*)mxGetData(ArrayList.weighted_mlem);
	if (MethodList.WeightedMean && MethodList.BSREM)
		ArrayList.ele_weighted_bsrem = (float*)mxGetData(ArrayList.weighted_bsrem);
	if (MethodList.WeightedMean && MethodList.MBSREM)
		ArrayList.ele_weighted_mbsrem = (float*)mxGetData(ArrayList.weighted_mbsrem);
	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
		ArrayList.ele_weighted_rosem = (float*)mxGetData(ArrayList.weighted_rosem);
	if (MethodList.WeightedMean && MethodList.RBIOSL)
		ArrayList.ele_weighted_rbi = (float*)mxGetData(ArrayList.weighted_rbi);
	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
		ArrayList.ele_weighted_cosem = (float*)mxGetData(ArrayList.weighted_cosem);

	if (MethodList.TV && MethodList.OSLOSEM)
		ArrayList.ele_TV_osem = (float*)mxGetData(ArrayList.TV_osem);
	if (MethodList.OSLMLEM && MethodList.TV)
		ArrayList.ele_TV_mlem = (float*)mxGetData(ArrayList.TV_mlem);
	if (MethodList.TV && MethodList.BSREM)
		ArrayList.ele_TV_bsrem = (float*)mxGetData(ArrayList.TV_bsrem);
	if (MethodList.TV && MethodList.MBSREM)
		ArrayList.ele_TV_mbsrem = (float*)mxGetData(ArrayList.TV_mbsrem);
	if (MethodList.TV && MethodList.ROSEMMAP)
		ArrayList.ele_TV_rosem = (float*)mxGetData(ArrayList.TV_rosem);
	if (MethodList.TV && MethodList.RBIOSL)
		ArrayList.ele_TV_rbi = (float*)mxGetData(ArrayList.TV_rbi);
	if (MethodList.TV && MethodList.OSLCOSEM > 0)
		ArrayList.ele_TV_cosem = (float*)mxGetData(ArrayList.TV_cosem);

	if (MethodList.AD && MethodList.OSLOSEM)
		ArrayList.ele_AD_osem = (float*)mxGetData(ArrayList.AD_osem);
	if (MethodList.OSLMLEM && MethodList.AD)
		ArrayList.ele_AD_mlem = (float*)mxGetData(ArrayList.AD_mlem);
	if (MethodList.AD && MethodList.BSREM)
		ArrayList.ele_AD_bsrem = (float*)mxGetData(ArrayList.AD_bsrem);
	if (MethodList.AD && MethodList.MBSREM)
		ArrayList.ele_AD_mbsrem = (float*)mxGetData(ArrayList.AD_mbsrem);
	if (MethodList.AD && MethodList.ROSEMMAP)
		ArrayList.ele_AD_rosem = (float*)mxGetData(ArrayList.AD_rosem);
	if (MethodList.AD && MethodList.RBIOSL)
		ArrayList.ele_AD_rbi = (float*)mxGetData(ArrayList.AD_rbi);
	if (MethodList.AD && MethodList.OSLCOSEM > 0)
		ArrayList.ele_AD_cosem = (float*)mxGetData(ArrayList.AD_cosem);

	if (MethodList.APLS && MethodList.OSLOSEM)
		ArrayList.ele_APLS_osem = (float*)mxGetData(ArrayList.APLS_osem);
	if (MethodList.OSLMLEM && MethodList.APLS)
		ArrayList.ele_APLS_mlem = (float*)mxGetData(ArrayList.APLS_mlem);
	if (MethodList.APLS && MethodList.BSREM)
		ArrayList.ele_APLS_bsrem = (float*)mxGetData(ArrayList.APLS_bsrem);
	if (MethodList.APLS && MethodList.MBSREM)
		ArrayList.ele_APLS_mbsrem = (float*)mxGetData(ArrayList.APLS_mbsrem);
	if (MethodList.APLS && MethodList.ROSEMMAP)
		ArrayList.ele_APLS_rosem = (float*)mxGetData(ArrayList.APLS_rosem);
	if (MethodList.APLS && MethodList.RBIOSL)
		ArrayList.ele_APLS_rbi = (float*)mxGetData(ArrayList.APLS_rbi);
	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
		ArrayList.ele_APLS_cosem = (float*)mxGetData(ArrayList.APLS_cosem);

	if (MethodList.OSLOSEM)
		ArrayList.ele_TGV_osem = (float*)mxGetData(ArrayList.TGV_osem);
	if (MethodList.OSLMLEM && MethodList.TGV)
		ArrayList.ele_TGV_mlem = (float*)mxGetData(ArrayList.TGV_mlem);
	if (MethodList.BSREM)
		ArrayList.ele_TGV_bsrem = (float*)mxGetData(ArrayList.TGV_bsrem);
	if (MethodList.MBSREM)
		ArrayList.ele_TGV_mbsrem = (float*)mxGetData(ArrayList.TGV_mbsrem);
	if (MethodList.ROSEMMAP)
		ArrayList.ele_TGV_rosem = (float*)mxGetData(ArrayList.TGV_rosem);
	if (MethodList.RBIOSL)
		ArrayList.ele_TGV_rbi = (float*)mxGetData(ArrayList.TGV_rbi);
	if (MethodList.OSLCOSEM > 0)
		ArrayList.ele_TGV_cosem = (float*)mxGetData(ArrayList.TGV_cosem);

	if (MethodList.OSLOSEM)
		ArrayList.ele_NLM_osem = (float*)mxGetData(ArrayList.NLM_osem);
	if (MethodList.OSLMLEM && MethodList.NLM)
		ArrayList.ele_NLM_mlem = (float*)mxGetData(ArrayList.NLM_mlem);
	if (MethodList.BSREM)
		ArrayList.ele_NLM_bsrem = (float*)mxGetData(ArrayList.NLM_bsrem);
	if (MethodList.MBSREM)
		ArrayList.ele_NLM_mbsrem = (float*)mxGetData(ArrayList.NLM_mbsrem);
	if (MethodList.ROSEMMAP)
		ArrayList.ele_NLM_rosem = (float*)mxGetData(ArrayList.NLM_rosem);
	if (MethodList.RBIOSL)
		ArrayList.ele_NLM_rbi = (float*)mxGetData(ArrayList.NLM_rbi);
	if (MethodList.OSLCOSEM > 0)
		ArrayList.ele_NLM_cosem = (float*)mxGetData(ArrayList.NLM_cosem);

	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			ArrayList.ele_custom_osem = (float*)mxGetData(ArrayList.custom_osem);
		}
		if (MethodList.OSLMLEM)
			ArrayList.ele_custom_mlem = (float*)mxGetData(ArrayList.custom_mlem);
		if (MethodList.BSREM)
			ArrayList.ele_custom_bsrem = (float*)mxGetData(ArrayList.custom_bsrem);
		if (MethodList.MBSREM)
			ArrayList.ele_custom_mbsrem = (float*)mxGetData(ArrayList.custom_mbsrem);
		if (MethodList.ROSEMMAP)
			ArrayList.ele_custom_rosem = (float*)mxGetData(ArrayList.custom_rosem);
		if (MethodList.RBIOSL)
			ArrayList.ele_custom_rbi = (float*)mxGetData(ArrayList.custom_rbi);
		if (MethodList.OSLCOSEM > 0) {
			ArrayList.ele_custom_cosem = (float*)mxGetData(ArrayList.custom_cosem);
			ArrayList.ele_c_osl_custom = (float*)mxGetData(ArrayList.c_osl_custom);
		}
	}


}

// Transfers the device data to host
// First transfer the ArrayFire arrays from the device to the host pointers pointing to the mxArrays
// Transfer the mxArrays to the cell
void device_to_host_cell(matlabArrays &ArrayList, const RecMethods &MethodList, AF_im_vectors & vec, uint32_t & oo, mxArray * cell, Weighting & w_vec)
{
	// Transfer data back to host
	if (MethodList.OSEM) {
		vec.OSEM.host(ArrayList.ele_os);
	}
	if (MethodList.MLEM)
		vec.MLEM.host(ArrayList.ele_ml);
	if (MethodList.MRAMLA)
		vec.MRAMLA.host(ArrayList.ele_ramlaM);
	if (MethodList.RAMLA)
		vec.RAMLA.host(ArrayList.ele_ramla);
	if (MethodList.ROSEM)
		vec.ROSEM.host(ArrayList.ele_rosem);
	if (MethodList.RBI)
		vec.RBI.host(ArrayList.ele_rbi);
	if (MethodList.DRAMA)
		vec.DRAMA.host(ArrayList.ele_drama);
	if (MethodList.COSEM) {
		vec.COSEM.host(ArrayList.ele_cosem);
	}
	if (MethodList.ECOSEM)
		vec.ECOSEM.host(ArrayList.ele_ecosem);
	if (MethodList.ACOSEM)
		vec.ACOSEM.host(ArrayList.ele_acosem);

	if (MethodList.MRP && MethodList.OSLOSEM)
		vec.MRP_OSEM.host(ArrayList.ele_mrp_osem);
	if (MethodList.MRP && MethodList.OSLMLEM)
		vec.MRP_MLEM.host(ArrayList.ele_mrp_mlem);
	if (MethodList.MRP && MethodList.BSREM)
		vec.MRP_BSREM.host(ArrayList.ele_mrp_bsrem);
	if (MethodList.MRP && MethodList.MBSREM)
		vec.MRP_MBSREM.host(ArrayList.ele_mrp_mbsrem);
	if (MethodList.MRP && MethodList.ROSEMMAP)
		vec.MRP_ROSEM.host(ArrayList.ele_mrp_rosem);
	if (MethodList.MRP && MethodList.RBIOSL)
		vec.MRP_RBI.host(ArrayList.ele_mrp_rbi);
	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
		vec.MRP_COSEM.host(ArrayList.ele_mrp_cosem);

	if (MethodList.Quad && MethodList.OSLOSEM)
		vec.Quad_OSEM.host(ArrayList.ele_quad_osem);
	if (MethodList.Quad && MethodList.OSLMLEM)
		vec.Quad_MLEM.host(ArrayList.ele_quad_mlem);
	if (MethodList.Quad && MethodList.BSREM)
		vec.Quad_BSREM.host(ArrayList.ele_quad_bsrem);
	if (MethodList.Quad && MethodList.MBSREM)
		vec.Quad_MBSREM.host(ArrayList.ele_quad_mbsrem);
	if (MethodList.Quad && MethodList.ROSEMMAP)
		vec.Quad_ROSEM.host(ArrayList.ele_quad_rosem);
	if (MethodList.Quad && MethodList.RBIOSL)
		vec.Quad_RBI.host(ArrayList.ele_quad_rbi);
	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
		vec.Quad_COSEM.host(ArrayList.ele_quad_cosem);

	if (MethodList.Huber && MethodList.OSLOSEM)
		vec.Huber_OSEM.host(ArrayList.ele_Huber_osem);
	if (MethodList.Huber && MethodList.OSLMLEM)
		vec.Huber_MLEM.host(ArrayList.ele_Huber_mlem);
	if (MethodList.Huber && MethodList.BSREM)
		vec.Huber_BSREM.host(ArrayList.ele_Huber_bsrem);
	if (MethodList.Huber && MethodList.MBSREM)
		vec.Huber_MBSREM.host(ArrayList.ele_Huber_mbsrem);
	if (MethodList.Huber && MethodList.ROSEMMAP)
		vec.Huber_ROSEM.host(ArrayList.ele_Huber_rosem);
	if (MethodList.Huber && MethodList.RBIOSL)
		vec.Huber_RBI.host(ArrayList.ele_Huber_rbi);
	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
		vec.Huber_COSEM.host(ArrayList.ele_Huber_cosem);

	if (MethodList.L && MethodList.OSLOSEM)
		vec.L_OSEM.host(ArrayList.ele_L_osem);
	if (MethodList.L && MethodList.OSLMLEM)
		vec.L_MLEM.host(ArrayList.ele_L_mlem);
	if (MethodList.L && MethodList.BSREM)
		vec.L_BSREM.host(ArrayList.ele_L_bsrem);
	if (MethodList.L && MethodList.MBSREM)
		vec.L_MBSREM.host(ArrayList.ele_L_mbsrem);
	if (MethodList.L && MethodList.ROSEMMAP)
		vec.L_ROSEM.host(ArrayList.ele_L_rosem);
	if (MethodList.L && MethodList.RBIOSL)
		vec.L_RBI.host(ArrayList.ele_L_rbi);
	if (MethodList.L && MethodList.OSLCOSEM > 0)
		vec.L_COSEM.host(ArrayList.ele_L_cosem);

	if (MethodList.FMH && MethodList.OSLOSEM)
		vec.FMH_OSEM.host(ArrayList.ele_fmh_osem);
	if (MethodList.FMH && MethodList.OSLMLEM)
		vec.FMH_MLEM.host(ArrayList.ele_fmh_mlem);
	if (MethodList.FMH && MethodList.BSREM)
		vec.FMH_BSREM.host(ArrayList.ele_fmh_bsrem);
	if (MethodList.FMH && MethodList.MBSREM)
		vec.FMH_MBSREM.host(ArrayList.ele_fmh_mbsrem);
	if (MethodList.FMH && MethodList.ROSEMMAP)
		vec.FMH_ROSEM.host(ArrayList.ele_fmh_rosem);
	if (MethodList.FMH && MethodList.RBIOSL)
		vec.FMH_RBI.host(ArrayList.ele_fmh_rbi);
	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
		vec.FMH_COSEM.host(ArrayList.ele_fmh_cosem);

	if (MethodList.WeightedMean && MethodList.OSLOSEM)
		vec.Weighted_OSEM.host(ArrayList.ele_weighted_osem);
	if (MethodList.WeightedMean && MethodList.OSLMLEM)
		vec.Weighted_MLEM.host(ArrayList.ele_weighted_mlem);
	if (MethodList.WeightedMean && MethodList.BSREM)
		vec.Weighted_BSREM.host(ArrayList.ele_weighted_bsrem);
	if (MethodList.WeightedMean && MethodList.MBSREM)
		vec.Weighted_MBSREM.host(ArrayList.ele_weighted_mbsrem);
	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
		vec.Weighted_ROSEM.host(ArrayList.ele_weighted_rosem);
	if (MethodList.WeightedMean && MethodList.RBIOSL)
		vec.Weighted_RBI.host(ArrayList.ele_weighted_rbi);
	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
		vec.Weighted_COSEM.host(ArrayList.ele_weighted_cosem);
	
	if (MethodList.TV && MethodList.OSLOSEM)
		vec.TV_OSEM.host(ArrayList.ele_TV_osem);
	if (MethodList.TV && MethodList.OSLMLEM)
		vec.TV_MLEM.host(ArrayList.ele_TV_mlem);
	if (MethodList.TV && MethodList.BSREM)
		vec.TV_BSREM.host(ArrayList.ele_TV_bsrem);
	if (MethodList.TV && MethodList.MBSREM)
		vec.TV_MBSREM.host(ArrayList.ele_TV_mbsrem);
	if (MethodList.TV && MethodList.ROSEMMAP)
		vec.TV_ROSEM.host(ArrayList.ele_TV_rosem);
	if (MethodList.TV && MethodList.RBIOSL)
		vec.TV_RBI.host(ArrayList.ele_TV_rbi);
	if (MethodList.TV && MethodList.OSLCOSEM > 0)
		vec.TV_COSEM.host(ArrayList.ele_TV_cosem);

	if (MethodList.AD && MethodList.OSLOSEM)
		vec.AD_OSEM.host(ArrayList.ele_AD_osem);
	if (MethodList.AD && MethodList.OSLMLEM)
		vec.AD_MLEM.host(ArrayList.ele_AD_mlem);
	if (MethodList.AD && MethodList.BSREM)
		vec.AD_BSREM.host(ArrayList.ele_AD_bsrem);
	if (MethodList.AD && MethodList.MBSREM)
		vec.AD_MBSREM.host(ArrayList.ele_AD_mbsrem);
	if (MethodList.AD && MethodList.ROSEMMAP)
		vec.AD_ROSEM.host(ArrayList.ele_AD_rosem);
	if (MethodList.AD && MethodList.RBIOSL)
		vec.AD_RBI.host(ArrayList.ele_AD_rbi);
	if (MethodList.AD && MethodList.OSLCOSEM > 0)
		vec.AD_COSEM.host(ArrayList.ele_AD_cosem);

	if (MethodList.APLS && MethodList.OSLOSEM)
		vec.APLS_OSEM.host(ArrayList.ele_APLS_osem);
	if (MethodList.APLS && MethodList.OSLMLEM)
		vec.APLS_MLEM.host(ArrayList.ele_APLS_mlem);
	if (MethodList.APLS && MethodList.BSREM)
		vec.APLS_BSREM.host(ArrayList.ele_APLS_bsrem);
	if (MethodList.APLS && MethodList.MBSREM)
		vec.APLS_MBSREM.host(ArrayList.ele_APLS_mbsrem);
	if (MethodList.APLS && MethodList.ROSEMMAP)
		vec.APLS_ROSEM.host(ArrayList.ele_APLS_rosem);
	if (MethodList.APLS && MethodList.RBIOSL)
		vec.APLS_RBI.host(ArrayList.ele_APLS_rbi);
	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
		vec.APLS_COSEM.host(ArrayList.ele_APLS_cosem);

	if (MethodList.TGV && MethodList.OSLOSEM)
		vec.TGV_OSEM.host(ArrayList.ele_TGV_osem);
	if (MethodList.TGV && MethodList.OSLMLEM)
		vec.TGV_MLEM.host(ArrayList.ele_TGV_mlem);
	if (MethodList.TGV && MethodList.BSREM)
		vec.TGV_BSREM.host(ArrayList.ele_TGV_bsrem);
	if (MethodList.TGV && MethodList.MBSREM)
		vec.TGV_MBSREM.host(ArrayList.ele_TGV_mbsrem);
	if (MethodList.TGV && MethodList.ROSEMMAP)
		vec.TGV_ROSEM.host(ArrayList.ele_TGV_rosem);
	if (MethodList.TGV && MethodList.RBIOSL)
		vec.TGV_RBI.host(ArrayList.ele_TGV_rbi);
	if (MethodList.TGV && MethodList.OSLCOSEM > 0)
		vec.TGV_COSEM.host(ArrayList.ele_TGV_cosem);

	if (MethodList.NLM && MethodList.OSLOSEM)
		vec.NLM_OSEM.host(ArrayList.ele_NLM_osem);
	if (MethodList.NLM && MethodList.OSLMLEM)
		vec.NLM_MLEM.host(ArrayList.ele_NLM_mlem);
	if (MethodList.NLM && MethodList.BSREM)
		vec.NLM_BSREM.host(ArrayList.ele_NLM_bsrem);
	if (MethodList.NLM && MethodList.MBSREM)
		vec.NLM_MBSREM.host(ArrayList.ele_NLM_mbsrem);
	if (MethodList.NLM && MethodList.ROSEMMAP)
		vec.NLM_ROSEM.host(ArrayList.ele_NLM_rosem);
	if (MethodList.NLM && MethodList.RBIOSL)
		vec.NLM_RBI.host(ArrayList.ele_NLM_rbi);
	if (MethodList.NLM && MethodList.OSLCOSEM > 0)
		vec.NLM_COSEM.host(ArrayList.ele_NLM_cosem);

	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM)
			vec.custom_OSEM.host(ArrayList.ele_custom_osem);
		if (MethodList.OSLMLEM)
			vec.custom_MLEM.host(ArrayList.ele_custom_mlem);
		if (MethodList.BSREM)
			vec.custom_BSREM.host(ArrayList.ele_custom_bsrem);
		if (MethodList.MBSREM)
			vec.custom_MBSREM.host(ArrayList.ele_custom_mbsrem);
		if (MethodList.ROSEMMAP)
			vec.custom_ROSEM.host(ArrayList.ele_custom_rosem);
		if (MethodList.RBIOSL)
			vec.custom_RBI.host(ArrayList.ele_custom_rbi);
		if (MethodList.OSLCOSEM > 0) {
			vec.custom_COSEM.host(ArrayList.ele_custom_cosem);
			vec.C_osl.host(ArrayList.ele_c_osl_custom);
		}
		if ((MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL))
			w_vec.D.host(ArrayList.ele_D_custom);
	}

	af::sync();

	if (MethodList.MLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mlem));
	oo++;
	if (MethodList.OSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.osem));
	oo++;
	if (MethodList.MRAMLA)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.ramlaM));
	oo++;
	if (MethodList.RAMLA)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.ramla));
	oo++;
	if (MethodList.ROSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.rosem));
	oo++;
	if (MethodList.RBI)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.rbi));
	oo++;
	if (MethodList.DRAMA)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.drama));
	oo++;
	if (MethodList.COSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.cosem));
	oo++;
	if (MethodList.ECOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.ecosem));
	oo++;
	if (MethodList.ACOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.acosem));
	oo++;

	if (MethodList.MRP && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_osem));
	oo++;
	if (MethodList.MRP && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_mlem));
	oo++;
	if (MethodList.MRP && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_bsrem));
	oo++;
	if (MethodList.MRP && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_mbsrem));
	oo++;
	if (MethodList.MRP && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_rosem));
	oo++;
	if (MethodList.MRP && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_rbi));
	oo++;
	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.mrp_cosem));
	oo++;

	if (MethodList.Quad && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_osem));
	oo++;
	if (MethodList.Quad && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_mlem));
	oo++;
	if (MethodList.Quad && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_bsrem));
	oo++;
	if (MethodList.Quad && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_mbsrem));
	oo++;
	if (MethodList.Quad && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_rosem));
	oo++;
	if (MethodList.Quad && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_rbi));
	oo++;
	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.quad_cosem));
	oo++;

	if (MethodList.Huber && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_osem));
	oo++;
	if (MethodList.Huber && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_mlem));
	oo++;
	if (MethodList.Huber && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_bsrem));
	oo++;
	if (MethodList.Huber && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_mbsrem));
	oo++;
	if (MethodList.Huber && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_rosem));
	oo++;
	if (MethodList.Huber && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_rbi));
	oo++;
	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.Huber_cosem));
	oo++;

	if (MethodList.L && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_osem));
	oo++;
	if (MethodList.L && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_mlem));
	oo++;
	if (MethodList.L && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_bsrem));
	oo++;
	if (MethodList.L && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_mbsrem));
	oo++;
	if (MethodList.L && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_rosem));
	oo++;
	if (MethodList.L && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_rbi));
	oo++;
	if (MethodList.L && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.L_cosem));
	oo++;

	if (MethodList.FMH && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_osem));
	oo++;
	if (MethodList.FMH && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_mlem));
	oo++;
	if (MethodList.FMH && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_bsrem));
	oo++;
	if (MethodList.FMH && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_mbsrem));
	oo++;
	if (MethodList.FMH && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_rosem));
	oo++;
	if (MethodList.FMH && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_rbi));
	oo++;
	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.fmh_cosem));
	oo++;

	if (MethodList.WeightedMean && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_osem));
	oo++;
	if (MethodList.WeightedMean && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_mlem));
	oo++;
	if (MethodList.WeightedMean && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_bsrem));
	oo++;
	if (MethodList.WeightedMean && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_mbsrem));
	oo++;
	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_rosem));
	oo++;
	if (MethodList.WeightedMean && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_rbi));
	oo++;
	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.weighted_cosem));
	oo++;

	if (MethodList.TV && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_osem));
	oo++;
	if (MethodList.TV && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_mlem));
	oo++;
	if (MethodList.TV && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_bsrem));
	oo++;
	if (MethodList.TV && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_mbsrem));
	oo++;
	if (MethodList.TV && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_rosem));
	oo++;
	if (MethodList.TV && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_rbi));
	oo++;
	if (MethodList.TV && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TV_cosem));
	oo++;

	if (MethodList.AD && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_osem));
	oo++;
	if (MethodList.AD && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_mlem));
	oo++;
	if (MethodList.AD && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_bsrem));
	oo++;
	if (MethodList.AD && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_mbsrem));
	oo++;
	if (MethodList.AD && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_rosem));
	oo++;
	if (MethodList.AD && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_rbi));
	oo++;
	if (MethodList.AD && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.AD_cosem));
	oo++;

	if (MethodList.APLS && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_osem));
	oo++;
	if (MethodList.APLS && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_mlem));
	oo++;
	if (MethodList.APLS && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_bsrem));
	oo++;
	if (MethodList.APLS && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_mbsrem));
	oo++;
	if (MethodList.APLS && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_rosem));
	oo++;
	if (MethodList.APLS && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_rbi));
	oo++;
	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.APLS_cosem));
	oo++;

	if (MethodList.TGV && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_osem));
	oo++;
	if (MethodList.TGV && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_mlem));
	oo++;
	if (MethodList.TGV && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_bsrem));
	oo++;
	if (MethodList.TGV && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_mbsrem));
	oo++;
	if (MethodList.TGV && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_rosem));
	oo++;
	if (MethodList.TGV && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_rbi));
	oo++;
	if (MethodList.TGV && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.TGV_cosem));
	oo++;

	if (MethodList.NLM && MethodList.OSLOSEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_osem));
	oo++;
	if (MethodList.NLM && MethodList.OSLMLEM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_mlem));
	oo++;
	if (MethodList.NLM && MethodList.BSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_bsrem));
	oo++;
	if (MethodList.NLM && MethodList.MBSREM)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_mbsrem));
	oo++;
	if (MethodList.NLM && MethodList.ROSEMMAP)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_rosem));
	oo++;
	if (MethodList.NLM && MethodList.RBIOSL)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_rbi));
	oo++;
	if (MethodList.NLM && MethodList.OSLCOSEM > 0)
		mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.NLM_cosem));
	oo++;


	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_osem));
		oo++;
		if (MethodList.OSLMLEM)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_mlem));
		oo++;
		if (MethodList.BSREM)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_bsrem));
		oo++;
		if (MethodList.MBSREM)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_mbsrem));
		oo++;
		if (MethodList.ROSEMMAP)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_rosem));
		oo++;
		if (MethodList.RBIOSL)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_rbi));
		oo++;
		if (MethodList.OSLCOSEM > 0)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.custom_cosem));
		oo++;
		if (MethodList.OSLCOSEM > 0)
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.c_osl_custom));
		oo++;
		if ((MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL))
			mxSetCell(cell, static_cast<mwIndex>(oo), mxDuplicateArray(ArrayList.D_custom));
		oo++;
	}
	else
		oo += 9u;
		//oo += 16u;
}

// Compute the epsilon value for MBSREM and MRAMLA
float MBSREM_epsilon(const af::array & Sino, const float epps, const uint32_t randoms_correction, const af::array& rand, const af::array& D)
{
	af::array P_Sino, apu;
	if (randoms_correction == 1u) {
		P_Sino = Sino(Sino > 0.f & rand == 0.f);
		apu = D + rand;
		apu = sum(Sino * log(apu) - apu);
	}
	else {
		P_Sino = Sino(Sino > 0.f);
		apu = af::sum(D * log(D) - D);
	}
	af::array hk_summa = Sino * af::log(Sino) - Sino;
	hk_summa(af::isNaN(hk_summa)) = 0.f;
	if (randoms_correction == 1u)
		hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Sino > 0.f & rand == 0.f), batchMinus);
	else
		hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Sino > 0.f), batchMinus);
	af::array epsilon = (af::min)(P_Sino, af::exp(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
	float eps;
	eps = af::min<float>(epsilon);
	eps = eps <= 0.f ? epps : eps;
	return eps;
}

// Various batch functions
af::array batchMinus(const af::array &lhs, const af::array &rhs) {
	return lhs - rhs;
}

af::array batchPlus(const af::array &lhs, const af::array &rhs) {
	return lhs + rhs;
}

af::array batchMul(const af::array &lhs, const af::array &rhs) {
	return lhs * rhs;
}

af::array batchDiv(const af::array &lhs, const af::array &rhs) {
	return lhs / rhs;
}

af::array batchNotEqual(const af::array &lhs, const af::array &rhs) {
	return lhs != rhs;
}

af::array padding(const af::array& im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz,
	const bool zero_pad)
{
	af::array padd;
	if (zero_pad == 1) {
		af::dtype type = im.type();
		if (Nz == 1) {
			if (im.dims(1) == 1)
				padd = moddims(im, Nx, Ny, Nz);
			else
				padd = im;
			padd = moddims(im, Nx, Ny, Nz);
			af::array out = af::constant(0, padd.dims(0) + 2 * Ndx, padd.dims(1) + 2 * Ndy, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(padd.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(padd.dims(1)))) = padd;
			padd = out;
		}
		else {
			if (im.dims(2) == 1)
				padd = moddims(im, Nx, Ny, Nz);
			else
				padd = im;
			padd = moddims(im, Nx, Ny, Nz);
			af::array out = af::constant(0, padd.dims(0) + 2 * Ndx, padd.dims(1) + 2 * Ndy, padd.dims(2) + 2 * Ndz, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(padd.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(padd.dims(1))),
				static_cast<double>(Ndz) + af::seq(static_cast<double>(padd.dims(2)))) = padd;
			padd = out;
		}
	}
	else {
		if (im.dims(1) == 1)
			padd = moddims(im, Nx, Ny, Nz);
		else
			padd = im;
		af::array out = af::join(0, af::flip(padd(af::seq(static_cast<double>(Ndx)), af::span, af::span), 0), padd, af::flip(padd(af::seq(static_cast<double>(padd.dims(0) - Ndx), static_cast<double>(padd.dims(0) - 1)), af::span, af::span), 0));
		out = af::join(1, af::flip(out(af::span, af::seq(static_cast<double>(Ndy)), af::span), 1), out, af::flip(out(af::span, af::seq(static_cast<double>(out.dims(1) - Ndy), static_cast<double>(out.dims(1) - 1)), af::span), 1));
		if (Nz == 1 || Ndz == 0) {
			//af::array out = af::constant(0.f, padd.dims(0) + 2 * Ndx, padd.dims(1) + 2 * Ndy);
			//out(static_cast<double>(Ndx) + af::seq(static_cast<double>(padd.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(padd.dims(1)))) = padd;
			//out(af::seq(static_cast<double>(Ndx - 1))
			//padd = out;
		}
		else {
			out = af::join(2, af::flip(out(af::span, af::span, af::seq(static_cast<double>(Ndz))), 2), out, af::flip(out(af::span, af::span, af::seq(static_cast<double>(out.dims(2) - Ndz), static_cast<double>(out.dims(2) - 1))), 2));
			//padd = moddims(im, Nx, Ny, Nz);
			//af::array out = af::constant(0.f, padd.dims(0) + 2 * Ndx, padd.dims(1) + 2 * Ndy, padd.dims(2) + 2 * Ndz);
			//out(static_cast<double>(Ndx) + af::seq(static_cast<double>(padd.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(padd.dims(1))),
			//	static_cast<double>(Ndz) + af::seq(static_cast<double>(padd.dims(2)))) = padd;
		}
		padd = out;
	}
	return padd;
}

// Compute the (isotropic) TV prior
af::array TVprior(const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const TVdata &S, const af::array& ima, const float epps, const uint32_t TVtype, 
	const Weighting & w_vec, const af::array& offsets) {
	af::array gradi;

	if (TVtype != 3) {
		af::array im = af::moddims(ima, Nx, Ny, Nz);
		// 1st order differentials
		af::array g = af::constant(0.f, Nx, Ny, Nz);
		af::array f = af::constant(0.f, Nx, Ny, Nz);
		af::array h = af::constant(0.f, Nx, Ny, Nz);
		f(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(im);
		f(Nx - 1u, af::span, af::span) = im(im.dims(0) - 1ULL, af::span, af::span) - im(0, af::span, af::span);
		g(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(im, 1);
		g(af::span, Ny - 1u, af::span) = im(af::span, im.dims(1) - 1ULL, af::span) - im(af::span, 0, af::span);
		h(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(im, 2);
		h(af::span, af::span, Nz - 1u) = im(af::span, af::span, im.dims(2) - 1ULL) - im(af::span, af::span, 0);

		af::array pval, apu1, apu2, apu3, apu4;
		g = af::flat(g);
		f = af::flat(f);
		h = af::flat(h);

		// If anatomical prior is used
		if (S.TV_use_anatomical) {
			if (TVtype == 1) {
				pval = af::sqrt(S.s1*af::pow(f, 2.) + S.s5*af::pow(g, 2.) + S.s9*af::pow(h, 2.) + S.s4*f*g + S.s7*f*h + S.s2*f*g + S.s8*h*g + S.s3*f*h + S.s6*h*g + S.TVsmoothing);
				apu1 = 0.5f * (2.f * S.s1*f + S.s4*g + S.s7*h + S.s2*g + S.s3*h) / pval;
				apu2 = 0.5f * (2.f * S.s5*g + S.s4*f + S.s2*f + S.s8*h + S.s6*h) / pval;
				apu3 = 0.5f * (2.f * S.s9*h + S.s8*g + S.s6*g + S.s7*f + S.s3*f) / pval;
				apu4 = 0.5f * (2.f * S.s1*f + 2.f * S.s5*g + 2.f * S.s9*h + S.s4*f + S.s2*f + S.s8*h + S.s6*h + S.s4*g + S.s7*h + S.s2*g + S.s3*h
					+ S.s8*g + S.s6*g + S.s7*f + S.s3*f) / pval;
			}
			else if (TVtype == 2) {
				af::array gp = af::constant(0.f, Nx, Ny, Nz);
				af::array fp = af::constant(0.f, Nx, Ny, Nz);
				af::array hp = af::constant(0.f, Nx, Ny, Nz);
				fp(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(S.reference_image);
				fp(Nx - 1u, af::span, af::span) = S.reference_image(im.dims(0) - 1ULL, af::span, af::span) - S.reference_image(0, af::span, af::span);
				gp(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(S.reference_image, 1);
				gp(af::span, Ny - 1u, af::span) = S.reference_image(af::span, im.dims(1) - 1ULL, af::span) - S.reference_image(af::span, 0, af::span);
				hp(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(S.reference_image, 2);
				hp(af::span, af::span, Nz - 1u) = S.reference_image(af::span, af::span, im.dims(2) - 1ULL) - S.reference_image(af::span, af::span, 0);

				pval = af::sqrt(af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) + S.T * (af::pow(fp, 2.) + af::pow(gp, 2.) + af::pow(hp, 2.)) + S.TVsmoothing);
				apu1 = f / pval;
				apu2 = g / pval;
				apu3 = h / pval;
				apu4 = (f + g + h) / pval;
			}
			// For APLS
			else if (TVtype == 4) {
				af::array gp = af::constant(0.f, Nx, Ny, Nz);
				af::array fp = af::constant(0.f, Nx, Ny, Nz);
				af::array hp = af::constant(0.f, Nx, Ny, Nz);
				fp(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(S.APLSReference);
				fp(Nx - 1u, af::span, af::span) = S.APLSReference(im.dims(0) - 1ULL, af::span, af::span) - S.APLSReference(0, af::span, af::span);
				gp(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(S.APLSReference, 1);
				gp(af::span, Ny - 1u, af::span) = S.APLSReference(af::span, im.dims(1) - 1ULL, af::span) - S.APLSReference(af::span, 0, af::span);
				hp(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(S.APLSReference, 2);
				hp(af::span, af::span, Nz - 1u) = S.APLSReference(af::span, af::span, im.dims(2) - 1ULL) - S.APLSReference(af::span, af::span, 0);

				fp = af::flat(fp);
				gp = af::flat(gp);
				hp = af::flat(hp);

				af::array epsilon = af::join(1, fp, gp, hp) / af::join(1, fp / af::sqrt(af::pow(fp, 2.) + S.eta * S.eta) + epps,
					gp / af::sqrt(af::pow(gp, 2.) + S.eta * S.eta) + epps, hp / af::sqrt(af::pow(hp, 2.) + S.eta * S.eta) + epps);
				af::array apu = sum(af::join(1, f, g, h) * epsilon, 1);

				pval = (af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) - af::pow(apu, 2.) + S.APLSsmoothing);
				pval(pval < 0) = 1.f;
				pval = af::sqrt(pval);
				apu1 = (f - (apu * epsilon(af::span, 0))) / pval;
				apu2 = (g - (apu * epsilon(af::span, 1))) / pval;
				apu3 = (h - (apu * epsilon(af::span, 2))) / pval;
				apu4 = (f - (apu * epsilon(af::span, 0)) + g - (apu * epsilon(af::span, 1)) + h - (apu * epsilon(af::span, 2))) / pval;
			}
		}
		// If anatomical prior is not used
		else {
			pval = af::sqrt(af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) + S.TVsmoothing);
			apu1 = f / pval;
			apu2 = g / pval;
			apu3 = h / pval;
			apu4 = (f + g + h) / pval;
		}
		apu1 = af::moddims(apu1, Nx, Ny, Nz);
		apu2 = af::moddims(apu2, Nx, Ny, Nz);
		apu3 = af::moddims(apu3, Nx, Ny, Nz);
		apu4 = af::moddims(apu4, Nx, Ny, Nz);
		// Derivatives
		apu1 = af::shift(apu1, 1);
		apu2 = af::shift(apu2, 0, 1);
		apu3 = af::shift(apu3, 0, 0, 1);
		gradi = apu4 - apu1 - apu2 - apu3;
		gradi = af::flat(gradi);
		//im = af::flat(im);
		gradi = batchFunc(gradi, 2.f*S.tau*((af::min))(ima), batchPlus);
	}
	else {
		if (S.TV_use_anatomical) {
			af::array padd = af::flat(padding(ima, Nx, Ny, Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd = padd(af::join(1, offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1)), offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL))));
			padd = af::moddims(padd, ima.dims(0), padd.dims(0) / ima.dims(0));
			af::array padd2 = af::flat(padding(S.reference_image, Nx, Ny, Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd2 = padd2(af::join(1, offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1)), offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL))));
			padd2 = af::moddims(padd2, ima.dims(0), padd2.dims(0) / ima.dims(0));
			gradi = af::sum(af::batchFunc(af::batchFunc(ima, padd, batchMinus) / std::pow(S.C, 2.f) * (1.f / af::sqrt(1.f + af::pow(af::batchFunc(ima, padd, batchMinus) / S.C, 2.) 
				+ af::pow(af::batchFunc(S.reference_image, padd2, batchMinus) / S.T, 2.))), w_vec.weights_TV.T(), batchMul), 1);
		}
		else {
			af::array padd = af::flat(padding(ima, Nx, Ny, Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd = padd(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
			padd = af::moddims(padd, ima.dims(0), padd.dims(0) / ima.dims(0));
			gradi = af::sum(af::batchFunc(af::batchFunc(ima, padd, batchMinus) / std::pow(S.C, 2.f) * (1.f / af::sqrt(1.f + af::pow(af::batchFunc(ima, padd, batchMinus) / S.C, 2.))), 
				w_vec.weights_TV.T(), batchMul), 1);
		}
	}

	return gradi;
}

//af::array MLEM(const af::array &im, const af::array &Summ, const af::array &rhs)
//{
//	//im = im / Summ * rhs;
//	return (im / Summ * rhs);
//}

//af::array OSL_MLEM(const af::array &im, const af::array &Summ, const af::array &rhs, const af::array &dU, const float beta)
//{
//	//im = im / (Summ + beta * dU) * rhs;
//	return (im / (Summ + beta * dU) * rhs);
//}

af::array EM(const af::array &im, const af::array &Summ, const af::array &rhs)
{
	//im = im / Summ * rhs;
	return (im / Summ * rhs);
}

af::array OSL(const af::array &Summ, const af::array &dU, const float beta, const float epps)
{
	af::array output = (Summ + beta * dU);
	output(output < epps) = epps;
	return output;
}

af::array MBSREM(const af::array & im, const af::array & rhs, const float U, const af::array & pj3, const float* lam, const uint32_t iter, const uint32_t im_dim,
	const float beta, const af::array &dU, const af::array & Summ, const float epps)
{
	af::array UU = af::constant(0.f, im_dim);
	af::array output;
	const af::array pp = im < (U / 2.f);
	UU(pp) = im(pp) / (pj3(pp));
	UU(!pp) = (U - im(!pp)) / (pj3(!pp));
	if (beta == 0.f)
		output = im + lam[iter] * UU * rhs;
	else
		output = im + lam[iter] * UU * (rhs - beta * dU);
	//if (beta == 0.f)
	//	output = im + lam[iter] * UU*(rhs - Summ);
	//else
	//	output = im + lam[iter] * UU*(rhs - Summ - beta * dU);
	output(output < epps) = epps;
	output(output >= U) = U - epps;
	return output;
}

af::array BSREM(const af::array & im, const af::array & rhs, const float * lam, const uint32_t iter)
{
	//im += lam[iter] * im * rhs;
	//im = (1.f - lam[iter]) * im + lam[iter] * im * rhs;
	return ((1.f - lam[iter]) * im + lam[iter] * im * rhs);
}

af::array ECOSEM(const af::array & im, const af::array & D, const af::array & OSEM_apu, const af::array & COSEM_apu, const float epps)
{
	float alpha_eco = 1.f;
	af::array output = alpha_eco * OSEM_apu + (1.f - alpha_eco)*COSEM_apu;
	float eco_s1 = af::sum<float>(D * (-COSEM_apu * af::log(im + epps) + im));
	float eco_s2 = af::sum<float>(D * (-COSEM_apu * af::log(output + epps) + output));
	while (alpha_eco > 0.0096f && eco_s1 < eco_s2) {
		alpha_eco *= 0.9f;
		output = alpha_eco * OSEM_apu + (1.f - alpha_eco)*COSEM_apu;
		eco_s2 = af::sum<float>(D*(-COSEM_apu * af::log(output + epps) + output));
	}
	if (alpha_eco <= 0.0096f)
		output = COSEM_apu;
	return output;
}

af::array ROSEM(const af::array & im, const af::array & Summ, const af::array & rhs, const float * lam, const uint32_t iter)
{
	//im = im + lam[iter] * im / Summ * rhs;
	//im = (1.f - lam[iter]) * im + lam[iter] * im / Summ * rhs;
	//return ((1.f - lam[iter]) * im + lam[iter] * im / Summ * rhs);
	return (im + lam[iter] * im / Summ * (rhs - Summ));
}

af::array RBI(const af::array & im, const af::array & Summ, const af::array & rhs, const af::array& D, const float beta, const af::array& dU)
{
	af::array output = im;
	if (beta == 0.f) {
		const float Summa = 1.f / af::max<float>(Summ / D);
		//output += Summa * (rhs - Summa * Summ);
		//output += (im / Summa) * (rhs - Summ);
		output += Summa * (im / D) * (rhs - Summ);
	}
	else {
		const float Summa = 1.f / af::max<float>((Summ + beta * dU) / (D + beta * dU));
		output += Summa * (im / (D + beta * dU)) * (rhs - Summ - beta * dU);
	}
	return output;
}

af::array DRAMA(const af::array & im, const af::array & Summ, const af::array & rhs, const float * lam, const uint32_t iter, const uint32_t sub_iter, const uint32_t subsets)
{
	return (im + lam[iter * subsets + sub_iter] * im / Summ * (rhs - Summ));
}

af::array MAP(const af::array & im, const float lam, const float beta, const af::array & dU, const float epps)
{
	af::array output = im - beta * lam * im * dU;
	output(output < epps) = epps;
	return output;
}

af::array COSEM(const af::array & im, const af::array & C_co, const af::array & D, const float h, const uint32_t COSEM_TYPE)
{
	af::array output;
	if (COSEM_TYPE == 1) {
		//if (beta == 0.f)
			output = af::pow(af::sum(C_co, 1) / D, h);
		//else
			//output = af::pow(af::sum(C_co, 1) / (D + beta * dU), h);
	}
	else {
		//if (beta == 0.f)
			output = (af::sum(C_co, 1) / D);
		//else
			//output = (af::sum(C_co, 1) / (D + beta * dU));
	}
	return output;
}

af::array MRP(const af::array & im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps,
	const af::array& offsets, const bool med_no_norm, const uint32_t im_dim)
{
	af::array padd = af::flat(padding(im, Nx, Ny, Nz, medx, medy, medz));
	//af::array grad = af::constant(0.f, im_dim, 1);
	//gfor(af::seq i, im_dim) {
	//	af::array apu = offsets(i, af::span);
	//	af::array testi = padd(apu);
	//	grad(i) = af::median(testi);
	//}
	af::array grad = af::median(af::moddims(padd(af::flat(offsets)), im_dim, offsets.dims(1)), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array Quadratic_prior(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t inffi,
	const af::array &offsets, const af::array &weights_quad, const uint32_t im_dim)
{
	//const af::array apu_pad = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));
	//af::array indeksi1 = offsets.col(inffi);
	//const af::array indeksi2 = af::join(1, offsets.cols(0, inffi - 1), offsets.cols(inffi + 1, af::end));
	//af::array grad = af::sum(af::batchFunc(af::batchFunc(apu_pad(indeksi1 + 0), af::moddims(apu_pad(indeksi2), im_dim, weights_quad.dims(0)), batchMinus), af::transpose(weights_quad), batchMul), 1);
	//int pituus = (Ndx * 2 + 1) * (Ndy * 2 + 1) * (Ndz * 2 + 1);
	const af::array apu_pad = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
	af::array grad;
	if (Ndz == 0 || Nz == 1) {
		grad = af::convolve2(apu_pad, weights_quad);
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::span);
	}
	else {
		grad = af::convolve3(apu_pad, weights_quad);
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::seq(Ndz, Nz + Ndz - 1));
	}
	grad = af::flat(grad);
	return grad;
}

af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t inffi,
	const af::array& offsets, const af::array& weights_huber, const uint32_t im_dim, const float delta)
{
	//const af::array apu_pad = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));
	//af::array indeksi1 = offsets.col(inffi);
	//const af::array indeksi2 = af::join(1, offsets.cols(0, inffi - 1), offsets.cols(inffi + 1, af::end));
	//af::array grad = af::sum(af::batchFunc(af::batchFunc(apu_pad(indeksi1 + 0), af::moddims(apu_pad(indeksi2), im_dim, weights_huber.dims(0)), batchMinus), af::transpose(weights_huber), batchMul), 1);
	af::array grad = Quadratic_prior(im, Ndx, Ndy, Ndz, Nx, Ny, Nz, inffi, offsets, weights_huber, im_dim);
	if (af::sum<dim_t>(delta >= af::abs(af::flat(grad))) == grad.elements() && af::sum<int>(af::flat(grad)) != 0)
		mexPrintf("Delta value of Huber prior larger than all the pixel difference values\n");
	grad(grad > delta) = delta;
	grad(grad < -delta) = -delta;
	return grad;
}

af::array FMH(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const uint32_t inffi,
	const af::array& offsets, const af::array& fmh_weights, const bool med_no_norm, const uint32_t alku_fmh, const uint32_t im_dim)
{
	af::array grad;
	af::array indeksi1;
	const af::array padd = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));
	uint32_t luup;
	if (Nz == 1 || Ndz == 0) {
		grad = af::constant(0.f, im_dim, 5);
		luup = 4;
	}
	else {
		grad = af::constant(0.f, im_dim, 14);
		luup = 13;
	}
	for (uint32_t ii = 0; ii < luup; ii++) {
		indeksi1 = af::flat(offsets(af::span, af::seq(Ndx * ii, offsets.dims(1) - Ndx * (ii) - 1, alku_fmh / Ndx - ii)));
		af::array apu_pad = af::moddims(padd(indeksi1 + 0), im_dim, fmh_weights.dims(0));
		grad(af::span, ii) = af::sum(af::batchFunc(apu_pad, af::transpose(fmh_weights(af::span, ii)), batchMul), 1);
	}
	indeksi1 = offsets.col(inffi);
	grad(af::span, af::end) = padd(indeksi1 + 0);
	grad = af::median(grad, 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array L_filter(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps,
	const af::array& offsets, const af::array& a_L, const bool med_no_norm, const uint32_t im_dim)
{
	af::array grad;
	af::array apu_pad = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));
	apu_pad = apu_pad(af::flat(offsets));
	apu_pad = af::sort(af::moddims(apu_pad, im_dim, a_L.dims(0)), 1);
	//grad = af::matmul(apu_pad, a_L);
	grad = af::sum(af::batchFunc(apu_pad, af::transpose(a_L), batchMul), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array Weighted_mean(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const af::array& offsets,
	const af::array& weighted_weights, const bool med_no_norm, const uint32_t im_dim, const uint32_t mean_type, const float w_sum)
{
	af::array grad;
	//af::array padd = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));
	//padd = af::moddims(padd(offsets), im_dim, offsets.dims(1));
	af::array padd = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
	if (mean_type == 1U) {
		if (Ndz == 0 || Nz == 1)
			grad = af::convolve2(padd, weighted_weights / w_sum);
		else
			grad = af::convolve3(padd, weighted_weights / w_sum);
		//grad = af::sum(af::batchFunc(padd, af::transpose(weighted_weights), batchMul), 1) / w_sum;
	}
	else if (mean_type == 2U) {
		if (Ndz == 0 || Nz == 1)
			grad = af::convolve2(1.f / padd, w_sum / weighted_weights);
		else
			grad = af::convolve3(1.f / padd, w_sum / weighted_weights);
		//grad = w_sum / af::sum(af::batchFunc(padd, af::transpose(weighted_weights), batchMul), 1);
	}
	else if (mean_type == 3U) {
		if (Ndz == 0 || Nz == 1)
			grad = af::exp(af::convolve2(af::log(padd), weighted_weights / w_sum));
		else
			grad = af::exp(af::convolve3(af::log(padd), weighted_weights / w_sum));
		//grad = af::exp(af::sum(af::batchFunc(af::log(padd), af::transpose(weighted_weights), batchMul), 1) / w_sum);
	}
	else {
		mexWarnMsgTxt("Unsupported mean type");
		grad = af::constant(0.f, Nx + Ndx * 2U, Ny + Ndy * 2U, Nz + Ndz * 2U);
	}
	if (Ndz == 0 || Nz == 1) {
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::span);
	}
	else {
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::seq(Ndz, Nz + Ndz - 1));
	}
	grad = af::flat(grad);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array AD(const af::array & im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const float TimeStepAD, const float KAD, const uint32_t NiterAD,
	const af_flux_function FluxType, const af_diffusion_eq DiffusionType, const bool med_no_norm)
{
	const af::array padd = af::moddims(im, Nx, Ny, Nz);
	af::array grad = af::anisotropicDiffusion(padd, TimeStepAD, KAD, NiterAD, FluxType, DiffusionType);
	grad = af::flat(grad);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array TGV(const af::array & im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t maxits, const float alpha, const float beta)
{
	af::array grad = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array u = grad;
	grad = af::flat(grad);
	const af::array imi = af::moddims(im, Nx, Ny, Nz);
	const float sigma = 1.f / 16.f;
	const float tau = 1.f / 8.f;
	af::array p1 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array p2 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array p3 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array v1 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array v2 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array v3 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array q1 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array q2 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array q3 = af::constant(0.f, Nx, Ny, Nz, f32);
	af::array q4 = af::constant(0.f, Nx, Ny, Nz, f32);

	af::array vb1 = v1;
	af::array vb2 = v2;
	af::array vb3 = v3;
	v1 = af::flat(v1);
	v2 = af::flat(v2);
	v3 = af::flat(v3);

	for (uint32_t kk = 0; kk < maxits; kk++) {
		af::array temp_u = u + imi;
		af::array eta1 = af::join(0, af::diff1(temp_u), temp_u(0, af::span, af::span) - temp_u(af::end, af::span, af::span)) - vb1;
		af::array eta2 = af::join(1, af::diff1(temp_u, 1), temp_u(af::span, 0, af::span) - temp_u(af::span, af::end, af::span)) - vb2;
		af::array eta3 = af::join(2, af::diff1(temp_u, 2), temp_u(af::span, af::span, 0) - temp_u(af::span, af::span, af::end)) - vb3;
		eta1 = p1 + sigma * eta1;
		eta2 = p2 + sigma * eta2;
		eta3 = p3 + sigma * eta3;
		af::array apu = (af::max)(1.f, af::sqrt(af::pow(eta1, 2.) + af::pow(eta2, 2.) + af::pow(eta3, 2.)) / beta);
		apu = af::flat(apu);
		p1 = af::flat(eta1) / (apu);
		p1 = af::moddims(p1, Nx, Ny, Nz);
		p2 = af::flat(eta2) / (apu);
		p2 = af::moddims(p2, Nx, Ny, Nz);
		p3 = af::flat(eta3) / (apu);
		p3 = af::moddims(p3, Nx, Ny, Nz);
		eta1 = af::join(0, af::diff1(vb1), vb1(0, af::span, af::span) - vb1(af::end, af::span, af::span));
		eta2 = af::join(1, af::diff1(vb2, 1), vb2(af::span, 0, af::span) - vb2(af::span, af::end, af::span));
		eta3 = af::join(2, af::diff1(vb3, 2), vb3(af::span, af::span, 0) - vb3(af::span, af::span, af::end));
		af::array eta4 = (af::join(0, af::diff1(vb2), vb2(0, af::span, af::span) - vb2(af::end, af::span, af::span)) + 
			af::join(1, af::diff1(vb1, 1), vb1(af::span, 0, af::span) - vb1(af::span, af::end, af::span)) +
			af::join(2, af::diff1(vb2, 2), vb2(af::span, af::span, 0) - vb2(af::span, af::span, af::end)) + 
			af::join(2, af::diff1(vb1, 2), vb1(af::span, af::span, 0) - vb1(af::span, af::span, af::end)) + 
			af::join(0, af::diff1(vb3), vb3(0, af::span, af::span) - vb3(af::end, af::span, af::span)) + 
			af::join(1, af::diff1(vb3, 1), vb3(af::span, 0, af::span) - vb3(af::span, af::end, af::span))) / 6.f;
		eta1 = q1 + sigma * eta1;
		eta2 = q2 + sigma * eta2;
		eta3 = q3 + sigma * eta3;
		eta4 = q4 + sigma * eta4;
		apu = (af::max)(1.f, af::sqrt(af::pow(eta1, 2.) + af::pow(eta2, 2.) + af::pow(eta3, 2.) + af::pow(eta4, 2.)) / alpha);
		apu = af::flat(apu);
		q1 = af::flat(eta1) / (apu);
		q1 = af::moddims(q1, Nx, Ny, Nz);
		q2 = af::flat(eta2) / (apu);
		q2 = af::moddims(q2, Nx, Ny, Nz);
		q3 = af::flat(eta3) / (apu);
		q3 = af::moddims(q3, Nx, Ny, Nz);
		q4 = af::flat(eta4) / (apu);
		q4 = af::moddims(q4, Nx, Ny, Nz);

		eta4 = af::join(0, p1(af::end, af::span, af::span) - p1(0, af::span, af::span), -af::diff1(p1)) +
			af::join(1, p2(af::span, af::end, af::span) - p2(af::span, 0, af::span), -af::diff1(p2, 1)) +
			af::join(2, p3(af::span, af::span, af::end) - p3(af::span, af::span, 0), -af::diff1(p3, 2));
		af::array uold = grad;
		grad -= tau * af::flat(eta4);
		eta1 = af::join(0, q1(af::end, af::span, af::span) - q1(0, af::span, af::span), -af::diff1(q1)) + 
			af::join(1, q4(af::span, af::end, af::span) - q4(af::span, 0, af::span), -af::diff1(q4, 1)) +
			af::join(2, q4(af::span, af::span, af::end) - q4(af::span, af::span, 0), -af::diff1(q4, 2)) - p1;
		eta2 = af::join(0, q4(af::end, af::span, af::span) - q4(0, af::span, af::span), -af::diff1(q4)) +
			af::join(1, q2(af::span, af::end, af::span) - q2(af::span, 0, af::span), -af::diff1(q2, 1)) +
			af::join(2, q4(af::span, af::span, af::end) - q4(af::span, af::span, 0), -af::diff1(q4, 2)) - p2;
		eta3 = af::join(0, q4(af::end, af::span, af::span) - q4(0, af::span, af::span), -af::diff1(q4)) +
			af::join(2, q3(af::span, af::span, af::end) - q3(af::span, af::span, 0), -af::diff1(q3, 2)) +
			af::join(1, q4(af::span, af::end, af::span) - q4(af::span, 0, af::span), -af::diff1(q4, 1)) - p3;
		af::array v1old = v1;
		af::array v2old = v2;
		af::array v3old = v3;
		v1 -= tau * af::flat(eta1);
		v2 -= tau * af::flat(eta2);
		v3 -= tau * af::flat(eta3);
		u = 2.f * grad - uold;
		u = af::moddims(u, Nx, Ny, Nz);
		vb1 = af::moddims(2.f * v1 - v1old, Nx, Ny, Nz);
		vb2 = af::moddims(2.f * v2 - v2old, Nx, Ny, Nz);
		vb3 = af::moddims(2.f * v3 - v3old, Nx, Ny, Nz);

	}

	return -grad;
}

//af::array im2col_3D(const af::array& A, const uint32_t blocksize1, const uint32_t blocksize2, const uint32_t blocksize3) {
//	const uint32_t Nx = A.dims(0);
//	const uint32_t Ny = A.dims(1);
//	const uint32_t Nz = A.dims(2);
//	af::array start_ind = af::flat(af::batchFunc(af::range(af::dim4(Nx - blocksize1 + 1), 0, u32) + 1, af::range(af::dim4(1, Ny - blocksize2 + 1), 1, u32) * Nx, batchPlus));
//	start_ind = af::reorder(af::transpose(af::batchFunc(start_ind, af::range(af::dim4(1, blocksize1), 1, u32), batchPlus)), 0, 2, 1);
//	start_ind = af::moddims(af::batchFunc(start_ind, af::range(af::dim4(1, blocksize2), 1, u32) * Nx, batchPlus), blocksize1 * blocksize2, start_ind.dims(2));
//	start_ind = af::batchFunc(start_ind, Nx * Ny * af::reorder(af::range(af::dim4(1, Nz - blocksize3 + 1), 1, u32), 0, 2, 1), batchPlus);
//	start_ind = af::tile(start_ind, blocksize3, 1, 1);
//	af::array apu = af::transpose(af::range(af::dim4(blocksize3), 0, u32) * Nx * Ny);
//	apu = af::flat(af::tile(apu, blocksize1 * blocksize2, 1));
//	start_ind = af::batchFunc(start_ind, apu, batchPlus);
//	return A(start_ind - 1);
//}
//
//af::array sparseSum(const af::array& W, const uint32_t s) {
//	af::array summa = af::constant(0.f, s, 1, f32);
//	af::array temp_row = af::sparseGetRowIdx(W);
//	af::array temp_val = af::sparseGetValues(W);
//	for (int kk = 0; kk < s; kk++) {
//		const int32_t* ind1 = temp_row(kk).host<int32_t>();
//		const int32_t* ind2 = temp_row(kk + 1).host<int32_t>();
//		if (*ind1 == *ind2)
//			continue;
//		summa(kk) = af::sum(temp_val(af::seq(static_cast<double>(*ind1), static_cast<double>(*ind2 - 1))));
//	}
//	return summa;
//}

af::array computeConvolution(const af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec, 
	const uint32_t n_rekos) {
	af::array apu = af::moddims(vec, Nx, Ny, Nz, n_rekos);
	//af::dim4 pad_dim1 = { w_vec.g_dim_x , w_vec.g_dim_y, w_vec.g_dim_z, 0 };
	//af::dim4 pad_dim2 = { w_vec.g_dim_x , w_vec.g_dim_y, w_vec.g_dim_z, 0 };
	//array apu = vec.im_os;
	for (int ff = 0; ff < n_rekos; ff++) {
		af::array apu2 = padding(apu(af::span, af::span, af::span, ff), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
		//af::array apu2 = apu(af::span, af::span, af::span, ff);
		//apu2 = af::pad(apu2, pad_dim1, pad_dim2, AF_PAD_SYM);
		apu2 = af::convolve3(apu2, g);
		apu2 = apu2(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		apu(af::span, af::span, af::span, ff) = apu2;
	}
	return af::flat(apu);
}

void deblur(af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec, const uint32_t iter, 
	const uint32_t subsets, const float epps) {
	//af::dim4 pad_dim1 = { w_vec.g_dim_x , w_vec.g_dim_y, w_vec.g_dim_z, 0 };
	//af::dim4 pad_dim2 = { w_vec.g_dim_x , w_vec.g_dim_y, w_vec.g_dim_z, 0 };
	af::array jelppi = padding(vec(af::span, iter + 1u), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	af::array apu = padding(vec(af::span, iter + 1u), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	//af::array jelppi = af::moddims(vec(af::span, iter + 1u), Nx, Ny, Nz);
	//af::array apu = af::moddims(vec(af::span, iter + 1u), Nx, Ny, Nz);
	//jelppi = af::pad(jelppi, pad_dim1, pad_dim2, AF_PAD_SYM);
	//apu = af::pad(apu, pad_dim1, pad_dim2, AF_PAD_SYM);
	for (int kk = 0; kk < subsets; kk++) {
		af::array apu2 = convolve3(jelppi, g) + epps;
		apu2 = apu2(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		apu2 = padding(apu2, Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
		jelppi *= af::convolve3(apu / apu2, g);
		jelppi = jelppi(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		if (kk < subsets - 1)
			jelppi = padding(jelppi, Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	}
	jelppi = af::flat(jelppi);
	vec(af::span, iter + 1u) = jelppi;
}

void computeDeblur(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps) {
	if (MethodList.OSEM) {
		deblur(vec.OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.MRAMLA) {
		deblur(vec.MRAMLA, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.RAMLA) {
		deblur(vec.RAMLA, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.ROSEM) {
		deblur(vec.ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.RBI) {
		deblur(vec.RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.DRAMA) {
		deblur(vec.DRAMA, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.COSEM) {
		deblur(vec.COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.ECOSEM) {
		deblur(vec.ECOSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.ACOSEM) {
		deblur(vec.ACOSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}

	if (MethodList.MRP) {
		if (MethodList.OSLOSEM) {
			deblur(vec.MRP_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.MRP_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.MRP_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.MRP_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.MRP_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.MRP_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}

	if (MethodList.Quad) {
		if (MethodList.OSLOSEM) {
			deblur(vec.Quad_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.Quad_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.Quad_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.Quad_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.Quad_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.Quad_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}

	if (MethodList.Huber) {
		if (MethodList.OSLOSEM) {
			deblur(vec.Huber_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.Huber_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.Huber_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.Huber_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.Huber_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.Huber_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}

	if (MethodList.L) {
		if (MethodList.OSLOSEM) {
			deblur(vec.L_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.L_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.L_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.L_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.L_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.L_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM) {
			deblur(vec.FMH_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.FMH_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.FMH_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.FMH_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.FMH_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.FMH_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM) {
			deblur(vec.Weighted_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.Weighted_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.Weighted_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.Weighted_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.Weighted_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.Weighted_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.TV) {
		if (MethodList.OSLOSEM) {
			deblur(vec.TV_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.TV_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.TV_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.TV_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.TV_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.TV_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.AD) {
		if (MethodList.OSLOSEM) {
			deblur(vec.AD_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.AD_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.AD_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.AD_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.AD_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.AD_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM) {
			deblur(vec.APLS_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.APLS_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.APLS_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.APLS_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.APLS_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.APLS_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.TGV) {
		if (MethodList.OSLOSEM) {
			deblur(vec.TGV_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.TGV_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.TGV_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.TGV_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.TGV_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.TGV_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.NLM) {
		if (MethodList.OSLOSEM) {
			deblur(vec.NLM_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.NLM_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.NLM_MBSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.NLM_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.NLM_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.NLM_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			deblur(vec.custom_OSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.BSREM) {
			deblur(vec.custom_BSREM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.MBSREM) {
			deblur(vec.ECOSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.ROSEMMAP) {
			deblur(vec.custom_ROSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.RBIOSL) {
			deblur(vec.custom_RBI, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.OSLCOSEM > 0) {
			deblur(vec.custom_COSEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
}


void computeDeblurMLEM(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps) {
	if (MethodList.MLEM) {
		deblur(vec.MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
	}
	if (MethodList.OSLMLEM) {
		if (MethodList.MRP) {
			deblur(vec.MRP_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.Quad) {
			deblur(vec.Quad_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.Huber) {
			deblur(vec.Huber_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.L) {
			deblur(vec.L_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.FMH) {
			deblur(vec.FMH_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.WeightedMean) {
			deblur(vec.Weighted_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.TV) {
			deblur(vec.TV_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.AD) {
			deblur(vec.AD_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.APLS) {
			deblur(vec.APLS_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.TGV) {
			deblur(vec.TGV_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
		if (MethodList.NLM) {
			deblur(vec.NLM_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
	if (MethodList.CUSTOM) {
		if (MethodList.OSLMLEM) {
			deblur(vec.custom_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps);
		}
	}
}