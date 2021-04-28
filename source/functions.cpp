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

#ifndef OPENCL
const char* getErrorString(CUresult error)
{
	const char* errstr;
	cuGetErrorString(error, &errstr);
	return errstr;
}
#endif

// Loads the input data and forms device data variables
void form_data_variables(AF_im_vectors & vec, std::vector<float> & beta, Weighting & w_vec, const mxArray *options, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t Niter, const af::array &x0, const uint32_t im_dim, const size_t koko_l, const RecMethods &MethodList, TVdata &data, 
	const uint32_t subsets, const uint32_t osa_iter0, const bool use_psf, const bool saveIter, const uint32_t Nt, const uint32_t iter0, const bool CT)
{
	// Load the number of priors, all MAP-algorithms, non-OS MAP algorithms, non-OS non-MAP algorithms, non-MAP algorithms and total number of algorithms
	w_vec.nPriors = getScalarUInt32(mxGetField(options, 0, "nPriors"), -2);
	if (MethodList.CUSTOM)
		w_vec.nPriorsTot = w_vec.nPriors + 1U;
	else
		w_vec.nPriorsTot = w_vec.nPriors;
	w_vec.nMAP = getScalarUInt32(mxGetField(options, 0, "nMAP"), -2);
	w_vec.nMAPML = getScalarUInt32(mxGetField(options, 0, "nMAPML"), -3);
	w_vec.nMAPOS = w_vec.nMAP - w_vec.nMAPML;
	w_vec.nMLEM = getScalarUInt32(mxGetField(options, 0, "nMLEM"), -4);
	w_vec.nOS = getScalarUInt32(mxGetField(options, 0, "nOS"), -5);
	w_vec.nTot = getScalarUInt32(mxGetField(options, 0, "nTot"), -6);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	w_vec.rekot = (uint32_t*)mxGetUint32s(mxGetField(options, 0, "rekoList"));
#else
	w_vec.rekot = (uint32_t*)mxGetData(mxGetField(options, 0, "rekoList"));
#endif
	w_vec.mIt.push_back(getScalarInt32(mxGetCell(mxGetField(options, 0, "mIt"), 0), -7));
	w_vec.mIt.push_back(getScalarInt32(mxGetCell(mxGetField(options, 0, "mIt"), 1), -8));
	if (DEBUG) {
		mexPrintf("nPriors = %u\n", w_vec.nPriors);
		mexPrintf("nMAP = %u\n", w_vec.nMAP);
		mexPrintf("nMAPML = %u\n", w_vec.nMAPML);
		mexPrintf("nTot = %u\n", w_vec.nTot);
		mexPrintf("nMLEM = %u\n", w_vec.nMLEM);
		mexEvalString("pause(.0001);");
	}
	uint32_t Ni = 1U;
	if (saveIter)
		Ni = Niter + 1U;
	// Load the necessary variables if the corresponding reconstruction method is used and set the initial value
	int yy = 0;
	// First non-MAP/prior-based algorithms (e.g. MLEM)
	if (MethodList.MLEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.OSEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.MRAMLA) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.RAMLA) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.ROSEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.RBI) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
	}
	
	if (MethodList.DRAMA) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
		
		// Relaxation parameter
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.lambda_DRAMA = (float*)mxGetSingles(mxGetField(options, 0, "lam_drama"));
#else
		w_vec.lambda_DRAMA = (float*)mxGetData(mxGetField(options, 0, "lam_drama"));
#endif
	}
	
	if (MethodList.COSEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
		
		// Complete data
		vec.C_co = af::constant(0.f, im_dim, subsets);
	}
	if (MethodList.ECOSEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
		
		if (!MethodList.COSEM) {
			
			// Complete data
			vec.C_co = af::constant(0.f, im_dim, subsets);
		}
	}
	if (MethodList.ACOSEM) {
		vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
		vec.imEstimates[yy](af::span, 0) = x0;
		yy++;
		
		// Complete data
		vec.C_aco = af::constant(0.f, im_dim, subsets);
	}

	if (MethodList.OSLCOSEM > 0)
		vec.C_osl = af::constant(0.f, im_dim, subsets);

	// Load the MAP/prior-based algorithms
	int tt = 0;
	for (uint32_t kk = 0; kk < w_vec.nPriors; kk++) {
		for (uint32_t ll = 0; ll < w_vec.nMAP; ll++) {
			const char* varChar = mxArrayToString(mxGetCell(mxGetField(options, 0, "varList"), tt));
			vec.imEstimates.push_back(af::constant(0.f, im_dim, Ni));
			vec.imEstimates[yy](af::span, 0) = x0;
			if (DEBUG) {
				mexPrintf("%s\n", varChar);
				mexEvalString("pause(.0001);");
			}
			// Load the regularization parameter as well if the prior is used
			beta.push_back(getScalarFloat(mxGetField(options, 0, varChar), -9));
			yy++;
			tt++;
		}
	}

	// CT-related variables such as number of projection images
	if (CT) {
		w_vec.size_y = getScalarUInt32(mxGetField(options, 0, "xSize"), -10);
		w_vec.nProjections = getScalarInt64(mxGetField(options, 0, "nProjections"), -11);
		w_vec.dPitch = getScalarFloat(mxGetField(options, 0, "dPitch"), -12);
	}
	// Load TV related input data
	if (MethodList.TV && MethodList.MAP) {
		// Is anatomical reference image used
		data.TV_use_anatomical = getScalarBool(mxGetField(options, 0, "TV_use_anatomical"), -13);
		// Tau-value
		data.tau = getScalarFloat(mxGetField(options, 0, "tau"), -14);
		// "Smoothing" parameter, prevents zero values in the square root
		data.TVsmoothing = getScalarFloat(mxGetField(options, 0, "TVsmoothing"), -15);
		// The type of TV prior used
		data.TVtype = getScalarUInt32(mxGetField(options, 0, "TVtype"), -16);
		// If anatomical prior is used, load the necessary coefficients
		if (data.TV_use_anatomical) {
			mxArray* TVdata_init = mxGetField(options, 0, "TVdata");
			if (data.TVtype == 1) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				data.s1 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s1")), afHost);
				data.s2 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s2")), afHost);
				data.s3 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s3")), afHost);
				data.s4 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s4")), afHost);
				data.s5 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s5")), afHost);
				data.s6 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s6")), afHost);
				data.s7 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s7")), afHost);
				data.s8 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s8")), afHost);
				data.s9 = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "s9")), afHost);
#else
				data.s1 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s1")), afHost);
				data.s2 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s2")), afHost);
				data.s3 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s3")), afHost);
				data.s4 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s4")), afHost);
				data.s5 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s5")), afHost);
				data.s6 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s6")), afHost);
				data.s7 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s7")), afHost);
				data.s8 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s8")), afHost);
				data.s9 = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "s9")), afHost);
#endif
			}
			else {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				data.reference_image = af::array(im_dim, (float*)mxGetSingles(mxGetField(TVdata_init, 0, "reference_image")), afHost);
#else
				data.reference_image = af::array(im_dim, (float*)mxGetData(mxGetField(TVdata_init, 0, "reference_image")), afHost);
#endif
			}
			data.T = getScalarFloat(mxGetField(options, 0, "T"), -17);
			data.C = getScalarFloat(mxGetField(options, 0, "C"), -18);
		}
		// Additional weights for the TV type 3
		if (data.TVtype == 3 && !MethodList.Quad) {
			w_vec.Ndx = getScalarUInt32(mxGetField(options, 0, "Ndx"), -19);
			w_vec.Ndy = getScalarUInt32(mxGetField(options, 0, "Ndy"), -20);
			w_vec.Ndz = getScalarUInt32(mxGetField(options, 0, "Ndz"), -21);
			w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			w_vec.weights_TV = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(mxGetField(options, 0, "weights_quad")), afHost);
#else
			w_vec.weights_TV = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_quad")), afHost);
#endif
		}
		if (data.TVtype == 3)
			data.C = getScalarFloat(mxGetField(options, 0, "C"), -22);
		if (data.TVtype == 4)
			data.SATVPhi = getScalarFloat(mxGetField(options, 0, "SATVPhi"), -23);
	}
	// General variables for neighborhood-based methods
	if ((MethodList.L || MethodList.FMH || MethodList.WeightedMean || MethodList.Quad || MethodList.Huber || MethodList.MRP || MethodList.NLM || (data.TVtype == 3 && MethodList.TV) || MethodList.RDP) && MethodList.MAP) {
		// Neighborhood size
		w_vec.Ndx = getScalarUInt32(mxGetField(options, 0, "Ndx"), -24);
		w_vec.Ndy = getScalarUInt32(mxGetField(options, 0, "Ndy"), -25);
		w_vec.Ndz = getScalarUInt32(mxGetField(options, 0, "Ndz"), -26);
		// Is normalization used in MRP, FMH, L, weighted mean or AD
		w_vec.med_no_norm = getScalarBool(mxGetField(options, 0, "med_no_norm"), -27);
		w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
	}
	if ((MethodList.L || MethodList.FMH || (data.TVtype == 3 && MethodList.TV) || MethodList.RDP) && MethodList.MAP) {
		// Index values for the neighborhood
//#ifdef OPENCL
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.tr_offsets = af::array(im_dim, w_vec.dimmu, (uint32_t*)mxGetUint32s(mxGetField(options, 0, "tr_offsets")), afHost);
#else
		w_vec.tr_offsets = af::array(im_dim, w_vec.dimmu, (uint32_t*)mxGetData(mxGetField(options, 0, "tr_offsets")), afHost);
#endif
//#else
//		w_vec.tr_offsets = af::array(im_dim, w_vec.dimmu, (uint32_t*)mxGetData(mxGetField(options, 0, "tr_offsets")), afHost);
//#endif
	}
	if (MethodList.FMH || MethodList.Quad || MethodList.Huber || MethodList.RDP)
		w_vec.inffi = getScalarUInt32(mxGetField(options, 0, "inffi"), -28);
	// Weights for the various priors
	if (MethodList.Quad && MethodList.MAP) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(mxGetField(options, 0, "weights_quad")), afHost);
#else
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_quad")), afHost);
#endif
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_quad * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_quad = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);

	}
	if (MethodList.RDP && MethodList.MAP) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.weights_RDP = af::array((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, (float*)mxGetSingles(mxGetField(options, 0, "weights_RDP")), afHost);
#else
		w_vec.weights_RDP = af::array((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, (float*)mxGetData(mxGetField(options, 0, "weights_RDP")), afHost);
#endif
		w_vec.weights_RDP = af::flat(w_vec.weights_RDP);
		w_vec.RDP_gamma = getScalarFloat(mxGetField(options, 0, "RDP_gamma"), -29);
	}
	if (MethodList.Huber && MethodList.MAP) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(mxGetField(options, 0, "weights_huber")), afHost);
#else
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, (float*)mxGetData(mxGetField(options, 0, "weights_huber")), afHost);
#endif
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_huber * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_huber = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
		w_vec.huber_delta = getScalarFloat(mxGetField(options, 0, "huber_delta"), -29);
	}
	if (MethodList.L && MethodList.MAP)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.a_L = af::array(w_vec.dimmu, (float*)mxGetSingles(mxGetField(options, 0, "a_L")), afHost);
#else
		w_vec.a_L = af::array(w_vec.dimmu, (float*)mxGetData(mxGetField(options, 0, "a_L")), afHost);
#endif
	if (MethodList.FMH && MethodList.MAP) {
		if (Nz == 1 || w_vec.Ndz == 0)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, (float*)mxGetSingles(mxGetField(options, 0, "fmh_weights")), afHost);
#else
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, (float*)mxGetData(mxGetField(options, 0, "fmh_weights")), afHost);
#endif
		else
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, (float*)mxGetSingles(mxGetField(options, 0, "fmh_weights")), afHost);
#else
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, (float*)mxGetData(mxGetField(options, 0, "fmh_weights")), afHost);
#endif
		w_vec.alku_fmh = getScalarUInt32(mxGetField(options, 0, "inffi"), -30);
	}
	if (MethodList.WeightedMean && MethodList.MAP) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, (float*)mxGetSingles(mxGetField(options, 0, "weighted_weights")), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
#else
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, (float*)mxGetData(mxGetField(options, 0, "weighted_weights")), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
#endif
		// Type of mean used (arithmetic, harmonic or geometric)
		w_vec.mean_type = getScalarInt32(mxGetField(options, 0, "mean_type"), -31);
		// Sum of the weights
		w_vec.w_sum = getScalarFloat(mxGetField(options, 0, "w_sum"), -31);
	}
	if (MethodList.AD && MethodList.MAP) {
		// Time-step value
		w_vec.TimeStepAD = getScalarFloat(mxGetField(options, 0, "TimeStepAD"), -32);
		// Conductance (edge value)
		w_vec.KAD = getScalarFloat(mxGetField(options, 0, "KAD"), -33);
		// Number of AD iterations
		w_vec.NiterAD = getScalarUInt32(mxGetField(options, 0, "NiterAD"), -34);
		// Flux type
		uint32_t Flux = getScalarUInt32(mxGetField(options, 0, "FluxType"), -35);
		// Diffusion type
		uint32_t Diffusion = getScalarUInt32(mxGetField(options, 0, "DiffusionType"), -36);
		if (Flux == 2U)
			w_vec.FluxType = AF_FLUX_QUADRATIC;
		else
			w_vec.FluxType = AF_FLUX_EXPONENTIAL;
		if (Diffusion == 2U)
			w_vec.DiffusionType = AF_DIFFUSION_MCDE;
		else
			w_vec.DiffusionType = AF_DIFFUSION_GRAD;
		if (!MethodList.L && !MethodList.FMH && !MethodList.WeightedMean && !MethodList.Quad && !MethodList.MRP)
			w_vec.med_no_norm = getScalarBool(mxGetField(options, 0, "med_no_norm"), -37);
	}
	if (MethodList.APLS && MethodList.MAP) {
		// Eta value
		data.eta = getScalarFloat(mxGetField(options, 0, "eta"), -38);
		// Tau-value
		if (!MethodList.TV)
			data.tau = getScalarFloat(mxGetField(options, 0, "tau"), -39);
		// Smoothing value
		data.APLSsmoothing = getScalarFloat(mxGetField(options, 0, "APLSsmoothing"), -40);
		// Anatomical reference image
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		data.APLSReference = af::array(Nx, Ny, Nz, (float*)mxGetSingles(mxGetField(options, 0, "APLS_ref_image")), afHost);
#else
		data.APLSReference = af::array(Nx, Ny, Nz, (float*)mxGetData(mxGetField(options, 0, "APLS_ref_image")), afHost);
#endif
	}
	if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.COSEM
		|| MethodList.ECOSEM || MethodList.ACOSEM || MethodList.PKMA || MethodList.OSLCOSEM > 0) {
		// Sum of the rows (measurements) of the system matrix
		//w_vec.D = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "D")), afHost);
		w_vec.D = af::constant(0.f, im_dim, 1);
		// For manual determination of the upper bound
		if ((MethodList.MRAMLA || MethodList.MBSREM) && Nt > 1U)
			w_vec.Amin = af::constant(0.f, koko_l, 1);
		else
			w_vec.Amin = af::constant(0.f, 1, 1);
		//w_vec.Amin = af::array(koko_l, (float*)mxGetData(mxGetField(options, 0, "Amin")), afHost);
		w_vec.MBSREM_prepass = getScalarBool(mxGetField(options, 0, "MBSREM_prepass"), -41);
	}
	if (MethodList.MRAMLA || MethodList.MBSREM) {
		// Relaxation parameter
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.lambda_MBSREM = (float*)mxGetSingles(mxGetField(options, 0, "lam_MBSREM"));
#else
		w_vec.lambda_MBSREM = (float*)mxGetData(mxGetField(options, 0, "lam_MBSREM"));
#endif
		// Upper bound
		w_vec.U = getScalarFloat(mxGetField(options, 0, "U"), -42);
	}
	if (DEBUG) {
		mexPrintf("w_vec.lambda_MBSREM = %f\n", w_vec.lambda_MBSREM);
		mexEvalString("pause(.0001);");
	}
	// Relaxation parameters
	if (MethodList.RAMLA || MethodList.BSREM)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.lambda_BSREM = (float*)mxGetSingles(mxGetField(options, 0, "lam"));
#else
		w_vec.lambda_BSREM = (float*)mxGetData(mxGetField(options, 0, "lam"));
#endif
	if (MethodList.ROSEM || MethodList.ROSEMMAP)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.lambda_ROSEM = (float*)mxGetSingles(mxGetField(options, 0, "lam_ROSEM"));
#else
		w_vec.lambda_ROSEM = (float*)mxGetData(mxGetField(options, 0, "lam_ROSEM"));
#endif
	if (MethodList.PKMA) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.lambda_PKMA = (float*)mxGetSingles(mxGetField(options, 0, "lam_PKMA"));
		w_vec.alpha_PKMA = (float*)mxGetSingles(mxGetField(options, 0, "alpha_PKMA"));
		w_vec.sigma_PKMA = (float*)mxGetSingles(mxGetField(options, 0, "sigma_PKMA"));
#else
		w_vec.lambda_PKMA = (float*)mxGetData(mxGetField(options, 0, "lam_PKMA"));
		w_vec.alpha_PKMA = (float*)mxGetData(mxGetField(options, 0, "alpha_PKMA"));
		w_vec.sigma_PKMA = (float*)mxGetData(mxGetField(options, 0, "sigma_PKMA"));
#endif
	}
	// Power factor for ACOSEM
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1)
		w_vec.h_ACOSEM = getScalarFloat(mxGetField(options, 0, "h"), -43);
	if (MethodList.TGV && MethodList.MAP) {
		data.TGVAlpha = getScalarFloat(mxGetField(options, 0, "alphaTGV"), -44);
		data.TGVBeta = getScalarFloat(mxGetField(options, 0, "betaTGV"), -45);
		data.NiterTGV = getScalarUInt32(mxGetField(options, 0, "NiterTGV"), -46);
	}
	if (MethodList.NLM && MethodList.MAP) {
		w_vec.NLM_anatomical = getScalarBool(mxGetField(options, 0, "NLM_use_anatomical"), -47);
		w_vec.NLTV = getScalarBool(mxGetField(options, 0, "NLTV"), -48);
		w_vec.NLM_MRP = getScalarBool(mxGetField(options, 0, "NLM_MRP"), -49);
		if (w_vec.NLM_anatomical)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			w_vec.NLM_ref = af::array(Nx, Ny, Nz, (float*)mxGetSingles(mxGetField(options, 0, "NLM_ref")), afHost);
#else
			w_vec.NLM_ref = af::array(Nx, Ny, Nz, (float*)mxGetData(mxGetField(options, 0, "NLM_ref")), afHost);
#endif
		w_vec.h2 = getScalarFloat(mxGetField(options, 0, "sigma"), -50);
		w_vec.h2 = w_vec.h2 * w_vec.h2;
		w_vec.Nlx = getScalarUInt32(mxGetField(options, 0, "Nlx"), -51);
		w_vec.Nly = getScalarUInt32(mxGetField(options, 0, "Nly"), -52);
		w_vec.Nlz = getScalarUInt32(mxGetField(options, 0, "Nlz"), -53);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), (float*)mxGetSingles(mxGetField(options, 0, "gaussianNLM")), afHost);
#else
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), (float*)mxGetData(mxGetField(options, 0, "gaussianNLM")), afHost);
#endif
	}
	if (MethodList.CUSTOM) {
		for (uint32_t kk = 0; kk < vec.imEstimates.size(); kk++) {
			const char* varTot = mxArrayToString(mxGetCell(mxGetField(options, 0, "varTot"), kk));
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			vec.imEstimates[kk](af::span, 0) = af::array(im_dim, (float*)mxGetSingles(mxGetField(mxGetField(options, 0, "im_vectors"), 0, varTot)), afHost);
#else
			vec.imEstimates[kk](af::span, 0) = af::array(im_dim, (float*)mxGetData(mxGetField(mxGetField(options, 0, "im_vectors"), 0, varTot)), afHost);
#endif
			if (DEBUG) {
				mexPrintf("%s\n", varTot);
				mexPrintf("vec.imEstimates[kk](af::span, 0) = %f\n", af::sum<float>(vec.imEstimates[kk](af::span, 0)));
				mexEvalString("pause(.0001);");
			}
		}
		tt = 0;
		for (uint32_t ll = 0; ll < w_vec.nMAP; ll++) {
			const char* varApu = mxArrayToString(mxGetCell(mxGetField(options, 0, "varApu"), tt));
			const char* varBeta = mxArrayToString(mxGetCell(mxGetField(options, 0, "varBeta"), tt));
			const char* varGrad = mxArrayToString(mxGetCell(mxGetField(options, 0, "varGrad"), tt));
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			vec.imEstimates.push_back(af::array(im_dim, (float*)mxGetSingles(mxGetField(options, 0, varApu)), afHost));
			w_vec.dU.push_back(af::array(im_dim, (float*)mxGetSingles(mxGetField(options, 0, varGrad)), afHost));
#else
			vec.imEstimates.push_back(af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, varApu)), afHost));
			w_vec.dU.push_back(af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, varGrad)), afHost));
#endif
			beta.push_back(getScalarFloat(mxGetField(options, 0, varBeta), -54));
			tt++;
		}
		if (MethodList.MBSREM) {
			if (iter0 > 0 || osa_iter0 > 0) {
				w_vec.U = getScalarFloat(mxGetField(options, 0, "U"), -55);
				w_vec.epsilon_mramla = getScalarFloat(mxGetField(options, 0, "epsilon_mramla"), -56);
			}
		}
		if (MethodList.OSLCOSEM > 0) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			vec.C_osl = af::array(im_dim, subsets, (float*)mxGetSingles(mxGetField(options, 0, "C_osl")), afHost);
#else
			vec.C_osl = af::array(im_dim, subsets, (float*)mxGetData(mxGetField(options, 0, "C_osl")), afHost);
#endif
		}
		if ((MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA))
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			w_vec.D = af::array(im_dim, (float*)mxGetSingles(mxGetField(options, 0, "D")), afHost);
#else
			w_vec.D = af::array(im_dim, (float*)mxGetData(mxGetField(options, 0, "D")), afHost);
#endif

		//uint32_t dd = n_rekos_mlem + n_rekos - 1U - nMAPOS;
		//vec.im_mlem(af::seq(n_rekos_mlem * im_dim - im_dim, n_rekos_mlem * im_dim - 1u)) = vec.imEstimates[dd];
		uint32_t dd = 0U;
		uint32_t yy = 0U;
		for (int kk = 0; kk < (w_vec.nMLEM + w_vec.nMAPML * w_vec.nPriorsTot); kk++) {
			vec.im_mlem(af::seq(yy, yy + im_dim - 1u)) = vec.imEstimates[dd];
			yy += im_dim;
			dd++;
		}
		//uint32_t dd = n_rekos_mlem + n_rekos - 1U;
		//uint32_t yy = n_rekos * im_dim;
		dd = 0U;
		yy = 0U;
		if (DEBUG) {
			mexPrintf("vec.imEstimates.size() = %d\n", vec.imEstimates.size());
			mexPrintf("[dd] = %u\n", dd);
			mexPrintf("[yy] = %u\n", yy);
			mexPrintf("vec.imEstimates[dd].dims(0) = %d\n", vec.imEstimates[dd].dims(0));
			mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
			mexEvalString("pause(.0001);");
		}
		for (int kk = 0; kk < (w_vec.nPriorsTot * w_vec.nMAPOS + w_vec.nOS); kk++) {
			vec.im_os(af::seq(yy, yy + im_dim - 1u)) = vec.imEstimates[dd];
			yy += im_dim;
			dd++;
		}
	}
	if (use_psf) {
		w_vec.g_dim_x = getScalarUInt32(mxGetField(options, 0, "g_dim_x"), -57);
		w_vec.g_dim_y = getScalarUInt32(mxGetField(options, 0, "g_dim_y"), -58);
		w_vec.g_dim_z = getScalarUInt32(mxGetField(options, 0, "g_dim_z"), -59);
		w_vec.deconvolution = getScalarBool(mxGetField(options, 0, "deblurring"), -60);
	}
}

// Obtain the reconstruction methods used
void get_rec_methods(const mxArray * options, RecMethods &MethodList)
{
	// Non-MAP/prior algorithms
	MethodList.MLEM = getScalarBool(mxGetField(options, 0, "MLEM"), -61);
	MethodList.OSEM = getScalarBool(mxGetField(options, 0, "OSEM"), -61);
	MethodList.RAMLA = getScalarBool(mxGetField(options, 0, "RAMLA"), -61);
	MethodList.MRAMLA = getScalarBool(mxGetField(options, 0, "MRAMLA"), -61);
	MethodList.ROSEM = getScalarBool(mxGetField(options, 0, "ROSEM"), -61);
	MethodList.RBI = getScalarBool(mxGetField(options, 0, "RBI"), -61);
	MethodList.DRAMA = getScalarBool(mxGetField(options, 0, "DRAMA"), -61);
	MethodList.COSEM = getScalarBool(mxGetField(options, 0, "COSEM"), -61);
	MethodList.ECOSEM = getScalarBool(mxGetField(options, 0, "ECOSEM"), -61);
	MethodList.ACOSEM = getScalarBool(mxGetField(options, 0, "ACOSEM"), -61);

	// Priors
	MethodList.MRP = getScalarBool(mxGetField(options, 0, "MRP"), -61);
	MethodList.Quad = getScalarBool(mxGetField(options, 0, "quad"), -61);
	MethodList.Huber = getScalarBool(mxGetField(options, 0, "Huber"), -61);
	MethodList.L = getScalarBool(mxGetField(options, 0, "L"), -61);
	MethodList.FMH = getScalarBool(mxGetField(options, 0, "FMH"), -61);
	MethodList.WeightedMean = getScalarBool(mxGetField(options, 0, "weighted_mean"), -61);
	MethodList.TV = getScalarBool(mxGetField(options, 0, "TV"), -61);
	MethodList.AD = getScalarBool(mxGetField(options, 0, "AD"), -61);
	MethodList.APLS = getScalarBool(mxGetField(options, 0, "APLS"), -61);
	MethodList.TGV = getScalarBool(mxGetField(options, 0, "TGV"), -61);
	MethodList.NLM = getScalarBool(mxGetField(options, 0, "NLM"), -61);
	MethodList.RDP = getScalarBool(mxGetField(options, 0, "RDP"), -61);

	// MAP/prior-based algorithms
	MethodList.OSLMLEM = getScalarBool(mxGetField(options, 0, "OSL_MLEM"), -61);
	MethodList.OSLOSEM = getScalarBool(mxGetField(options, 0, "OSL_OSEM"), -61);
	MethodList.BSREM = getScalarBool(mxGetField(options, 0, "BSREM"), -61);
	MethodList.MBSREM = getScalarBool(mxGetField(options, 0, "MBSREM"), -61);
	MethodList.ROSEMMAP = getScalarBool(mxGetField(options, 0, "ROSEM_MAP"), -61);
	MethodList.RBIOSL = getScalarBool(mxGetField(options, 0, "OSL_RBI"), -61);
	MethodList.OSLCOSEM = getScalarUInt32(mxGetField(options, 0, "OSL_COSEM"), -61);
	MethodList.PKMA = getScalarBool(mxGetField(options, 0, "PKMA"), -61);

	// Whether MAP/prior-based algorithms are used
	MethodList.MAP = getScalarBool(mxGetField(options, 0, "MAP"), -61);

	// Custom prior
	MethodList.CUSTOM = getScalarBool(mxGetField(options, 0, "custom"), -61);
}

// Create the MATLAB output
// Creates the mxArrays that are later placed in the cell array
// Creates the array pointers to the mxArrays
//void create_matlab_output(matlabArrays &ArrayList, const mwSize *dimmi, const RecMethods &MethodList, const uint32_t dim_n)
//{
//	//mwSize dimu[2] = { 128 * 128 * 63, 9 };
//	// Output arrays for the MATLAB cell array
//	if (MethodList.MLEM)
//		ArrayList.mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.OSEM)
//		ArrayList.osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRAMLA)
//		ArrayList.ramlaM = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.ramlaM = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.RAMLA)
//		ArrayList.ramla = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.ramla = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.ROSEM)
//		ArrayList.rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.RBI)
//		ArrayList.rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.DRAMA)
//		ArrayList.drama = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.drama = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.COSEM)
//		ArrayList.cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.ECOSEM)
//		ArrayList.ecosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.ecosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.ACOSEM)
//		ArrayList.acosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.acosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.MRP && MethodList.OSLOSEM)
//		ArrayList.mrp_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.OSLMLEM)
//		ArrayList.mrp_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.BSREM)
//		ArrayList.mrp_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.MBSREM)
//		ArrayList.mrp_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.ROSEMMAP)
//		ArrayList.mrp_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.RBIOSL)
//		ArrayList.mrp_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
//		ArrayList.mrp_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.mrp_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.Quad && MethodList.OSLOSEM)
//		ArrayList.quad_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.OSLMLEM)
//		ArrayList.quad_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.BSREM)
//		ArrayList.quad_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.MBSREM)
//		ArrayList.quad_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.ROSEMMAP)
//		ArrayList.quad_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.RBIOSL)
//		ArrayList.quad_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
//		ArrayList.quad_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.quad_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.Huber && MethodList.OSLOSEM)
//		ArrayList.Huber_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.OSLMLEM)
//		ArrayList.Huber_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.BSREM)
//		ArrayList.Huber_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.MBSREM)
//		ArrayList.Huber_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.ROSEMMAP)
//		ArrayList.Huber_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.RBIOSL)
//		ArrayList.Huber_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
//		ArrayList.Huber_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.Huber_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.L && MethodList.OSLOSEM)
//		ArrayList.L_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.OSLMLEM)
//		ArrayList.L_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.BSREM)
//		ArrayList.L_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.MBSREM)
//		ArrayList.L_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.ROSEMMAP)
//		ArrayList.L_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.RBIOSL)
//		ArrayList.L_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.L && MethodList.OSLCOSEM > 0)
//		ArrayList.L_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.L_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.FMH && MethodList.OSLOSEM)
//		ArrayList.fmh_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.OSLMLEM)
//		ArrayList.fmh_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.BSREM)
//		ArrayList.fmh_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.MBSREM)
//		ArrayList.fmh_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.ROSEMMAP)
//		ArrayList.fmh_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.RBIOSL)
//		ArrayList.fmh_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
//		ArrayList.fmh_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.fmh_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.WeightedMean && MethodList.OSLOSEM)
//		ArrayList.weighted_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.OSLMLEM)
//		ArrayList.weighted_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.BSREM)
//		ArrayList.weighted_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.MBSREM)
//		ArrayList.weighted_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
//		ArrayList.weighted_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.RBIOSL)
//		ArrayList.weighted_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
//		ArrayList.weighted_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.weighted_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.TV && MethodList.OSLOSEM)
//		ArrayList.TV_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.OSLMLEM)
//		ArrayList.TV_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.BSREM)
//		ArrayList.TV_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.MBSREM)
//		ArrayList.TV_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.ROSEMMAP)
//		ArrayList.TV_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.RBIOSL)
//		ArrayList.TV_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.TV && MethodList.OSLCOSEM > 0)
//		ArrayList.TV_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TV_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.AD && MethodList.OSLOSEM)
//		ArrayList.AD_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.OSLMLEM)
//		ArrayList.AD_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.BSREM)
//		ArrayList.AD_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.MBSREM)
//		ArrayList.AD_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.ROSEMMAP)
//		ArrayList.AD_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.RBIOSL)
//		ArrayList.AD_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.AD && MethodList.OSLCOSEM > 0)
//		ArrayList.AD_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.AD_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.APLS && MethodList.OSLOSEM)
//		ArrayList.APLS_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.OSLMLEM)
//		ArrayList.APLS_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.BSREM)
//		ArrayList.APLS_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.MBSREM)
//		ArrayList.APLS_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.ROSEMMAP)
//		ArrayList.APLS_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.RBIOSL)
//		ArrayList.APLS_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
//		ArrayList.APLS_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.APLS_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.OSLOSEM)
//		ArrayList.TGV_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.OSLMLEM)
//		ArrayList.TGV_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.BSREM)
//		ArrayList.TGV_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MBSREM)
//		ArrayList.TGV_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.ROSEMMAP)
//		ArrayList.TGV_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.RBIOSL)
//		ArrayList.TGV_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.OSLCOSEM > 0)
//		ArrayList.TGV_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.TGV_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//
//	if (MethodList.OSLOSEM)
//		ArrayList.NLM_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.OSLMLEM)
//		ArrayList.NLM_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.BSREM)
//		ArrayList.NLM_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.MBSREM)
//		ArrayList.NLM_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.ROSEMMAP)
//		ArrayList.NLM_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.RBIOSL)
//		ArrayList.NLM_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.OSLCOSEM > 0)
//		ArrayList.NLM_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//	else
//		ArrayList.NLM_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	if (MethodList.CUSTOM) {
//		if (MethodList.OSLOSEM) {
//			ArrayList.custom_osem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		}
//		else
//			ArrayList.custom_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.OSLMLEM)
//			ArrayList.custom_mlem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		else
//			ArrayList.custom_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.BSREM)
//			ArrayList.custom_bsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		else
//			ArrayList.custom_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.MBSREM)
//			ArrayList.custom_mbsrem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		else
//			ArrayList.custom_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.ROSEMMAP)
//			ArrayList.custom_rosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		else
//			ArrayList.custom_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.RBIOSL)
//			ArrayList.custom_rbi = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		else
//			ArrayList.custom_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		if (MethodList.OSLCOSEM > 0) {
//			ArrayList.custom_cosem = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//			ArrayList.c_osl_custom = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
//		}
//		else {
//			ArrayList.custom_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//			ArrayList.c_osl_custom = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//		}
//		if (MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI) {
//			mwSize dimmiD[3] = { static_cast<mwSize>(0), static_cast<mwSize>(0), static_cast<mwSize>(0)};
//			for (int ii = 0; ii < 3; ii++)
//				dimmiD[ii] = dimmi[ii];
//			ArrayList.D_custom = mxCreateNumericArray(3, dimmiD, mxSINGLE_CLASS, mxREAL);
//		}
//		else
//			ArrayList.D_custom = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
//	}
//
//
//	if (MethodList.MLEM)
//		ArrayList.ele_ml = (float*)mxGetData(ArrayList.mlem);
//	if (MethodList.OSEM)
//		ArrayList.ele_os = (float*)mxGetData(ArrayList.osem);
//	if (MethodList.MRAMLA)
//		ArrayList.ele_ramlaM = (float*)mxGetData(ArrayList.ramlaM);
//	if (MethodList.RAMLA)
//		ArrayList.ele_ramla = (float*)mxGetData(ArrayList.ramla);
//	if (MethodList.ROSEM)
//		ArrayList.ele_rosem = (float*)mxGetData(ArrayList.rosem);
//	if (MethodList.RBI)
//		ArrayList.ele_rbi = (float*)mxGetData(ArrayList.rbi);
//	if (MethodList.DRAMA)
//		ArrayList.ele_drama = (float*)mxGetData(ArrayList.drama);
//	if (MethodList.COSEM)
//		ArrayList.ele_cosem = (float*)mxGetData(ArrayList.cosem);
//	if (MethodList.ECOSEM)
//		ArrayList.ele_ecosem = (float*)mxGetData(ArrayList.ecosem);
//	if (MethodList.ACOSEM)
//		ArrayList.ele_acosem = (float*)mxGetData(ArrayList.acosem);
//
//	if (MethodList.OSLOSEM && MethodList.MRP)
//		ArrayList.ele_mrp_osem = (float*)mxGetData(ArrayList.mrp_osem);
//	if (MethodList.OSLMLEM && MethodList.MRP)
//		ArrayList.ele_mrp_mlem = (float*)mxGetData(ArrayList.mrp_mlem);
//	if (MethodList.MRP && MethodList.BSREM)
//		ArrayList.ele_mrp_bsrem = (float*)mxGetData(ArrayList.mrp_bsrem);
//	if (MethodList.MRP && MethodList.MBSREM)
//		ArrayList.ele_mrp_mbsrem = (float*)mxGetData(ArrayList.mrp_mbsrem);
//	if (MethodList.MRP && MethodList.ROSEMMAP)
//		ArrayList.ele_mrp_rosem = (float*)mxGetData(ArrayList.mrp_rosem);
//	if (MethodList.MRP && MethodList.RBIOSL)
//		ArrayList.ele_mrp_rbi = (float*)mxGetData(ArrayList.mrp_rbi);
//	if (MethodList.MRP && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_mrp_cosem = (float*)mxGetData(ArrayList.mrp_cosem);
//
//	if (MethodList.Quad && MethodList.OSLOSEM)
//		ArrayList.ele_quad_osem = (float*)mxGetData(ArrayList.quad_osem);
//	if (MethodList.OSLMLEM && MethodList.Quad)
//		ArrayList.ele_quad_mlem = (float*)mxGetData(ArrayList.quad_mlem);
//	if (MethodList.Quad && MethodList.BSREM)
//		ArrayList.ele_quad_bsrem = (float*)mxGetData(ArrayList.quad_bsrem);
//	if (MethodList.Quad && MethodList.MBSREM)
//		ArrayList.ele_quad_mbsrem = (float*)mxGetData(ArrayList.quad_mbsrem);
//	if (MethodList.Quad && MethodList.ROSEMMAP)
//		ArrayList.ele_quad_rosem = (float*)mxGetData(ArrayList.quad_rosem);
//	if (MethodList.Quad && MethodList.RBIOSL)
//		ArrayList.ele_quad_rbi = (float*)mxGetData(ArrayList.quad_rbi);
//	if (MethodList.Quad && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_quad_cosem = (float*)mxGetData(ArrayList.quad_cosem);
//
//	if (MethodList.Huber && MethodList.OSLOSEM)
//		ArrayList.ele_Huber_osem = (float*)mxGetData(ArrayList.Huber_osem);
//	if (MethodList.OSLMLEM && MethodList.Huber)
//		ArrayList.ele_Huber_mlem = (float*)mxGetData(ArrayList.Huber_mlem);
//	if (MethodList.Huber && MethodList.BSREM)
//		ArrayList.ele_Huber_bsrem = (float*)mxGetData(ArrayList.Huber_bsrem);
//	if (MethodList.Huber && MethodList.MBSREM)
//		ArrayList.ele_Huber_mbsrem = (float*)mxGetData(ArrayList.Huber_mbsrem);
//	if (MethodList.Huber && MethodList.ROSEMMAP)
//		ArrayList.ele_Huber_rosem = (float*)mxGetData(ArrayList.Huber_rosem);
//	if (MethodList.Huber && MethodList.RBIOSL)
//		ArrayList.ele_Huber_rbi = (float*)mxGetData(ArrayList.Huber_rbi);
//	if (MethodList.Huber && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_Huber_cosem = (float*)mxGetData(ArrayList.Huber_cosem);
//
//	if (MethodList.L && MethodList.OSLOSEM)
//		ArrayList.ele_L_osem = (float*)mxGetData(ArrayList.L_osem);
//	if (MethodList.OSLMLEM && MethodList.L)
//		ArrayList.ele_L_mlem = (float*)mxGetData(ArrayList.L_mlem);
//	if (MethodList.L && MethodList.BSREM)
//		ArrayList.ele_L_bsrem = (float*)mxGetData(ArrayList.L_bsrem);
//	if (MethodList.L && MethodList.MBSREM)
//		ArrayList.ele_L_mbsrem = (float*)mxGetData(ArrayList.L_mbsrem);
//	if (MethodList.L && MethodList.ROSEMMAP)
//		ArrayList.ele_L_rosem = (float*)mxGetData(ArrayList.L_rosem);
//	if (MethodList.L && MethodList.RBIOSL)
//		ArrayList.ele_L_rbi = (float*)mxGetData(ArrayList.L_rbi);
//	if (MethodList.L && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_L_cosem = (float*)mxGetData(ArrayList.L_cosem);
//
//	if (MethodList.FMH && MethodList.OSLOSEM)
//		ArrayList.ele_fmh_osem = (float*)mxGetData(ArrayList.fmh_osem);
//	if (MethodList.OSLMLEM && MethodList.FMH)
//		ArrayList.ele_fmh_mlem = (float*)mxGetData(ArrayList.fmh_mlem);
//	if (MethodList.FMH && MethodList.BSREM)
//		ArrayList.ele_fmh_bsrem = (float*)mxGetData(ArrayList.fmh_bsrem);
//	if (MethodList.FMH && MethodList.MBSREM)
//		ArrayList.ele_fmh_mbsrem = (float*)mxGetData(ArrayList.fmh_mbsrem);
//	if (MethodList.FMH && MethodList.ROSEMMAP)
//		ArrayList.ele_fmh_rosem = (float*)mxGetData(ArrayList.fmh_rosem);
//	if (MethodList.FMH && MethodList.RBIOSL)
//		ArrayList.ele_fmh_rbi = (float*)mxGetData(ArrayList.fmh_rbi);
//	if (MethodList.FMH && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_fmh_cosem = (float*)mxGetData(ArrayList.fmh_cosem);
//
//	if (MethodList.WeightedMean && MethodList.OSLOSEM)
//		ArrayList.ele_weighted_osem = (float*)mxGetData(ArrayList.weighted_osem);
//	if (MethodList.OSLMLEM && MethodList.WeightedMean)
//		ArrayList.ele_weighted_mlem = (float*)mxGetData(ArrayList.weighted_mlem);
//	if (MethodList.WeightedMean && MethodList.BSREM)
//		ArrayList.ele_weighted_bsrem = (float*)mxGetData(ArrayList.weighted_bsrem);
//	if (MethodList.WeightedMean && MethodList.MBSREM)
//		ArrayList.ele_weighted_mbsrem = (float*)mxGetData(ArrayList.weighted_mbsrem);
//	if (MethodList.WeightedMean && MethodList.ROSEMMAP)
//		ArrayList.ele_weighted_rosem = (float*)mxGetData(ArrayList.weighted_rosem);
//	if (MethodList.WeightedMean && MethodList.RBIOSL)
//		ArrayList.ele_weighted_rbi = (float*)mxGetData(ArrayList.weighted_rbi);
//	if (MethodList.WeightedMean && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_weighted_cosem = (float*)mxGetData(ArrayList.weighted_cosem);
//
//	if (MethodList.TV && MethodList.OSLOSEM)
//		ArrayList.ele_TV_osem = (float*)mxGetData(ArrayList.TV_osem);
//	if (MethodList.OSLMLEM && MethodList.TV)
//		ArrayList.ele_TV_mlem = (float*)mxGetData(ArrayList.TV_mlem);
//	if (MethodList.TV && MethodList.BSREM)
//		ArrayList.ele_TV_bsrem = (float*)mxGetData(ArrayList.TV_bsrem);
//	if (MethodList.TV && MethodList.MBSREM)
//		ArrayList.ele_TV_mbsrem = (float*)mxGetData(ArrayList.TV_mbsrem);
//	if (MethodList.TV && MethodList.ROSEMMAP)
//		ArrayList.ele_TV_rosem = (float*)mxGetData(ArrayList.TV_rosem);
//	if (MethodList.TV && MethodList.RBIOSL)
//		ArrayList.ele_TV_rbi = (float*)mxGetData(ArrayList.TV_rbi);
//	if (MethodList.TV && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_TV_cosem = (float*)mxGetData(ArrayList.TV_cosem);
//
//	if (MethodList.AD && MethodList.OSLOSEM)
//		ArrayList.ele_AD_osem = (float*)mxGetData(ArrayList.AD_osem);
//	if (MethodList.OSLMLEM && MethodList.AD)
//		ArrayList.ele_AD_mlem = (float*)mxGetData(ArrayList.AD_mlem);
//	if (MethodList.AD && MethodList.BSREM)
//		ArrayList.ele_AD_bsrem = (float*)mxGetData(ArrayList.AD_bsrem);
//	if (MethodList.AD && MethodList.MBSREM)
//		ArrayList.ele_AD_mbsrem = (float*)mxGetData(ArrayList.AD_mbsrem);
//	if (MethodList.AD && MethodList.ROSEMMAP)
//		ArrayList.ele_AD_rosem = (float*)mxGetData(ArrayList.AD_rosem);
//	if (MethodList.AD && MethodList.RBIOSL)
//		ArrayList.ele_AD_rbi = (float*)mxGetData(ArrayList.AD_rbi);
//	if (MethodList.AD && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_AD_cosem = (float*)mxGetData(ArrayList.AD_cosem);
//
//	if (MethodList.APLS && MethodList.OSLOSEM)
//		ArrayList.ele_APLS_osem = (float*)mxGetData(ArrayList.APLS_osem);
//	if (MethodList.OSLMLEM && MethodList.APLS)
//		ArrayList.ele_APLS_mlem = (float*)mxGetData(ArrayList.APLS_mlem);
//	if (MethodList.APLS && MethodList.BSREM)
//		ArrayList.ele_APLS_bsrem = (float*)mxGetData(ArrayList.APLS_bsrem);
//	if (MethodList.APLS && MethodList.MBSREM)
//		ArrayList.ele_APLS_mbsrem = (float*)mxGetData(ArrayList.APLS_mbsrem);
//	if (MethodList.APLS && MethodList.ROSEMMAP)
//		ArrayList.ele_APLS_rosem = (float*)mxGetData(ArrayList.APLS_rosem);
//	if (MethodList.APLS && MethodList.RBIOSL)
//		ArrayList.ele_APLS_rbi = (float*)mxGetData(ArrayList.APLS_rbi);
//	if (MethodList.APLS && MethodList.OSLCOSEM > 0)
//		ArrayList.ele_APLS_cosem = (float*)mxGetData(ArrayList.APLS_cosem);
//
//	if (MethodList.OSLOSEM)
//		ArrayList.ele_TGV_osem = (float*)mxGetData(ArrayList.TGV_osem);
//	if (MethodList.OSLMLEM && MethodList.TGV)
//		ArrayList.ele_TGV_mlem = (float*)mxGetData(ArrayList.TGV_mlem);
//	if (MethodList.BSREM)
//		ArrayList.ele_TGV_bsrem = (float*)mxGetData(ArrayList.TGV_bsrem);
//	if (MethodList.MBSREM)
//		ArrayList.ele_TGV_mbsrem = (float*)mxGetData(ArrayList.TGV_mbsrem);
//	if (MethodList.ROSEMMAP)
//		ArrayList.ele_TGV_rosem = (float*)mxGetData(ArrayList.TGV_rosem);
//	if (MethodList.RBIOSL)
//		ArrayList.ele_TGV_rbi = (float*)mxGetData(ArrayList.TGV_rbi);
//	if (MethodList.OSLCOSEM > 0)
//		ArrayList.ele_TGV_cosem = (float*)mxGetData(ArrayList.TGV_cosem);
//
//	if (MethodList.OSLOSEM)
//		ArrayList.ele_NLM_osem = (float*)mxGetData(ArrayList.NLM_osem);
//	if (MethodList.OSLMLEM && MethodList.NLM)
//		ArrayList.ele_NLM_mlem = (float*)mxGetData(ArrayList.NLM_mlem);
//	if (MethodList.BSREM)
//		ArrayList.ele_NLM_bsrem = (float*)mxGetData(ArrayList.NLM_bsrem);
//	if (MethodList.MBSREM)
//		ArrayList.ele_NLM_mbsrem = (float*)mxGetData(ArrayList.NLM_mbsrem);
//	if (MethodList.ROSEMMAP)
//		ArrayList.ele_NLM_rosem = (float*)mxGetData(ArrayList.NLM_rosem);
//	if (MethodList.RBIOSL)
//		ArrayList.ele_NLM_rbi = (float*)mxGetData(ArrayList.NLM_rbi);
//	if (MethodList.OSLCOSEM > 0)
//		ArrayList.ele_NLM_cosem = (float*)mxGetData(ArrayList.NLM_cosem);
//
//	if (MethodList.CUSTOM) {
//		if (MethodList.OSLOSEM) {
//			ArrayList.ele_custom_osem = (float*)mxGetData(ArrayList.custom_osem);
//		}
//		if (MethodList.OSLMLEM)
//			ArrayList.ele_custom_mlem = (float*)mxGetData(ArrayList.custom_mlem);
//		if (MethodList.BSREM)
//			ArrayList.ele_custom_bsrem = (float*)mxGetData(ArrayList.custom_bsrem);
//		if (MethodList.MBSREM)
//			ArrayList.ele_custom_mbsrem = (float*)mxGetData(ArrayList.custom_mbsrem);
//		if (MethodList.ROSEMMAP)
//			ArrayList.ele_custom_rosem = (float*)mxGetData(ArrayList.custom_rosem);
//		if (MethodList.RBIOSL)
//			ArrayList.ele_custom_rbi = (float*)mxGetData(ArrayList.custom_rbi);
//		if (MethodList.OSLCOSEM > 0) {
//			ArrayList.ele_custom_cosem = (float*)mxGetData(ArrayList.custom_cosem);
//			ArrayList.ele_c_osl_custom = (float*)mxGetData(ArrayList.c_osl_custom);
//		}
//		if (MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI) {
//			ArrayList.ele_D_custom = (float*)mxGetData(ArrayList.D_custom);
//		}
//	}
//
//
//}

// Transfers the device data to host
// First transfer the ArrayFire arrays from the device to the host pointers pointing to the mxArrays
// Transfer the mxArrays to the cell
void device_to_host_cell(const RecMethods &MethodList, AF_im_vectors & vec, uint32_t & oo, mxArray * cell, Weighting & w_vec,
	const mwSize* dimmi, const uint32_t dim_n)
{
	uint32_t kk;
	// Transfer data back to host
	for (kk = 0; kk < w_vec.nTot; kk++) {
		mxArray* apu = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		float* apuF = (float*)mxGetSingles(apu);
#else
		float* apuF = (float*)mxGetData(apu);
#endif
		vec.imEstimates[kk].host(apuF);
		af::sync();
		mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), mxDuplicateArray(apu));
	}
	if (MethodList.CUSTOM) {
		kk = w_vec.nTot;
		if (MethodList.OSLCOSEM > 0) {
			mxArray* apu = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apuF = (float*)mxGetSingles(apu);
#else
			float* apuF = (float*)mxGetData(apu);
#endif
			vec.C_osl.host(apuF);
			af::sync();
			mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), mxDuplicateArray(apu));
		}
		kk++;
		if (MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBI || MethodList.RBIOSL || MethodList.PKMA) {
			mwSize dimmiD[3] = { static_cast<mwSize>(0), static_cast<mwSize>(0), static_cast<mwSize>(0) };
			for (int ii = 0; ii < 3; ii++)
				dimmiD[ii] = dimmi[ii];
			mxArray* apu = mxCreateNumericArray(3, dimmiD, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apuF = (float*)mxGetSingles(apu);
#else
			float* apuF = (float*)mxGetData(apu);
#endif
			w_vec.D.host(apuF);
			af::sync();
			mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), mxDuplicateArray(apu));
		}
		kk++;
		if (MethodList.MBSREM) {
			mxArray* epsilon = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* epsilon1 = (float*)mxGetSingles(epsilon);
#else
			float* epsilon1 = (float*)mxGetData(epsilon);
#endif
			epsilon1[0] = w_vec.epsilon_mramla;
			mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), (epsilon));
		}
		kk++;
		if (MethodList.MBSREM) {
			mxArray* UU = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* U = (float*)mxGetSingles(UU);
#else
			float* U = (float*)mxGetData(UU);
#endif
			U[0] = w_vec.U;
			mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), (UU));
		}
	}
}

// Compute the epsilon value for MBSREM and MRAMLA
float MBSREM_epsilon(const af::array & Sino, const float epps, const uint32_t randoms_correction, const af::array& rand, const af::array& D, 
	const bool TOF, const int64_t nBins, const bool CT)
{
	float eps;
	if (CT) {
		af::array hk_summa = -af::exp(-Sino) / Sino - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (randoms_correction == 1u) {
			Iind = (Sino > 0.f & rand == 0.f);
			if (af::sum<float>(Iind) == 0.f)
				return 1e8f;
			P_Sino = Sino(Iind);
			apu = D + rand;
			apu = af::sum(-af::exp(-apu) / Sino - apu);
		}
		else {
			Iind = (Sino > 0.f);
			P_Sino = Sino(Iind);
			apu = af::sum(-af::exp(-D) / Sino - D);
		}
		if (randoms_correction == 1u)
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		else
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		//if (DEBUG) {
		//	mexPrintf("exp = %f\n", af::min<float>(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		//	mexEvalString("pause(.0001);");
		//}
		af::array epsilon = (af::min)(P_Sino, af::log(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
	else {
		af::array hk_summa = Sino * af::log(Sino) - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (TOF && randoms_correction) {
			af::array rInd = rand == 0.f;
			P_Sino = Sino(Sino > 0.f & af::tile(rInd, nBins));
			apu = D + af::tile(rand, nBins);
			apu = af::sum(Sino * af::log(apu) - apu);
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Sino > 0.f & af::tile(rInd, nBins)), batchMinus);
		}
		else {
			if (randoms_correction == 1u) {
				Iind = (Sino > 0.f & rand == 0.f);
				if (af::sum<float>(Iind) == 0.f)
					return 1e8f;
				P_Sino = Sino(Iind);
				apu = D + rand;
				apu = af::sum(Sino * af::log(apu) - apu);
			}
			else {
				Iind = (Sino > 0.f);
				P_Sino = Sino(Iind);
				apu = af::sum(Sino * af::log(D) - D);
			}
			if (randoms_correction == 1u)
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
			else
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		}
		//if (DEBUG) {
		//	mexPrintf("exp = %f\n", af::min<float>(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		//	mexEvalString("pause(.0001);");
		//}
		af::array epsilon = (af::min)(P_Sino, af::exp(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
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
		af::array out = padd;
		if (Ndx > 0)
			out = af::join(0, af::flip(padd(af::seq(static_cast<double>(Ndx)), af::span, af::span), 0), padd, af::flip(padd(af::seq(static_cast<double>(padd.dims(0) - Ndx), static_cast<double>(padd.dims(0) - 1)), af::span, af::span), 0));
		if (Ndy > 0)
			out = af::join(1, af::flip(out(af::span, af::seq(static_cast<double>(Ndy)), af::span), 1), out, af::flip(out(af::span, af::seq(static_cast<double>(out.dims(1) - Ndy), static_cast<double>(out.dims(1) - 1)), af::span), 1));
		if (Nz == 1 || Ndz == 0) {
		}
		else {
			out = af::join(2, af::flip(out(af::span, af::span, af::seq(static_cast<double>(Ndz))), 2), out, af::flip(out(af::span, af::span, af::seq(static_cast<double>(out.dims(2) - Ndz), static_cast<double>(out.dims(2) - 1))), 2));
		}
		padd = out;
	}
	return padd;
}

af::array EM(const af::array &im, const af::array &Summ, const af::array &rhs)
{
	return (im / Summ * rhs);
}

af::array OSL(const af::array &Summ, const af::array &dU, const float beta, const float epps)
{
	af::array output = (Summ + beta * dU + epps);
	//output(output < epps) = epps;
	//if (DEBUG) {
	//	mexPrintf("beta = %f\n", beta);
	//}
	return output;
}

af::array MBSREM(const af::array & im, const af::array & rhs, const float U, const af::array & pj3, const float* lam, const uint32_t iter, const uint32_t im_dim,
	const float beta, const af::array &dU, const af::array & Summ, const float epps)
{
	//af::array UU = af::constant(0.f, im_dim);
	af::array output;
	const af::array pp = im < (U / 2.f);
	//UU(pp) = im(pp);
	af::array UU = im / pj3;
	UU(!pp) = (U - im(!pp)) / (pj3(!pp));
	if (beta == 0.f)
		output = im + lam[iter] * UU * (rhs - Summ);
	else
		output = im + lam[iter] * UU * (rhs - beta * dU - Summ);
	output(output < epps) = epps;
	output(output >= U) = U - epps;
	//if (DEBUG) {
	//	mexPrintf("U = %f\n", U);
	//	mexPrintf("lam[iter] = %f\n", lam[iter]);
	//	mexPrintf("beta = %f\n", beta);
	//	mexPrintf("Summ = %f\n", af::sum<float>(Summ));
	//	mexPrintf("output = %f\n", af::sum<float>(lam[iter] * UU * rhs));
	//	mexEvalString("pause(.0001);");
	//}
	return output;
}

af::array BSREM(const af::array & im, const af::array & rhs, const float * lam, const uint32_t iter, const af::array& Summ)
{
	return ((1.f - lam[iter] * Summ) * im + lam[iter] * im * rhs);
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
	return (im + lam[iter] * im / Summ * (rhs - Summ));
}

af::array RBI(const af::array & im, const af::array & Summ, const af::array & rhs, const af::array& D, const float beta, const af::array& dU)
{
	af::array output = im;
	if (beta == 0.f) {
		const float Summa = 1.f / af::max<float>(Summ / D);
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
		output = af::pow(af::sum(C_co, 1) / D, h);
	}
	else {
		output = (af::sum(C_co, 1) / D);
	}
	return output;
}

af::array PKMA(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const float* alpha, const float* sigma, const af::array& D, 
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const float epps, const float beta, const af::array& dU)
{
	const af::array S = (im + epps)  / D;
	const af::array im_ = im;
	af::array im_apu = im - lam[iter] * S * (Summ - rhs + beta * dU);
	im_apu(im_apu < epps) = epps;
	uint32_t ind = iter * subsets + osa_iter;
	im_apu = (1.f - alpha[ind]) * im_ + alpha[ind] * (sigma[ind] * im_apu);
	return im_apu;
}
af::array MRP(const af::array& im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps,
	const af::array& offsets, const bool med_no_norm, const uint32_t im_dim, const kernelStruct& OpenCLStruct)
{
	af::array padd = padding(im, Nx, Ny, Nz, medx, medy, medz);
#ifdef OPENCL
	cl::Kernel kernelMed(OpenCLStruct.kernelMed);
	uint32_t kernelIndMed = 0U;
	const af::dim4 dimmi(padd.dims(0), padd.dims(1), padd.dims(2));
	cl::NDRange global_size(padd.dims(0), padd.dims(1), padd.dims(2));
	//cl::NDRange local_size(64U);
	af::array grad = af::constant(0.f, dimmi);
	padd = af::flat(padd);
	grad = af::flat(grad);
	cl::Buffer d_grad = cl::Buffer(*grad.device<cl_mem>(), true);
	cl::Buffer d_padd = cl::Buffer(*padd.device<cl_mem>(), true);
	if (DEBUG) {
		mexPrintf("padd = %f\n", af::sum<float>(padd));
	}
	af::sync();
	(*OpenCLStruct.af_queue).finish();
	kernelMed.setArg(kernelIndMed++, d_padd);
	kernelMed.setArg(kernelIndMed++, d_grad);
	kernelMed.setArg(kernelIndMed++, Nx);
	kernelMed.setArg(kernelIndMed++, Ny);
	kernelMed.setArg(kernelIndMed++, Nz);
	cl_int status = (*OpenCLStruct.af_queue).enqueueNDRangeKernel(kernelMed, cl::NullRange, global_size, cl::NullRange);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to launch the Median filter kernel\n");
		mexEvalString("pause(.0001);");
	}
	else if (DEBUG) {
		mexPrintf("Median kernel launched successfully\n");
		mexEvalString("pause(.0001);");
	}
	status = (*OpenCLStruct.af_queue).finish();
	grad.unlock();
	padd.unlock();
	af::sync();
#else
	CUresult status = CUDA_SUCCESS;
	const af::dim4 dimmi(padd.dims(0), padd.dims(1), padd.dims(2));
	af::array grad = af::constant(0.f, dimmi);
	padd = af::flat(padd);
	grad = af::flat(grad);
	CUdeviceptr* d_grad = grad.device<CUdeviceptr>();
	CUdeviceptr* d_padd = padd.device<CUdeviceptr>();
	af::sync();
	void* args[] = { reinterpret_cast<void*>(&d_padd), reinterpret_cast<void*>(&d_grad), (void*)&Nx, (void*)&Ny , (void*)&Nz};
	status = cuLaunchKernel(OpenCLStruct.kernelMed, dimmi[0], dimmi[1], dimmi[2], 1, 1, 1, 0, *OpenCLStruct.af_cuda_stream, &args[0], 0);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to launch the median filter kernel\n");
		mexEvalString("pause(.0001);");
	}
	status = cuCtxSynchronize();
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Queue finish failed after kernel\n");
		mexEvalString("pause(.0001);");
	}
	grad.unlock();
	padd.unlock();
	if (DEBUG) {
		mexPrintf("dims(0) = %d\n", dimmi[0]);
		mexPrintf("dims(1) = %d\n", dimmi[1]);
		mexPrintf("dims(2) = %d\n", dimmi[2]);
		mexPrintf("grad = %f\n", af::sum<float>(grad));
		mexPrintf("padd = %f\n", af::sum<float>(padd));
	}
	//padd = af::flat(padd);
	//af::array grad = af::median(af::moddims(padd(af::flat(offsets)), im_dim, offsets.dims(1)), 1);
#endif
	grad = af::moddims(grad, dimmi);
	grad = grad(af::seq(medx, Nx + medx - 1), af::seq(medy, Ny + medy - 1), af::seq(medz, Nz + medz - 1));
	grad = af::flat(grad) + epps;
	//if (DEBUG) {
	//	mexPrintf("grad = %f\n", af::sum<float>(grad));
	//	mexPrintf("im = %f\n", af::sum<float>(im));
	//	mexPrintf("erotus = %f\n", af::sum<float>(im - grad));
	//}
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - (grad)) / (grad);
	af::sync();
	if (DEBUG) {
		mexPrintf("grad2 = %f\n", af::sum<float>(grad));
	}
	return grad;
}

af::array Quadratic_prior(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t inffi,
	const af::array &offsets, const af::array &weights_quad, const uint32_t im_dim)
{
	const af::array apu_pad = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
	af::array grad;
	af::array weights = weights_quad;
	//weights(weights.elements() / 2) *= -1.f;
	//af::eval(weights);
	//weights = af::moddims(weights, weights_quad.dims(0), weights_quad.dims(1), weights_quad.dims(2));
	//if (DEBUG) {
	//	const char* R = af::toString("weights", weights);
	//	mexPrintf("weights = %s\n", R);
	//	mexPrintf("weights.elements() / 2 = %d\n", weights.elements() / 2);
	//}
	if (Ndz == 0 || Nz == 1) {
		grad = af::convolve2(apu_pad, weights);
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::span);
	}
	else {
		grad = af::convolve3(apu_pad, weights);
		grad = grad(af::seq(Ndx, Nx + Ndx - 1), af::seq(Ndy, Ny + Ndy - 1), af::seq(Ndz, Nz + Ndz - 1));
	}
	grad = af::flat(grad);
	return grad;
}

af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t inffi,
	const af::array& offsets, const af::array& weights_huber, const uint32_t im_dim, const float delta)
{
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
		//grad(af::span, ii) = af::sum(af::batchFunc(apu_pad, af::transpose(fmh_weights(af::span, ii)), batchMul), 1);
		grad(af::span, ii) = af::matmul(apu_pad, fmh_weights(af::span, ii));
	}
	indeksi1 = offsets.col(alku_fmh);
	grad(af::span, af::end) = padd(indeksi1 + 0U);
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
	grad = af::sum(af::batchFunc(apu_pad, af::transpose(a_L), batchMul), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + epps);
	return grad;
}

af::array Weighted_mean(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, 
	const af::array& weighted_weights, const bool med_no_norm, const uint32_t im_dim, const uint32_t mean_type, const float w_sum)
{
	af::array grad = af::constant(0.f, Nx, Ny, Nz);
	const float wsum = af::sum<float>(af::flat(weighted_weights));
	if (mean_type == 1U) {
		af::array padd = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || Nz == 1)
			grad = af::convolve2(padd, weighted_weights / wsum);
		else
			grad = af::convolve3(padd, weighted_weights / wsum);
	}
	else if (mean_type == 2U) {
		af::array padd = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || Nz == 1)
			grad = 1.f / af::convolve2(1.f / padd, weighted_weights / wsum);
		else
			grad = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
	}
	else if (mean_type == 3U) {
		af::array padd = padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || Nz == 1)
			grad = af::exp(af::convolve2(af::log(padd), weighted_weights / wsum));
		else
			grad = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
	}
	else if (mean_type == 4U) {
		grad = af::constant(0.f, im.dims(0));
		af::array padd = padding(im, Nx, Ny, Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::convolve3(padd, weighted_weights / wsum);
		//mexPrintf("Nx + Ndx * 4 - 1 = %u\n", Nx + Ndx * 4 - 1);
		if (Ndz == 0 || Nz == 1) {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::seq(Ndz, Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + Nx - 1), af::seq(kk, kk + Ny - 1), af::seq(ll, ll + Nz - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		grad = batchFunc(im, mm, batchDiv);
		grad(grad < 1e-3f) = 1e-3f;
		grad = batchFunc(af::log(grad), af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 5U) {
		af::array padd = padding(im, Nx, Ny, Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
		if (Ndz == 0 || Nz == 1) {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::seq(Ndz, Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + Nx - 1), af::seq(kk, kk + Ny - 1), af::seq(ll, ll + Nz - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-1f) = 1e-1f;
		grad = batchFunc(batchFunc(im_t, mm, batchMinus), im_t, batchDiv);
		grad = grad - af::pow(batchFunc((batchFunc(im_t, mm, batchMinus)), std::sqrt(2.f) * (im_t), batchDiv), 2.);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 6U) {
		af::array padd = padding(im, Nx, Ny, Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
		if (Ndz == 0 || Nz == 1) {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, Nx + Ndx * 3 - 1), af::seq(Ndy, Ny + Ndy * 3 - 1), af::seq(Ndz, Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + Nx - 1), af::seq(kk, kk + Ny - 1), af::seq(ll, ll + Nz - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-1f) = 1e-1f;
		grad = 1.f - batchFunc(mm, im_t, batchDiv);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else {
		mexWarnMsgTxt("Unsupported mean type");
	}
	if (mean_type <= 3U) {
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
	}
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


// Compute the TV prior
af::array TVprior(const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const TVdata& S, const af::array& ima, const float epps, const uint32_t TVtype,
	const Weighting& w_vec, const af::array& offsets) {
	af::array gradi;

	if (TVtype != 3U) {
		const af::array im = af::moddims(ima, Nx, Ny, Nz);
		// 1st order differentials
		af::array g = af::constant(0.f, Nx, Ny, Nz);
		af::array f = af::constant(0.f, Nx, Ny, Nz);
		af::array h = af::constant(0.f, Nx, Ny, Nz);
		f(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(im);
		//f(Nx - 1u, af::span, af::span) = im(im.dims(0) - 1ULL, af::span, af::span) - im(0, af::span, af::span);
		f(af::end, af::span, af::span) = f(Nx - 2u, af::span, af::span) * -1.f;
		g(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(im, 1);
		//g(af::span, Ny - 1u, af::span) = im(af::span, im.dims(1) - 1ULL, af::span) - im(af::span, 0, af::span);
		g(af::span, af::end, af::span) = g(af::span, Ny - 2u, af::span) * -1.f;
		h(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(im, 2);
		//h(af::span, af::span, Nz - 1u) = im(af::span, af::span, im.dims(2) - 1ULL) - im(af::span, af::span, 0);
		h(af::span, af::span, af::end) = h(af::span, af::span, Nz - 2u) * -1.f;

		af::array pval, apu1, apu2, apu3, apu4;
		g = af::flat(g);
		f = af::flat(f);
		h = af::flat(h);

			// If anatomical prior is used
			if (S.TV_use_anatomical || TVtype == 5U) {
				if (TVtype == 1U) {
					pval = af::sqrt(S.s1 * af::pow(f, 2.) + S.s5 * af::pow(g, 2.) + S.s9 * af::pow(h, 2.) + S.s4 * f * g + S.s7 * f * h + S.s2 * f * g + S.s8 * h * g + S.s3 * f * h + S.s6 * h * g + S.TVsmoothing);
					apu1 = 0.5f * (2.f * S.s1 * f + S.s4 * g + S.s7 * h + S.s2 * g + S.s3 * h) / pval;
					apu2 = 0.5f * (2.f * S.s5 * g + S.s4 * f + S.s2 * f + S.s8 * h + S.s6 * h) / pval;
					apu3 = 0.5f * (2.f * S.s9 * h + S.s8 * g + S.s6 * g + S.s7 * f + S.s3 * f) / pval;
					apu4 = 0.5f * (2.f * S.s1 * f + 2.f * S.s5 * g + 2.f * S.s9 * h + S.s4 * f + S.s2 * f + S.s8 * h + S.s6 * h + S.s4 * g + S.s7 * h + S.s2 * g + S.s3 * h
						+ S.s8 * g + S.s6 * g + S.s7 * f + S.s3 * f) / pval;
				}
				else if (TVtype == 2U) {
					const af::array reference_image = af::moddims(S.reference_image, Nx, Ny, Nz);
					af::array gp = af::constant(0.f, Nx, Ny, Nz);
					af::array fp = af::constant(0.f, Nx, Ny, Nz);
					af::array hp = af::constant(0.f, Nx, Ny, Nz);
					fp(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(reference_image);
					//fp(Nx - 1u, af::span, af::span) = S.reference_image(im.dims(0) - 1ULL, af::span, af::span) - S.reference_image(0, af::span, af::span);
					fp(af::end, af::span, af::span) = fp(Nx - 2u, af::span, af::span) * -1.f;
					gp(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(reference_image, 1);
					//gp(af::span, Ny - 1u, af::span) = S.reference_image(af::span, im.dims(1) - 1ULL, af::span) - S.reference_image(af::span, 0, af::span);
					gp(af::span, af::end, af::span) = gp(af::span, Ny - 2u, af::span) * -1.f;
					hp(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(reference_image, 2);
					//hp(af::span, af::span, Nz - 1u) = S.reference_image(af::span, af::span, im.dims(2) - 1ULL) - S.reference_image(af::span, af::span, 0);
					hp(af::span, af::span, af::end) = hp(af::span, af::span, Nz - 2u) * -1.f;

					gp = af::flat(gp);
					fp = af::flat(fp);
					hp = af::flat(hp);

					pval = af::sqrt(af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) + S.T * (af::pow(fp, 2.) + af::pow(gp, 2.) + af::pow(hp, 2.)) + S.TVsmoothing);
					apu1 = f / pval;
					apu2 = g / pval;
					apu3 = h / pval;
					apu4 = (f + g + h) / pval;
				}
				// For APLS
				else if (TVtype == 5U) {
					af::array gp = af::constant(0.f, Nx, Ny, Nz);
					af::array fp = af::constant(0.f, Nx, Ny, Nz);
					af::array hp = af::constant(0.f, Nx, Ny, Nz);
					fp(af::seq(0, Nx - 2u), af::span, af::span) = -af::diff1(S.APLSReference);
					//fp(Nx - 1u, af::span, af::span) = S.APLSReference(im.dims(0) - 1ULL, af::span, af::span) - S.APLSReference(0, af::span, af::span);
					fp(af::end, af::span, af::span) = fp(Nx - 2u, af::span, af::span) * -1.f;
					gp(af::span, af::seq(0, Ny - 2u), af::span) = -af::diff1(S.APLSReference, 1);
					//gp(af::span, Ny - 1u, af::span) = S.APLSReference(af::span, im.dims(1) - 1ULL, af::span) - S.APLSReference(af::span, 0, af::span);
					gp(af::span, af::end, af::span) = gp(af::span, Ny - 2u, af::span) * -1.f;
					hp(af::span, af::span, af::seq(0, Nz - 2u)) = -af::diff1(S.APLSReference, 2);
					//hp(af::span, af::span, Nz - 1u) = S.APLSReference(af::span, af::span, im.dims(2) - 1ULL) - S.APLSReference(af::span, af::span, 0);
					hp(af::span, af::span, af::end) = hp(af::span, af::span, Nz - 2u) * -1.f;

					fp = af::flat(fp) + epps;
					gp = af::flat(gp) + epps;
					hp = af::flat(hp) + epps;

					//af::array epsilon = af::join(1, fp, gp, hp) / af::join(1, fp / af::sqrt(af::pow(fp, 2.) + S.eta * S.eta) + epps,
					//	gp / af::sqrt(af::pow(gp, 2.) + S.eta * S.eta) + epps, hp / af::sqrt(af::pow(hp, 2.) + S.eta * S.eta) + epps);
					const af::array epsilon = af::batchFunc(af::join(1, fp, gp, hp), af::sqrt(fp * fp + gp * gp + hp * hp + S.eta * S.eta), batchDiv);
					const af::array apu = af::sum(af::join(1, f, g, h) * epsilon, 1);

					pval = (f * f + g * g + h * h - apu * apu + S.APLSsmoothing);
					pval(pval <= 0.f) = S.APLSsmoothing;
					pval = af::sqrt(pval);
					apu1 = (f - (apu * epsilon(af::span, 0))) / pval;
					apu2 = (g - (apu * epsilon(af::span, 1))) / pval;
					apu3 = (h - (apu * epsilon(af::span, 2))) / pval;
					apu4 = (f - (apu * epsilon(af::span, 0)) + g - (apu * epsilon(af::span, 1)) + h - (apu * epsilon(af::span, 2))) / pval;
				}
			}
			// If anatomical prior is not used
			else {
				if (TVtype == 4U) {
					//af::array ff = af::constant(0.f, f.dims(0), f32);
					//ff(f > 0) = 1.f;
					//ff(f < 0) = -1.f;
					//af::array gg = af::constant(0.f, g.dims(0), f32);
					//gg(g > 0) = 1.f;
					//gg(g < 0) = -1.f;
					//af::array hh = af::constant(0.f, h.dims(0), f32);
					//hh(h > 0) = 1.f;
					//hh(h < 0) = -1.f;
					const af::array ff = f / af::abs(f);
					const af::array gg = g / af::abs(g);
					const af::array hh = h / af::abs(h);
					if (S.SATVPhi == 0.f) {
						apu1 = ff;
						apu2 = gg;
						apu3 = hh;
					}
					else {
						apu1 = ff - ff / (af::abs(f) / S.SATVPhi + 1.f);
						apu2 = gg - gg / (af::abs(g) / S.SATVPhi + 1.f);
						apu3 = hh - hh / (af::abs(h) / S.SATVPhi + 1.f);
					}
				}
				else {
					pval = af::sqrt(af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) + S.TVsmoothing);
					apu1 = f / pval;
					apu2 = g / pval;
					apu3 = h / pval;
				}
				apu4 = apu1 + apu2 + apu3;
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
			gradi = gradi + 2.f * S.tau * af::min<float>(af::flat(ima));
	}
	else {
		if (S.TV_use_anatomical) {
			af::array padd = af::flat(padding(ima, Nx, Ny, Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd = padd(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
			padd = af::moddims(padd, ima.dims(0), padd.dims(0) / ima.dims(0));
			af::array padd2 = af::flat(padding(S.reference_image, Nx, Ny, Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd2 = padd2(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
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
		af::array apu = (af::max)(1.f, af::sqrt(eta1 * eta1 + eta2 * eta2 + eta3 * eta3) / beta);
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
		apu = (af::max)(1.f, af::sqrt(eta1 * eta1 + eta2 * eta2 + eta3 * eta3 + eta4 * eta4) / alpha);
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

af::array RDP(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const af::array& weights_RDP, 
	const uint32_t im_dim, const float gamma, const af::array& offsets, const uint32_t inffi)
{
	const af::array im_apu = af::flat(padding(im, Nx, Ny, Nz, Ndx, Ndy, Ndz));

	const af::array indeksi2 = af::join(1, offsets.cols(0, inffi - 1), offsets.cols(inffi + 1, af::end));
	af::array grad = af::batchFunc(im, af::moddims(im_apu(indeksi2), im_dim, weights_RDP.dims(0)), batchDiv);
	grad = af::matmul((((grad - 1.f) * (gamma * af::abs(grad - 1.f) + grad + 3.f))  / af::pow2(grad + 1.f + gamma * af::abs(grad - 1.f))), weights_RDP);
	grad = af::flat(grad);
	return grad;
}


af::array computeConvolution(const af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec, 
	const uint32_t n_rekos) {
	af::array apu = af::moddims(vec, Nx, Ny, Nz, n_rekos);
	for (uint32_t ff = 0U; ff < n_rekos; ff++) {
		af::array apu2 = padding(apu(af::span, af::span, af::span, ff), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
		apu2 = af::convolve3(apu2, g);
		apu2 = apu2(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		apu(af::span, af::span, af::span, ff) = apu2;
	}
	return af::flat(apu);
}

void deblur(af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec, const uint32_t iter, 
	const uint32_t subsets, const float epps, const bool saveIter) {
	uint32_t it = 0U;
	if (saveIter)
		it = iter + 1U;
	af::array jelppi = padding(vec(af::span, it), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	af::array apu = padding(vec(af::span, it), Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	for (uint32_t kk = 0U; kk < subsets; kk++) {
		af::array apu2 = convolve3(jelppi, g) + epps;
		apu2 = apu2(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		apu2 = padding(apu2, Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
		jelppi *= af::convolve3(apu / apu2, g);
		jelppi = jelppi(af::seq(w_vec.g_dim_x + 1, Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, Ny + w_vec.g_dim_y), af::seq(w_vec.g_dim_z + 1, Nz + w_vec.g_dim_z));
		if (kk < subsets - 1)
			jelppi = padding(jelppi, Nx, Ny, Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	}
	jelppi = af::flat(jelppi);
	vec(af::span, it) = jelppi;
}

void computeDeblur(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps, const bool saveIter) {

	for (uint32_t kk = 0U; kk < w_vec.nTot; kk++)
		deblur(vec.imEstimates[kk], g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
}


//void computeDeblurMLEM(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
//	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps, const bool saveIter) {
//	if (MethodList.MLEM) {
//		deblur(vec.MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//	}
//	if (MethodList.OSLMLEM) {
//		if (MethodList.MRP) {
//			deblur(vec.MRP_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.Quad) {
//			deblur(vec.Quad_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.Huber) {
//			deblur(vec.Huber_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.L) {
//			deblur(vec.L_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.FMH) {
//			deblur(vec.FMH_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.WeightedMean) {
//			deblur(vec.Weighted_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.TV) {
//			deblur(vec.TV_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.AD) {
//			deblur(vec.AD_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.APLS) {
//			deblur(vec.APLS_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.TGV) {
//			deblur(vec.TGV_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//		if (MethodList.NLM) {
//			deblur(vec.NLM_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//	}
//	if (MethodList.CUSTOM) {
//		if (MethodList.OSLMLEM) {
//			deblur(vec.custom_MLEM, g, Nx, Ny, Nz, w_vec, iter, subsets, epps, saveIter);
//		}
//	}
//}