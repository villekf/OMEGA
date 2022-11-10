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
void form_data_variables(AF_im_vectors & vec, std::vector<float> & beta, Weighting & w_vec, const mxArray *options, scalarStruct& inputScalars, 
	const RecMethods &MethodList, TVdata &data, const uint32_t Nt, const uint32_t iter0, std::vector<std::vector<float*>>& imEstimates)
{
	// Load the number of priors, all MAP-algorithms, non-OS MAP algorithms, non-OS non-MAP algorithms, non-MAP algorithms and total number of algorithms
	w_vec.nPriors = getScalarUInt32(getField(options, 0, "nPriors"), -2);
	if (MethodList.CUSTOM)
		w_vec.nPriorsTot = w_vec.nPriors + 1U;
	else
		w_vec.nPriorsTot = w_vec.nPriors;
	w_vec.nMAP = getScalarUInt32(getField(options, 0, "nMAP"), -2);
	w_vec.nMAPOS = w_vec.nMAP - w_vec.nMAPML;
	w_vec.nOS = getScalarUInt32(getField(options, 0, "nOS"), -5);
	w_vec.nTot = getScalarUInt32(getField(options, 0, "nTot"), -6);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	w_vec.rekot = (uint32_t*)mxGetUint32s(getField(options, 0, "rekoList"));
#else
	w_vec.rekot = (uint32_t*)mxGetData(getField(options, 0, "rekoList"));
#endif
	w_vec.mIt.push_back(getScalarInt32(mxGetCell(getField(options, 0, "mIt"), 0), -7));
	w_vec.mIt.push_back(getScalarInt32(mxGetCell(getField(options, 0, "mIt"), 1), -8));
	if (DEBUG) {
		mexPrintf("nPriors = %u\n", w_vec.nPriors);
		mexPrintf("nMAP = %u\n", w_vec.nMAP);
		mexPrintf("nMAPML = %u\n", w_vec.nMAPML);
		mexPrintf("nTot = %u\n", w_vec.nTot);
		mexEvalString("pause(.0001);");
	}
	uint32_t Ni = 1U;
	if (inputScalars.saveIter)
		Ni = inputScalars.Niter + 1U;
	// Load the necessary variables if the corresponding reconstruction method is used and set the initial value
	int yy = 0;
	// First non-MAP/prior-based algorithms (e.g. MLEM)
	//if (MethodList.OSEM) {
	//	imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
	//	if (saveIter)
	//		imEstimates[yy](af::span, 0) = x0;
	//	yy++;
	//}
	//
	//if (MethodList.MRAMLA) {
	//	imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
	//	imEstimates[yy](af::span, 0) = x0;
	//	yy++;
	//}
	//
	//if (MethodList.RAMLA) {
	//	imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
	//	imEstimates[yy](af::span, 0) = x0;
	//	yy++;
	//}
	//
	//if (MethodList.ROSEM) {
	//	imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
	//	imEstimates[yy](af::span, 0) = x0;
	//	yy++;
	//}
	//
	//if (MethodList.RBI) {
	//	imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
	//	imEstimates[yy](af::span, 0) = x0;
	//	yy++;
	//}
	
	if (MethodList.DRAMA) {
		//imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
		//imEstimates[yy](af::span, 0) = x0;
		//yy++;
		
		// Relaxation parameter
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.lambda_DRAMA = (float*)mxGetSingles(getField(options, 0, "lam_drama"));
#else
		w_vec.lambda_DRAMA = (float*)mxGetData(getField(options, 0, "lam_drama"));
#endif
	}
	
	if (MethodList.COSEM) {
		//imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
		//imEstimates[yy](af::span, 0) = x0;
		//yy++;
		
		// Complete data
		vec.C_co = af::constant(0.f, inputScalars.im_dim, inputScalars.subsets);
	}
	if (MethodList.ECOSEM) {
		//imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
		//imEstimates[yy](af::span, 0) = x0;
		//yy++;
		
		if (!MethodList.COSEM) {
			
			// Complete data
			vec.C_co = af::constant(0.f, inputScalars.im_dim, inputScalars.subsets);
		}
	}
	if (MethodList.ACOSEM) {
		//imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
		//imEstimates[yy](af::span, 0) = x0;
		//yy++;
		
		// Complete data
		vec.C_aco = af::constant(0.f, inputScalars.im_dim, inputScalars.subsets);
	}

	if (MethodList.OSLCOSEM > 0)
		vec.C_osl = af::constant(0.f, inputScalars.im_dim, inputScalars.subsets);

	// Load the MAP/prior-based algorithms
	int tt = 0;
	for (uint32_t kk = 0; kk < w_vec.nPriors; kk++) {
		for (uint32_t ll = 0; ll < w_vec.nMAP; ll++) {
			const char* varChar = mxArrayToString(mxGetCell(getField(options, 0, "varList"), tt));
			//imEstimates.push_back(af::constant(0.f, inputScalars.im_dim, Ni));
			//imEstimates[yy](af::span, 0) = x0;
			if (DEBUG) {
				mexPrintf("%s\n", varChar);
				mexEvalString("pause(.0001);");
			}
			// Load the regularization parameter as well if the prior is used
			beta.push_back(getScalarFloat(getField(options, 0, varChar), -9));
			yy++;
			tt++;
		}
	}

	if (inputScalars.maskFP)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.maskFP = (uint8_t*)mxGetUint8s(getField(options, 0, "maskFP"));
#else
		w_vec.maskFP = (uint8_t*)mxGetData(getField(options, 0, "maskFP"));
#endif
	if (inputScalars.maskBP && (inputScalars.projector_type == 5 || inputScalars.projector_type == 4 || inputScalars.projector_type == 14)) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.maskBP = (uint8_t*)mxGetUint8s(getField(options, 0, "maskBP"));
#else
		w_vec.maskBP = (uint8_t*)mxGetData(getField(options, 0, "maskBP"));
#endif
		const size_t nMask = mxGetNumberOfElements(getField(options, 0, "maskBP"));
		if (DEBUG) {
			mexPrintf("nMask = %d\n", nMask);
			mexEvalString("pause(.0001);");
		}
	}
	// CT-related variables such as number of projection images
	if (inputScalars.CT) {
		w_vec.size_y = getScalarUInt32(getField(options, 0, "xSize"), -10);
		w_vec.size_x = getScalarUInt32(getField(options, 0, "ySize"), -10);
		w_vec.nProjections = getScalarInt64(getField(options, 0, "nProjections"), -11);
		//w_vec.dPitch = getScalarFloat(getField(options, 0, "dPitch"), -12);
		w_vec.dPitchX = (float)mxGetScalar(getField(options, 0, "dPitchX"));
		w_vec.dPitchY = (float)mxGetScalar(getField(options, 0, "dPitchY"));
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		w_vec.uv = (float*)mxGetSingles(getField(options, 0, "uV"));
//#else
//		w_vec.uv = (float*)mxGetData(getField(options, 0, "uV"));
//#endif
	}
	//else if (inputScalars.PET) {
	//	w_vec.nProjections = (int64_t)mxGetScalar(getField(options, 0, "nProjections"));
	//	w_vec.size_y = (uint32_t)mxGetScalar(getField(options, 0, "Nang"));
	//	w_vec.size_x = (uint32_t)mxGetScalar(getField(options, 0, "Ndist"));
	//	w_vec.dPitchX = (float)mxGetScalar(getField(options, 0, "cr_p"));
	//	w_vec.dPitchY = (float)mxGetScalar(getField(options, 0, "cr_pz"));
	//}
	else {
		w_vec.size_y = (uint32_t)mxGetScalar(getField(options, 0, "Nang"));
		w_vec.size_x = (uint32_t)mxGetScalar(getField(options, 0, "Ndist"));
		w_vec.nProjections = (int64_t)mxGetScalar(getField(options, 0, "nProjections"));
		// Detector pitch
		w_vec.dPitchX = (float)mxGetScalar(getField(options, 0, "cr_p"));
		w_vec.dPitchY = (float)mxGetScalar(getField(options, 0, "cr_pz"));
	}
	if (inputScalars.projector_type == 4U || inputScalars.projector_type == 14 || inputScalars.projector_type == 41)
		w_vec.kerroin4 = (float)mxGetScalar(getField(options, 0, "kerroin"));

	if (inputScalars.projector_type == 6U) {
		const mwSize* ind = mxGetDimensions(getField(options, 0, "gFilter"));
		if (DEBUG) {
			mexPrintf("indX = %d\n", ind[0]);
			mexPrintf("indY = %d\n", ind[1]);
			mexEvalString("pause(.0001);");
		}
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.angles = (float*)mxGetSingles(getField(options, 0, "angles"));
		w_vec.gFilter = af::array(ind[0], ind[1], ind[2], ind[3], (float*)mxGetSingles(getField(options, 0, "gFilter")));
		w_vec.distInt = (uint32_t*)mxGetUint32s(getField(options, 0, "blurPlanes"));
#else
		w_vec.angles = (float*)mxGetData(getField(options, 0, "angles"));
		w_vec.gFilter = af::array(ind[0], ind[1], ind[2], ind[3], (float*)mxGetData(getField(options, 0, "angles")));
		w_vec.distInt = (uint32_t*)mxGetData(getField(options, 0, "blurPlanes"));
#endif
		if (DEBUG) {
			mexPrintf("w_vec.gFilter.dims(0) = %d\n", w_vec.gFilter.dims(0));
			mexPrintf("w_vec.gFilter.dims(1) = %d\n", w_vec.gFilter.dims(1));
			mexPrintf("w_vec.gFilter.dims(2) = %d\n", w_vec.gFilter.dims(2));
			mexPrintf("w_vec.gFilter.dims(3) = %d\n", w_vec.gFilter.dims(3));
			mexPrintf("w_vec.distInt[0] = %d\n", w_vec.distInt[0]);
			mexEvalString("pause(.0001);");
		}
	}

	// Load TV related input data
	if (MethodList.TV && MethodList.MAP) {
		// Is anatomical reference image used
		data.TV_use_anatomical = getScalarBool(getField(options, 0, "TV_use_anatomical"), -13);
		// Tau-value
		data.tau = getScalarFloat(getField(options, 0, "tau"), -14);
		// "Smoothing" parameter, prevents zero values in the square root
		data.TVsmoothing = getScalarFloat(getField(options, 0, "TVsmoothing"), -15);
		// The type of TV prior used
		data.TVtype = getScalarUInt32(getField(options, 0, "TVtype"), -16);
		// If anatomical prior is used, load the necessary coefficients
		if (data.TV_use_anatomical) {
			mxArray* TVdata_init = getField(options, 0, "TVdata");
			if (data.TVtype == 1) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
				data.s1 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s1")), afHost);
				data.s2 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s2")), afHost);
				data.s3 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s3")), afHost);
				data.s4 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s4")), afHost);
				data.s5 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s5")), afHost);
				data.s6 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s6")), afHost);
				data.s7 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s7")), afHost);
				data.s8 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s8")), afHost);
				data.s9 = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "s9")), afHost);
#else
				data.s1 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s1")), afHost);
				data.s2 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s2")), afHost);
				data.s3 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s3")), afHost);
				data.s4 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s4")), afHost);
				data.s5 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s5")), afHost);
				data.s6 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s6")), afHost);
				data.s7 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s7")), afHost);
				data.s8 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s8")), afHost);
				data.s9 = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "s9")), afHost);
#endif
			}
			else {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
				data.reference_image = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(TVdata_init, 0, "reference_image")), afHost);
#else
				data.reference_image = af::array(inputScalars.im_dim, (float*)mxGetData(getField(TVdata_init, 0, "reference_image")), afHost);
#endif
			}
			data.T = getScalarFloat(getField(options, 0, "T"), -17);
			data.C = getScalarFloat(getField(options, 0, "C"), -18);
		}
		// Additional weights for the TV type 3
		if (data.TVtype == 3 && !MethodList.Quad) {
			w_vec.Ndx = getScalarUInt32(getField(options, 0, "Ndx"), -19);
			w_vec.Ndy = getScalarUInt32(getField(options, 0, "Ndy"), -20);
			w_vec.Ndz = getScalarUInt32(getField(options, 0, "Ndz"), -21);
			w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			w_vec.weights_TV = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(getField(options, 0, "weights_quad")), afHost);
#else
			w_vec.weights_TV = af::array(w_vec.dimmu - 1, (float*)mxGetData(getField(options, 0, "weights_quad")), afHost);
#endif
		}
		if (data.TVtype == 3)
			data.C = getScalarFloat(getField(options, 0, "C"), -22);
		if (data.TVtype == 4)
			data.SATVPhi = getScalarFloat(getField(options, 0, "SATVPhi"), -23);
	}
	// General variables for neighborhood-based methods
	if ((MethodList.L || MethodList.FMH || MethodList.WeightedMean || MethodList.Quad || MethodList.Huber || MethodList.MRP || MethodList.NLM || (data.TVtype == 3 && MethodList.TV) || MethodList.RDP) && MethodList.MAP) {
		// Neighborhood size
		w_vec.Ndx = getScalarUInt32(getField(options, 0, "Ndx"), -24);
		w_vec.Ndy = getScalarUInt32(getField(options, 0, "Ndy"), -25);
		w_vec.Ndz = getScalarUInt32(getField(options, 0, "Ndz"), -26);
		// Is normalization used in MRP, FMH, L, weighted mean or AD
		w_vec.med_no_norm = getScalarBool(getField(options, 0, "med_no_norm"), -27);
		w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
	}
	if ((MethodList.L || MethodList.FMH || (data.TVtype == 3 && MethodList.TV)) && MethodList.MAP) {
		// Index values for the neighborhood
//#ifdef OPENCL
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.tr_offsets = af::array(inputScalars.im_dim, w_vec.dimmu, (uint32_t*)mxGetUint32s(getField(options, 0, "tr_offsets")), afHost);
#else
		w_vec.tr_offsets = af::array(inputScalars.im_dim, w_vec.dimmu, (uint32_t*)mxGetData(getField(options, 0, "tr_offsets")), afHost);
#endif
//#else
//		w_vec.tr_offsets = af::array(inputScalars.im_dim, w_vec.dimmu, (uint32_t*)mxGetData(getField(options, 0, "tr_offsets")), afHost);
//#endif
	}
	if (MethodList.FMH || MethodList.Quad || MethodList.Huber)
		w_vec.inffi = getScalarUInt32(getField(options, 0, "inffi"), -28);
	// Weights for the various priors
	if (MethodList.Quad && MethodList.MAP) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(getField(options, 0, "weights_quad")), afHost);
#else
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, (float*)mxGetData(getField(options, 0, "weights_quad")), afHost);
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
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		w_vec.weights_RDP = af::array((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, (float*)mxGetSingles(getField(options, 0, "weights_RDP")), afHost);
//#else
//		w_vec.weights_RDP = af::array((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, (float*)mxGetData(getField(options, 0, "weights_RDP")), afHost);
//#endif
//		w_vec.weights_RDP = af::flat(w_vec.weights_RDP);
		w_vec.RDP_gamma = getScalarFloat(getField(options, 0, "RDP_gamma"), -29);
	}
	if (MethodList.Huber && MethodList.MAP) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, (float*)mxGetSingles(getField(options, 0, "weights_huber")), afHost);
#else
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, (float*)mxGetData(getField(options, 0, "weights_huber")), afHost);
#endif
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_huber * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_huber = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
		w_vec.huber_delta = getScalarFloat(getField(options, 0, "huber_delta"), -29);
	}
	if (MethodList.L && MethodList.MAP)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.a_L = af::array(w_vec.dimmu, (float*)mxGetSingles(getField(options, 0, "a_L")), afHost);
#else
		w_vec.a_L = af::array(w_vec.dimmu, (float*)mxGetData(getField(options, 0, "a_L")), afHost);
#endif
	if (MethodList.FMH && MethodList.MAP) {
		if (inputScalars.Nz == 1 || w_vec.Ndz == 0)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, (float*)mxGetSingles(getField(options, 0, "fmh_weights")), afHost);
#else
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, (float*)mxGetData(getField(options, 0, "fmh_weights")), afHost);
#endif
		else
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, (float*)mxGetSingles(getField(options, 0, "fmh_weights")), afHost);
#else
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, (float*)mxGetData(getField(options, 0, "fmh_weights")), afHost);
#endif
		w_vec.alku_fmh = getScalarUInt32(getField(options, 0, "inffi"), -30);
	}
	if (MethodList.WeightedMean && MethodList.MAP) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, (float*)mxGetSingles(getField(options, 0, "weighted_weights")), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
#else
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, (float*)mxGetData(getField(options, 0, "weighted_weights")), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
#endif
		// Type of mean used (arithmetic, harmonic or geometric)
		w_vec.mean_type = getScalarInt32(getField(options, 0, "mean_type"), -31);
		// Sum of the weights
		w_vec.w_sum = getScalarFloat(getField(options, 0, "w_sum"), -31);
	}
	if (MethodList.AD && MethodList.MAP) {
		// Time-step value
		w_vec.TimeStepAD = getScalarFloat(getField(options, 0, "TimeStepAD"), -32);
		// Conductance (edge value)
		w_vec.KAD = getScalarFloat(getField(options, 0, "KAD"), -33);
		// Number of AD iterations
		w_vec.NiterAD = getScalarUInt32(getField(options, 0, "NiterAD"), -34);
		// Flux type
		uint32_t Flux = getScalarUInt32(getField(options, 0, "FluxType"), -35);
		// Diffusion type
		uint32_t Diffusion = getScalarUInt32(getField(options, 0, "DiffusionType"), -36);
		if (Flux == 2U)
			w_vec.FluxType = AF_FLUX_QUADRATIC;
		else
			w_vec.FluxType = AF_FLUX_EXPONENTIAL;
		if (Diffusion == 2U)
			w_vec.DiffusionType = AF_DIFFUSION_MCDE;
		else
			w_vec.DiffusionType = AF_DIFFUSION_GRAD;
		if (!MethodList.L && !MethodList.FMH && !MethodList.WeightedMean && !MethodList.Quad && !MethodList.MRP)
			w_vec.med_no_norm = getScalarBool(getField(options, 0, "med_no_norm"), -37);
	}
	if (MethodList.APLS && MethodList.MAP) {
		// Eta value
		data.eta = getScalarFloat(getField(options, 0, "eta"), -38);
		// Tau-value
		if (!MethodList.TV)
			data.tau = getScalarFloat(getField(options, 0, "tau"), -39);
		// Smoothing value
		data.APLSsmoothing = getScalarFloat(getField(options, 0, "APLSsmoothing"), -40);
		// Anatomical reference image
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		data.APLSReference = af::array(inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, (float*)mxGetSingles(getField(options, 0, "APLS_ref_image")), afHost);
#else
		data.APLSReference = af::array(inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, (float*)mxGetData(getField(options, 0, "APLS_ref_image")), afHost);
#endif
	}
	if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.COSEM
		|| MethodList.ECOSEM || MethodList.ACOSEM || MethodList.PKMA || MethodList.OSLCOSEM > 0) {
		// Sum of the rows (measurements) of the system matrix
		//w_vec.D = af::array(inputScalars.im_dim, (float*)mxGetData(getField(options, 0, "D")), afHost);
		w_vec.D = af::constant(0.f, inputScalars.im_dim, 1);
		// For manual determination of the upper bound
		if ((MethodList.MRAMLA || MethodList.MBSREM) && Nt > 1U)
			w_vec.Amin = af::constant(0.f, inputScalars.koko, 1);
		else
			w_vec.Amin = af::constant(0.f, 1, 1);
		//w_vec.Amin = af::array(koko_l, (float*)mxGetData(getField(options, 0, "Amin")), afHost);
		w_vec.MBSREM_prepass = getScalarBool(getField(options, 0, "MBSREM_prepass"), -41);
	}
	if (MethodList.MRAMLA || MethodList.MBSREM) {
		// Relaxation parameter
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.lambda_MBSREM = (float*)mxGetSingles(getField(options, 0, "lam_MBSREM"));
#else
		w_vec.lambda_MBSREM = (float*)mxGetData(getField(options, 0, "lam_MBSREM"));
#endif
		// Upper bound
		w_vec.U = getScalarFloat(getField(options, 0, "U"), -42);
	}
	if (DEBUG) {
		mexPrintf("w_vec.lambda_MBSREM = %f\n", w_vec.lambda_MBSREM);
		mexEvalString("pause(.0001);");
	}
	// Relaxation parameters
	if (MethodList.RAMLA || MethodList.BSREM)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.lambda_BSREM = (float*)mxGetSingles(getField(options, 0, "lam"));
#else
		w_vec.lambda_BSREM = (float*)mxGetData(getField(options, 0, "lam"));
#endif
	if (MethodList.ROSEM || MethodList.ROSEMMAP)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.lambda_ROSEM = (float*)mxGetSingles(getField(options, 0, "lam_ROSEM"));
#else
		w_vec.lambda_ROSEM = (float*)mxGetData(getField(options, 0, "lam_ROSEM"));
#endif
	if (MethodList.PKMA) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.lambda_PKMA = (float*)mxGetSingles(getField(options, 0, "lam_PKMA"));
		w_vec.alpha_PKMA = (float*)mxGetSingles(getField(options, 0, "alpha_PKMA"));
		w_vec.sigma_PKMA = (float*)mxGetSingles(getField(options, 0, "sigma_PKMA"));
#else
		w_vec.lambda_PKMA = (float*)mxGetData(getField(options, 0, "lam_PKMA"));
		w_vec.alpha_PKMA = (float*)mxGetData(getField(options, 0, "alpha_PKMA"));
		w_vec.sigma_PKMA = (float*)mxGetData(getField(options, 0, "sigma_PKMA"));
#endif
	}
	// Power factor for ACOSEM
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1)
		w_vec.h_ACOSEM = getScalarFloat(getField(options, 0, "h"), -43);
	w_vec.h_ACOSEM_2 = 1.f / w_vec.h_ACOSEM;
	if (MethodList.TGV && MethodList.MAP) {
		data.TGVAlpha = getScalarFloat(getField(options, 0, "alphaTGV"), -44);
		data.TGVBeta = getScalarFloat(getField(options, 0, "betaTGV"), -45);
		data.NiterTGV = getScalarUInt32(getField(options, 0, "NiterTGV"), -46);
	}
	if (MethodList.NLM && MethodList.MAP) {
		w_vec.NLM_anatomical = getScalarBool(getField(options, 0, "NLM_use_anatomical"), -47);
		w_vec.NLTV = getScalarBool(getField(options, 0, "NLTV"), -48);
		w_vec.NLM_MRP = getScalarBool(getField(options, 0, "NLM_MRP"), -49);
		if (w_vec.NLM_anatomical)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			w_vec.NLM_ref = (float*)mxGetSingles(getField(options, 0, "NLM_ref"));
//			w_vec.NLM_ref = af::array(inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, (float*)mxGetSingles(getField(options, 0, "NLM_ref")), afHost);
#else
			w_vec.NLM_ref = (float*)mxGetData(getField(options, 0, "NLM_ref"));
//			w_vec.NLM_ref = af::array(inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, (float*)mxGetData(getField(options, 0, "NLM_ref")), afHost);
#endif
		w_vec.h2 = getScalarFloat(getField(options, 0, "sigma"), -50);
		w_vec.h2 = w_vec.h2 * w_vec.h2;
		w_vec.Nlx = getScalarUInt32(getField(options, 0, "Nlx"), -51);
		w_vec.Nly = getScalarUInt32(getField(options, 0, "Nly"), -52);
		w_vec.Nlz = getScalarUInt32(getField(options, 0, "Nlz"), -53);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), (float*)mxGetSingles(getField(options, 0, "gaussianNLM")), afHost);
#else
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), (float*)mxGetData(getField(options, 0, "gaussianNLM")), afHost);
#endif
	}
	if (MethodList.CUSTOM) {
		for (uint32_t kk = 0; kk < imEstimates.size(); kk++) {
			const char* varTot = mxArrayToString(mxGetCell(getField(options, 0, "varTot"), kk));
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			imEstimates[kk][0] = (float*)mxGetSingles(getField(getField(options, 0, "im_vectors"), 0, varTot));
#else
			imEstimates[kk][0] = (float*)mxGetData(getField(getField(options, 0, "im_vectors"), 0, varTot));
#endif
			if (DEBUG) {
				mexPrintf("%s\n", varTot);
				//mexPrintf("imEstimates[kk](af::span, 0) = %f\n", af::sum<float>(imEstimates[kk](af::span, 0)));
				mexEvalString("pause(.0001);");
			}
		}
		tt = 0;
		for (uint32_t ll = 0; ll < w_vec.nMAP; ll++) {
			const char* varApu = mxArrayToString(mxGetCell(getField(options, 0, "varApu"), tt));
			const char* varBeta = mxArrayToString(mxGetCell(getField(options, 0, "varBeta"), tt));
			const char* varGrad = mxArrayToString(mxGetCell(getField(options, 0, "varGrad"), tt));
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			//imEstimates.push_back(af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(options, 0, varApu)), afHost));
			imEstimates.push_back(std::vector<float*>(imEstimates[0].size(), (float*)mxGetSingles(getField(options, 0, varApu))));
			w_vec.dU.push_back(af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(options, 0, varGrad)), afHost));
#else
			imEstimates.push_back(std::vector<float*>(imEstimates[0].size(), (float*)mxGetData(getField(options, 0, varApu))));
			//imEstimates.push_back(af::array(inputScalars.im_dim, (float*)mxGetData(getField(options, 0, varApu)), afHost));
			w_vec.dU.push_back(af::array(inputScalars.im_dim, (float*)mxGetData(getField(options, 0, varGrad)), afHost));
#endif
			beta.push_back(getScalarFloat(getField(options, 0, varBeta), -54));
			tt++;
		}
		if (MethodList.MBSREM) {
			if (iter0 > 0 || inputScalars.osa_iter0 > 0) {
				w_vec.U = getScalarFloat(getField(options, 0, "U"), -55);
				w_vec.epsilon_mramla = getScalarFloat(getField(options, 0, "epsilon_mramla"), -56);
			}
		}
		if (MethodList.OSLCOSEM > 0) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			vec.C_osl = af::array(inputScalars.im_dim, inputScalars.subsets, (float*)mxGetSingles(getField(options, 0, "C_osl")), afHost);
#else
			vec.C_osl = af::array(inputScalars.im_dim, inputScalars.subsets, (float*)mxGetData(getField(options, 0, "C_osl")), afHost);
#endif
		}
		if ((MethodList.OSLCOSEM > 0 || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA))
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			w_vec.D = af::array(inputScalars.im_dim, (float*)mxGetSingles(getField(options, 0, "D")), afHost);
#else
			w_vec.D = af::array(inputScalars.im_dim, (float*)mxGetData(getField(options, 0, "D")), afHost);
#endif

		//uint32_t dd = n_rekos_mlem + n_rekos - 1U - nMAPOS;
		//vec.im_mlem(af::seq(n_rekos_mlem * inputScalars.im_dim - inputScalars.im_dim, n_rekos_mlem * inputScalars.im_dim - 1u)) = imEstimates[dd];
		uint32_t dd = 0U;
		uint32_t yy = 0U;
		//uint32_t dd = n_rekos_mlem + n_rekos - 1U;
		//uint32_t yy = n_rekos * inputScalars.im_dim;
		if (DEBUG) {
			//mexPrintf("imEstimates.size() = %d\n", imEstimates.size());
			mexPrintf("[dd] = %u\n", dd);
			mexPrintf("[yy] = %u\n", yy);
			//mexPrintf("imEstimates[dd].dims(0) = %d\n", imEstimates[dd].dims(0));
			mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
			mexEvalString("pause(.0001);");
		}
		for (uint32_t kk = 0; kk < (w_vec.nPriorsTot * w_vec.nMAPOS + w_vec.nOS); kk++) {
			vec.im_os(af::seq(yy, yy + inputScalars.im_dim - 1u)) = af::array(inputScalars.im_dim, imEstimates[dd][0], afHost);
			yy += inputScalars.im_dim;
			dd++;
		}
	}
	if (inputScalars.use_psf) {
		w_vec.g_dim_x = getScalarUInt32(getField(options, 0, "g_dim_x"), -57);
		w_vec.g_dim_y = getScalarUInt32(getField(options, 0, "g_dim_y"), -58);
		w_vec.g_dim_z = getScalarUInt32(getField(options, 0, "g_dim_z"), -59);
		w_vec.deconvolution = getScalarBool(getField(options, 0, "deblurring"), -60);
	}
	if (inputScalars.use_psf && w_vec.deconvolution) {
		w_vec.deblur_iterations = (uint32_t)mxGetScalar(getField(options, 0, "deblur_iterations"));
	}
	if (MethodList.CPLS || MethodList.CPTV || MethodList.CPTVKL || MethodList.CPLSKL) {
		w_vec.tauCP = getScalarFloat(getField(options, 0, "tauCP"), -61);
		w_vec.sigmaCP = getScalarFloat(getField(options, 0, "sigmaCP"), -62);
		w_vec.sigma2CP = getScalarFloat(getField(options, 0, "sigmaCP"), -62);
		w_vec.thetaCP = getScalarFloat(getField(options, 0, "thetaCP"), -63);
		w_vec.powerIterations = getScalarUInt32(getField(options, 0, "poowerIterations"), -63);
		w_vec.alphaCPTV = getScalarFloat(getField(options, 0, "beta"), -63);
	}
	w_vec.derivType = getScalarUInt32(getField(options, 0, "derivativeType"), -63);
}

// Obtain the reconstruction methods used
void get_rec_methods(const mxArray * options, RecMethods &MethodList) {
	// Non-MAP/prior or single prior algorithms
	//MethodList.MLEM = getScalarBool(getField(options, 0, "MLEM"), -61);
	MethodList.OSEM = getScalarBool(getField(options, 0, "OSEM"), -61);
	MethodList.RAMLA = getScalarBool(getField(options, 0, "RAMLA"), -61);
	MethodList.MRAMLA = getScalarBool(getField(options, 0, "MRAMLA"), -61);
	MethodList.ROSEM = getScalarBool(getField(options, 0, "ROSEM"), -61);
	MethodList.RBI = getScalarBool(getField(options, 0, "RBI"), -61);
	MethodList.DRAMA = getScalarBool(getField(options, 0, "DRAMA"), -61);
	MethodList.COSEM = getScalarBool(getField(options, 0, "COSEM"), -61);
	MethodList.ECOSEM = getScalarBool(getField(options, 0, "ECOSEM"), -61);
	MethodList.ACOSEM = getScalarBool(getField(options, 0, "ACOSEM"), -61);
	MethodList.LSQR = getScalarBool(getField(options, 0, "LSQR"), -61);
	MethodList.CGLS = getScalarBool(getField(options, 0, "CGLS"), -61);
	MethodList.CPLS = getScalarBool(getField(options, 0, "CPLS"), -61);
	MethodList.CPTV = getScalarBool(getField(options, 0, "CPTV"), -61);
	MethodList.CPTVKL = getScalarBool(getField(options, 0, "CPTVKL"), -61);
	if (MethodList.LSQR || MethodList.CGLS)
		MethodList.initAlg = true;

	// Priors
	MethodList.MRP = getScalarBool(getField(options, 0, "MRP"), -61);
	MethodList.Quad = getScalarBool(getField(options, 0, "quad"), -61);
	MethodList.Huber = getScalarBool(getField(options, 0, "Huber"), -61);
	MethodList.L = getScalarBool(getField(options, 0, "L"), -61);
	MethodList.FMH = getScalarBool(getField(options, 0, "FMH"), -61);
	MethodList.WeightedMean = getScalarBool(getField(options, 0, "weighted_mean"), -61);
	MethodList.TV = getScalarBool(getField(options, 0, "TV"), -61);
	MethodList.AD = getScalarBool(getField(options, 0, "AD"), -61);
	MethodList.APLS = getScalarBool(getField(options, 0, "APLS"), -61);
	MethodList.TGV = getScalarBool(getField(options, 0, "TGV"), -61);
	MethodList.NLM = getScalarBool(getField(options, 0, "NLM"), -61);
	MethodList.RDP = getScalarBool(getField(options, 0, "RDP"), -61);

	// MAP/prior-based algorithms
	//MethodList.OSLMLEM = getScalarBool(getField(options, 0, "OSL_MLEM"), -61);
	MethodList.OSLOSEM = getScalarBool(getField(options, 0, "OSL_OSEM"), -61);
	MethodList.BSREM = getScalarBool(getField(options, 0, "BSREM"), -61);
	MethodList.MBSREM = getScalarBool(getField(options, 0, "MBSREM"), -61);
	MethodList.ROSEMMAP = getScalarBool(getField(options, 0, "ROSEM_MAP"), -61);
	MethodList.RBIOSL = getScalarBool(getField(options, 0, "OSL_RBI"), -61);
	MethodList.OSLCOSEM = getScalarUInt32(getField(options, 0, "OSL_COSEM"), -61);
	MethodList.PKMA = getScalarBool(getField(options, 0, "PKMA"), -61);
	//MethodList.CPTV = getScalarBool(getField(options, 0, "CPTV"), -61);

	// Whether MAP/prior-based algorithms are used
	MethodList.MAP = getScalarBool(getField(options, 0, "MAP"), -61);

	// Custom prior
	MethodList.CUSTOM = getScalarBool(getField(options, 0, "custom"), -61);
}

// Transfers the device data to host
// First transfer the ArrayFire arrays from the device to the host pointers pointing to the mxArrays
// Transfer the mxArrays to the cell
void device_to_host_cell(const RecMethods &MethodList, AF_im_vectors & vec, uint32_t & oo, mxArray * cell, Weighting & w_vec,
	const mwSize* dimmi, const uint32_t dim_n, const std::vector<std::vector<float*>>& imEstimates, const scalarStruct& inputScalars)
{
	uint32_t kk;
	uint64_t yy = 0u;
	if (DEBUG) {
		mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
		mexPrintf("vec.im_os.dims(1) = %d\n", vec.im_os.dims(1));
		mexEvalString("pause(.0001);");
	}
	// Transfer data back to host
	for (kk = 0; kk < w_vec.nTot; kk++) {
		mxArray* apu = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		float* apuF = (float*)mxGetSingles(apu);
#else
		float* apuF = (float*)mxGetData(apu);
#endif
		//imEstimates[kk].host(apuF);
		if (inputScalars.saveIter) {
			std::copy(imEstimates[kk].begin(), imEstimates[kk].end(), &apuF);
		}
		else
			vec.im_os(af::seq(yy, yy + inputScalars.im_dim - 1u)).host(apuF);
		af::sync();
		mxSetCell(cell, static_cast<mwIndex>(w_vec.rekot[kk]), mxDuplicateArray(apu));
		yy += inputScalars.im_dim;
	}
	if (MethodList.CUSTOM) {
		kk = w_vec.nTot;
		if (MethodList.OSLCOSEM > 0) {
			mxArray* apu = mxCreateNumericArray(dim_n, dimmi, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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
	const bool zero_pad, const uint32_t Nw)
{
	af::array out = im;
	if (zero_pad == 1) {
		af::dtype type = out.type();
		if (Nz == 1) {
			if (out.dims(1) == 1)
				out = moddims(out, Nx, Ny, Nz, Nw);
			//padd = moddims(out, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
			out = af::constant(0, out.dims(0) + 2 * Ndx, out.dims(1) + 2 * Ndy, 1, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))), af::span, af::span) = out;
			//out = out;
		}
		else {
			if (out.dims(2) == 1)
				out = moddims(out, Nx, Ny, Nz, Nw);
			af::array out = af::constant(0, out.dims(0) + 2 * Ndx, out.dims(1) + 2 * Ndy, out.dims(2) + 2 * Ndz, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))),
				static_cast<double>(Ndz) + af::seq(static_cast<double>(out.dims(2))), af::span) = out;
			//out = out;
		}
	}
	else {
		if (out.dims(1) == 1)
			out = moddims(out, Nx, Ny, Nz, Nw);
		if (Ndx > 0)
			out = af::join(0, af::flip(out(af::seq(static_cast<double>(Ndx)), af::span, af::span, af::span), 0), out,
				af::flip(out(af::seq(static_cast<double>(out.dims(0) - Ndx), static_cast<double>(out.dims(0) - 1)), af::span, af::span, af::span), 0));
		if (Ndy > 0)
			out = af::join(1, af::flip(out(af::span, af::seq(static_cast<double>(Ndy)), af::span, af::span), 1), out,
				af::flip(out(af::span, af::seq(static_cast<double>(out.dims(1) - Ndy), static_cast<double>(out.dims(1) - 1)), af::span, af::span), 1));
		if (Nz == 1 || Ndz == 0) {
		}
		else {
			out = af::join(2, af::flip(out(af::span, af::span, af::seq(static_cast<double>(Ndz)), af::span), 2), out,
				af::flip(out(af::span, af::span, af::seq(static_cast<double>(out.dims(2) - Ndz), static_cast<double>(out.dims(2) - 1)), af::span), 2));
		}
		//out = out;
	}
	return out;
}

af::array EM(const af::array &im, const af::array &Summ, const af::array &rhs)
{
	//const af::array joku = im / (Summ);
	//mexPrintf("joku = %f\n", af::sum<float>(joku));
	//mexEvalString("pause(.0001);");
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

af::array MBSREM(const af::array & im, const af::array & rhs, const float U, const af::array & D, const float* lam, const uint32_t iter, 
	const float beta, const af::array &dU, const af::array & Summ, const scalarStruct inputScalars)
{
	//af::array UU = af::constant(0.f, inputScalars.im_dim);
	af::array output;
	const af::array pp = im < (U / 2.f);
	//UU(pp) = im(pp);
	af::array UU = im / D;
	UU(!pp) = (U - im(!pp)) / (D(!pp));
	if (beta == 0.f)
		output = im + lam[iter] * UU * (rhs - Summ);
	else
		output = im + lam[iter] * UU * (rhs - beta * dU - Summ);
	output(output < inputScalars.epps) = inputScalars.epps;
	output(output >= U) = U - inputScalars.epps;
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
	//return ((1.f - lam[iter] * Summ) * im + lam[iter] * im * rhs);
	return (im + lam[iter] * im * rhs);
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
		output += Summa * (im / D) * (rhs);
	}
	else {
		const float Summa = 1.f / af::max<float>((Summ + beta * dU) / (D + beta * dU));
		output += Summa * (im / (D + beta * dU)) * (rhs - beta * dU);
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

af::array PKMA(const af::array& im, const af::array& Summ, const af::array& rhs, Weighting& w_vec, const af::array& D, 
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const float epps, const float beta, const af::array& dU) {
	af::array S;
	if (w_vec.IEM)
		S = af::max(af::join(1, af::constant(1.f, im.dims(0)) * w_vec.nuIEM, w_vec.preRef, im), 1) / D;
	else
		S = (im + epps)  / D;
	const af::array im_ = im;
	af::array im_apu = im - w_vec.lambda_PKMA[iter] * S * (rhs + beta * dU);
	im_apu(im_apu < epps) = epps;
	uint32_t ind = iter * subsets + osa_iter;
	im_apu = (1.f - w_vec.alpha_PKMA[ind]) * im_ + w_vec.alpha_PKMA[ind] * (w_vec.sigma_PKMA[ind] * im_apu);
	return im_apu;
}

void LSQR(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t iter, AF_im_vectors& vec) {
	if (iter == 0)
		vec.wLSQR = im;
	im = rhs - w_vec.betaLSQR * im;
	w_vec.alphaLSQR = af::norm(im);
	im = im / w_vec.alphaLSQR;
	const float rho_ = sqrt(w_vec.rhoLSQR * w_vec.rhoLSQR + w_vec.betaLSQR * w_vec.betaLSQR);
	const float c = w_vec.rhoLSQR / rho_;
	const float s = w_vec.betaLSQR / rho_;
	w_vec.thetaLSQR = s * w_vec.alphaLSQR;
	w_vec.rhoLSQR = -c * w_vec.alphaLSQR;
	const float phi_ = c * w_vec.phiLSQR;
	w_vec.phiLSQR = s * w_vec.phiLSQR;
	af::array fApu = (phi_ / rho_) * vec.wLSQR + vec.fLSQR;
	vec.fLSQR = fApu.copy();
	af::array wApu = im - (w_vec.thetaLSQR / rho_) * vec.wLSQR;
	vec.wLSQR = wApu.copy();
	if (iter == inputScalars.Niter - 1)
		im = vec.fLSQR;
}

void CGLS(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t iter, AF_im_vectors& vec) {
	//float gamma_ = af::norm(rhs);
	//gamma_ *= gamma_;
	const float gamma_ = af::sum<float>(rhs * rhs);
	const float beta = gamma_ / w_vec.gammaCGLS;
	af::array fApu = vec.fCGLS + w_vec.alphaCGLS * im;
	vec.fCGLS = fApu.copy();
	if (iter == inputScalars.Niter - 1)
		im = vec.fCGLS;
	else
		im = rhs + beta * im;
	w_vec.gammaCGLS = gamma_;
}

void CPLS(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec) {
	const af::array uPrev = vec.uCP;
	const af::array uApu = vec.uCP - w_vec.tauCP * rhs;
	vec.uCP = uApu.copy();
	im = vec.uCP + w_vec.thetaCP * (vec.uCP - uPrev);
}

void computeGradient(const af::array& im, const scalarStruct& inputScalars, af::array& f, af::array& g, af::array& h, const int type) {
	//im = af::moddims(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//f = af::moddims(f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//g = af::moddims(g, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//h = af::moddims(h, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array f = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array g = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array h = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	// 1st order forward differences
	if (type == 0) {
		f(af::seq(0, af::end - 1LL), af::span, af::span) = af::diff1(im);
		f(af::end, af::span, af::span) = -1.f * im(af::end, af::span, af::span);
		g(af::span, af::seq(0, af::end - 1LL), af::span) = af::diff1(im, 1);
		g(af::span, af::end, af::span) = -1.f * im(af::span, af::end, af::span);
		h(af::span, af::span, af::seq(0, af::end - 1LL)) = af::diff1(im, 2);
		h(af::span, af::span, af::end) = -1.f * im(af::span, af::span, af::end);
	}
	// 1st order bacward differences
	else if (type == 1) {
		f(af::seq(1, af::end), af::span, af::span) = -af::diff1(im);
		f(0, af::span, af::span) = -1.f * im(0, af::span, af::span);
		g(af::span, af::seq(1, af::end), af::span) = -af::diff1(im, 1);
		g(af::span, 0, af::span) = -1.f * im(af::span, 0, af::span);
		h(af::span, af::span, af::seq(1, af::end)) = -af::diff1(im, 2);
		h(af::span, af::span, 0) = -1.f * im(af::span, af::span, 0);
	}
	// 1st order central differences
	else {
		f = (shift(im, -1) - shift(im, 1)) * .5f;
		f(0, af::span, af::span) = im(1, af::span, af::span) - im(0, af::span, af::span);
		f(af::end, af::span, af::span) = im(af::end, af::span, af::span) - im(af::end - 1, af::span, af::span);
		g = (shift(im, 0, -1) - shift(im, 0, 1)) * .5f;
		g(af::span, 0, af::span) = im(af::span, 1, af::span) - im(af::span, 0, af::span);
		g(af::span, af::end, af::span) = im(af::span, af::end, af::span) - im(af::span, af::end - 1, af::span);
		h = (shift(im, 0, 0, -1) - shift(im, 0, 0, 1)) * .5f;
		h(af::span, af::span, 0) = im(af::span, af::span, 1) - im(af::span, af::span, 0);
		h(af::span, af::span, af::end) = im(af::span, af::span, af::end) - im(af::span, af::span, af::end - 1);
	}
	f = af::flat(f);
	g = af::flat(g);
	h = af::flat(h);
	//return af::join(0, f, g, h);
}

void computeDivergence(af::array rhs, const af::array& im, const scalarStruct& inputScalars, const int type) {
	rhs = af::moddims(rhs, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array output = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array f = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array g = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	//af::array h = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	// 1st order transpose of forward differences
	if (type == 0) {
		rhs(af::seq(1, af::end), af::span, af::span) += -af::diff1(im(af::span, af::span, af::span, 0));
		rhs(0, af::span, af::span) += -1.f * im(0, af::span, af::span, 0);
		rhs(af::span, af::seq(1, af::end), af::span) += -af::diff1(im(af::span, af::span, af::span, 1), 1);
		rhs(af::span, 0, af::span) += -1.f * im(af::span, 0, af::span, 1);
		rhs(af::span, af::span, af::seq(1, af::end)) += -af::diff1(im(af::span, af::span, af::span, 2), 2);
		rhs(af::span, af::span, 0) += -1.f * im(af::span, af::span, 0, 2);
		//f(af::seq(1, af::end), af::span, af::span) = -af::diff1(im(af::span, af::span, af::span, 0));
		//f(0, af::span, af::span) = -1.f * im(0, af::span, af::span, 0);
		//g(af::span, af::seq(1, af::end), af::span) = -af::diff1(im(af::span, af::span, af::span, 1), 1);
		//g(af::span, 0, af::span) = -1.f * im(af::span, 0, af::span, 1);
		//h(af::span, af::span, af::seq(1, af::end)) = -af::diff1(im(af::span, af::span, af::span, 2), 2);
		//h(af::span, af::span, 0) = -1.f * im(af::span, af::span, 0, 2);
	}
	// 1st order transpose of bacward differences
	else if (type == 1) {
		//f(af::seq(0, inputScalars.Nx - 2u), af::span, af::span) = af::diff1(im(af::span, af::span, af::span, 0));
		//f(af::end, af::span, af::span) = -1.f * im(af::end, af::span, af::span, 0);
		//g(af::span, af::seq(0, inputScalars.Ny - 2u), af::span) = af::diff1(im(af::span, af::span, af::span, 1), 1);
		//g(af::span, af::end, af::span) = -1.f * im(af::span, af::end, af::span, 1);
		//h(af::span, af::span, af::seq(0, inputScalars.Nz - 2u)) = af::diff1(im(af::span, af::span, af::span, 2), 2);
		//h(af::span, af::span, af::end) = -1.f * im(af::span, af::span, af::end, 2);
	}
	// 1st order transpose of central differences
	else {
		//f = (shift(im(af::span, af::span, af::span, 0), 1) - shift(im(af::span, af::span, af::span, 0), -1)) * .5f;
		//f(0, af::span, af::span) = im(0, af::span, af::span, 0) - im(1, af::span, af::span, 0);
		//f(af::end, af::span, af::span) = im(af::end - 1, af::span, af::span, 0) - im(af::end, af::span, af::span, 0);
		//g = (shift(im(af::span, af::span, af::span, 1), 0, 1) - shift(im(af::span, af::span, af::span, 1), 0, -1)) * .5f;
		//g(af::span, 0, af::span) = im(af::span, 0, af::span, 1) - im(af::span, 1, af::span, 1);
		//g(af::span, af::end, af::span) = im(af::span, af::end - 1, af::span, 1) - im(af::span, af::end, af::span, 1);
		//h = (shift(im(af::span, af::span, af::span, 2), 0, 0, 1) - shift(im(af::span, af::span, af::span, 2), 0, 0, -1)) * .5f;
		//h(af::span, af::span, 0) = im(af::span, af::span, 0, 2) - im(af::span, af::span, 1, 2);
		//h(af::span, af::span, af::end) = im(af::span, af::span, af::end - 1, 2) - im(af::span, af::span, af::end, 2);
	}
	//f = af::flat(f);
	//g = af::flat(g);
	//h = af::flat(h);
	//return af::flat(f + g + h);
	//return af::flat(f);
	//return af::flat(output);
}

void computeSecondOrderGradient(af::array& f, af::array& g, af::array& h, const int type) {
	// Second order derivatives
	if (type == 0) {
		f = af::shift(f, 1);
		g = af::shift(g, 0, 1);
		h = af::shift(h, 0, 0, 1);
	}
	else if (type == 1) {
		f = af::shift(f, -1);
		g = af::shift(g, 0, -1);
		h = af::shift(h, 0, 0, -1);
	}
	else {
		f = (shift(f, -1) - shift(f, 1)) * .5f;
		g = (shift(g, 0, -1) - shift(g, 0, 1)) * .5f;
		h = (shift(h, 0, 0, -1) - shift(h, 0, 0, 1)) * .5f;
	}
}

void CPTV(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj) {
	proj.CPTVHelperq(vec.qCPTV, w_vec.alphaCPTV);
	proj.CPTVDiv(vec.qCPTV, rhs, inputScalars);
	CPLS(im, rhs, inputScalars, w_vec, vec);
	proj.CPTVGrad(im, vec.qCPTV, inputScalars, w_vec.sigma2CP);
}

af::array MRP(const af::array& im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const scalarStruct& inputScalars,
	const af::array& offsets, const bool med_no_norm, ProjectorClass& proj) {
	af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, medx, medy, medz);
	const af::dim4 dimmi(padd.dims(0), padd.dims(1), padd.dims(2));
	af::array grad = af::constant(0.f, dimmi);
	padd = af::flat(padd);
	grad = af::flat(grad);
	proj.computeMRP(padd, grad, inputScalars);
	grad = af::moddims(grad, dimmi);
	grad = grad(af::seq(medx, inputScalars.Nx + medx - 1), af::seq(medy, inputScalars.Ny + medy - 1), 
		af::seq(medz, inputScalars.Nz + medz - 1));
	grad = af::flat(grad) + inputScalars.epps;
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

af::array Quadratic_prior(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array &offsets, const af::array &weights_quad)
{
	const af::array apu_pad = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz);
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
	if (Ndz == 0 || inputScalars.Nz == 1) {
		grad = af::convolve2(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx + Ndx - 1), af::seq(Ndy, inputScalars.Ny + Ndy - 1), af::span);
	}
	else {
		grad = af::convolve3(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx + Ndx - 1), af::seq(Ndy, inputScalars.Ny + Ndy - 1), af::seq(Ndz, inputScalars.Nz + Ndz - 1));
	}
	grad = af::flat(grad);
	return grad;
}

af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array& offsets, const af::array& weights_huber, const float delta)
{
	af::array grad = Quadratic_prior(im, Ndx, Ndy, Ndz, inputScalars, inffi, offsets, weights_huber);
	if (af::sum<dim_t>(delta >= af::abs(af::flat(grad))) == grad.elements() && af::sum<int>(af::flat(grad)) != 0)
		mexPrintf("Delta value of Huber prior larger than all the pixel difference values\n");
	grad(grad > delta) = delta;
	grad(grad < -delta) = -delta;
	return grad;
}

af::array FMH(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array& offsets, const af::array& fmh_weights, const bool med_no_norm, const uint32_t alku_fmh)
{
	af::array grad;
	af::array indeksi1;
	const af::array padd = af::flat(padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz));
	uint32_t luup;
	if (inputScalars.Nz == 1 || Ndz == 0) {
		grad = af::constant(0.f, inputScalars.im_dim, 5);
		luup = 4;
	}
	else {
		grad = af::constant(0.f, inputScalars.im_dim, 14);
		luup = 13;
	}
	for (uint32_t ii = 0; ii < luup; ii++) {
		indeksi1 = af::flat(offsets(af::span, af::seq(Ndx * ii, offsets.dims(1) - Ndx * (ii) - 1, alku_fmh / Ndx - ii)));
		af::array apu_pad = af::moddims(padd(indeksi1 + 0), inputScalars.im_dim, fmh_weights.dims(0));
		//grad(af::span, ii) = af::sum(af::batchFunc(apu_pad, af::transpose(fmh_weights(af::span, ii)), batchMul), 1);
		grad(af::span, ii) = af::matmul(apu_pad, fmh_weights(af::span, ii));
	}
	indeksi1 = offsets.col(alku_fmh);
	grad(af::span, af::end) = padd(indeksi1 + 0U);
	grad = af::median(grad, 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

af::array L_filter(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& offsets, const af::array& a_L, const bool med_no_norm)
{
	af::array grad;
	af::array apu_pad = af::flat(padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz));
	apu_pad = apu_pad(af::flat(offsets));
	apu_pad = af::sort(af::moddims(apu_pad, inputScalars.im_dim, a_L.dims(0)), 1);
	grad = af::sum(af::batchFunc(apu_pad, af::transpose(a_L), batchMul), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

af::array Weighted_mean(const af::array & im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& weighted_weights, const bool med_no_norm, const uint32_t mean_type, const float w_sum)
{
	af::array grad = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	const float wsum = af::sum<float>(af::flat(weighted_weights));
	if (mean_type == 1U) {
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz == 1)
			grad = af::convolve2(padd, weighted_weights / wsum);
		else
			grad = af::convolve3(padd, weighted_weights / wsum);
	}
	else if (mean_type == 2U) {
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz == 1)
			grad = 1.f / af::convolve2(1.f / padd, weighted_weights / wsum);
		else
			grad = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
	}
	else if (mean_type == 3U) {
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz == 1)
			grad = af::exp(af::convolve2(af::log(padd), weighted_weights / wsum));
		else
			grad = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
	}
	else if (mean_type == 4U) {
		grad = af::constant(0.f, im.dims(0));
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::convolve3(padd, weighted_weights / wsum);
		//mexPrintf("Nx + Ndx * 4 - 1 = %u\n", Nx + Ndx * 4 - 1);
		if (Ndz == 0 || inputScalars.Nz == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx - 1), af::seq(kk, kk + inputScalars.Ny - 1), af::seq(ll, ll + inputScalars.Nz - 1));
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
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
		if (Ndz == 0 || inputScalars.Nz == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx - 1), af::seq(kk, kk + inputScalars.Ny - 1), af::seq(ll, ll + inputScalars.Nz - 1));
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
		af::array padd = padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
		if (Ndz == 0 || inputScalars.Nz == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx - 1), af::seq(kk, kk + inputScalars.Ny - 1), af::seq(ll, ll + inputScalars.Nz - 1));
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
		if (Ndz == 0 || inputScalars.Nz == 1) {
			grad = grad(af::seq(Ndx, inputScalars.Nx + Ndx - 1), af::seq(Ndy, inputScalars.Ny + Ndy - 1), af::span);
		}
		else {
			grad = grad(af::seq(Ndx, inputScalars.Nx + Ndx - 1), af::seq(Ndy, inputScalars.Ny + Ndy - 1), af::seq(Ndz, inputScalars.Nz + Ndz - 1));
		}
		grad = af::flat(grad);
		if (med_no_norm)
			grad = im - grad;
		else
			grad = (im - grad) / (grad + inputScalars.epps);
	}
	return grad;
}

af::array AD(const af::array & im, const scalarStruct& inputScalars, const float TimeStepAD, const float KAD, const uint32_t NiterAD,
	const af_flux_function FluxType, const af_diffusion_eq DiffusionType, const bool med_no_norm)
{
	const af::array padd = af::moddims(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	af::array grad = af::anisotropicDiffusion(padd, TimeStepAD, KAD, NiterAD, FluxType, DiffusionType);
	grad = af::flat(grad);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}


// Compute the TV prior
af::array TVprior(const scalarStruct& inputScalars, const TVdata& S, const af::array& ima, const uint32_t TVtype,
	const Weighting& w_vec, const af::array& offsets) {
	af::array gradi;

	if (TVtype != 3U) {
		const af::array im = af::moddims(ima, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		// 1st order differentials
		af::array g = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		af::array f = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		af::array h = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		//f(af::seq(0, inputScalars.Nx - 2u), af::span, af::span) = -af::diff1(im);
		////f(Nx - 1u, af::span, af::span) = im(im.dims(0) - 1ULL, af::span, af::span) - im(0, af::span, af::span);
		//f(af::end, af::span, af::span) = f(inputScalars.Nx - 2u, af::span, af::span) * -1.f;
		//g(af::span, af::seq(0, inputScalars.Ny - 2u), af::span) = -af::diff1(im, 1);
		////g(af::span, Ny - 1u, af::span) = im(af::span, im.dims(1) - 1ULL, af::span) - im(af::span, 0, af::span);
		//g(af::span, af::end, af::span) = g(af::span, inputScalars.Ny - 2u, af::span) * -1.f;
		//h(af::span, af::span, af::seq(0, inputScalars.Nz - 2u)) = -af::diff1(im, 2);
		////h(af::span, af::span, Nz - 1u) = im(af::span, af::span, im.dims(2) - 1ULL) - im(af::span, af::span, 0);
		//h(af::span, af::span, af::end) = h(af::span, af::span, inputScalars.Nz - 2u) * -1.f;

		//g = af::flat(g);
		//f = af::flat(f);
		//h = af::flat(h);
		computeGradient(im, inputScalars, f, g, h, w_vec.derivType);
		af::array pval, apu1, apu2, apu3;

			// If anatomical prior is used
			if (S.TV_use_anatomical || TVtype == 5U) {
				if (TVtype == 1U) {
					pval = af::sqrt(S.s1 * af::pow(f, 2.) + S.s5 * af::pow(g, 2.) + S.s9 * af::pow(h, 2.) + S.s4 * f * g + S.s7 * f * h + S.s2 * f * g + S.s8 * h * g + S.s3 * f * h + S.s6 * h * g + S.TVsmoothing);
					apu1 = 0.5f * (2.f * S.s1 * f + S.s4 * g + S.s7 * h + S.s2 * g + S.s3 * h) / pval;
					apu2 = 0.5f * (2.f * S.s5 * g + S.s4 * f + S.s2 * f + S.s8 * h + S.s6 * h) / pval;
					apu3 = 0.5f * (2.f * S.s9 * h + S.s8 * g + S.s6 * g + S.s7 * f + S.s3 * f) / pval;
					gradi = 0.5f * (2.f * S.s1 * f + 2.f * S.s5 * g + 2.f * S.s9 * h + S.s4 * f + S.s2 * f + S.s8 * h + S.s6 * h + S.s4 * g + S.s7 * h + S.s2 * g + S.s3 * h
						+ S.s8 * g + S.s6 * g + S.s7 * f + S.s3 * f) / pval;
				}
				else if (TVtype == 2U) {
					const af::array reference_image = af::moddims(S.reference_image, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					af::array gp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					af::array fp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					af::array hp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					computeGradient(reference_image, inputScalars, fp, gp, hp, w_vec.derivType);
					//fp(af::seq(0, inputScalars.Nx - 2u), af::span, af::span) = -af::diff1(reference_image);
					////fp(Nx - 1u, af::span, af::span) = S.reference_image(im.dims(0) - 1ULL, af::span, af::span) - S.reference_image(0, af::span, af::span);
					//fp(af::end, af::span, af::span) = fp(inputScalars.Nx - 2u, af::span, af::span) * -1.f;
					//gp(af::span, af::seq(0, inputScalars.Ny - 2u), af::span) = -af::diff1(reference_image, 1);
					////gp(af::span, Ny - 1u, af::span) = S.reference_image(af::span, im.dims(1) - 1ULL, af::span) - S.reference_image(af::span, 0, af::span);
					//gp(af::span, af::end, af::span) = gp(af::span, inputScalars.Ny - 2u, af::span) * -1.f;
					//hp(af::span, af::span, af::seq(0, inputScalars.Nz - 2u)) = -af::diff1(reference_image, 2);
					////hp(af::span, af::span, Nz - 1u) = S.reference_image(af::span, af::span, im.dims(2) - 1ULL) - S.reference_image(af::span, af::span, 0);
					//hp(af::span, af::span, af::end) = hp(af::span, af::span, inputScalars.Nz - 2u) * -1.f;

					//gp = af::flat(gp);
					//fp = af::flat(fp);
					//hp = af::flat(hp);

					pval = af::sqrt(af::pow(f, 2.) + af::pow(g, 2.) + af::pow(h, 2.) + S.T * (af::pow(fp, 2.) + af::pow(gp, 2.) + af::pow(hp, 2.)) + S.TVsmoothing);
					apu1 = f / pval;
					apu2 = g / pval;
					apu3 = h / pval;
					gradi = (f + g + h) / pval;
				}
				// For APLS
				else if (TVtype == 5U) {
					af::array gp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					af::array fp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					af::array hp = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
					computeGradient(S.APLSReference, inputScalars, fp, gp, hp, w_vec.derivType);
					//fp(af::seq(0, inputScalars.Nx - 2u), af::span, af::span) = -af::diff1(S.APLSReference);
					////fp(Nx - 1u, af::span, af::span) = S.APLSReference(im.dims(0) - 1ULL, af::span, af::span) - S.APLSReference(0, af::span, af::span);
					//fp(af::end, af::span, af::span) = fp(inputScalars.Nx - 2u, af::span, af::span) * -1.f;
					//gp(af::span, af::seq(0, inputScalars.Ny - 2u), af::span) = -af::diff1(S.APLSReference, 1);
					////gp(af::span, Ny - 1u, af::span) = S.APLSReference(af::span, im.dims(1) - 1ULL, af::span) - S.APLSReference(af::span, 0, af::span);
					//gp(af::span, af::end, af::span) = gp(af::span, inputScalars.Ny - 2u, af::span) * -1.f;
					//hp(af::span, af::span, af::seq(0, inputScalars.Nz - 2u)) = -af::diff1(S.APLSReference, 2);
					////hp(af::span, af::span, Nz - 1u) = S.APLSReference(af::span, af::span, im.dims(2) - 1ULL) - S.APLSReference(af::span, af::span, 0);
					//hp(af::span, af::span, af::end) = hp(af::span, af::span, inputScalars.Nz - 2u) * -1.f;

					fp = af::flat(fp) + inputScalars.epps;
					gp = af::flat(gp) + inputScalars.epps;
					hp = af::flat(hp) + inputScalars.epps;

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
					gradi = (f - (apu * epsilon(af::span, 0)) + g - (apu * epsilon(af::span, 1)) + h - (apu * epsilon(af::span, 2))) / pval;
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
				gradi = apu1 + apu2 + apu3;
			}
			apu1 = af::moddims(apu1, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
			apu2 = af::moddims(apu2, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
			apu3 = af::moddims(apu3, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
			gradi = af::moddims(gradi, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
			computeSecondOrderGradient(apu1, apu2, apu3, w_vec.derivType);
			if (w_vec.derivType < 2)
				gradi = apu1 + apu2 + apu3 - gradi;
			else
				gradi = apu1 + apu2 + apu3;
			gradi = af::flat(gradi);
			gradi = gradi + 2.f * S.tau * af::min<float>(af::flat(ima));
	}
	else {
		if (S.TV_use_anatomical) {
			af::array padd = af::flat(padding(ima, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd = padd(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
			padd = af::moddims(padd, ima.dims(0), padd.dims(0) / ima.dims(0));
			af::array padd2 = af::flat(padding(S.reference_image, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd2 = padd2(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
			padd2 = af::moddims(padd2, ima.dims(0), padd2.dims(0) / ima.dims(0));
			gradi = af::sum(af::batchFunc(af::batchFunc(ima, padd, batchMinus) / std::pow(S.C, 2.f) * (1.f / af::sqrt(1.f + af::pow(af::batchFunc(ima, padd, batchMinus) / S.C, 2.)
				+ af::pow(af::batchFunc(S.reference_image, padd2, batchMinus) / S.T, 2.))), w_vec.weights_TV.T(), batchMul), 1);
		}
		else {
			af::array padd = af::flat(padding(ima, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz));
			padd = padd(af::join(1, af::flat(offsets(af::span, af::seq(0, (w_vec.dimmu - 1) / 2 - 1))), af::flat(offsets(af::span, af::seq((w_vec.dimmu - 1) / 2 + 1, offsets.dims(1) - 1ULL)))));
			padd = af::moddims(padd, ima.dims(0), padd.dims(0) / ima.dims(0));
			gradi = af::sum(af::batchFunc(af::batchFunc(ima, padd, batchMinus) / std::pow(S.C, 2.f) * (1.f / af::sqrt(1.f + af::pow(af::batchFunc(ima, padd, batchMinus) / S.C, 2.))),
				w_vec.weights_TV.T(), batchMul), 1);
		}
	}

	return gradi;
}

af::array TGV(const af::array & im, const scalarStruct& inputScalars, const uint32_t maxits, const float alpha, const float beta)
{
	af::array grad = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array u = grad;
	grad = af::flat(grad);
	const af::array imi = af::moddims(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	const float sigma = 1.f / 16.f;
	const float tau = 1.f / 8.f;
	af::array p1 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array p2 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array p3 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array v1 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array v2 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array v3 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array q1 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array q2 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array q3 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);
	af::array q4 = af::constant(0.f, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, f32);

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
		p1 = af::moddims(p1, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		p2 = af::flat(eta2) / (apu);
		p2 = af::moddims(p2, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		p3 = af::flat(eta3) / (apu);
		p3 = af::moddims(p3, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
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
		q1 = af::moddims(q1, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		q2 = af::flat(eta2) / (apu);
		q2 = af::moddims(q2, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		q3 = af::flat(eta3) / (apu);
		q3 = af::moddims(q3, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		q4 = af::flat(eta4) / (apu);
		q4 = af::moddims(q4, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);

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
		u = af::moddims(u, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		vb1 = af::moddims(2.f * v1 - v1old, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		vb2 = af::moddims(2.f * v2 - v2old, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
		vb3 = af::moddims(2.f * v3 - v3old, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);

	}

	return -grad;
}

af::array RDP(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_RDP,
	const float gamma, const af::array& offsets, const uint32_t inffi)
{
	const af::array im_apu = af::flat(padding(im, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, Ndx, Ndy, Ndz));

	const af::array indeksi2 = af::join(1, offsets.cols(0, inffi - 1), offsets.cols(inffi + 1, af::end));
	af::array grad = af::batchFunc(im, af::moddims(im_apu(indeksi2), inputScalars.im_dim, weights_RDP.dims(0)), batchDiv);
	af::array apu = grad + 1.f + gamma * af::abs(grad - 1.f);
	grad = af::matmul((((grad - 1.f) * (gamma * af::abs(grad - 1.f) + grad + 3.f))  / (apu * apu)), weights_RDP);
	grad = af::flat(grad);
	return grad;
}


af::array computeConvolution(const af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec,
	const uint32_t nRekos) {
	af::array apu = af::moddims(vec, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, nRekos);
	//for (uint32_t ff = 0U; ff < inputScalars.nRekos2; ff++) {
		padding(apu, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1, 
			false, inputScalars.nRekos2);
		apu = af::convolve3(apu, g);
	//	apu2 = apu2(af::seq(w_vec.g_dim_x + 1, inputScalars.Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, inputScalars.Ny + w_vec.g_dim_y), 
	//		af::seq(w_vec.g_dim_z + 1, inputScalars.Nz + w_vec.g_dim_z));
	//	apu(af::span, af::span, af::span, ff) = apu2;
	//}
	return af::flat(apu);
}

void deblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec) {
	//uint32_t it = 0U;
	//if (inputScalars.saveIter)
	//	it = iter + 1U;
	af::array jelppi = padding(vec, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	af::array apu = padding(vec, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	for (uint32_t kk = 0U; kk < w_vec.deblur_iterations; kk++) {
		af::array apu2 = convolve3(jelppi, g) + inputScalars.epps;
		apu2 = apu2(af::seq(w_vec.g_dim_x + 1, inputScalars.Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, inputScalars.Ny + w_vec.g_dim_y), 
			af::seq(w_vec.g_dim_z + 1, inputScalars.Nz + w_vec.g_dim_z));
		apu2 = padding(apu2, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
		jelppi *= af::convolve3(apu / apu2, g);
		jelppi = jelppi(af::seq(w_vec.g_dim_x + 1, inputScalars.Nx + w_vec.g_dim_x), af::seq(w_vec.g_dim_y + 1, inputScalars.Ny + w_vec.g_dim_y), 
			af::seq(w_vec.g_dim_z + 1, inputScalars.Nz + w_vec.g_dim_z));
		if (kk < w_vec.deblur_iterations - 1)
			jelppi = padding(jelppi, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, w_vec.g_dim_x + 1, w_vec.g_dim_y + 1, w_vec.g_dim_z + 1);
	}
	vec = af::flat(jelppi);
	//vec = jelppi;
}

//void computeDeblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec,
//	const RecMethods& MethodList, const uint32_t iter) {
//
//	deblur(vec, g, inputScalars, w_vec);
//}



af::array NLM(ProjectorClass& proj, const af::array& im, Weighting& w_vec, const scalarStruct& inputScalars)
{
	int32_t type = 0;
	if (w_vec.NLTV)
		type = 1;
	else if (w_vec.NLM_MRP)
		type = 2;
	else
		type = 0;
	af::array grad = af::constant(0.f, im.elements(), 1);
	af::sync();
	proj.computeNLM(grad, im, inputScalars, w_vec, type);
	grad.unlock();
	im.unlock();
	w_vec.gaussianNLM.unlock();

	if (w_vec.NLM_MRP) {
		grad = im - grad;
	}
	//af::sync();
	return grad;
}

void initializeRHS(AF_im_vectors& vec, const scalarStruct& inputScalars) {
	if (inputScalars.atomic_64bit)
		vec.rhs_os[0] = (af::constant(0LL, static_cast<size_t>(inputScalars.im_dim) * static_cast<size_t>(inputScalars.nRekos), 1, s64));
	else if (inputScalars.atomic_32bit)
		vec.rhs_os[0] = (af::constant(0, static_cast<size_t>(inputScalars.im_dim) * static_cast<size_t>(inputScalars.nRekos), 1, s32));
	else
		vec.rhs_os[0] = (af::constant(0.f, static_cast<size_t>(inputScalars.im_dim) * static_cast<size_t>(inputScalars.nRekos), 1));
}

int initializationStep(Weighting& w_vec, af::array& mData, AF_im_vectors& vec, ProjectorClass& proj, scalarStruct& inputScalars,
	std::vector<int64_t> length, uint64_t m_size, uint64_t st, const RecMethods& MethodList, uint32_t curIter, const af::array& g, 
	af::array& meanBP, std::vector<af::array>& Summ) {

	if (curIter == 0) {
		af::sync();
		int status = 0;
		uint64_t yy = 0u;
		if (MethodList.LSQR) {
			//w_vec.wLSQR = af::constant(0.f, vec.im_os.dims(0));
			//w_vec.vhLSQR = af::constant(0.f, vec.im_os.dims(0));
			vec.fLSQR = vec.im_os.copy();
			w_vec.betaLSQR = af::norm(mData);
			mData = mData / w_vec.betaLSQR;
			//if (DEBUG) {
			//	mexPrintf("!!!!!!!!!!!!!!!!!!!!!!!mData = %f\n", af::sum<float>(mData));
			//	mexEvalString("pause(.0001);");
			//}
			if (inputScalars.projector_type == 6)
				backprojectionSPECT(mData, Summ, w_vec, vec, inputScalars, length[0], 0, 0, 0, 0, 0);
			else
				status = proj.backwardProjection(vec, inputScalars, w_vec, mData, 0, length, st, m_size, meanBP);
			if (status != 0) {
				vec.rhs_os[0].unlock();
				mData.unlock();
				return -1;
			}
			af::sync();
			vec.rhs_os[0].unlock();
			mData.unlock();
			af::sync();
			if (inputScalars.atomic_64bit)
				vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH;
			else if (inputScalars.atomic_32bit)
				vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH32;
			if (DEBUG) {
				af_print_mem_info("mem info", -1);
				mexPrintf("!!!!!!!!!!!!!!!!!!!!!!!vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os[0]));
				mexEvalString("pause(.0001);");
			}
			w_vec.alphaLSQR = af::norm(vec.rhs_os[0]);
			vec.im_os = vec.rhs_os[0] / w_vec.alphaLSQR;
			if (inputScalars.use_psf) {
				vec.im_os_blurred = computeConvolution(vec.rhs_os[0], g, inputScalars, w_vec, inputScalars.nRekos2);
				w_vec.alphaLSQR = af::norm(vec.im_os_blurred);
				vec.im_os_blurred = vec.im_os_blurred / w_vec.alphaLSQR;
				vec.wLSQR = vec.im_os_blurred.copy();
			}
			else {
				vec.wLSQR = vec.im_os.copy();
			}
			if (DEBUG) {
				mexPrintf("!!!!!!!!!!!!!!!!!!!!!!!vec.im_os = %f\n", af::sum<float>(vec.im_os));
				mexEvalString("pause(.0001);");
			}
			w_vec.phiLSQR = w_vec.betaLSQR;
			w_vec.rhoLSQR = w_vec.alphaLSQR;
			af::sync();
		}
		else if (MethodList.CGLS) {
			vec.rCGLS = mData;
			vec.fCGLS = vec.im_os.copy();
			if (inputScalars.projector_type == 6)
				backprojectionSPECT(mData, Summ, w_vec, vec, inputScalars, length[0], 0, 0, 0, 0, 0);
			else
				status = proj.backwardProjection(vec, inputScalars, w_vec, mData, 0, length, st, m_size, meanBP);
			if (status != 0) {
				vec.rhs_os[0].unlock();
				mData.unlock();
				return -1;
			}
			af::sync();
			vec.rhs_os[0].unlock();
			mData.unlock();
			af::sync();
			if (inputScalars.atomic_64bit)
				vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH;
			else if (inputScalars.atomic_32bit)
				vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH32;
			vec.im_os = vec.rhs_os[0];
			if (inputScalars.use_psf)
				vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
			if (DEBUG) {
				af_print_mem_info("mem info", -1);
				mexEvalString("pause(.0001);");
			}
			w_vec.gammaCGLS = af::sum<float>(vec.rhs_os[0] * vec.rhs_os[0]);
			//w_vec.gammaCGLS = af::norm(vec.rhs_os);
			//w_vec.gammaCGLS *= w_vec.gammaCGLS;
		}
		else if (inputScalars.use_psf)
			vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
		if (MethodList.CPLS || MethodList.CPTV || MethodList.CPTVKL || MethodList.CPLSKL) {
			vec.pCP = af::constant(0.f, m_size);
			vec.uCP = af::constant(0.f, inputScalars.im_dim);
			if (MethodList.CPTV || MethodList.CPTVKL)
				vec.qCPTV = af::constant(0.f, static_cast<dim_t>(inputScalars.im_dim) * 3LL);
		}
		if (MethodList.initAlg)
			vec.rhs_os.clear();
			//initializeRHS(vec, inputScalars);
	}
	else if (inputScalars.use_psf) {
		vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
	}
	af::sync();
	af::deviceGC();
	return 0;
}

void forwardProjectionSPECT(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars, 
	const int64_t length, const uint32_t uu) {
	uint32_t u1 = uu;
	const af::array apuArr = af::moddims(vec.im_os, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz);
	for (int kk = 0; kk < length; kk++) {
		af::array kuvaRot = af::rotate(apuArr, w_vec.angles[u1], true, AF_INTERP_BILINEAR);
		kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
		kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
		kuvaRot = kuvaRot(af::span, af::span, af::seq(w_vec.distInt[u1], af::end));
		kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
		kuvaRot = af::reorder(af::sum(kuvaRot, 0), 1, 2, 0);
		fProj(af::span, af::span, kk) = kuvaRot;
		u1++;
	}
	fProj = af::flat(fProj);
	fProj(fProj < inputScalars.epps) = inputScalars.epps;
}

void backprojectionSPECT(af::array& fProj, std::vector<af::array>& Summ, const Weighting& w_vec, AF_im_vectors& vec, 
	const scalarStruct& inputScalars, const int64_t length, const uint32_t uu, const uint32_t osa_iter, const uint32_t iter, 
	const uint8_t compute_norm_matrix, const uint32_t iter0) {
	fProj = af::moddims(fProj, inputScalars.size_x, w_vec.size_y, length);
	af::array apuBP2 = af::constant(0.f, inputScalars.size_x * w_vec.size_y * inputScalars.size_x, length);
	uint32_t u1 = uu;
	for (int kk = 0; kk < length; kk++) {
		af::array apuBP = af::constant(0.f, inputScalars.size_x, inputScalars.size_x, w_vec.size_y);
		af::array kuvaRot = fProj(af::span, af::span, kk);
		kuvaRot = af::reorder(kuvaRot, 1, 0, 2);
		af::eval(kuvaRot);
		kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
		kuvaRot = kuvaRot(af::span, af::span, af::seq(w_vec.distInt[u1], af::end));
		kuvaRot = reorder(kuvaRot, 2, 1, 0);
		apuBP(af::seq(w_vec.distInt[u1], af::end), af::span, af::span) = kuvaRot;
		apuBP = af::rotate(apuBP, -w_vec.angles[u1], true, AF_INTERP_BILINEAR);
		apuBP2(af::span, kk) = af::flat(apuBP);
		u1++;
	}
	vec.rhs_os[0] = af::sum(apuBP2, 1);
	vec.rhs_os[0](vec.rhs_os[0] < inputScalars.epps && vec.rhs_os[0] >= 0.f) = inputScalars.epps;
	if ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1) {
		apuBP2 = af::constant(0.f, inputScalars.size_x * w_vec.size_y * inputScalars.size_x, length);
		u1 = uu;
		for (int kk = 0; kk < length; kk++) {
			af::array apuSumm = af::constant(0.f, inputScalars.size_x, inputScalars.size_x, w_vec.size_y);
			af::array kuvaRot = af::constant(1.f, w_vec.size_y, inputScalars.size_x);
			kuvaRot = af::convolve2(kuvaRot, w_vec.gFilter(af::span, af::span, af::span, u1));
			kuvaRot = kuvaRot(af::span, af::span, af::seq(w_vec.distInt[u1], af::end));
			kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
			apuSumm(af::seq(w_vec.distInt[u1], af::end), af::span, af::span) = kuvaRot;
			apuSumm = af::rotate(apuSumm, -w_vec.angles[u1], true, AF_INTERP_BILINEAR);
			apuBP2(af::span, kk) = af::flat(apuSumm);
			u1++;
		}
		if (compute_norm_matrix == 2) {
			Summ[osa_iter] = af::sum(apuBP2, 1);
			Summ[osa_iter](Summ[osa_iter] < inputScalars.epps) = 1.f;
		}
		else {
			Summ[0] = af::sum(apuBP2, 1);
			Summ[0](Summ[osa_iter] < inputScalars.epps) = 1.f;
		}
	}
}

int computeACOSEMWeight(const scalarStruct& inputScalars, const std::vector<int64_t>& length, float& uu, const uint32_t osa_iter, const af::array& mData, 
	const uint64_t m_size, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t subSum) {
	int status = 0;
	if (inputScalars.projector_type < 4) {
		uu = proj.computeSum(length, osa_iter);
	}
	else
		uu = af::sum<float>(mData);
	af::array outputFP;
	if (inputScalars.projector_type != 6) {
		outputFP = af::constant(0.f, m_size);
		status = proj.update_opencl_inputs(vec, inputScalars);
		if (status != 0)
			return -1;
		af::sync();
		status = proj.forwardProjection(vec, inputScalars, w_vec, outputFP, osa_iter, length, 0, m_size);
		af::sync();
		if (inputScalars.use_psf)
			vec.im_os_blurred.unlock();
		else
			vec.im_os.unlock();
		outputFP.unlock();
		if (status != 0)
			return -1;
	}
	else {
		outputFP = af::constant(0.f, inputScalars.size_x, w_vec.size_y, length[osa_iter]);
		forwardProjectionSPECT(outputFP, w_vec, vec, inputScalars, length[osa_iter], subSum);
	}
	if (inputScalars.CT)
		w_vec.ACOSEM_rhs = af::sum<float>(af::exp(-outputFP));
	else
		w_vec.ACOSEM_rhs = af::sum<float>(outputFP);
	return 0;
}

int powerMethod(scalarStruct& inputScalars, Weighting& w_vec, const std::vector<int64_t>& length, const uint64_t m_size, ProjectorClass& proj, 
	AF_im_vectors& vec, const af::array& g, const RecMethods& MethodList) {
	int status = 0;
	std::vector<af::array> Summ;
	af::array meanBP;
	w_vec.tauCP = 0.f;
	vec.im_os = af::randn(inputScalars.im_dim);
	vec.im_os = vec.im_os / af::norm(vec.im_os);
	af::array outputFP = af::constant(0.f, m_size);
	for (int kk = 0; kk < w_vec.powerIterations; kk++) {
		af_print_mem_info("mem info", -1);
		if (inputScalars.use_psf)
			vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
		af::sync();
		if (inputScalars.projector_type == 6)
			forwardProjectionSPECT(outputFP, w_vec, vec, inputScalars, length[0], 0);
		else
			status = proj.forwardProjection(vec, inputScalars, w_vec, outputFP, 0, length, 0, m_size);
		af::sync();
		if (inputScalars.use_psf)
			vec.im_os_blurred.unlock();
		else
			vec.im_os.unlock();
		outputFP.unlock();
		if (status != 0)
			return -1;
		if (DEBUG) {
			mexPrintf("Power forward projection complete\n");
			mexEvalString("pause(.0001);");
		}
		af::sync();
		if (inputScalars.projector_type == 6)
			backprojectionSPECT(outputFP, Summ, w_vec, vec, inputScalars, length[0], 0, 0, 0, 0, 0);
		else
			status = proj.backwardProjection(vec, inputScalars, w_vec, outputFP, 0, length, 0, m_size, meanBP);
		af::sync();
		vec.rhs_os[0].unlock();
		outputFP.unlock();
		if (status != 0)
			return -1;
		if (inputScalars.use_psf) {
			vec.rhs_os[0] = computeConvolution(vec.rhs_os[0], g, inputScalars, w_vec);
		}
		if (MethodList.CPTV || MethodList.CPTVKL) {

		}
		w_vec.tauCP = af::dot<float>(vec.im_os, vec.rhs_os[0]);
		vec.im_os = vec.rhs_os[0];
		vec.im_os = vec.im_os / af::norm(vec.im_os);
		if (inputScalars.verbose == 2 || DEBUG) {
			mexPrintf("Largest eigenvalue at iteration %d is %f\n", kk, w_vec.tauCP);
			mexEvalString("pause(.0001);");
		}
	}
	w_vec.sigmaCP = 1.f;
	w_vec.sigma2CP = 1.f;
	mexPrintf("Largest eigenvalue is %f\n", w_vec.tauCP);
	mexEvalString("pause(.0001);");
	return 0;
}