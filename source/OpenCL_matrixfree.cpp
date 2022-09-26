/**************************************************************************
* Matrix free computations for OMEGA (Implementation 2).
* This is the main file for the matrix-free computations in OMEGA. In this
* file the MATLAB variables are loaded and either the matrix-free 
* reconstructions are computed or the custom prior reconstruction.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly inaccurate.
* Both functions use ArrayFire functions and thus require the installation
* of the ArrayFire library.
* 
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include "AF_opencl_functions.hpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 40)
		mexErrMsgTxt("Too few input arguments.  There must be at least 40.");
	else if (nrhs > 74)
		mexErrMsgTxt("Too many input arguments.  There can be at most 74.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	int ind = 0;
	scalarStruct inputScalars;
	// Load the input arguments
	// These are mainly same as with implementation 3 (check the comments in OpenCL_matrixfree_multi_gpu.cpp)
	const char *k_path = mxArrayToString(prhs[ind]);
	ind++;

	inputScalars.Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.Nx = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.dx = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.dz = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.by = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.bx = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.bz = getScalarFloat(prhs[ind], ind);
	ind++;

	// Coordinates of the detectors in z-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* z_det = (float*)mxGetSingles(prhs[ind]);
#else
	const float* z_det = (float*)mxGetData(prhs[ind]);
#endif
	const size_t size_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in x-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* x = (float*)mxGetSingles(prhs[ind]);
#else
	const float* x = (float*)mxGetData(prhs[ind]);
#endif
	const size_t size_of_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in y-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* y = (float*)mxGetSingles(prhs[ind]);
#else
	const float* y = (float*)mxGetData(prhs[ind]);
#endif
	ind++;

	inputScalars.dy = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.maxyy = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.maxxx = getScalarFloat(prhs[ind], ind);
	ind++;

	const uint32_t NSinos = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.NSlices = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.size_x = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.zmax = getScalarFloat(prhs[ind], ind);
	ind++;

	const uint32_t TotSinos = getScalarUInt32(prhs[ind], ind);
	ind++;

	const bool verbose = getScalarBool(prhs[ind], ind);
	ind++;
	
	// Detector pair numbers, for raw list-mode data
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const uint16_t* L = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* L = (uint16_t*)mxGetData(prhs[ind]);
#endif
	const size_t numRows = mxGetM(prhs[ind]);
	ind++;

	// Location (ring numbers) of pseudo rings, if present
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const uint32_t* pseudos = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
	const uint32_t* pseudos = (uint32_t*)mxGetData(prhs[ind]);
#endif
	const uint32_t pRows = (uint32_t)mxGetNumberOfElements(prhs[ind]);
	ind++;

	inputScalars.det_per_ring = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Is TOF data used?
	inputScalars.TOF = getScalarBool(prhs[ind], ind);
	ind++;

	// Size of single TOF-subset
	const int64_t TOFSize = getScalarInt64(prhs[ind], ind);
	ind++;

	// Variance of the Gaussian TOF
	inputScalars.sigma_x = getScalarFloat(prhs[ind], ind);
	ind++;

	// Centers of the TOF-bins
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* TOFCenter = (float*)mxGetSingles(prhs[ind]);
#else
	const float* TOFCenter = (float*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Index offset for TOF subsets
	inputScalars.nBins = getScalarInt64(prhs[ind], ind);
	ind++;


	inputScalars.dec = getScalarUInt32(prhs[ind], ind);
	ind++;

	// The device used
	const uint32_t device = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.raw = getScalarUInt8(prhs[ind], ind);
	ind++;

	//af::setDevice(device);

	const char* fileName = mxArrayToString(prhs[ind]);
	ind++;

	const uint32_t type = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.use_psf = getScalarBool(prhs[ind], ind);
	ind++;

	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[ind]);
	ind++;

	inputScalars.bzb = inputScalars.bz + static_cast<float>(inputScalars.Nz) * inputScalars.dz;

	if (type < 2) {

		if (nrhs != 74)
			mexErrMsgTxt("Invalid number of input arguments. There must be 74.");

		// attenuation values
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* atten = (float*)mxGetSingles(prhs[ind]);
#else
		const float* atten = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_atten = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Normalization coefficients
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* norm = (float*)mxGetSingles(prhs[ind]);
#else
		const float* norm = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_norm = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Number of measurements/LORs
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const int64_t* pituus = (int64_t*)mxGetInt64s(prhs[ind]);
#else
		const int64_t* pituus = (int64_t*)mxGetData(prhs[ind]);
#endif
		const size_t nPituus = mxGetNumberOfElements(prhs[ind]);
		ind++;

		inputScalars.attenuation_correction = getScalarUInt32(prhs[ind], ind);
		ind++;

		inputScalars.normalization_correction = getScalarUInt32(prhs[ind], ind);
		ind++;

		const uint32_t Niter = getScalarUInt32(prhs[ind], ind);
		ind++;

		inputScalars.subsets = getScalarUInt32(prhs[ind], ind);
		ind++;

		//const uint8_t* rekot = (uint8_t*)mxGetData(prhs[37]);
		const size_t size_reko = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float epps = getScalarFloat(prhs[ind], ind);
		ind++;

		// Number of voxels the current LOR/ray traverses
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint16_t* lor1 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// XY-indices of the detector coordinates of each LOR
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint32_t* xy_index = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Z-indices of the detector coordinates of each LOR
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint16_t* z_index = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Is any OS-method used
		const bool osem_bool = getScalarBool(prhs[ind], ind);
		ind++;

		inputScalars.tube_width = getScalarFloat(prhs[ind], ind);
		ind++;

		//const float crystal_size_z = getScalarFloat(prhs[ind], ind);
		//ind++;

		// Center coordinates of voxels in the X-dimension
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* x_center = (float*)mxGetSingles(prhs[ind]);
#else
		const float* x_center = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_center_x = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Center coordinates of voxels in the Y-dimension
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* y_center = (float*)mxGetSingles(prhs[ind]);
#else
		const float* y_center = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_center_y = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Center coordinates of voxels in the Z-dimension
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* z_center = (float*)mxGetSingles(prhs[ind]);
#else
		const float* z_center = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_center_z = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Randoms
		const mxArray* sc_ra = prhs[ind];
		ind++;

		// Randoms corrections
		inputScalars.randoms_correction = getScalarUInt32(prhs[ind], ind);
		ind++;

		// The type of projector used (Siddon or orthogonal)
		inputScalars.projector_type = getScalarUInt32(prhs[ind], ind);
		ind++;

		// If true, then the precomputed LOR voxels counts are used
		inputScalars.precompute = getScalarBool(prhs[ind], ind);
		ind++;

		// Number of rays in Siddon
		inputScalars.n_rays = getScalarUInt16(prhs[ind], ind);
		ind++;

		// Number of rays in Siddon (axial)
		inputScalars.n_rays3D = getScalarUInt16(prhs[ind], ind);
		ind++;

		// Crystal pitch in z-direction
		inputScalars.cr_pz = getScalarFloat(prhs[ind], ind);
		ind++;

		const mxArray* options = prhs[ind];
		ind++;

		const bool saveIter = getScalarBool(getField(options, 0, "save_iter"), ind);
		size_t Ni = 0ULL;
		if (saveIter)
			Ni = static_cast<size_t>(Niter);
		const size_t outSize = static_cast<size_t>(inputScalars.Nx) * static_cast<size_t>(inputScalars.Ny) * static_cast<size_t>(inputScalars.Nz);
		const size_t outSize2 = Ni + 1ULL;

		// Implementation 2
		//if (type == 0) {

		//Cell array containing the measurements
		const mxArray* Sin = prhs[ind];
		ind++;

		// Number of time steps
		const uint32_t Nt = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Use 64-bit integer atomic functions if possible
		const bool use_64bit_atomics = getScalarBool(prhs[ind], ind);
		ind++;

		// Number of OS-reconstruction algorithms (including priors)
		const uint32_t n_rekos = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Number of MLEM algorithms (including priors)
		const uint32_t n_rekos_mlem = getScalarUInt32(prhs[ind], ind);
		ind++;

		// What type of reconstruction, needed for the OpenCL kernel
		// E.g. 2 means COSEM (different computations)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint8_t* reko_type = (uint8_t*)mxGetUint8s(prhs[ind]);
#else
		const uint8_t* reko_type = (uint8_t*)mxGetData(prhs[ind]);
#endif
		ind++;


#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint8_t* reko_type_mlem = (uint8_t*)mxGetUint8s(prhs[ind]);
#else
		const uint8_t* reko_type_mlem = (uint8_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Global correction factor
		inputScalars.global_factor = getScalarFloat(prhs[ind], ind);
		ind++;

		inputScalars.bmin = getScalarFloat(prhs[ind], ind);
		ind++;

		inputScalars.bmax = getScalarFloat(prhs[ind], ind);
		ind++;

		inputScalars.Vmax = getScalarFloat(prhs[ind], ind);
		ind++;

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* V = (float*)mxGetSingles(prhs[ind]);
#else
		const float* V = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_V = mxGetNumberOfElements(prhs[ind]);
		ind++;

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const float* gaussian = (float*)mxGetSingles(prhs[ind]);
#else
		const float* gaussian = (float*)mxGetData(prhs[ind]);
#endif
		const size_t size_gauss = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Create the output cell array
		mxArray* cell_array_ptr = mxCreateCellMatrix(size_reko + 1ULL, Nt);

		// Output dimensions
		const mwSize dim[4] = { static_cast<mwSize>(inputScalars.Nx), static_cast<mwSize>(inputScalars.Ny), static_cast<mwSize>(inputScalars.Nz), static_cast<mwSize>(outSize2) };


		inputScalars.PET = getScalarBool(getField(options, 0, "PET"), ind);
		inputScalars.CT = getScalarBool(getField(options, 0, "CT"), ind);
		inputScalars.SPECT = getScalarBool(getField(options, 0, "SPECT"), ind);
		inputScalars.PITCH = getScalarBool(getField(options, 0, "PITCH"), ind);
		inputScalars.listmode = getScalarBool(getField(options, 0, "listmode"), ind);
		if (inputScalars.projector_type == 5) {
			inputScalars.meanFP = getScalarBool(getField(options, 0, "meanFP"), ind);
			inputScalars.meanBP = getScalarBool(getField(options, 0, "meanBP"), ind);
		}
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			inputScalars.maskBP = getScalarBool(getField(options, 0, "useMaskBP"), ind);
		}
		inputScalars.maskFP = getScalarBool(getField(options, 0, "useMaskFP"), ind);
		if (inputScalars.projector_type == 2 || inputScalars.projector_type == 3) {
			inputScalars.orthXY = getScalarBool(getField(options, 0, "orthTransaxial"), ind);
			inputScalars.orthZ = getScalarBool(getField(options, 0, "orthAxial"), ind);
		}
		inputScalars.nProjections = getScalarInt64(getField(options, 0, "nProjections"), ind);
		inputScalars.subsetType = getScalarUInt32(getField(options, 0, "subset_type"), ind);
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			inputScalars.dL = getScalarFloat(getField(options, 0, "dL"), ind);
			inputScalars.d_Scale.s[0] = getScalarFloat(getField(options, 0, "dScaleX"), ind);
			inputScalars.d_Scale.s[1] = getScalarFloat(getField(options, 0, "dScaleY"), ind);
			inputScalars.d_Scale.s[2] = getScalarFloat(getField(options, 0, "dScaleZ"), ind);
			//inputScalars.d_Scale.s[3] = 0.f;
			inputScalars.dSize.s[0] = getScalarFloat(getField(options, 0, "dSizeX"), ind);
			inputScalars.dSize.s[1] = getScalarFloat(getField(options, 0, "dSizeY"), ind);
			inputScalars.dSize.s[2] = getScalarFloat(getField(options, 0, "dSizeZ"), ind);
			if (inputScalars.projector_type == 5) {
				inputScalars.dSizeBP.s[0] = getScalarFloat(getField(options, 0, "dSizeXBP"), ind);
				inputScalars.dSizeBP.s[1] = getScalarFloat(getField(options, 0, "dSizeZBP"), ind);
			}
		}

		size_t mDim = mxGetNumberOfElements(mxGetCell(Sin, 0));
		size_t koko = 0ULL;
		if (inputScalars.raw)
			koko = numRows / 2;
		else {
			koko = mDim / inputScalars.nBins;
			//for (int uu = 0; uu < nPituus; uu++)
			//	if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			//		koko += pituus[uu] * mDim;
			//	else
			//		koko += pituus[uu];
		}
		if (DEBUG) {
			mexPrintf("koko = %u\n", koko);
			mexPrintf("size_z = %u\n", size_z);
			mexPrintf("inputScalars.maskBP = %u\n", inputScalars.maskBP);
			mexPrintf("inputScalars.maskFP = %u\n", inputScalars.maskFP);
			mexPrintf("inputScalars.projector_type = %u\n", inputScalars.projector_type);
			mexEvalString("pause(.0001);");
		}
		//const uint8_t listmode = (uint8_t)mxGetScalar(getField(options, 0, "listmode"));
		//const bool CT = (bool)mxGetScalar(getField(options, 0, "CT"));

		//for (int tt = 0; tt < 6; tt++) {
		//	const char* varChar = mxArrayToString(mxGetCell(mxGetField(options, 0, "varList"), tt));
		//	if (DEBUG) {
		//		mexPrintf("%s\n", varChar);
		//		mexEvalString("pause(.0001);");
		//	}
		//}
		try {

			reconstruction_AF_matrixfree(koko, lor1, z_det, x, y, Sin, sc_ra, inputScalars, Niter, options, pituus, xy_index, 
				z_index, TotSinos, cell_array_ptr, dim, verbose, atten, size_atten, norm, size_norm, epps, k_path, Nt, pseudos, pRows, L, size_z, osem_bool,
				fileName, x_center, y_center, z_center, size_center_x, size_center_y, size_of_x, size_center_z,
				header_directory, device, use_64bit_atomics, n_rekos, n_rekos_mlem, reko_type, reko_type_mlem,
				V, size_V, gaussian, size_gauss, saveIter, TOFSize, TOFCenter);
			//}


			plhs[0] = cell_array_ptr;

			// Clear ArrayFire memory
			af::deviceGC();
		}
		catch (const std::exception& e) {
			af::deviceGC();
			mexErrMsgTxt(e.what());
		}
	}
	// Compute the number of voxels each LOR traverses (using AF device)
//	else if (type == 2) {
//
//		if (nrhs != 40)
//			mexErrMsgTxt("Invalid number of input arguments. There must be 40.");
//
//		// Starting block
//		const uint32_t block1 = getScalarUInt32(prhs[ind], ind);
//		ind++;
//
//		// Ending block
//		const uint32_t blocks = getScalarUInt32(prhs[ind], ind);
//		ind++;
//
//		// Number of sinograms
//		const uint32_t NSinos2 = getScalarUInt32(prhs[ind], ind);
//		ind++;
//
//		// Total number of sinograms
//		const uint16_t TotSinos2 = getScalarUInt16(prhs[ind], ind);
//		ind++;
//
//		size_t loop_var_par = 1ULL;
//
//		if (raw)
//			loop_var_par = numRows / 2ULL;
//		else
//			loop_var_par = static_cast<size_t>(NSinos2) * static_cast<size_t>(size_x);
//
//		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);
//
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		uint16_t* lor = (uint16_t*)mxGetUint16s(plhs[0]);
//#else
//		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);
//#endif
//
//		try {
//			find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos2, verbose, loop_var_par,
//				k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, size_of_x, header_directory);
//			af::deviceGC();
//		}
//		catch (const std::exception& e) {
//			af::deviceGC();
//			mexErrMsgTxt(e.what());
//		}
//
//	}

	return;
}