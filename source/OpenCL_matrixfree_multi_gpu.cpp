/**************************************************************************
* Matrix free computations for OMEGA for the multi-GPU/device case.
* Supports heterogeneous computing, but only on the same platform (i.e. if 
* you have a CPU and an integrated GPU from the same vendor and OpenCL
* runtime for both you can utilize both of them at the same time). Mixing
* CPU and GPU requires the optimization of the GPU to CPU value.
* This code can also be run for single device.
* This code is very similar to the other matrix-free code, but this one
* can also be run without installing ArrayFire, i.e. this uses pure OpenCL.
* For purely OSEM or MLEM reconstructions this code should be the fastest,
* regardless of the amount of devices used.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
* 
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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
#include "functions_multigpu.hpp"



void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 38)
		mexErrMsgTxt("Too few input arguments.  There must be at least 38.");
	else if (nrhs > 72)
		mexErrMsgTxt("Too many input arguments.  There can be at most 72.");

	if (nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be at least one.");
	else if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.  There can be at most one.");

	scalarStruct inputScalars;

	// Load the input arguments
	int ind = 0;

	// Path to the kernel (.cl) files
	const char* k_path = mxArrayToString(prhs[ind]);
	ind++;

	// Image size in y-direction
	inputScalars.Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in x-direction
	inputScalars.Nx = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in z-direction
	inputScalars.Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in x-direction
	inputScalars.dx = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in z-direction
	inputScalars.dz = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	inputScalars.by = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	inputScalars.bx = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	inputScalars.bz = getScalarFloat(prhs[ind], ind);
	ind++;

	// Coordinates of the detectors in z-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* z_det = (float*)mxGetSingles(prhs[ind]);
#else
	const float* z_det = (float*)mxGetData(prhs[ind]);
#endif
	inputScalars.size_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in x-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* x = (float*)mxGetSingles(prhs[ind]);
#else
	const float* x = (float*)mxGetData(prhs[ind]);
#endif
	const size_t numel_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in y-direction
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* y = (float*)mxGetSingles(prhs[ind]);
#else
	const float* y = (float*)mxGetData(prhs[ind]);
#endif
	inputScalars.numelY = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		//mexPrintf("y[363] = %f\n", y[363]);
		//mexPrintf("y[362] = %f\n", y[362]);
		//mexPrintf("y[364] = %f\n", y[364]);
		mexPrintf("numel_y = %u\n", inputScalars.numelY);
		mexEvalString("pause(.0001);");
	}

	// Distance between adjacent pixels in y-direction
	inputScalars.dy = getScalarFloat(prhs[ind], ind);
	ind++;

	// The maximum elements of the pixel space in both x- and y-directions
	inputScalars.maxyy = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.maxxx = getScalarFloat(prhs[ind], ind);
	ind++;

	// Number of slices included
	inputScalars.NSlices = getScalarFloat(prhs[ind], ind);
	ind++;

	// Number of detector indices
	inputScalars.size_x = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Maximum value of the z-direction detector coordinates
	inputScalars.zmax = getScalarFloat(prhs[ind], ind);
	ind++;

	// Are status messages displayed
	const bool verbose = getScalarBool(prhs[ind], ind);
	ind++;

	// Detector pair numbers, for raw data
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
	inputScalars.pRows = (uint32_t)mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Number of detectors per ring
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

	// Platform used
	const uint32_t device = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Filename of the current kernel (.cl) file
	const char* fileName = mxArrayToString(prhs[ind]);
	ind++;

	// Is raw data used
	inputScalars.raw = getScalarUInt8(prhs[ind], ind);
	ind++;

	// Coefficient that determines how much more measurements the other device (usually GPU) has compared to the other device
	const float kerroin = getScalarFloat(prhs[ind], ind);
	ind++;

	// Which of the below sections are used (i.e. is this a precomputation phase, implementation 3 or forward-backwards projection)
	const uint32_t type = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[ind]);
	ind++;
	if (DEBUG) {
		mexPrintf("ind = %u\n", ind);
		mexEvalString("pause(.0001);");
	}

	const mxArray* options = prhs[ind];
	ind++;

	inputScalars.bzb = inputScalars.bz + static_cast<float>(inputScalars.Nz) * inputScalars.dz;

	std::string s(fileName);

	s += (std::to_string(device) + ".bin");

	fileName = s.c_str();

	// Fixed local size
	size_t local_size[2] = { 64ULL, 0ULL };
	//local_size[0] = 64ULL;
	//local_size[1] = 0ULL;

	inputScalars.PET = getScalarBool(getField(options, 0, "PET"), ind);
	if (DEBUG) {
		mexPrintf("inputScalars.PET = %u\n", inputScalars.PET);
		mexEvalString("pause(.0001);");
	}
	inputScalars.SPECT = getScalarBool(getField(options, 0, "SPECT"), ind);
	if (DEBUG) {
		mexPrintf("inputScalars.SPECT = %u\n", inputScalars.SPECT);
		mexEvalString("pause(.0001);");
	}
	inputScalars.PITCH = getScalarBool(getField(options, 0, "PITCH"), ind);
	inputScalars.nProjections = getScalarInt64(getField(options, 0, "nProjections"), ind);
	inputScalars.subsetType = getScalarUInt32(getField(options, 0, "subset_type"), ind);
	if ((bool)mxGetScalar(getField(options, 0, "useSubsets")))
		inputScalars.subsets = 2U;
	const uint8_t listmode = (uint8_t)mxGetScalar(getField(options, 0, "listmode"));
	const bool CT = (bool)mxGetScalar(getField(options, 0, "CT"));
	if (DEBUG) {
		mexPrintf("inputScalars.PITCH = %u\n", inputScalars.PITCH);
		mexPrintf("CT = %u\n", CT);
		mexPrintf("inputScalars.subsets = %u\n", inputScalars.subsets);
		mexPrintf("inputScalars.size_z = %u\n", inputScalars.size_z);
		mexPrintf("inputScalars.attenuation_correction = %u\n", inputScalars.attenuation_correction);
		mexPrintf("inputScalars.normalization_correction = %u\n", inputScalars.normalization_correction);
		mexPrintf("inputScalars.subsetType = %u\n", inputScalars.subsetType);
		mexEvalString("pause(.0001);");
	}


	if (type == 2) {

		if (nrhs != 38)
			mexErrMsgTxt("Incorrect number of input arguments. There has to be 38.");

		// Starting ring
		const uint32_t block1 = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Ending ring
		const uint32_t blocks = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Number of sinograms used
		const uint32_t NSinos = getScalarUInt32(prhs[ind], ind);
		ind++;

		// The total number of sinograms
		const uint16_t TotSinos = getScalarUInt32(prhs[ind], ind);
		ind++;

		size_t loop_var_par = 1ULL;

		if (inputScalars.raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(inputScalars.size_x);

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		uint16_t* lor = (uint16_t*)mxGetUint16s(plhs[0]);
#else
		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);
#endif

		// Find the number of voxels each LOR traverses
		find_LORs(lor, z_det, x, y, TotSinos, verbose, loop_var_par, inputScalars, k_path, pseudos, L, fileName, device, numel_x, header_directory, local_size);

	}
	else if (type < 2) {

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
		ind++;

		// Is the attenuation correction included
		inputScalars.attenuation_correction = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Is the normalization correction included
		inputScalars.normalization_correction = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Number of voxels that each LOR/ray traverses
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
		const size_t koko_l = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Z-indices of the detector coordinates of each LOR
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		const uint16_t* z_index = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// For 2D orthogonal, the "strip width"
		//inputScalars.tube_width = getScalarFloat(prhs[ind], ind);
		//ind++;

		// For orthogonal, either the "tube width" or strip width
		inputScalars.tube_width = getScalarFloat(prhs[ind], ind);
		ind++;

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

		// Number of rays in Siddon (transaxial)
		inputScalars.n_rays = getScalarUInt16(prhs[ind], ind);
		ind++;

		// Number of rays in Siddon (axial)
		inputScalars.n_rays3D = getScalarUInt16(prhs[ind], ind);
		ind++;

		// Crystal pitch in z-direction
		inputScalars.cr_pz = getScalarFloat(prhs[ind], ind);
		ind++;

		// Measurement data
		const mxArray* Sin = prhs[ind];
		ind++;

		// Use 64-bit integer atomic functions if possible
		const bool use_64bit_atomics = getScalarBool(prhs[ind], ind);
		ind++;

		size_t koko;
		if (inputScalars.raw)
			koko = numRows;
		else
			koko = koko_l;

		if (inputScalars.projector_type > 0 && inputScalars.projector_type < 4)
			local_size[0] = 128ULL;
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || ((inputScalars.PET || inputScalars.SPECT || CT) && listmode == 0)) {
			local_size[0] = 8ULL;
			local_size[1] = 8ULL;
		}
		if (inputScalars.projector_type == 5) {
			inputScalars.meanFP = getScalarBool(getField(options, 0, "meanFP"), ind);
			inputScalars.meanBP = getScalarBool(getField(options, 0, "meanBP"), ind);
		}
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			inputScalars.maskFP = getScalarBool(getField(options, 0, "useMaskFP"), ind);
			inputScalars.maskBP = getScalarBool(getField(options, 0, "useMaskBP"), ind);
		}
		if (inputScalars.projector_type == 2 || inputScalars.projector_type == 3) {
			inputScalars.orthXY = getScalarBool(getField(options, 0, "orthTransaxial"), ind);
			inputScalars.orthZ = getScalarBool(getField(options, 0, "orthAxial"), ind);
		}
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			//ind += 6;
//				//if (DEBUG) {
//				//	mexPrintf("ind = %u\n", ind);
//				//	mexEvalString("pause(.0001);");
//				//}
			inputScalars.dL = getScalarFloat(getField(options, 0, "dL"), ind);
			//inputScalars.detY = getScalarFloat(mxGetField(options, 0, "detY"), ind);
			//inputScalars.tStart = getScalarFloat(mxGetField(options, 0, "tStart"), ind);
			//inputScalars.tStep = getScalarFloat(mxGetField(options, 0, "tStep"), ind);
			//inputScalars.d_Scale = { { getScalarFloat(mxGetField(options, 0, "dScaleX"), ind), getScalarFloat(mxGetField(options, 0, "dScaleY"), ind),
			//	getScalarFloat(mxGetField(options, 0, "dScaleZ"), ind), 0.f } };
			inputScalars.d_Scale.s[0] = getScalarFloat(getField(options, 0, "dScaleX"), ind);
			inputScalars.d_Scale.s[1] = getScalarFloat(getField(options, 0, "dScaleY"), ind);
			inputScalars.d_Scale.s[2] = getScalarFloat(getField(options, 0, "dScaleZ"), ind);
			//inputScalars.d_Scale.s[3] = 0.f;
			inputScalars.dSize.s[0] = getScalarFloat(getField(options, 0, "dSizeX"), ind);
			inputScalars.dSize.s[1] = getScalarFloat(getField(options, 0, "dSizeY"), ind);
			inputScalars.dSize.s[2] = getScalarFloat(getField(options, 0, "dSizeZ"), ind);
			//if (inputScalars.SPECT) {
			//	inputScalars.cThickness = getScalarFloat(mxGetField(options, 0, "cThickness"), ind);
			//	inputScalars.cSizeX = getScalarUInt32(mxGetField(options, 0, "cSizeX"), ind);
			//	inputScalars.cSizeY = getScalarUInt32(mxGetField(options, 0, "cSizeY"), ind);
			//}
			if (inputScalars.projector_type == 5) {
				inputScalars.dSizeBP.s[0] = getScalarFloat(getField(options, 0, "dSizeXBP"), ind);
				inputScalars.dSizeBP.s[1] = getScalarFloat(getField(options, 0, "dSizeZBP"), ind);
				//inputScalars.meanV = getScalarFloat(mxGetField(options, 0, "meanV"), ind);
			}
			//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			//				offsetH = (float*)mxGetSingles(mxGetField(options, 0, "horizontalOffset"));
			//#else
			//				offsetH = (float*)mxGetData(mxGetField(options, 0, "horizontalOffset"));
			//#endif
			//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			//				offsetB = (float*)mxGetSingles(mxGetField(options, 0, "bedOffset"));
			//#else
			//				offsetB = (float*)mxGetData(mxGetField(options, 0, "bedOffset"));
			//#endif
		}

		if (type == 0) {

			if (nrhs != 64)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 64.");

			// Right hand side for forward or backprojections
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			const float* rhs = (float*)mxGetSingles(prhs[ind]);
#else
			const float* rhs = (float*)mxGetData(prhs[ind]);
#endif
			const size_t size_rhs = mxGetNumberOfElements(prhs[ind]);
			ind++;

			// Is the normalization constant computed (sum(A))
			const cl_uchar no_norm = getScalarUInt8(prhs[ind], ind);
			ind++;

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

			// Use PSF in Siddon
			inputScalars.use_psf = getScalarBool(prhs[ind], ind);
			ind++;

			//float SD = 0.f, SO = 0.f;
			//float* offsetH = nullptr, *offsetB = nullptr;
			if (DEBUG) {
				mexPrintf("ind = %u\n", ind);
				mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
				mexPrintf("local_size[0] = %u\n", local_size[0]);
				mexEvalString("pause(.0001);");
			}
			if (DEBUG) {
				mexPrintf("Variables loaded\n");
				mexEvalString("pause(.0001);");
			}

			const size_t outSize3 = static_cast<size_t>(inputScalars.Nx) * static_cast<size_t>(inputScalars.Ny) * static_cast<size_t>(inputScalars.Nz);
			size_t outSize;
			if (size_rhs == outSize3)
				outSize = pituus[0];
			else
				outSize = outSize3;
			const size_t outSize2 = 1;

			const uint16_t TotSinos = inputScalars.size_z / 2ULL;

			mxArray* output;
			output = mxCreateCellMatrix(2, 1);
			if (DEBUG) {
				mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
				mexEvalString("pause(.0001);");
			}

			// Forward/backward projection
			reconstruction_f_b_proj(koko, lor1, z_det, x, y, rhs, sc_ra, inputScalars, pituus, koko_l, xy_index, z_index, TotSinos, verbose, 
				atten, size_atten, norm, size_norm, k_path, pseudos, L, fileName, device, kerroin, output, size_rhs, no_norm, numel_x, x_center, 
				y_center, z_center, size_center_x, size_center_y, size_center_z, header_directory,	Sin, use_64bit_atomics, V, size_V, local_size, 
				options, TOFSize, TOFCenter);
			plhs[0] = output;

		}
		else if (type == 1) {

			if (nrhs != 72)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 72.");

			// Number of sinograms used
			const uint32_t NSinos = getScalarUInt32(prhs[ind], ind);
			ind++;

			// Total number of sinograms
			const uint16_t TotSinos = getScalarUInt16(prhs[ind], ind);
			ind++;

			// Number of iterations
			inputScalars.Niter = getScalarUInt32(prhs[ind], ind);
			ind++;

			// Number of subsets
			inputScalars.subsets = getScalarUInt32(prhs[ind], ind);
			ind++;

			// Which reconstruction methods are used
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			const uint8_t* rekot = (uint8_t*)mxGetSingles(prhs[ind]);
#else
			const uint8_t* rekot = (uint8_t*)mxGetData(prhs[ind]);
#endif
			const size_t size_reko = mxGetNumberOfElements(prhs[ind]);
			ind++;

			// Epsilon value
			inputScalars.epps = getScalarFloat(prhs[ind], ind);
			ind++;

			// Number of time steps
			inputScalars.Nt = getScalarUInt32(prhs[ind], ind);
			ind++;

			// Is OSEM used
			const bool osem_bool = getScalarBool(prhs[ind], ind);
			ind++;

			// Use PSF in Siddon
			inputScalars.use_psf = getScalarBool(prhs[ind], ind);
			ind++;

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

			const size_t outSize = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;
			size_t Ni = 0U;
			const bool saveIter = getScalarBool(getField(prhs[ind], 0, "save_iter"), -1);
			if (saveIter)
				Ni = static_cast<size_t>(inputScalars.Niter);
			const size_t outSize2 = Ni + 1ULL;

			// Create the output cell array
			mxArray* cell_array_ptr;
			cell_array_ptr = mxCreateCellMatrix(1, inputScalars.Nt);

			// Implementation 3 (multi-device OpenCL)
			reconstruction_multigpu(koko, lor1, z_det, x, y, Sin, sc_ra, prhs[ind], inputScalars, pituus, koko_l, xy_index, z_index, TotSinos, 
				cell_array_ptr, verbose, atten, size_atten, norm, size_norm, rekot, k_path, size_reko, pseudos, L, osem_bool, fileName, device, 
				kerroin, numel_x, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, header_directory, use_64bit_atomics, 
				V, size_V, local_size, gaussian, size_gauss, TOFSize, TOFCenter);

			plhs[0] = cell_array_ptr;
		}
	}

	return;
}