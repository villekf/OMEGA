/**************************************************************************
* Matrix free computations for OMEGA for implementations 3 and 5.
* This code is very similar to the other matrix-free code, but this one
* can also be run without installing ArrayFire.
* Currently this supports host inputs as well as MATLAB gpuArray inputs
* (CUDA only at the moment). Since Apple hardware uses shared memory
* the host data is essentially device data on Metal side.
* 
* Copyright(C) 2020-2026 Ville-Veikko Wettenhovi
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
#ifdef MATLAB
#include "mfunctions.h"
#endif
#if defined(MATLABGPU)
#include "gpu/mxGPUArray.h"
#endif
#ifndef METAL
#include "ProjectorClass.h"
#endif
#include "multi_gpu_reconstruction.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
#if defined(MATLABGPU)
	// Initialize the MATLAB GPU API
	// This also makes MATLAB's CUDA context current, which the
	// projector then attaches to
	if (mxInitGPU() != MX_GPU_SUCCESS)
		mexErrMsgTxt("mxInitGPU failed to initialize the MATLAB GPU API.");
#endif
	// Check for the number of input and output arguments
	if (nrhs < 53)
		mexErrMsgTxt("Too few input arguments. There must be at least 53.");
	else if (nrhs > 53)
		mexErrMsgTxt("Too many input arguments. There can be at most 53.");

	if (nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There must be at least one.");
	else if (nlhs > 2)
		mexErrMsgTxt("Too many output arguments. There can be at most two.");

	int ind = 0;
	scalarStruct inputScalars;
	// Load the input arguments

	if (DEBUG) {
		mexPrintBase("ind0 = %u\n", ind);
		mexEval();
	}

	// The number of x-voxels in the estimated image
	size_t sX = mxGetNumberOfElements(prhs[ind]);
	uint32_t* Nx = getUint32s(prhs[ind], "solu");
	inputScalars.Nx = std::vector<uint32_t>(Nx, Nx + sX);
	ind++;

	// The number of y-voxels in the estimated image
	sX = mxGetNumberOfElements(prhs[ind]);
	uint32_t* Ny = getUint32s(prhs[ind], "solu");
	inputScalars.Ny = std::vector<uint32_t>(Ny, Ny + sX);
	ind++;

	// The number of z-voxels in the estimated image
	sX = mxGetNumberOfElements(prhs[ind]);
	uint32_t* Nz = getUint32s(prhs[ind], "solu");
	inputScalars.Nz = std::vector<uint32_t>(Nz, Nz + sX);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind1 = %u\n", ind);
		mexEval();
	}

	// The size of x-voxels in the estimated image
	sX = mxGetNumberOfElements(prhs[ind]);
	float* dx = getSingles(prhs[ind], "solu");
	inputScalars.dx = std::vector<float>(dx, dx + sX);
	ind++;

	// The size of y-voxels in the estimated image
	sX = mxGetNumberOfElements(prhs[ind]);
	float* dy = getSingles(prhs[ind], "solu");
	inputScalars.dy = std::vector<float>(dy, dy + sX);
	ind++;

	// The size of z-voxels in the estimated image
	sX = mxGetNumberOfElements(prhs[ind]);
	float* dz = getSingles(prhs[ind], "solu");
	inputScalars.dz = std::vector<float>(dz, dz + sX);
	ind++;

	// The distance from the origin to the corner of the image (x-direction)
	sX = mxGetNumberOfElements(prhs[ind]);
	float* bx = getSingles(prhs[ind], "solu");
	inputScalars.bx = std::vector<float>(bx, bx + sX);
	ind++;

	// The distance from the origin to the corner of the image (y-direction)
	sX = mxGetNumberOfElements(prhs[ind]);
	float* by = getSingles(prhs[ind], "solu");
	inputScalars.by = std::vector<float>(by, by + sX);
	ind++;

	// The distance from the origin to the corner of the image (z-direction)
	sX = mxGetNumberOfElements(prhs[ind]);
	float* bz = getSingles(prhs[ind], "solu");
	inputScalars.bz = std::vector<float>(bz, bz + sX);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind2 = %u\n", ind);
		mexEval();
	}

	// Coordinates of the detectors in z-direction (PET/SPECT) or the directional vectors for the detector panel pixels (CT)
	const float* z_det = getSingles(prhs[ind], "solu");
	inputScalars.size_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in x/y-directions
	float* x = getSingles(prhs[ind], "solu");
	inputScalars.size_of_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind3 = %u\n", ind);
		mexEval();
	}

	// The size of the first dimension in the input sinogram/projection
	inputScalars.nRowsD = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.verbose = getScalarInt8(prhs[ind], ind);
	ind++;

	// Detector pair numbers, for raw data
	const uint16_t* L = getUint16s(prhs[ind], "solu");
	const size_t numRows = mxGetM(prhs[ind]);
	inputScalars.sizeL = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind4 = %u\n", ind);
		mexEval();
	}

	// Is TOF data used?
	inputScalars.TOF = getScalarBool(prhs[ind], ind);
	ind++;

	// Size of single TOF-subset
	inputScalars.TOFSize = getScalarInt64(prhs[ind], ind);
	ind++;

	// Variance of the Gaussian TOF
	inputScalars.sigma_x = getScalarFloat(prhs[ind], ind);
	ind++;

	// Centers of the TOF-bins
	inputScalars.TOFCenter = getSingles(prhs[ind], "solu");
	ind++;

	// Index offset for TOF subsets
	inputScalars.nBins = getScalarInt64(prhs[ind], ind);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind5 = %u\n", ind);
		mexEval();
	}

	// The device used
	inputScalars.platform = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.raw = getScalarUInt8(prhs[ind], ind);
	ind++;

	inputScalars.use_psf = getScalarBool(prhs[ind], ind);
	ind++;

	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind6 = %u\n", ind);
		mexEval();
	}

	// attenuation values
	const float* atten = getSingles(prhs[ind], "solu");
	inputScalars.size_atten = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Normalization coefficients
	const float* norm = getSingles(prhs[ind], "solu");
	inputScalars.size_norm = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind7 = %u\n", ind);
		mexEval();
	}

	// Number of measurements/LORs
	const int64_t* pituus = getInt64s(prhs[ind], "solu");
	const size_t nPituus = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("nPituus = %u\n", nPituus);
		mexEval();
	}

	// Is the attenuation correction included
	inputScalars.attenuation_correction = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Is the normalization correction included
	inputScalars.normalization_correction = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.Niter = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.subsets = getScalarUInt32(prhs[ind], ind);
	ind++;

	inputScalars.epps = getScalarFloat(prhs[ind], ind);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind8 = %u\n", ind);
		mexEval();
	}

	// XY-indices of the detector coordinates of each LOR
	const uint32_t* xy_index = getUint32s(prhs[ind], "solu");
	inputScalars.sizeXY = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Z-indices of the detector coordinates of each LOR
	const uint16_t* z_index = getUint16s(prhs[ind], "solu");
	inputScalars.sizeZ = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind9 = %u\n", ind);
		mexEval();
	}

	inputScalars.tube_width = getScalarFloat(prhs[ind], ind);
	ind++;

	// Center coordinates of voxels in the X-dimension
	inputScalars.x_center = getSingles(prhs[ind], "solu");
	inputScalars.size_center_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Center coordinates of voxels in the Y-dimension
	inputScalars.y_center = getSingles(prhs[ind], "solu");
	inputScalars.size_center_y = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Center coordinates of voxels in the Z-dimension
	inputScalars.z_center = getSingles(prhs[ind], "solu");
	inputScalars.size_center_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind10 = %u\n", ind);
		mexEval();
	}

	// Randoms
	const mxArray* sc_ra = prhs[ind];
	ind++;

	// Randoms corrections
	inputScalars.randoms_correction = getScalarUInt32(prhs[ind], ind);
	ind++;

	// The type of projector used (Siddon or orthogonal)
	inputScalars.projector_type = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Number of rays in Siddon
	inputScalars.n_rays = getScalarUInt16(prhs[ind], ind);
	ind++;

	// Number of rays in Siddon (axial)
	inputScalars.n_rays3D = getScalarUInt16(prhs[ind], ind);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind11 = %u\n", ind);
		mexEval();
	}

	const mxArray* options = prhs[ind];
	ind++;

	//Cell array containing the measurements
	const mxArray* Sin = prhs[ind];
	ind++;

	// Number of time steps
	inputScalars.Nt = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Use 64-bit integer atomic functions if possible
	inputScalars.atomic_64bit = getScalarBool(prhs[ind], ind);
	ind++;

	if (DEBUG) {
		mexPrintBase("ind12 = %u\n", ind);
		mexEval();
	}

	inputScalars.bmin = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.bmax = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.Vmax = getScalarFloat(prhs[ind], ind);
	ind++;

	inputScalars.V = getSingles(prhs[ind], "solu");
	inputScalars.size_V = mxGetNumberOfElements(prhs[ind]);
	ind++;

	inputScalars.gaussian = getSingles(prhs[ind], "solu");
	const size_t size_gauss = mxGetNumberOfElements(prhs[ind]);
	ind++;

	const int type = getScalarInt32(prhs[ind], ind);
	ind++;

	const int no_norm = getScalarInt32(prhs[ind], ind);
	ind++;

	if (inputScalars.verbose >= 3) {
		mexPrint("Loaded MEX inputs");
	}

	inputScalars.saveIter = getScalarBool(options, 0, "save_iter");
	inputScalars.saveNIter = getUint32s(options, "saveNIter");
	inputScalars.saveIterationsMiddle = mxGetNumberOfElements(getField(options, 0, "saveNIter"));
	size_t Ni = 0ULL;
	if (inputScalars.saveIter)
		Ni = static_cast<size_t>(inputScalars.Niter);
	else if (inputScalars.saveIterationsMiddle > 0)
		Ni = inputScalars.saveIterationsMiddle;
	const size_t outSize = static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * static_cast<size_t>(inputScalars.Nz[0]);
	const size_t outSize2 = Ni + 1ULL;

	// Output dimensions
	const mwSize dim[4] = { static_cast<mwSize>(inputScalars.Nx[0]), static_cast<mwSize>(inputScalars.Ny[0]), static_cast<mwSize>(inputScalars.Nz[0]), static_cast<mwSize>(outSize2) };

	loadInput(inputScalars, options, type);
	inputScalars.subsetsUsed = getScalarUInt32(getField(options, 0, "subsets"));
	inputScalars.timestepsUsed = inputScalars.Nt;
	if (type > 0) {
		inputScalars.osa_iter0 = getScalarUInt32(getField(options, 0, "currentSubset"));
		inputScalars.subsetsUsed = inputScalars.osa_iter0 + 1;
		inputScalars.timestep0 = getScalarUInt32(getField(options, 0, "currentTimestep"));
		inputScalars.timestepsUsed = inputScalars.timestep0 + 1;
	}

	// Detect MATLAB gpuArray inputs (MATLABGPU/CUDA build only)
	// This depends on the "type", that is whether forward or backprojection is used
	// as well as whether host or gpuArray data is used
	// If gpuArray data is present, then the type uses gpuArrays for both inputs
	// and outputs, otherwise host data for both
	bool imOnDevice = false;
	bool measOnDevice = false;
#if defined(MATLABGPU)
	mxGPUArray const* gIm = nullptr;
	mxGPUArray const* gMeas = nullptr;
	mxGPUArray* outGPU = nullptr;
	mxGPUArray* sensGPU = nullptr;
	const mxArray* x0Field = getField(options, 0, "x0");
	{
		const bool x0IsGPU = mxIsGPUArray(x0Field);
		const bool SinIsGPU = mxIsGPUArray(Sin);
		if (type == 0 && (x0IsGPU || SinIsGPU))
			mexErrMsgTxt("Implementation 3 (type == 0) does not support gpuArray inputs.");
		imOnDevice = (type == 1) && x0IsGPU;
		measOnDevice = (type == 2) && SinIsGPU;
	}
#endif

	if (DEBUG) {
		mexPrint("Set output vector");
	}

	mwSize mDim[1] = { 1 };
	mwSize d[1] = { 1 };

	if (type == 1) {
		mDim[0] = pituus[inputScalars.osa_iter0 + 1] - pituus[inputScalars.osa_iter0];
		if (DEBUG) {
			mexPrintBase("mDim = %u\n", mDim[0]);
			mexPrintBase("inputScalars.osa_iter0 = %u\n", inputScalars.osa_iter0);
			mexEval();
		}
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			mDim[0] = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * mDim[0];
		mDim[0] *= inputScalars.nBins;
		if (DEBUG) {
			mexPrintBase("mDim = %u\n", mDim[0]);
			mexEval();
		}
	}
	else if (type == 2 || type == 0)
		mDim[0] = std::accumulate(inputScalars.im_dim.begin(), inputScalars.im_dim.end(), (int64_t)0);

	if (DEBUG) {
		mexPrintBase("type = %u\n", type);
		mexPrintBase("mDim = %u\n", mDim[0]);
		mexEval();
	}
	mxArray* array_ptr = nullptr;
	mxArray* sens_ptr = nullptr;
	// Previously the host outputs were created here in the conditional
	// However, now that the gpuArray support is included as well it is
	// codewise more efficient to define only the sizes and class types here
	// These can then be shared in the output creation, whether a host or
	// gpuArray output is used
	mwSize const* arrDims = mDim;
	mxClassID arrClassID = mxSINGLE_CLASS;
	mwSize const* sensDims = d;
	mxClassID sensClassID = mxSINGLE_CLASS;
	if (type == 1) {
		arrDims = mDim; 
		arrClassID = mxSINGLE_CLASS;
		sensDims = d; 
		sensClassID = mxSINGLE_CLASS;
	}
	else if (type == 2 && (inputScalars.atomic_32bit || inputScalars.atomic_64bit)) {
		if (no_norm == 0)
			if (inputScalars.atomic_32bit) {
				sensDims = mDim; 
				sensClassID = mxINT32_CLASS;
			}
			else {
				sensDims = mDim; 
				sensClassID = mxINT64_CLASS;
			}
		else
			if (inputScalars.atomic_32bit) {
				sensDims = d; 
				sensClassID = mxINT32_CLASS;
			}
			else {
				sensDims = d; 
				sensClassID = mxINT64_CLASS;
			}
		if (inputScalars.atomic_32bit) {
			arrDims = mDim; 
			arrClassID = mxINT32_CLASS;
		}
		else {
			arrDims = mDim; 
			arrClassID = mxINT64_CLASS;
		}
	}
	else {
		arrDims = mDim; arrClassID = mxSINGLE_CLASS;
		if (no_norm == 0 && type == 2) {
			sensDims = mDim; 
			sensClassID = mxSINGLE_CLASS;
		}
		else
			if (type == 0 && inputScalars.atomic_32bit) {
				sensDims = d; 
				sensClassID = mxINT32_CLASS;
			}
			else if (type == 0 && inputScalars.atomic_64bit) {
				sensDims = mDim; 
				sensClassID = mxINT64_CLASS;
			}
			else {
				sensDims = d; 
				sensClassID = mxSINGLE_CLASS;
			}
	}

#if defined(MATLABGPU)
	if ((type == 1 && imOnDevice) || (type == 2 && measOnDevice)) {
		// gpuArray outputs
		outGPU = mxGPUCreateGPUArray(1, arrDims, arrClassID, mxREAL, MX_GPU_INITIALIZE_VALUES);
		sensGPU = mxGPUCreateGPUArray(1, sensDims, sensClassID, mxREAL, MX_GPU_INITIALIZE_VALUES);
	}
	else
#endif
	{
		// Host outputs
		array_ptr = mxCreateNumericArray(1, arrDims, arrClassID, mxREAL);
		sens_ptr = mxCreateNumericArray(1, sensDims, sensClassID, mxREAL);
	}

	if (DEBUG) {
		mexPrint("Output vector set");
	}

	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;

	// Struct containing the necessary variables for the priors
	Weighting w_vec;

	if (type == 0) {
		// Obtain the reconstruction methods used
		get_rec_methods(options, MethodList);
	}
	// Load the necessary data from the MATLAB input (options) and create the necessary variables
	form_data_variables(w_vec, options, inputScalars, MethodList);

	if (DEBUG) {
		mexPrint("Reconstruction methods obtained");
	}
	// Add index-based reconstruction support
	if (inputScalars.listmode) {
		if (inputScalars.indexBased) {
			w_vec.trIndex = getUint16s(options, "trIndex", 0);
			w_vec.axIndex = getUint16s(options, "axIndex", 0);
		}
		else {
			w_vec.listCoord = getSingles(options, "x", 0);
		}
		if (inputScalars.TOF)
			w_vec.TOFIndices = getUint8s(options, "TOFIndices", 0);
	}

	const float* Sino = nullptr;
#if defined(MATLABGPU)
	if (measOnDevice) {
		gMeas = mxGPUCreateFromMxArray(Sin);
		if (mxGPUGetClassID(gMeas) != mxSINGLE_CLASS || mxGPUGetComplexity(gMeas) != mxREAL) {
			mxGPUDestroyGPUArray(gMeas);
			mexErrMsgTxt("The measurement gpuArray input must be real, single precision.");
		}
		Sino = (const float*)mxGPUGetDataReadOnly(gMeas);
		inputScalars.size_meas = mxGPUGetNumberOfElements(gMeas);
	}
	else
#endif
	{
		Sino = getSingles(Sin, "solu");
		inputScalars.size_meas = mxGetNumberOfElements(Sin);
	}
	const float* randoms = getSingles(sc_ra, "solu");
	const float* extraCorr = getSingles(options, "ScatterC", 0);
	const float* x0 = nullptr;
#if defined(MATLABGPU)
	if (imOnDevice) {
		gIm = mxGPUCreateFromMxArray(x0Field);
		if (mxGPUGetClassID(gIm) != mxSINGLE_CLASS || mxGPUGetComplexity(gIm) != mxREAL) {
			mxGPUDestroyGPUArray(gIm);
			mexErrMsgTxt("The initial image (options.x0) gpuArray input must be real, single precision.");
		}
		x0 = (const float*)mxGPUGetDataReadOnly(gIm);
	}
	else
#endif
	{
		x0 = getSingles(options, "x0");
	}
	if (DEBUG && !imOnDevice) {
		mexPrintBase("x0[0] = %f\n", x0[0]);
		mexPrintBase("x0.dim = %u\n", mxGetNumberOfElements(mxGetField(options, 0, "x0")));
		mexEval();
	}

	if (DEBUG) {
		mexPrint("Pointers set");
	}
#if defined(MATLABGPU)
	// Fill in the device-resident inputs/outputs
	// Members are left null (i.e. host arrays are used) unless this call's own input for this type is 
	// a gpuArray; the outputs (outGPU/sensGPU) were only allocated above under that same condition, so 
	// this stays consistent with the "output lives where the input lives"
	deviceIO devIO;
	if (imOnDevice)
		devIO.im = x0;
	if (measOnDevice)
		devIO.meas = Sino;
	if (outGPU)
		devIO.output = mxGPUGetData(outGPU);
	if (sensGPU)
		devIO.sensIm = mxGPUGetData(sensGPU);
#endif
	if (inputScalars.atomic_32bit && (type == 2)) {
		int32_t* output = array_ptr ? getInt32s(array_ptr, "solu") : nullptr;
		int32_t* sensIm = sens_ptr ? getInt32s(sens_ptr, "solu") : nullptr;
		reconstruction_multigpu(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L
#if defined(MATLABGPU)
			, devIO
#endif
		);
#if defined(MATLABGPU)
		if (outGPU) {
			plhs[0] = mxGPUCreateMxArrayOnGPU(outGPU);
			mxGPUDestroyGPUArray(outGPU);
			if (nlhs > 1)
				plhs[1] = mxGPUCreateMxArrayOnGPU(sensGPU);
			mxGPUDestroyGPUArray(sensGPU);
		}
		else
#endif
		{
			plhs[0] = array_ptr;
			if (nlhs > 1)
				plhs[1] = sens_ptr;
			else
				mxDestroyArray(sens_ptr);
		}
	}
	else if (inputScalars.atomic_64bit && (type == 2)) {
		int64_t* output = array_ptr ? getInt64s(array_ptr, "solu") : nullptr;
		int64_t* sensIm = sens_ptr ? getInt64s(sens_ptr, "solu") : nullptr;
		reconstruction_multigpu(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L
#if defined(MATLABGPU)
			, devIO
#endif
		);
#if defined(MATLABGPU)
		if (outGPU) {
			plhs[0] = mxGPUCreateMxArrayOnGPU(outGPU);
			mxGPUDestroyGPUArray(outGPU);
			if (nlhs > 1)
				plhs[1] = mxGPUCreateMxArrayOnGPU(sensGPU);
			mxGPUDestroyGPUArray(sensGPU);
		}
		else
#endif
		{
			plhs[0] = array_ptr;
			if (nlhs > 1)
				plhs[1] = sens_ptr;
			else
				mxDestroyArray(sens_ptr);
		}
	}
	else if (inputScalars.atomic_64bit && (type == 0)) {
		float* output = array_ptr ? getSingles(array_ptr, "solu") : nullptr;
		int64_t* sensIm = sens_ptr ? getInt64s(sens_ptr, "solu") : nullptr;
		reconstruction_multigpu(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L
#if defined(MATLABGPU)
			, devIO
#endif
		);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else if (inputScalars.atomic_32bit && (type == 0)) {
		float* output = array_ptr ? getSingles(array_ptr, "solu") : nullptr;
		int32_t* sensIm = sens_ptr ? getInt32s(sens_ptr, "solu") : nullptr;
		reconstruction_multigpu(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L
#if defined(MATLABGPU)
			, devIO
#endif
		);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else {
		float* output = array_ptr ? getSingles(array_ptr, "solu") : nullptr;
		float* sensIm = sens_ptr ? getSingles(sens_ptr, "solu") : nullptr;
		reconstruction_multigpu(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L
#if defined(MATLABGPU)
			, devIO
#endif
		);
#if defined(MATLABGPU)
		if (outGPU) {
			plhs[0] = mxGPUCreateMxArrayOnGPU(outGPU);
			mxGPUDestroyGPUArray(outGPU);
			if (nlhs > 1)
				plhs[1] = mxGPUCreateMxArrayOnGPU(sensGPU);
			mxGPUDestroyGPUArray(sensGPU);
		}
		else
#endif
		{
			plhs[0] = array_ptr;
			if (nlhs > 1)
				plhs[1] = sens_ptr;
			else
				mxDestroyArray(sens_ptr);
		}
	}

#if defined(MATLABGPU)
	if (gIm)
		mxGPUDestroyGPUArray(gIm);
	if (gMeas)
		mxGPUDestroyGPUArray(gMeas);
#endif

	return;
}