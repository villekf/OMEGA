/*************************************************************************************************************************************************
* Matrix free computations for OMEGA (Implementation 2).
* This is the main file for the matrix-free computations in OMEGA. In this file the MATLAB variables are loaded and either the matrix-free 
* reconstructions are computed or the custom prior reconstruction (custom priors are not supported anymore).
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus can be slightly inaccurate.
* This uses ArrayFire functions and thus requires the installation of the ArrayFire library.
* Note that in the current version, this same file is used for both OpenCL and CUDA backends despite the name suggesting otherwise.
* 
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*************************************************************************************************************************************************/
#include "reconstructionAF.h"

/// <summary>
/// The input mexFunction that uses the MATLAB/Octave input and gives the output produced by the implementation 2
/// </summary>
/// <param name="nlhs number of output parameters (lhs = left hand side)"></param>
/// <param name="plhs mxArray containing all the output variables (the image estimates and (optionally) the forward projections"></param>
/// <param name="nrhs number of input parameters (rhs = right hand side)"></param>
/// <param name="prhs mxArray containing all the input variables"></param>
/// <returns></returns>
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 51)
		mexErrMsgTxt("Too few input arguments. There must be exactly 51.");
	else if (nrhs > 51)
		mexErrMsgTxt("Too many input arguments. There must be exactly 51.");
	if (nlhs < 2)
		mexErrMsgTxt("Invalid number of output arguments. There must be at least two.");
	else if (nlhs > 3)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most three.");

	int ind = 0;
	scalarStruct inputScalars;
	// Load the input arguments
	// These are mainly same as with implementation 3 (check the comments in OpenCL_matrixfree_multi_gpu.cpp)

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
		mexPrintBase("size_z = %u\n", inputScalars.size_z);
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
	const uint32_t device = getScalarUInt32(prhs[ind], ind);
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
	ind++;

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
	const mwSize dim[5] = { static_cast<mwSize>(inputScalars.Nx[0]), static_cast<mwSize>(inputScalars.Ny[0]), static_cast<mwSize>(inputScalars.Nz[0]), static_cast<mwSize>(outSize2), static_cast<mwSize>(inputScalars.Nt) }; 

	loadInput(inputScalars, options);

	inputScalars.subsetsUsed = inputScalars.subsets;

	if (inputScalars.verbose >= 3) {
		mexPrint("Loaded MEX inputs");
	}

	size_t mDim = mxGetNumberOfElements(Sin) / static_cast<size_t>(inputScalars.Nt);

	mxArray* cell_array_ptr, * FPptr, *resPtr;
	if (CELL)
		cell_array_ptr = mxCreateCellMatrix(static_cast<mwSize>(inputScalars.nMultiVolumes) + 1, 1);
	else
		cell_array_ptr = mxCreateNumericArray(5, dim, mxSINGLE_CLASS, mxREAL);
	if (inputScalars.raw)
		inputScalars.koko = numRows / 2;
	else {
		inputScalars.koko = mDim / inputScalars.nBins;
	}
	if (inputScalars.storeFP) {
		FPptr = mxCreateCellMatrix(static_cast<mwSize>(inputScalars.subsets * inputScalars.Niter), 1);
	}
	else
		FPptr = mxCreateCellMatrix(1,1);

	if (inputScalars.storeResidual) {
		resPtr = mxCreateNumericMatrix(static_cast<mwSize>(inputScalars.subsets * inputScalars.Niter), 1, mxSINGLE_CLASS, mxREAL);
	}
	else
		resPtr = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	af::setDevice(device);

	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;

	// Struct containing the necessary variables for the priors
	Weighting w_vec;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);

	if (inputScalars.listmode) {
		if (inputScalars.indexBased) {
			w_vec.trIndex = getUint16s(options, "trIndex", 0);
			w_vec.axIndex = getUint16s(options, "axIndex", 0);
		}
		else {
			w_vec.listCoord = getSingles(options, "x", 0);
			//inputScalars.size_of_x = mxGetNumberOfElements(getField(options, 0, "x"));
		}
		if (inputScalars.TOF)
			w_vec.TOFIndices = getUint8s(options, "TOFIndices", 0);
	}

	if (DEBUG) {
		mexPrintBase("ind = %u\n", ind);
		mexPrintBase("koko = %u\n", inputScalars.koko);
		mexPrintBase("size_z = %u\n", inputScalars.size_z);
		mexPrintBase("inputScalars.maskBP = %u\n", inputScalars.maskBP);
		mexPrintBase("inputScalars.maskFP = %u\n", inputScalars.maskFP);
		mexPrintBase("inputScalars.offset = %u\n", inputScalars.offset);
		mexPrintBase("inputScalars.projector_type = %u\n", inputScalars.projector_type);
		mexPrintBase("inputScalars.FPType = %u\n", inputScalars.FPType);
		mexPrintBase("inputScalars.BPType = %u\n", inputScalars.BPType);
		mexPrintBase("inputScalars.useExtendedFOV = %u\n", inputScalars.useExtendedFOV);
		mexPrintBase("inputScalars.eFOV = %u\n", inputScalars.eFOV);
		mexPrintBase("inputScalars.TGV2D = %u\n", inputScalars.TGV2D);
		mexPrintBase("inputScalars.NxPrior = %u\n", inputScalars.NxPrior);
		mexPrintBase("inputScalars.NyPrior = %u\n", inputScalars.NyPrior);
		mexPrintBase("inputScalars.NzPrior = %u\n", inputScalars.NzPrior);
		mexPrintBase("inputScalars.im_dim = %u\n", inputScalars.im_dim[0]);
		mexPrintBase("inputScalars.Nx = %u\n", inputScalars.Nx[0]);
		mexPrintBase("inputScalars.Ny = %u\n", inputScalars.Ny[0]);
		mexPrintBase("inputScalars.Nz = %u\n", inputScalars.Nz[0]);
		mexPrintBase("inputScalars.Nf = %u\n", inputScalars.Nf);
		mexPrintBase("inputScalars.nColsD = %u\n", inputScalars.nColsD);
		mexPrintBase("inputScalars.nRowsD = %u\n", inputScalars.nRowsD);
		mexPrintBase("inputScalars.bmin = %f\n", inputScalars.bmin);
		mexPrintBase("inputScalars.bmax = %f\n", inputScalars.bmax);
		mexPrintBase("inputScalars.Vmax = %f\n", inputScalars.Vmax);
		mexPrintBase("inputScalars.size_V = %u\n", inputScalars.size_V);
		mexEval();
	}

	inputScalars.Nxy = inputScalars.Nx[0] * inputScalars.Ny[0];
	inputScalars.im_dim[0] = static_cast<int64_t>(inputScalars.Nxy) * static_cast<int64_t>(inputScalars.Nz[0]);
	if (inputScalars.multiResolution) {
		for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++)
			inputScalars.im_dim[ii] = static_cast<int64_t>(inputScalars.Nx[ii]) * static_cast<int64_t>(inputScalars.Ny[ii]) * static_cast<int64_t>(inputScalars.Nz[ii]);
	}

	// Load the necessary data from the MATLAB input (options) and create the necessary variables
	form_data_variables(w_vec, options, inputScalars, MethodList);

#if !defined MTYPE 
	const float* Sino = getSingles(Sin, "solu");
	const float* randoms = getSingles(sc_ra, "solu");
#else
	const uint16_t* Sino = getUint16s(Sin, "solu");
	const uint16_t* randoms = getUint16s(sc_ra, "solu");
#endif
	const float* extraCorr = getSingles(options, "ScatterC", 0);
	const float* x0 = getSingles(options, "x0");
	float* residual = getSingles(resPtr, "solu");

	if (inputScalars.verbose >= 3) {
		mexPrint("Loaded struct values. Starting reconstruction itself...");
	}
	try {

		int status = reconstructionAF(z_det, x, Sino, randoms, inputScalars, device, pituus, w_vec, MethodList, header_directory, x0,
			cell_array_ptr, FPptr, atten, norm, extraCorr, size_gauss, xy_index, z_index, residual, L);


		plhs[0] = cell_array_ptr;
		plhs[1] = FPptr;
		if (nlhs > 2)
			plhs[2] = resPtr;

		// Clear ArrayFire memory
		af::deviceGC();
		if ((inputScalars.verbose >= 3 || DEBUG) && status == 0)
			mexPrint("Reconstruction completed successfully!");
		else if (status != 0)
			mexPrint("Reconstruction failed!");
	}
	catch (const std::exception& e) {
		af::deviceGC();
		mexErrMsgTxt(e.what());
	}
	return;
}