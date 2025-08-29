#ifdef MATLAB
#include "mfunctions.h"
#endif
#include "reconstruction_metal.mm"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
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
	if (type > 0) {
		inputScalars.osa_iter0 = getScalarUInt32(getField(options, 0, "currentSubset"));
		inputScalars.subsetsUsed = inputScalars.osa_iter0 + 1;
	}


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
	if (type == 1) {
		array_ptr = mxCreateNumericArray(1, mDim, mxSINGLE_CLASS, mxREAL);
		sens_ptr = mxCreateNumericArray(1, d, mxSINGLE_CLASS, mxREAL);
	}
	else if (type == 2 && (inputScalars.atomic_32bit || inputScalars.atomic_64bit)) {
		if (no_norm == 0)
			if (inputScalars.atomic_32bit)
				sens_ptr = mxCreateNumericArray(1, mDim, mxINT32_CLASS, mxREAL);
			else
				sens_ptr = mxCreateNumericArray(1, mDim, mxINT64_CLASS, mxREAL);
		else
			if (inputScalars.atomic_32bit)
				sens_ptr = mxCreateNumericArray(1, d, mxINT32_CLASS, mxREAL);
			else
				sens_ptr = mxCreateNumericArray(1, d, mxINT64_CLASS, mxREAL);
		if (inputScalars.atomic_32bit)
			array_ptr = mxCreateNumericArray(1, mDim, mxINT32_CLASS, mxREAL);
		else
			array_ptr = mxCreateNumericArray(1, mDim, mxINT64_CLASS, mxREAL);
	}
	else {
		array_ptr = mxCreateNumericArray(1, mDim, mxSINGLE_CLASS, mxREAL);
		if (no_norm == 0 && type == 2)
			sens_ptr = mxCreateNumericArray(1, mDim, mxSINGLE_CLASS, mxREAL);
		else
			if (type == 0 && inputScalars.atomic_32bit)
				sens_ptr = mxCreateNumericArray(1, d, mxINT32_CLASS, mxREAL);
			else if (type == 0 && inputScalars.atomic_64bit)
				sens_ptr = mxCreateNumericArray(1, mDim, mxINT64_CLASS, mxREAL);
			else
				sens_ptr = mxCreateNumericArray(1, d, mxSINGLE_CLASS, mxREAL);
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
	if (inputScalars.listmode)
		w_vec.listCoord = getSingles(options, "x", 0);

	const float* Sino = getSingles(Sin, "solu");
	const float* randoms = getSingles(sc_ra, "solu");
	const float* extraCorr = getSingles(options, "ScatterC", 0);
	const float* x0 = getSingles(options, "x0");
	if (DEBUG) {
		mexPrintBase("x0[0] = %f\n", x0[0]);
		mexPrintBase("x0.dim = %u\n", mxGetNumberOfElements(mxGetField(options, 0, "x0")));
		mexEval();
	}

	if (DEBUG) {
		mexPrint("Pointers set");
	}
	if (inputScalars.atomic_32bit && (type == 2)) {
		int32_t* output = getInt32s(array_ptr, "solu");
		int32_t* sensIm = getInt32s(sens_ptr, "solu");
		reconstruction_metal(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else if (inputScalars.atomic_64bit && (type == 2)) {
		int64_t* output = getInt64s(array_ptr, "solu");
		int64_t* sensIm = getInt64s(sens_ptr, "solu");
		reconstruction_metal(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else if (inputScalars.atomic_64bit && (type == 0)) {
		float* output = getSingles(array_ptr, "solu");
		int64_t* sensIm = getInt64s(sens_ptr, "solu");
		reconstruction_metal(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else if (inputScalars.atomic_32bit && (type == 0)) {
		float* output = getSingles(array_ptr, "solu");
		int32_t* sensIm = getInt32s(sens_ptr, "solu");
		reconstruction_metal(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}
	else {
		float* output = getSingles(array_ptr, "solu");
		float* sensIm = getSingles(sens_ptr, "solu");
		reconstruction_metal(z_det, x, inputScalars, w_vec, MethodList, pituus, header_directory, Sino, x0, output, sensIm, type, no_norm, randoms, atten, norm, extraCorr, size_gauss, xy_index, z_index, L);
		plhs[0] = array_ptr;
		if (nlhs > 1)
			plhs[1] = sens_ptr;
		else
			mxDestroyArray(sens_ptr);
	}

	return;
}