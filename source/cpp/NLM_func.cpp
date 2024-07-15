

#include "NLM.h"
#include "mexFunktio.h"

using namespace std;


void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[])

{
	// Check for the number of input and output arguments
	if (nrhs < 19)
		mexErrMsgTxt("Too few input arguments. There must be at least 19.");
	else if (nrhs > 19)
		mexErrMsgTxt("Too many input arguments. There can be at most 19.");

	if (nlhs > 1 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most one.");

	int ind = 0;
	// Load the input arguments

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const double* u_ref = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* u_ref = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const double* u = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* u = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const double* gaussian = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* gaussian = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	const int32_t search_window_x = getScalarInt32(prhs[ind], ind);
	ind++;

	const int32_t search_window_y = getScalarInt32(prhs[ind], ind);
	ind++;

	const int32_t search_window_z = getScalarInt32(prhs[ind], ind);
	ind++;

	const int32_t patch_window_x = getScalarInt32(prhs[ind], ind);
	ind++;

	const int32_t patch_window_y = getScalarInt32(prhs[ind], ind);
	ind++;

	const int32_t patch_window_z = getScalarInt32(prhs[ind], ind);
	ind++;

	const uint32_t Nx = getScalarUInt32(prhs[ind], ind);
	ind++;

	const uint32_t Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	const uint32_t Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	const double h = getScalarDouble(prhs[ind], ind);
	ind++;

	const int32_t type = getScalarInt32(prhs[ind], ind);
	ind++;

	const double gamma = getScalarDouble(prhs[ind], ind);
	ind++;

	const double epps = getScalarDouble(prhs[ind], ind);
	ind++;

	const double p = getScalarDouble(prhs[ind], ind);
	ind++;

	const double q = getScalarDouble(prhs[ind], ind);
	ind++;

	const double c = getScalarDouble(prhs[ind], ind);
	ind++;

	const uint32_t N = Nx * Ny * Nz;
	const int32_t Nxy = static_cast<int32_t>(Nx * Ny);

	plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	double* grad = (double*)mxGetDoubles(plhs[0]);
#else
	double* grad = (double*)mxGetData(plhs[0]);
#endif

	NLMFunc(grad, u_ref, u, gaussian, search_window_x, search_window_y, search_window_z, patch_window_x, patch_window_y, patch_window_z, Nx,
		Ny, Nz, Nxy, h, type, gamma, epps, p, q, c);

	return;
}