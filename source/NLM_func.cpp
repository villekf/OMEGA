

#include "NLM.h"
#include "mexFunktio.h"

using namespace std;


void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[])

{
	// Check for the number of input and output arguments
	if (nrhs < 15)
		mexErrMsgTxt("Too few input arguments. There must be at least 15.");
	else if (nrhs > 15)
		mexErrMsgTxt("Too many input arguments. There can be at most 15.");

	if (nlhs > 1 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most one.");

	int ind = 0;
	// Load the input arguments

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* u_ref = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* u_ref = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* u = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* u = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
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

	const double epps = getScalarDouble(prhs[ind], ind);
	ind++;

	const uint32_t N = Nx * Ny * Nz;
	const int32_t Nxy = static_cast<int32_t>(Nx * Ny);

	plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	double* grad = (double*)mxGetDoubles(plhs[0]);
#else
	double* grad = (double*)mxGetData(plhs[0]);
#endif

	NLM(grad, u_ref, u, gaussian, search_window_x, search_window_y, search_window_z, patch_window_x, patch_window_y, patch_window_z, Nx,
		Ny, Nz, Nxy, h, type, epps);

	return;
}