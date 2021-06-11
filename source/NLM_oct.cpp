

#include "NLM.h"
#include <octave/oct.h>

using namespace std;


DEFUN_DLD(NLM_oct, prhs, nargout, "NLM") {

	int ind = 0;
	// Load the input arguments

	const NDArray u_ref_ = prhs(ind).array_value();
	ind++;

	const NDArray u_ = prhs(ind).array_value();
	ind++;

	const NDArray gaussian_ = prhs(ind).array_value();
	ind++;

	const int32_t search_window_x = prhs(ind).int32_scalar_value();
	ind++;

	const int32_t search_window_y = prhs(ind).int32_scalar_value();
	ind++;

	const int32_t search_window_z = prhs(ind).int32_scalar_value();
	ind++;

	const int32_t patch_window_x = prhs(ind).int32_scalar_value();
	ind++;

	const int32_t patch_window_y = prhs(ind).int32_scalar_value();
	ind++;

	const int32_t patch_window_z = prhs(ind).int32_scalar_value();
	ind++;

	const uint32_t Nx = prhs(ind).uint32_scalar_value();
	ind++;

	const uint32_t Ny = prhs(ind).uint32_scalar_value();
	ind++;

	const uint32_t Nz = prhs(ind).uint32_scalar_value();
	ind++;

	const double h = prhs(ind).scalar_value();
	ind++;

	const int32_t type = prhs(ind).int32_scalar_value();
	ind++;

	const double epps = prhs(ind).scalar_value();
	ind++;

	const uint32_t N = Nx * Ny * Nz;
	const int32_t Nxy = static_cast<int32_t>(Nx * Ny);

	NDArray grad_(dim_vector(N, 1));

	double* grad = grad_.fortran_vec();

	const double* u_ref = u_ref_.fortran_vec();
	const double* gaussian = gaussian_.fortran_vec();
	const double* u = u_.fortran_vec();

	NLM(grad, u_ref, u, gaussian, search_window_x, search_window_y, search_window_z, patch_window_x, patch_window_y, patch_window_z, Nx,
		Ny, Nz, Nxy, h, type, epps);


	octave_value_list retval(nargout);

	retval(0) = octave_value(grad_);

	return retval;
}