/**************************************************************************
* These functions handle the device selection, queue creation, program
* building and kernel creation, as well as the output data and kernel 
* release.
*
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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
#include "functions_multigpu.hpp"

using namespace std;

// Main reconstruction function for implementation 3
void reconstruction_multigpu(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, 
	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, 
	const float dy, const float dz, const float bx, const float by, const float bz,	const float bzb, const float maxxx, const float maxyy, const float zmax, 
	const float NSlices, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, 
	const uint32_t TotSinos, mxArray* cell, const bool verbose, const uint32_t randoms_correction,  const uint32_t attenuation_correction, 
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets, 
	const float epps, const uint8_t* rekot,	const char* k_path, const size_t size_rekot, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, 
	const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool use_psf, 
	const uint32_t device, float kerroin, const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, 
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, 
	const char* header_directory, const bool precompute, const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, 
	const bool use_64bit_atomics, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V, const size_t size_V, 
	const size_t local_size, const float* gaussian, const size_t size_gauss, const bool TOF, const int64_t TOFSize, const float sigma_x, const float* TOFCenter, 
	const int64_t nBins) {

	// Total number of voxels
	const uint32_t im_dim = Nx * Ny * Nz;

	cl_int status = CL_SUCCESS;
	cl::Context context;
	cl_uint num_devices_context;
	cl::Kernel kernel;
	cl::Kernel kernel_sum;
	cl::Kernel kernel_mlem;
	cl::Kernel kernel_3Dconvolution;
	cl::Kernel kernel_3Dconvolution_f;
	cl::Kernel kernel_vectorMult, kernel_vectorDiv;
	int cpu_device = -1;
	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	const uint32_t Nxyz = Nx * Ny * Nz;
	const uint8_t fp = 0;

	size_t size;

	// Maximum of 16 devices

	cl::vector<cl::Device> devices;

	// Get the OpenCL context
	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, subsets,
		raw);

	if (status != CL_SUCCESS) {
		mexPrintf("Failed to get platforms\n");
		return;
	}

	// Create the same number of command queues as there are devices
	cl::Program program;
	std::vector<cl::CommandQueue> commandQueues;


	const uint32_t scatter = static_cast<uint32_t>((bool)mxGetScalar(mxGetField(options, 0, "scatter")));

	const bool listmode = (bool)mxGetScalar(mxGetField(options, 0, "listmode"));

	// Build the program and get the command queues
	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, 
		projector_type, header_directory, crystal_size_z, precompute, raw, attenuation_correction, normalization, dec, fp, local_size, n_rays, n_rays3D, 
		false, cr_pz, dx, use_psf, scatter, randoms_correction, TOF, nBins, listmode);

	if (status != CL_SUCCESS) {
		mexPrintf("Failed to build programs\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Program created\n");
		mexEvalString("pause(.0001);");
	}

	// Create kernels
	// Orthogonal distance based
	if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1)))) {
		kernel = cl::Kernel(program, "kernel_multi", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create orthogonal OpenCL kernel\n");
			return;
		}
	}
	// Improved Siddon's ray tracer
	else {
		kernel = cl::Kernel(program, "siddon_multi", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create Siddon OpenCL kernel\n");
			return;
		}
	}

	if (use_psf) {
		kernel_3Dconvolution = cl::Kernel(program, "Convolution3D", &status);
		kernel_3Dconvolution_f = cl::Kernel(program, "Convolution3D_f", &status);
		kernel_vectorMult = cl::Kernel(program, "vectorMult", &status);
		kernel_vectorDiv = cl::Kernel(program, "vectorDiv", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create convolution OpenCL kernel\n");
			return;
		}
	}

	// Kernel for the summing (combination) of data from different devices
	kernel_sum = cl::Kernel(program, "summa", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	// MLEM/OSEM kernel
	kernel_mlem = cl::Kernel(program, "mlem", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OSEM/MLEM OpenCL kernel\n");
		return;
	}
	if (DEBUG) {
		mexPrintf("OpenCL kernels successfully built\n");
		mexEvalString("pause(.0001);");
	}

	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	// Compute the estimates
	OSEM_MLEM(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, options, dx, dy, dz, bx, 
		by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, 
		normalization, atten, size_atten, norm, size_norm, subsets, epps, Nt, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel, kernel_sum, 
		kernel_mlem, kernel_3Dconvolution, kernel_3Dconvolution_f, kernel_vectorMult, kernel_vectorDiv, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, atomic_64bit,
		compute_norm_matrix, precompute, dec, projector_type, n_rays, n_rays3D, cr_pz, cell, osem_bool, global_factor, bmin, bmax, Vmax, V, size_V, local_size, 
		use_psf, gaussian, size_gauss, scatter, TOF, TOFSize, sigma_x, TOFCenter, nBins, devices);


	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}
	return;
}

// Main reconstruction function, forward/backward projection
void reconstruction_f_b_proj(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const int64_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const char* k_path, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, 
	const char* fileName, const uint32_t device, float kerroin, mxArray* output, const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x,
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, const mxArray* Sin, const bool use_64bit_atomics, const float global_factor, 
	const float bmin, const float bmax, const float Vmax, const float* V, const size_t size_V, const size_t local_size, const bool use_psf, const mxArray* options, 
	const bool TOF, const int64_t TOFSize, const float sigma_x, const float* TOFCenter, const int64_t nBins) {
	// This functions very similarly to the above function

	const uint32_t im_dim = Nx * Ny * Nz;

	cl_int status = CL_SUCCESS;
	cl::Context context;
	cl_uint num_devices_context;
	cl::Kernel kernel, kernel_sum;
	int cpu_device = -1;
	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	uint32_t Nxyz = Nx * Ny * Nz;

	// If 1, then the forward projection is computed
	uint8_t fp = size_rhs == im_dim;
	if (fp == 0)
		fp = 2u;
	if (atomic_64bit && fp == 1)
		atomic_64bit = false;

	size_t size;

	cl::vector<cl::Device> devices;

	const bool listmode = (bool)mxGetScalar(mxGetField(options, 0, "listmode"));

	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, 1u,
		raw);

	if (status != CL_SUCCESS) {
		return;
	}


	const uint32_t scatter = static_cast<uint32_t>((bool)mxGetScalar(mxGetField(options, 0, "scatter")));

	cl::Program program;
	std::vector<cl::CommandQueue> commandQueues;

	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, projector_type, header_directory, crystal_size_z, 
		precompute, raw, attenuation_correction, normalization, dec, fp, local_size, n_rays, n_rays3D, false, cr_pz, dx, use_psf, scatter, randoms_correction, TOF, nBins, listmode);

	if (status != CL_SUCCESS) {
		mexPrintf("Failed to build programs\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Program created\n");
		mexEvalString("pause(.0001);");
	}

	if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1)))) {
		kernel = cl::Kernel(program, "kernel_multi", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}
	else {
		kernel = cl::Kernel(program, "siddon_multi", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create Siddon OpenCL kernel\n");
			return;
		}
	}

	kernel_sum = cl::Kernel(program, "summa", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create sum OpenCL kernel\n");
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	f_b_project(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, rhs, sc_ra, Nx, Ny, Nz, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, 
		normalization, atten, size_atten, norm, size_norm, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel_sum, kernel, output,  
		size_rhs, no_norm, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, precompute, dec, 
		projector_type, n_rays, n_rays3D, cr_pz, Sin, atomic_64bit, global_factor, bmin, bmax, Vmax, V, size_V, fp, local_size, options, scatter, TOF, 
		TOFSize, sigma_x, TOFCenter, nBins);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}
	return;
}


// Main reconstruction function, find the number of voxels each LOR traverses
void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz, 
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t size_x, 
	const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, const char* k_path, const uint32_t* pseudos,
	const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const char* fileName, 
	const uint32_t device, const size_t numel_x, const char* header_directory, const size_t local_size) {

	// Number of voxels
	const uint32_t im_dim = Nx * Ny * Nz;
	// Never use 64 bit atomics here
	bool atomic_64bit = false;

	cl_int status = CL_SUCCESS;
	cl::Context context;
	cl_uint num_devices_context;
	cl::Kernel kernel;

	cl::vector<cl::Device> devices;

	// Get the context for a single device
	status = clGetPlatformsContextSingle(device, context, num_devices_context, devices);

	if (status != CL_SUCCESS) {
		mexPrintf("Error while getting platforms\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Platforms obtained\n");
		mexEvalString("pause(.0001);");
	}

	// Get queues
	cl::Program program;
	std::vector<cl::CommandQueue> commandQueues;

	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, 50u, header_directory, 0.f, false, raw, 
		0, 0, 0, 0, local_size, 0, 0, true, 0.f, dx, false, 0, 0, false, 1LL);

	if (status != CL_SUCCESS) {
		mexPrintf("Error while building or getting queues\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Queues obtained\n");
		mexEvalString("pause(.0001);");
	}

	kernel = cl::Kernel(program, "siddon_precomp", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Kernel created\n");
		mexEvalString("pause(.0001);");
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	// Compute the voxel count
	precomp_siddon(num_devices_context, context, commandQueues, lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, 
		kernel, numel_x, local_size);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	return;
}
