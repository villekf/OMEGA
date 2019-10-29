/**************************************************************************
* These functions handle the device selection, queue creation, program
* building and kernel creation, as well as the output data and kernel 
* release.
*
* Copyright(C) 2019  Ville - Veikko Wettenhovi
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
	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, 
	const float dz, const float bx, const float by, const float bz,	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, 
	const uint32_t* pituus, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, 
	mxArray* cell, const bool verbose, const uint32_t randoms_correction,  const uint32_t attenuation_correction, 
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets, 
	const float epps, const uint8_t* rekot,	const char* k_path, const size_t size_rekot, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, 
	const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, 
	const uint32_t device, float kerroin, const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, 
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, 
	const char* header_directory, const bool precompute, const int32_t dec, const uint16_t n_rays, const float cr_pz, const bool use_64bit_atomics) {

	// Total number of voxels
	const uint32_t im_dim = Nx * Ny * Nz;

	cl_int status = CL_SUCCESS;
	cl_context context = NULL;
	cl_uint num_devices_context;
	cl_kernel kernel;
	cl_kernel kernel_sum;
	cl_kernel kernel_mlem;
	int cpu_device = -1;
	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	uint32_t Nxyz = Nx * Ny * Nz;

	size_t size;

	// Maximum of 16 devices
	cl_device_id * devices = (cl_device_id*)alloca(16 * sizeof(cl_device_id));

	// Get the OpenCL context
	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, subsets,
		raw);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		mexPrintf("Failed to get platforms\n");
		return;
	}

	// Create the same number of command queues as there are devices
	cl_program program = NULL;
	cl_command_queue *commandQueues = (cl_command_queue*)alloca(sizeof(cl_command_queue) * num_devices_context);

	// Build the program and get the command queues
	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, projector_type, header_directory);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		//status = clReleaseProgram(program);
		for (size_t i = 0ULL; i < num_devices_context; i++) {
			//status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		mexPrintf("Failed to build programs\n");
		return;
	}

	// Create kernels
	// Orthogonal distance based
	if (projector_type == 2u) {
		kernel = clCreateKernel(program, "orth_multi", &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			for (size_t i = 0ULL; i < num_devices_context; i++) {
				status = clReleaseCommandQueue(commandQueues[i]);
				status = clReleaseDevice(devices[i]);
			}
			status = clReleaseContext(context);
			status = clReleaseProgram(program);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}
	// Improved Siddon's ray tracer
	else {
		kernel = clCreateKernel(program, "siddon_multi", &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			for (size_t i = 0ULL; i < num_devices_context; i++) {
				status = clReleaseCommandQueue(commandQueues[i]);
				status = clReleaseDevice(devices[i]);
			}
			status = clReleaseContext(context);
			status = clReleaseProgram(program);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}

	// Kernel for the summing (combination) of data from different devices
	kernel_sum = clCreateKernel(program, "summa", &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		for (size_t i = 0ULL; i < num_devices_context; i++) {
			status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		status = clReleaseContext(context);
		status = clReleaseProgram(program);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	// MLEM/OSEM kernel
	kernel_mlem = clCreateKernel(program, "mlem", &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		for (size_t i = 0ULL; i < num_devices_context; i++) {
			status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		status = clReleaseContext(context);
		status = clReleaseProgram(program);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}
	if (verbose) {
		mexPrintf("OpenCL kernels successfully built\n");
		mexEvalString("pause(.0001);");
	}

	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	// Compute the estimates
	OSEM_MLEM(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, options, dx, dy, dz, bx, 
		by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, 
		normalization, atten, size_atten, norm, size_norm, subsets, epps, Nt, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel, kernel_sum, 
		kernel_mlem, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, atomic_64bit, 
		compute_norm_matrix, precompute, dec, projector_type, n_rays, cr_pz, cell, osem_bool);


	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	// Release queues and devices
	for (size_t i = 0ULL; i < num_devices_context; i++) {
		status = clReleaseCommandQueue(commandQueues[i]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
		}
		status = clReleaseDevice(devices[i]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
		}
	}

	// Release context
	status = clReleaseContext(context);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}

	// Release program
	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}

	// Release kernels
	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
	status = clReleaseKernel(kernel_sum);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
	status = clReleaseKernel(kernel_mlem);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
}

// Main reconstruction function, forward/backward projection
void reconstruction_f_b_proj(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const char* k_path, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, 
	const char* fileName, const uint32_t device, float kerroin, float* output, float* normalizer, const size_t size_rhs, const bool no_norm, const size_t numel_x, 
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const int32_t dec, const uint16_t n_rays, const float cr_pz) {
	// This functions very similarly to the above function

	const uint32_t im_dim = Nx * Ny * Nz;

	cl_int status = CL_SUCCESS;
	cl_context context = NULL;
	cl_uint num_devices_context;
	cl_kernel kernel, kernel_sum;
	int cpu_device = -1;
	bool atomic_64bit = false;
	cl_uchar compute_norm_matrix = 1u;
	uint32_t Nxyz = Nx * Ny * Nz;

	size_t size;

	cl_device_id* devices = (cl_device_id*)alloca(16 * sizeof(cl_device_id));

	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, 1u,
		raw);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		return;
	}

	cl_program program = NULL;
	cl_command_queue* commandQueues = (cl_command_queue*)alloca(sizeof(cl_command_queue) * num_devices_context);

	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, projector_type, header_directory);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		status = clReleaseProgram(program);
		for (size_t i = 0; i < num_devices_context; i++) {
			status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		return;
	}

	if (projector_type == 2u) {
		kernel = clCreateKernel(program, "f_b_project_orth", &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			for (size_t i = 0ULL; i < num_devices_context; i++) {
				status = clReleaseCommandQueue(commandQueues[i]);
				status = clReleaseDevice(devices[i]);
			}
			status = clReleaseContext(context);
			status = clReleaseProgram(program);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}
	else {
		kernel = clCreateKernel(program, "f_b_project_siddon", &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			for (size_t i = 0ULL; i < num_devices_context; i++) {
				status = clReleaseCommandQueue(commandQueues[i]);
				status = clReleaseDevice(devices[i]);
			}
			status = clReleaseContext(context);
			status = clReleaseProgram(program);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}

	kernel_sum = clCreateKernel(program, "summa", &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		for (size_t i = 0ULL; i < num_devices_context; i++) {
			status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		status = clReleaseContext(context);
		status = clReleaseProgram(program);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	f_b_project(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, rhs, sc_ra, Nx, Ny, Nz, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, 
		normalization, atten, size_atten, norm, size_norm, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel_sum, kernel, output, normalizer, 
		size_rhs, no_norm, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, precompute, dec, 
		projector_type, n_rays, cr_pz);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	for (size_t i = 0; i < num_devices_context; i++)
		status = clReleaseCommandQueue(commandQueues[i]);

	status = clReleaseContext(context);

	status = clReleaseProgram(program);

	status = clReleaseKernel(kernel_sum);
	status = clReleaseKernel(kernel);
}


// Main reconstruction function, find the number of voxels each LOR traverses
void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, 
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, 
	const float NSlices, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par, const char* k_path, 
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const char* fileName, 
	const uint32_t device, float kerroin, const size_t numel_x, const char* header_directory) {

	// Number of voxels
	const uint32_t im_dim = Nx * Ny * Nz;
	// Never use 64 bit atomics here
	bool atomic_64bit = false;

	cl_int status = CL_SUCCESS;
	cl_context context = NULL;
	cl_uint num_devices_context;
	cl_kernel kernel;

	cl_device_id * devices = (cl_device_id*)alloca(sizeof(cl_device_id));

	// Get the context for a single device
	status = clGetPlatformsContextSingle(device, context, num_devices_context, devices);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		mexPrintf("Error while getting platforms\n");
		return;
	}

	// Get queues
	cl_program program = NULL;
	cl_command_queue *commandQueues = (cl_command_queue*)alloca(sizeof(cl_command_queue) * num_devices_context);

	status = ClBuildProgramGetQueues(program, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, 50u, header_directory);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		//status = clReleaseProgram(program);
		for (size_t i = 0; i < num_devices_context; i++) {
			//status = clReleaseCommandQueue(commandQueues[i]);
			status = clReleaseDevice(devices[i]);
		}
		mexPrintf("Error while building or getting queues\n");
		return;
	}

	kernel = clCreateKernel(program, "siddon_precomp", &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	// Compute the voxel count
	precomp_siddon(num_devices_context, context, commandQueues, lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, 
		kernel, numel_x);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	for (size_t i = 0; i < num_devices_context; i++) {
		status = clReleaseCommandQueue(commandQueues[i]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
		}
		status = clReleaseDevice(devices[i]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
		}
	}

	status = clReleaseContext(context);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}

	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}

	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
}
