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
	const mxArray* sc_ra, const mxArray* options, scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index, 
	const uint16_t* z_index, const uint32_t TotSinos, mxArray* cell, const bool verbose, const float* atten, const size_t size_atten, const float* norm, 
	const size_t size_norm, const uint8_t* rekot, const char* k_path, const size_t size_rekot, const uint32_t* pseudos, const uint16_t* L, 
	const bool osem_bool, const char* fileName, const uint32_t device, float kerroin, const size_t numel_x, const float* x_center, const float* y_center, 
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const char* header_directory, 
	const bool use_64bit_atomics, const float* V, const size_t size_V, size_t local_size[], const float* gaussian, const size_t size_gauss, 
	const int64_t TOFSize, const float* TOFCenter) {

	// Total number of voxels
	const uint32_t im_dim = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;

	cl_int status = CL_SUCCESS;
	cl::Context context;
	cl_uint num_devices_context;
	cl::Kernel kernel;
	cl::Kernel kernel_sum;
	cl::Kernel kernel_mlem;
	cl::Kernel kernel_3Dconvolution;
	cl::Kernel kernel_3Dconvolution_f;
	cl::Kernel kernel_vectorMult, kernel_vectorDiv;
	cl::Kernel kernelBP;
	int cpu_device = -1;
	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	const uint32_t Nxyz = im_dim;

	size_t size;

	// Maximum of 16 devices

	cl::vector<cl::Device> devices;

	// Get the OpenCL context
	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, inputScalars);

	if (status != CL_SUCCESS) {
		mexPrintf("Failed to get platforms\n");
		return;
	}

	std::string deviceName = devices[0].getInfo<CL_DEVICE_VENDOR>(&status);
	std::string NV("NVIDIA Corporation");
	if (NV.compare(deviceName) == 0 && inputScalars.projector_type == 1 && local_size[1] == 0ULL)
		local_size[0] = 32ULL;

	// Create the same number of command queues as there are devices
	cl::Program program, programAux;
	std::vector<cl::CommandQueue> commandQueues;


	inputScalars.scatter = static_cast<uint32_t>((bool)mxGetScalar(getField(options, 0, "scatter")));

	const uint8_t listmode = (uint8_t)mxGetScalar(getField(options, 0, "listmode"));
	const bool computeSensImag = (bool)mxGetScalar(getField(options, 0, "compute_sensitivity_image"));
	const bool CT = (bool)mxGetScalar(getField(options, 0, "CT"));
	const bool atomic_32bit = (bool)mxGetScalar(getField(options, 0, "use_32bit_atomics"));

	if (listmode == 1 && computeSensImag)
		compute_norm_matrix = 0u;

	// Build the program and get the command queues
	status = ClBuildProgramGetQueues(program, programAux, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, atomic_32bit,
		inputScalars, header_directory, local_size,	false, listmode, CT);

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
	if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u || 
		(inputScalars.projector_type == 1u && ((inputScalars.precompute || (inputScalars.n_rays * inputScalars.n_rays3D) == 1)))) {
		kernel = cl::Kernel(program, "proj123SiddonSingleRay", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create orthogonal OpenCL kernel\n");
			return;
		}
	}
	else if (inputScalars.projector_type == 4) {
		//if (inputScalars.fp == 1 && !inputScalars.SPECT)
			kernel = cl::Kernel(program, "projectorType4Forward", &status);
		//else
			kernelBP = cl::Kernel(program, "projectorType4Backward", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create projector type 4 OpenCL kernels\n");
			return;
		}
	}
	else if (inputScalars.projector_type == 5) {
		//if (inputScalars.fp == 1 && !inputScalars.SPECT)
			kernel = cl::Kernel(program, "projectorType5Forward", &status);
		//else
			kernelBP = cl::Kernel(program, "projectorType5Backward", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create projector type 5 OpenCL kernels\n");
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

	if (inputScalars.use_psf) {
		kernel_3Dconvolution = cl::Kernel(programAux, "Convolution3D", &status);
		kernel_3Dconvolution_f = cl::Kernel(programAux, "Convolution3D_f", &status);
		kernel_vectorMult = cl::Kernel(programAux, "vectorMult", &status);
		kernel_vectorDiv = cl::Kernel(programAux, "vectorDiv", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create convolution OpenCL kernel\n");
			return;
		}
	}

	// Kernel for the summing (combination) of data from different devices
	kernel_sum = cl::Kernel(programAux, "summa", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	// MLEM/OSEM kernel
	kernel_mlem = cl::Kernel(programAux, "mlem", &status);
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
	OSEM_MLEM(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, Sin, sc_ra, inputScalars, options, 
		pituus, koko_l, xy_index, z_index, TotSinos, verbose, atten, size_atten, norm, size_norm, pseudos, L, im_dim, kernel, kernel_sum, 
		kernel_mlem, kernel_3Dconvolution, kernel_3Dconvolution_f, kernel_vectorMult, kernel_vectorDiv, kernelBP, numel_x, x_center, y_center, 
		z_center, size_center_x, size_center_y, size_center_z, atomic_64bit, atomic_32bit, compute_norm_matrix, cell, osem_bool, V, size_V, 
		local_size, gaussian, size_gauss, TOFSize, TOFCenter, devices);


	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}
	return;
}

// Main reconstruction function, forward/backward projection
void reconstruction_f_b_proj(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, 
	const mxArray* sc_ra, scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, 
	const uint32_t TotSinos, const bool verbose, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const char* k_path, const uint32_t* pseudos, const uint16_t* L, const char* fileName, const uint32_t device, float kerroin, mxArray* output, 
	const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x, const float* x_center, const float* y_center, const float* z_center, 
	const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const char* header_directory, const mxArray* Sin, 
	const bool use_64bit_atomics, const float* V, const size_t size_V, size_t local_size[], const mxArray* options, const int64_t TOFSize, 
	const float* TOFCenter) {
	// This functions very similarly to the above function

	const uint32_t im_dim = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;
	if (DEBUG) {
		mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
		mexPrintf("inputScalars.Nx = %u\n", inputScalars.Nx);
		mexPrintf("inputScalars.Ny = %u\n", inputScalars.Ny);
		mexPrintf("inputScalars.Nz = %u\n", inputScalars.Nz);
		mexEvalString("pause(.0001);");
	}

	cl_int status = CL_SUCCESS;
	cl::Context context;
	cl_uint num_devices_context;
	cl::Kernel kernel, kernel_sum;
	int cpu_device = -1;
	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	uint32_t Nxyz = im_dim;

	// If 1, then the forward projection is computed
	inputScalars.fp = size_rhs == im_dim;
	if (inputScalars.fp == 0)
		inputScalars.fp = 2u;
	if (atomic_64bit && inputScalars.fp == 1)
		atomic_64bit = false;

	size_t size;

	cl::vector<cl::Device> devices;

	const bool listmode = (bool)mxGetScalar(getField(options, 0, "listmode"));
	const bool CT = (bool)mxGetScalar(getField(options, 0, "CT"));
	const bool atomic_32bit = (bool)mxGetScalar(getField(options, 0, "use_32bit_atomics"));

	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, 
		Nxyz, inputScalars);

	std::string deviceName = devices[0].getInfo<CL_DEVICE_VENDOR>(&status);
	std::string NV("NVIDIA Corporation");
	if (NV.compare(deviceName) == 0 && inputScalars.projector_type == 1 && local_size[1] == 0ULL)
		local_size[0] = 32ULL;

	if (status != CL_SUCCESS) {
		return;
	}

	inputScalars.scatter = static_cast<uint32_t>((bool)mxGetScalar(getField(options, 0, "scatter")));

	cl::Program program, programAux;
	std::vector<cl::CommandQueue> commandQueues;

	status = ClBuildProgramGetQueues(program, programAux, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, 
		atomic_32bit, inputScalars, header_directory, local_size, false, listmode, CT);

	if (status != CL_SUCCESS) {
		mexPrintf("Failed to build programs\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Program created\n");
		mexEvalString("pause(.0001);");
	}
	if (DEBUG) {
		mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
		mexEvalString("pause(.0001);");
	}

	if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u || 
		(inputScalars.projector_type == 1u && ((inputScalars.precompute || (inputScalars.n_rays * inputScalars.n_rays3D) == 1)))) {
		kernel = cl::Kernel(program, "proj123SiddonSingleRay", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create OpenCL kernel\n");
			return;
		}
	}
	else if (inputScalars.projector_type == 4) {
		if (inputScalars.fp == 1 && !inputScalars.SPECT)
			kernel = cl::Kernel(program, "projectorType4Forward", &status);
		else
			kernel = cl::Kernel(program, "projectorType4Backward", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create projector type 4 OpenCL kernel\n");
			return;
		}
	}
	else if (inputScalars.projector_type == 5) {
		if (inputScalars.fp == 1 && !inputScalars.SPECT)
			kernel = cl::Kernel(program, "projectorType5Forward", &status);
		else
			kernel = cl::Kernel(program, "projectorType5Backward", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create projector type 5 OpenCL kernel\n");
			return;
		}
	}
	else {
		kernel = cl::Kernel(program, "proj1SiddonMultiRay", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create Siddon OpenCL kernel\n");
			return;
		}
	}

	kernel_sum = cl::Kernel(programAux, "summa", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create sum OpenCL kernel\n");
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}
	if (DEBUG) {
		mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
		mexEvalString("pause(.0001);");
	}

	f_b_project(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, rhs, sc_ra, inputScalars, pituus, 
		koko_l, xy_index, z_index, TotSinos, verbose, atten, size_atten, norm, size_norm, pseudos, L, im_dim, kernel_sum, kernel, output,  
		size_rhs, no_norm, numel_x, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, Sin, atomic_64bit, atomic_32bit, 
		V, size_V, local_size, options, TOFSize, TOFCenter);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}
	return;
}


// Main reconstruction function, find the number of voxels each LOR traverses
void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, 
	scalarStruct inputScalars, const char* k_path, const uint32_t* pseudos,	const uint16_t* L,const char* fileName, const uint32_t device, 
	const size_t numel_x, const char* header_directory, const size_t local_size[]) {

	// Number of voxels
	const uint32_t im_dim = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;
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
	cl::Program program, programAux;
	std::vector<cl::CommandQueue> commandQueues;

	status = ClBuildProgramGetQueues(program, programAux, k_path, context, num_devices_context, devices, verbose, commandQueues, atomic_64bit, 
		false, inputScalars, header_directory, 0, local_size, true);

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
	precomp_siddon(num_devices_context, context, commandQueues, lor, z_det, x, y, inputScalars, TotSinos, verbose, loop_var_par, pseudos, L, im_dim,
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
