/**************************************************************************
* Matrix free computations for OMEGA.
* In this file the OpenCL buffers are created, calls to other necessary
* functions are made and the OpenCL kernels are launched. This file
* contains the code for the matrix-free reconstructions in OMEGA using the
* implementation 2.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
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
#include "AF_opencl_functions.hpp"

// Use ArrayFire namespace for convenience
using namespace af;

void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx,
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
	const float NSlices, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par, const char* k_path,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z,
	const char* fileName, const uint32_t device, const size_t numel_x, const char* header_directory) {


	const uint32_t im_dim = Nx * Ny * Nz;
	bool atomic_64bit = false;
	const uint64_t local_size = 64ULL;

	cl_int status = CL_SUCCESS;
	cl_kernel kernel;
	cl_uint num_devices_context = 1u;

	cl_context context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue commandQueue = afcl::getQueue();

	cl_program program = NULL;
	cl_command_queue* commandQueues = (cl_command_queue*)alloca(sizeof(cl_command_queue) * num_devices_context);

	commandQueues[0] = commandQueue;

	std::string kernel_path;

	kernel_path = k_path;
	kernel_path += ".cl";

	std::fstream sourceFile(kernel_path.c_str());
	std::string content((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
	const char* sourceCode = new char[content.size()];
	sourceCode = content.c_str();
	program = clCreateProgramWithSource(context, 1, (const char**)& sourceCode, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	std::string options = header_directory;
	if (raw == 1)
		options += " -DRAW";
	options += " -DFIND_LORS";
	options += (" -DLOCAL_SIZE=" + std::to_string(local_size));
	options += " -DCAST=float";
	status = clBuildProgram(program, num_devices_context, &af_device_id, options.c_str(), NULL, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to build OpenCL program. Build log: \n");
		size_t len;
		char* buffer;
		clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		buffer = (char*)calloc(len, sizeof(size_t));
		clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
		mexPrintf("%s\n", buffer);
		return;
	}

	kernel = clCreateKernel(program, "siddon_precomp", &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	precomp_siddon(num_devices_context, context, commandQueues, lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, 
		size_x, TotSinos, verbose, loop_var_par, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel, numel_x, local_size);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		getErrorString(status);
	}

	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		getErrorString(status);
	}
	return;
}
