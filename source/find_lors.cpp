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
	const float NSlices, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, const char* k_path,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z,
	const char* fileName, const uint32_t device, const size_t numel_x, const char* header_directory) {

	af::setDevice(device);

	if (DEBUG) {
		mexPrintf("Started\n");
		mexEvalString("pause(.0001);");
	}

	const uint32_t im_dim = Nx * Ny * Nz;
	bool atomic_64bit = false;
	const uint64_t local_size = 64ULL;

	cl_int status = CL_SUCCESS;
	cl::Kernel kernel;
	cl_uint num_devices_context = 1u;

	cl::Context context(afcl::getContext(true));
	std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>(&status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	else if (DEBUG) {
		mexPrintf("Context created\n");
		mexPrintf("devices.size() = %u\n", devices.size());
		mexEvalString("pause(.0001);");
	}
	cl::Device af_device_id = devices[0];
	cl::CommandQueue af_queue(afcl::getQueue(true));

	if (DEBUG) {
		mexPrintf("Queue created\n");
		mexEvalString("pause(.0001);");
	}

	cl::Program program;
	std::vector<cl::CommandQueue> commandQueues;

	commandQueues.push_back(af_queue);


	std::string kernelFile = header_directory;
	std::string kernel_path;

	kernel_path = k_path;
	kernel_path += ".cl";
	// Load the source text file
	std::ifstream sourceFile(kernel_path.c_str());
	std::string contentF((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
	// Load the header text file
	std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
	std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
	std::string content = contentHeader + contentF;
	std::vector<std::string> testi;
	testi.push_back(content);
	cl::Program::Sources source(testi);
	program = cl::Program(context, source);
	std::string options = "-DFIND_LORS";
	if (raw == 1)
		options += " -DRAW";
	options += (" -DLOCAL_SIZE=" + std::to_string(local_size));
	options += " -DCAST=float"; 
	if (DEBUG) {
		mexPrintf("%s\n", options.c_str());
		mexEvalString("pause(.0001);");
	}
	//try {
		status = program.build(options.c_str());
		if (status == CL_SUCCESS) {
			mexPrintf("OpenCL program built\n");
			mexEvalString("pause(.0001);");
		}
	//}
	//catch (cl::Error& e) {
		else {
			mexPrintf("Failed to build OpenCL program.\n");
			std::vector<cl::Device> dev;
			context.getInfo(CL_CONTEXT_DEVICES, &dev);
			for (int ll = 0; ll < dev.size(); ll++) {
				cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
				if (status != CL_BUILD_ERROR)
					continue;
				std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
				std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
				mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
			}
			return;
		}
		//mexPrintf("%s\n", e.what());
	//}

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

	precomp_siddon(num_devices_context, context, commandQueues, lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, 
		size_x, TotSinos, verbose, loop_var_par, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel, numel_x, local_size);


	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	return;
}
