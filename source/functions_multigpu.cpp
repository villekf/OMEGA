/**************************************************************************
* All the functions needed for the matrix-free OpenCL image reconstruction
* in the multi-GPU/device case. There are functions for context creation 
* for both multi-device and single device case (clGetPlatformsContextSingle) 
* and a function for program building and command queue creation 
* (ClBuildProgramGetQueues).
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

// Get the OpenCL context for the current platform
cl_int clGetPlatformsContext(const uint32_t device, const float kerroin, cl::Context& context, size_t& size, int& cpu_device, 
	cl_uint& num_devices_context, std::vector<cl::Device> &devices, bool& atomic_64bit, cl_uchar& compute_norm_matrix, const uint32_t Nxyz,
	const uint32_t subsets, const uint8_t raw) {
	cl_int status = CL_SUCCESS;
	cl_uint num_platforms;
	cl_float mem_portions;
	// These variables are used to determine if the normalization constant is memory-wise possible to store
	if (raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	cl_float image_bytes = static_cast<cl_float>(Nxyz) * 8.f;

	// Get the number of platforms 
	std::vector<cl::Platform> platforms;
	status = cl::Platform::get(&platforms);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}

	if (device > platforms.size() - 1) {
		std::cerr << "The specified platform number is greater than the available platform numbers!" << std::endl;
		status = -1;
		return status;
	}

	// Get context properties from the chosen platform
	cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, reinterpret_cast <cl_context_properties>(platforms[device]()), 0 };

	// Create context from the chosen platform
	if (kerroin == 0.f) {
		// If a single device was selected (options.cpu_to_gpu_factor = 0), use GPU if possible
		context = cl::Context(CL_DEVICE_TYPE_GPU, properties, NULL, NULL, &status);
		if (status != CL_SUCCESS) {
			// Otherwise CPU
			context = cl::Context(CL_DEVICE_TYPE_CPU, properties, NULL, NULL, &status);
		}
	}
	// Use all devices if options.cpu_to_gpu_factor > 0
	else
		context = cl::Context(CL_DEVICE_TYPE_ALL, properties, NULL, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	// Get device IDs
	std::vector<cl::Device> devices2;
	status = context.getInfo(CL_CONTEXT_DEVICES, &devices2);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}

	// Ignore devices with less than 2GB of memory
	// Also check global memory size and the feasibility of storing the normalization constants
	cl_int n_ignores = 0;
	cl_int n_gpus = 0;
	std::vector<cl_int> ignores(devices2.size());
	if (devices2.size() > 1ULL) {
		cl_ulong mem_max = 0ULL;
		for (size_t i = 0ULL; i < devices2.size(); i++)
		{
			cl_device_type type;
			cl_ulong mem;
			status = devices2[i].getInfo(CL_DEVICE_TYPE, &type);
			switch (type) {
			case CL_DEVICE_TYPE_GPU:
				status = devices2[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &mem);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				if ((mem / (1024ULL * 1024ULL)) < 2000ULL || (kerroin == 0.f && mem < mem_max)) {
					ignores[i] = 1;
					n_ignores++;
				}
				else {
					ignores[i] = 0;
					mem_max = mem;
					if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
						compute_norm_matrix = 0u;
					else
						compute_norm_matrix = 1u;
				}
				n_gpus++;
				break;
			case CL_DEVICE_TYPE_CPU:
				status = devices2[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &mem);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				cpu_device = static_cast<int>(i);
				if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
					compute_norm_matrix = 0u;
				else
					compute_norm_matrix = 1u;
				break;
			case CL_DEVICE_TYPE_ACCELERATOR:
				break;
			}
		}
	}
	else {
		cl_ulong mem;
		status = devices2[0].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &mem);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
			compute_norm_matrix = 0u;
		else
			compute_norm_matrix = 1u;
	}

	// Get the number of devices
	status = context.getInfo(CL_CONTEXT_NUM_DEVICES, &num_devices_context);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}

	// Remove devices that were ignored above
	uint32_t ll = 0u;
	for (size_t i = 0; i < num_devices_context; i++) {
		if (ignores[i] == 1) {
			num_devices_context--;
		}
		else {
			devices.push_back(devices2[i]);
		}
	}

	return status;
}

// Get context for a single device
cl_int clGetPlatformsContextSingle(const uint32_t device, cl::Context& context, cl_uint& num_devices_context, std::vector<cl::Device>& devices) {
	cl_int status = CL_SUCCESS;
	cl_uint num_platforms = 0;

	cl_uint num_devices = 0;
	std::vector<cl::Platform> platforms;
	status = cl::Platform::get(&platforms);

	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}

	if (device > platforms.size() - 1 || device < 0) {
		std::cerr << "The specified platform number is greater or smaller than the available platform numbers!" << std::endl;
		return -1;
	}

	platforms[device].getInfo(CL_DEVICE_TYPE_ALL, &num_devices);

	std::vector<cl::Device> devices2;
	status = platforms[device].getDevices(CL_DEVICE_TYPE_ALL, &devices2);

	// Choose the GPU with highest amount of memory
	// If no GPU, use CPU
	if (num_devices > 1) {
		cl_ulong mem_prev = 0;
		for (size_t i = 0; i < num_devices; i++) {
			cl_device_type type;
			status = devices2[i].getInfo(CL_DEVICE_TYPE, &type);
			switch (type) {
			case CL_DEVICE_TYPE_GPU:
				cl_ulong mem;
				status = devices2[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &mem);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				if (mem_prev < mem) {
					devices[0] = devices2[i];
					mem_prev = mem;
				}
				break;
			case CL_DEVICE_TYPE_CPU:
				if (mem_prev == 0)
					devices[0] = devices2[i];
				break;
			case CL_DEVICE_TYPE_ACCELERATOR:
				break;
			}
		}
	}
	else
		devices[0] = devices2[0];


	cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, reinterpret_cast <cl_context_properties>(platforms[device]()), 0 };

	context = cl::Context(devices, properties, NULL, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}

	num_devices_context = 1;

	return CL_SUCCESS;
}

// Build the programs and get the command queues
cl_int ClBuildProgramGetQueues(cl::Program& program, const char* k_path, const cl::Context context, const cl_uint num_devices_context,
	const std::vector<cl::Device>& devices, const bool verbose, std::vector<cl::CommandQueue>& commandQueues, bool& atomic_64bit, const uint32_t projector_type, const char* header_directory,
	const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction, const uint32_t normalization_correction,
	const int32_t dec, const uint8_t fp, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D, const bool find_lors, const float cr_pz,
	const float dx, const bool use_psf, const uint32_t scatter, const uint32_t randoms_correction) {
	cl_int status = CL_SUCCESS;


	//cl_device_id* devices2 = (cl_device_id*)alloca(num_devices_context * sizeof(cl_device_id));
	std::vector<cl::Device> devices2;
	for (size_t i = 0; i < num_devices_context; i++) {
		devices2.push_back(devices[i]);
	}

	std::string options = header_directory;
	options += " -cl-single-precision-constant";
	if (crystal_size_z == 0.f && projector_type == 2u)
		options += " -DCRYST";
	if ((crystal_size_z > 0.f && projector_type == 2u) || projector_type == 3u)
		options += " -DCRYSTZ";
	if (projector_type == 3u)
		options += " -DVOL";
	if (precompute)
		options += " -DPRECOMPUTE";
	if (raw == 1)
		options += " -DRAW";
	if (projector_type == 1u)
		options += " -DSIDDON";
	else if (projector_type == 2u || projector_type == 3u)
		options += " -DORTH";
	if (attenuation_correction == 1u)
		options += " -DATN";
	if (normalization_correction == 1u)
		options += " -DNORM";
	if ((projector_type == 2 || projector_type == 3u || (projector_type == 1u && use_psf && (precompute || (n_rays * n_rays3D) == 1))) && dec > 0)
		options += (" -DDEC=" + std::to_string(dec));
	if (fp < 2)
		options += " -DFP";
	if (fp != 1)
		options += " -DBP";
	if (projector_type == 1u && !precompute && (n_rays * n_rays3D) > 1) {
		options += (" -DN_RAYS=" + std::to_string(n_rays * n_rays3D));
		options += (" -DN_RAYS2D=" + std::to_string(n_rays));
		options += (" -DN_RAYS3D=" + std::to_string(n_rays3D));
	}
	if (find_lors)
		options += " -DFIND_LORS";
	if (use_psf)
		options += " -DPSF";
	if (scatter == 1u)
		options += " -DSCATTER";
	if (randoms_correction == 1u)
		options += " -DRANDOMS";
	//if (projector_type == 1u && use_psf && (precompute || (n_rays * n_rays3D) == 1)) {
	//	options += " -DORTH";
	//	options += " -DCRYSTZ";
	//	uint32_t limit = static_cast<uint32_t>(std::ceil(cr_pz / dx));
	//	options += (" -DPSF_LIMIT=" + std::to_string(limit));
	//	options += (" -DX=" + std::to_string(dx));
	//	options += (" -DSIGMA=" + std::to_string(cr_pz));
	//}
	options += (" -DLOCAL_SIZE=" + std::to_string(local_size));
	size_t pituus;
	if (atomic_64bit) {
		pituus = options.length();
		options += " -DCAST=ulong";
		options += " -DATOMIC";
		options += (" -DTH=" + std::to_string(TH));
	}
	else
		options += " -DCAST=float";
	//mexPrintf("%s\n", options.c_str());
	// If integer atomic 64-bit operations are enabled, check if they are supported by the device(s)
	if (atomic_64bit) {
		std::string kernel_path_atom;

		kernel_path_atom = k_path;
		//kernel_path_atom += "_64atom.cl";
		kernel_path_atom += ".cl";
		// Load the source text file
		std::fstream sourceFile_atom(kernel_path_atom.c_str());
		std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
		std::vector<std::string> testi;
		testi.push_back(content_atom);
		cl::Program::Sources source(testi);
		// Create the program from the source
		// Build the program
		program = cl::Program(context, source);
		//try {
			status = program.build(options.c_str());
			if (status != CL_SUCCESS) {
				mexPrintf("Failed to build 64-bit atomics program.\n");
			}
			else {
				mexPrintf("OpenCL program (64-bit atomics) built\n");
			}
		//}
		//catch (cl::Error& e) {
		//	//mexPrintf("%s\n", e.what());
		//	options.erase(pituus, options.size() + 1);
		//	options += " -DCAST=float";
		//	mexPrintf("Failed to build 64-bit atomics program.\n");
		//	status = -1;
		//}
	}
	else
		status = -1;
	// If not, use 32-bit atomic add (float)
	if (status != CL_SUCCESS) {
		status = CL_SUCCESS;
		atomic_64bit = false;

		std::string kernel_path;

		kernel_path = k_path;
		kernel_path += ".cl";
		std::fstream sourceFile(kernel_path.c_str());
		std::string content((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
		std::vector<std::string> testi;
		testi.push_back(content);
		cl::Program::Sources source(testi);
		program = cl::Program(context, source);
		//try {
			status = program.build(options.c_str());
			if (status == CL_SUCCESS) {
				mexPrintf("OpenCL program built\n");
			}
			else {
				getErrorString(status);
				mexPrintf("Failed to build OpenCL program.\n");
				std::vector<cl::Device> dev;
				context.getInfo(CL_CONTEXT_DEVICES, &dev);
				for (int ll = 0; ll < dev.size(); ll++) {
					cl_build_status b_status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
					if (b_status != CL_BUILD_ERROR)
						continue;
					std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
					std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
					mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
				}
				return status;

			}
		//}
		//catch (cl::Error& e) {
		//	mexPrintf("%s\n", e.what());
		//	status = -1;
		//}
	}

	// Create the command queues
	// Enable out of order execution (devices can compute kernels at the same time)
	for (size_t i = 0; i < num_devices_context; i++) {
		commandQueues.push_back(cl::CommandQueue(context, devices[i], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status));
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	return status;
}