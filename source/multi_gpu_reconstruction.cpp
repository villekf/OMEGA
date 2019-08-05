#include "functions_multigpu.hpp"

using namespace std;

// Main reconstruction function
void reconstruction_multigpu(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, 
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint8_t* rekot,
	const char* k_path, const size_t size_rekot, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const uint32_t device, float kerroin, 
	const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, 
	const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, const int32_t dec) {


	const uint32_t im_dim = Nx * Ny * Nz;

	mxArray *mlem = mxCreateNumericArray(2, dimmi, mxSINGLE_CLASS, mxREAL);
	float *ele_ml = (float*)mxGetData(mlem);

	cl_int status = CL_SUCCESS;
	cl_context context = NULL;
	cl_uint num_devices_context;
	cl_kernel kernel;
	cl_kernel kernel_sum;
	cl_kernel kernel_mlem;
	int cpu_device = -1;
	bool atomic_64bit = false;
	cl_uchar compute_norm_matrix = 1u;
	uint32_t Nxyz = Nx * Ny * Nz;

	size_t size;

	cl_device_id * devices = (cl_device_id*)alloca(16 * sizeof(cl_device_id));

	status = clGetPlatformsContext(device, kerroin, context, size, cpu_device, num_devices_context, devices, atomic_64bit, compute_norm_matrix, Nxyz, subsets,
		raw);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		mexPrintf("Failed to get platforms\n");
		return;
	}

	cl_program program = NULL;
	cl_command_queue *commandQueues = (cl_command_queue*)alloca(sizeof(cl_command_queue) * num_devices_context);

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

	//size_t size_testi = 1u;
	//cl_kernel kernel_testi = clCreateKernel(program, "test", &status);
	//if (status != CL_SUCCESS) {
	//	std::cerr << getErrorString(status) << std::endl;
	//	mexPrintf("Failed to create OpenCL kernel\n");
	//	return;
	//}
	//cl_mem d_var = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ulong), NULL, &status);
	//cl_ulong atomic = 0u;
	//for (size_t i = 0; i < num_devices_context; i++) {
	//	status = clEnqueueWriteBuffer(commandQueues[i], d_var, CL_TRUE, 0, sizeof(cl_ulong), &atomic, 0, NULL, NULL);
	//	clFinish(commandQueues[i]);
	//	clSetKernelArg(kernel_testi, 0, sizeof(cl_mem), &d_var);
	//	status = clEnqueueNDRangeKernel(commandQueues[i], kernel_testi, 1, NULL, &size_testi, NULL, 0, NULL, NULL);
	//	clFinish(commandQueues[i]);
	//	status = clEnqueueReadBuffer(commandQueues[i], d_var, CL_TRUE, 0, sizeof(cl_ulong), &atomic, 0, NULL, NULL);
	//	clFinish(commandQueues[i]);
	//	mexPrintf("atomic = %u\n", atomic);
	//	if (atomic == 0)
	//		atomic_64bit = false;
	//	mexPrintf("atomic_64bit = %u\n", atomic_64bit);
	//}

	mexPrintf("atomic_64bit = %u\n", atomic_64bit);
	mexPrintf("compute_norm_matrix = %u\n", compute_norm_matrix);
	mexPrintf("randoms_correction = %u\n", randoms_correction);
	mexPrintf("normalization = %u\n", normalization);
	mexPrintf("projector_type = %u\n", projector_type);

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
		mexPrintf("size_center_x = %u\n", size_center_x);
	}
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

	OSEM_MLEM(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, options, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, normalization, atten, size_atten, norm, size_norm, 
		subsets, epps, Nt, pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel, kernel_sum, kernel_mlem, ele_ml, numel_x, tube_width, crystal_size_z, 
		x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, atomic_64bit, compute_norm_matrix, precompute, dec, projector_type);


	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

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
	status = clReleaseKernel(kernel_sum);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
	status = clReleaseKernel(kernel_mlem);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
	//clReleaseKernel(kernel_testi);
	//if (status != CL_SUCCESS) {
	//	std::cerr << getErrorString(status) << std::endl;
	//}
	//clReleaseMemObject(d_var);


	if (osem_bool)
		mxSetCell(cell, 1, mlem);
	else
		mxSetCell(cell, 0, mlem);
}

// Main reconstruction function
void reconstruction_f_b_proj(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const char* k_path, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, 
	const char* fileName, const uint32_t device, float kerroin, float* output, float* normalizer, const size_t size_rhs, const bool no_norm, const size_t numel_x, 
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const int32_t dec) {


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
		mexPrintf("size_center_x = %u\n", size_center_x);
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
		mexPrintf("Created OpenCL kernel\n");
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

	size_t outSize;
	if (size_rhs == Nx * Ny * Nz)
		outSize = pituus[0];
	else
		outSize = Nx * Ny * Nz;
	const size_t outSize2 = 1;
	const size_t outSize3 = Nx * Ny * Nz;

	f_b_project(num_devices_context, kerroin, cpu_device, context, commandQueues, koko, lor1, z_det, x, y, rhs, sc_ra, Nx, Ny, Nz, dx, dy, dz, bx, by,
		bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, normalization, atten, size_atten, norm, size_norm,
		pseudos, det_per_ring, prows, L, raw, size_z, im_dim, kernel_sum, kernel, output, normalizer, size_rhs, no_norm, numel_x, tube_width, crystal_size_z,
		x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, precompute, dec, projector_type);


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


// Main reconstruction function
void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, 
	const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par, 
	const char* k_path, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const char* fileName, const uint32_t device, float kerroin, const size_t numel_x, const char* header_directory) {


	const uint32_t im_dim = Nx * Ny * Nz;
	bool atomic_64bit = false;

	cl_int status = CL_SUCCESS;
	cl_context context = NULL;
	cl_uint num_devices_context;
	cl_kernel kernel;

	cl_device_id * devices = (cl_device_id*)alloca(sizeof(cl_device_id));

	status = clGetPlatformsContextSingle(device, context, num_devices_context, devices);

	if (status != CL_SUCCESS) {
		status = clReleaseContext(context);
		mexPrintf("Error while getting platforms\n");
		return;
	}

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
