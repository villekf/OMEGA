#include "mex.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <CL/cl.hpp>

using std::vector;

// this error code was taken from: https://stackoverflow.com/questions/24326432/convenient-way-to-show-opencl-error-codes
const char *getErrorString(cl_int error)
{
	switch (error) {
		// run-time and JIT compiler errors
	case 0: return "CL_SUCCESS";
	case -1: return "CL_DEVICE_NOT_FOUND";
	case -2: return "CL_DEVICE_NOT_AVAILABLE";
	case -3: return "CL_COMPILER_NOT_AVAILABLE";
	case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case -5: return "CL_OUT_OF_RESOURCES";
	case -6: return "CL_OUT_OF_HOST_MEMORY";
	case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
	case -8: return "CL_MEM_COPY_OVERLAP";
	case -9: return "CL_IMAGE_FORMAT_MISMATCH";
	case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case -11: return "CL_BUILD_PROGRAM_FAILURE";
	case -12: return "CL_MAP_FAILURE";
	case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
	case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
	case -15: return "CL_COMPILE_PROGRAM_FAILURE";
	case -16: return "CL_LINKER_NOT_AVAILABLE";
	case -17: return "CL_LINK_PROGRAM_FAILURE";
	case -18: return "CL_DEVICE_PARTITION_FAILED";
	case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

		// compile-time errors
	case -30: return "CL_INVALID_VALUE";
	case -31: return "CL_INVALID_DEVICE_TYPE";
	case -32: return "CL_INVALID_PLATFORM";
	case -33: return "CL_INVALID_DEVICE";
	case -34: return "CL_INVALID_CONTEXT";
	case -35: return "CL_INVALID_QUEUE_PROPERTIES";
	case -36: return "CL_INVALID_COMMAND_QUEUE";
	case -37: return "CL_INVALID_HOST_PTR";
	case -38: return "CL_INVALID_MEM_OBJECT";
	case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case -40: return "CL_INVALID_IMAGE_SIZE";
	case -41: return "CL_INVALID_SAMPLER";
	case -42: return "CL_INVALID_BINARY";
	case -43: return "CL_INVALID_BUILD_OPTIONS";
	case -44: return "CL_INVALID_PROGRAM";
	case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
	case -46: return "CL_INVALID_KERNEL_NAME";
	case -47: return "CL_INVALID_KERNEL_DEFINITION";
	case -48: return "CL_INVALID_KERNEL";
	case -49: return "CL_INVALID_ARG_INDEX";
	case -50: return "CL_INVALID_ARG_VALUE";
	case -51: return "CL_INVALID_ARG_SIZE";
	case -52: return "CL_INVALID_KERNEL_ARGS";
	case -53: return "CL_INVALID_WORK_DIMENSION";
	case -54: return "CL_INVALID_WORK_GROUP_SIZE";
	case -55: return "CL_INVALID_WORK_ITEM_SIZE";
	case -56: return "CL_INVALID_GLOBAL_OFFSET";
	case -57: return "CL_INVALID_EVENT_WAIT_LIST";
	case -58: return "CL_INVALID_EVENT";
	case -59: return "CL_INVALID_OPERATION";
	case -60: return "CL_INVALID_GL_OBJECT";
	case -61: return "CL_INVALID_BUFFER_SIZE";
	case -62: return "CL_INVALID_MIP_LEVEL";
	case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
	case -64: return "CL_INVALID_PROPERTY";
	case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
	case -66: return "CL_INVALID_COMPILER_OPTIONS";
	case -67: return "CL_INVALID_LINKER_OPTIONS";
	case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";
	case -69: return "CL_INVALID_PIPE_SIZE";
	case -70: return "CL_INVALID_DEVICE_QUEUE";
	case -71: return "CL_INVALID_SPEC_ID";
	case -72: return "CL_MAX_SIZE_RESTRICTION_EXCEEDED";

		// extension errors
	case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
	case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
	case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
	case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
	case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
	case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
	case -1006: return "CL_INVALID_D3D11_DEVICE_KHR";
	case -1007: return "CL_INVALID_D3D11_RESOURCE_KHR";
	case -1008: return "CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR";
	case -1009: return "CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR";
	case -1010: return "CL_INVALID_DX9_MEDIA_ADAPTER_KHR";
	case -1011: return "CL_INVALID_DX9_MEDIA_SURFACE_KHR";
	case -1012: return "CL_DX9_MEDIA_SURFACE_ALREADY_ACQUIRED_KHR";
	case -1013: return "CL_DX9_MEDIA_SURFACE_NOT_ACQUIRED_KHR";
	case -1057: return "CL_DEVICE_PARTITION_FAILED_EXT";
	case -1058: return "CL_INVALID_PARTITION_COUNT_EXT";
	case -1059: return "CL_INVALID_PARTITION_NAME_EXT";
	case -1093: return "CL_INVALID_EGL_OBJECT_KHR";
	case -1092: return "CL_EGL_RESOURCE_NOT_ACQUIRED_KHR";
	case -1094: return "CL_INVALID_ACCELERATOR_INTEL";
	case -1095: return "CL_INVALID_ACCELERATOR_TYPE_INTEL";
	case -1096: return "CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL";
	case -1097: return "CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL";
	case -1098: return "CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL";
	case -1099: return "CL_INVALID_VA_API_MEDIA_SURFACE_INTEL";
	case -1100: return "CL_VA_API_MEDIA_SURFACE_ALREADY_ACQUIRED_INTEL";
	case -1101: return "CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL";
	default: return "Unknown OpenCL error";
	}
}

void reconstruction(const size_t koko, const uint16_t* lor1, const double* z_det, const double* x, const double* y, const int Nx, const int Ny, 
	const int Nz, const double d, const double dz, const double bxf, const double byf, const double bzf, const double bxb, 
	const double byb, const double bzb, const double bx, const double by, const double bz, const double maxxx, const double maxyy, const double minxx, const double minyy, 
	const double zmax, const double xxv, const double yyv, const double xxa, const double yya, const double NSlices, const int pituus, const uint64_t summa, 
	const size_t koko_l, const uint32_t* xy_index, const uint32_t* z_index, const int size_x, const int TotSinos, double* ele, int* indices,
	const bool verbose, const int attenuation_correction, const double* atten, const int size_atten, const char* k_path, const uint32_t* lor2, 
	const size_t koko2, const int device) {

	//double tottime;

	//try {

	cl_int status;
	cl_uint numPlatforms = 0;
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	// Allocate enough space for each platform
	cl_platform_id *platforms = NULL;
	platforms = (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
	// Fill in the platforms
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	// Retrieve the number of devices
	cl_uint numDevices = 0;
	status = clGetDeviceIDs(platforms[1], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	// Allocate enough space for each device
	cl_device_id *devices;
	devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
	// Fill in the devices
	status = clGetDeviceIDs(platforms[1], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);

	cl_context af_context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);

	cl_device_id af_device_id = devices[0];

	// Only create a command-queue for the first device
	cl_command_queue af_queue = clCreateCommandQueue(af_context, af_device_id, 0, &status);


	status = CL_SUCCESS;

	std::string kernel_path;

	//cl_program program = clCreateProgramWithSource(af_context, 1, (const char **)&programSource, NULL, &status);

	kernel_path = k_path;
	std::fstream sourceFile(kernel_path.c_str());
	std::string content((std::istreambuf_iterator<char>(sourceFile)),std::istreambuf_iterator<char>());
	const char* sourceCode = new char[content.size()];
	sourceCode = content.c_str();
	cl_program program = clCreateProgramWithSource(af_context, 1, (const char **)&sourceCode, NULL, &status);

	//mexPrintf("path: %s\n", kernel_path.c_str());
	//mexEvalString("pause(.0001);");

	//mexPrintf(sourceCode);
	//mexPrintf("%d\n", content.size());

	//std::ifstream sourceFile("siddon_kernel.cl");
	//std::string sourceCode(std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));
	//cl::Program::Sources programSource(1, std::make_pair(sourceCode.c_str(), sourceCode.length()));

	//cl_program program = clCreateProgramWithSource(af_context, 1, (const char **)&programSource, NULL, &status);

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to create OpenCL program\n");
		return;
	}
	else if (verbose) {
		mexPrintf("OpenCL program successfully created\n");
		mexEvalString("pause(.0001);");
	}


	status = clBuildProgram(program, 1, &af_device_id, NULL, NULL, NULL);

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to build OpenCL program. Build log: \n");
		size_t len;
		char *buffer;
		clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		buffer = (char*)calloc(len, sizeof(size_t));
		clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
		mexPrintf("%s\n", buffer);
		return;
	}
	else if (verbose) {
		mexPrintf("OpenCL program successfully built\n");
		mexEvalString("pause(.0001);");
	}

	cl_kernel kernel = clCreateKernel(program, "siddon", &status);

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexErrMsgTxt("Failed to create OpenCL kernel\n");
		return;
	}
	else if (verbose) {
		mexPrintf("OpenCL kernel successfully created\n");
		mexEvalString("pause(.0001);");
	}

	cl_mem d_Nx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_Ny = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_Nz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_dz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_d = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_by = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bxf = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_byf = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bzf = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bxb = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_byb = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_bzb = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_maxxx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_maxyy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_minxx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_minyy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_zmax = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_NSlices = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_xxv = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_yyv = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_xxa = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_yya = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double), NULL, &status);
	cl_mem d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double) * TotSinos * 2, NULL, &status);
	cl_mem d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double) * size_x * 2, NULL, &status);
	cl_mem d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double) * size_x * 2, NULL, &status);
	cl_mem d_atten = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double) * size_atten, NULL, &status);
	cl_mem d_size_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_TotSinos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_attenuation_correction = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_xyindex = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int) * koko, NULL, &status);
	cl_mem d_zindex = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int) * koko, NULL, &status);
	cl_mem d_row = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int) * koko2, NULL, &status);
	cl_mem d_lor1 = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko_l, NULL, &status);

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		status = clReleaseProgram(program);
		status = clReleaseKernel(kernel);
		mexErrMsgTxt("Buffer creation failed\n");
	}
	else if (verbose){
		mexPrintf("Buffers created\n");
		mexEvalString("pause(.0001);");
	}
	//cl_mem d_z_det = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double)*mxGetNumberOfElements(prhs[9]), NULL, &status);
	//cl_mem d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double)*mxGetNumberOfElements(prhs[10]), NULL, &status);
	//cl_mem d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(double)*mxGetNumberOfElements(prhs[11]), NULL, &status);

	status = clEnqueueWriteBuffer(af_queue, d_Nx, CL_FALSE, 0, sizeof(int), &Nx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_Ny, CL_FALSE, 0, sizeof(int), &Ny, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_Nz, CL_FALSE, 0, sizeof(int), &Nz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_dz, CL_FALSE, 0, sizeof(double), &dz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_d, CL_FALSE, 0, sizeof(double), &d, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bx, CL_FALSE, 0, sizeof(double), &bx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_by, CL_FALSE, 0, sizeof(double), &by, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bz, CL_FALSE, 0, sizeof(double), &bz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bxf, CL_FALSE, 0, sizeof(double), &bxf, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_byf, CL_FALSE, 0, sizeof(double), &byf, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bzf, CL_FALSE, 0, sizeof(double), &bzf, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bxb, CL_FALSE, 0, sizeof(double), &bxb, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_byb, CL_FALSE, 0, sizeof(double), &byb, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bzb, CL_FALSE, 0, sizeof(double), &bzb, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_maxxx, CL_FALSE, 0, sizeof(double), &maxxx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_maxyy, CL_FALSE, 0, sizeof(double), &maxyy, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_minxx, CL_FALSE, 0, sizeof(double), &minxx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_minyy, CL_FALSE, 0, sizeof(double), &minyy, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_zmax, CL_FALSE, 0, sizeof(double), &zmax, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_NSlices, CL_FALSE, 0, sizeof(double), &NSlices, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_xxv, CL_FALSE, 0, sizeof(double), &xxv, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_yyv, CL_FALSE, 0, sizeof(double), &yyv, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_xxa, CL_FALSE, 0, sizeof(double), &xxa, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_yya, CL_FALSE, 0, sizeof(double), &yya, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(double) * size_x * 2, x, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(double) * size_x * 2, y, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(double) * TotSinos * 2, z_det, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_atten, CL_FALSE, 0, sizeof(double) * size_atten, atten, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_size_x, CL_FALSE, 0, sizeof(int), &size_x, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_TotSinos, CL_FALSE, 0, sizeof(int), &TotSinos, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_attenuation_correction, CL_FALSE, 0, sizeof(int), &attenuation_correction, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_xyindex, CL_FALSE, 0, sizeof(int) * koko, xy_index, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_zindex, CL_FALSE, 0, sizeof(int) * koko, z_index, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_row, CL_FALSE, 0, sizeof(int) * koko2, lor2, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_lor1, CL_FALSE, 0, sizeof(uint16_t) * koko_l, lor1, 0, NULL, NULL);

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		status = clReleaseProgram(program);
		status = clReleaseKernel(kernel);
		mexErrMsgTxt("Buffer write failed\n");
	}
	else if (verbose) {
		mexPrintf("Buffers written\n");
		mexEvalString("pause(.0001);");
	}

	//int osa_iter = 0;

	//pz_osem_apu = pz_osem(span, 0);

	//int32_t zero = 0;
	//double zerof = 0.f;

	int im_dim = Nx * Ny * Nz;


			size_t length = pituus;

			size_t outputSize1 = sizeof(int32_t) * summa;
			size_t outputSize2 = sizeof(double) * summa;


			cl_mem d_indices = clCreateBuffer(af_context, CL_MEM_READ_WRITE, outputSize1, NULL, &status);

			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				status = clReleaseProgram(program);
				status = clReleaseKernel(kernel);
				mexErrMsgTxt("Failed to create buffer for system matrix column indices. Check that you have sufficient memory available\n");
			}
			else if (verbose) {
				mexPrintf("Buffer for system matrix column indices created\n");
				mexEvalString("pause(.0001);");
			}

			cl_mem d_elements = clCreateBuffer(af_context, CL_MEM_READ_WRITE, outputSize2, NULL, &status);

			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				//status = clReleaseMemObject(d_indices);
				status = clReleaseProgram(program);
				status = clReleaseKernel(kernel);
				mexErrMsgTxt("Failed to create buffer for system matrix elements. Check that you have sufficient memory available\n");
			}
			else if (verbose) {
				mexPrintf("Buffer for system matrix elements created\n");
				mexEvalString("pause(.0001);");
			}

			clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_indices);
			clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_elements);
			clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_lor1);
			clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_Nx);
			clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_Ny);
			clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_Nz);
			clSetKernelArg(kernel, 6, sizeof(cl_mem), &d_dz);
			clSetKernelArg(kernel, 7, sizeof(cl_mem), &d_d);
			clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_bz);
			clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_bx);
			clSetKernelArg(kernel, 10, sizeof(cl_mem), &d_by);
			clSetKernelArg(kernel, 11, sizeof(cl_mem), &d_bxf);
			clSetKernelArg(kernel, 12, sizeof(cl_mem), &d_byf);
			clSetKernelArg(kernel, 13, sizeof(cl_mem), &d_bzf);
			clSetKernelArg(kernel, 14, sizeof(cl_mem), &d_bxb);
			clSetKernelArg(kernel, 15, sizeof(cl_mem), &d_byb);
			clSetKernelArg(kernel, 16, sizeof(cl_mem), &d_bzb);
			clSetKernelArg(kernel, 17, sizeof(cl_mem), &d_maxxx);
			clSetKernelArg(kernel, 18, sizeof(cl_mem), &d_maxyy);
			clSetKernelArg(kernel, 19, sizeof(cl_mem), &d_minxx);
			clSetKernelArg(kernel, 20, sizeof(cl_mem), &d_minyy);
			clSetKernelArg(kernel, 21, sizeof(cl_mem), &d_zmax);
			clSetKernelArg(kernel, 22, sizeof(cl_mem), &d_NSlices);
			clSetKernelArg(kernel, 23, sizeof(cl_mem), &d_x);
			clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_y);
			clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_z);
			clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_size_x);
			clSetKernelArg(kernel, 27, sizeof(cl_mem), &d_row);
			clSetKernelArg(kernel, 28, sizeof(cl_mem), &d_xxv);
			clSetKernelArg(kernel, 29, sizeof(cl_mem), &d_yyv);
			clSetKernelArg(kernel, 30, sizeof(cl_mem), &d_xxa);
			clSetKernelArg(kernel, 31, sizeof(cl_mem), &d_yya);
			clSetKernelArg(kernel, 32, sizeof(cl_mem), &d_xyindex);
			clSetKernelArg(kernel, 33, sizeof(cl_mem), &d_zindex);
			clSetKernelArg(kernel, 34, sizeof(cl_mem), &d_TotSinos);
			clSetKernelArg(kernel, 35, sizeof(cl_mem), &d_attenuation_correction);
			clSetKernelArg(kernel, 36, sizeof(cl_mem), &d_atten);
			////clFinish(af_queue);
			status = clEnqueueNDRangeKernel(af_queue, kernel, 1, NULL, &length, NULL, 0, NULL, NULL);

			////clFinish(af_queue);
			////sync();

			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				clReleaseCommandQueue(af_queue);
				clReleaseProgram(program);
				clReleaseKernel(kernel);
				clReleaseMemObject(d_indices);
				clReleaseMemObject(d_elements);
				clReleaseContext(af_context);
				mexErrMsgTxt("Failed to execute the OpenCL kernel\n");
			}
			else if (verbose) {
				mexPrintf("OpenCL kernel executed successfully\n");
				mexEvalString("pause(.0001);");
			}

			//sync();
			//if (osa_iter == 2)
			//	(s_indices).host(ele);

			////clFinish(af_queue);
			//if (osa_iter == 0) {
			//	status = clEnqueueReadBuffer(af_queue, d_indices, CL_TRUE, 0, outputSize1, ele, 0, NULL, NULL);
			//}

			//if (status != CL_SUCCESS) {
			//	std::cerr << getErrorString(status) << std::endl;
			//	mexPrintf("Failed to copy the output indices to host\n");
			//}
			//else if (verbose)
			//	mexPrintf("Output indices successfully copied to host\n");

			//status = clEnqueueReadBuffer(af_queue, d_elements, CL_TRUE, 0, outputSize, ele, 0, NULL, NULL);

			//if (status != CL_SUCCESS) {
			//	std::cerr << getErrorString(status) << std::endl;
			//	mexPrintf("Failed to copy the output elements to host\n");
			//}
			//else if (verbose)
			//	mexPrintf("Output elements successfully copied to host\n");

			clFinish(af_queue);

			status = clEnqueueReadBuffer(af_queue, d_indices, CL_TRUE, 0, outputSize1, indices, 0, NULL, NULL);
			status = clEnqueueReadBuffer(af_queue, d_elements, CL_TRUE, 0, outputSize2, ele, 0, NULL, NULL);




			//if (verbose) {
			//	mexPrintf("Memory released\n");
			//	mexEvalString("pause(.0001);");
			//}



	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexErrMsgTxt("Failed to release program\n");
	}
	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexErrMsgTxt("Failed to release kernel\n");
	}
	clReleaseCommandQueue(af_queue);
	status = clReleaseMemObject(d_indices);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexErrMsgTxt("Failed to release the memory of column indices\n");
	}
	status = clReleaseMemObject(d_elements);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexErrMsgTxt("Failed to release the memory of elements\n");
	}
	clReleaseMemObject(d_Nx);
	clReleaseMemObject(d_Ny);
	clReleaseMemObject(d_Nz);
	clReleaseMemObject(d_dz);
	clReleaseMemObject(d_d);
	clReleaseMemObject(d_bz);
	clReleaseMemObject(d_bx);
	clReleaseMemObject(d_by);
	clReleaseMemObject(d_bxb);
	clReleaseMemObject(d_byb);
	clReleaseMemObject(d_bzb);
	clReleaseMemObject(d_bxf);
	clReleaseMemObject(d_byf);
	clReleaseMemObject(d_bzf);
	clReleaseMemObject(d_maxxx);
	clReleaseMemObject(d_minxx);
	clReleaseMemObject(d_maxyy);
	clReleaseMemObject(d_minyy);
	clReleaseMemObject(d_zmax);
	clReleaseMemObject(d_NSlices);
	clReleaseMemObject(d_xxv);
	clReleaseMemObject(d_yyv);
	clReleaseMemObject(d_xxa);
	clReleaseMemObject(d_yya);
	clReleaseMemObject(d_z);
	clReleaseMemObject(d_x);
	clReleaseMemObject(d_y);
	clReleaseMemObject(d_atten);
	clReleaseMemObject(d_size_x);
	clReleaseMemObject(d_TotSinos);
	clReleaseMemObject(d_attenuation_correction);
	clReleaseMemObject(d_xyindex);
	clReleaseMemObject(d_zindex);
	clReleaseMemObject(d_row);
	clReleaseMemObject(d_lor1);
	clReleaseContext(af_context);

	//int32_t* output = (int32_t*)std::malloc(outputSize);
	//double* output = (double*)std::malloc(outputSize);


	//for (int kk = 0; kk < summa[0]; kk++) {
	//	if (output[kk] > 0)
	//		mexPrintf("output %d\n", output[kk]);
	//}

	//free(sourceCode);

	//}
	//return indices;
	//return output;
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	if (nrhs != 32)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 32.");

	//if (nlhs != 1)
	//	mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	//if (nlhs != 2)
	//	mexErrMsgTxt("Invalid number of output arguments.  There must be exactly two.");


	//int pixelsize = (int)mxGetScalar(prhs[0]);

	int Ny = (int)mxGetScalar(prhs[1]);

	int Nx = (int)mxGetScalar(prhs[2]);

	int Nz = (int)mxGetScalar(prhs[3]);

	double d = (double)mxGetScalar(prhs[4]);

	double dz = (double)mxGetScalar(prhs[5]);

	double by = (double)mxGetScalar(prhs[6]);

	double bx = (double)mxGetScalar(prhs[7]);

	double bz = (double)mxGetScalar(prhs[8]);

	double *z_det = (double*)mxGetData(prhs[9]);

	double *x = (double*)mxGetData(prhs[10]);

	double *y = (double*)mxGetData(prhs[11]);

	double *iij = (double*)mxGetData(prhs[12]);

	double *jji = (double*)mxGetData(prhs[13]);

	double *kkj = (double*)mxGetData(prhs[14]);

	double *yy = (double*)mxGetData(prhs[15]);

	double *xx = (double*)mxGetData(prhs[16]);

	int NSinos = (int)mxGetScalar(prhs[17]);

	double NSlices = (double)mxGetScalar(prhs[18]);

	//vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[12]));

	//vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[13]));

	//vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[14]));

	//vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[15]) - 1);

	//vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[16]) - 1);

	double yyv = yy[1] - yy[0];

	double xxv = xx[1] - xx[0];

	double yya = yy[0];

	double xxa = xx[0];

	int size_x = (int)mxGetScalar(prhs[19]);

	double zmax = (double)mxGetScalar(prhs[20]);

	int TotSinos = (int)mxGetScalar(prhs[21]);

	double *atten = (double*)mxGetData(prhs[22]);
	int size_atten = mxGetNumberOfElements(prhs[22]);

	//vector<int> indeksit(index, index + mxGetNumberOfElements(prhs[24]) - 1);

	uint32_t pituus = (uint32_t)mxGetScalar(prhs[23]);

	int attenuation_correction = (int)mxGetScalar(prhs[24]);

	uint16_t* lor1 = (uint16_t*)mxGetData(prhs[25]);
	size_t koko_l = mxGetNumberOfElements(prhs[25]);

	uint64_t summa = (uint64_t)mxGetScalar(prhs[26]);

	uint32_t* xy_index = (uint32_t*)mxGetData(prhs[27]);
	size_t koko = mxGetNumberOfElements(prhs[27]);

	uint32_t* z_index = (uint32_t*)mxGetData(prhs[28]);

	uint32_t* lor2 = (uint32_t*)mxGetData(prhs[29]);
	size_t koko2 = mxGetNumberOfElements(prhs[29]);

	bool verbose = (bool)mxGetScalar(prhs[30]);

	int device = (int)mxGetScalar(prhs[31]);

	//int device = 2;
	//af::setDevice(1);
	//af::setDevice(device);

	//mexPrintf("size_atten %d\n", size_atten);
	//mexEvalString("pause(.0001);");

	//int loop_var_par = 1;

	//vector<int> lor;

	//vector<int> indices;

	//vector<double> elements;

	//clock_t time = clock();

	//vector<vector<int32_t>> rows;

	//rows.assign(subsets, vector<int32_t>(ind_size, 0));

	//for (int ii = 0; ii < subsets; ii++) {
	//	int oo = 1;
	//	for (int kk = pituus[ii]; kk < pituus[ii + 1]; kk++) {
	//		int32_t apu = static_cast<int32_t>(lor1[kk]);
	//		rows[ii][oo] = (oo > 1) ? apu + rows[ii][oo - 1] : apu;
	//		oo++;
	//	}
	//}

	//time = clock() - time;

	//mexPrintf("Function elapsed time is %f seconds\n", ((double)time) / CLOCKS_PER_SEC);

	//char *argv = (char *)mxCalloc(1, sizeof(char *));
	char *k_path = mxArrayToString(prhs[0]);

	double maxyy = *std::max_element(yy, yy + Ny + 1);
	double minyy = *std::min_element(yy, yy + Ny + 1);

	double maxxx = *std::max_element(xx, xx + Nx + 1);
	double minxx = *std::min_element(xx, xx + Nx + 1);

	double bxf = bx + iij[0] * d;
	double byf = by + jji[0] * d;
	double bzf = bz + kkj[0] * dz;
	double bxb = bx + iij[mxGetNumberOfElements(prhs[12]) - 1] * d;
	double byb = by + jji[mxGetNumberOfElements(prhs[13]) - 1] * d;
	double bzb = bz + kkj[mxGetNumberOfElements(prhs[14]) - 1] * dz;


	//mexPrintf("pituus %d\n", pituus);
	//mexPrintf("koko_l %f\n", iij[mxGetNumberOfElements(prhs[12]) - 1]);

	//size_t outSize = summa[0];
	//size_t outSize2 = 1;
	//size_t outSize = pituus[3] - pituus[2];
	//size_t outSize2 = 1;
	mwSize N = Nx * Ny * Nz;
	mwSize nzmax = summa;
	mwSize rows = pituus;

	plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);
	double* elements = (double*)mxGetData(plhs[0]);
	mwIndex* indices2 = mxGetIr(plhs[0]);
	mwIndex* lor = mxGetJc(plhs[0]);

	int* indices = new int[summa];

	//plhs[1] = mxCreateNumericMatrix(summa, 1, mxDOUBLE_CLASS, mxREAL);

	//double* elements = (double*)mxGetData(plhs[1]);

	//plhs[1] = mxCreateNumericMatrix(summa, 1, mxINT32_CLASS, mxREAL);

	//int* indices = (int*)mxGetData(plhs[1]);

	reconstruction(koko, lor1, z_det, x, y, Nx, Ny, Nz, d, dz, bxf, byf, bzf, bxb, byb, bzb, bx, by, bz, maxxx, maxyy, minxx, minyy, zmax, xxv, 
		yyv, xxa, yya, NSlices, pituus, summa, koko_l, xy_index, z_index, size_x, TotSinos, elements, indices, verbose, attenuation_correction, atten, size_atten,
		k_path, lor2, koko2, device);

	for (int kk = 0; kk < summa; kk++)
		indices2[kk] = indices[kk];

	for (int kk = 0; kk <= pituus; kk++)
		lor[kk] = lor2[kk];

	//mxSetJc(plhs[0], (mwIndex*)mxGetData(prhs[29]));

	delete[] indices;

	return;
}