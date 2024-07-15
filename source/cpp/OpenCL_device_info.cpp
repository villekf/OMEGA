#include "mex.h"
#include "opencl_error.hpp"
#include <string>
#include <iostream>


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	cl_int status;
	cl_uint num_platforms;
	size_t size;
	cl_context context = NULL;

	std::string output_string;

	status = clGetPlatformIDs(0, NULL, &num_platforms);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	output_string = "No. platforms: " + std::to_string(num_platforms) + "\n";
	//mexPrintf("No. platforms: %u\n", num_platforms);

	cl_platform_id platforms[64];

	status = clGetPlatformIDs(num_platforms, platforms, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	cl_uint num_devices;
	for (int kk = 0; kk < num_platforms; kk++) {
		status = clGetPlatformInfo(platforms[kk], CL_PLATFORM_NAME, 0, NULL, &size);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		char * name = (char *)alloca(sizeof(char) * size);
		status = clGetPlatformInfo(platforms[kk], CL_PLATFORM_NAME, size, name, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		//mexPrintf("\nPlatform %d: %s\n", kk, name);
		output_string += "\nPlatform " + std::to_string(kk) + ": " + name + "\n";
		status = clGetDeviceIDs(platforms[kk], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		output_string += "No. devices:  " + std::to_string(num_devices) + "\n";
		//mexPrintf("No. devices: %u\n", num_devices);
		cl_device_id *devices = new cl_device_id[num_devices];
		status = clGetDeviceIDs(platforms[kk], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);
		if (status != CL_SUCCESS) {
			delete[] devices;
			getErrorString(status);
			return;
		}

		for (int ii = 0; ii < num_devices; ii++) {
			cl_ulong mem;
			size_t size2;
			status = clGetDeviceInfo(devices[ii], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
			if (status != CL_SUCCESS) {
				delete[] devices;
				getErrorString(status);
				return;
			}
			status = clGetDeviceInfo(devices[ii], CL_DEVICE_NAME, 0, NULL, &size2);
			if (status != CL_SUCCESS) {
				delete[] devices;
				getErrorString(status);
				return;
			}
			char * name2 = (char *)alloca(sizeof(char) * size2);
			status = clGetDeviceInfo(devices[ii], CL_DEVICE_NAME, size2, name2, NULL);
			if (status != CL_SUCCESS) {
				delete[] devices;
				getErrorString(status);
				return;
			}
			mem = mem / (1024 * 1024);
			//mexPrintf("Device %d: %s with %u MB of memory\n", ii, name2, mem);
			output_string += "Device  " + std::to_string(ii) + ": " + name2 + " with " + std::to_string(mem) + " MB of memory\n";
		}
		delete[] devices;
	}
	const char* c = output_string.c_str();

	//mexPrintf("\n%s\n", c);

	plhs[0] = mxCreateString(c);
}