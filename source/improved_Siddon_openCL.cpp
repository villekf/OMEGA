/**************************************************************************
* Computes OSEM in OpenCL for OMEGA.
* This is the matrix version of the OSEM. First the improved Siddon is used
* to calculate the system matrix. ArrayFire sparse matrix is formed from 
* the output row and column indices (CSR format) and values. ArrayFire 
* functions are then used to compute the OSEM estimate.
* Output is the OSEM estimate.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
* This function is also most likely slower than the matrix free version. 
* This function was chronologically first of the OpenCL version though and 
* thus remains as an alternative method. This function, however, won't
* receive any updates besides bug fixes.
* Currently this code is slow on (at least) CUDA devices.
* This code has been deprecated and might be removed in a future release.
* 
* Copyright (C) 2019 Ville-Veikko Wettenhovi
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
#include "functions.hpp"
#include "AF_opencl_functions.hpp"

using namespace af;


void reconstruction(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* Sin, const int Nx, const int Ny, 
	const int Nz, const int Niter, const float* x0, const float dx, const float dy, const float dz, const float bzb, 
	const float bx, const float by, const float bz, const float maxxx, const float maxyy, const float zmax, const float NSlices, const int* pituus, 
	const uint64_t* summa, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const int size_x, const int TotSinos, float* ele,
	const bool verbose, const int attenuation_correction, const float* atten, const int size_atten, const int subsets, const float epps, const char* k_path, 
	const int* pseudos, const int det_per_ring, const int prows, const uint16_t* L, const uint8_t raw, const size_t size_z) {

	const array lor(koko_l, lor1, afHost);

	array zindex, xindex, LL;

	if (raw) {
		LL = array(koko, L, afHost);
		zindex = constant(0, 1, 1, u16);
		xindex = constant(0, 1, 1, u32);
	}
	else {
		zindex = array(koko, z_index, afHost);
		xindex = array(koko, xy_index, afHost);
		LL = constant(0, 1, 1, u16);
	}

	const array Sino(koko_l, Sin, afHost);

	array pz_osem = constant(1.f, Ny*Nx*Nz, Niter + 1);

	const array x00(Nx*Ny*Nz, x0, afHost);

	pz_osem(span, 0) = x00;

	array pz_osem_apu;

	static cl_context af_context = afcl::getContext();
	static cl_device_id af_device_id = afcl::getDeviceId();
	static cl_command_queue af_queue = afcl::getQueue();

	int status = CL_SUCCESS;

	std::string kernel_path;

	kernel_path = k_path;
	std::fstream sourceFile(kernel_path.c_str());
	std::string content((std::istreambuf_iterator<char>(sourceFile)),std::istreambuf_iterator<char>());
	const char* sourceCode = new char[content.size()];
	sourceCode = content.c_str();
	cl_program program = clCreateProgramWithSource(af_context, 1, (const char **)&sourceCode, NULL, &status);

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to create OpenCL program\n");
		return;
	}
	else if (verbose) {
		mexPrintf("OpenCL program successfully created\n");
		mexEvalString("pause(.0001);");
	}


	status = clBuildProgram(program, 1, &af_device_id, NULL, NULL, NULL);

	if (status != CL_SUCCESS) {
		getErrorString(status);
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
		getErrorString(status);
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}
	else if (verbose) {
		mexPrintf("OpenCL kernel successfully created\n");
		mexEvalString("pause(.0001);");
	}

	cl_mem d_x, d_y, d_z;

	cl_mem d_Nx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_Ny = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_Nz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_dz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_dx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_dy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_bz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_bx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_by = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_bzb = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_maxxx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_maxyy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_zmax = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	cl_mem d_NSlices = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	if (raw) {
		d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
		d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x, NULL, &status);
		d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x, NULL, &status);
	}
	else {
		d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * TotSinos * 2, NULL, &status);
		d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x * 2, NULL, &status);
		d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x * 2, NULL, &status);
	}
	cl_mem d_atten = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
	cl_mem d_size_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_TotSinos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_attenuation_correction = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_det_per_ring = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	cl_mem d_pseudos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int) * prows, NULL, &status);
	cl_mem d_raw = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &status);

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer creation failed\n");
	}
	else if (verbose){
		mexPrintf("Buffers created\n");
		mexEvalString("pause(.0001);");
	}

	status = clEnqueueWriteBuffer(af_queue, d_Nx, CL_FALSE, 0, sizeof(int), &Nx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_Ny, CL_FALSE, 0, sizeof(int), &Ny, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_Nz, CL_FALSE, 0, sizeof(int), &Nz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_dz, CL_FALSE, 0, sizeof(float), &dz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_dx, CL_FALSE, 0, sizeof(float), &dx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_dy, CL_FALSE, 0, sizeof(float), &dy, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bx, CL_FALSE, 0, sizeof(float), &bx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_by, CL_FALSE, 0, sizeof(float), &by, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bz, CL_FALSE, 0, sizeof(float), &bz, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_bzb, CL_FALSE, 0, sizeof(float), &bzb, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_maxxx, CL_FALSE, 0, sizeof(float), &maxxx, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_maxyy, CL_FALSE, 0, sizeof(float), &maxyy, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_zmax, CL_FALSE, 0, sizeof(float), &zmax, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_NSlices, CL_FALSE, 0, sizeof(float), &NSlices, 0, NULL, NULL);
	if (raw) {
		status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(float) * size_x, x, 0, NULL, NULL);
		status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(float) * size_x, y, 0, NULL, NULL);
		status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
	}
	else {
		status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(float) * size_x * 2, x, 0, NULL, NULL);
		status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(float) * size_x * 2, y, 0, NULL, NULL);
		status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(float) * TotSinos * 2, z_det, 0, NULL, NULL);
	}
	status = clEnqueueWriteBuffer(af_queue, d_atten, CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_size_x, CL_FALSE, 0, sizeof(int), &size_x, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_TotSinos, CL_FALSE, 0, sizeof(int), &TotSinos, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_attenuation_correction, CL_FALSE, 0, sizeof(int), &attenuation_correction, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_det_per_ring, CL_FALSE, 0, sizeof(int), &det_per_ring, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_pseudos, CL_FALSE, 0, sizeof(int) * prows, pseudos, 0, NULL, NULL);
	status = clEnqueueWriteBuffer(af_queue, d_raw, CL_FALSE, 0, sizeof(uint8_t), &raw, 0, NULL, NULL);

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer write failed\n");
	}
	else if (verbose) {
		mexPrintf("Buffers written\n");
		mexEvalString("pause(.0001);");
	}

	const int im_dim = Nx * Ny * Nz;

	for (int iter = 0; iter < Niter; iter++) {

		pz_osem_apu = pz_osem(span, iter);

		for (int osa_iter = 0; osa_iter < subsets; osa_iter++) {

			//clFinish(af_queue);

			//if (verbose) {
			//	mexPrintf("Next iteration started\n");
			//	mexEvalString("pause(.0001);");
			//}

			size_t length = pituus[osa_iter + 1] - pituus[osa_iter];

			size_t outputSize1 = sizeof(int) * summa[osa_iter];
			size_t outputSize2 = sizeof(float) * summa[osa_iter];

			array sub_index_array = range(dim4(length), 0, s32) + pituus[osa_iter];

			array apu_lor = lor(sub_index_array);

			array row_csr = join(0, constant(0, 1, s32), accum(apu_lor.as(s32)));

			//row_csr = join(0, constant(0, 1, s32), row_csr);

			array Lo, xo, zo;

			if (raw) {
				Lo = LL(range(dim4(length * 2), 0, u32) + 2 * pituus[osa_iter]);
			}
			else {
				xo = xindex(sub_index_array);
				zo = zindex(sub_index_array);
			}

			cl_mem d_indices = clCreateBuffer(af_context, CL_MEM_READ_WRITE, outputSize1, NULL, &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create buffer for system matrix indices. Check that you have sufficient memory available\n");
			}
			//else if (verbose) {
			//	mexPrintf("Buffer for system matrix indices created\n");
			//	mexEvalString("pause(.0001);");
			//}

			cl_mem d_elements = clCreateBuffer(af_context, CL_MEM_READ_WRITE, outputSize2, NULL, &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create buffer for system matrix elements. Check that you have sufficient memory available\n");
			}
			//else if (verbose) {
			//	mexPrintf("Buffer for system matrix elements created\n");
			//	mexEvalString("pause(.0001);");
			//}

			cl_mem * d_lor = apu_lor.device<cl_mem>();
			cl_mem * d_row = row_csr.device<cl_mem>();

			cl_mem *d_xyindex, *d_zindex, *d_L;

			if (raw) {
				d_xyindex = xindex.device<cl_mem>();
				d_zindex = zindex.device<cl_mem>();
				d_L = Lo.device<cl_mem>();
			}
			else {
				d_xyindex = xo.device<cl_mem>();
				d_zindex = zo.device<cl_mem>();
				d_L = LL.device<cl_mem>();
			}

			clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_indices);
			clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_elements);
			clSetKernelArg(kernel, 2, sizeof(cl_mem), d_lor);
			clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_Nx);
			clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_Ny);
			clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_Nz);
			clSetKernelArg(kernel, 6, sizeof(cl_mem), &d_dz);
			clSetKernelArg(kernel, 7, sizeof(cl_mem), &d_dx);
			clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_dy);
			clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_bz);
			clSetKernelArg(kernel, 10, sizeof(cl_mem), &d_bx);
			clSetKernelArg(kernel, 11, sizeof(cl_mem), &d_by);
			clSetKernelArg(kernel, 12, sizeof(cl_mem), &d_bzb);
			clSetKernelArg(kernel, 13, sizeof(cl_mem), &d_maxxx);
			clSetKernelArg(kernel, 14, sizeof(cl_mem), &d_maxyy);
			clSetKernelArg(kernel, 15, sizeof(cl_mem), &d_zmax);
			clSetKernelArg(kernel, 16, sizeof(cl_mem), &d_NSlices);
			clSetKernelArg(kernel, 17, sizeof(cl_mem), &d_x);
			clSetKernelArg(kernel, 18, sizeof(cl_mem), &d_y);
			clSetKernelArg(kernel, 19, sizeof(cl_mem), &d_z);
			clSetKernelArg(kernel, 20, sizeof(cl_mem), &d_size_x);
			clSetKernelArg(kernel, 21, sizeof(cl_mem), d_row);
			clSetKernelArg(kernel, 22, sizeof(cl_mem), d_xyindex);
			clSetKernelArg(kernel, 23, sizeof(cl_mem), d_zindex);
			clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_TotSinos);
			clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_attenuation_correction);
			clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_atten);
			clSetKernelArg(kernel, 27, sizeof(cl_mem), d_L);
			clSetKernelArg(kernel, 28, sizeof(cl_mem), &d_det_per_ring);
			clSetKernelArg(kernel, 29, sizeof(cl_mem), &d_pseudos);
			clSetKernelArg(kernel, 30, sizeof(int), &prows);
			clSetKernelArg(kernel, 31, sizeof(cl_mem), &d_raw);
			status = clEnqueueNDRangeKernel(af_queue, kernel, 1, NULL, &length, NULL, 0, NULL, NULL);

			////clFinish(af_queue);
			////sync();

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to execute the OpenCL kernel\n");
			}
			//else if (verbose) {
			//	mexPrintf("OpenCL kernel executed successfully\n");
			//	mexEvalString("pause(.0001);");
			//}

			if (raw) {
				zindex.unlock();
				xindex.unlock();
				Lo.unlock();
			}
			else {
				zo.unlock();
				xo.unlock();
				LL.unlock();
			}
			apu_lor.unlock();
			row_csr.unlock();

			////clFinish(af_queue);

			array s_elements = afcl::array(summa[osa_iter], d_elements, f32, true);
			array s_indices = afcl::array(summa[osa_iter], d_indices, s32, true);

			//////sync();

			array H = sparse(row_csr.dims(0) - 1, im_dim, s_elements, row_csr, s_indices);
			array Summ = matmul(H, constant(1.f, H.dims(0), 1, f32), AF_MAT_TRANS);
			pz_osem_apu = (pz_osem_apu / (Summ + epps)) * (matmul(H, (Sino(sub_index_array) / (matmul(H, pz_osem_apu) + epps)), AF_MAT_TRANS));
			//eval(Summ);
			//sync();

			mexPrintf("OSEM sub-iteration %d complete\n", osa_iter + 1);
			mexEvalString("pause(.0001);");

			clFinish(af_queue);

			status = clReleaseMemObject(d_indices);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to release the memory of indices\n");
			}
			status = clReleaseMemObject(d_elements);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to release the memory of elements\n");
			}


			//if (verbose) {
			//	mexPrintf("Memory released\n");
			//	mexEvalString("pause(.0001);");
			//}

		}

		pz_osem(span, iter + 1) = pz_osem_apu;
		mexPrintf("Iteration %d complete\n", iter + 1);
		mexEvalString("pause(.0001);");
	}

	pz_osem.host(ele);

	clReleaseMemObject(d_Nx);
	clReleaseMemObject(d_Ny);
	clReleaseMemObject(d_Nz);
	clReleaseMemObject(d_dz);
	clReleaseMemObject(d_dx);
	clReleaseMemObject(d_dy);
	clReleaseMemObject(d_bz);
	clReleaseMemObject(d_bx);
	clReleaseMemObject(d_by);
	clReleaseMemObject(d_bzb);
	clReleaseMemObject(d_maxxx);
	clReleaseMemObject(d_maxyy);
	clReleaseMemObject(d_zmax);
	clReleaseMemObject(d_NSlices);
	clReleaseMemObject(d_z);
	clReleaseMemObject(d_x);
	clReleaseMemObject(d_y);
	clReleaseMemObject(d_atten);
	clReleaseMemObject(d_size_x);
	clReleaseMemObject(d_TotSinos);
	clReleaseMemObject(d_attenuation_correction);
	clReleaseMemObject(d_raw);
	clReleaseMemObject(d_det_per_ring);
	clReleaseMemObject(d_pseudos);

	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to release program\n");
	}
	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to release kernel\n");
	}

}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	if (nrhs != 39)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 39.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	const int Ny = (int)mxGetScalar(prhs[1]);

	const int Nx = (int)mxGetScalar(prhs[2]);

	const int Nz = (int)mxGetScalar(prhs[3]);

	const float d = (float)mxGetScalar(prhs[4]);

	const float dz = (float)mxGetScalar(prhs[5]);

	const float by = (float)mxGetScalar(prhs[6]);

	const float bx = (float)mxGetScalar(prhs[7]);

	const float bz = (float)mxGetScalar(prhs[8]);

	const float *z_det = (float*)mxGetData(prhs[9]);
	const size_t size_z = mxGetNumberOfElements(prhs[9]);

	const float *x = (float*)mxGetData(prhs[10]);

	const float *y = (float*)mxGetData(prhs[11]);

	const float dy = (float)mxGetScalar(prhs[12]);

	const float maxyy = (float)mxGetScalar(prhs[13]);

	const float maxxx = (float)mxGetScalar(prhs[14]);

	const int NSinos = (int)mxGetScalar(prhs[15]);

	const float NSlices = (float)mxGetScalar(prhs[16]);

	const int size_x = (int)mxGetScalar(prhs[17]);

	const float zmax = (float)mxGetScalar(prhs[18]);

	const int TotSinos = (int)mxGetScalar(prhs[19]);

	const float *atten = (float*)mxGetData(prhs[20]);
	const int size_atten = mxGetNumberOfElements(prhs[20]);

	const int *pituus = (int*)mxGetData(prhs[21]);

	const int attenuation_correction = (int)mxGetScalar(prhs[22]);

	const int Niter = (int)mxGetScalar(prhs[23]);

	const int subsets = (int)mxGetScalar(prhs[24]);

	const bool* rekot = (bool*)mxGetData(prhs[25]);

	const float epps = (float)mxGetScalar(prhs[26]);

	const float* Sin = (float*)mxGetData(prhs[27]);

	const float* x0 = (float*)mxGetData(prhs[28]);

	const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[29]);

	const uint64_t* summa = (uint64_t*)mxGetData(prhs[30]);

	const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[31]);

	const uint16_t* z_index = (uint16_t*)mxGetData(prhs[32]);

	const uint16_t *L = (uint16_t*)mxGetData(prhs[33]);

	const int *pseudos = (int*)mxGetData(prhs[34]);
	const int pRows = (int)mxGetM(prhs[34]);

	const int det_per_ring = (int)mxGetScalar(prhs[35]);

	const uint8_t raw = (uint8_t)mxGetScalar(prhs[36]);

	const bool verbose = (bool)mxGetScalar(prhs[37]);

	const int device = (int)mxGetScalar(prhs[38]);

	af::setDevice(device);

	const char *k_path = mxArrayToString(prhs[0]);

	const float bzb = bz + static_cast<float>(Nz) * dz;

	size_t koko, koko_l;
	koko_l = mxGetNumberOfElements(prhs[27]);
	if (raw) {
		koko = mxGetNumberOfElements(prhs[33]);
	}
	else {
		koko = mxGetNumberOfElements(prhs[31]);
	}

	//size_t outSize = summa[0];
	//size_t outSize2 = 1;
	//size_t outSize = pituus[3] - pituus[2];
	//size_t outSize2 = 1;
	const size_t outSize = Nx * Ny * Nz;
	const size_t outSize2 = Niter + 1;

	plhs[0] = mxCreateNumericMatrix(outSize, outSize2, mxSINGLE_CLASS, mxREAL);

	float* elements = (float*)mxGetData(plhs[0]);

	reconstruction(koko, lor1, z_det, x, y, Sin, Nx, Ny, Nz, Niter, x0, d, dy, dz, bzb, bx, by, bz, maxxx, maxyy, zmax, 
		NSlices, pituus, summa, koko_l, xy_index, z_index, size_x, TotSinos, elements, verbose, attenuation_correction, atten, size_atten, subsets, epps, 
		k_path, pseudos, det_per_ring, pRows, L, raw, size_z);


	sync();
	deviceGC();

	return;
}