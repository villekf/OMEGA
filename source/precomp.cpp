/**************************************************************************
* This function computes the precomputation phase in OpenCL.
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

void precomp_siddon(const cl_uint& num_devices_context, const cl_context& context, const cl_command_queue* commandQueues, uint16_t* lor1, const float* z_det, 
	const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, 
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t size_x, 
	const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, 
	const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim, const cl_kernel& kernel, const size_t numel_x, const size_t local_size) {

	cl_int status = CL_SUCCESS;
	const uint32_t Nxy = Nx * Ny;

	cl_ushort zero = 0;

	size_t osa_length = loop_var_par;
	std::vector<size_t> cumsum((num_devices_context + 1), 0);
	std::vector<size_t> length(num_devices_context, 0);
	for (cl_uint i = 0; i < num_devices_context; i++) {
		length[i] = osa_length;
		cumsum[i + 1] = cumsum[i] + length[i];
	}

	cl_mem* d_z = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_x = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_y = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_pseudos = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_L = (cl_mem*)malloc((num_devices_context) * sizeof(cl_mem));
	cl_mem* d_lor = (cl_mem*)malloc((num_devices_context) * sizeof(cl_mem));

	// Create the necessary buffers
	for (cl_uint i = 0; i < num_devices_context; i++) {
		d_z[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_x[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_y[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_pseudos[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_lor[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		if (raw) {
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = clEnqueueWriteBuffer(commandQueues[i], d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}

		//status = clEnqueueFillBuffer(commandQueues[i], d_lor[i], &zero, sizeof(cl_ushort), 0, sizeof(cl_ushort) * length[num_devices_context + i], 0, NULL, NULL);

		if (raw) {
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}

		status = clFlush(commandQueues[i]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	clSetKernelArg(kernel, 0, sizeof(uint32_t), &Nxy);
	clSetKernelArg(kernel, 1, sizeof(uint32_t), &im_dim);
	clSetKernelArg(kernel, 2, sizeof(uint32_t), &Nx);
	clSetKernelArg(kernel, 3, sizeof(uint32_t), &Ny);
	clSetKernelArg(kernel, 4, sizeof(uint32_t), &Nz);
	clSetKernelArg(kernel, 5, sizeof(float), &dz);
	clSetKernelArg(kernel, 6, sizeof(float), &dx);
	clSetKernelArg(kernel, 7, sizeof(float), &dy);
	clSetKernelArg(kernel, 8, sizeof(float), &bz);
	clSetKernelArg(kernel, 9, sizeof(float), &bx);
	clSetKernelArg(kernel, 10, sizeof(float), &by);
	clSetKernelArg(kernel, 11, sizeof(float), &bzb);
	clSetKernelArg(kernel, 12, sizeof(float), &maxxx);
	clSetKernelArg(kernel, 13, sizeof(float), &maxyy);
	clSetKernelArg(kernel, 14, sizeof(float), &zmax);
	clSetKernelArg(kernel, 15, sizeof(float), &NSlices);
	clSetKernelArg(kernel, 16, sizeof(uint32_t), &size_x);
	clSetKernelArg(kernel, 17, sizeof(uint16_t), &TotSinos);
	clSetKernelArg(kernel, 18, sizeof(uint32_t), &det_per_ring);
	clSetKernelArg(kernel, 19, sizeof(uint8_t), &raw);
	clSetKernelArg(kernel, 20, sizeof(uint32_t), &prows);

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_event event1;

	for (cl_uint i = 0; i < num_devices_context; i++) {
		size_t erotus = length[i] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		const size_t global_size = length[i] + erotus;
		const uint64_t m_size = static_cast<uint64_t>(length[i]);

		clSetKernelArg(kernel, 21, sizeof(cl_mem), &d_pseudos[i]);
		clSetKernelArg(kernel, 22, sizeof(cl_mem), &d_x[i]);
		clSetKernelArg(kernel, 23, sizeof(cl_mem), &d_y[i]);
		clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_z[i]);
		clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_lor[i]);
		clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_L[i]);
		clSetKernelArg(kernel, 27, sizeof(uint64_t), &m_size);
		status = clEnqueueNDRangeKernel(commandQueues[i], kernel, 1, NULL, &global_size, &local_size, 0, NULL, &event1);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}

	//for (cl_uint i = 0; i < num_devices_context; i++) {
	//	clFinish(commandQueues[i]);
	//}

	status = clEnqueueReadBuffer(commandQueues[0], d_lor[0], CL_TRUE, 0, sizeof(uint16_t) * osa_length, lor1, 1, &event1, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	clReleaseEvent(event1);

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clReleaseMemObject(d_z[i]);
		clReleaseMemObject(d_x[i]);
		clReleaseMemObject(d_y[i]);
		clReleaseMemObject(d_pseudos[i]);
		clReleaseMemObject(d_L[i]);
		clReleaseMemObject(d_lor[i]);
	}
}