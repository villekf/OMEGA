/**************************************************************************
* This function computes the precomputation phase in OpenCL.
*
* Copyright(C) 2020 Ville - Veikko Wettenhovi
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

void precomp_siddon(const cl_uint& num_devices_context, const cl::Context& context, const std::vector<cl::CommandQueue>& commandQueues, uint16_t* lor1, const float* z_det,
	const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, 
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t size_x, 
	const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows,
	const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim, const cl::Kernel& kernel, const size_t numel_x, const size_t local_size) {

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

	std::vector<cl::Buffer> d_z(num_devices_context);
	std::vector<cl::Buffer> d_x(num_devices_context);
	std::vector<cl::Buffer> d_y(num_devices_context);
	std::vector<cl::Buffer> d_pseudos(num_devices_context);
	std::vector<cl::Buffer> d_L(num_devices_context);
	std::vector<cl::Buffer> d_lor(num_devices_context);

	// Create the necessary buffers
	for (cl_uint i = 0; i < num_devices_context; i++) {
		d_z[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_x[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_y[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_pseudos[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_lor[i] = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (raw) {
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = commandQueues[i].enqueueWriteBuffer(d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[i].enqueueWriteBuffer(d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[i].enqueueWriteBuffer(d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[i].enqueueWriteBuffer(d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}

		//status = clEnqueueFillBuffer(commandQueues[i], d_lor[i], &zero, sizeof(cl_ushort), 0, sizeof(cl_ushort) * length[num_devices_context + i], 0, NULL, NULL);

		if (raw) {
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t), L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}

		status = commandQueues[i].flush();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	cl::Kernel kernel_ = kernel;
	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	kernel_.setArg(0, Nxy);
	kernel_.setArg(1, im_dim);
	kernel_.setArg(2, Nx);
	kernel_.setArg(3, Ny);
	kernel_.setArg(4, Nz);
	kernel_.setArg(5, dz);
	kernel_.setArg(6, dx);
	kernel_.setArg(7, dy);
	kernel_.setArg(8, bz);
	kernel_.setArg(9, bx);
	kernel_.setArg(10, by);
	kernel_.setArg(11, bzb);
	kernel_.setArg(12, maxxx);
	kernel_.setArg(13, maxyy);
	kernel_.setArg(14, zmax);
	kernel_.setArg(15, NSlices);
	kernel_.setArg(16, size_x);
	kernel_.setArg(17, TotSinos);
	kernel_.setArg(18, det_per_ring);
	kernel_.setArg(19, raw);
	kernel_.setArg(20, prows);

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	std::vector<cl::Event> event1(1);

	for (cl_uint i = 0; i < num_devices_context; i++) {
		size_t erotus = length[i] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		const size_t global_size = length[i] + erotus;
		const uint64_t m_size = static_cast<uint64_t>(length[i]);
		cl::NDRange global(global_size);
		cl::NDRange local(local_size);

		if (DEBUG) {
			mexPrintf("global_size = %u\n", global_size);
			mexPrintf("local_size = %u\n", local_size);
			mexPrintf("m_size = %u\n", m_size);
			mexEvalString("pause(.0001);");
		}

		kernel_.setArg(21, d_pseudos[i]);
		kernel_.setArg(22, d_x[i]);
		kernel_.setArg(23, d_y[i]);
		kernel_.setArg(24, d_z[i]);
		kernel_.setArg(25, d_lor[i]);
		kernel_.setArg(26, d_L[i]);
		kernel_.setArg(27, m_size);
		status = commandQueues[i].enqueueNDRangeKernel(kernel_, 0, global, local, NULL, &event1[0]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		else if (DEBUG) {
			mexPrintf("Kernel complete\n");
			mexEvalString("pause(.0001);");
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	status = commandQueues[0].enqueueReadBuffer(d_lor[0], CL_TRUE, 0, sizeof(uint16_t) * osa_length, lor1, &event1);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
}