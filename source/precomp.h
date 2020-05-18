#pragma once
#include <cstdint>
#include <stdio.h>
#include "opencl_error.hpp"

void precomp_siddon(const cl_uint& num_devices_context, const cl_context& context, const cl_command_queue* commandQueues, uint16_t* lor1, const float* z_det,
	const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx,
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t size_x,
	const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows,
	const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim, const cl_kernel& kernel, const size_t numel_x, const size_t local_size);