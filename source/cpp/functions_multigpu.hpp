/**************************************************************************
* Header for matrix-free OpenCL functions for the multi-GPU/device case
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
#pragma once
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <vector>
#include "precomp.h"
#include <string>
#include <cmath>
#include "ProjectorClass.h"

#define DEBUG true

void OSEM_MLEM(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl::Context& context,
	const std::vector<cl::CommandQueue>& commandQueues, const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, 
	const mxArray* Sin, const mxArray* sc_ra, scalarStruct inputScalars, const mxArray* options, const int64_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const bool verbose, const float* atten, const size_t size_atten,
	const float* norm, const size_t size_norm, const uint16_t* L, const uint32_t im_dim, const cl::Kernel& kernel,
	const cl::Kernel& kernel_sum, const cl::Kernel& kernel_mlem, const cl::Kernel& kernel_convolution, const cl::Kernel& kernel_convolution_f,
	const cl::Kernel& kernel_vectorMult, const cl::Kernel& kernel_vectorDiv, cl::Kernel& kernelBP, const size_t numel_x, const float* x_center,
	const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z,
	const bool atomic_64bit, const bool atomic_32bit, const cl_uchar compute_norm_matrix, mxArray* cell, const bool osem_bool, const float* V,
	const size_t size_V, const size_t local_size[], const float* gaussian, const size_t size_gauss, const int64_t TOFSize, const float* TOFCenter,
	const cl::vector<cl::Device> devices);

void f_b_project(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl::Context& context, 
	const std::vector<cl::CommandQueue> & commandQueues, const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, 
	const float* rhs, const mxArray* sc_ra, scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index,
	const uint16_t* z_index, const bool verbose, const float* atten, const size_t size_atten, const float* norm, 
	const size_t size_norm,	const uint16_t* L, const uint32_t im_dim, const cl::Kernel& kernel_sum, const cl::Kernel& kernel, 
	mxArray* output, const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x, const float* x_center, const float* y_center, 
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const mxArray* Sin, 
	const bool atomic_64bit, const bool atomic_32bit, const float* V, const size_t size_V, const size_t local_size[], const mxArray* options, 
	const int64_t TOFSize, const float* TOFCenter);

cl_int clGetPlatformsContext(const uint32_t device, const float kerroin, cl::Context& context, size_t& size, int& cpu_device,
	cl_uint& num_devices_context, cl::vector<cl::Device> & devices, bool& atomic_64bit, cl_uchar& compute_norm_matrix, const uint32_t Nxyz, const scalarStruct inputScalars);

cl_int clGetPlatformsContextSingle(const uint32_t device, cl::Context& context, cl_uint& num_devices_context, cl::vector<cl::Device> & devices);

cl_int ClBuildProgramGetQueues(cl::Program& program, cl::Program& programAux, const char* k_path, const cl::Context context, const cl_uint num_devices_context,
	const cl::vector<cl::Device> & devices, const bool verbose, std::vector<cl::CommandQueue> & commandQueues, bool& atomic_64bit, const bool atomic_32bit, 
	const scalarStruct inputScalars, const char* header_directory, const size_t local_size[], const bool find_lors, 
	const uint8_t listmode = 0, const bool CT = false);

void reconstruction_multigpu(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const mxArray* Sin, const mxArray* sc_ra, 
	const mxArray* options, scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, 
	mxArray* cell, const bool verbose, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint8_t* rekot, const char* k_path, 
	const size_t size_rekot, const uint16_t* L, const bool osem_bool, const char* fileName, const uint32_t device, float kerroin, 
	const size_t numel_x, const float* x_center, const float* y_center,	const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_center_z, const char* header_directory, const bool use_64bit_atomics, const float* V, const size_t size_V, size_t local_size[], 
	const float* gaussian, const size_t size_gauss, const int64_t TOFSize, const float* TOFCenter);

void reconstruction_f_b_proj(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* rhs, const mxArray* sc_ra,
	scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const bool verbose, 
	const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const char* k_path, const uint16_t* L, 
	const char* fileName, const uint32_t device, float kerroin, mxArray* output, const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x,
	const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, 
	const char* header_directory, const mxArray* Sin, const bool use_64bit_atomics, const float* V, const size_t size_V, size_t local_size[], const mxArray* options, 
	const int64_t TOFSize, const float* TOFCenter);

void find_LORs(uint16_t* lor, const float* z_det, const float* x, const bool verbose, const size_t loop_var_par, 
	scalarStruct inputScalars, const char* k_path, const uint16_t* L, const char* fileName, const uint32_t device, 
	const size_t numel_x, const char* header_directory, const size_t local_size[]);

void precomp_siddon(const cl_uint& num_devices_context, const cl::Context& context, const std::vector<cl::CommandQueue>& commandQueues, uint16_t* lor1, const float* z_det,
	const float* x, scalarStruct inputScalars, const bool verbose, const size_t loop_var_par, const uint16_t* L, const uint32_t im_dim, const cl::Kernel& kernel, 
	const size_t numel_x, const size_t local_size[]);