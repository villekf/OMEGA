/**************************************************************************
* Header for matrix-free CUDA functions
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
#pragma once
#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "functions.hpp"
//#include <driver_types.h>

#pragma pack(1)

#define TH 100000000000.f
#define DEBUG false

//#undef max
//#undef min

const char* getErrorString(CUresult error);

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(CUresult code, const char* file, int line, bool abort = true)
{
	if (code != CUDA_SUCCESS)
	{
		const char* errstr;
		cuGetErrorString(code, &errstr);
		mexPrintf("GPUassert: %s %s %d\n", errstr, file, line);
	}
}

// Struct for the various estimates used in OpenCL kernels
typedef struct _CUDA_im_vectors {
	CUdeviceptr* d_im_mlem, *d_rhs_mlem, *d_im_os, *d_rhs_os;
} CUDA_im_vectors;

// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
//typedef struct _RecMethodsOpenCL {
//	cl_char MLEM, OSEM, MRAMLA, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
//	cl_char MRP, Quad, L, FMH, WeightedMean, TV, AD, APLS, TGV, NLM;
//	cl_char OSLMLEM, OSLOSEM, MBSREM, BSREM, ROSEMMAP, RBIMAP;
//	cl_char OSLCOSEM;
//} RecMethodsOpenCL;

// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
typedef struct _RecMethodsOpenCL {
	char MLEM, OSEM, MRAMLA, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
	char MRP, Quad, Huber, L, FMH, WeightedMean, TV, AD, APLS, TGV, NLM;
	char OSLMLEM, OSLOSEM, MBSREM, BSREM, ROSEMMAP, RBIOSL, OSLCOSEM;
} RecMethodsOpenCL;

// Update the OpenCL inputs for the current iteration/subset
void update_cuda_inputs(AF_im_vectors & vec, CUDA_im_vectors &vec_opencl, const bool mlem, const uint32_t im_dim, const uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const RecMethods MethodList, const bool atomic_64bit, const bool use_psf);

// For OpenCL
void OpenCLRecMethods(const RecMethods& MethodList, RecMethodsOpenCL& MethodListOpenCL);

CUresult createAndWriteBuffers(CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_L,
	std::vector<CUdeviceptr>& d_zindex, std::vector<CUdeviceptr>& d_xyindex, std::vector<CUdeviceptr>& d_Sino, std::vector<CUdeviceptr>& d_sc_ra,
	const uint32_t size_x, const size_t size_z, const uint32_t TotSinos, const size_t size_atten, const size_t size_norm, const size_t size_scat, const uint32_t prows,
	std::vector<size_t>& length, const float* x, const float* y, const float* z_det, const uint32_t* xy_index, const uint16_t* z_index,
	const uint16_t* lor1, const uint16_t* L, const float* Sin, const uint8_t raw, const uint32_t subsets, const int64_t* pituus, const float* atten,
	const float* norm, const float* scat, const uint32_t* pseudos, const float* V, CUdeviceptr& d_atten, std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat, CUdeviceptr& d_pseudos, CUdeviceptr& d_V,
	CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, CUdeviceptr& d_zcenter, const float* x_center, const float* y_center, const float* z_center,
	const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const size_t size_of_x, const size_t size_V,
	const bool randoms_correction, const mxArray* sc_ra, const bool precompute, CUdeviceptr& d_lor_mlem, CUdeviceptr& d_L_mlem, CUdeviceptr& d_zindex_mlem,
	CUdeviceptr& d_xyindex_mlem, CUdeviceptr& d_Sino_mlem, CUdeviceptr& d_sc_ra_mlem, CUdeviceptr& d_reko_type, CUdeviceptr& d_reko_type_mlem, const bool osem_bool,
	const bool mlem_bool, const size_t koko, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const uint32_t n_rekos, const uint32_t n_rekos_mlem,
	CUdeviceptr& d_norm_mlem, CUdeviceptr& d_scat_mlem, const bool TOF, const int64_t nBins, const bool loadTOF, CUdeviceptr& d_TOFCenter, const float* TOFCenter);

// Prepass phase for MRAMLA, MBSREM, COSEM, ACOSEM, ECOSEM, RBI
void MRAMLA_prepass_CUDA(const uint32_t subsets, const uint32_t im_dim, const int64_t* pituus, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_zindex,
	std::vector<CUdeviceptr>& d_xyindex, Weighting& w_vec, std::vector<af::array>& Summ, std::vector<CUdeviceptr>& d_Sino, size_t koko_l, af::array& cosem,
	af::array& C_co, af::array& C_aco, af::array& C_osl, uint32_t alku, std::vector<CUdeviceptr>& d_L, uint8_t raw, RecMethodsOpenCL& MethodList,
	std::vector<size_t> length, uint8_t compute_norm_matrix, std::vector<CUdeviceptr>& d_sc_ra, af::array& E, const uint32_t det_per_ring, CUdeviceptr& d_pseudos,
	const uint32_t prows, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dz, const float dx, const float dy, const float bz, const float bx,
	const float by, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z, const uint32_t size_x,
	const uint32_t TotSinos, CUdeviceptr& d_atten,
	std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat, const float epps, const uint32_t Nxy, const float tube_width, const float crystal_size_z, const float bmin, const float bmax,
	const float Vmax, CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, CUdeviceptr& d_zcenter, CUdeviceptr& d_V, const float dc_z, const uint16_t n_rays, const uint16_t n_rays3D,
	const bool precompute, const uint32_t projector_type, const CUstream& af_cuda_stream, const float global_factor, CUdeviceptr& d_reko_type, CUfunction& kernel_mbsrem,
	const bool atomic_64bit, const bool use_psf, const af::array& g, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins,
	const bool randoms_correction, const float sigma_x, CUdeviceptr& d_TOFCenter, const uint32_t Nt = 1U);

nvrtcResult createProgramCUDA(const bool verbose, const char* k_path, const char* fileName,
	nvrtcProgram& program_os, nvrtcProgram& program_ml, nvrtcProgram& program_mbsrem, bool& atomic_64bit, const char* header_directory,
	const uint32_t projector_type, const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction,
	const uint32_t normalization_correction, const int32_t dec, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D,
	const RecMethods MethodList, const bool osem_bool, const bool mlem_bool, const uint32_t n_rekos, const uint32_t n_rekos_mlem,
	const Weighting& w_vec, const uint32_t osa_iter0, const float cr_pz, const float dx, const bool use_psf, const uint32_t scatter, const uint32_t randoms_correction, 
	const bool TOF, const int64_t nBins);

nvrtcResult buildProgramCUDA(const bool verbose, const char* k_path, nvrtcProgram& program, bool& atomic_64bit, std::vector<const char*> &options);

nvrtcResult createKernelsCUDA(const bool verbose, nvrtcProgram& program_os, nvrtcProgram& program_ml, nvrtcProgram& program_mbsrem, CUfunction& kernel_os, CUfunction& kernel_ml,
	CUfunction& kernel_mbsrem, CUfunction& kernelNLM, const bool osem_bool, const bool mlem_bool, const RecMethods& MethodList, const Weighting& w_vec, const bool precompute, const uint32_t projector_type,
	const uint16_t n_rays, const uint16_t n_rays3D);

void computeOSEstimatesCUDA(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, RecMethodsOpenCL& MethodListOpenCL, const uint32_t im_dim,
	af::array* testi, const float epps, const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const TVdata& data, std::vector<size_t>& length, std::vector<CUdeviceptr>& d_Sino, bool& break_iter, af::array& pj3, const uint32_t n_rekos2,
	const int64_t* pituus, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_zindex, std::vector<CUdeviceptr>& d_xyindex, const CUstream& af_cuda_stream,
	std::vector<af::array>& Summ, CUfunction& kernel_mramla, std::vector<CUdeviceptr>& d_L, const uint8_t raw, const size_t koko, const bool atomic_64bit,
	const uint8_t compute_norm_matrix, std::vector<CUdeviceptr>& d_sc_ra, af::array& E, std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat,
	const uint32_t det_per_ring, CUdeviceptr& d_pseudos, const uint32_t prows, const float dz, const float dx, const float dy, const float bz, const float bx,
	const float by, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z,
	const uint32_t size_x, const uint32_t TotSinos, CUdeviceptr& d_atten, const uint32_t Nxy, const float tube_width, const float crystal_size_z, const float bmin,
	const float bmax, const float Vmax, CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, CUdeviceptr& d_zcenter, CUdeviceptr& d_V, const float dc_z, const uint16_t n_rays,
	const uint16_t n_rays3D, const bool precompute, const uint32_t projector_type, const float global_factor, CUdeviceptr& d_reko_type, CUfunction& kernel_mbsrem,
	const bool use_psf, const af::array& g, const kernelStruct& CUDAStruct, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins,
	const bool randoms_correction, const float sigma_x, CUdeviceptr& d_TOFCenter);