/**************************************************************************
* Header for matrix-free OpenCL functions
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
#include "functions.hpp"

#pragma pack(1)

//#undef max
//#undef min


// Struct for the various estimates used in OpenCL kernels
typedef struct _OpenCL_im_vectors {
	cl::Buffer d_im_mlem, d_rhs_mlem, d_im_os, d_rhs_os;
} OpenCL_im_vectors;

// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
typedef struct _RecMethodsOpenCL {
	cl_char MLEM, OSEM, MRAMLA, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
	cl_char MRP, Quad, Huber, L, FMH, WeightedMean, TV, AD, APLS, TGV, NLM;
	cl_char OSLMLEM, OSLOSEM, MBSREM, BSREM, ROSEMMAP, RBIOSL, OSLCOSEM;
} RecMethodsOpenCL;

// Update the OpenCL inputs for the current iteration/subset
void update_opencl_inputs(AF_im_vectors& vec, OpenCL_im_vectors& vec_opencl, const bool mlem, const uint32_t im_dim, const uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const RecMethods MethodList, const bool atomic_64bit, const bool use_psf);

// For OpenCL
void OpenCLRecMethods(const RecMethods &MethodList, RecMethodsOpenCL &MethodListOpenCL);

// Save the OpenCL program binary
//cl_int SaveProgramBinary(const bool verbose, const char* k_path, cl_context af_context, cl_device_id af_device_id, const char* fileName, cl_program &program, 
//	bool& atomic_64bit, const uint32_t device, const char* header_directory, const bool force_build, const uint32_t projector_type, const float crystal_size_z,
//	const bool precompute);

// Load the OpenCL binary and create an OpenCL program from it
//cl_int CreateProgramFromBinary(cl_context af_context, cl_device_id af_device_id, FILE *fp, cl_program &program);

cl_int createKernels(cl::Kernel& kernel_ml, cl::Kernel& kernel, cl::Kernel& kernel_mramla, cl::Kernel& kernelNLM, const bool osem_bool, const cl::Program& program_os, const cl::Program& program_ml,
	const cl::Program& program_mbsrem, const RecMethods MethodList, const Weighting w_vec, const uint32_t projector_type, const bool mlem_bool, const bool precompute,
	const uint16_t n_rays, const uint16_t n_rays3D);

cl_int createAndWriteBuffers(cl::Buffer& d_x, cl::Buffer& d_y, cl::Buffer& d_z, std::vector<cl::Buffer>& d_lor, std::vector<cl::Buffer>& d_L, std::vector<cl::Buffer>& d_zindex,
	std::vector<cl::Buffer>& d_xyindex, std::vector<cl::Buffer>& d_Sino, std::vector<cl::Buffer>& d_sc_ra, const uint32_t size_x, const size_t size_z,
	const uint32_t TotSinos, const size_t size_atten, const size_t size_norm, const size_t size_scat, const uint32_t prows, std::vector<size_t>& length, const float* x, const float* y,
	const float* z_det, const uint32_t* xy_index, const uint16_t* z_index, const uint16_t* lor1, const uint16_t* L, const float* Sin, const uint8_t raw,
	cl::Context& af_context, const uint32_t subsets, const int64_t* pituus, const float* atten, const float* norm, const float* scat, const uint32_t* pseudos, const float* V,
	cl::CommandQueue& af_queue, cl::Buffer& d_atten, std::vector<cl::Buffer>& d_norm, std::vector<cl::Buffer>& d_scat, cl::Buffer& d_pseudos, cl::Buffer& d_V, cl::Buffer& d_xcenter, 
	cl::Buffer& d_ycenter, cl::Buffer& d_zcenter, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_center_z, const size_t size_of_x, const size_t size_V, const bool atomic_64bit, const bool randoms_correction, const mxArray* sc_ra, const bool precompute, 
	cl::Buffer& d_lor_mlem,	cl::Buffer& d_L_mlem, cl::Buffer& d_zindex_mlem, cl::Buffer& d_xyindex_mlem, cl::Buffer& d_Sino_mlem, cl::Buffer& d_sc_ra_mlem, cl::Buffer& d_reko_type, 
	cl::Buffer& d_reko_type_mlem, const bool osem_bool,	const bool mlem_bool, const size_t koko, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const uint32_t n_rekos, 
	const uint32_t n_rekos_mlem, cl::Buffer& d_norm_mlem, cl::Buffer& d_scat_mlem, const bool TOF, const int64_t nBins, const bool loadTOF, cl::Buffer& d_TOFCenter, const float* TOFCenter);

// Prepass phase for MRAMLA, MBSREM, COSEM, ACOSEM, ECOSEM, RBI
void MRAMLA_prepass(const uint32_t subsets, const uint32_t im_dim, const int64_t* pituus, const std::vector<cl::Buffer>& lor, const std::vector<cl::Buffer>& zindex,
	const std::vector<cl::Buffer>& xindex, cl::Program program, const cl::CommandQueue& af_queue, const cl::Context af_context, Weighting& w_vec,
	std::vector<af::array>& Summ, std::vector<cl::Buffer>& d_Sino, const size_t koko_l, af::array& cosem, af::array& C_co,
	af::array& C_aco, af::array& C_osl, const uint32_t alku, cl::Kernel& kernel_mramla, const std::vector<cl::Buffer>& L, const uint8_t raw,
	const RecMethodsOpenCL MethodListOpenCL, const std::vector<size_t> length, const bool atomic_64bit, const cl_uchar compute_norm_matrix,
	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
	const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins,
	const size_t koko, const bool randoms_correction, const uint32_t Nt = 1U);

void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dx,
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
	const float NSlices, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const size_t loop_var_par, const char* k_path,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z,
	const char* fileName, const uint32_t device, const size_t numel_x, const char* header_directory);

cl_int createProgram(const bool verbose, const char* k_path, cl::Context& af_context, cl::Device& af_device_id, const char* fileName,
	cl::Program& program_os, cl::Program& program_ml, cl::Program& program_mbsrem, bool& atomic_64bit, const uint32_t device, const char* header_directory,
	const uint32_t projector_type, const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction,
	const uint32_t normalization_correction, const int32_t dec, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D,
	const bool find_lors, const RecMethods MethodList, const bool osem_bool, const bool mlem_bool, const uint32_t n_rekos, const uint32_t n_rekos_mlem,
	const Weighting& w_vec, const uint32_t osa_iter0, const float cr_pz, const float dx, const bool use_psf, const uint32_t scatter, const uint32_t randoms_correction,
	const bool TOF, const int64_t nBins);

cl_int buildProgram(const bool verbose, std::string content, cl::Context& af_context, cl::Device& af_device_id, cl::Program& program,
	bool& atomic_64bit, std::string options);
//cl_int buildProgram(const bool verbose, const char* k_path, cl::Context& af_context, cl::Device& af_device_id, cl::Program& program,
//	bool& atomic_64bit, std::string options);

void computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, af::array* testi, const float epps,
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const TVdata& data, std::vector<size_t>& length, std::vector<cl::Buffer>& d_Sino, bool& break_iter, af::array& pj3, const uint32_t n_rekos2, const int64_t* pituus,
	const std::vector<cl::Buffer>& d_lor, const std::vector<cl::Buffer>& d_zindex, const std::vector<cl::Buffer>& d_xyindex, cl::Program& program_mbsrem, const cl::CommandQueue& af_queue,
	const cl::Context& af_context, std::vector<af::array>& Summ, cl::Kernel& kernel_mramla, const std::vector<cl::Buffer>& d_L, const uint8_t raw,
	const RecMethodsOpenCL& MethodListOpenCL, const size_t koko, const bool atomic_64bit, const cl_uchar compute_norm_matrix, cl::Kernel& kernelNLM,
	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
	const af::array& g, const kernelStruct& OpenCLStruct, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, const bool randoms_correction);