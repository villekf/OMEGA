/**************************************************************************
* These functions handle the device selection, queue creation, program
* building and kernel creation, as well as the output data and kernel 
* release.
*
* Copyright(C) 2020-2024 Ville-Veikko Wettenhovi
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
//#include "functions_multigpu.hpp"
#include "ProjectorClass.h"

// Main reconstruction function for implementation 3
template <typename T>
inline void reconstruction_multigpu(const float* z_det, const float* x, scalarStruct& inputScalars, Weighting& w_vec, RecMethods& MethodList, const int64_t* pituus,
	const char* header_directory, const float* meas, const float* im, T* output, T* sensIm, const int type = 0, const int no_norm = 1, const float* rand = nullptr, const float* atten = nullptr, 
	const float* norm = nullptr, const float* extraCorr = nullptr, const size_t size_gauss = 0, const uint32_t* xy_index = nullptr,
	const uint16_t* z_index = nullptr, const uint16_t* L = nullptr) {

	// Number of measurements in each subset
	std::vector<int64_t> length(inputScalars.subsetsUsed);


	if (DEBUG) {
		mexPrintBase("inputScalars.subsets = %u\n", inputScalars.subsets);
		mexPrintBase("inputScalars.subsetsUsed = %u\n", inputScalars.subsetsUsed);
		mexEval();
	}

	for (uint32_t kk = 0; kk < inputScalars.subsetsUsed; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];
	uint64_t m_size = length[inputScalars.osa_iter0];

	cl_int status = CL_SUCCESS;

	if (DEBUG) {
		mexPrint("Adding projector");
	}
	ProjectorClass proj;
	status = proj.addProjector(inputScalars, w_vec, MethodList, header_directory, type);
	if (status != 0)
		return;
	proj.no_norm = no_norm;

	// Create OpenCL buffers, CUDA arrays or OneAPI buffers
	status = proj.createBuffers(inputScalars, w_vec, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, length, MethodList, type);
	if (status != 0)
		return;

	// Input constant data to the kernels
	status = proj.initializeKernel(inputScalars, w_vec);
	if (status != 0)
		return;
	status = proj.setDynamicKernelData(inputScalars, w_vec);
	if (status != 0)
		return;
	cl::detail::size_t_array region = { { 0, 0, 0 } };
	int64_t imTot = 0ULL;

	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[inputScalars.osa_iter0];
	if (type == 1) {
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("inputScalars.osa_iter0 = %u\n", inputScalars.osa_iter0);
			mexPrintBase("im[0] = %f\n", im[0]);
			mexEval();
		}
		proj.d_output = cl::Buffer(proj.CLContext, CL_MEM_WRITE_ONLY, sizeof(float) * m_size * inputScalars.nBins, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
			imTot += (inputScalars.Ny[ii] + 1) * (inputScalars.Nz[ii] + 1) * inputScalars.Nx[ii];
		//cl::detail::size_t_array region = { { 0, 0, 0 } };
		//size_t uu = 0;
		//for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		//	region[0] = inputScalars.Nx[ii];
		//	region[1] = inputScalars.Ny[ii];
		//	region[2] = inputScalars.Nz[ii];
		//	proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
		//	status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os, CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return;
		//	}
		//	uu += inputScalars.im_dim[ii];
		//}
		status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_output, 0.f, 0, sizeof(float) * m_size * inputScalars.nBins);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else if (type == 2) {
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("proj.no_norm = %u\n", proj.no_norm);
			mexPrintBase("inputScalars.nMultiVolumes = %u\n", inputScalars.nMultiVolumes);
			mexEval();
		}
		if (inputScalars.BPType == 5)
			proj.d_output = cl::Buffer(proj.CLContext, CL_MEM_READ_ONLY, sizeof(float) * static_cast<uint64_t>(inputScalars.nRowsD + 1) * static_cast<uint64_t>(inputScalars.nColsD + 1) * length[inputScalars.osa_iter0], NULL, &status);
		else
			proj.d_output = cl::Buffer(proj.CLContext, CL_MEM_READ_ONLY, sizeof(float) * m_size * inputScalars.nBins, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (DEBUG) {
				mexPrintBase("inputScalars.im_dim[ii] = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("sizeof(T) = %u\n", sizeof(T));
				mexPrintBase("sizeof(T) * inputScalars.im_dim[ii] = %u\n", sizeof(T) * inputScalars.im_dim[ii]);
				mexEval();
			}
			proj.vec_opencl.d_rhs_os.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_WRITE_ONLY, sizeof(T) * inputScalars.im_dim[ii], NULL, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (proj.no_norm == 0) {
				proj.d_Summ.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_WRITE_ONLY, sizeof(T) * inputScalars.im_dim[ii], NULL, &status));
				status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_Summ[ii], (T)0, 0, sizeof(T) * inputScalars.im_dim[ii]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else
				proj.d_Summ.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_WRITE_ONLY, sizeof(T), NULL, &status));
			status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.vec_opencl.d_rhs_os[ii], (T)0, 0, sizeof(T) * inputScalars.im_dim[ii]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (inputScalars.BPType == 5)
			status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_output, CL_FALSE, 0, sizeof(float) * static_cast<uint64_t>(inputScalars.nRowsD + 1) * static_cast<uint64_t>(inputScalars.nColsD + 1) * length[inputScalars.osa_iter0], meas);
		else
			status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_output, CL_FALSE, 0, sizeof(float) * m_size * inputScalars.nBins, meas);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else if (type == 0) {
		size_t uu = 0;
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			region[0] = inputScalars.Nx[ii];
			region[1] = inputScalars.Ny[ii];
			region[2] = inputScalars.Nz[ii];
			//status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os[ii], CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
			//if (status != CL_SUCCESS) {
			//	getErrorString(status);
			//	return;
			//}
			//uu += inputScalars.im_dim[ii];
			proj.d_imFinal.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * inputScalars.im_dim[ii], NULL, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (inputScalars.use_psf) {
				proj.d_imTemp.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * inputScalars.im_dim[ii], NULL, &status));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			proj.vec_opencl.d_rhs_os.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(T)* inputScalars.im_dim[ii], NULL, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (proj.no_norm == 0) {
				proj.d_Summ.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(T) * inputScalars.im_dim[ii], NULL, &status));
				status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_Summ[ii], (T)0, 0, sizeof(T) * inputScalars.im_dim[ii]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.vec_opencl.d_rhs_os[ii], (T)0, 0, sizeof(T) * inputScalars.im_dim[ii]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_imFinal[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &im[uu]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (inputScalars.use_psf) {
				status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_imTemp[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &im[uu]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			uu += inputScalars.im_dim[ii];
			if (inputScalars.use_psf) {
				proj.d_g = cl::Buffer(proj.CLContext, CL_MEM_READ_ONLY, sizeof(float) * size_gauss, NULL, &status);
				status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_g, CL_FALSE, 0, sizeof(float) * size_gauss, inputScalars.gaussian);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
		uu = 0;
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
			m_size = length[osa_iter];
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];
			proj.d_meas.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_ONLY, sizeof(float) * m_size * inputScalars.nBins, NULL, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_meas[osa_iter], CL_FALSE, 0, sizeof(float) * m_size * inputScalars.nBins, &meas[uu]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (inputScalars.randoms_correction) {
				proj.d_rand.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_ONLY, sizeof(float) * m_size, NULL, &status));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = proj.CLCommandQueue[0].enqueueWriteBuffer(proj.d_rand[osa_iter], CL_FALSE, 0, sizeof(float) * m_size, &rand[uu]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			uu += m_size * inputScalars.nBins;
		}
	}
	for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
	for (uint32_t iter = 0; iter < inputScalars.Niter; iter++) {
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
			uint32_t subIter = osa_iter;
			m_size = length[osa_iter];
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];
			if (type == 0) {
				proj.d_output = cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * m_size * inputScalars.nBins, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_output, 0.f, 0, sizeof(float) * m_size * inputScalars.nBins);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (inputScalars.CT) {
					proj.d_outputCT = cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * m_size, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_outputCT, 0.f, 0, sizeof(float) * m_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					region[0] = inputScalars.Nx[ii];
					region[1] = inputScalars.Ny[ii];
					region[2] = inputScalars.Nz[ii];
					proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
					if (inputScalars.use_psf) {
						status = proj.computeConvolution(inputScalars, ii);
						if (status != CL_SUCCESS) {
							return;
						}
						status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(proj.d_imTemp[ii], proj.vec_opencl.d_image_os, 0, proj.origin, region);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrint("Image copy failed\n");
							return;
						}
					}
					else {
						status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(proj.d_imFinal[ii], proj.vec_opencl.d_image_os, 0, proj.origin, region);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrint("Image copy failed\n");
							return;
						}
					}
					for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
						proj.CLCommandQueue[i].finish();
					}
					status = proj.forwardProjection(inputScalars, w_vec, osa_iter, length, m_size, ii);
					if (status != CL_SUCCESS) {
						return;
					}
				}
				status = proj.computeForward(inputScalars, length, osa_iter);
				if (status != CL_SUCCESS) {
					return;
				}
			}
			else if (type == 1) {
				size_t uu = 0;
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					region[0] = inputScalars.Nx[ii];
					region[1] = inputScalars.Ny[ii];
					region[2] = inputScalars.Nz[ii];
					if (inputScalars.FPType == 5) {
						region[0] = inputScalars.Ny[ii] + 1;
						region[1] = inputScalars.Nz[ii] + 1;
						region[2] = inputScalars.Nx[ii];
						if (DEBUG) {
							mexPrintBase("region[0] = %u\n", region[0]);
							mexPrintBase("region[1] = %u\n", region[1]);
							mexPrintBase("region[2] = %u\n", region[2]);
							mexEval();
						}
						proj.vec_opencl.d_image_os_int = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os_int, CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						region[0] = inputScalars.Nx[ii] + 1;
						region[1] = inputScalars.Nz[ii] + 1;
						region[2] = inputScalars.Ny[ii];
						uu += imTot;
					}
					if (DEBUG) {
						mexPrintBase("uu = %u\n", uu);
						mexEval();
					}
					proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os, CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					status = proj.forwardProjection(inputScalars, w_vec, osa_iter, length, m_size, ii);
					if (status != CL_SUCCESS) {
						return;
					}
					if (inputScalars.FPType == 5)
						uu -= imTot;
					uu += inputScalars.im_dim[ii];
				}
			}
			if (type == 2 || type == 0) {
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					status = proj.backwardProjection(inputScalars, w_vec, osa_iter, length, m_size, false, ii, ii);
					if (status != CL_SUCCESS) {
						return;
					}
					if (type == 0) {
						status = proj.computeEstimate(inputScalars, ii);
						if (status != CL_SUCCESS) {
							return;
						}
					}
				}
			}
		}
	}
	for (cl_uint i = 0; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
	if (type == 1) {
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
			mexEval();
		}
		status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_output, CL_FALSE, 0, sizeof(float) * m_size * inputScalars.nBins, output);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else if (type == 2) {
		size_t uu = 0;
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (inputScalars.atomic_64bit)
				status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.vec_opencl.d_rhs_os[ii], CL_FALSE, 0, sizeof(cl_long) * inputScalars.im_dim[ii], &output[uu]);
			else if (inputScalars.atomic_32bit)
				status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.vec_opencl.d_rhs_os[ii], CL_FALSE, 0, sizeof(cl_int) * inputScalars.im_dim[ii], &output[uu]);
			else
				status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.vec_opencl.d_rhs_os[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &output[uu]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (proj.no_norm == 0) {
				if (inputScalars.atomic_64bit)
					status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_Summ[ii], CL_FALSE, 0, sizeof(cl_long) * inputScalars.im_dim[ii], &sensIm[uu]);
				else if (inputScalars.atomic_32bit)
					status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_Summ[ii], CL_FALSE, 0, sizeof(cl_int) * inputScalars.im_dim[ii], &sensIm[uu]);
				else
					status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_Summ[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &sensIm[uu]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			uu += inputScalars.im_dim[ii];
		}
	}
	else if (type == 0) {
		size_t uu = 0;
		int ii = 0;
		//for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_imFinal[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &output[uu]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		//	uu += inputScalars.im_dim[ii];
		//}
	}

	for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
	return;
}
