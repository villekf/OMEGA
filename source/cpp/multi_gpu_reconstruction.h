/**************************************************************************
* This is a hybrid Objective-C++ / C++ function for the device selection,
* queue creation, program building and kernel creation, as well as the
* output data and kernel release.
* The file is pure C++ for OpenCL (implementations 3 and 5).
* For Metal (METAL preprocessor directive declared) the file contains
* Objective-C syntax.
*
* Copyright(C) 2020-2026 Ville-Veikko Wettenhovi, Niilo Saarlemo
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
#include "ProjectorClass.h"

// Main reconstruction function for implementations 3 and 5
template <typename T, typename C>
inline void reconstruction_multigpu(const float* z_det, const float* x, scalarStruct& inputScalars, Weighting& w_vec, RecMethods& MethodList, const int64_t* pituus,
	const char* header_directory, const float* meas, const float* im, T* output, C* sensIm, const int type = 0, const int no_norm = 1, const float* rand = nullptr, const float* atten = nullptr,
	const float* norm = nullptr, const float* extraCorr = nullptr, const size_t size_gauss = 0, const uint32_t* xy_index = nullptr,
	const uint16_t* z_index = nullptr, const uint16_t* L = nullptr) {

	const C tyyppi = (C)0;
	const size_t nLength = static_cast<size_t>(inputScalars.subsets) * static_cast<size_t>(inputScalars.Nt);
	std::vector<int64_t> length(nLength); // Number of measurements in each subset

	if (DEBUG) {
		mexPrintBase("inputScalars.subsets = %u\n", inputScalars.subsets);
		mexPrintBase("inputScalars.subsetsUsed = %u\n", inputScalars.subsetsUsed);
		mexPrintBase("inputScalars.osa_iter0 = %u\n", inputScalars.osa_iter0);
		mexEval();
	}

	for (size_t kk = 0; kk < nLength; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];
	// Index of the single subset/timestep pair this call processes
	const size_t indD0 = static_cast<size_t>(inputScalars.osa_iter0) + static_cast<size_t>(inputScalars.timestep0) * static_cast<size_t>(inputScalars.subsets);
	uint64_t m_size = length[indD0];
	if (DEBUG) mexPrint("Adding projector");
	STATUS_t status = SUCCESS_VALUE;

	ProjectorClass proj;
	status = proj.addProjector(inputScalars, w_vec, MethodList, header_directory, type);

	if (status != 0)
		return;
	proj.no_norm = no_norm;

	// Check for correct sizes before writing the buffers
	{
		const size_t vecSize = ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			? static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) : 1ULL;
		const size_t lastMeas = static_cast<size_t>(pituus[inputScalars.subsetsUsed]) * vecSize;
		const bool normalizationIndexedData = inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads;
		bool bad = false;
		auto checkSize = [&](const char* name, const size_t have, const size_t need) {
			if (have < need) {
				mexPrintBase("%s: host array has %llu elements but ", name, static_cast<unsigned long long>(have));
				mexPrintBase("%llu are required\n", static_cast<unsigned long long>(need));
				mexEval();
				bad = true;
			}
		};
		if (inputScalars.maskFP)
			checkSize("maskFP", inputScalars.size_maskFP, static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) *
				(inputScalars.maskFPZ > 1 ? static_cast<size_t>(inputScalars.maskFPZ) : 1ULL));
		if (inputScalars.maskBP)
			checkSize("maskBP", inputScalars.size_maskBP, static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) *
				static_cast<size_t>(inputScalars.maskBPZ));
		if (inputScalars.normalization_correction && inputScalars.size_norm > 1ULL)
			checkSize("normalization", inputScalars.size_norm, normalizationIndexedData ? static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) * static_cast<size_t>(inputScalars.nHeads) : lastMeas);
		if (inputScalars.attenuation_correction) {
			if (inputScalars.CTAttenuation)
				checkSize("attenuation image", inputScalars.size_atten, static_cast<size_t>(inputScalars.im_dim[0]));
			else
				checkSize("attenuation", inputScalars.size_atten, lastMeas);
		}
		// Check "measurement" size when needed (backprojection and implementation 3)
		// BDD uses integral image and thus requires a separate one
		if (type == 2) {
			const size_t needMeas = inputScalars.BPType == 5
				? static_cast<size_t>(inputScalars.nRowsD + 1) * static_cast<size_t>(inputScalars.nColsD + 1) * static_cast<size_t>(length[indD0])
				: static_cast<size_t>(length[indD0]) * vecSize * static_cast<size_t>(inputScalars.nBins);
			checkSize("measurements", inputScalars.size_meas, needMeas);
		}
		else if (type == 0) {
			size_t needMeas = 0ULL;
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
				needMeas += static_cast<size_t>(length[kk]) * vecSize * static_cast<size_t>(inputScalars.nBins);
			checkSize("measurements", inputScalars.size_meas, needMeas);
		}
		if (bad) {
			mexPrint("Aborting: one or more inputs are smaller than the reconstruction geometry requires");
			return;
		}
	}

	// Create OpenCL buffers, CUDA arrays or OneAPI buffers (in the future)
	status = proj.createBuffers(inputScalars, w_vec, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, length, MethodList, type);
	if (status != 0)
		return;

	// Input constant data to the kernels
	status = proj.initializeKernel(inputScalars, w_vec);
	if (status != 0)
		return;

#ifndef METAL
	cl::detail::size_t_array region = { { 0, 0, 0 } };
#else
    NS::SharedPtr<MTL::TextureDescriptor> pTextureDesc = NS::TransferPtr(MTL::TextureDescriptor::alloc()->init());
    pTextureDesc->setTextureType(MTL::TextureType::TextureType3D);
    pTextureDesc->setPixelFormat(MTL::PixelFormat::PixelFormatR32Float);
    pTextureDesc->setWidth(0);
    pTextureDesc->setHeight(0);
    pTextureDesc->setDepth(0);
    std::array<NS::UInteger, 3> region = { 0, 0, 0 };
#endif
	int64_t imTot = 0ULL;

	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[indD0];

	// type 0 = Implementation 3
	// type 1 = Implementation 5 forward projection
	// type 2 = Implementation 5 backprojection
	if (type == 1) { // Forward projection A*x
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("inputScalars.osa_iter0 = %u\n", inputScalars.osa_iter0);
			mexPrintBase("im[0] = %f\n", im[0]);
			mexEval();
		}
		proj.d_output = proj.makeDeviceBuffer(sizeof(float) * m_size * inputScalars.nBins, CL_MEM_READ_WRITE, status);
		CHECK(status, "\n", );
		status = proj.fillDeviceBuffer(proj.d_output, 0.f, sizeof(float) * m_size * inputScalars.nBins);
		CHECK(status, "\n", );
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			imTot += (inputScalars.Ny[ii] + 1) * (inputScalars.Nz[ii] + 1) * inputScalars.Nx[ii];
			proj.d_Summ.emplace_back(proj.makeDeviceBuffer(sizeof(T), CL_MEM_READ_WRITE, status));
			CHECK(status, "\n", );
		}
	}
	if (type == 2) { // Backward projection A'*y
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("proj.no_norm = %u\n", proj.no_norm);
			mexPrintBase("inputScalars.nMultiVolumes = %u\n", inputScalars.nMultiVolumes);
			mexEval();
		}
		if (inputScalars.BPType == 5)
			proj.d_output = proj.makeDeviceBuffer(sizeof(float) * static_cast<uint64_t>(inputScalars.nRowsD + 1) * static_cast<uint64_t>(inputScalars.nColsD + 1) * length[indD0], CL_MEM_READ_ONLY, status);
		else
			proj.d_output = proj.makeDeviceBuffer(sizeof(float) * m_size * inputScalars.nBins, CL_MEM_READ_ONLY, status);
		CHECK(status, "\n", );

		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (DEBUG) {
				mexPrintBase("inputScalars.im_dim[ii] = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("sizeof(T) = %u\n", sizeof(T));
				mexPrintBase("sizeof(T) * inputScalars.im_dim[ii] = %u\n", sizeof(T) * inputScalars.im_dim[ii]);
				mexEval();
			}

			proj.vec_opencl.d_rhs_os.emplace_back(proj.makeDeviceBuffer(sizeof(T) * inputScalars.im_dim[ii], CL_MEM_READ_WRITE, status));
			CHECK(status, "\n", );
			if (proj.no_norm == 0) {
				proj.d_Summ.emplace_back(proj.makeDeviceBuffer(sizeof(T) * inputScalars.im_dim[ii], CL_MEM_READ_WRITE, status));
				CHECK(status, "\n", );
				status = proj.fillDeviceBuffer(proj.d_Summ[ii], (T)0, sizeof(T) * inputScalars.im_dim[ii]);
				CHECK(status, "\n", );
			} else {
				proj.d_Summ.emplace_back(proj.makeDeviceBuffer(sizeof(T), CL_MEM_READ_WRITE, status));
				CHECK(status, "\n", );
			}
			status = proj.fillDeviceBuffer(proj.vec_opencl.d_rhs_os[ii], (T)0, sizeof(T) * inputScalars.im_dim[ii]);
			CHECK(status, "\n", );

		}
		if (inputScalars.BPType == 5)
			status = proj.writeDeviceBuffer(proj.d_output, meas, sizeof(float) * static_cast<uint64_t>(inputScalars.nRowsD + 1) * static_cast<uint64_t>(inputScalars.nColsD + 1) * length[indD0]);
		else
			status = proj.writeDeviceBuffer(proj.d_output, meas, sizeof(float) * m_size * inputScalars.nBins);
		CHECK(status, "\n", );
	}
#ifndef METAL // Metal has no support for implementation 3
	else if (type == 0) {
		size_t uu = 0;
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			region[0] = inputScalars.Nx[ii];
			region[1] = inputScalars.Ny[ii];
			region[2] = inputScalars.Nz[ii];
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
			proj.vec_opencl.d_rhs_os.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(C) * inputScalars.im_dim[ii], NULL, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			//if (proj.no_norm == 0) {
			//	proj.d_Summ.emplace_back(cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(C) * inputScalars.im_dim[ii], NULL, &status));
			//	status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_Summ[ii], (C)0, 0, sizeof(C) * inputScalars.im_dim[ii]);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return;
			//	}
			//}
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
		const size_t nSummPerSubset = static_cast<size_t>(inputScalars.nMultiVolumes) + 1ULL;
		proj.d_meas.resize(inputScalars.subsetsUsed);
		if (inputScalars.randoms_correction)
			proj.d_rand.resize(inputScalars.subsetsUsed);
		proj.d_Summ.resize(static_cast<size_t>(inputScalars.subsetsUsed) * nSummPerSubset);
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
			m_size = length[osa_iter];
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];
			proj.d_meas[osa_iter] = proj.makeDeviceBuffer(sizeof(float) * m_size * inputScalars.nBins, CL_MEM_READ_ONLY, status);
			CHECK(status, "\n", );
			status = proj.writeDeviceBuffer(proj.d_meas[osa_iter], &meas[uu], sizeof(float) * m_size * inputScalars.nBins);
			CHECK(status, "\n", );
			if (inputScalars.randoms_correction) {
				proj.d_rand[osa_iter] = proj.makeDeviceBuffer(sizeof(float) * m_size, CL_MEM_READ_ONLY, status);
				CHECK(status, "\n", );
				status = proj.writeDeviceBuffer(proj.d_rand[osa_iter], &rand[uu], sizeof(float) * m_size);
				CHECK(status, "\n", );
			}
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				const size_t sInd = static_cast<size_t>(ii) + static_cast<size_t>(osa_iter) * nSummPerSubset;
				if (proj.no_norm == 0) {
					proj.d_Summ[sInd] = proj.makeDeviceBuffer(sizeof(C) * inputScalars.im_dim[ii], CL_MEM_READ_WRITE, status);
					CHECK(status, "\n", );
					status = proj.fillDeviceBuffer(proj.d_Summ[sInd], (C)0, sizeof(C) * inputScalars.im_dim[ii]);
					CHECK(status, "\n", );
				} else {
					proj.d_Summ[sInd] = proj.makeDeviceBuffer(sizeof(C), CL_MEM_READ_WRITE, status);
					CHECK(status, "\n", );
				}
			}
			uu += m_size * inputScalars.nBins;
		}
	}
	for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
#endif
	for (uint32_t iter = 0; iter < inputScalars.Niter; iter++) {
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
            for (uint32_t timestep = inputScalars.timestep0; timestep < inputScalars.timestepsUsed; timestep++) {
                // Same [osa_iter + timestep * subsets] convention ProjectorClass uses internally
                const size_t indD = static_cast<size_t>(osa_iter) + static_cast<size_t>(timestep) * static_cast<size_t>(inputScalars.subsets);
                m_size = length[indD];
                if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
                    m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[indD];
#if defined(OPENCL)
                if (type == 0) {
                    proj.d_output = cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * m_size * inputScalars.nBins, NULL, &status);
                    CHECK(status, "\n", );
                    status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_output, 0.f, 0, sizeof(float) * m_size * inputScalars.nBins);
                    CHECK(status, "\n", );
                    if (inputScalars.CT) {
                        proj.d_outputCT = cl::Buffer(proj.CLContext, CL_MEM_READ_WRITE, sizeof(float) * m_size, NULL, &status);
                        CHECK(status, "\n", );
                        status = proj.CLCommandQueue[0].enqueueFillBuffer(proj.d_outputCT, 0.f, 0, sizeof(float) * m_size);
                        CHECK(status, "\n", );
                    }
                    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
                        region[0] = inputScalars.Nx[ii];
                        region[1] = inputScalars.Ny[ii];
                        region[2] = inputScalars.Nz[ii];
                        proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
                        if (inputScalars.use_psf) {
                            status = proj.computeConvolutionF(inputScalars, ii);
                            CHECK(status, "\n", );
                            status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(proj.d_imTemp[ii], proj.vec_opencl.d_image_os, 0, proj.origin, region);
                            CHECK(status, "\n", );
                        }
                        else {
                            status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(proj.d_imFinal[ii], proj.vec_opencl.d_image_os, 0, proj.origin, region);
                            CHECK(status, "\n", );
                        }
                        for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
                            proj.CLCommandQueue[i].finish();
                        }
                        status = proj.forwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, ii);
                        CHECK(status, "\n", );
                    }
                    status = proj.computeForward(inputScalars, length, osa_iter);
                    CHECK(status, "\n", );

                }
    #endif
                if (type == 1) {
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
    #ifndef METAL
                            proj.vec_opencl.d_image_os_int = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
                            CHECK(status, "\n", );
                            status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os_int, CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
                            CHECK(status, "\n", );
    #else
	                            pTextureDesc->setWidth(region[0]);
	                            pTextureDesc->setHeight(region[1]);
	                            pTextureDesc->setDepth(region[2]);
	                            proj.vec_opencl.d_image_os_int = NS::TransferPtr(proj.mtlDevice->newTexture(pTextureDesc.get()));
	                            status = proj.vec_opencl.d_image_os_int ? SUCCESS_VALUE : -1;
	                            CHECK(status, "\n", );

	                            MTL::Region mtlRegion = MTL::Region(
	                                0, 0, 0,
	                                region[0], region[1], region[2]
                            );
                            NS::UInteger bytesPerRow = region[0] * 4;
                            NS::UInteger bytesPerImage = bytesPerRow * region[1];
                            proj.vec_opencl.d_image_os_int->replaceRegion(mtlRegion, 0, 0, &im[uu], bytesPerRow, bytesPerImage);
    #endif
                            region[0] = inputScalars.Nx[ii] + 1;
                            region[1] = inputScalars.Nz[ii] + 1;
                            region[2] = inputScalars.Ny[ii];
                            uu += imTot;
                        }

                        if (DEBUG) {
                            mexPrintBase("uu = %u\n", uu);
                            mexEval();
                        }
                        if (inputScalars.useBuffers) {
                            proj.vec_opencl.d_im = proj.makeDeviceBuffer(sizeof(float) * inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], CL_MEM_READ_ONLY, status);
                            CHECK(status, "\n", );
                            status = proj.writeDeviceBuffer(proj.vec_opencl.d_im, &im[uu], sizeof(float) * inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii]);
                            CHECK(status, "\n", );
                        } else {
    #ifndef METAL
                            proj.vec_opencl.d_image_os = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
                            CHECK(status, "\n", );
                            status = proj.CLCommandQueue[0].enqueueWriteImage(proj.vec_opencl.d_image_os, CL_FALSE, proj.origin, region, 0, 0, &im[uu]);
                            CHECK(status, "\n", );
    #else
                            // Standalone image-mode input is uploaded directly to a Metal
                            // texture. Clear any staging buffer so forwardProjection keeps
                            // this texture instead of trying to refresh it from stale data.
                            proj.vec_opencl.d_im.reset();
	                            pTextureDesc->setWidth(region[0]);
	                            pTextureDesc->setHeight(region[1]);
	                            pTextureDesc->setDepth(region[2]);
	                            proj.vec_opencl.d_image_os = NS::TransferPtr(proj.mtlDevice->newTexture(pTextureDesc.get()));
	                            status = proj.vec_opencl.d_image_os ? SUCCESS_VALUE : -1;
	                            CHECK(status, "\n", );

	                            MTL::Region mtlRegion = MTL::Region(
	                                0, 0, 0,
	                                region[0], region[1], region[2]
                            );
                            NS::UInteger bytesPerRow = region[0] * 4;
                            NS::UInteger bytesPerImage = bytesPerRow * region[1];
                            proj.vec_opencl.d_image_os->replaceRegion(mtlRegion, 0, 0, &im[uu], bytesPerRow, bytesPerImage);
    #endif
                        }

                        status = proj.forwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, ii);
                        CHECK(status, "\n", );
                        if (inputScalars.FPType == 5)
                            uu -= imTot;

                        uu += inputScalars.im_dim[ii];
                    }
                }
                if (type == 2 || type == 0) {
                    for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
                        int uu = ii;
                        if (type == 0) {
                            uu += osa_iter * (inputScalars.nMultiVolumes + 1);
                            status = proj.fillDeviceBuffer(proj.vec_opencl.d_rhs_os[ii], (C)0, sizeof(C) * inputScalars.im_dim[ii]);
                            CHECK(status, "\n", );
                            status = proj.backwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, MethodList, false, ii, ii, uu);

                        } else {
                            status = proj.backwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, MethodList, false, ii, uu);
                        }
                        CHECK(status, "\n", );
    #ifndef METAL // Metal has no support for implementation 3
                        if (type == 0) {
                            if (inputScalars.use_psf) {
                                status = proj.computeConvolution(inputScalars, proj.vec_opencl.d_rhs_os[ii], ii, tyyppi);
                                if (status != CL_SUCCESS) {
                                    return;
                                }
                                if (proj.no_norm == 0) {
                                    status = proj.computeConvolution(inputScalars, proj.d_Summ[uu], ii, tyyppi);
                                    if (status != CL_SUCCESS) {
                                        return;
                                    }
                                }
                            }
                            status = proj.computeEstimate(inputScalars, ii, uu);
                            if (status != CL_SUCCESS) {
                                return;
                            }
                        }
    #endif
                    }
                }
            }
        }
		if (type == 0)
			proj.no_norm = 1;
	}
#ifndef METAL
	for (cl_uint i = 0; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
#endif
	if (type == 1) {
		if (DEBUG) {
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
			mexEval();
		}
		status = proj.readDeviceBuffer(proj.d_output, output, sizeof(float) * m_size * inputScalars.nBins);
		CHECK(status, "\n", );
	} else if (type == 2) {
		size_t uu = 0;
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (inputScalars.atomic_64bit)
				status = proj.readDeviceBuffer(proj.vec_opencl.d_rhs_os[ii], &output[uu], sizeof(INT64_t) * inputScalars.im_dim[ii]);
			else if (inputScalars.atomic_32bit)
				status = proj.readDeviceBuffer(proj.vec_opencl.d_rhs_os[ii], &output[uu], sizeof(INT32_t) * inputScalars.im_dim[ii]);
			else
				status = proj.readDeviceBuffer(proj.vec_opencl.d_rhs_os[ii], &output[uu], sizeof(float) * inputScalars.im_dim[ii]);
			CHECK(status, "\n", );
			if (proj.no_norm == 0) {
				if (inputScalars.atomic_64bit)
					status = proj.readDeviceBuffer(proj.d_Summ[ii], &sensIm[uu], sizeof(INT64_t) * inputScalars.im_dim[ii]);
				else if (inputScalars.atomic_32bit)
					status = proj.readDeviceBuffer(proj.d_Summ[ii], &sensIm[uu], sizeof(INT32_t) * inputScalars.im_dim[ii]);
				else
					status = proj.readDeviceBuffer(proj.d_Summ[ii], &sensIm[uu], sizeof(float) * inputScalars.im_dim[ii]);
				CHECK(status, "\n", );
			}
			uu += inputScalars.im_dim[ii];
		}
	}
#ifndef METAL // Metal has no support for implementation 3
	else if (type == 0) {
		size_t uu = 0;
		int ii = 0;
		status = proj.CLCommandQueue[0].enqueueReadBuffer(proj.d_imFinal[ii], CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[ii], &output[uu]);
		CHECK(status, "\n", );
	}

	for (cl_uint i = 0ULL; i < proj.CLCommandQueue.size(); i++) {
		proj.CLCommandQueue[i].finish();
	}
#endif
	return;
}
