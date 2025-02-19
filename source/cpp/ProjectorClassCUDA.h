/*******************************************************************************************************************************************
* Class object for forward and backward projections.
*
* Copyright (C) 2022-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License  along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/
#pragma once
#include "structs.h"
/// <summary>
/// Class object for forward and backward projections. CUDA version
/// </summary>
class ProjectorClass {
	//private:
		// Local size
	size_t local_size[3];
	size_t local_sizePrior[3];
	// Kernel input indices
	unsigned int kernelInd_MRAMLA = 0;
	unsigned int kernelIndFP = 0;
	unsigned int kernelIndBP = 0;
	unsigned int kernelIndFPSubIter = 0;
	unsigned int kernelIndBPSubIter = 0;
	// Crystal pitch
	float2 dPitch;
	// Image dimensions
	int3 d_NOrig, d_NPrior;
	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusPrior[3];
	size_t erotusPriorEFOV[3];
	size_t erotusSens[3];
	size_t kSize = 0ULL;
	// Local and global sizes
	unsigned int local[3];
	unsigned int global[3];
	unsigned int localPrior[3];
	unsigned int globalPrior[3];
	unsigned int globalPriorEFOV[3];
	struct CUDAMemAlloc {
		bool xC = false;
		bool yC = false;
		bool zC = false;
		bool V = false;
		bool TOF = false;
		bool eFOV = false;
		bool GGMRF = false;
		bool maskFP = false;
		bool maskBP = false;
		bool atten = false;
		bool attenM = false;
		bool norm = false;
		bool extra = false;
		bool raw = false;
		bool subInd = false;
		bool priorMask = false;
		int NLMRef = 0;
		bool BPIm = false;
		bool proj5Im = false;
		bool auxMod = false;
		bool FPMod = false;
		bool BPMod = false;
		bool SensMod = false;
		bool xFull = false;
		bool zFull = false;
		bool offsetT = false;
		bool indexBased = false;
		bool TOFIndex = false;
		bool angle = false;
		bool rayShifts = false;
		int zType = -1;
		int xSteps = -1;
		int zSteps = -1;
		int nSteps = 0;
		int eSteps = 0;
		int lSteps = 0;
		int iSteps = 0;
		int TOFSteps = 0;
		int aSteps = 0;
		int oSteps = 0;
	};
	CUDAMemAlloc memAlloc;
	bool useBuffers = true;

	template <typename K, typename T>
	inline K make_vec3(T a, T b, T c) {
		K apu;
		apu.x = a;
		apu.y = b;
		apu.z = c;
		return apu;
	}

	/// <summary>
	/// This function creates the CUDA programs for the forward and backward projections and for NLM/MRP/RDP/TV
	/// </summary>
	/// <param name="programFP the program to store forward projection program"></param>
	/// <param name="programBP the program to store backprojection program"></param>
	/// <param name="programAux the program to store auxliary (such as priors) programs"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="w_vec specifies some of the special options used"></param>
	/// <param name="local_size the local size"></param>
	/// <returns></returns>
	inline nvrtcResult createProgram(CUmodule& programFP, CUmodule& programBP,
		CUmodule& programAux, const char* header_directory, scalarStruct& inputScalars, const RecMethods MethodList,
		const Weighting& w_vec, const size_t local_size[], const int type = -1) {

		int compMajor = 0, compMinor = 0;
		cuDeviceGetAttribute(&compMajor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, CUDeviceID[0]);
		cuDeviceGetAttribute(&compMinor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, CUDeviceID[0]);

		nvrtcResult status = NVRTC_SUCCESS;


		std::string kernelFile = header_directory;
		std::string kernel_path, kernel_pathBP;
		std::string contentFP, contentBP;
		std::string contentAux;
		std::vector<const char*> options;
		int uu = 0;
		char buffer0[35];
		char buffer1[30];
		char buffer2[30];
		char buffer3[30];
		char buffer4[30];
		char buffer5[30];
		char buffer6[30];
		char buffer7[30];
		char buffer8[30];
		char buffer9[30];
		char buffer10[30];
		char buffer11[30];
		//char spectBuffer1[30];
		//char spectBuffer2[30];
		//char spectBuffer3[30];
		//char spectBuffer4[30];
		//char spectBuffer5[30];
		//char spectBuffer6[30];
		//char spectBuffer7[30];
		//char spectBuffer8[30];

		std::snprintf(buffer0, 35, "--gpu-architecture=compute_%d%d", compMajor, compMinor);
		options.push_back(buffer0);
		options.push_back("-DCUDA");
		if (inputScalars.useMAD) {
			options.push_back("--use_fast_math");
			options.push_back("-DUSEMAD");
		}
		if ((inputScalars.useImages && inputScalars.FPType != 4 && inputScalars.FPType != 5 && inputScalars.BPType != 5) || (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 5)) {
			options.push_back("-DUSEIMAGES");
			inputScalars.useBuffers = false;
			useBuffers = false;
		}
		std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
		// Load the header text file
		std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
		// Load orthogonal/volume of intersection headers if applicable
		if (inputScalars.FPType == 2 || inputScalars.BPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 3) {
			if (inputScalars.orthXY)
				options.push_back("-DCRYSTXY");
			if (inputScalars.orthZ)
				options.push_back("-DCRYSTZ");
			//std::ifstream sourceHeader1(kernelFile + "general_orth_opencl_functions.h");
			//std::string contentHeader1((std::istreambuf_iterator<char>(sourceHeader1)), std::istreambuf_iterator<char>());
			std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
			std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
			//contentHeader += contentHeader1 + contentHeader3;
			contentHeader += contentHeader3;
		}

		kernel_path = kernelFile;
		kernel_pathBP = kernelFile;
		if (inputScalars.FPType > 0) {
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				//if (!inputScalars.precompute && (inputScalars.n_rays * inputScalars.n_rays3D) > 1)
				//	kernel_path += "multidevice_siddon_no_precomp.cu");
				//else
				kernel_path += "projectorType123.cl";
			}
			else if (inputScalars.FPType == 4)
				kernel_path += "projectorType4.cl";
			else if (inputScalars.FPType == 5)
				kernel_path += "projectorType5.cl";
			std::ifstream sourceFile(kernel_path.c_str());
			std::string contentFFP((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
			contentFP = contentHeader + contentFFP;
		}
		if (inputScalars.BPType > 0) {
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				//if (!inputScalars.precompute && (inputScalars.n_rays * inputScalars.n_rays3D) > 1)
				//	kernel_pathBP += "multidevice_siddon_no_precomp.cu");
				//else
				kernel_pathBP += "projectorType123.cl";
			}
			else if (inputScalars.BPType == 4)
				kernel_pathBP += "projectorType4.cl";
			else if (inputScalars.BPType == 5)
				kernel_pathBP += "projectorType5.cl";
			std::ifstream sourceFileBP(kernel_pathBP.c_str());
			std::string contentFBP((std::istreambuf_iterator<char>(sourceFileBP)), std::istreambuf_iterator<char>());
			contentBP = contentHeader + contentFBP;
		}

		// Load the source text file
		// Set all preprocessor definitions
		const bool siddonVal = (inputScalars.FPType == 1 || inputScalars.BPType == 1 || inputScalars.FPType == 4 || inputScalars.BPType == 4) ? true : false;
		//if (inputScalars.FPType == 3 || inputScalars.BPType == 3)
		//	options.push_back("-DVOL");
		if (inputScalars.raw == 1)
			options.push_back("-DRAW");
		if (inputScalars.maskFP) {
			options.push_back("-DMASKFP");
			if (inputScalars.maskFPZ > 1)
				options.push_back("-DMASKFP3D");
		}
		if (inputScalars.useTotLength)
			options.push_back("-DTOTLENGTH");
		if (inputScalars.maskBP) {
			options.push_back("-DMASKBP");
			if (inputScalars.maskBPZ > 1)
				options.push_back("-DMASKBP3D");
		}
		if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights)
			options.push_back("-DFDK");
		if (inputScalars.offset)
			options.push_back("-DOFFSET");
		//if (inputScalars.FPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 2 || inputScalars.BPType == 3)
		//	options.push_back("-DORTH");
		if (inputScalars.attenuation_correction == 1u && inputScalars.CTAttenuation)
			options.push_back("-DATN");
		else if (inputScalars.attenuation_correction == 1u && !inputScalars.CTAttenuation)
			options.push_back("-DATNM");
		if (inputScalars.normalization_correction == 1u)
			options.push_back("-DNORM");
		if (inputScalars.scatter == 1u)
			options.push_back("-DSCATTER");
		if (inputScalars.randoms_correction == 1u)
			options.push_back("-DRANDOMS");
		if (inputScalars.nLayers > 1U) {
			if (inputScalars.listmode > 0 && inputScalars.indexBased)
				std::snprintf(buffer11, 30, "-DNLAYERS=%d", static_cast<int32_t>(inputScalars.nLayers));
			else
				std::snprintf(buffer11, 30, "-DNLAYERS=%d", static_cast<int32_t>(inputScalars.nProjections / (inputScalars.nLayers * inputScalars.nLayers)));
			options.push_back(buffer11);
		}
		if (inputScalars.TOF) {
			options.push_back("-DTOF");
		}
		if (inputScalars.CT)
			options.push_back("-DCT");
		else if (inputScalars.PET)
			options.push_back("-DPET");
		else if (inputScalars.SPECT) {
			options.push_back("-DSPECT");
			std::snprintf(buffer2, 30, "-DN_RAYS=%d", static_cast<int32_t>(inputScalars.n_rays * inputScalars.n_rays3D));
			options.push_back(buffer2);
			std::snprintf(buffer3, 30, "-DN_RAYS2D=%d", static_cast<int32_t>(inputScalars.n_rays));
			options.push_back(buffer3);
			std::snprintf(buffer4, 30, "-DN_RAYS3D=%d", static_cast<int32_t>(inputScalars.n_rays3D));
			options.push_back(buffer4);
			/*std::snprintf(spectBuffer1, 30, "-DCOL_D=%f", inputScalars.colD);
			options.push_back(spectBuffer1);
			std::snprintf(spectBuffer2, 30, "-DCOL_L=%f", inputScalars.colL);
			options.push_back(spectBuffer2);
			std::snprintf(spectBuffer3, 30, "-DDSEPTAL=%f", inputScalars.dSeptal);
			options.push_back(spectBuffer3);
			std::snprintf(spectBuffer4, 30, "-DHEXORIENTATION=%u", static_cast<uint8_t>(inputScalars.hexOrientation));
			options.push_back(spectBuffer4);
			std::snprintf(spectBuffer5, 30, "-DCONEMETHOD=%u", static_cast<uint8_t>(inputScalars.coneMethod));
			options.push_back(spectBuffer5);

			if (inputScalars.coneMethod == 3) {
				inputScalars.nRaySPECT = std::pow(std::ceil(std::sqrt(inputScalars.nRaySPECT)), 2);
			}
			std::snprintf(spectBuffer6, 30, "-DNRAYSPECT=%u", static_cast<uint16_t>(inputScalars.nRaySPECT));
			options.push_back(spectBuffer6);


			uint32_t nHexSPECT;
			if (inputScalars.coneMethod != 1) {
				std::snprintf(spectBuffer7, 30, "-DN_RAYS=%u", static_cast<uint16_t>(inputScalars.nRaySPECT));
				options.push_back(spectBuffer7);
				options.push_back("-DN_RAYS2D=1");
				options.push_back("-DN_RAYS3D=1");
				nHexSPECT = 1;
			} else {
				options.push_back("-DN_RAYS=1");
				options.push_back("-DN_RAYS2D=1");
				options.push_back("-DN_RAYS3D=1");
				nHexSPECT = std::pow(std::ceil(w_vec.dPitchX / inputScalars.colD), 2);
			}

			std::snprintf(spectBuffer8, 30, "-DNHEXSPECT=%u", static_cast<uint16_t>(nHexSPECT));
			options.push_back(spectBuffer8);*/
		}

		std::snprintf(buffer1, 30, "-DNBINS=%d", static_cast<int32_t>(inputScalars.nBins));
		options.push_back(buffer1);
		if (inputScalars.listmode == 1)
			options.push_back("-DLISTMODE");
		else if (inputScalars.listmode == 2)
			options.push_back("-DLISTMODE2");
		if (inputScalars.listmode > 0 && inputScalars.indexBased)
			options.push_back("-DINDEXBASED");
		if (siddonVal && (inputScalars.n_rays * inputScalars.n_rays3D) > 1) {
			std::snprintf(buffer2, 30, "-DN_RAYS=%d", static_cast<int32_t>(inputScalars.n_rays * inputScalars.n_rays3D));
			options.push_back(buffer2);
			std::snprintf(buffer3, 30, "-DN_RAYS2D=%d", static_cast<int32_t>(inputScalars.n_rays));
			options.push_back(buffer3);
			std::snprintf(buffer4, 30, "-DN_RAYS3D=%d", static_cast<int32_t>(inputScalars.n_rays3D));
			options.push_back(buffer4);
		}
		if (inputScalars.pitch)
			options.push_back("-DPITCH");
		if (((inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
			options.push_back("-DSUBSETS");
		if (local_size[1] > 0ULL) {
			std::snprintf(buffer5, 30, "-DLOCAL_SIZE=%d", static_cast<int32_t>(local_size[0]));
			options.push_back(buffer5);
			std::snprintf(buffer6, 30, "-DLOCAL_SIZE2=%d", static_cast<int32_t>(local_size[1]));
			options.push_back(buffer6);
		}
		else {
			std::snprintf(buffer5, 30, "-DLOCAL_SIZE=%d", static_cast<int32_t>(local_size[0]));
			options.push_back(buffer5);
		}
		if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
			std::snprintf(buffer7, 30, "-DSTYPE=%d", static_cast<int32_t>(inputScalars.subsetType));
			options.push_back(buffer7);
			std::snprintf(buffer8, 30, "-DNSUBSETS=%d", static_cast<int32_t>(inputScalars.subsets));
			options.push_back(buffer8);
		}
		if (DEBUG) {
			mexPrintBase("path = %s\n", kernel_path.c_str());
			mexPrintBase("pathBP = %s\n", kernel_pathBP.c_str());
			mexPrintBase("file = %s\n", kernelFile.c_str());
			//mexPrintBase("contentFP = %s\n", contentFP.c_str());
			mexPrintBase("inputScalars.BPType = %u\n", inputScalars.BPType);
			mexPrintBase("inputScalars.FPType = %u\n", inputScalars.FPType);
			mexEval();
		}
		// Build projector program
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			std::vector<const char*> os_options = options;
			os_options.push_back("-DAF");
			//if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.TOF) && inputScalars.dec > 0)
			//	os_options += (" -DDEC=" + std::to_string(inputScalars.dec));
			//os_options += (" -DN_REKOS=" + std::to_string(inputScalars.nRekos));
			//if (inputScalars.nRekos == 1)
			//	os_options.push_back("-DNREKOS1");
			//else if (inputScalars.nRekos == 2)
			//	os_options.push_back("-DNREKOS2");
			os_options.push_back("-DSIDDON");
			os_options.push_back("-DATOMICF");
			std::vector<const char*> os_optionsFP = os_options;
			os_optionsFP.push_back("-DFP");
			if (inputScalars.FPType == 3)
				os_optionsFP.push_back("-DVOL");
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3)
				os_optionsFP.push_back("-DORTH");
			//if (inputScalars.projector_type < 4) {
			//	os_optionsFP.push_back("-DBP");
			//	if (MethodList.MRAMLA || MethodList.MBSREM)
			//		os_optionsFP.push_back("-DMRAMLA");
			//	if (MethodList.COSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0u || MethodList.ECOSEM)
			//		os_optionsFP.push_back("-DCOSEM");
			//}

			//if (inputScalars.projector_type < 4)
			//	status = buildProgram(inputScalars.verbose, contentBP, CLContext, CUDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
			//else {
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build FP 1-3 program\n");
				}
				status = buildProgram(inputScalars.verbose, contentFP, programFP, os_optionsFP);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("FP 1-3 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.FPMod = true;
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build BP 1-3 program\n");
				}
				os_options.push_back("-DBP");
				if (inputScalars.BPType == 3)
					os_options.push_back("-DVOL");
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
					os_options.push_back("-DORTH");
				status = buildProgram(inputScalars.verbose, contentBP, programBP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("BP 1-3 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.BPMod = true;
			}
		}
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
			std::vector<const char*> os_options = options;
			if (inputScalars.FPType == 4)
				os_options.push_back("-DFP");
			if (inputScalars.BPType == 4 && inputScalars.CT)
				os_options.push_back("-DBP");
			//if (inputScalars.subsets > 1)
			//	os_options += (" -DSTYPE=" + std::to_string(inputScalars.subsetType));
			//if (inputScalars.projector_type == 41)
			//	os_options.push_back("-DPTYPE41");
			os_options.push_back("-DPTYPE4");
			if (!inputScalars.largeDim) {
				std::snprintf(buffer9, 30, "-DNVOXELS=%d", static_cast<int32_t>(NVOXELS));
				os_options.push_back(buffer9);
			}
			if (inputScalars.FPType == 4) {
				status = buildProgram(inputScalars.verbose, contentFP, programFP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("FP 4 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.FPMod = true;
			}
			if (!inputScalars.CT && inputScalars.BPType == 4) {
				os_options = options;
				os_options.push_back("-DPTYPE4");
				os_options.push_back("-DBP");
				os_options.push_back("-DATOMICF");
				status = buildProgram(inputScalars.verbose, contentBP, programBP,os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("BP 4 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.BPMod = true;
			}
			else if (inputScalars.CT && inputScalars.BPType == 4 && inputScalars.FPType != 4) {
				status = buildProgram(inputScalars.verbose, contentBP, programBP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("BP 4 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.BPMod = true;
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
			std::vector<const char*> os_options = options;
			os_options.push_back("-DPROJ5");
			if (inputScalars.meanFP)
				os_options.push_back("-DMEANDISTANCEFP");
			else if (inputScalars.meanBP)
				os_options.push_back("-DMEANDISTANCEBP");
			if (inputScalars.FPType == 5)
				os_options.push_back("-DFP");
			if (inputScalars.BPType == 5)
				os_options.push_back("-DBP");
			if (inputScalars.pitch) {
				std::snprintf(buffer9, 30, "-DNVOXELS5=%d", static_cast<int32_t>(1));
				os_options.push_back(buffer9);
			}
			else {
				std::snprintf(buffer9, 30, "-DNVOXELS5=%d", static_cast<int32_t>(NVOXELS5));
				os_options.push_back(buffer9);
			}
			std::snprintf(buffer10, 30, "-DNVOXELSFP=%d", static_cast<int32_t>(NVOXELSFP));
			os_options.push_back(buffer10);
			if (inputScalars.FPType == 5) {
				status = buildProgram(inputScalars.verbose, contentFP, programFP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("FP 5 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.FPMod = true;
			}
			else {
				status = buildProgram(inputScalars.verbose, contentBP, programBP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("BP 5 program built\n");
				}
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.BPMod = true;
			}
		}
		if (inputScalars.computeSensImag) {
			std::vector<const char*> os_options = options;
			os_options.push_back("-DBP");
			os_options.push_back("-DATOMICF");
			os_options.push_back("-DSENS");
			if (inputScalars.BPType == 4) {
				os_options.push_back("-DPTYPE4");
				os_options.push_back(buffer9);
			}
			else
				os_options.push_back("-DSIDDON");
			status = buildProgram(inputScalars.verbose, contentBP, programSens, os_options);
			if (status == NVRTC_SUCCESS)
				memAlloc.SensMod = true;
		}
		// Build prior programs
		if (MethodList.NLM || MethodList.MRP || MethodList.RDP || w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]
			|| MethodList.TV || MethodList.APLS || MethodList.hyperbolic || MethodList.ProxTV || MethodList.ProxTGV || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM ||
			MethodList.CPType || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF || type == 0) {
			if (DEBUG) {
				mexPrint("Building aux programs\n");
			}
			std::vector<const char*> optionsAux;
			int uu = 0;
			char buffer11[30];
			char buffer12[30];
			char buffer13[30];
			char buffer14[30];
			char buffer15[30];
			char buffer16[30];
			char buffer17[30];
			char buffer18[30];
			char buffer19[30];
			char buffer20[30];
			optionsAux.push_back(buffer0);
			optionsAux.push_back("-DCUDA");
			if (inputScalars.useMAD) {
				optionsAux.push_back("--use_fast_math");
				optionsAux.push_back("-DUSEMAD");
			}
			std::string auxKernelPath = kernelFile + "auxKernels.cl";
			std::ifstream sourceFileAux(auxKernelPath.c_str());
			std::string contentAAux((std::istreambuf_iterator<char>(sourceFileAux)), std::istreambuf_iterator<char>());
			contentAux = contentHeader + contentAAux;
			if (inputScalars.useImages)
				optionsAux.push_back("-DUSEIMAGES");
			if (inputScalars.useExtendedFOV)
				optionsAux.push_back("-DEFOV");
			if (inputScalars.use64BitIndices) {
				optionsAux.push_back("-DLTYPE=long long");
				optionsAux.push_back("-DLTYPE3=long3");
			}
			if (type != 0)
				optionsAux.push_back("-DAF");
			else {
				if (inputScalars.CT)
					optionsAux.push_back("-DCT");
				if (inputScalars.randoms_correction)
					optionsAux.push_back("-DRANDOMS");
				if (inputScalars.use_psf)
					optionsAux.push_back("-DPSF");
			}
			if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
				optionsAux.push_back("-DMASKPRIOR");
			if (inputScalars.eFOV)
				optionsAux.push_back("-DEFOVZ");
			if (MethodList.MRP) {
				optionsAux.push_back("-DMEDIAN");
				std::snprintf(buffer11, 30, "-DSEARCH_WINDOW_X=%d", static_cast<int32_t>(w_vec.Ndx));
				optionsAux.push_back(buffer11);
				std::snprintf(buffer12, 30, "-DSEARCH_WINDOW_Y=%d", static_cast<int32_t>(w_vec.Ndy));
				optionsAux.push_back(buffer12);
				std::snprintf(buffer13, 30, "-DSEARCH_WINDOW_Z=%d", static_cast<int32_t>(w_vec.Ndz));
				optionsAux.push_back(buffer13);
			}
			if (MethodList.NLM || MethodList.ProxNLM) {
				//optionsAux.push_back("-DPROXNLM");
				optionsAux.push_back("-DNLM_");
				if (w_vec.NLM_MRP) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(2));
					optionsAux.push_back(buffer11);
				}
				else if (w_vec.NLTV) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(1));
					optionsAux.push_back(buffer11);
				}
				else if (w_vec.NLRD) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(3));
					optionsAux.push_back(buffer11);
				}
				else if (w_vec.NLLange) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(4));
					optionsAux.push_back(buffer11);
				}
				else if (w_vec.NLLangeFiltered) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(5));
					optionsAux.push_back(buffer11);
				}
				else if (w_vec.NLGGMRF) {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(6));
					optionsAux.push_back(buffer11);
				}
				else {
					std::snprintf(buffer11, 30, "-DTYPE=%d", static_cast<int32_t>(0));
					optionsAux.push_back(buffer11);
				}
				if (w_vec.NLAdaptive)
					optionsAux.push_back("-DNLMADAPTIVE");
				if (w_vec.NLM_anatomical)
					optionsAux.push_back("-DNLMREF");
				std::snprintf(buffer12, 30, "-DSWINDOWX=%d", static_cast<int32_t>(w_vec.Ndx));
				optionsAux.push_back(buffer12);
				std::snprintf(buffer13, 30, "-DSWINDOWY=%d", static_cast<int32_t>(w_vec.Ndy));
				optionsAux.push_back(buffer13);
				std::snprintf(buffer14, 30, "-DSWINDOWZ=%d", static_cast<int32_t>(w_vec.Ndz));
				optionsAux.push_back(buffer14);
				std::snprintf(buffer15, 30, "-DPWINDOWX=%d", static_cast<int32_t>(w_vec.Nlx));
				optionsAux.push_back(buffer15);
				std::snprintf(buffer16, 30, "-DPWINDOWY=%d", static_cast<int32_t>(w_vec.Nly));
				optionsAux.push_back(buffer16);
				std::snprintf(buffer17, 30, "-DPWINDOWZ=%d", static_cast<int32_t>(w_vec.Nlz));
				optionsAux.push_back(buffer17);
			}
			if (MethodList.GGMRF) {
				optionsAux.push_back("-DGGMRF");
				std::snprintf(buffer11, 30, "-DSWINDOWX=%d", static_cast<int32_t>(w_vec.Ndx));
				optionsAux.push_back(buffer11);
				std::snprintf(buffer12, 30, "-DSWINDOWY=%d", static_cast<int32_t>(w_vec.Ndy));
				optionsAux.push_back(buffer12);
				std::snprintf(buffer13, 30, "-DSWINDOWZ=%d", static_cast<int32_t>(w_vec.Ndz));
				optionsAux.push_back(buffer13);
			}
			if (MethodList.hyperbolic) {
				optionsAux.push_back("-DHYPER");
				std::snprintf(buffer11, 30, "-DSWINDOWX=%d", static_cast<int32_t>(w_vec.Ndx));
				optionsAux.push_back(buffer11);
				std::snprintf(buffer12, 30, "-DSWINDOWY=%d", static_cast<int32_t>(w_vec.Ndy));
				optionsAux.push_back(buffer12);
				std::snprintf(buffer13, 30, "-DSWINDOWZ=%d", static_cast<int32_t>(w_vec.Ndz));
				optionsAux.push_back(buffer13);
			}
			if (MethodList.RDP) {
				optionsAux.push_back("-DRDP");
				if (w_vec.RDPLargeNeighbor) {
					optionsAux.push_back("-DRDPCORNERS");
					std::snprintf(buffer11, 30, "-DSWINDOWX=%d", static_cast<int32_t>(w_vec.Ndx));
					optionsAux.push_back(buffer11);
					std::snprintf(buffer12, 30, "-DSWINDOWY=%d", static_cast<int32_t>(w_vec.Ndy));
					optionsAux.push_back(buffer12);
					std::snprintf(buffer13, 30, "-DSWINDOWZ=%d", static_cast<int32_t>(w_vec.Ndz));
					optionsAux.push_back(buffer13);
				}
				if (w_vec.RDP_anatomical)
					optionsAux.push_back("-DRDPREF");
			}
			if (MethodList.ProxRDP && w_vec.RDPLargeNeighbor)
				optionsAux.push_back("-DRDPCORNERS");
			if (MethodList.TV && !w_vec.data.TV_use_anatomical) {
				optionsAux.push_back("-DTVGRAD");
				if (w_vec.data.TVtype == 6)
					optionsAux.push_back("-DTVW1");
				else if (w_vec.data.TVtype == 4)
					optionsAux.push_back("-DSATV");
				else if (w_vec.data.TVtype == 3)
					optionsAux.push_back("-DJPTV");
				if (w_vec.derivType > 0) {
					std::snprintf(buffer11, 30, "-DDIFFTYPE=%d", static_cast<int32_t>(w_vec.derivType));
					optionsAux.push_back(buffer11);
				}
			}
			else if ((MethodList.TV && w_vec.data.TV_use_anatomical) || MethodList.APLS) {
				optionsAux.push_back("-DTVGRAD");
				if (w_vec.data.TVtype == 1)
					optionsAux.push_back("-DANATOMICAL1");
				else if (w_vec.data.TVtype == 2)
					optionsAux.push_back("-DANATOMICAL2");
				else if (w_vec.data.TVtype == 5)
					optionsAux.push_back("-DANATOMICAL3");
				if (w_vec.derivType > 0) {
					std::snprintf(buffer11, 30, "-DDIFFTYPE=%d", static_cast<int32_t>(w_vec.derivType));
					optionsAux.push_back(buffer11);
				}
			}
			if (MethodList.ProxTV) {
				optionsAux.push_back("-DPROXTV");
				if (w_vec.UseL2Ball)
					optionsAux.push_back("-DL2");
				if (w_vec.derivType > 0) {
					std::snprintf(buffer11, 30, "-DDIFFTYPE=%d", static_cast<int32_t>(w_vec.derivType));
					optionsAux.push_back(buffer11);
				}
			}
			if (MethodList.ProxTGV || MethodList.TGV) {
				optionsAux.push_back("-DPROXTV");
				optionsAux.push_back("-DPROXTGV");
				if (w_vec.UseL2Ball)
					optionsAux.push_back("-DL2");
				if (w_vec.derivType > 0) {
					std::snprintf(buffer11, 30, "-DDIFFTYPE=%d", static_cast<int32_t>(w_vec.derivType));
					optionsAux.push_back(buffer11);
				}
				if (!inputScalars.TGV2D)
					optionsAux.push_back("-DTGVZ");
			}
			if (MethodList.ProxRDP)
				optionsAux.push_back("-DPROXRDP");
			if (local_sizePrior[1] > 0ULL) {
				std::snprintf(buffer18, 30, "-DLOCAL_SIZE=%d", static_cast<int32_t>(local_sizePrior[0]));
				optionsAux.push_back(buffer18);
				std::snprintf(buffer19, 30, "-DLOCAL_SIZE2=%d", static_cast<int32_t>(local_sizePrior[1]));
				optionsAux.push_back(buffer19);
				std::snprintf(buffer20, 30, "-DLOCAL_SIZE3=%d", static_cast<int32_t>(local_sizePrior[2]));
				optionsAux.push_back(buffer20);
			}
			else {
				std::snprintf(buffer18, 30, "-DLOCAL_SIZE=%d", static_cast<int32_t>(local_sizePrior[0]));
				optionsAux.push_back(buffer18);
			}
			if (MethodList.PKMA)
				optionsAux.push_back("-DPKMA");
			else if (MethodList.MBSREM)
				optionsAux.push_back("-DMBSREM");
			else if (MethodList.BSREM)
				optionsAux.push_back("-DBSREM");
			else if (MethodList.CPType) {
				optionsAux.push_back("-DPDHG");
				if (inputScalars.subsets > 1)
					optionsAux.push_back("-DSUBSETS");
			}
			status = buildProgram(inputScalars.verbose, contentAux, programAux, optionsAux);
			if (status == NVRTC_SUCCESS && DEBUG) {
				mexPrint("Aux program built\n");
			}
			else if (status != NVRTC_SUCCESS)
				return status;
			if (status == NVRTC_SUCCESS)
				memAlloc.auxMod = true;
		}
		if (DEBUG) {
			mexPrintBase("status = %u\n", status);
			mexEval();
		}
		return status;
	}

	/// <summary>
	/// Builds the input CUDA program
	/// </summary>
	/// <param name="verbose the level of verbosity"></param>
	/// <param name="contentFP program code"></param>
	/// <param name="program the program where to store the built program"></param>
	/// <param name="options preprocessor values for the build"></param>
	/// <returns></returns>
	inline nvrtcResult buildProgram(const int8_t verbose, std::string& content, CUmodule& module, std::vector<const char*>& options) {
		nvrtcResult status = NVRTC_SUCCESS;
		CUresult status2 = CUDA_SUCCESS;
		nvrtcProgram program;
		if (DEBUG || verbose >= 3) {
			for (int ll = 0; ll < options.size(); ll++)
				mexPrintBase("%s ", options[ll]);
			mexPrintBase("%s\n", "");
		}
		//const char* sourceCode = new char[content.size()];
		//sourceCode = content.c_str();
		status = nvrtcCreateProgram(&program, content.c_str(), "32bit", 0, NULL, NULL);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		//mexPrintBase("%s\n", content.c_str());

		// Build the program
		status = nvrtcCompileProgram(program, options.size(), options.data());
		// Build log in case of failure
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrint("Failed to build CUDA program. Build log: \n");
			size_t len;
			char* buffer;
			nvrtcGetProgramLogSize(program, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			nvrtcGetProgramLog(program, buffer);
			mexPrintBase("%s\n", buffer);
			free(buffer);
			nvrtcDestroyProgram(&program);
			return status;
		}
		else if (verbose > 1)
			mexPrint("CUDA program built\n");
		size_t ptxSize;
		status = nvrtcGetPTXSize(program, &ptxSize);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		char* ptx = new char[ptxSize];
		status = nvrtcGetPTX(program, ptx);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		status2 = cuModuleLoadDataEx(&module, ptx, 0, 0, 0);
		if (status2 != CUDA_SUCCESS) {
			getErrorString(status2);
			return NVRTC_ERROR_BUILTIN_OPERATION_FAILURE;
		}
		if (DEBUG) {
			mexPrintBase("ptxSize = %u\n", ptxSize);
		}
		// Destroy the program.
		status = nvrtcDestroyProgram(&program);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		//delete[] sourceCode;
		delete[] ptx;
		return status;
	}

	/// <summary>
	/// Creates the necessary CUDA kernels from the input programs
	/// </summary>
	/// <param name="kernelFP forward projection kernel"></param>
	/// <param name="kernelBP backprojection kernel"></param>
	/// <param name="kernelNLM NLM kernel"></param>
	/// <param name="kernelMed MRP kernel"></param>
	/// <param name="kernelRDP RDP kernel"></param>
	/// <param name="programFP program containing forward projection"></param>
	/// <param name="programBP program containing backprojection"></param>
	/// <param name="programAux program containing NLM/MRP/RDP"></param>
	/// <param name="MethodList reconstruction algorithms selected"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <returns></returns>
	inline CUresult createKernels(CUfunction& kernelFP, CUfunction& kernelBP, CUfunction& kernelNLM, CUfunction& kernelMed,
		CUfunction& kernelRDP, CUfunction& kernelGGMRF, const CUmodule& programFP, const CUmodule& programBP, const CUmodule& programAux,
		const RecMethods& MethodList, const Weighting& w_vec, const scalarStruct& inputScalars, const int type = -1) {
		CUresult status = CUDA_SUCCESS;
		// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
			if (inputScalars.FPType == 4) {
				status = cuModuleGetFunction(&kernelFP, programFP, "projectorType4Forward");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 4 FP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 4 FP successfully created\n");
				}
			}
			if (inputScalars.BPType == 4) {
				if (inputScalars.FPType == 4 && inputScalars.CT)
					status = cuModuleGetFunction(&kernelBP, programFP, "projectorType4Backward");
				else if (!inputScalars.CT)
					status = cuModuleGetFunction(&kernelBP, programBP, "projectorType4Forward");
				else
					status = cuModuleGetFunction(&kernelBP, programBP, "projectorType4Backward");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 4 BP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 4 BP successfully created\n");
				}
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
			if (inputScalars.FPType == 5) {
				status = cuModuleGetFunction(&kernelFP, programFP, "projectorType5Forward");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 5 FP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 5 FP successfully created\n");
				}
			}
			if (inputScalars.BPType == 5) {
				if (inputScalars.FPType == 5)
					status = cuModuleGetFunction(&kernelBP, programFP, "projectorType5Backward");
				else
					status = cuModuleGetFunction(&kernelBP, programBP, "projectorType5Backward");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 5 BP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 5 BP successfully created\n");
				}
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
				status = cuModuleGetFunction(&kernelFP, programFP, "projectorType123");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 1-3 FP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 1-3 FP successfully created\n");
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				status = cuModuleGetFunction(&kernelBP, programBP, "projectorType123");
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 1-3 BP kernel\n");
					return status;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("CUDA kernel for projector type 1-3 BP successfully created\n");
				}
			}

		}

		if (MethodList.NLM) {
			status = cuModuleGetFunction(&kernelNLM, programAux, "NLM");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create NLM kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("NLM kernel successfully created\n");
			}
		}
		if (MethodList.MRP) {
			status = cuModuleGetFunction(&kernelMed, programAux, "medianFilter3D");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create Median kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Median kernel successfully created\n");
			}
		}
		if (MethodList.RDP) {
			status = cuModuleGetFunction(&kernelRDP, programAux, "RDPKernel");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create RDP kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("RDP kernel successfully created\n");
			}
		}
		if (MethodList.GGMRF) {
			status = cuModuleGetFunction(&kernelGGMRF, programAux, "GGMRFKernel");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create GGMRF kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("GGMRF kernel successfully created\n");
			}
		}
		if (MethodList.TV) {
			status = cuModuleGetFunction(&kernelTV, programAux, "TVKernel");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create TV kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("TV kernel successfully created\n");
			}
		}
		if (MethodList.hyperbolic) {
			status = cuModuleGetFunction(&kernelHyper, programAux, "hyperbolicKernel");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create hyperbolic prior kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Hyperbolic prior kernel successfully created\n");
			}
		}
		if (MethodList.PKMA || MethodList.BSREM || MethodList.MBSREM) {
			status = cuModuleGetFunction(&kernelPoisson, programAux, "PoissonUpdate");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create Poisson Update kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Poisson Update kernel successfully created\n");
			}
		}
		if (MethodList.CPType) {
			status = cuModuleGetFunction(&kernelPDHG, programAux, "PDHGUpdate");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create PDHG Update kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("PDHG Update kernel successfully created\n");
			}
		}
		if (MethodList.ProxTV) {
			status = cuModuleGetFunction(&kernelProxTVq, programAux, "ProxTVq");
			status = cuModuleGetFunction(&kernelProxTVDiv, programAux, "ProxTVDivergence");
			status = cuModuleGetFunction(&kernelProxTVGrad, programAux, "ProxTVGradient");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create CPTV kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("CPTV kernel successfully created\n");
			}
		}
		if (MethodList.ProxRDP) {
			status = cuModuleGetFunction(&kernelProxq, programAux, "Proxq");
			status = cuModuleGetFunction(&kernelProxRDP, programAux, "ProxRDP");
			status = cuModuleGetFunction(&kernelProxTrans, programAux, "ProxTrans");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create proximal RDP kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Proximal RDP kernel successfully created\n");
			}
		}
		if (MethodList.ProxNLM) {
			status = cuModuleGetFunction(&kernelProxq, programAux, "Proxq");
			status = cuModuleGetFunction(&kernelProxNLM, programAux, "ProxNLM");
			status = cuModuleGetFunction(&kernelProxTrans, programAux, "ProxTrans");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create proximal NLM kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Proximal NLM kernel successfully created\n");
			}
		}
		if (MethodList.ProxTGV) {
			status = cuModuleGetFunction(&kernelProxTVq, programAux, "ProxTVq");
			status = cuModuleGetFunction(&kernelProxTGVq, programAux, "ProxTGVq");
			status = cuModuleGetFunction(&kernelProxTVDiv, programAux, "ProxTVDivergence");
			status = cuModuleGetFunction(&kernelProxTVGrad, programAux, "ProxTVGradient");
			status = cuModuleGetFunction(&kernelProxTGVDiv, programAux, "ProxTGVDivergence");
			status = cuModuleGetFunction(&kernelProxTGVSymmDeriv, programAux, "ProxTGVSymmDeriv");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create CPTGV kernel\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("CPTGV kernel successfully created\n");
			}
		}
		if (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) {
			status = cuModuleGetFunction(&kernelElementMultiply, programAux, "vectorElementMultiply");
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create element-wise kernels\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Element-wise kernels successfully created\n");
			}
			status = cuModuleGetFunction(&kernelElementDivision, programAux, "vectorElementDivision");

			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create element-wise kernels\n");
				return status;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Element-wise kernels successfully created\n");
			}
		}
		if (inputScalars.computeSensImag) {
			if (inputScalars.BPType == 4)
				status = cuModuleGetFunction(&kernelSensList, programSens, "projectorType4Forward");
			else
				status = cuModuleGetFunction(&kernelSensList, programSens, "projectorType123");
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create sensitivity image kernels\n");
				return status;
			}
		}
		//if (type == 0) {
		//	kernelsumma = cuModuleGetFunction(programAux, "summa", &status);
		//	if (status != CUDA_SUCCESS) {
		//		getErrorString(status);
		//		mexPrint("Failed to create implementation 3 kernels\n");
		//		return status;
		//	}
		//	kernelEstimate = cuModuleGetFunction(programAux, "computeEstimate", &status);
		//	if (status != CUDA_SUCCESS) {
		//		getErrorString(status);
		//		mexPrint("Failed to create implementation 3 kernels\n");
		//		return status;
		//	}
		//	kernelForward = cuModuleGetFunction(programAux, "forward", &status);
		//	if (inputScalars.use_psf)
		//		kernelPSFf = cuModuleGetFunction(programAux, "Convolution3D_f", &status);
		//	//kernelDiv = cuModuleGetFunction(programAux, "vectorDiv", &status);
		//	//kernelMult = cuModuleGetFunction(programAux, "vectorMult", &status);

		//	if (status != CUDA_SUCCESS) {
		//		getErrorString(status);
		//		mexPrint("Failed to create implementation 3 kernels\n");
		//		return status;
		//	}
		//	else if (DEBUG || inputScalars.verbose >= 2) {
		//		mexPrint("Implementation 3 kernels successfully created\n");
		//	}
		//}
		return status;
	}
public:
	//CUcontext CLContext;
	std::vector<CUdevice> CUDeviceID;
	std::vector<CUstream> CLCommandQueue;
	CUfunction kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, kernelElementMultiply, kernelElementDivision,
		kernelTV, kernelProxTGVSymmDeriv, kernelProxTGVDiv, kernelProxTGVq, kernelPoisson, kernelPDHG, kernelProxRDP, kernelProxq, kernelProxTrans, kernelProxNLM, kernelGGMRF,
		kernelsumma, kernelEstimate, kernelPSF, kernelPSFf, kernelDiv, kernelMult, kernelForward, kernelSensList, kernelApu, kernelHyper;
	CUmodule programFP, programBP, programAux, programSens;
	CUdeviceptr d_angle, d_xcenter, d_ycenter, d_zcenter, d_V, d_TOFCenter, *d_output, *d_meanBP, *d_meanFP, d_eFOVIndices, d_weights, *d_inputB, *d_W, *d_gaussianNLM;
	CUtexObject d_maskFP, d_maskBP, d_maskPrior;
	CUtexObject d_inputImage, d_imageX, d_imageY, d_attenIm, d_urefIm, d_inputI, d_RDPrefI;
	CUarray atArray, uRefArray, maskArrayBP, maskArrayPrior, BPArray, FPArray, integArrayXY, imArray;
	CUdeviceptr *d_qX, *d_qY, *d_qZ;
	CUdeviceptr *d_rX, *d_rY, *d_rXY, *d_rZ, *d_rXZ, *d_rYZ;
	CUdeviceptr *d_vX, *d_vY, *d_vZ;
	CUdeviceptr *d_vector, *d_input;
	CUdeviceptr* d_im, *d_rhs, *d_U, d_g, d_uref, *d_refIm, d_attenB, d_maskFPB, d_maskBPB, *d_RDPref;
	CUdeviceptr d_rayShiftsDetector, d_rayShiftsSource; // SPECT
	//CUdeviceptr d_outputCT;
	std::vector<void*> FPArgs, BPArgs, SensArgs;
	CUDA_im_vectors vec_opencl;
	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
	std::vector<float3> b, d, bmax;
	// Image dimensions
	std::vector<int3> d_N;

	std::vector<CUtexObject> d_maskFP3;
	std::vector<CUdeviceptr> d_maskFPB;
	std::vector<CUarray> maskArrayFP;
	std::vector<CUdeviceptr> d_normFull, d_scatFull, d_xFull, d_zFull;
	std::vector<CUdeviceptr> d_L;
	std::vector<CUdeviceptr*> d_Summ;
	std::vector<CUdeviceptr> d_zindex;
	std::vector<CUdeviceptr> d_xyindex;
	std::vector<CUdeviceptr> d_norm;
	std::vector<CUdeviceptr> d_scat;
	std::vector<CUdeviceptr> d_x;
	std::vector<CUdeviceptr> d_z;
	std::vector<CUdeviceptr> d_atten;
	std::vector<CUdeviceptr> d_T;
	std::vector<CUdeviceptr> d_trIndex;
	std::vector<CUdeviceptr> d_axIndex;
	std::vector<CUdeviceptr> d_TOFIndex;
	CUDA_TEXTURE_DESC texDesc;
	CUDA_ARRAY_DESCRIPTOR arr2DDesc;
	CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
	CUDA_RESOURCE_DESC resDesc;
	CUDA_RESOURCE_VIEW_DESC viewDesc;
	unsigned char no_norm = 0;
	int proj6 = 1;
	size_t memSize = 0ULL;
	std::vector<std::vector<size_t>> erotusBP, erotusPDHG;
	~ProjectorClass() {
		CUresult status = CUDA_SUCCESS;
		if (memAlloc.FPMod)
			getErrorString(cuModuleUnload(programFP));
		if (memAlloc.BPMod)
			getErrorString(cuModuleUnload(programBP));
		if (memAlloc.auxMod)
			getErrorString(cuModuleUnload(programAux));
		if (memAlloc.SensMod)
			getErrorString(cuModuleUnload(programSens));
		//if (memAlloc.xC)
		//	getErrorString(cuMemFree(d_xcenter));
		//if (memAlloc.yC)
		//	getErrorString(cuMemFree(d_ycenter));
		//if (memAlloc.zC)
		//	getErrorString(cuMemFree(d_zcenter));
		if (memAlloc.attenM) {
			for (int kk = 0; kk < memAlloc.aSteps; kk++) {
				getErrorString(cuMemFree(d_atten[kk]));
			}
		}
		if (memAlloc.V)
			getErrorString(cuMemFree(d_V));
		if (memAlloc.xSteps >= 0) {
			for (int kk = 0; kk <= memAlloc.xSteps; kk++) {
				getErrorString(cuMemFree(d_x[kk]));
			}
		}
		if (memAlloc.zType == 0) {
			getErrorString(cuMemFree(d_z[memAlloc.zSteps]));
		}
		else if (memAlloc.zType == 1) {
			for (int kk = 0; kk <= memAlloc.zSteps; kk++) {
				getErrorString(cuMemFree(d_z[kk]));
			}
		}
		if (memAlloc.offsetT) {
			for (int kk = 0; kk < memAlloc.oSteps; kk++) {
				getErrorString(cuMemFree(d_T[kk]));
			}
		}
		if (memAlloc.TOF) {
			getErrorString(cuMemFree(d_TOFCenter));
		}
		if (memAlloc.eFOV) {
			getErrorString(cuMemFree(d_eFOVIndices));
		}
		if (memAlloc.rayShifts) {
			getErrorString(cuMemFree(d_rayShiftsDetector));
			getErrorString(cuMemFree(d_rayShiftsSource));
		}
		if (memAlloc.GGMRF) {
			getErrorString(cuMemFree(d_weights));
		}
		if (memAlloc.norm) {
			for (int kk = 0; kk < memAlloc.nSteps; kk++) {
				getErrorString(cuMemFree(d_norm[kk]));
			}
		}
		if (memAlloc.extra) {
			for (int kk = 0; kk < memAlloc.eSteps; kk++) {
				getErrorString(cuMemFree(d_scat[kk]));
			}
		}
		if (memAlloc.indexBased) {
			for (int kk = 0; kk < memAlloc.iSteps; kk++) {
				getErrorString(cuMemFree(d_trIndex[kk]));
				getErrorString(cuMemFree(d_axIndex[kk]));
			}
		}
		if (memAlloc.TOFIndex) {
			for (int kk = 0; kk < memAlloc.TOFSteps; kk++) {
				getErrorString(cuMemFree(d_TOFIndex[kk]));
			}
		}
		if (memAlloc.angle)
			getErrorString(cuMemFree(d_angle));
		if (memAlloc.xFull)
			getErrorString(cuMemFree(d_xFull[0]));
		if (memAlloc.zFull)
			getErrorString(cuMemFree(d_zFull[0]));
		if (memAlloc.maskFP) {
			if (useBuffers) {
				for (int ll = 0; ll < d_maskFPB.size(); ll++)
					getErrorString(cuMemFree(d_maskFPB[ll]));
			}
			else {
				if (d_maskFP3.size() > 0) {
					for (int ll = 0; ll < d_maskFP3.size(); ll++)
						getErrorString(cuTexObjectDestroy(d_maskFP3[ll]));
				}
				else {
					getErrorString(cuTexObjectDestroy(d_maskFP));
				}
				for (int ll = 0; ll < maskArrayFP.size(); ll++)
					getErrorString(cuArrayDestroy(maskArrayFP[ll]));
			}
		}
		if (memAlloc.maskBP) {
			if (useBuffers) {
				getErrorString(cuMemFree(d_maskBPB));
			}
			else {
				getErrorString(cuTexObjectDestroy(d_maskBP));
				getErrorString(cuArrayDestroy(maskArrayBP));
			}
		}
		if (memAlloc.atten) {
			if (useBuffers) {
				getErrorString(cuMemFree(d_attenB));
			}
			else {
				getErrorString(cuTexObjectDestroy(d_attenIm));
				getErrorString(cuArrayDestroy(atArray));
			}
		}
		if (memAlloc.priorMask) {
			getErrorString(cuTexObjectDestroy(d_maskPrior));
			getErrorString(cuArrayDestroy(maskArrayPrior));
		}
		if (memAlloc.NLMRef == 1) {
			getErrorString(cuTexObjectDestroy(d_urefIm));
			getErrorString(cuArrayDestroy(uRefArray));
		}
		else if (memAlloc.NLMRef == 2) {
			getErrorString(cuMemFree(d_uref));
		}
	}
	/// <summary>
	/// This function creates the projector class object
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <returns></returns>
	inline int addProjector(scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, const char* header_directory, const int type = -1) {
		// Set-up the local group size
		local_size[0] = 32U;
		local_size[1] = 1U;
		local_size[2] = 1U;
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || (inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT)))
			local_size[0] = 128ULL;
		if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
			if (inputScalars.nColsD > 1 && !(inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT))) {
				local_size[0] = 16U;
				local_size[1] = 16U;
			}
		}
		if (DEBUG) {
			mexPrintBase("inputScalars.nColsD = %u\n", inputScalars.nColsD);
			mexPrintBase("inputScalars.nRowsD = %u\n", inputScalars.nRowsD);
			mexPrintBase("local_size[0] = %u\n", local_size[0]);
			mexPrintBase("local_size[1] = %u\n", local_size[1]);
			mexEval();
		}
		// Local group for priors
		local_sizePrior[0] = 16U;
		local_sizePrior[1] = 16U;
		local_sizePrior[2] = 1U;
		CUresult status = CUDA_SUCCESS;
		nvrtcResult status2 = NVRTC_SUCCESS;
		proj6 = 0;


		// Create the CUDA context and command queue and assign the device
//#ifdef AF
		//CLContext = afcu::getContext(true);
		//std::vector<CUdevice> devices = CLContext.getInfo<CL_CONTEXT_DEVICES>(&status);
		//if (status != CUDA_SUCCESS) {
		//	getErrorString(status);
		//	return status;
		//}
		int af_id = af::getDevice();
		CUDeviceID.push_back(afcu::getNativeId(af_id));
		//CUstream testi = afcu::getStream(CUDeviceID[0]);
		CLCommandQueue.push_back(afcu::getStream(CUDeviceID[0]));
//#else
//		status = clGetPlatformsContext(inputScalars.platform, CLContext, CLCommandQueue, inputScalars.usedDevices, CUDeviceID);
//#endif
		// For NVIDIA cards, 32 local size seems more optimal with 1D kernelFP
		//std::string deviceName = CUDeviceID[0].getInfo<CL_DEVICE_VENDOR>(&status);
		//std::string NV("NVIDIA Corporation");
		//if (NV.compare(deviceName) == 0 && (inputScalars.projector_type == 1 || inputScalars.projector_type == 11) && local_size[1] == 1ULL)
		//	local_size[0] = 32ULL;
		//if (DEBUG) {
		//	unsigned long long apu = CUDeviceID[0].getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(&status);
		//	unsigned int apu2 = CUDeviceID[0].getInfo<CL_DEVICE_ADDRESS_BITS>(&status);
		//	mexPrintBase("CL_DEVICE_MAX_MEM_ALLOC_SIZE = %llu\n", apu);
		//	mexPrintBase("CL_DEVICE_ADDRESS_BITS = %u\n", apu2);
		//	mexEval();
		//}

		status2 = createProgram(programFP, programBP, programAux, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		if (status2 != NVRTC_SUCCESS) {
			std::cerr << "Error while creating program" << std::endl;
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 2) {
			mexPrint("CUDA programs successfully created\n");
		}

		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, MethodList, w_vec, inputScalars, type);
		if (status != CUDA_SUCCESS) {
			mexPrint("Failed to create kernels\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 2) {
			mexPrint("CUDA kernels successfully created\n");
		}
		//format.image_channel_order = CL_A;
		//format.image_channel_data_type = CL_FLOAT;
		//formatMask.image_channel_order = CL_A;
		//formatMask.image_channel_data_type = CL_UNSIGNED_INT8;

		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			erotus[0] = inputScalars.nRowsD % local_size[0];
			if (inputScalars.FPType == 5)
				erotus[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP) % local_size[1];
			else
				erotus[1] = inputScalars.nColsD % local_size[1];
			if (erotus[1] > 0)
				erotus[1] = (local_size[1] - erotus[1]);
			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
		}

		if ((MethodList.ProxTGV || MethodList.ProxTV || MethodList.ProxRDP)) {
			erotusPriorEFOV[0] = inputScalars.NxPrior % local_sizePrior[0];
			erotusPriorEFOV[1] = inputScalars.NyPrior % local_sizePrior[1];
			erotusPriorEFOV[2] = inputScalars.NzPrior % local_sizePrior[2];
			if (erotusPriorEFOV[0] > 0)
				erotusPriorEFOV[0] = (local_sizePrior[0] - erotusPriorEFOV[0]);
			if (erotusPriorEFOV[1] > 0)
				erotusPriorEFOV[1] = (local_sizePrior[1] - erotusPriorEFOV[1]);
			if (erotusPriorEFOV[2] > 0)
				erotusPriorEFOV[2] = (local_sizePrior[1] - erotusPriorEFOV[2]);
			globalPriorEFOV[0] = (inputScalars.NxPrior + erotusPriorEFOV[0]) / local_sizePrior[0];
			globalPriorEFOV[1] = (inputScalars.NyPrior + erotusPriorEFOV[1]) / local_sizePrior[1];
			globalPriorEFOV[2] = (inputScalars.NzPrior + erotusPriorEFOV[2]) / local_sizePrior[2];
		}

		erotusBP.resize(2);
		erotusPDHG.resize(2);
		if (MethodList.CPType || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM) {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				erotusPDHG[0].emplace_back(inputScalars.Nx[ii] % local_sizePrior[0]);
				erotusPDHG[1].emplace_back(inputScalars.Ny[ii] % local_sizePrior[1]);
				if (erotusPDHG[0][ii] > 0)
					erotusPDHG[0][ii] = (local_sizePrior[0] - erotusPDHG[0][ii]);
				if (erotusPDHG[1][ii] > 0)
					erotusPDHG[1][ii] = (local_sizePrior[1] - erotusPDHG[1][ii]);
			}
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			erotusBP[0].emplace_back(inputScalars.Nx[ii] % local_size[0]);
			erotusBP[1].emplace_back(inputScalars.Ny[ii] % local_size[1]);
			if (erotusBP[0][ii] > 0)
				erotusBP[0][ii] = (local_size[0] - erotusBP[0][ii]);
			if (erotusBP[1][ii] > 0)
				erotusBP[1][ii] = (local_size[1] - erotusBP[1][ii]);
		}
		local[0] = local_size[0];
		local[1] = local_size[1];
		local[2] = 1;
		localPrior[0] = local_sizePrior[0];
		localPrior[1] = local_sizePrior[1];
		localPrior[2] = local_sizePrior[2];
		erotusPrior[0] = inputScalars.Nx[0] % local_sizePrior[0];
		erotusPrior[1] = inputScalars.Ny[0] % local_sizePrior[1];
		erotusPrior[2] = inputScalars.Nz[0] % local_sizePrior[2];
		if (erotusPrior[0] > 0)
			erotusPrior[0] = (local_sizePrior[0] - erotusPrior[0]);
		if (erotusPrior[1] > 0)
			erotusPrior[1] = (local_sizePrior[1] - erotusPrior[1]);
		if (erotusPrior[2] > 0)
			erotusPrior[2] = (local_sizePrior[1] - erotusPrior[2]);
		globalPrior[0] = (inputScalars.Nx[0] + erotusPrior[0]) / localPrior[0];
		globalPrior[1] = (inputScalars.Ny[0] + erotusPrior[1]) / localPrior[1];
		globalPrior[2] = (inputScalars.Nz[0] + erotusPrior[2]) / localPrior[2];
		d_NOrig = make_vec3<int3>(static_cast<int>(inputScalars.NxOrig), static_cast<int>(inputScalars.NyOrig), static_cast<int>(inputScalars.NzOrig));
		d_NPrior = make_vec3<int3>(static_cast<int>(inputScalars.NxPrior), static_cast<int>(inputScalars.NyPrior), static_cast<int>(inputScalars.NzPrior));
		dPitch = { w_vec.dPitchX, w_vec.dPitchY };
		b.resize(inputScalars.nMultiVolumes + 1);
		d.resize(inputScalars.nMultiVolumes + 1);
		d_N.resize(inputScalars.nMultiVolumes + 1);
		bmax.resize(inputScalars.nMultiVolumes + 1);
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			b[ii] = make_vec3<float3>(inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii]);
			d[ii] = make_vec3<float3>(inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii]);
			d_N[ii] = make_vec3<int3>(static_cast<int>(inputScalars.Nx[ii]), static_cast<int>(inputScalars.Ny[ii]), static_cast<int>(inputScalars.Nz[ii]));
			bmax[ii] = make_vec3<float3>(static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
				static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii],
				static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii]);
		}
		if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
			erotusSens[0] = inputScalars.det_per_ring % local_size[0];
			erotusSens[1] = inputScalars.det_per_ring % local_size[1];
			if (erotusSens[1] > 0)
				erotusSens[1] = (local_size[1] - erotusSens[1]);
			if (erotusSens[0] > 0)
				erotusSens[0] = (local_size[0] - erotusSens[0]);
			d_xFull.resize(1);
			d_zFull.resize(1);
		}
		if (DEBUG)
			mexPrint("Luuppi valmis\n");
		//vec_opencl.d_rhs_os.resize(1);
		d_Summ.resize(1);
		d_Summ[0] = nullptr;
		return 0;
	}

	/// <summary>
	/// This function first creates the necessary CUDA buffers and then writes the data to them
	/// </summary>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="x the x/y/z coordinates for the detectors (PET and SPECT) or source and detector (CT). z-coordinate applies only for CT"></param>
	/// <param name="z_det the z coordinates for the detectors (PET and SPECT) or the directional vectors for the detector panel pixels (CT)"></param>
	/// <param name="xy_index subset indices for subsets types &lt; 8, x/y dimensions"></param>
	/// <param name="z_index same as above but for z dimension"></param>
	/// <param name="L raw data detector indices"></param>
	/// <param name="pituus cumulative sum of length"></param>
	/// <param name="atten attenuation image"></param>
	/// <param name="norm normalization matrix"></param>
	/// <param name="extraCorr scatter data (for multiplicative scatter correction)"></param>
	/// <param name="x_center x-coordinates of the voxel centers"></param>
	/// <param name="y_center y-coordinates of the voxel centers"></param>
	/// <param name="z_center z-coordinates of the voxel centers"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="Sin measurement data (sinograms or projections)"></param>
	/// <param name="sc_ra randoms and/or scatter data (for additive scatter correction or for randoms correction)"></param>
	/// <returns></returns>
	inline int createAndWriteBuffers(const std::vector<int64_t>& length, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* L, const int64_t* pituus, const float* atten, const float* norm, const float* extraCorr,
		const scalarStruct& inputScalars, const Weighting& w_vec, const RecMethods& MethodList) {
		CUresult status = CUDA_SUCCESS;
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		// Create the necessary buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			status = cuMemAlloc(&d_weights, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			memAlloc.GGMRF = true;
		}
		//else if (MethodList.hyperbolic) {
		//	status = cuMemAlloc(&d_weights, sizeof(float) * ((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1));
		//	if (status != CUDA_SUCCESS) {
		//		getErrorString(status);
		//		return -1;
		//	}
		//	memAlloc.GGMRF = true;
		//}
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (inputScalars.useImages) {
				std::memset(&texDesc, 0, sizeof(texDesc));
				std::memset(&resDesc, 0, sizeof(resDesc));
				std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
				std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
				std::memset(&viewDesc, 0, sizeof(viewDesc));
				arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
				arr3DDesc.NumChannels = 1;
				arr3DDesc.Height = inputScalars.Nx[0];
				arr3DDesc.Width = inputScalars.Ny[0];
				arr3DDesc.Depth = inputScalars.Nz[0];
				status = cuArray3DCreate(&uRefArray, &arr3DDesc);
				CUDA_MEMCPY3D cpy3d;
				std::memset(&cpy3d, 0, sizeof(cpy3d));
				cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
				cpy3d.srcHost = w_vec.NLM_ref;
				cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(float);
				cpy3d.srcHeight = inputScalars.Nx[0];
				cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
				cpy3d.dstArray = uRefArray;
				cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(float);
				cpy3d.Height = inputScalars.Nx[0];
				cpy3d.Depth = inputScalars.Nz[0];
				status = cuMemcpy3D(&cpy3d);
				resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
				resDesc.res.array.hArray = uRefArray;
				texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
				viewDesc.height = inputScalars.Nx[0];
				viewDesc.width = inputScalars.Ny[0];
				viewDesc.depth = inputScalars.Nz[0];
				viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
				status = cuTexObjectCreate(&d_urefIm, &resDesc, &texDesc, &viewDesc);
			}
			else
				status = cuMemAlloc(&d_uref, sizeof(float) * inputScalars.im_dim[0]);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.useImages)
				memAlloc.NLMRef = 1;
			else
				memAlloc.NLMRef = 2;
		}
		if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
			std::memset(&texDesc, 0, sizeof(texDesc));
			std::memset(&resDesc, 0, sizeof(resDesc));
			std::memset(&viewDesc, 0, sizeof(viewDesc));
			if (inputScalars.maskBPZ > 1) {
				std::memset(&arr2DDesc, 0, sizeof(arr3DDesc));
				arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
				arr3DDesc.NumChannels = 1;
				arr3DDesc.Height = inputScalars.Nx[0];
				arr3DDesc.Width = inputScalars.Ny[0];
				arr3DDesc.Depth = inputScalars.Nz[0];
				status = cuArray3DCreate(&maskArrayPrior, &arr3DDesc);
				CUDA_MEMCPY3D cpy3d;
				std::memset(&cpy3d, 0, sizeof(cpy3d));
				cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
				cpy3d.srcHost = w_vec.maskPrior;
				cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(uint8_t);
				cpy3d.srcHeight = inputScalars.Nx[0];
				cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
				cpy3d.dstArray = maskArrayPrior;
				cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(uint8_t);
				cpy3d.Height = inputScalars.Nx[0];
				cpy3d.Depth = inputScalars.maskBPZ;
				status = cuMemcpy3D(&cpy3d);
			}
			else {
				std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
				arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
				arr2DDesc.NumChannels = 1;
				arr2DDesc.Height = inputScalars.Nx[0];
				arr2DDesc.Width = inputScalars.Ny[0];
				status = cuArrayCreate(&maskArrayPrior, &arr2DDesc);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				CUDA_MEMCPY2D cpy2d;
				std::memset(&cpy2d, 0, sizeof(cpy2d));
				cpy2d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
				cpy2d.srcHost = w_vec.maskPrior;
				cpy2d.srcPitch = inputScalars.Ny[0] * sizeof(uint8_t);
				cpy2d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
				cpy2d.dstArray = maskArrayPrior;
				cpy2d.WidthInBytes = inputScalars.Ny[0] * sizeof(uint8_t);
				cpy2d.Height = inputScalars.Nx[0];
				status = cuMemcpy2D(&cpy2d);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
			resDesc.res.array.hArray = maskArrayPrior;
			texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			//texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
			texDesc.flags = CU_TRSF_READ_AS_INTEGER;
			viewDesc.height = inputScalars.Nx[0];
			viewDesc.width = inputScalars.Ny[0];
			if (inputScalars.maskBPZ > 1)
				viewDesc.depth = inputScalars.maskBPZ;
			viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
			status = cuTexObjectCreate(&d_maskPrior, &resDesc, &texDesc, &viewDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			memAlloc.priorMask = true;
		}
		if (inputScalars.projector_type != 6) {
			status = cuMemAlloc(&d_V, sizeof(float) * inputScalars.size_V);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			memAlloc.V = true;
			// Detector coordinates
			if ((!inputScalars.CT && inputScalars.listmode == 0) || inputScalars.indexBased) {
				status = cuMemAlloc(&d_x[0], sizeof(float) * inputScalars.size_of_x);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.xSteps++;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				status = cuMemAlloc(&d_xFull[0], sizeof(float) * inputScalars.size_of_x);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.xFull = true;
			}
			//status = cuMemAlloc(&d_xcenter, sizeof(float) * inputScalars.size_center_x);;
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			//memAlloc.xC = true;
			//status = cuMemAlloc(&d_ycenter, sizeof(float) * inputScalars.size_center_y);
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			//memAlloc.yC = true;
			//status = cuMemAlloc(&d_zcenter, sizeof(float) * inputScalars.size_center_z);
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			//memAlloc.zC = true;
			// Attenuation data for image-based attenuation
			if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.useBuffers)
					status = cuMemAlloc(&d_attenB, sizeof(float) * inputScalars.im_dim[0]);
				else {
					std::memset(&texDesc, 0, sizeof(texDesc));
					std::memset(&resDesc, 0, sizeof(resDesc));
					std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
					std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
					std::memset(&viewDesc, 0, sizeof(viewDesc));
					arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
					arr3DDesc.NumChannels = 1;
					arr3DDesc.Height = inputScalars.Nx[0];
					arr3DDesc.Width = inputScalars.Ny[0];
					arr3DDesc.Depth = inputScalars.Nz[0];
					status = cuArray3DCreate(&atArray, &arr3DDesc);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					CUDA_MEMCPY3D cpy3d;
					std::memset(&cpy3d, 0, sizeof(cpy3d));
					cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
					cpy3d.srcHost = atten;
					cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(float);
					cpy3d.srcHeight = inputScalars.Nx[0];
					cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
					cpy3d.dstArray = atArray;
					cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(float);
					cpy3d.Height = inputScalars.Nx[0];
					cpy3d.Depth = inputScalars.Nz[0];
					status = cuMemcpy3D(&cpy3d);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
					resDesc.res.array.hArray = atArray;
					texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
					viewDesc.height = inputScalars.Nx[0];
					viewDesc.width = inputScalars.Ny[0];
					viewDesc.depth = inputScalars.Nz[0];
					viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
					status = cuTexObjectCreate(&d_attenIm, &resDesc, &texDesc, &viewDesc);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				memAlloc.atten = true;
			}
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						d_maskFPB.resize(inputScalars.subsetsUsed);
						for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
							status = cuMemAlloc(&d_maskFPB[kk], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk]);
					}
					else {
						std::memset(&texDesc, 0, sizeof(texDesc));
						std::memset(&resDesc, 0, sizeof(resDesc));
						std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
						std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
						std::memset(&viewDesc, 0, sizeof(viewDesc));
						if (inputScalars.maskFPZ > 1) {
							maskArrayFP.resize(inputScalars.subsetsUsed);
							arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr3DDesc.NumChannels = 1;
							arr3DDesc.Height = inputScalars.nRowsD;
							arr3DDesc.Width = inputScalars.nColsD;
							CUDA_MEMCPY3D cpy3d;
							std::memset(&cpy3d, 0, sizeof(cpy3d));
							cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
							cpy3d.srcPitch = inputScalars.nColsD * sizeof(uint8_t);
							cpy3d.srcHeight = inputScalars.nRowsD;
							cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
							cpy3d.WidthInBytes = inputScalars.nColsD * sizeof(uint8_t);
							cpy3d.Height = inputScalars.nRowsD;
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
								arr3DDesc.Depth = length[kk];
								status = cuArray3DCreate(&maskArrayFP[kk], &arr3DDesc);
								if (status != CUDA_SUCCESS) {
									getErrorString(status);
									return -1;
								}
								cpy3d.Depth = length[kk];
								cpy3d.srcHost = &w_vec.maskFP[pituus[kk] * vecSize];
								cpy3d.dstArray = maskArrayFP[kk];
								status = cuMemcpy3D(&cpy3d);
								if (status != CUDA_SUCCESS) {
									getErrorString(status);
									return -1;
								}
							}
						}
						else {
							maskArrayFP.resize(1);
							arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr2DDesc.NumChannels = 1;
							arr2DDesc.Height = inputScalars.nRowsD;
							arr2DDesc.Width = inputScalars.nColsD;
							status = cuArrayCreate(&maskArrayFP[0], &arr2DDesc);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
							CUDA_MEMCPY2D cpy2d;
							std::memset(&cpy2d, 0, sizeof(cpy2d));
							cpy2d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
							cpy2d.srcHost = w_vec.maskFP;
							cpy2d.srcPitch = inputScalars.nColsD * sizeof(uint8_t);
							cpy2d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
							cpy2d.dstArray = maskArrayFP[0];
							cpy2d.WidthInBytes = inputScalars.nColsD * sizeof(uint8_t);
							cpy2d.Height = inputScalars.nRowsD;
							status = cuMemcpy2D(&cpy2d);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
						}
						resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						texDesc.flags = CU_TRSF_READ_AS_INTEGER;
						viewDesc.height = inputScalars.nRowsD;
						viewDesc.width = inputScalars.nColsD;
						//viewDesc.depth = 1;
						viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
						if (inputScalars.maskFPZ > 1) {
							texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							d_maskFP3.resize(inputScalars.subsetsUsed);
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
								viewDesc.depth = length[kk];
								resDesc.res.array.hArray = maskArrayFP[kk];
								status = cuTexObjectCreate(&d_maskFP3[kk], &resDesc, &texDesc, &viewDesc);
							}
						}
						else {
							resDesc.res.array.hArray = maskArrayFP[0];
							status = cuTexObjectCreate(&d_maskFP, &resDesc, &texDesc, &viewDesc);
						}
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
					memAlloc.maskFP = true;
				}
				if (inputScalars.maskBP) {
					if (inputScalars.useBuffers)
						status = cuMemAlloc(&d_maskBPB, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ);
					else {
						std::memset(&texDesc, 0, sizeof(texDesc));
						std::memset(&resDesc, 0, sizeof(resDesc));
						std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
						std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
						std::memset(&viewDesc, 0, sizeof(viewDesc));
						if (inputScalars.maskBPZ > 1) {
							arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr3DDesc.NumChannels = 1;
							arr3DDesc.Height = inputScalars.Nx[0];
							arr3DDesc.Width = inputScalars.Ny[0];
							arr3DDesc.Depth = inputScalars.Nz[0];
							status = cuArray3DCreate(&maskArrayBP, &arr3DDesc);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
							CUDA_MEMCPY3D cpy3d;
							std::memset(&cpy3d, 0, sizeof(cpy3d));
							cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
							cpy3d.srcHost = w_vec.maskBP;
							cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(uint8_t);
							cpy3d.srcHeight = inputScalars.Nx[0];
							cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
							cpy3d.dstArray = maskArrayBP;
							cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(uint8_t);
							cpy3d.Height = inputScalars.Nx[0];
							cpy3d.Depth = inputScalars.Nz[0];
							status = cuMemcpy3D(&cpy3d);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
						}
						else {
							arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr2DDesc.NumChannels = 1;
							arr2DDesc.Height = inputScalars.Nx[0];
							arr2DDesc.Width = inputScalars.Ny[0];
							status = cuArrayCreate(&maskArrayBP, &arr2DDesc);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
							CUDA_MEMCPY2D cpy2d;
							std::memset(&cpy2d, 0, sizeof(cpy2d));
							cpy2d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
							cpy2d.srcHost = w_vec.maskBP;
							cpy2d.srcPitch = inputScalars.Ny[0] * sizeof(uint8_t);
							cpy2d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
							cpy2d.dstArray = maskArrayBP;
							cpy2d.WidthInBytes = inputScalars.Ny[0] * sizeof(uint8_t);
							cpy2d.Height = inputScalars.Nx[0];
							status = cuMemcpy2D(&cpy2d);
							if (status != CUDA_SUCCESS) {
								getErrorString(status);
								return -1;
							}
						}
						resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
						resDesc.res.array.hArray = maskArrayBP;
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						//texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						texDesc.flags = CU_TRSF_READ_AS_INTEGER;
						viewDesc.height = inputScalars.Nx[0];
						viewDesc.width = inputScalars.Ny[0];
						//viewDesc.depth = 1;
						viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
						if (inputScalars.maskBPZ > 1) {
							viewDesc.depth = inputScalars.maskBPZ;
							texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						}
						status = cuTexObjectCreate(&d_maskBP, &resDesc, &texDesc, &viewDesc);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
					memAlloc.maskBP = true;
				}
			}
			if (inputScalars.eFOV) {
				status = cuMemAlloc(&d_eFOVIndices, sizeof(uint8_t) * inputScalars.Nz[0]);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.eFOV = true;
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				status = cuMemAlloc(&d_angle, sizeof(float) * inputScalars.nProjections);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.angle = true;
			}
			if (inputScalars.TOF) {
				// TOF bin centers
				status = cuMemAlloc(&d_TOFCenter, sizeof(float) * inputScalars.nBins);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.TOF = true;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				status = cuMemAlloc(&d_zFull[0], sizeof(float) * inputScalars.size_z);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.zFull = true;
			}
			if (inputScalars.SPECT) {
				status = cuMemAlloc(&d_rayShiftsDetector, sizeof(float) * 2 * inputScalars.n_rays);
				status = cuMemAlloc(&d_rayShiftsSource, sizeof(float) * 2 * inputScalars.n_rays);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memAlloc.rayShifts = true;
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if (DEBUG) {
					mexPrintBase("length[kk] = %u\n", length[kk]);
					mexPrintBase("kk = %u\n", kk);
					mexPrintBase("vecSize = %u\n", vecSize);
					mexEval();
				}
				if (inputScalars.CT && inputScalars.listmode != 1) {
					if (inputScalars.pitch) {
						status = cuMemAlloc(&d_z[kk], sizeof(float) * length[kk] * 6);
					}
					else
						status = cuMemAlloc(&d_z[kk], sizeof(float) * length[kk] * 2);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.zType = 1;
					memAlloc.zSteps++;
				}
				else {
					if (inputScalars.PET && inputScalars.listmode == 0) {
						if (inputScalars.nLayers > 1)
							status = cuMemAlloc(&d_z[kk], sizeof(float) * length[kk] * 3);
						else
							status = cuMemAlloc(&d_z[kk], sizeof(float) * length[kk] * 2);
						memAlloc.zType = 1;
						memAlloc.zSteps++;
					}
					else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
						status = cuMemAlloc(&d_z[kk], sizeof(float) * inputScalars.size_z);
						memAlloc.zType = 0;
						memAlloc.zSteps = kk;
					}
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					status = cuMemAlloc(&d_T[kk], sizeof(float) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.offsetT = true;
					memAlloc.oSteps++;
				}
				if (inputScalars.CT || (inputScalars.listmode > 0 && !inputScalars.indexBased)) {
					if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF || (inputScalars.CT && inputScalars.listmode == 0)) {
						status = cuMemAlloc(&d_x[kk], sizeof(float) * length[kk] * 6);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						memAlloc.xSteps++;
					}
				}
				if (inputScalars.normalization_correction) {
					status = cuMemAlloc(&d_norm[kk], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.norm = true;
					memAlloc.nSteps++;
				}
				if (inputScalars.scatter == 1U) {
					status = cuMemAlloc(&d_scat[kk], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.extra = true;
					memAlloc.eSteps++;
				}
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					status = cuMemAlloc(&d_atten[kk], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.attenM = true;
					memAlloc.aSteps++;
				}
				// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
				if (inputScalars.raw && inputScalars.listmode != 1) {
					status = cuMemAlloc(&d_L[kk], sizeof(uint16_t) * length[kk] * 2);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.raw = true;
					memAlloc.lSteps++;
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					status = cuMemAlloc(&d_xyindex[kk], sizeof(uint32_t) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = cuMemAlloc(&d_zindex[kk], sizeof(uint16_t) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memAlloc.subInd = true;
					memAlloc.iSteps++;
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						status = cuMemAlloc(&d_trIndex[kk], sizeof(uint16_t) * length[kk] * 2);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						status = cuMemAlloc(&d_axIndex[kk], sizeof(uint16_t) * length[kk] * 2);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						memAlloc.indexBased = true;
						memAlloc.iSteps++;
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.TOF) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						status = cuMemAlloc(&d_TOFIndex[kk], sizeof(uint8_t) * length[kk]);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						memAlloc.TOFIndex = true;
						memAlloc.TOFSteps++;
					}
				}
			}
		}

		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Buffer creation failed\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer creation succeeded\n");
		}


		// assign values to the buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			status = cuMemcpyHtoD(d_weights, w_vec.weights, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		//else if (MethodList.hyperbolic) {
		//	status = cuMemcpyHtoD(d_weights, w_vec.weights, sizeof(float) * ((w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1));
		//	if (status != CUDA_SUCCESS) {
		//		getErrorString(status);
		//		return -1;
		//	}
		//}
		if (inputScalars.projector_type != 6) {
			status = cuMemcpyHtoD(d_V, inputScalars.V, sizeof(float) * inputScalars.size_V);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((!inputScalars.CT && inputScalars.listmode == 0) || inputScalars.indexBased) {
				status = cuMemcpyHtoD(d_x[0], x, sizeof(float) * inputScalars.size_of_x);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				status = cuMemcpyHtoD(d_xFull[0], x, sizeof(float) * inputScalars.size_of_x);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = cuMemcpyHtoD(d_zFull[0], z_det, sizeof(float) * inputScalars.size_z);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			//status = cuMemcpyHtoD(d_xcenter, inputScalars.x_center, sizeof(float) * inputScalars.size_center_x);
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			//status = cuMemcpyHtoD(d_ycenter, inputScalars.y_center, sizeof(float) * inputScalars.size_center_y);
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			//status = cuMemcpyHtoD(d_zcenter, inputScalars.z_center, sizeof(float) * inputScalars.size_center_z);
			//if (status != CUDA_SUCCESS) {
			//	getErrorString(status);
			//	return -1;
			//}
			if ((inputScalars.maskFP || inputScalars.maskBP) && inputScalars.useBuffers) {
				if (inputScalars.maskFP) {
					if (inputScalars.maskFPZ > 1)
						for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
							status = cuMemcpyHtoD(d_maskFPB[kk], &w_vec.maskFP[pituus[kk] * vecSize], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk]);
					else
						status = cuMemcpyHtoD(d_maskFPB[0], w_vec.maskFP, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else if (inputScalars.maskBP) {
					status = cuMemcpyHtoD(d_maskBPB, w_vec.maskBP, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
			}
			if (inputScalars.attenuation_correction && inputScalars.CTAttenuation && inputScalars.useBuffers) {
				status = cuMemcpyHtoD(d_attenB, atten, sizeof(float) * inputScalars.im_dim[0]);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
				if (!inputScalars.useImages)
					status = cuMemcpyHtoD(d_uref, w_vec.NLM_ref, sizeof(float) * inputScalars.im_dim[0]);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				status = cuMemcpyHtoD(d_angle, w_vec.angles, sizeof(float) * inputScalars.nProjections);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.TOF) {
				status = cuMemcpyHtoD(d_TOFCenter, inputScalars.TOFCenter, sizeof(float) * inputScalars.nBins);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.eFOV) {
				status = cuMemcpyHtoD(d_eFOVIndices, w_vec.eFOVIndices, sizeof(uint8_t) * inputScalars.Nz[0]);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.SPECT) {
				status = cuMemcpyHtoD(d_rayShiftsDetector, w_vec.rayShiftsDetector, sizeof(float) * 2 * inputScalars.n_rays);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = cuMemcpyHtoD(d_rayShiftsSource, w_vec.rayShiftsSource, sizeof(float) * 2 * inputScalars.n_rays);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if (inputScalars.CT && inputScalars.listmode == 0) {
					if (inputScalars.pitch)
						status = cuMemcpyHtoD(d_z[kk], &z_det[pituus[kk] * 6], sizeof(float) * length[kk] * 6);
					else
						status = cuMemcpyHtoD(d_z[kk], &z_det[pituus[kk] * 2], sizeof(float) * length[kk] * 2);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else {
					if (inputScalars.PET && inputScalars.listmode == 0) {
						if (inputScalars.nLayers > 1)
							status = cuMemcpyHtoD(d_z[kk], &z_det[pituus[kk] * 3], sizeof(float) * length[kk] * 3);
						else
							status = cuMemcpyHtoD(d_z[kk], &z_det[pituus[kk] * 2], sizeof(float) * length[kk] * 2);
					}
					else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased))
						status = cuMemcpyHtoD(d_z[kk], z_det, sizeof(float) * inputScalars.size_z);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						//return -1;
					}
				}
				if (inputScalars.CT && inputScalars.listmode == 0) {
					status = cuMemcpyHtoD(d_x[kk], &x[pituus[kk] * 6], sizeof(float) * length[kk] * 6);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
					if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF) {
						status = cuMemcpyHtoD(d_x[kk], &w_vec.listCoord[pituus[kk] * 6], sizeof(float) * length[kk] * 6);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					status = cuMemcpyHtoD(d_T[kk], &inputScalars.T[pituus[kk]], sizeof(float) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.raw && inputScalars.listmode != 1) {
					status = cuMemcpyHtoD(d_L[kk], &L[pituus[kk] * 2], sizeof(uint16_t) * length[kk] * 2);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					status = cuMemcpyHtoD(d_zindex[kk], &z_index[pituus[kk]], sizeof(uint16_t) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = cuMemcpyHtoD(d_xyindex[kk], &xy_index[pituus[kk]], sizeof(uint32_t) * length[kk]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						status = cuMemcpyHtoD(d_trIndex[kk], &w_vec.trIndex[pituus[kk] * 2], sizeof(uint16_t) * length[kk] * 2);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						status = cuMemcpyHtoD(d_axIndex[kk], &w_vec.axIndex[pituus[kk] * 2], sizeof(uint16_t) * length[kk] * 2);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.TOF) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						status = cuMemcpyHtoD(d_TOFIndex[kk], &w_vec.TOFIndices[pituus[kk]], sizeof(uint8_t) * length[kk]);
						if (status != CUDA_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
				if (inputScalars.normalization_correction) {
					status = cuMemcpyHtoD(d_norm[kk], &norm[pituus[kk] * vecSize], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.scatter == 1U) {
					status = cuMemcpyHtoD(d_scat[kk], &extraCorr[pituus[kk] * vecSize], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					status = cuMemcpyHtoD(d_atten[kk], &atten[pituus[kk] * vecSize], sizeof(float) * length[kk] * vecSize);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (DEBUG) {
					mexPrintBase("length[kk] = %d\n", length[kk]);
					if (inputScalars.pitch && kk > 0) {
						mexPrintBase("z_det[pituus[kk] * 6] = %f\n", z_det[pituus[kk] * 6 - 1]);
						mexPrintBase("pituus[kk] * 6 = %d\n", pituus[kk] * 6 - 1);
					}
					mexPrintBase("pituus[kk] = %d\n", pituus[kk]);
					mexPrintBase("erotus = %d\n", pituus[kk + 1] - pituus[kk]);
					mexPrintBase("kk = %d\n", kk);
					mexEval();
				}
				status = cuCtxSynchronize();
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}

		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Buffer write failed\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer write succeeded\n");
		}
		return 0;
	}

	/// <summary>
	/// Resizes required vectors and then calls the function to create and write buffers. Also creates two necessary images
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="x the x/y/z coordinates for the detectors (PET and SPECT) or source and detector (CT). z-coordinate applies only for CT"></param>
	/// <param name="z_det the z coordinates for the detectors (PET and SPECT) or the directional vectors for the detector panel pixels (CT)"></param>
	/// <param name="xy_index subset indices for subsets types &lt; 8, x/y dimensions"></param>
	/// <param name="z_index same as above but for z dimension"></param>
	/// <param name="lor1 LORs to be discarded, i.e. the ray does not intersect the image"></param>
	/// <param name="L raw data detector indices"></param>
	/// <param name="pituus cumulative sum of length"></param>
	/// <param name="atten attenuation image"></param>
	/// <param name="norm normalization matrix"></param>
	/// <param name="extraCorr scatter data (for multiplicative scatter correction)"></param>
	/// <param name="V precomputed volume values for the volume of intersection based projector"></param>
	/// <param name="x_center x-coordinates of the voxel centers"></param>
	/// <param name="y_center y-coordinates of the voxel centers"></param>
	/// <param name="z_center z-coordinates of the voxel centers"></param>
	/// <param name="sc_ra randoms and/or scatter data (for additive scatter correction or for randoms correction)"></param>
	/// <param name="TOFCenter TOF bin center values"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="Sin measurement data (sinograms or projections)"></param>
	/// <param name="reko_type for reconstruction algorithms requiring unique operations in FP or BP"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <returns></returns>
	inline int createBuffers(scalarStruct& inputScalars, Weighting& w_vec, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* L, const int64_t* pituus, const float* atten, const float* norm, const float* extraCorr,
		const std::vector<int64_t>& length, const RecMethods& MethodList, const int type = 0) {
		int status = 0;
		//if (inputScalars.precompute)
		//	d_lor.resize(inputScalars.subsets);
		if (inputScalars.raw)
			d_L.resize(inputScalars.subsetsUsed);
		if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1) {
			d_xyindex.resize(inputScalars.subsetsUsed);
			d_zindex.resize(inputScalars.subsetsUsed);
		}
		if (inputScalars.listmode > 0 && inputScalars.indexBased) {
			d_trIndex.resize(inputScalars.subsetsUsed);
			d_axIndex.resize(inputScalars.subsetsUsed);
		}
		if (inputScalars.listmode > 0 && inputScalars.TOF) {
			d_TOFIndex.resize(inputScalars.subsetsUsed);
		}
		//if (inputScalars.randoms_correction)
		//	d_sc_ra.resize(inputScalars.subsets);
		if (inputScalars.normalization_correction)
			d_norm.resize(inputScalars.subsetsUsed);
		if (inputScalars.scatter)
			d_scat.resize(inputScalars.subsetsUsed);
		if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
			d_atten.resize(inputScalars.subsetsUsed);
		d_x.resize(inputScalars.subsetsUsed);
		d_z.resize(inputScalars.subsetsUsed);
		if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
			d_T.resize(inputScalars.subsetsUsed);
		//d_Sino.resize(inputScalars.TOFsubsets);
		//if (type < 2) {
		//	for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		//		cl::size_type imX = inputScalars.Nx[ii];
		//		cl::size_type imY = inputScalars.Ny[ii];
		//		cl::size_type imZ = inputScalars.Nz[ii];
		//		if (inputScalars.FPType == 5) {
		//			vec_opencl.d_image_os_int.emplace_back(CUtexObject(CLContext, CL_MEM_READ_ONLY, format, imY + 1, imZ + 1, imX, 0, 0, NULL, &status));
		//			imX++;
		//			cl::size_type aY = imY;
		//			imY = imZ + 1;
		//			imZ = aY;
		//		}
		//		vec_opencl.d_image_os.emplace_back(CUtexObject(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status));
		//		if (status != CUDA_SUCCESS) {
		//			getErrorString(status);
		//			mexPrint("Failed to create input images\n");
		//			return status;
		//		}
		//	}
		//}

		status = createAndWriteBuffers(length, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, inputScalars, w_vec, MethodList);
		if (status != 0) {
			return status;
		}
		return 0;
	}

	/// <summary>
	/// Inputs constant values to the kernels, i.e. values that do not change in each time step or iteration
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <returns></returns>
	inline int initializeKernel(scalarStruct& inputScalars, Weighting& w_vec) {
		int status = 0;

		if (inputScalars.FPType == 4 || inputScalars.FPType == 5) {
			FPArgs.emplace_back(&inputScalars.nRowsD);
			FPArgs.emplace_back(&inputScalars.nColsD);
			FPArgs.emplace_back(&dPitch);
		}

		if (inputScalars.BPType == 4 || inputScalars.BPType == 5) {
			BPArgs.emplace_back(&inputScalars.nRowsD);
			BPArgs.emplace_back(&inputScalars.nColsD);
			BPArgs.emplace_back(&dPitch);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				SensArgs.emplace_back(&inputScalars.nRowsD);
				SensArgs.emplace_back(&inputScalars.nColsD);
				SensArgs.emplace_back(&dPitch);
			}
		}
		if (inputScalars.FPType == 4) {
			FPArgs.emplace_back(&inputScalars.dL);
			FPArgs.emplace_back(&inputScalars.global_factor);
		}
		if (inputScalars.BPType == 4 && !inputScalars.CT) {
			BPArgs.emplace_back(&inputScalars.dL);
			BPArgs.emplace_back(&inputScalars.global_factor);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				SensArgs.emplace_back(&inputScalars.dL);
				SensArgs.emplace_back(&inputScalars.global_factor);
			}
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			// Set the kernelFP parameters that do not change
			FPArgs.emplace_back(&inputScalars.global_factor);
			FPArgs.emplace_back(&inputScalars.epps);
			FPArgs.emplace_back(&inputScalars.nRowsD);
			FPArgs.emplace_back(&inputScalars.det_per_ring);
			//FPArgs.emplace_back(&inputScalars.Nxy);
			FPArgs.emplace_back(&inputScalars.sigma_x);

			if (inputScalars.SPECT) {
				FPArgs.emplace_back(&d_rayShiftsDetector);
				FPArgs.emplace_back(&d_rayShiftsSource);
			}

			FPArgs.emplace_back(&dPitch);
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (inputScalars.FPType == 2)
					FPArgs.emplace_back(&inputScalars.tube_width);
				else
					FPArgs.emplace_back(&inputScalars.cylRadiusProj3);
				FPArgs.emplace_back(&inputScalars.bmin);
				FPArgs.emplace_back(&inputScalars.bmax);
				FPArgs.emplace_back(&inputScalars.Vmax);
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			BPArgs.emplace_back(&inputScalars.global_factor);
			BPArgs.emplace_back(&inputScalars.epps);
			BPArgs.emplace_back(&inputScalars.nRowsD);
			BPArgs.emplace_back(&inputScalars.det_per_ring);
			//BPArgs.emplace_back(&inputScalars.Nxy);
			BPArgs.emplace_back(&inputScalars.sigma_x);

			if (inputScalars.SPECT) {
				BPArgs.emplace_back(&d_rayShiftsDetector);
				BPArgs.emplace_back(&d_rayShiftsSource);
			}

			BPArgs.emplace_back(&dPitch);
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (inputScalars.FPType == 2)
					BPArgs.emplace_back(&inputScalars.tube_width);
				else
					BPArgs.emplace_back(&inputScalars.cylRadiusProj3);
				BPArgs.emplace_back(&inputScalars.bmin);
				BPArgs.emplace_back(&inputScalars.bmax);
				BPArgs.emplace_back(&inputScalars.Vmax);
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				SensArgs.emplace_back(&inputScalars.global_factor);
				SensArgs.emplace_back(&inputScalars.epps);
				SensArgs.emplace_back(&inputScalars.nRowsD);
				SensArgs.emplace_back(&inputScalars.det_per_ring);
				//SensArgs.emplace_back(&inputScalars.Nxy);
				SensArgs.emplace_back(&inputScalars.sigma_x);
				SensArgs.emplace_back(&dPitch);
				if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
					if (inputScalars.FPType == 2)
						SensArgs.emplace_back(&inputScalars.tube_width);
					else
						SensArgs.emplace_back(&inputScalars.cylRadiusProj3);
					SensArgs.emplace_back(&inputScalars.bmin);
					SensArgs.emplace_back(&inputScalars.bmax);
					SensArgs.emplace_back(&inputScalars.Vmax);
				}
			}
		}
		//if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
		//	BPArgs.emplace_back(&inputScalars.T);
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if (inputScalars.TOF) {
				FPArgs.emplace_back(&d_TOFCenter);
				if (DEBUG) {
					mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
					mexEval();
				}
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				//FPArgs.emplace_back(&d_xcenter);
				//FPArgs.emplace_back(&d_ycenter);
				//FPArgs.emplace_back(&d_zcenter);
				FPArgs.emplace_back(&d_V);
			}
			FPArgs.emplace_back(&inputScalars.nColsD);
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if (inputScalars.TOF) {
				BPArgs.emplace_back(&d_TOFCenter);
				if (DEBUG) {
					mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
					mexEval();
				}
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				//BPArgs.emplace_back(&d_xcenter);
				//BPArgs.emplace_back(&d_ycenter);
				//BPArgs.emplace_back(&d_zcenter);
				BPArgs.emplace_back(&d_V);
			}
			BPArgs.emplace_back(&inputScalars.nColsD);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.TOF) {
					SensArgs.emplace_back(&d_TOFCenter);
				}
				if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
					SensArgs.emplace_back(&d_V);
				}
				SensArgs.emplace_back(&inputScalars.nColsD);
			}
		}
		if ((inputScalars.BPType == 4 || inputScalars.FPType == 4) && !inputScalars.CT && inputScalars.TOF) {
			if (inputScalars.FPType == 4) {
				FPArgs.emplace_back(&d_TOFCenter);
				FPArgs.emplace_back(&inputScalars.sigma_x);
			}
			if (inputScalars.BPType == 4) {
				BPArgs.emplace_back(&d_TOFCenter);
				BPArgs.emplace_back(&inputScalars.sigma_x);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					SensArgs.emplace_back(&d_TOFCenter);
					SensArgs.emplace_back(&inputScalars.sigma_x);
				}
			}
		}
		if (DEBUG) {
			mexPrintBase("kernelIndFP = %u\n", FPArgs.size());
			mexPrintBase("kernelIndBP = %u\n", BPArgs.size());
			mexEval();
		}
		return 0;
	}

	/// <summary>
	/// Loads "dynamic" data, i.e. if the reconstruction is dynamic (time-varying) then this loads the measurement data/randoms/scatter for the next time step
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="Sino measurement data (sinograms or projections)"></param>
	/// <param name="randomsData randoms and/or scatter data (for additive scatter correction or for randoms correction)"></param>
	/// <param name="extraCorr scatter data (for multiplicative scatter correction)"></param>
	/// <param name="pituus cumulative sum of length"></param>
	/// <returns></returns>
	inline int loadDynamicData(scalarStruct& inputScalars, const std::vector<int64_t>& length, const float* extraCorr, const int64_t* pituus, const uint32_t tt) {

		CUresult status = CUDA_SUCCESS;
		for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
			if (inputScalars.scatter == 1u) {
				status = cuMemcpyHtoD(d_scat[kk], &extraCorr[pituus[kk] + inputScalars.koko * tt], sizeof(float) * length[kk]);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		return 0;
	}

	/// <summary>
	/// Sets kernel parameters that do not change per iteration but only per time step
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <returns></returns>
	inline int setDynamicKernelData(scalarStruct& inputScalars, Weighting& w_vec) {
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.attenuation_correction && !inputScalars.CT && inputScalars.CTAttenuation) {
			if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 || inputScalars.FPType == 4)) {
				if (inputScalars.useBuffers)
					FPArgs.emplace_back(&d_attenB);
				else
					FPArgs.emplace_back(&d_attenIm);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.BPType == 4) {
				if (inputScalars.useBuffers)
					BPArgs.emplace_back(&d_attenB);
				else
					BPArgs.emplace_back(&d_attenIm);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					if (inputScalars.useBuffers)
						SensArgs.emplace_back(&d_attenB);
					else
						SensArgs.emplace_back(&d_attenIm);
				}
			}
		}
		return 0;
	}

	template <typename T>
	inline int loadCoord(scalarStruct& inputScalars, const int64_t length, const T* listCoord, const T* listCoordAx = nullptr, const uint8_t* TOFIndices = nullptr) {
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.indexBased) {
			getErrorString(cuMemFree(d_trIndex[0]));
			getErrorString(cuMemFree(d_axIndex[0]));
			status = cuMemAlloc(&d_trIndex[0], sizeof(uint16_t) * length * 2);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = cuMemAlloc(&d_axIndex[0], sizeof(uint16_t) * length * 2);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = cuMemcpyHtoD(d_trIndex[0], listCoord, sizeof(uint16_t) * length * 2);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = cuMemcpyHtoD(d_axIndex[0], listCoordAx, sizeof(uint16_t) * length * 2);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		else {
			getErrorString(cuMemFree(d_x[0]));
			status = cuMemAlloc(&d_x[0], sizeof(float) * length * 6);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = cuMemcpyHtoD(d_x[0], listCoord, sizeof(float) * length * 6);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.TOF) {
			getErrorString(cuMemFree(d_TOFIndex[0]));
			status = cuMemAlloc(&d_TOFIndex[0], sizeof(uint8_t) * length);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = cuMemcpyHtoD(d_TOFIndex[0], TOFIndices, sizeof(uint8_t) * length);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		return 0;
	}

	/// <summary>
	/// Compute the forward projection for the selected projector type
	/// </summary>
	/// <param name="vec image estimates and backprojection"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="outputFP the output forward projection array"></param>
	/// <param name="osa_iter current subset (sub-iteration)"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="m_size for projector types 1-3, the total number of LORs"></param>
	/// <returns></returns>
	inline int forwardProjection(scalarStruct& inputScalars, Weighting& w_vec, uint32_t osa_iter, const std::vector<int64_t>& length, uint64_t m_size, int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kTemp = FPArgs;
		if (inputScalars.FPType == 5) {
			global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
			global[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1]) / local[1];
			global[2] = length[osa_iter];
		}
		else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
			global[1] = (inputScalars.nColsD  + erotus[1]) / local[1];
			global[2] = length[osa_iter];
		}
		else {
			erotus[0] = length[osa_iter] % local_size[0];

			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
			global[0] = (length[osa_iter] + erotus[0]) / local[0];
			global[1] = 1;
			global[2] = 1;
		}

		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("local[0] = %u\n", local[0]);
			mexPrintBase("local[1] = %u\n", local[1]);
			mexPrintBase("local[2] = %u\n", local[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("d[ii].s0 = %f\n", d[ii].x);
			mexPrintBase("d[ii].s1 = %f\n", d[ii].y);
			mexPrintBase("d[ii].s2 = %f\n", d[ii].z);
			mexPrintBase("b[ii].s0 = %f\n", b[ii].x);
			mexPrintBase("b[ii].s1 = %f\n", b[ii].y);
			mexPrintBase("b[ii].s2 = %f\n", b[ii].z);
			mexPrintBase("bmax[ii].s0 = %f\n", bmax[ii].x);
			mexPrintBase("bmax[ii].s1 = %f\n", bmax[ii].y);
			mexPrintBase("bmax[ii].s2 = %f\n", bmax[ii].z);
			mexPrintBase("kernelIndFPSubIter = %u\n", kernelIndFPSubIter);
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
			mexPrintBase("size_y = %u\n", inputScalars.nColsD);
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexPrintBase("maskBP = %u\n", inputScalars.maskBP);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("ii = %u\n", ii);
			mexPrintBase("NVOXELS = %u\n", NVOXELS);
			mexPrintBase("NVOXELS5 = %u\n", NVOXELS5);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexEval();
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.FPType == 4) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
				kTemp.emplace_back(&d_atten[osa_iter]);
		}
		if (inputScalars.FPType == 5 || inputScalars.FPType == 4) {
			kTemp.emplace_back(&d_N[ii]);
			kTemp.emplace_back(&b[ii]);
			if (inputScalars.FPType == 5) {
				kTemp.emplace_back(&inputScalars.dSize[ii]);
				kTemp.emplace_back(&d[ii]);
				kTemp.emplace_back(&inputScalars.d_Scale[ii]);
			}
			else {
				kTemp.emplace_back(&bmax[ii]);
				kTemp.emplace_back(&inputScalars.d_Scale4[ii]);
			}
		}
		//mexPrint("1!!!!\n");
		if (inputScalars.FPType == 4) {
			kTemp.emplace_back(&vec_opencl.d_image_os);
			kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
			//mexPrint("2!!!!\n");
			if ((inputScalars.listmode == 0 && !inputScalars.CT) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
				kTemp.emplace_back(&d_x[0]);
			else
				kTemp.emplace_back(&d_x[osa_iter]);
			if ((inputScalars.CT || inputScalars.PET || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
				kTemp.emplace_back(&d_z[osa_iter]);
			else
				kTemp.emplace_back(&d_z[inputScalars.osa_iter0]);
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						kTemp.emplace_back(&d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1)
							kTemp.emplace_back(&d_maskFP3[osa_iter]);
						else
							kTemp.emplace_back(&d_maskFP);
				}
			}
			kTemp.emplace_back((void*)&length[osa_iter]);
			//mexPrint("3!!!!\n");
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				kTemp.emplace_back(&d_xyindex[osa_iter]);
				kTemp.emplace_back(&d_zindex[osa_iter]);
			}
			if (inputScalars.raw) {
				kTemp.emplace_back(&d_L[osa_iter]);
				kTemp.emplace_back(&inputScalars.det_per_ring);
			}
			if (inputScalars.normalization_correction)
				kTemp.emplace_back(&d_norm[osa_iter]);
			if (inputScalars.scatter)
				kTemp.emplace_back(&d_scat[osa_iter]);
			//mexPrint("4!!!!\n");
			kTemp.emplace_back(&no_norm);
			kTemp.emplace_back(&m_size);
			kTemp.emplace_back(&osa_iter);
			kTemp.emplace_back(&ii);
			//mexPrint("5!!!!\n");
		}
		else if (inputScalars.FPType == 5) {
			if (!inputScalars.loadTOF && inputScalars.listmode > 0)
				kTemp.emplace_back(&d_x[0]);
			else
				kTemp.emplace_back(&d_x[osa_iter]);
			kTemp.emplace_back(&d_z[osa_iter]);
			kTemp.emplace_back(&vec_opencl.d_image_os);
			kTemp.emplace_back(&vec_opencl.d_image_os_int);
			kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
			if (inputScalars.meanFP) {

			}
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						kTemp.emplace_back(&d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1)
							kTemp.emplace_back(&d_maskFP3[osa_iter]);
						else
							kTemp.emplace_back(&d_maskFP);
				}
			}
			kTemp.emplace_back((void*)&length[osa_iter]);
		}
		else if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
				kTemp.emplace_back(&d_atten[osa_iter]);
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						kTemp.emplace_back(&d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1)
							kTemp.emplace_back(&d_maskFP3[osa_iter]);
						else
							kTemp.emplace_back(&d_maskFP);
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0) {
				kTemp.emplace_back((void*)&length[osa_iter]);
			}
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !inputScalars.CT) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
				kTemp.emplace_back(&d_x[0]);
			else
				kTemp.emplace_back(&d_x[osa_iter]);
			if ((inputScalars.CT || inputScalars.PET || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
				kTemp.emplace_back(&d_z[osa_iter]);
			else
				kTemp.emplace_back(&d_z[inputScalars.osa_iter0]);
			if (inputScalars.normalization_correction)
				kTemp.emplace_back(&d_norm[osa_iter]);
			if (inputScalars.scatter)
				kTemp.emplace_back(&d_scat[osa_iter]);
			kTemp.emplace_back(reinterpret_cast<void*>(&d_Summ[uu]));
			kTemp.emplace_back(&d_N[ii]);
			kTemp.emplace_back(&d[ii]);
			kTemp.emplace_back(&b[ii]);
			kTemp.emplace_back(&bmax[ii]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				kTemp.emplace_back(&d_xyindex[osa_iter]);
				kTemp.emplace_back(&d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				if (!inputScalars.loadTOF) {
					kTemp.emplace_back(&d_trIndex[0]);
					kTemp.emplace_back(&d_axIndex[0]);
				}
				else {
					kTemp.emplace_back(&d_trIndex[osa_iter]);
					kTemp.emplace_back(&d_axIndex[osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				if (!inputScalars.loadTOF) {
					kTemp.emplace_back(&d_TOFIndex[0]);
				}
				else {
					kTemp.emplace_back(&d_TOFIndex[osa_iter]);
				}
			}
			if (inputScalars.raw) {
				kTemp.emplace_back(&d_L[osa_iter]);
			}
			if (inputScalars.useBuffers)
				kTemp.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_im));
			else
				kTemp.emplace_back(&vec_opencl.d_image_os);
			kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
			kTemp.emplace_back(&no_norm);
			kTemp.emplace_back(&m_size);
			kTemp.emplace_back(&osa_iter);
			kTemp.emplace_back(&ii);
		}
		status = cuLaunchKernel(kernelFP, global[0], global[1], global[2], local[0], local[1], local[2], 0, CLCommandQueue[0], kTemp.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward projection kernel launched successfully\n");
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (!inputScalars.useBuffers) {
			status = cuTexObjectDestroy(vec_opencl.d_image_os);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(FPArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (inputScalars.FPType == 5) {
			status = cuTexObjectDestroy(vec_opencl.d_image_os_int);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(integArrayXY);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Forward projection completed");
		return 0;
	}

	/// <summary>
	/// Compute the backprojection for the selected projector type
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="osa_iter current subset (sub-iteration)"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="m_size for projector types 1-3, the total number of LORs"></param>
	/// <param name="compSens if true, computes the sensitivity image as well"></param>
	/// <returns></returns>
	inline int backwardProjection(scalarStruct& inputScalars, Weighting& w_vec, uint32_t osa_iter, std::vector<int64_t>& length, uint64_t m_size, const bool compSens = false, int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kTemp = BPArgs;
		if (inputScalars.listmode > 0 && compSens) {
			kernelApu = kernelBP;
			kernelBP = kernelSensList;
			kTemp = SensArgs;
		}

		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
				global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
				global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
				global[2] = length[osa_iter];
			}
			else if (inputScalars.listmode > 0 && compSens) {
				global[0] = static_cast<size_t>(inputScalars.det_per_ring + erotusSens[0]) / local[0];
				global[1] = (static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1]) / local[1];
				global[2] = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.nLayers);
			}
			else {
				erotus[0] = length[osa_iter] % local_size[0];

				if (erotus[0] > 0)
					erotus[0] = (local_size[0] - erotus[0]);
				global[0] = (length[osa_iter] + erotus[0]) / local[0];
				global[1] = 1;
				global[2] = 1;
			}

			if (DEBUG) {
				mexPrintBase("global[0] = %u\n", global[0]);
				mexPrintBase("local[0] = %u\n", local[0]);
				mexPrintBase("local[1] = %u\n", local[1]);
				mexPrintBase("global[1] = %u\n", global[1]);
				mexPrintBase("global[2] = %u\n", global[2]);
				if (inputScalars.listmode > 0 && compSens) {
					mexPrintBase("erotusSens[0] = %u\n", erotusSens[0]);
					mexPrintBase("erotusSens[1] = %u\n", erotusSens[1]);
				}
				else {
				mexPrintBase("erotus[0] = %u\n", erotus[0]);
				mexPrintBase("erotus[1] = %u\n", erotus[1]);
				}
				mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
				mexPrintBase("m_size = %u\n", m_size);
				mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
				mexPrintBase("size_y = %u\n", inputScalars.nColsD);
				mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
				mexPrintBase("listmode = %u\n", inputScalars.listmode);
				mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("no_norm = %u\n", no_norm);
				mexPrintBase("compSens = %u\n", compSens);
				mexEval();
				//mexEvalString("pause(2);");
			}

			// Set kernelBP arguments
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
				kTemp.emplace_back(&d_atten[osa_iter]);
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						kTemp.emplace_back(&d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1)
							kTemp.emplace_back(&d_maskFP3[osa_iter]);
						else
							kTemp.emplace_back(&d_maskFP);
				}
				if (inputScalars.maskBP) {
					if (inputScalars.useBuffers)
						kTemp.emplace_back(&d_maskBPB);
					else
						kTemp.emplace_back(&d_maskBP);
					if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
						if (inputScalars.useBuffers)
							kTemp.emplace_back(&d_maskBPB);
						else
							kTemp.emplace_back(&d_maskBP);
					}
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0)
				kTemp.emplace_back(&length[osa_iter]);
			if (compSens) {
				kTemp.emplace_back(&d_xFull[0]);
				kTemp.emplace_back(&d_zFull[0]);
				kTemp.emplace_back(&inputScalars.rings);
			}
			else {
				if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !inputScalars.CT) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
					kTemp.emplace_back(&d_x[0]);
				else
					kTemp.emplace_back(&d_x[osa_iter]);
				if ((inputScalars.CT || inputScalars.PET || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
					kTemp.emplace_back(&d_z[osa_iter]);
				else if (inputScalars.indexBased && inputScalars.listmode > 0)
					kTemp.emplace_back(&d_z[0]);
				else
					kTemp.emplace_back(&d_z[inputScalars.osa_iter0]);
			}
			if (compSens) {
				if (inputScalars.normalization_correction)
					kTemp.emplace_back(&d_normFull[0]);
				if (inputScalars.scatter)
					kTemp.emplace_back(&d_scatFull[0]);
			}
			else {
				if (inputScalars.normalization_correction)
					kTemp.emplace_back(&d_norm[osa_iter]);
				if (inputScalars.scatter)
					kTemp.emplace_back(&d_scat[osa_iter]);
			}
			kTemp.emplace_back(reinterpret_cast<void*>(&d_Summ[uu]));
			kTemp.emplace_back(&d_N[ii]);
			kTemp.emplace_back(&d[ii]);
			kTemp.emplace_back(&b[ii]);
			kTemp.emplace_back(&bmax[ii]);
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
					kTemp.emplace_back(&d_xyindex[osa_iter]);
					kTemp.emplace_back(&d_zindex[osa_iter]);
				}
			if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
				if (!inputScalars.loadTOF) {
					kTemp.emplace_back(&d_trIndex[0]);
					kTemp.emplace_back(&d_axIndex[0]);
				}
				else {
					kTemp.emplace_back(&d_trIndex[osa_iter]);
					kTemp.emplace_back(&d_axIndex[osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				if (!inputScalars.loadTOF) {
					kTemp.emplace_back(&d_TOFIndex[0]);
				}
				else {
					kTemp.emplace_back(&d_TOFIndex[osa_iter]);
				}
			}
			if (inputScalars.raw)
				kTemp.emplace_back(&d_L[osa_iter]);
			kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
			kTemp.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_rhs_os[uu]));
			kTemp.emplace_back(&no_norm);
			kTemp.emplace_back(&m_size);
			kTemp.emplace_back(&osa_iter);
			kTemp.emplace_back(&ii);
		}
		else {
			if (inputScalars.CT) {

				//arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
				//arr3DDesc.NumChannels = 1;
				//arr3DDesc.Height = inputScalars.Nx[0];
				//arr3DDesc.Width = inputScalars.Ny[0];
				//arr3DDesc.Depth = inputScalars.Nz[0];
				//status = cuArray3DCreate(&uRefArray, &arr3DDesc);
				//CUDA_MEMCPY3D cpy3d;
				//std::memset(&cpy3d, 0, sizeof(cpy3d));
				//cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
				//cpy3d.srcHost = w_vec.NLM_ref;
				//cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(float);
				//cpy3d.srcHeight = inputScalars.Nx[0];
				//cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
				//cpy3d.dstArray = uRefArray;
				//cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(float);
				//cpy3d.Height = inputScalars.Nx[0];
				//cpy3d.Depth = inputScalars.Nz[0];
				//status = cuMemcpy3D(&cpy3d);
				if (!inputScalars.useBuffers) {
					std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
					arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
					arr3DDesc.NumChannels = 1;
					arr3DDesc.Height = inputScalars.nColsD;
					arr3DDesc.Width = inputScalars.nRowsD;
					arr3DDesc.Depth = length[osa_iter];
					if (inputScalars.BPType == 5) {
						arr3DDesc.Height++;
						arr3DDesc.Width++;
					}
					if (DEBUG) {
						mexPrintBase("arr3DDesc.NumChannels= %u\n", arr3DDesc.NumChannels);
						mexPrintBase("arr3DDesc.Height = %u\n", arr3DDesc.Height);
						mexPrintBase("arr3DDesc.Width = %u\n", arr3DDesc.Width);
						mexPrintBase("arr3DDesc.Depth = %u\n", arr3DDesc.Depth);
						mexEval();
						//mexEvalString("pause(2);");
					}
					status = cuArray3DCreate(&BPArray, &arr3DDesc);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						mexPrint("Array creation failed\n");
						return -1;
					}
					else if (DEBUG)
						mexPrint("Array creation succeeded\n");
					CUDA_MEMCPY3D cpy3d;
					std::memset(&cpy3d, 0, sizeof(cpy3d));
					cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
					cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(d_output);
					cpy3d.srcPitch = inputScalars.nRowsD * sizeof(float);
					cpy3d.srcHeight = inputScalars.nColsD;
					cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
					cpy3d.dstArray = BPArray;
					cpy3d.WidthInBytes = inputScalars.nRowsD * sizeof(float);
					cpy3d.Height = inputScalars.nColsD;
					cpy3d.Depth = length[osa_iter];
					if (inputScalars.BPType == 5) {
						cpy3d.srcPitch += sizeof(float);
						cpy3d.srcHeight++;
						cpy3d.WidthInBytes += sizeof(float);
						cpy3d.Height++;
					}
					status = cuMemcpy3D(&cpy3d);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						mexPrint("Array mem copy failed\n");
						return -1;
					}
					else if (DEBUG)
						mexPrint("Array mem copy succeeded\n");
					CUDA_RESOURCE_DESC resDescIm;
					std::memset(&resDescIm, 0, sizeof(resDescIm));
					std::memset(&texDesc, 0, sizeof(texDesc));
					std::memset(&viewDesc, 0, sizeof(viewDesc));
					viewDesc.height = inputScalars.nColsD;
					viewDesc.width = inputScalars.nRowsD;
					viewDesc.depth = length[osa_iter];
					viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
					resDescIm.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
					resDescIm.res.array.hArray = BPArray;
					if (inputScalars.BPType == 4) {
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
					}
					else {
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
						texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
						viewDesc.height++;
						viewDesc.width++;
					}
					status = cuTexObjectCreate(&d_inputImage, &resDescIm, &texDesc, &viewDesc);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						mexPrint("Image creation failed\n");
						return -1;
					}
					status = cuCtxSynchronize();
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
						mexPrint("Queue finish failed after image copy\n");
						return -1;
					}
				}
				if (inputScalars.BPType == 4) {
					global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
					global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
					if (!inputScalars.largeDim)
						global[2] = (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS;
					else
						global[2] = inputScalars.Nz[ii];
				}
				else if (inputScalars.BPType == 5) {
					if (inputScalars.pitch) {
						global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
						global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
						global[2] = inputScalars.Nz[ii];
					}
					else {
						global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
						global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
						global[2] = (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5;
					}
				}
				else {
					global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
					global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
					global[2] = inputScalars.Nz[ii];
				}

				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("local[0] = %u\n", local[0]);
					mexPrintBase("local[1] = %u\n", local[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
					mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
					mexPrintBase("size_y = %u\n", inputScalars.nColsD);
					mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexEval();
					//mexEvalString("pause(2);");
				}
				if (inputScalars.offset)
					kTemp.emplace_back(&d_T[osa_iter]);
				if (inputScalars.BPType == 5 || inputScalars.BPType == 4) {
					kTemp.emplace_back(&d_N[ii]);
					kTemp.emplace_back(&b[ii]);
					kTemp.emplace_back(&d[ii]);
					if (inputScalars.BPType == 5) {
						kTemp.emplace_back(&inputScalars.d_Scale[ii]);
						kTemp.emplace_back(&inputScalars.dSizeBP);
					}
					else {
						kTemp.emplace_back(&w_vec.kerroin4[ii]);
					}
				}
				//mexPrint("1!!!!\n");
				if (inputScalars.BPType == 4) {
					if (inputScalars.useBuffers)
						kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
					else
						kTemp.emplace_back(&d_inputImage);
					if (inputScalars.CT && inputScalars.DSC > 0.f) {
						kTemp.emplace_back(&d_angle);
						kTemp.emplace_back(&inputScalars.DSC);
					}
					//mexPrint("2!!!!\n");
					kTemp.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_rhs_os[uu]));
					//mexPrint("3!!!!\n");
					if (compSens)
						kTemp.emplace_back(&d_xFull[0]);
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0)
							kTemp.emplace_back(&d_x[0]);
						else
							kTemp.emplace_back(&d_x[osa_iter]);
					if (compSens)
						kTemp.emplace_back(&d_zFull[0]);
					else
						kTemp.emplace_back(&d_z[osa_iter]);
					//mexPrint("4!!!!\n");
					kTemp.emplace_back(reinterpret_cast<void*>(&d_Summ[uu]));
					//mexPrint("5!!!!\n");
				}
				else {
					if (compSens)
						kTemp.emplace_back(&d_xFull[0]);
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0)
							kTemp.emplace_back(&d_x[0]);
						else
							kTemp.emplace_back(&d_x[osa_iter]);
					if (compSens)
						kTemp.emplace_back(&d_zFull[0]);
					else
						kTemp.emplace_back(&d_z[osa_iter]);
					kTemp.emplace_back(&d_inputImage);
					kTemp.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_rhs_os[uu]));
					kTemp.emplace_back(reinterpret_cast<void*>(&d_Summ[uu]));
					if (inputScalars.meanBP) {
						kTemp.emplace_back(&d_meanBP);
					}
				}
			}
			else {
				if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
					global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
					global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
					global[2] = length[osa_iter];
				}
				else if (inputScalars.listmode > 0 && compSens) {
					global[0] = static_cast<size_t>(inputScalars.det_per_ring + erotusSens[0]) / local[0];
					global[1] = (static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1]) / local[1];
					global[2] = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings);
				}
				else {
					erotus[0] = length[osa_iter] % local_size[0];

					if (erotus[0] > 0)
						erotus[0] = (local_size[0] - erotus[0]);
					global[0] = (length[osa_iter] + erotus[0]) / local[0];
					global[1] = 1;
					global[2] = 1;
				}

				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("local[0] = %u\n", local[0]);
					mexPrintBase("local[1] = %u\n", local[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					mexPrintBase("erotus[0] = %u\n", erotus[0]);
					mexPrintBase("erotus[1] = %u\n", erotus[1]);
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
					mexPrintBase("size_y = %u\n", inputScalars.nColsD);
					mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("osa_iter = %u\n", osa_iter);
					mexEval();
					//mexEvalString("pause(2);");
				}
				kTemp.emplace_back(&d_N[ii]);
				kTemp.emplace_back(&b[ii]);
				kTemp.emplace_back(&bmax[ii]);
				kTemp.emplace_back(&inputScalars.d_Scale4[ii]);
				kTemp.emplace_back(reinterpret_cast<void*>(&d_output));
				kTemp.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_rhs_os[uu]));
				if (compSens) {
					kTemp.emplace_back(&d_xFull[0]);
					kTemp.emplace_back(&d_zFull[0]);
					kTemp.emplace_back(&inputScalars.rings);
					kTemp.emplace_back(&inputScalars.det_per_ring);
				}
				else {
					if ((inputScalars.listmode == 0 && !inputScalars.CT) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
						kTemp.emplace_back(&d_x[0]);
					else
						kTemp.emplace_back(&d_x[osa_iter]);
					if ((inputScalars.CT || inputScalars.PET || inputScalars.listmode > 0))
						kTemp.emplace_back(&d_z[osa_iter]);
					else
						kTemp.emplace_back(&d_z[inputScalars.osa_iter0]);
				}
				kTemp.emplace_back(&length[osa_iter]);
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
					kTemp.emplace_back(&d_xyindex[osa_iter]);
					kTemp.emplace_back(&d_zindex[osa_iter]);
				}
				if (inputScalars.raw) {
					kTemp.emplace_back(&d_L[osa_iter]);
					kTemp.emplace_back(&inputScalars.det_per_ring);
				}
				if (inputScalars.normalization_correction)
					kTemp.emplace_back(&d_norm[osa_iter]);
				if (inputScalars.scatter)
					kTemp.emplace_back(&d_scat[osa_iter]);
				kTemp.emplace_back(reinterpret_cast<void*>(&d_Summ[uu]));
			}
			kTemp.emplace_back(&no_norm);
			//mexPrint("6!!!!\n");
			if (inputScalars.maskBP) {
				if (inputScalars.useBuffers)
					kTemp.emplace_back(&d_maskBPB);
				else
					kTemp.emplace_back(&d_maskBP);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					if (inputScalars.useBuffers)
						kTemp.emplace_back(&d_maskBPB);
					else
						kTemp.emplace_back(&d_maskBP);
				}
			}
			if (inputScalars.CT)
				kTemp.emplace_back(&length[osa_iter]);
			else {
				kTemp.emplace_back(&m_size);
				kTemp.emplace_back(&osa_iter);
			}
			kTemp.emplace_back(&ii);
			//mexPrint("7!!!!\n");
		}
		status = cuLaunchKernel(kernelBP, global[0], global[1], global[2], local[0], local[1], local[2], 0, CLCommandQueue[0], kTemp.data(), 0);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Backprojection kernel launched successfully\n");
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5) {
			if (!inputScalars.useBuffers) {
				status = cuTexObjectDestroy(d_inputImage);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
				status = cuArrayDestroy(BPArray);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
			}
		}
		if (inputScalars.listmode > 0 && compSens) {
			kernelBP = kernelApu;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Backprojection computed");
		return 0;
	}

	/// <summary>
	/// Release buffers needed only by the initial computation of the sensitivity image including all measurements (e.g. image-based preconditioners 2-3)
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <returns></returns>
	inline void releaseBuffer(const scalarStruct& inputScalars) {
		d_xFull.clear();
		d_zFull.clear();
		if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
			d_normFull.clear();
		}
		if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) {
			d_scatFull.clear();
		}
		//if (inputScalars.precompute) {
		//	d_lorFull.clear();
		//}
	}

	/// <summary>
	/// Get the total global memory of the selected device
	/// </summary>
	inline int64_t getGlobalMem() {
		CUresult status = CUDA_SUCCESS;
		size_t mem;
		size_t memF;
		unsigned long long mem_loc;
		status = cuMemGetInfo(&memF, &mem);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		int apu = 0;
		cuDeviceGetAttribute(&apu, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, CUDeviceID[0]);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (DEBUG) {
			mexPrintBase("mem_loc = %u\n", apu);
			mexPrintBase("memFree = %u\n", memF);
		}
		return mem;
	}

	/// <summary>
	/// Compute median root prior (MRP)
	/// </summary>
	/// <param name="padd the padded input array (current estimate)"></param>
	/// <param name="grad the output gradient array"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <returns></returns>
	inline int computeMRP(const scalarStruct& inputScalars, const uint64_t global_size[]) {
		std::vector<void*> kArgs;
		CUresult status = CUDA_SUCCESS;
		unsigned int gSize[3];
		unsigned int erotus[2];
		erotus[0] = localPrior[0] - (global_size[0] % localPrior[0]);
		erotus[1] = localPrior[1] - (global_size[1] % localPrior[1]);
		gSize[0] = (global_size[0] + erotus[0]) / localPrior[0];
		gSize[1] = (global_size[1] + erotus[1]) / localPrior[1];
		gSize[2] = global_size[2];
		status = cuCtxSynchronize();
		kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NOrig);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		status = cuLaunchKernel(kernelMed, gSize[0], gSize[1], gSize[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Median filter kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Median kernel launched successfully\n");
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after MRP kernel\n");
			return -1;
		}
		return 0;
	}

	/// <summary>
	/// Non-local means (NLM) prior
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <returns></returns>
	inline int computeNLM(const scalarStruct& inputScalars, Weighting& w_vec, float beta) {
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA NLM gradient computation");
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
		if (inputScalars.largeDim)
			globalPrior[2] = inputScalars.Nz[0];
		status = cuCtxSynchronize();
		const int3 searchWindow = { static_cast<int>(w_vec.Ndx) , static_cast<int>(w_vec.Ndy) , static_cast<int>(w_vec.Ndz) };
		const int3 patchWindow = { static_cast<int>(w_vec.Nlx) , static_cast<int>(w_vec.Nly) , static_cast<int>(w_vec.Nlz) };
		if (DEBUG) {
			mexPrintBase("w_vec.Ndx = %u\n", w_vec.Ndx);
			mexPrintBase("w_vec.Ndy = %u\n", w_vec.Ndy);
			mexPrintBase("w_vec.Ndz = %u\n", w_vec.Ndz);
			mexPrintBase("w_vec.Nlx = %u\n", w_vec.Nlx);
			mexPrintBase("w_vec.Nly = %u\n", w_vec.Nly);
			mexPrintBase("w_vec.Nlz = %u\n", w_vec.Nlz);
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPrior[0] = %u\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %u\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %u\n", globalPrior[2]);
			mexPrintBase("localPrior[0] = %u\n", localPrior[0]);
			mexPrintBase("localPrior[1] = %u\n", localPrior[1]);
			mexPrintBase("localPrior[2] = %u\n", localPrior[2]);
			mexPrintBase("w_vec.h2 = %f\n", w_vec.h2);
			mexPrintBase("w_vec.RDP_gamma = %f\n", w_vec.RDP_gamma);
			mexPrintBase("useImages = %d\n", inputScalars.useImages);
			mexEval();
			//mexEvalString("pause(2);");
		}
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		if (inputScalars.useImages) {
			kArgs.emplace_back(&d_inputI);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		}
		kArgs.emplace_back(&d_gaussianNLM);
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NOrig);
		kArgs.emplace_back(&w_vec.h2);
		kArgs.emplace_back(&apu);
		kArgs.emplace_back(&beta);
		if (w_vec.NLRD || w_vec.NLLange || w_vec.NLGGMRF)
			kArgs.emplace_back(&w_vec.RDP_gamma);
		if (w_vec.NLGGMRF) {
			kArgs.emplace_back(&w_vec.GGMRF_p);
			kArgs.emplace_back(&w_vec.GGMRF_q);
			kArgs.emplace_back(&w_vec.GGMRF_c);
		}
		if (w_vec.NLAdaptive)
			kArgs.emplace_back(&w_vec.NLAdaptiveConstant);
		if (w_vec.NLM_anatomical)
			if (inputScalars.useImages)
				kArgs.emplace_back(&d_urefIm);
			else
				kArgs.emplace_back(&d_uref);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		//Compute the kernel
		status = cuLaunchKernel(kernelNLM, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the NLM kernel\n");
			return status;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after NLM kernel\n");
			return status;
		}
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA NLM gradient computed");
		return 0;
	}

	/// <summary>
	/// Compute relative difference prior (RDP)
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="gamma controls the shape of the prior"></param>
	/// <param name="weights_RDP (UNUSED) the voxel weights for RDP"></param>
	/// <returns></returns>
	inline int computeRDP(const scalarStruct& inputScalars, float gamma, float beta, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA RDP gradient computation");
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior[2] = inputScalars.Nz[0];
		status = cuCtxSynchronize();
		if (DEBUG) {
			mexPrintBase("inputScalars.epps = %.9f\n", inputScalars.epps);
			mexPrintBase("gamma = %f\n", gamma);
			mexPrintBase("inputScalars.Nx = %d\n", inputScalars.Nx[0]);
			mexPrintBase("inputScalars.Ny = %d\n", inputScalars.Ny[0]);
			mexPrintBase("inputScalars.Nz * inputScalars.nRekos = %d\n", inputScalars.Nz[0] * inputScalars.nRekos);
			mexPrintBase("globalPrior[0] = %d\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %d\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %d\n", globalPrior[2]);
			mexEval();
		}
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		if (inputScalars.useImages) {
			kArgs.emplace_back(&d_inputI);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		}
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NOrig);
		kArgs.emplace_back(&gamma);
		kArgs.emplace_back(&apu);
		kArgs.emplace_back(&beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		if (RDPLargeNeighbor) {
			kArgs.emplace_back(&d_weights);
			if (useRDPRef)
				if (inputScalars.useImages)
					kArgs.emplace_back(&d_RDPrefI);
				else
					kArgs.emplace_back(reinterpret_cast<void*>(&d_RDPref));
		}
		// Compute the kernel
		status = cuLaunchKernel(kernelRDP, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the RDP kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after RDP kernel\n");
			return -1;
		}
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			if (RDPLargeNeighbor && useRDPRef) {
				status = cuTexObjectDestroy(d_RDPrefI);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
			}
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA RDP gradient computed");
		return 0;
	}

	/// <summary>
	/// Compute relative generalized Gaussian Markov random field prior (GGMRF)
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="p constant controlling the powers near from the origin"></param>
	/// <param name="q constant controlling the powers distant from the origin"></param>
	/// <param name="c constant controlling the approximate threshold of transition between low and high contrast regions"></param>
	/// <param name="beta regularization parameter"></param>
	/// <returns></returns>
	inline int computeGGMRF(const scalarStruct& inputScalars, float p, float q, float c, float pqc, float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA GGMRF gradient computation");
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior[2] = inputScalars.Nz[0];
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
		if (DEBUG) {
			mexPrintBase("p = %f\n", p);
			mexPrintBase("q = %f\n", q);
			mexPrintBase("c = %f\n", c);
			mexPrintBase("pqc = %f\n", pqc);
			mexPrintBase("inputScalars.Nx = %d\n", inputScalars.Nx[0]);
			mexPrintBase("inputScalars.Ny = %d\n", inputScalars.Ny[0]);
			mexPrintBase("inputScalars.Nz * inputScalars.nRekos = %d\n", inputScalars.Nz[0] * inputScalars.nRekos);
			mexPrintBase("globalPrior[0] = %d\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %d\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %d\n", globalPrior[2]);
			mexEval();
		}
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		if (inputScalars.useImages) {
			kArgs.emplace_back(&d_inputI);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		}
		kArgs.emplace_back(&d_weights);
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&p);
		kArgs.emplace_back(&q);
		kArgs.emplace_back(&c);
		kArgs.emplace_back(&pqc);
		kArgs.emplace_back(&beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		// Compute the kernel
		status = cuLaunchKernel(kernelGGMRF, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the GGMRF kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after GGMRF kernel\n");
			return -1;
		}
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA GGMRF gradient computed");
		return 0;
	}


	inline int ProxHelperQ(float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
		unsigned int kernelIndProxRDP = 0ULL;
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qX));
		kArgs.emplace_back(&alpha);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal RDP helper kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after proximal RDP helper kernel\n");
			return -1;
		}
		return 0;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TV prior
	/// </summary>
	/// <param name="q the input TV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <returns></returns>
	inline int ProxTVHelperQ(float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qY));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qZ));
		kArgs.emplace_back(&alpha);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after kernel\n");
			return -1;
		}
		return 0;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TGV prior
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
	inline int ProxTGVHelperQ(const scalarStruct& inputScalars, float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		status = cuCtxSynchronize();
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rY));
		if (!inputScalars.TGV2D)
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rZ));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rXY));
		if (!inputScalars.TGV2D) {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rYZ));
		}
		kArgs.emplace_back(&alpha);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTGVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after kernel\n");
			return -1;
		}
		return 0;
	}

	/// <summary>
	/// Divergence of the TV prior
	/// </summary>
	/// <param name="im the input array from where the divergence is computed"></param>
	/// <param name="input the backprojection, to which the divergence is added"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
	inline int ProxTVDiv(const scalarStruct& inputScalars) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TV divergence");
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
		if (DEBUG) {
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
			mexEval();
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed before divergence kernel\n");
			return -1;
		}
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NPrior);
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qY));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qZ));
		kArgs.emplace_back(reinterpret_cast<void*>(&vec_opencl.d_rhs_os[0]));
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV divergence kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after divergence kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Proximal TV divergence computed");
		return 0;
	}

	/// <summary>
	/// TV prior (gradient)
	/// </summary>
	/// <param name="im the input image (from which the gradient/TV is computed)"></param>
	/// <param name="input the output TV"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <param name="sigma adjustable constant for some of the priors"></param>
	/// <param name="v divergence of the symmetric derivative for TGV"></param>
	/// <returns></returns>
	inline int ProxTVGrad(const scalarStruct& inputScalars, float sigma2, const size_t vSize) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TV gradient");
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPrior[0]);
			mexPrintBase("global[1] = %u\n", globalPrior[1]);
			mexPrintBase("global[2] = %u\n", globalPrior[2]);
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
			mexPrintBase("vSize = %u\n", vSize);
			mexEval();
		}
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NPrior);
		kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qY));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qZ));
		kArgs.emplace_back(&sigma2);
		if (vSize > 0) {
			if (!inputScalars.TGV2D)
				kArgs.emplace_back(reinterpret_cast<void*>(&d_vZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vX));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vY));
		}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTVGrad, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV gradient kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after gradient kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Proximal TV gradient computed");
		return 0;
	}

	/// <summary>
	/// Symmetric derivative for TGV
	/// </summary>
	/// <param name="v input array"></param>
	/// <param name="q the output symmetric derivative array"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma2 the sigma value of CP/PDHG (1 for PKMA)"></param>
	/// <returns></returns>
	inline int ProxTGVSymmDeriv(const scalarStruct& inputScalars, float sigma2) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TGV symmetric derivative");
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
		unsigned int kernelIndCPTGV = 0ULL;
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
			mexEval();
		}
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NPrior);
		kArgs.emplace_back(reinterpret_cast<void*>(&d_vX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_vY));
		if (!inputScalars.TGV2D)
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vZ));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rY));
		if (!inputScalars.TGV2D) {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXY));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rYZ));
		}
		else
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXY));
		kArgs.emplace_back(&sigma2);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTGVSymmDeriv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV symmetric derivative kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after symmetric derivative kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Proximal TGV symmetric derivative computed");
		return 0;
	}

	/// <summary>
	/// Divergence for TGV
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="v output of the divergence"></param>
	/// <param name="p the TV gradient"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="theta theta value of CP/PDHG or the momentum parameter for PKMA"></param>
	/// <param name="tau tau value of CP/PDHG (1 for PKMA)"></param>
	/// <returns></returns>
	inline int ProxTGVDiv(const scalarStruct& inputScalars, float theta, float tau) {
		if (inputScalars.verbose >= 3) {
			mexPrint("Starting Proximal TGV divergence");
		}
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		status = cuCtxSynchronize();
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NPrior);
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rY));
		if (!inputScalars.TGV2D) {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXY));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXZ));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rYZ));
		}
		else
			kArgs.emplace_back(reinterpret_cast<void*>(&d_rXY));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_vX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_vY));
		if (!inputScalars.TGV2D)
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vZ));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qX));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qY));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_qZ));
		kArgs.emplace_back(&theta);
		kArgs.emplace_back(&tau);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		// Compute the kernel
		status = cuLaunchKernel(kernelProxTGVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV divergence kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after divergence kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Proximal TGV divergence complete");
		return 0;
	}

	/// <summary>
	/// In-place element-wise computations, both multiplication and division supported, for either 1D or 2D arrays
	/// </summary>
	/// <param name="vector input array"></param>
	/// <param name="input input and output array"></param>
	/// <param name="mult if true, performs multiplication, otherwise division"></param>
	/// <param name="D2 if true, assumes 2D case, otherwise 1D"></param>
	/// <returns></returns>
	inline int elementWiseComp(const bool mult, const uint64_t size[], bool D2 = false) {
		const unsigned int gSize[3] = {static_cast<unsigned int>(size[0]), static_cast<unsigned int>(size[1]), static_cast<unsigned int>(size[2])};
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		unsigned char D = static_cast<unsigned char>(D2);
		if (DEBUG) {
			mexPrintBase("gSize[0] = %u\n", gSize[0]);
			mexPrintBase("gSize[1] = %u\n", gSize[1]);
			mexPrintBase("gSize[2] = %u\n", gSize[2]);
			mexEval();
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to synchronize before element-wise kernel\n");
			return -1;
		}
		if (mult) {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vector));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_input));
			kArgs.emplace_back(&D);
			// Compute the kernel
			status = cuLaunchKernel(kernelElementMultiply, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_vector));
			kArgs.emplace_back(reinterpret_cast<void*>(&d_input));
			// Compute the kernel
			status = cuLaunchKernel(kernelElementDivision, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		}
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the element-wise kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after element-wise kernel\n");
			return -1;
		}
		return 0;
	}

	/// <summary>
	/// The gradient of hyperbolic prior
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma adjustable weighting parameter"></param>
	/// <param name="beta regularization parameter"></param>
	/// <returns></returns>
	inline int hyperGradient(const scalarStruct& inputScalars, float sigma, float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA hyperbolic prior gradient computation");
		CUresult status = CUDA_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior[2] = inputScalars.Nz[0];
		status = cuCtxSynchronize();
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		std::vector<void*> kArgs;
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		if (inputScalars.useImages) {
			kArgs.emplace_back(&d_inputI);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		}
		float smooth = inputScalars.epps;
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NOrig);
		kArgs.emplace_back(&sigma);
		kArgs.emplace_back(&smooth);
		kArgs.emplace_back(&beta);
		kArgs.emplace_back(&d_weights);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		// Compute the kernel
		status = cuLaunchKernel(kernelHyper, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the hyperbolic prior gradient kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after hyperbolic prior gradient kernel\n");
			return -1;
		}
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA hyperbolic prior gradient computed");
		return 0;
	}

	/// <summary>
	/// The gradient of TV prior
	/// </summary>
	/// <param name="grad output gradient array"></param>
	/// <param name="im input image (from which the gradient is computed)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma various adjustable parameters for some of the priors"></param>
	/// <param name="smooth smoothing value that allows differentiation"></param>
	/// <returns></returns>
	inline int TVGradient(const scalarStruct& inputScalars, float sigma, float smooth, float beta, float C = 0.f, const int type = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA TV gradient computation");
		CUresult status = CUDA_SUCCESS;
		status = cuCtxSynchronize();
		if (inputScalars.largeDim)
			globalPrior[2] =  inputScalars.Nz[0];
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("smooth = %f\n", smooth);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		std::vector<void*> kArgs;
		kArgs.emplace_back(reinterpret_cast<void*>(&d_W));
		if (inputScalars.useImages) {
			kArgs.emplace_back(&d_inputI);
		}
		else {
			kArgs.emplace_back(reinterpret_cast<void*>(&d_inputB));
		}
		kArgs.emplace_back(&d_N[0]);
		kArgs.emplace_back(&d_NOrig);
		kArgs.emplace_back(&sigma);
		kArgs.emplace_back(&smooth);
		kArgs.emplace_back(&beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			kArgs.emplace_back(&d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kArgs.emplace_back(&d_eFOVIndices);
		if (type == 2 || type == 3)
			kArgs.emplace_back(&C);
		if (type > 0)
			kArgs.emplace_back(reinterpret_cast<void*>(&d_refIm));
		// Compute the kernel
		status = cuLaunchKernel(kernelTV, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the TV gradient kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after TV gradient kernel\n");
			return -1;
		}
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA TV gradient computed");
		return 0;
	}


	inline int PoissonUpdate(const scalarStruct& inputScalars, float lambda, float epps, float alpha, const int ii = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA Poisson update (PKMA/MBSREM/BSREM) computation");
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		status = cuCtxSynchronize();
		global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localPrior[0];
		global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localPrior[1];
		global[2] = inputScalars.Nz[ii];
		bool apu = inputScalars.enforcePositivity;
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].z);
			mexPrintBase("lambda = %.8f\n", lambda);
			mexPrintBase("alpha = %f\n", alpha);
			mexEval();
		}
		kArgs.emplace_back(reinterpret_cast<void*>(&d_im));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rhs));
		kArgs.emplace_back(&d_N[ii]);
		kArgs.emplace_back(&lambda);
		kArgs.emplace_back(&epps);
		kArgs.emplace_back(&alpha);
		kArgs.emplace_back(&apu);
		// Compute the kernel
		status = cuLaunchKernel(kernelPoisson, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Poisson update kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after Poisson update kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA Poisson update computed");
		return 0;
	}

	inline int PDHGUpdate(const scalarStruct& inputScalars, float epps, float theta, float tau, const int ii = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting CUDA PDHG update computation");
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localPrior[0];
		global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localPrior[1];
		global[2] = inputScalars.Nz[ii];
		bool apu = inputScalars.enforcePositivity;
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].z);
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		kArgs.emplace_back(reinterpret_cast<void*>(&d_im));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_rhs));
		kArgs.emplace_back(reinterpret_cast<void*>(&d_U));
		kArgs.emplace_back(&d_N[ii]);
		kArgs.emplace_back(&epps);
		kArgs.emplace_back(&theta);
		kArgs.emplace_back(&tau);
		kArgs.emplace_back(&apu);
		// Compute the kernel
		status = cuLaunchKernel(kernelPDHG, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the PDHG update kernel\n");
			return -1;
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after PDHG update kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("CUDA PDHG update computed");
		return 0;
	}

	inline int transferTex(const scalarStruct& inputScalars, CUdeviceptr* input, const bool RDP = false) {

		CUresult status = CUDA_SUCCESS;
		CUDA_TEXTURE_DESC texDesc;
		CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
		CUDA_RESOURCE_DESC resDesc;
		CUDA_RESOURCE_VIEW_DESC viewDesc;
		std::memset(&texDesc, 0, sizeof(texDesc));
		std::memset(&resDesc, 0, sizeof(resDesc));
		std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
		std::memset(&viewDesc, 0, sizeof(viewDesc));
		arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
		arr3DDesc.NumChannels = 1;
		arr3DDesc.Height = inputScalars.Nx[0];
		arr3DDesc.Width = inputScalars.Ny[0];
		arr3DDesc.Depth = inputScalars.Nz[0];
		status = cuArray3DCreate(&imArray, &arr3DDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to create NLM image array\n");
			return -1;
		}
		CUDA_MEMCPY3D cpy3d;
		std::memset(&cpy3d, 0, sizeof(cpy3d));
		cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
		cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(input);
		cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(float);
		cpy3d.srcHeight = inputScalars.Nx[0];
		cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		cpy3d.dstArray = imArray;
		cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(float);
		cpy3d.Height = inputScalars.Nx[0];
		cpy3d.Depth = inputScalars.Nz[0];
		status = cuMemcpy3D(&cpy3d);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to copy NLM image array\n");
			return -1;
		}
		resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
		resDesc.res.array.hArray = imArray;
		texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
		viewDesc.height = inputScalars.Nx[0];
		viewDesc.width = inputScalars.Ny[0];
		viewDesc.depth = inputScalars.Nz[0];
		viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
		if (RDP)
			status = cuTexObjectCreate(&d_RDPrefI, &resDesc, &texDesc, &viewDesc);
		else
			status = cuTexObjectCreate(&d_inputI, &resDesc, &texDesc, &viewDesc);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("NLM image copy failed\n");
			return -1;
		}
		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Synchronization failed\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Synchronization completed\n");
		return 0;
	}

};