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
/// Class object for forward and backward projections. OpenCL version
/// </summary>
class ProjectorClass {
//private:
	// Local size
	size_t local_size[3];
	size_t local_sizePrior[3];
	// Kernel input indices
	cl_uint kernelInd_MRAMLA = 0;
	cl_uint kernelIndFP = 0;
	cl_uint kernelIndBP = 0;
	cl_uint kernelIndFPSubIter = 0;
	cl_uint kernelIndBPSubIter = 0;
	cl_uint kernelIndSens = 0;
	// Crystal pitch
	cl_float2 dPitch;
	// Image dimensions
	cl_int3 d_NOrig, d_NPrior;
	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusPrior[3];
	size_t erotusPriorEFOV[3];
	size_t erotusSens[3];
	// Local and global sizes
	cl::NDRange local, global, localPrior, globalPrior, globalPriorEFOV;

	bool constantBuffer = false;

	// Get the OpenCL context for the current platform
	cl_int clGetPlatformsContext(const uint32_t platform, cl::Context& context, std::vector<cl::CommandQueue>& commandQueues, const std::vector<uint32_t>& usedDevices, std::vector<cl::Device>& devices) {
		cl_int status = CL_SUCCESS;

		// Get the number of platforms 
		std::vector<cl::Platform> platforms;
		status = cl::Platform::get(&platforms);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (DEBUG) {
			mexPrintBase("platforms.size() = %u\n", platforms.size());
			mexEval();
		}

		if (platforms.size() == 0) {
			std::cerr << "No platforms available!" << std::endl;
			status = -1;
			return status;
		}
		if (platform >= platforms.size()) {
			std::cerr << "The specified platform number is greater than the available platform numbers!" << std::endl;
			status = -1;
			return status;
		}
		if (DEBUG) {
			mexPrintBase("platform = %u\n", platform);
			mexEval();
		}

		// Get context properties from the chosen platform
		cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, reinterpret_cast <cl_context_properties>(platforms[platform]()), 0 };

		// Create context from the chosen platform
		// If a single device was selected (options.cpu_to_gpu_factor = 0), use GPU if possible
		context = cl::Context(CL_DEVICE_TYPE_ALL, properties, NULL, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		// Get device IDs
		std::vector<cl::Device> devices2;
		status = context.getInfo(CL_CONTEXT_DEVICES, &devices2);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		devices.push_back(devices2[usedDevices[0]]);
		if (DEBUG) {
			mexPrintBase("devices.size() = %u\n", devices.size());
			mexEval();
		}

		// Create the command queues
		// Enable out of order execution (devices can compute kernels at the same time)
		for (size_t i = 0; i < devices.size(); i++) {
			commandQueues.push_back(cl::CommandQueue(context, devices[i], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
		if (DEBUG) {
			mexPrintBase("commandQueues.size() = %u\n", commandQueues.size());
			mexEval();
		}

		for (cl_uint i = 0; i < commandQueues.size(); i++) {
			commandQueues[i].finish();
		}

		return status;
	}

	/// <summary>
	/// This function creates the OpenCL programs for the forward and backward projections and for NLM/MRP/RDP/TV
	/// </summary>
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
	/// <param name="programFP the program to store forward projection program"></param>
	/// <param name="programBP the program to store backprojection program"></param>
	/// <param name="programAux the program to store auxliary (such as priors) programs"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="w_vec specifies some of the special options used"></param>
	/// <param name="local_size the local size"></param>
	/// <returns></returns>
	inline cl_int createProgram(cl::Context& CLContext, cl::Device& CLDeviceID, cl::Program& programFP, cl::Program& programBP,
		cl::Program& programAux, cl::Program& programSens, const char* header_directory, scalarStruct& inputScalars, const RecMethods MethodList,
		const Weighting& w_vec, const size_t local_size[], const int type = -1) {

		cl_int status = CL_SUCCESS;

		std::string kernelFile = header_directory;
		std::string kernel_path, kernel_pathBP;
		std::string contentFP, contentBP;
		std::string contentAux;
		std::string options = "-cl-single-precision-constant";
		options += " -DOPENCL";
		if (inputScalars.useMAD) {
			options += " -cl-fast-relaxed-math";
			options += " -DUSEMAD";
		}
		if ((inputScalars.useImages && inputScalars.FPType != 4 && inputScalars.FPType != 5 && inputScalars.BPType != 5) || (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 5)) {
			options += " -DUSEIMAGES";
			inputScalars.useBuffers = false;
		}
		std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
		// Load the header text file
		std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
		// Load orthogonal/volume of intersection headers if applicable
		if (inputScalars.FPType == 2 || inputScalars.BPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 3) {
			if (inputScalars.orthXY)
				options += " -DCRYSTXY";
			if (inputScalars.orthZ)
				options += " -DCRYSTZ";
			std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
			std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
			contentHeader += contentHeader3;
		}

		kernel_path = kernelFile;
		kernel_pathBP = kernelFile;
		if (inputScalars.FPType > 0 && inputScalars.FPType != 6) {
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
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
		if (inputScalars.BPType > 0 && inputScalars.BPType != 6) {
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
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
		if (constantBuffer || (inputScalars.listmode > 0 && !inputScalars.indexBased))
			options += " -DUSEGLOBAL";
		if (inputScalars.raw == 1)
			options += " -DRAW";
		if (inputScalars.maskFP) {
			options += " -DMASKFP";
			if (inputScalars.maskFPZ > 1)
				options += " -DMASKFP3D";
		}
		if (inputScalars.maskBP) {
			options += " -DMASKBP";
			if (inputScalars.maskBPZ > 1)
				options += " -DMASKBP3D";
		}
		if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights)
			options += " -DFDK";
		if (inputScalars.offset)
			options += " -DOFFSET";
		if (inputScalars.attenuation_correction == 1u && inputScalars.CTAttenuation)
			options += " -DATN";
		else if (inputScalars.attenuation_correction == 1u && !inputScalars.CTAttenuation)
			options += " -DATNM";
		if (inputScalars.normalization_correction == 1u)
			options += " -DNORM";
		if (inputScalars.scatter == 1u)
			options += " -DSCATTER";
		if (inputScalars.randoms_correction == 1u)
			options += " -DRANDOMS";
		if (inputScalars.nLayers > 1U)
			if (inputScalars.listmode > 0 && inputScalars.indexBased)
				options += (" -DNLAYERS=" + std::to_string(inputScalars.nLayers));
			else
			options += (" -DNLAYERS=" + std::to_string(inputScalars.nProjections / (inputScalars.nLayers * inputScalars.nLayers)));
		if (inputScalars.TOF) {
			options += " -DTOF";
		}
		if (inputScalars.CT)
			options += " -DCT";
		else if (inputScalars.PET)
			options += " -DPET";
		else if (inputScalars.SPECT) {
			options += " -DSPECT";
			options += (" -DN_RAYS=" + std::to_string(inputScalars.n_rays * inputScalars.n_rays3D));
			options += (" -DN_RAYS2D=" + std::to_string(inputScalars.n_rays));
			options += (" -DN_RAYS3D=" + std::to_string(inputScalars.n_rays3D));
		}

		options += (" -DNBINS=" + std::to_string(inputScalars.nBins));
		if (inputScalars.listmode == 1)
			options += " -DLISTMODE";
		else if (inputScalars.listmode == 2)
			options += " -DLISTMODE2";
		if (inputScalars.listmode > 0 && inputScalars.indexBased)
			options += " -DINDEXBASED";
		if (siddonVal && ((inputScalars.n_rays * inputScalars.n_rays3D) > 1)) {
			options += (" -DN_RAYS=" + std::to_string(inputScalars.n_rays * inputScalars.n_rays3D));
			options += (" -DN_RAYS2D=" + std::to_string(inputScalars.n_rays));
			options += (" -DN_RAYS3D=" + std::to_string(inputScalars.n_rays3D));
		}
		if (inputScalars.pitch)
			options += " -DPITCH";
		if (((inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
			options += " -DSUBSETS";
		if (local_size[1] > 0ULL) {
			options += (" -DLOCAL_SIZE=" + std::to_string(local_size[0]));
			options += (" -DLOCAL_SIZE2=" + std::to_string(local_size[1]));
		}
		else
			options += (" -DLOCAL_SIZE=" + std::to_string(local_size[0]));
		if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
			options += (" -DSTYPE=" + std::to_string(inputScalars.subsetType));
			options += (" -DNSUBSETS=" + std::to_string(inputScalars.subsets));
		}
		if (DEBUG) {
			mexPrintBase("path = %s\n", kernel_path.c_str());
			mexPrintBase("pathBP = %s\n", kernel_pathBP.c_str());
			mexPrintBase("file = %s\n", kernelFile.c_str());
			mexPrintBase("inputScalars.BPType = %u\n", inputScalars.BPType);
			mexPrintBase("inputScalars.FPType = %u\n", inputScalars.FPType);
			mexEval();
		}
		// Build projector program
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			std::string os_options = options;
			os_options += " -DSIDDON";
			os_options += " -DATOMICF";
			std::string os_optionsFP = os_options;
			os_optionsFP += " -DFP";
			if (inputScalars.FPType == 3)
				os_optionsFP += " -DVOL";
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3)
				os_optionsFP += " -DORTH";
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build FP 1-3 program\n");
				}
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
				if (status == CL_SUCCESS && DEBUG) {
					mexPrint("FP 1-3 program built\n");
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build BP 1-3 program\n");
				}
				os_options += " -DBP";
				if (inputScalars.BPType == 3)
					os_options += " -DVOL";
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
					os_options += " -DORTH";
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
				if (status == CL_SUCCESS && DEBUG) {
					mexPrint("BP 1-3 program built\n");
				}
			}
		}
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
			std::string os_options = options;
			if (inputScalars.FPType == 4)
				os_options += " -DFP";
			if (inputScalars.BPType == 4 && inputScalars.CT)
				os_options += " -DBP";
			if (inputScalars.CT) {
				inputScalars.atomic_64bit = false;
				inputScalars.atomic_32bit = false;
			}
			os_options += " -DPTYPE4";
			if (!inputScalars.largeDim)
				os_options += (" -DNVOXELS=" + std::to_string(NVOXELS));
			if (inputScalars.FPType == 4)
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			if (!inputScalars.CT && inputScalars.BPType == 4) {
				os_options = options;
				os_options += " -DPTYPE4";
				os_options += " -DBP";
				os_options += " -DATOMICF";
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			}
			else if (inputScalars.CT && inputScalars.BPType == 4 && inputScalars.FPType != 4)
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
			std::string os_options = options;
			os_options += " -DPROJ5";
			if (inputScalars.meanFP)
				os_options += " -DMEANDISTANCEFP";
			else if (inputScalars.meanBP)
				os_options += " -DMEANDISTANCEBP";
			if (inputScalars.FPType == 5)
				os_options += " -DFP";
			if (inputScalars.BPType == 5)
				os_options += " -DBP";
			if (inputScalars.pitch)
				os_options += (" -DNVOXELS5=" + std::to_string(1));
			else
				os_options += (" -DNVOXELS5=" + std::to_string(NVOXELS5));
			os_options += (" -DNVOXELSFP=" + std::to_string(NVOXELSFP));
			if (inputScalars.FPType == 5)
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			else
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
		}
		if (inputScalars.computeSensImag) {
			std::string os_options = options;
			os_options += " -DBP";
			os_options += " -DATOMICF";
			os_options += " -DSENS";
			if (inputScalars.BPType == 4) {
				os_options += " -DPTYPE4";
				os_options += (" -DNVOXELS=" + std::to_string(NVOXELS));
			}
			else
				os_options += " -DSIDDON";
			status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programSens, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
		}
		// Build prior programs
		if (MethodList.NLM || MethodList.MRP || MethodList.RDP || w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]
			|| MethodList.TV || MethodList.APLS || MethodList.hyperbolic || MethodList.ProxTV || MethodList.ProxTGV || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM ||
			MethodList.CPType || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF || type == 0) {
			if (DEBUG) {
				mexPrint("Building aux programs\n");
			}
			std::string auxKernelPath = kernelFile + "auxKernels.cl";
			std::ifstream sourceFileAux(auxKernelPath.c_str());
			std::string contentAAux((std::istreambuf_iterator<char>(sourceFileAux)), std::istreambuf_iterator<char>());
			contentAux = contentHeader + contentAAux;
			options = "-cl-single-precision-constant";
			options += " -DOPENCL";
			if (inputScalars.useMAD) {
				options += " -cl-fast-relaxed-math";
				options += " -DUSEMAD";
			}
			
			if (inputScalars.useImages)
				options += " -DUSEIMAGES";
			if (inputScalars.useExtendedFOV)
				options += " -DEFOV";
			if (inputScalars.use64BitIndices) {
				options += " -DLTYPE=long";
				options += " -DLTYPE3=long3";
			}
			if (type == 2) {
				if (inputScalars.use_psf)
					options += " -DPSF";
			}
			else if (type == 0) {
				if (inputScalars.CT)
					options += " -DCT";
				if (inputScalars.randoms_correction)
					options += " -DRANDOMS";
				if (inputScalars.use_psf)
					options += " -DPSF";
			}
			else
				options += " -DAF";
			if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
				options += " -DMASKPRIOR";
				if (inputScalars.maskBPZ > 1)
					options += " -DMASKBP3D";
			}
			if (inputScalars.eFOV)
				options += " -DEFOVZ";
			if (MethodList.MRP) {
				options += " -DMEDIAN";
				options += (" -DSEARCH_WINDOW_X=" + std::to_string(w_vec.Ndx));
				options += (" -DSEARCH_WINDOW_Y=" + std::to_string(w_vec.Ndy));
				options += (" -DSEARCH_WINDOW_Z=" + std::to_string(w_vec.Ndz));
			}
			if (MethodList.NLM) {
				options += " -DNLM_";
				if (w_vec.NLM_MRP)
					options += (" -DNLTYPE=" + std::to_string(2));
				else if (w_vec.NLTV)
					options += (" -DNLTYPE=" + std::to_string(1));
				else if (w_vec.NLRD)
					options += (" -DNLTYPE=" + std::to_string(3));
				else if (w_vec.NLLange)
					options += (" -DNLTYPE=" + std::to_string(4));
				else if (w_vec.NLLangeFiltered)
					options += (" -DNLTYPE=" + std::to_string(5));
				else if (w_vec.NLGGMRF)
					options += (" -DNLTYPE=" + std::to_string(6));
				else
					options += (" -DNLTYPE=" + std::to_string(0));
				if (w_vec.NLAdaptive)
					options += " -DNLMADAPTIVE";
				if (w_vec.NLM_anatomical)
					options += " -DNLMREF";
				options += (" -DSWINDOWX=" + std::to_string(w_vec.Ndx));
				options += (" -DSWINDOWY=" + std::to_string(w_vec.Ndy));
				options += (" -DSWINDOWZ=" + std::to_string(w_vec.Ndz));
				options += (" -DPWINDOWX=" + std::to_string(w_vec.Nlx));
				options += (" -DPWINDOWY=" + std::to_string(w_vec.Nly));
				options += (" -DPWINDOWZ=" + std::to_string(w_vec.Nlz));
			}
			if (MethodList.GGMRF) {
				options += " -DGGMRF";
				options += (" -DSWINDOWX=" + std::to_string(w_vec.Ndx));
				options += (" -DSWINDOWY=" + std::to_string(w_vec.Ndy));
				options += (" -DSWINDOWZ=" + std::to_string(w_vec.Ndz));
			}
			if (MethodList.hyperbolic) {
				options += " -DHYPER";
				options += (" -DSWINDOWX=" + std::to_string(w_vec.Ndx));
				options += (" -DSWINDOWY=" + std::to_string(w_vec.Ndy));
				options += (" -DSWINDOWZ=" + std::to_string(w_vec.Ndz));
			}
			if (MethodList.RDP) {
				options += " -DRDP";
				if (w_vec.RDPLargeNeighbor) {
					options += " -DRDPCORNERS";
					options += (" -DSWINDOWX=" + std::to_string(w_vec.Ndx));
					options += (" -DSWINDOWY=" + std::to_string(w_vec.Ndy));
					options += (" -DSWINDOWZ=" + std::to_string(w_vec.Ndz));
				}
				if (w_vec.RDP_anatomical)
					options += " -DRDPREF";
			}
			if (MethodList.ProxRDP && w_vec.RDPLargeNeighbor)
				options += " -DRDPCORNERS";
			if (MethodList.TV && !w_vec.data.TV_use_anatomical) {
				options += " -DTVGRAD";
				if (w_vec.data.TVtype == 6)
					options += " -DTVW1";
				else if (w_vec.data.TVtype == 4)
					options += " -DSATV";
				else if (w_vec.data.TVtype == 2)
					options += " -DJPTV";
				if (w_vec.derivType > 0)
					options += (" -DDIFFTYPE=" + std::to_string(w_vec.derivType));
			}
			else if ((MethodList.TV && w_vec.data.TV_use_anatomical) || MethodList.APLS) {
				options += " -DTVGRAD";
				if (w_vec.data.TVtype == 1)
					options += " -DANATOMICAL1";
				else if (w_vec.data.TVtype == 2)
					options += " -DANATOMICAL2";
				else if (w_vec.data.TVtype == 5)
					options += " -DANATOMICAL3";
				if (w_vec.derivType > 0)
					options += (" -DDIFFTYPE=" + std::to_string(w_vec.derivType));
			}
			if (MethodList.ProxTV) {
				options += " -DPROXTV";
				if (w_vec.UseL2Ball)
					options += " -DL2";
				if (w_vec.derivType > 0)
					options += (" -DDIFFTYPE=" + std::to_string(w_vec.derivType));
			}
			if (MethodList.ProxTGV || MethodList.TGV) {
				options += " -DPROXTV";
				options += " -DPROXTGV";
				if (w_vec.UseL2Ball)
					options += " -DL2";
				if (w_vec.derivType > 0)
					options += (" -DDIFFTYPE=" + std::to_string(w_vec.derivType));
				if (!inputScalars.TGV2D)
					options += " -DTGVZ";
			}
			if (MethodList.ProxRDP)
				options += " -DPROXRDP";
			if (MethodList.ProxNLM) {
				options += " -DPROXNLM";
				if (w_vec.NLM_MRP)
					options += (" -DNLTYPE=" + std::to_string(2));
				else if (w_vec.NLTV)
					options += (" -DNLTYPE=" + std::to_string(1));
				else if (w_vec.NLRD)
					options += (" -DNLTYPE=" + std::to_string(3));
				else if (w_vec.NLLange)
					options += (" -DNLTYPE=" + std::to_string(4));
				else if (w_vec.NLLangeFiltered)
					options += (" -DNLTYPE=" + std::to_string(5));
				else
					options += (" -DNLTYPE=" + std::to_string(0));
				if (w_vec.NLM_anatomical)
					options += " -DNLMREF";
				options += (" -DSWINDOWX=" + std::to_string(w_vec.Ndx));
				options += (" -DSWINDOWY=" + std::to_string(w_vec.Ndy));
				options += (" -DSWINDOWZ=" + std::to_string(w_vec.Ndz));
				options += (" -DPWINDOWX=" + std::to_string(w_vec.Nlx));
				options += (" -DPWINDOWY=" + std::to_string(w_vec.Nly));
				options += (" -DPWINDOWZ=" + std::to_string(w_vec.Nlz));
			}
			if (local_sizePrior[1] > 0ULL) {
				options += (" -DLOCAL_SIZE=" + std::to_string(local_sizePrior[0]));
				options += (" -DLOCAL_SIZE2=" + std::to_string(local_sizePrior[1]));
				options += (" -DLOCAL_SIZE3=" + std::to_string(local_sizePrior[2]));
			}
			else
				options += (" -DLOCAL_SIZE=" + std::to_string(local_sizePrior[0]));
			if (MethodList.PKMA)
				options += " -DPKMA";
			else if (MethodList.MBSREM || MethodList.MRAMLA)
				options += " -DMBSREM";
			else if (MethodList.BSREM || MethodList.RAMLA)
				options += " -DBSREM";
			else if (MethodList.CPType) {
				options += " -DPDHG";
				if (inputScalars.subsets > 1)
					options += " -DSUBSETS";
			}
			status = buildProgram(inputScalars.verbose, contentAux, CLContext, CLDeviceID, programAux, inputScalars.atomic_64bit, inputScalars.atomic_32bit, options);
		}
		if (DEBUG) {
			mexPrintBase("status = %u\n", status);
			mexPrintBase("w_vec.NLM_MRP = %u\n", w_vec.NLM_MRP);
			mexPrintBase("w_vec.NLTV = %u\n", w_vec.NLTV);
			mexPrintBase("w_vec.NLRD = %u\n", w_vec.NLRD);
			mexPrintBase("w_vec.NLLange = %u\n", w_vec.NLLange);
			mexPrintBase("w_vec.NLLangeFiltered = %u\n", w_vec.NLLangeFiltered);
			mexEval();
		}
		return status;
	}

	/// <summary>
	/// Builds the input OpenCL program
	/// </summary>
	/// <param name="verbose the level of verbosity"></param>
	/// <param name="contentFP program code"></param>
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
	/// <param name="program the program where to store the built program"></param>
	/// <param name="atomic_64bit are 64-bit (int64) atomics used"></param>
	/// <param name="atomic_32bit are 32-bit (int) atomics used"></param>
	/// <param name="options preprocessor values for the build"></param>
	/// <returns></returns>
	inline cl_int buildProgram(const int8_t verbose, std::string contentFP, cl::Context& CLContext, cl::Device& CLDeviceID, cl::Program& program,
		bool& atomic_64bit, const bool atomic_32bit, std::string options) {
		cl_int status = CL_SUCCESS;
		size_t pituus;
		if (atomic_64bit) {
			pituus = options.length();
			options += " -DCAST=long";
			options += " -DATOMIC";
			options += (" -DTH=" + std::to_string(TH));
		}
		else if (atomic_32bit) {
			options += " -DCAST=int";
			options += " -DATOMIC32";
			options += (" -DTH=" + std::to_string(TH32));
		}
		else {
			options += " -DCAST=float";
		}
		if (DEBUG || verbose >= 2)
			mexPrintBase("%s\n", options.c_str());
		if (atomic_64bit) {
			cl::string apu = CLDeviceID.getInfo<CL_DEVICE_EXTENSIONS>();
			cl::string apu2 = "cl_khr_int64_base_atomics";
			size_t var = apu.find(apu2);
			if (var < 0) {
				options.erase(pituus, options.size() + 1);
				options += " -DCAST=float";
				status = -1;
			}
			else {
				std::vector<std::string> testi;
				testi.push_back(contentFP);
				cl::Program::Sources source(testi);
				program = cl::Program(CLContext, source);
				status = program.build(options.c_str());
				if (status == CL_SUCCESS && (DEBUG || verbose >= 2)) {
					mexPrint("OpenCL program (64-bit atomics) built\n");
				}
				else if (status != CL_SUCCESS) {
					mexPrint("Failed to build 64-bit atomics program.\n");
					if (DEBUG) {
						getErrorString(status);
						std::vector<cl::Device> dev;
						CLContext.getInfo(CL_CONTEXT_DEVICES, &dev);
						for (int ll = 0; ll < dev.size(); ll++) {
							cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
							if (status != CL_BUILD_ERROR)
								continue;
							std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
							std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
							mexPrintBase("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
						}
					}
					options.erase(pituus, options.size() + 1);
					options += " -DCAST=float";
				}
			}
		}
		else
			status = -1;
		// If not, use 32-bit atomic add (float)
		if (status != CL_SUCCESS) {
			status = CL_SUCCESS;
			atomic_64bit = false;
			std::vector<std::string> testi;
			testi.push_back(contentFP);
			cl::Program::Sources source(testi);
			program = cl::Program(CLContext, source);
			status = program.build(options.c_str());
			if (status == CL_SUCCESS && (DEBUG || verbose >= 2)) {
				mexPrint("OpenCL program built\n");
			}
			else if (status != CL_SUCCESS) {
				mexPrint("Failed to build OpenCL program.\n");
				getErrorString(status);
				std::vector<cl::Device> dev;
				CLContext.getInfo(CL_CONTEXT_DEVICES, &dev);
				for (int ll = 0; ll < dev.size(); ll++) {
					cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
					if (status != CL_BUILD_ERROR)
						continue;
					std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
					std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
					mexPrintBase("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
				}
			}
		}
		return status;
	}

	/// <summary>
	/// Creates the necessary OpenCL kernels from the input programs
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
	inline cl_int createKernels(cl::Kernel& kernelFP, cl::Kernel& kernelBP, cl::Kernel& kernelNLM, cl::Kernel& kernelMed,
		cl::Kernel& kernelRDP, cl::Kernel& kernelGGMRF, const cl::Program& programFP, const cl::Program& programBP, const cl::Program& programAux, 
		const cl::Program& programSens, const RecMethods& MethodList, const Weighting& w_vec, const scalarStruct& inputScalars, const int type = -1) {
		cl_int status = CL_SUCCESS;
		// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
			if (inputScalars.FPType == 4) {
				kernelFP = cl::Kernel(programFP, "projectorType4Forward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 4 FP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("OpenCL kernel for projector type 4 FP successfully created\n");
				}
			}
			if (inputScalars.BPType == 4) {
				if (inputScalars.FPType == 4 && inputScalars.CT)
					kernelBP = cl::Kernel(programFP, "projectorType4Backward", &status);
				else if (!inputScalars.CT)
					kernelBP = cl::Kernel(programBP, "projectorType4Forward", &status);
				else
					kernelBP = cl::Kernel(programBP, "projectorType4Backward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 4 BP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("OpenCL kernel for projector type 4 BP successfully created\n");
				}
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
			if (inputScalars.FPType == 5) {
				kernelFP = cl::Kernel(programFP, "projectorType5Forward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 5 FP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("OpenCL kernel for projector type 5 FP successfully created\n");
				}
			}
			if (inputScalars.BPType == 5) {
				if (inputScalars.FPType == 5)
					kernelBP = cl::Kernel(programFP, "projectorType5Backward", &status);
				else //if (inputScalars.BPType == 5)
					kernelBP = cl::Kernel(programBP, "projectorType5Backward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create projector type 5 BP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose >= 2) {
					mexPrint("OpenCL kernel for projector type 5 BP successfully created\n");
				}
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3))
				kernelFP = cl::Kernel(programFP, "projectorType123", &status);
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3)
				kernelBP = cl::Kernel(programBP, "projectorType123", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create OS-methods kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("OpenCL kernel successfully created\n");
			}
		}

		if (MethodList.NLM) {
			kernelNLM = cl::Kernel(programAux, "NLM", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create NLM kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("NLM kernel successfully created\n");
			}
		}
		if (MethodList.MRP) {
			kernelMed = cl::Kernel(programAux, "medianFilter3D", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create Median kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Median kernel successfully created\n");
			}
		}
		if (MethodList.RDP) {
			kernelRDP = cl::Kernel(programAux, "RDPKernel", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create RDP kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("RDP kernel successfully created\n");
			}
		}
		if (MethodList.GGMRF) {
			kernelGGMRF = cl::Kernel(programAux, "GGMRFKernel", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create GGMRF kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("GGMRF kernel successfully created\n");
			}
		}
		if (MethodList.TV || MethodList.APLS) {
			kernelTV = cl::Kernel(programAux, "TVKernel", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create TV kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("TV kernel successfully created\n");
			}
		}
		if (MethodList.hyperbolic) {
			kernelHyper = cl::Kernel(programAux, "hyperbolicKernel", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create hyperbolic prior kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Hyperbolic prior kernel successfully created\n");
			}
		}
		if (MethodList.PKMA || MethodList.BSREM || MethodList.MBSREM || MethodList.MRAMLA || MethodList.RAMLA) {
			kernelPoisson = cl::Kernel(programAux, "PoissonUpdate", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create Poisson Update kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Poisson Update kernel successfully created\n");
			}
		}
		if (MethodList.CPType) {
			kernelPDHG = cl::Kernel(programAux, "PDHGUpdate", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create PDHG Update kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("PDHG Update kernel successfully created\n");
			}
		}
		if (MethodList.ProxTV) {
			kernelProxTVq = cl::Kernel(programAux, "ProxTVq", &status);
			kernelProxTVDiv = cl::Kernel(programAux, "ProxTVDivergence", &status);
			kernelProxTVGrad = cl::Kernel(programAux, "ProxTVGradient", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create CPTV kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("CPTV kernel successfully created\n");
			}
		}
		if (MethodList.ProxRDP) {
			kernelProxq = cl::Kernel(programAux, "Proxq", &status);
			kernelProxRDP = cl::Kernel(programAux, "ProxRDP", &status);
			kernelProxTrans = cl::Kernel(programAux, "ProxTrans", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create proximal RDP kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Proximal RDP kernel successfully created\n");
			}
		}
		if (MethodList.ProxNLM) {
			kernelProxq = cl::Kernel(programAux, "Proxq", &status);
			kernelProxNLM = cl::Kernel(programAux, "ProxNLM", &status);
			kernelProxTrans = cl::Kernel(programAux, "ProxTrans", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create proximal NLM kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Proximal NLM kernel successfully created\n");
			}
		}
		if (MethodList.ProxTGV) {
			kernelProxTVq = cl::Kernel(programAux, "ProxTVq", &status);
			kernelProxTGVq = cl::Kernel(programAux, "ProxTGVq", &status);
			kernelProxTVDiv = cl::Kernel(programAux, "ProxTVDivergence", &status);
			kernelProxTVGrad = cl::Kernel(programAux, "ProxTVGradient", &status);
			kernelProxTGVDiv = cl::Kernel(programAux, "ProxTGVDivergence", &status);
			kernelProxTGVSymmDeriv = cl::Kernel(programAux, "ProxTGVSymmDeriv", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create CPTGV kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("CPTGV kernel successfully created\n");
			}
		}
		if (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) {
			kernelElementMultiply = cl::Kernel(programAux, "vectorElementMultiply", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create element-wise kernels\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Element-wise kernels successfully created\n");
			}
			kernelElementDivision = cl::Kernel(programAux, "vectorElementDivision", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create element-wise kernels\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Element-wise kernels successfully created\n");
			}
		}
		if (type == 0) {
			kernelsumma = cl::Kernel(programAux, "summa", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create implementation 3 kernels\n");
				return -1;
			}
			kernelEstimate = cl::Kernel(programAux, "computeEstimate", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create implementation 3 kernels\n");
				return -1;
			}
			kernelForward = cl::Kernel(programAux, "forward", &status);
			if (inputScalars.use_psf) {
				kernelPSFf = cl::Kernel(programAux, "Convolution3D_f", &status);
				kernelPSF = cl::Kernel(programAux, "Convolution3D", &status);
			}

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create implementation 3 kernels\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose >= 2) {
				mexPrint("Implementation 3 kernels successfully created\n");
			}
		}
		//else if (type == 2) {
		//	if (inputScalars.use_psf) {
		//		kernelPSF = cl::Kernel(programAux, "Convolution3D", &status);
		//		if (status != CL_SUCCESS) {
		//			getErrorString(status);
		//			mexPrint("Failed to create PSF kernel\n");
		//			return -1;
		//		}
		//	}
		//}
		if (inputScalars.computeSensImag) {
			if (inputScalars.BPType == 4)
				kernelSensList = cl::Kernel(programSens, "projectorType4Forward", &status);
			else
				kernelSensList = cl::Kernel(programSens, "projectorType123", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create sensitivity image kernels\n");
				return -1;
			}
		}
		return status;
	}
public:
	cl::Context CLContext;
	std::vector<cl::Device> CLDeviceID;
	std::vector<cl::CommandQueue> CLCommandQueue;
	OpenCL_im_vectors vec_opencl;
	cl::Kernel kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, kernelElementMultiply, kernelElementDivision, 
		kernelTV, kernelProxTGVSymmDeriv, kernelProxTGVDiv, kernelProxTGVq, kernelPoisson, kernelPDHG, kernelProxRDP, kernelProxq, kernelProxTrans, kernelProxNLM, kernelGGMRF,
		kernelsumma, kernelEstimate, kernelPSF, kernelPSFf, kernelDiv, kernelMult, kernelForward, kernelSensList, kernelApu, kernelHyper;
	cl::Buffer d_xcenter, d_ycenter, d_zcenter, d_V, d_TOFCenter, d_output, d_meanBP, d_meanFP, d_eFOVIndices, d_weights, d_inputB, d_W, d_gaussianNLM;
	cl::Image2D d_maskFP, d_maskBP, d_maskPrior;
	cl::Image3D d_maskBP3, d_maskPrior3;
	cl::Image3D d_inputImage, d_attenIm, d_urefIm, d_inputI, d_RDPrefI;
	cl::Buffer d_qX, d_qY, d_qZ;
	cl::Buffer d_rX, d_rY, d_rXY, d_rZ, d_rXZ, d_rYZ;
	cl::Buffer d_vX, d_vY, d_vZ;
	cl::Buffer d_angle;
	cl::Buffer d_vector, d_input;
	cl::Buffer d_im, d_rhs, d_U, d_g, d_uref, d_refIm, d_RDPref;
	cl::Buffer d_outputCT, d_maskBPB, d_attenB;
	cl::Buffer d_rayShiftsDetector, d_rayShiftsSource; // SPECT
	size_t memSize = 0ULL;
	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
	std::vector<cl_float3> b, d, bmax;
	std::vector<cl_int3> d_N;

	std::vector<cl::Image3D> d_maskFP3;
	std::vector<cl::Buffer> d_maskFPB;
	std::vector<cl::Buffer> d_LFull, d_zindexFull, d_xyindexFull, d_normFull, d_scatFull, d_xFull, d_zFull;
	std::vector<cl::Buffer> d_L;
	std::vector<cl::Buffer> d_Summ;
	std::vector<cl::Buffer> d_meas;
	std::vector<cl::Buffer> d_rand;
	std::vector<cl::Buffer> d_imTemp;
	std::vector<cl::Buffer> d_imFinal;
	std::vector<cl::Buffer> d_zindex;
	std::vector<cl::Buffer> d_xyindex;
	std::vector<cl::Buffer> d_trIndex;
	std::vector<cl::Buffer> d_axIndex;
	std::vector<cl::Buffer> d_norm;
	std::vector<cl::Buffer> d_scat;
	std::vector<cl::Buffer> d_x;
	std::vector<cl::Buffer> d_z;
	std::vector<cl::Buffer> d_atten;
	std::vector<cl::Buffer> d_T;
	// Image origin
	cl::detail::size_t_array origin = { { 0, 0, 0 } }; 
	cl::detail::size_t_array region = { { 0, 0, 0 } };
	// Image format
	cl::ImageFormat format;
	cl::ImageFormat formatMask;
	std::vector<std::vector<size_t>> erotusBP, erotusPDHG;
	cl_uchar no_norm = 0;
	int proj6 = 1;
	~ProjectorClass() {	}
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
		local_size[0] = 64ULL;
		local_size[1] = 1ULL;
		local_size[2] = 1ULL;
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || (inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT)))
			local_size[0] = 128ULL;
		if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
			if (inputScalars.nColsD > 1 && !(inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT))) {
				local_size[0] = 16ULL;
				local_size[1] = 16ULL;
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
		local_sizePrior[0] = 16ULL;
		local_sizePrior[1] = 16ULL;
		local_sizePrior[2] = 1ULL;
		cl_int status = CL_SUCCESS;
		proj6 = 0;

		// Create the OpenCL context and command queue and assign the device
#ifdef AF
		CLContext = afcl::getContext(true);
		std::vector<cl::Device> devices = CLContext.getInfo<CL_CONTEXT_DEVICES>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		CLDeviceID.push_back(devices[0]);
		CLCommandQueue.push_back(cl::CommandQueue(afcl::getQueue(true), true));
#else
		status = clGetPlatformsContext(inputScalars.platform, CLContext, CLCommandQueue, inputScalars.usedDevices, CLDeviceID);
#endif
		// For NVIDIA cards, 32 local size seems more optimal with 1D kernelFP
		std::string deviceName = CLDeviceID[0].getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		if (NV.compare(deviceName) == 0 && (inputScalars.projector_type == 1 || inputScalars.projector_type == 11) && local_size[1] == 1ULL)
			local_size[0] = 32ULL;
		if (DEBUG) {
			std::string deviceName2 = CLDeviceID[0].getInfo<CL_DEVICE_NAME>(&status);
			cl_ulong apu = CLDeviceID[0].getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(&status);
			cl_uint apu2 = CLDeviceID[0].getInfo<CL_DEVICE_ADDRESS_BITS>(&status);
			mexPrintBase("CL_DEVICE_MAX_MEM_ALLOC_SIZE = %llu\n", apu);
			mexPrintBase("CL_DEVICE_ADDRESS_BITS = %u\n", apu2);
			mexPrint(deviceName.c_str());
			mexPrint(deviceName2.c_str());
			mexEval();
		}
		cl_ulong constantBufferSize = CLDeviceID[0].getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>(&status);

		if ((inputScalars.size_of_x + inputScalars.size_z) * sizeof(float) >= constantBufferSize)
			constantBuffer = true;

		cl::Program programFP, programBP, programAux, programSens;

		status = createProgram(CLContext, CLDeviceID[0], programFP, programBP, programAux, programSens, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		if (status != CL_SUCCESS) {
			std::cerr << "Error while creating program" << std::endl;
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 2) {
			mexPrint("OpenCL programs successfully created\n");
		}

		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, programSens, MethodList, w_vec, inputScalars, type);
		if (status != CL_SUCCESS) {
			mexPrint("Failed to create kernels\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 2) {
			mexPrint("OpenCL kernels successfully created\n");
		}
		format.image_channel_order = CL_A;
		format.image_channel_data_type = CL_FLOAT;
		formatMask.image_channel_order = CL_A;
		formatMask.image_channel_data_type = CL_UNSIGNED_INT8;

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
			globalPriorEFOV = { inputScalars.NxPrior + erotusPriorEFOV[0], inputScalars.NyPrior + erotusPriorEFOV[1], inputScalars.NzPrior + erotusPriorEFOV[2] };
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
		local = { local_size[0] , local_size[1] };
		localPrior = { local_sizePrior[0] , local_sizePrior[1], local_sizePrior[2] };
		erotusPrior[0] = inputScalars.Nx[0] % local_sizePrior[0];
		erotusPrior[1] = inputScalars.Ny[0] % local_sizePrior[1];
		erotusPrior[2] = inputScalars.Nz[0] % local_sizePrior[2];
		if (erotusPrior[0] > 0)
			erotusPrior[0] = (local_sizePrior[0] - erotusPrior[0]);
		if (erotusPrior[1] > 0)
			erotusPrior[1] = (local_sizePrior[1] - erotusPrior[1]);
		if (erotusPrior[2] > 0)
			erotusPrior[2] = (local_sizePrior[1] - erotusPrior[2]);
		globalPrior = { inputScalars.Nx[0] + erotusPrior[0], inputScalars.Ny[0] + erotusPrior[1], inputScalars.Nz[0] + erotusPrior[2] };
		d_NOrig = { static_cast<cl_int>(inputScalars.NxOrig), static_cast<cl_int>(inputScalars.NyOrig), static_cast<cl_int>(inputScalars.NzOrig) };
		d_NPrior = { static_cast<cl_int>(inputScalars.NxPrior), static_cast<cl_int>(inputScalars.NyPrior), static_cast<cl_int>(inputScalars.NzPrior) };
		dPitch = { w_vec.dPitchX, w_vec.dPitchY };
		b.resize(inputScalars.nMultiVolumes + 1);
		d.resize(inputScalars.nMultiVolumes + 1);
		d_N.resize(inputScalars.nMultiVolumes + 1);
		bmax.resize(inputScalars.nMultiVolumes + 1);
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			b[ii] = { inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii] };
			d[ii] = { inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii] };
			d_N[ii] = { static_cast<cl_int>(inputScalars.Nx[ii]), static_cast<cl_int>(inputScalars.Ny[ii]), static_cast<cl_int>(inputScalars.Nz[ii]) };
			bmax[ii] = { static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
				static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii],
				static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii] };
		}
		if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
			erotusSens[0] = inputScalars.det_per_ring % local_size[0];
			erotusSens[1] = inputScalars.det_per_ring % local_size[1];
			if (erotusSens[1] > 0)
				erotusSens[1] = (local_size[1] - erotusSens[1]);
			if (erotusSens[0] > 0)
				erotusSens[0] = (local_size[0] - erotusSens[0]);
		}
		region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
		return 0;
	}

	/// <summary>
	/// This function first creates the necessary OpenCL buffers and then writes the data to them
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
	inline cl_int createAndWriteBuffers(const std::vector<int64_t>& length, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* L,const int64_t* pituus, const float* atten, const float* norm, const float* extraCorr, 
		const scalarStruct& inputScalars, const Weighting& w_vec, const RecMethods& MethodList) {
		cl_int status = CL_SUCCESS;
		size_t vecSize = 1;
		cl::size_type imX = inputScalars.Nx[0];
		cl::size_type imY = inputScalars.Ny[0];
		cl::size_type imZ = inputScalars.Nz[0];
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (inputScalars.useImages)
				d_urefIm = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			else
				d_uref = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (MethodList.NLM || MethodList.RDP || (MethodList.TV && !w_vec.data.TV_use_anatomical) || MethodList.GGMRF || MethodList.hyperbolic) {
			if (inputScalars.useImages) {
				d_inputI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
				if (status != 0) {
					getErrorString(status);
					mexPrint("Failed to create prior image\n");
					return -1;
				}
			}
		}
		if (MethodList.RDP && w_vec.RDPLargeNeighbor && w_vec.RDP_anatomical) {
			if (inputScalars.useImages) {
				d_RDPrefI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
				if (status != 0) {
					getErrorString(status);
					mexPrint("Failed to create RDP reference image\n");
					return -1;
				}
			}
		}
		// Create the necessary buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			d_weights = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
			imX = inputScalars.Nx[0];
			imY = inputScalars.Ny[0];
			imZ = inputScalars.maskBPZ;
			if (DEBUG) {
				mexPrintBase("imX = %u\n", imX);
				mexPrintBase("imY = %u\n", imY);
				mexPrintBase("imZ = %u\n", imZ);
				mexEval();
			}
			if (imZ > 1)
				d_maskPrior3 = cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, imZ, 0, 0, NULL, &status);
			else
				d_maskPrior = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				d_V = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_V, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			// Detector coordinates
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				d_x[0] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.nLayers > 1)
					d_xFull.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x, NULL, &status));
				//d_xFull.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x / 2, NULL, &status));
				else
					d_xFull.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x, NULL, &status));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			// Attenuation data for image-based attenuation
			if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				imZ = inputScalars.Nz[0];
				if (inputScalars.useBuffers)
					d_attenB = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0], NULL, &status);
				else
					d_attenIm = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.useBuffers) {
					if (inputScalars.maskFP) {
						if (inputScalars.maskFPZ > 1) {
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
								d_maskFPB.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk], NULL, &status));
						}
						else
							d_maskFPB.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint8_t)* inputScalars.nRowsD * inputScalars.nColsD, NULL, &status));
					}
					if (inputScalars.maskBP)
						d_maskBPB = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ, NULL, &status);
				}
				else {
					if (inputScalars.maskFP) {
						imX = inputScalars.nRowsD;
						imY = inputScalars.nColsD;
						imZ = inputScalars.maskFPZ;
						if (imZ > 1) {
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
								d_maskFP3.emplace_back(cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, length[kk], 0, 0, NULL, &status));
						}
						else
							d_maskFP = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
					}
					if (inputScalars.maskBP) {
						imX = inputScalars.Nx[0];
						imY = inputScalars.Ny[0];
						imZ = inputScalars.maskBPZ;
						if (DEBUG) {
							mexPrintBase("imX = %u\n", imX);
							mexPrintBase("imY = %u\n", imY);
							mexPrintBase("imZ = %u\n", imZ);
							mexEval();
						}
						if (imZ > 1)
							d_maskBP3 = cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, imZ, 0, 0, NULL, &status);
						else
							d_maskBP = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
					}
				}
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.SPECT) {
				d_rayShiftsDetector = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				d_rayShiftsSource = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.eFOV) {
				d_eFOVIndices = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nz[0], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				d_angle = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			// TOF bin centers
			if (inputScalars.TOF) {
				d_TOFCenter = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nBins, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.nLayers > 1)
					d_zFull.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z, NULL, &status));
				else
					d_zFull.emplace_back(cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z, NULL, &status));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if (DEBUG) {
					mexPrintBase("length[kk] = %u\n", length[kk]);
					mexPrintBase("kk = %u\n", kk);
					mexPrintBase("vecSize = %u\n", vecSize);
					mexEval();
				}
				if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode != 1) {
					if (inputScalars.pitch)
						d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 6, NULL, &status);
					else
						d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else {
					if (inputScalars.PET && inputScalars.listmode == 0)
						if (inputScalars.nLayers > 1)
							d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 3, NULL, &status);
						else
							d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2, NULL, &status);
					else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased))
						d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z, NULL, &status);
					else
						d_z[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					d_T[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if ((inputScalars.CT || inputScalars.SPECT) || (inputScalars.listmode > 0 && !inputScalars.indexBased)) {
					if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF || ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0))
						d_x[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 6, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
					d_norm[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) {
					d_scat[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					d_atten[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
				if (inputScalars.raw && inputScalars.listmode != 1) {
					d_L[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					d_xyindex[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk], NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					d_zindex[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						d_trIndex[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						d_axIndex[kk] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
			}
		}

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Buffer creation failed\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer creation succeeded\n");
		}


		// assign values to the buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			status = CLCommandQueue[0].enqueueWriteBuffer(d_weights, CL_FALSE, 0, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, w_vec.weights);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			memSize += (sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1)) / 1048576ULL;
		}
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			cl::detail::size_t_array region = { { 0, 0, 0 } };
			region[0] = inputScalars.Nx[0];
			region[1] = inputScalars.Ny[0];
			region[2] = inputScalars.Nz[0];
			if (inputScalars.useImages)
				status = CLCommandQueue[0].enqueueWriteImage(d_urefIm, CL_FALSE, origin, region, 0, 0, w_vec.NLM_ref);
			else
				status = CLCommandQueue[0].enqueueWriteBuffer(d_uref, CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[0], w_vec.NLM_ref);
			memSize += (sizeof(float) * inputScalars.im_dim[0]) / 1048576ULL;
		}
		if (inputScalars.maskFP || inputScalars.maskBP || inputScalars.useExtendedFOV) {
			//mexPrintBase("inputScalars.useBuffers = %u\n", inputScalars.useBuffers);
			//mexEval();
			if (inputScalars.useBuffers) {
				if (inputScalars.maskFP) {
					if (inputScalars.maskFPZ > 1) {
						for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
							status = CLCommandQueue[0].enqueueWriteBuffer(d_maskFPB[kk], CL_FALSE, 0, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk], &w_vec.maskFP[pituus[kk] * vecSize]);
					}
					else
						status = CLCommandQueue[0].enqueueWriteBuffer(d_maskFPB[0], CL_FALSE, 0, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD, w_vec.maskFP);
				}
				if (inputScalars.maskBP)
					status = CLCommandQueue[0].enqueueWriteBuffer(d_maskBPB, CL_FALSE, 0, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ, w_vec.maskBP);
			}
			else {
				cl::detail::size_t_array region = { { 1, 1, 1 } };
				if (inputScalars.maskFP) {
					region[0] = inputScalars.nRowsD;
					region[1] = inputScalars.nColsD;
					region[2] = inputScalars.maskFPZ;
					if (inputScalars.maskFPZ > 1) {
						for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
							region[2] = length[kk];
							status = CLCommandQueue[0].enqueueWriteImage(d_maskFP3[kk], CL_FALSE, origin, region, 0, 0, &w_vec.maskFP[pituus[kk] * vecSize]);
						}
					}
					else
						status = CLCommandQueue[0].enqueueWriteImage(d_maskFP, CL_FALSE, origin, region, 0, 0, w_vec.maskFP);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(bool) * inputScalars.nRowsD * inputScalars.nColsD) / 1048576ULL;
				}
				if (inputScalars.maskBP) {
					region[0] = inputScalars.Nx[0];
					region[1] = inputScalars.Ny[0];
					region[2] = inputScalars.maskBPZ;
					if (DEBUG) {
						mexPrintBase("region[0] = %u\n", region[0]);
						mexPrintBase("region[1] = %u\n", region[1]);
						mexPrintBase("region[2] = %u\n", region[2]);
						mexEval();
					}
					if (inputScalars.maskBPZ > 1)
						status = CLCommandQueue[0].enqueueWriteImage(d_maskBP3, CL_FALSE, origin, region, 0, 0, w_vec.maskBP);
					else
						status = CLCommandQueue[0].enqueueWriteImage(d_maskBP, CL_FALSE, origin, region, 0, 0, w_vec.maskBP);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(bool) * inputScalars.Nx[0] * inputScalars.Ny[0]) / 1048576ULL;
				}
			}
			if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
				cl::detail::size_t_array region = { { 1, 1, 1 } };
				region[0] = inputScalars.Nx[0];
				region[1] = inputScalars.Ny[0];
				region[2] = inputScalars.maskBPZ;
				if (inputScalars.maskBPZ > 1)
					status = CLCommandQueue[0].enqueueWriteImage(d_maskPrior3, CL_FALSE, origin, region, 0, 0, w_vec.maskPrior);
				else
					status = CLCommandQueue[0].enqueueWriteImage(d_maskPrior, CL_FALSE, origin, region, 0, 0, w_vec.maskPrior);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(bool) * inputScalars.Nx[0] * inputScalars.Ny[0]) / 1048576ULL;
			}
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_V, CL_FALSE, 0, sizeof(float) * inputScalars.size_V, inputScalars.V);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * inputScalars.size_V) / 1048576ULL;
			}
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_x[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_of_x, x);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * inputScalars.size_of_x) / 1048576ULL;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.nLayers > 1) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_xFull[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_of_x, x);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = CLCommandQueue[0].enqueueWriteBuffer(d_zFull[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_z, z_det);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_xFull[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_of_x, x);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = CLCommandQueue[0].enqueueWriteBuffer(d_zFull[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_z, z_det);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
			}
			if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.useBuffers)
					status = CLCommandQueue[0].enqueueWriteBuffer(d_attenB, CL_FALSE, 0, sizeof(float) * inputScalars.im_dim[0], atten);
				else {
					cl::detail::size_t_array region = { { 0, 0, 0 } };
					region[0] = inputScalars.Nx[0];
					region[1] = inputScalars.Ny[0];
					region[2] = inputScalars.Nz[0];
					status = CLCommandQueue[0].enqueueWriteImage(d_attenIm, CL_FALSE, origin, region, 0, 0, atten);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				memSize += (sizeof(float) * inputScalars.im_dim[0]) / 1048576ULL;
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_angle, CL_FALSE, 0, sizeof(float) * inputScalars.nProjections, w_vec.angles);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * inputScalars.nProjections) / 1048576ULL;
			}
			if (inputScalars.TOF) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_TOFCenter, CL_FALSE, 0, sizeof(float) * inputScalars.nBins, inputScalars.TOFCenter);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * inputScalars.nBins) / 1048576ULL;
			}
			status = CLCommandQueue[0].finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.eFOV) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_eFOVIndices, CL_FALSE, 0, sizeof(uint8_t) * inputScalars.Nz[0], w_vec.eFOVIndices);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(uint8_t) * inputScalars.Nz[0]) / 1048576ULL;
			}
			if (inputScalars.SPECT) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_rayShiftsDetector, CL_FALSE, 0, sizeof(float) * 2 * inputScalars.n_rays, w_vec.rayShiftsDetector);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * 2 * inputScalars.n_rays) / 1048576ULL;

				status = CLCommandQueue[0].enqueueWriteBuffer(d_rayShiftsSource, CL_FALSE, 0, sizeof(float) * 2 * inputScalars.n_rays, w_vec.rayShiftsSource);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * 2 * inputScalars.n_rays) / 1048576ULL;
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
					int64_t kerroin = 2;
					if (inputScalars.pitch)
						kerroin = 6;
					status = CLCommandQueue[0].enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * kerroin, &z_det[pituus[kk] * kerroin]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * kerroin) / 1048576ULL;
				}
				else {
					if (inputScalars.PET && inputScalars.listmode == 0) {
						int64_t kerroin = 2;
						if (inputScalars.nLayers > 1)
							int64_t kerroin = 3;
						status = CLCommandQueue[0].enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * kerroin, &z_det[pituus[kk] * kerroin]);
						memSize += (sizeof(float) * length[kk] * kerroin) / 1048576ULL;
					}
					else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
						status = CLCommandQueue[0].enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * inputScalars.size_z, z_det);
						memSize += (sizeof(float) * inputScalars.size_z) / 1048576ULL;
					}
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_x[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 6, &x[pituus[kk] * 6]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * 6) / 1048576ULL;
				}
				else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
					if (kk < inputScalars.TOFsubsets || inputScalars.loadTOF)
						status = CLCommandQueue[0].enqueueWriteBuffer(d_x[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 6, &w_vec.listCoord[pituus[kk] * 6]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * 6) / 1048576ULL;
				}
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_T[kk], CL_FALSE, 0, sizeof(float) * length[kk], &inputScalars.T[pituus[kk]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.raw) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(uint16_t) * length[kk] * 2) / 1048576ULL;
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &z_index[pituus[kk]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = CLCommandQueue[0].enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(uint32_t) * length[kk] + sizeof(uint16_t) * length[kk]) / 1048576ULL;
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased) {
					if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF)) {
						status = CLCommandQueue[0].enqueueWriteBuffer(d_trIndex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &w_vec.trIndex[pituus[kk] * 2]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						memSize += (sizeof(uint16_t) * length[kk] * 2) / 1048576ULL;
						status = CLCommandQueue[0].enqueueWriteBuffer(d_axIndex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &w_vec.axIndex[pituus[kk] * 2]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						memSize += (sizeof(uint16_t) * length[kk] * 2) / 1048576ULL;
					}
				}
				status = CLCommandQueue[0].finish();
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_norm[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &norm[pituus[kk] * vecSize]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
				}
				if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_scat[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &extraCorr[pituus[kk] * vecSize]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
				}
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					status = CLCommandQueue[0].enqueueWriteBuffer(d_atten[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &atten[pituus[kk] * vecSize]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
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
					mexPrintBase("memSize = %u\n", memSize);
					mexEval();
				}
				status = CLCommandQueue[0].finish();
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}

		if (status != CL_SUCCESS) {
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
		cl_int status = CL_SUCCESS;
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
		if (inputScalars.normalization_correction)
			d_norm.resize(inputScalars.subsetsUsed);
		if (inputScalars.scatter)
			d_scat.resize(inputScalars.subsetsUsed);
		if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
			d_atten.resize(inputScalars.subsetsUsed);
		if (inputScalars.projector_type != 6) {
			d_x.resize(inputScalars.subsetsUsed);
			d_z.resize(inputScalars.subsetsUsed);
		}
		if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
			d_T.resize(inputScalars.subsetsUsed);

		status = createAndWriteBuffers(length, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, inputScalars, w_vec, MethodList);
		if (status != CL_SUCCESS) {
			return -1;
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
		cl_int status = CL_SUCCESS;

		if (inputScalars.FPType == 4 || inputScalars.FPType == 5) {
			kernelFP.setArg(kernelIndFP++, inputScalars.nRowsD);
			kernelFP.setArg(kernelIndFP++, inputScalars.nColsD);
			status = kernelFP.setArg(kernelIndFP++, dPitch);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}

		if (inputScalars.BPType == 4 || inputScalars.BPType == 5) {
			kernelBP.setArg(kernelIndBP++, inputScalars.nRowsD);
			kernelBP.setArg(kernelIndBP++, inputScalars.nColsD);
			kernelBP.setArg(kernelIndBP++, dPitch);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				kernelSensList.setArg(kernelIndSens++, inputScalars.nRowsD);
				kernelSensList.setArg(kernelIndSens++, inputScalars.nColsD);
				kernelSensList.setArg(kernelIndSens++, dPitch);
			}
		}
		if (inputScalars.FPType == 4) {
			kernelFP.setArg(kernelIndFP++, inputScalars.dL);
			status = kernelFP.setArg(kernelIndFP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.BPType == 4 && !inputScalars.CT) {
			kernelBP.setArg(kernelIndBP++, inputScalars.dL);
			status = kernelBP.setArg(kernelIndBP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				kernelSensList.setArg(kernelIndSens++, inputScalars.dL);
				status = kernelSensList.setArg(kernelIndSens++, inputScalars.global_factor);
			}
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			// Set the kernelFP parameters that do not change
			status = kernelFP.setArg(kernelIndFP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.epps));
			getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.nRowsD));
			getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.det_per_ring));
			getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.sigma_x));
			if (inputScalars.SPECT) {
				getErrorString(kernelFP.setArg(kernelIndFP++, d_rayShiftsDetector));
				getErrorString(kernelFP.setArg(kernelIndFP++, d_rayShiftsSource));
			}

			getErrorString(kernelFP.setArg(kernelIndFP++, dPitch));
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (inputScalars.FPType == 2) {
					getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.tube_width));
				}
				else {
					getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.cylRadiusProj3));
				}
				getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.bmin));
				getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.bmax));
				getErrorString(kernelFP.setArg(kernelIndFP++, inputScalars.Vmax));
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			status = kernelBP.setArg(kernelIndBP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.epps));
			getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.nRowsD));
			getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.det_per_ring));
			getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.sigma_x));
			if (inputScalars.SPECT) {
				getErrorString(kernelBP.setArg(kernelIndBP++, d_rayShiftsDetector));
				getErrorString(kernelBP.setArg(kernelIndBP++, d_rayShiftsSource));
			}
			getErrorString(kernelBP.setArg(kernelIndBP++, dPitch));
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (inputScalars.BPType == 2) {
					getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.tube_width));
				}
				else {
					getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.cylRadiusProj3));
				}
				getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.bmin));
				getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.bmax));
				getErrorString(kernelBP.setArg(kernelIndBP++, inputScalars.Vmax));
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				status = kernelSensList.setArg(kernelIndSens++, inputScalars.global_factor);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				kernelSensList.setArg(kernelIndSens++, inputScalars.epps);
				kernelSensList.setArg(kernelIndSens++, inputScalars.nRowsD);
				kernelSensList.setArg(kernelIndSens++, inputScalars.det_per_ring);
				kernelSensList.setArg(kernelIndSens++, inputScalars.sigma_x);
				kernelSensList.setArg(kernelIndSens++, dPitch);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					if (inputScalars.BPType == 2)
						kernelSensList.setArg(kernelIndSens++, inputScalars.tube_width);
					else
						kernelSensList.setArg(kernelIndSens++, inputScalars.cylRadiusProj3);
					kernelSensList.setArg(kernelIndSens++, inputScalars.bmin);
					kernelSensList.setArg(kernelIndSens++, inputScalars.bmax);
					kernelSensList.setArg(kernelIndSens++, inputScalars.Vmax);
				}
			}
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if (DEBUG) {
				mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
				mexEval();
			}
			if (inputScalars.TOF) {
				status = kernelFP.setArg(kernelIndFP++, d_TOFCenter);
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				status = kernelFP.setArg(kernelIndFP++, d_V);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFP++, inputScalars.nColsD);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if (inputScalars.TOF)
				kernelBP.setArg(kernelIndBP++, d_TOFCenter);
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				kernelBP.setArg(kernelIndBP++, d_V);
			}
			status = kernelBP.setArg(kernelIndBP++, inputScalars.nColsD);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.TOF)
					kernelSensList.setArg(kernelIndSens++, d_TOFCenter);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					kernelSensList.setArg(kernelIndSens++, d_V);
				}
				status = kernelSensList.setArg(kernelIndSens++, inputScalars.nColsD);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		if ((inputScalars.BPType == 4 || inputScalars.FPType == 4) && !inputScalars.CT && inputScalars.TOF) {
			if (inputScalars.FPType == 4) {
				kernelFP.setArg(kernelIndFP++, d_TOFCenter);
				kernelFP.setArg(kernelIndFP++, inputScalars.sigma_x);
			}
			if (inputScalars.BPType == 4) {
				kernelBP.setArg(kernelIndBP++, d_TOFCenter);
				kernelBP.setArg(kernelIndBP++, inputScalars.sigma_x);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					kernelSensList.setArg(kernelIndSens++, d_TOFCenter);
					kernelSensList.setArg(kernelIndSens++, inputScalars.sigma_x);
				}
			}
		}
		if (DEBUG) {
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("kernelIndBP = %u\n", kernelIndBP);
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

		cl_int status = CL_SUCCESS;
		for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
			if (inputScalars.scatter == 1u) {
				status = CLCommandQueue[0].enqueueWriteBuffer(d_scat[kk], CL_TRUE, 0, sizeof(float) * length[kk], &extraCorr[pituus[kk] + inputScalars.koko * tt]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				memSize += (sizeof(float) * length[kk]) / 1048576ULL;
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
		cl_int status = CL_SUCCESS;
		if (inputScalars.attenuation_correction && !inputScalars.CT && inputScalars.CTAttenuation) {
			if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 || inputScalars.FPType == 4)) {
				if (inputScalars.useBuffers)
					status = kernelFP.setArg(kernelIndFP++, d_attenB);
				else
					status = kernelFP.setArg(kernelIndFP++, d_attenIm);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.BPType == 4) {
				if (inputScalars.useBuffers)
					status = kernelBP.setArg(kernelIndBP++, d_attenB);
				else
					status = kernelBP.setArg(kernelIndBP++, d_attenIm);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					if (inputScalars.useBuffers)
						kernelSensList.setArg(kernelIndBP++, d_attenB);
					else
						kernelSensList.setArg(kernelIndBP++, d_attenIm);
				}
			}
		}
		return status;
	}
	template <typename T>
	inline int loadCoord(scalarStruct& inputScalars, const int64_t length, const T* listCoord, const T* listCoordAx = nullptr) {
		cl_int status = CL_SUCCESS;
		if (inputScalars.indexBased) {
			d_trIndex[0] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			d_axIndex[0] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = CLCommandQueue[0].enqueueWriteBuffer(d_trIndex[0], CL_FALSE, 0, sizeof(uint16_t) * length * 2, listCoord);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = CLCommandQueue[0].enqueueWriteBuffer(d_axIndex[0], CL_FALSE, 0, sizeof(uint16_t) * length * 2, listCoordAx);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		else {
		d_x[0] = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * length * 6, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		status = CLCommandQueue[0].enqueueWriteBuffer(d_x[0], CL_FALSE, 0, sizeof(float) * length * 6, listCoord);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
			}
		}
		return 0;
	}

	int computeConvolutionF(const scalarStruct& inputScalars, const int ii = 0) {
		int status = CL_SUCCESS;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		cl_uint kernelInd = 0U;
		kernelPSFf.setArg(kernelInd++, d_imFinal[ii]);
		kernelPSFf.setArg(kernelInd++, d_imTemp[ii]);
		kernelPSFf.setArg(kernelInd++, d_g);
		kernelPSFf.setArg(kernelInd++, inputScalars.g_dim_x);
		kernelPSFf.setArg(kernelInd++, inputScalars.g_dim_y);
		kernelPSFf.setArg(kernelInd++, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSFf, cl::NDRange(), globalC, localPrior, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");

		return status;
	}

	template <typename T>
	int computeConvolution(const scalarStruct& inputScalars, cl::Buffer& input, const int ii = 0, const T& cType = 0) {
		int status = CL_SUCCESS;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		cl::Buffer d_BPApu = cl::Buffer(CLContext, CL_MEM_READ_WRITE, sizeof(T) * inputScalars.im_dim[ii], NULL, &status);
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		cl_uint kernelInd = 0U;
		kernelPSF.setArg(kernelInd++, input);
		kernelPSF.setArg(kernelInd++, d_BPApu);
		kernelPSF.setArg(kernelInd++, d_g);
		kernelPSF.setArg(kernelInd++, inputScalars.g_dim_x);
		kernelPSF.setArg(kernelInd++, inputScalars.g_dim_y);
		kernelPSF.setArg(kernelInd++, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSF, cl::NDRange(), globalC, localPrior, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");
		status = CLCommandQueue[0].enqueueCopyBuffer(d_BPApu, input, 0, 0, sizeof(T) * inputScalars.im_dim[ii]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}

		return status;
	}

	int computeForward(const scalarStruct& inputScalars, const std::vector<int64_t>& length, const uint32_t osa_iter) {
		int status = CL_SUCCESS;
		cl::NDRange localF = { 64, 1, 1 };
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { inputScalars.nRowsD * inputScalars.nColsD * static_cast<size_t>(length[osa_iter]) * inputScalars.nBins, 1, 1 };
		else
			global = { static_cast<cl::size_type>(length[osa_iter]) * inputScalars.nBins, 1, 1 };
		
		size_t erotusF = global[0] % localF[0];
		if (erotusF > 0)
			erotusF = localF[0] - erotusF;
		global = { static_cast<cl::size_type>(global[0] + erotusF), 1, 1 };
		//else {
			//erotus[0] = length[osa_iter] % local_size[0];

			//if (erotus[0] > 0)
			//	erotus[0] = (local_size[0] - erotus[0]);
			//global = { static_cast<cl::size_type>(length[osa_iter] + erotus[0]), 1, 1 };
			//global = { static_cast<cl::size_type>(length[osa_iter]), 1, 1 };
		//}
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			//mexPrintBase("local[0] = %u\n", local[0]);
			//mexPrintBase("local[1] = %u\n", local[1]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("global.dimensions() = %u\n", global.dimensions());
			//mexPrintBase("local.dimensions() = %u\n", local.dimensions());
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexEval();
		}
		cl_uint kernelInd = 0U;
		kernelForward.setArg(kernelInd++, d_output);
		kernelForward.setArg(kernelInd++, d_meas[osa_iter]);
		if (inputScalars.CT)
			kernelForward.setArg(kernelInd++, d_outputCT);
		if (inputScalars.randoms_correction)
			kernelForward.setArg(kernelInd++, d_rand[osa_iter]);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelForward, cl::NDRange(), global, localF, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward step kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Forward step computed");

		return status;
	}

	int computeEstimate(const scalarStruct& inputScalars, const int ii = 0, const int uu = 0) {
		int status = CL_SUCCESS;
		global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		cl_uint kernelInd = 0U;
		kernelEstimate.setArg(kernelInd++, d_Summ[uu]);
		kernelEstimate.setArg(kernelInd++, vec_opencl.d_rhs_os[ii]);
		kernelEstimate.setArg(kernelInd++, d_imFinal[ii]);
		kernelEstimate.setArg(kernelInd++, inputScalars.epps);
		kernelEstimate.setArg(kernelInd++, d_N[ii]);
		kernelEstimate.setArg(kernelInd++, no_norm);
		//if (inputScalars.use_psf) {
		//	kernelEstimate.setArg(kernelInd++, d_g);
		//}
		if (inputScalars.CT)
			kernelEstimate.setArg(kernelInd++, inputScalars.flat);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelEstimate, cl::NDRange(), global, localPrior, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Estimate step kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Estimate step computed");

		return status;
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
	inline int forwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const std::vector<int64_t>&length, const uint64_t m_size, const int32_t ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
		kernelIndFPSubIter = kernelIndFP;
		cl_int status = CL_SUCCESS;
		if (inputScalars.FPType == 5)
			global = { inputScalars.nRowsD + erotus[0], (inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1], static_cast<size_t>(length[osa_iter]) };
		else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[osa_iter]) };
		else {
			erotus[0] = length[osa_iter] % local_size[0];

			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
			global = { static_cast<cl::size_type>(length[osa_iter] + erotus[0]), 1, 1 };
		}

		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("local[0] = %u\n", local[0]);
			mexPrintBase("local[1] = %u\n", local[1]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("d_N[ii].s0 = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N[ii].s1 = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N[ii].s2 = %u\n", d_N[ii].s[2]);
			mexPrintBase("d[ii].s0 = %f\n", d[ii].s[0]);
			mexPrintBase("d[ii].s1 = %f\n", d[ii].s[1]);
			mexPrintBase("d[ii].s2 = %f\n", d[ii].s[2]);
			mexPrintBase("b[ii].s0 = %f\n", b[ii].s[0]);
			mexPrintBase("b[ii].s1 = %f\n", b[ii].s[1]);
			mexPrintBase("b[ii].s2 = %f\n", b[ii].s[2]);
			mexPrintBase("bmax[ii].s0 = %f\n", bmax[ii].s[0]);
			mexPrintBase("bmax[ii].s1 = %f\n", bmax[ii].s[1]);
			mexPrintBase("bmax[ii].s2 = %f\n", bmax[ii].s[2]);
			mexPrintBase("global.dimensions() = %u\n", global.dimensions());
			mexPrintBase("local.dimensions() = %u\n", local.dimensions());
			mexPrintBase("kernelIndFPSubIter = %u\n", kernelIndFPSubIter);
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("FPType = %u\n", inputScalars.FPType);
			if (inputScalars.FPType == 4) {
				mexPrintBase("dL = %f\n", inputScalars.dL);
				mexPrintBase("d_Scale4[ii].s[0] = %f\n", inputScalars.d_Scale4[ii].s[0]);
				mexPrintBase("d_Scale4[ii].s[1] = %f\n", inputScalars.d_Scale4[ii].s[1]);
				mexPrintBase("d_Scale4[ii].s[2] = %f\n", inputScalars.d_Scale4[ii].s[2]);
			}
			mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
			mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexPrintBase("maskBP = %u\n", inputScalars.maskBP);
			mexPrintBase("maskFP = %u\n", inputScalars.maskFP);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("ii = %u\n", ii);
			mexPrintBase("NVOXELS = %u\n", NVOXELS);
			mexPrintBase("NVOXELS5 = %u\n", NVOXELS5);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexPrintBase("memSize = %u\n", memSize);
			mexPrintBase("subsetType = %u\n", inputScalars.subsetType);
			mexEval();
		}

		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.FPType == 4) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
				status = kernelFP.setArg(kernelIndFPSubIter++, d_atten[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.FPType == 4) {
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_N[ii]));
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, b[ii]));
			if (inputScalars.FPType == 5) {
				status = kernelFP.setArg(kernelIndFPSubIter++, inputScalars.dSize[ii]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelFP.setArg(kernelIndFPSubIter++, d[ii]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelFP.setArg(kernelIndFPSubIter++, inputScalars.d_Scale[ii]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else {
				getErrorString(kernelFP.setArg(kernelIndFPSubIter++, bmax[ii]));
				getErrorString(kernelFP.setArg(kernelIndFPSubIter++, inputScalars.d_Scale4[ii]));
			}
		}
		if (inputScalars.FPType == 4) {
			status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, d_output);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.listmode == 0 && !inputScalars.CT) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[0]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.CT || inputScalars.PET || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[inputScalars.osa_iter0]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1)
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP3[osa_iter]);
					else
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_xyindex[osa_iter]));
				getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_zindex[osa_iter]));
			}
			if (inputScalars.raw) {
				kernelFP.setArg(kernelIndFPSubIter++, d_L[osa_iter]);
				kernelFP.setArg(kernelIndFPSubIter++, inputScalars.det_per_ring);
			}
			if (inputScalars.normalization_correction)
				status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[osa_iter]);
			if (inputScalars.scatter)
				status = kernelFP.setArg(kernelIndFPSubIter++, d_scat[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, no_norm);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_ulong>(m_size));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, osa_iter);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, ii));
		}
		else if (inputScalars.FPType == 5) {
			if (!inputScalars.loadTOF && inputScalars.listmode > 0)
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[0]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os_int);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, d_output);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.meanFP) {
				status = kernelBP.setArg(kernelIndBPSubIter++, d_meanFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1)
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP3[osa_iter]);
					else
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]));
		}
		else if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
				status = kernelFP.setArg(kernelIndFPSubIter++, d_atten[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1)
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP3[osa_iter]);
					else
						status = kernelFP.setArg(kernelIndFPSubIter++, d_maskFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0) {
				status = kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[0]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[inputScalars.osa_iter0]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.normalization_correction)
				kernelFP.setArg(kernelIndFPSubIter++, d_norm[osa_iter]);
			if (inputScalars.scatter)
				kernelFP.setArg(kernelIndFPSubIter++, d_scat[osa_iter]);
			status = kernelFP.setArg(kernelIndFPSubIter++, d_Summ[uu]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_N[ii]));
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d[ii]));
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, b[ii]));
			status = kernelFP.setArg(kernelIndFPSubIter++, bmax[ii]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				status = kernelFP.setArg(kernelIndFPSubIter++, d_xyindex[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelFP.setArg(kernelIndFPSubIter++, d_zindex[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				if (!inputScalars.loadTOF) {
					getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_trIndex[0]));
					getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_axIndex[0]));
				}
				else {
					getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_trIndex[osa_iter]));
					getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_axIndex[osa_iter]));
				}
			}
			if (inputScalars.raw)
				kernelFP.setArg(kernelIndFPSubIter++, d_L[osa_iter]);
			if (inputScalars.useBuffers)
				status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_im);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, d_output));
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, no_norm));
			status = kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_ulong>(m_size));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelFP.setArg(kernelIndFPSubIter++, osa_iter);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelFP.setArg(kernelIndFPSubIter++, ii));
		}
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelFP, cl::NDRange(), global, local, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward projection kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.verbose >= 3)
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
	inline int backwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const std::vector<int64_t>&length, const uint64_t m_size, const bool compSens = false, const int32_t ii = 0, const int uu = 0, 
		int ee = -1) {
		if (inputScalars.verbose >= 3)
			mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
		if (ee < 0)
			ee = uu;
		cl_int status = CL_SUCCESS;
		kernelIndBPSubIter = kernelIndBP;
		if (inputScalars.listmode > 0 && compSens) {
			kernelApu = kernelBP;
			kernelBP = kernelSensList;
			kernelIndBPSubIter = kernelIndSens;
		}

		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[osa_iter]) };
			else if (inputScalars.listmode > 0 && compSens)
				if (inputScalars.nLayers > 1)
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], 
					static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.nLayers)};
				else
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) };
			else {
				erotus[0] = length[osa_iter] % local_size[0];

				if (erotus[0] > 0)
					erotus[0] = (local_size[0] - erotus[0]);
				global = { static_cast<cl::size_type>(length[osa_iter] + erotus[0]), 1, 1 };
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
				mexPrintBase("global.dimensions() = %u\n", global.dimensions());
				mexPrintBase("local.dimensions() = %u\n", local.dimensions());
				mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
				mexPrintBase("m_size = %u\n", m_size);
				mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
				mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
				mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
				mexPrintBase("listmode = %u\n", inputScalars.listmode);
				mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("no_norm = %u\n", no_norm);
				mexPrintBase("ii = %u\n", ii);
				mexPrintBase("memSize = %u\n", memSize);
				mexPrintBase("compSens = %u\n", compSens);
				mexEval();
			}

			// Set kernelBP arguments
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
				status = kernelBP.setArg(kernelIndBPSubIter++, d_atten[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						status = kernelBP.setArg(kernelIndBPSubIter++, d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1)
							status = kernelBP.setArg(kernelIndBPSubIter++, d_maskFP3[osa_iter]);
						else
							status = kernelBP.setArg(kernelIndBPSubIter++, d_maskFP);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				if (inputScalars.maskBP) {
					if (inputScalars.useBuffers)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBPB);
					else
						if (inputScalars.maskBPZ > 1)
							status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP3);
						else
							status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
						if (inputScalars.useBuffers)
							status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBPB);
						else
							if (inputScalars.maskBPZ > 1)
								status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBP3);
							else
								status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBP);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0)
				status = kernelBP.setArg(kernelIndBPSubIter++, length[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (compSens) {
				status = kernelBP.setArg(kernelIndBPSubIter++, d_xFull[0]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_zFull[0]));
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, inputScalars.rings));
			}
			else {
				if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
					status = kernelBP.setArg(kernelIndBPSubIter++, d_x[0]);
				else
					status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased)))
					status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
				else if (inputScalars.indexBased && inputScalars.listmode > 0)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_z[0]);
				else
					status = kernelBP.setArg(kernelIndBPSubIter++, d_z[inputScalars.osa_iter0]);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (compSens) {
				if (inputScalars.normalization_correction)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_normFull[0]);
				if (inputScalars.scatter)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_scatFull[0]);
			}
			else {
				if (inputScalars.normalization_correction)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[osa_iter]);
				if (inputScalars.scatter)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_scat[osa_iter]);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ[ee]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_N[ii]));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d[ii]));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, b[ii]));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, bmax[ii]));
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
					getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_xyindex[osa_iter]));
					getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_zindex[osa_iter]));
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
					if (!inputScalars.loadTOF) {
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_trIndex[0]));
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_axIndex[0]));
					}
					else {
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_trIndex[osa_iter]));
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_axIndex[osa_iter]));
					}
				}
				if (inputScalars.raw)
					kernelBP.setArg(kernelIndBPSubIter++, d_L[osa_iter]);
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_output));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os[uu]));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, no_norm));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_ulong>(m_size)));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, osa_iter));
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, ii));
		}
		else {
			if (inputScalars.CT) {

				if (!inputScalars.useBuffers) {
					cl::size_type imX = inputScalars.nRowsD;
					cl::size_type imY = inputScalars.nColsD;
					cl::size_type imZ = length[osa_iter];
					if (inputScalars.BPType == 5) {
						imX++;
						imY++;
					}
					cl::detail::size_t_array region = { imX, imY, imZ };
					d_inputImage = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrint("Image creation failed\n");
						return -1;
					}

					status = CLCommandQueue[0].enqueueCopyBufferToImage(d_output, d_inputImage, 0, origin, region);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrint("Image copy failed\n");
						return -1;
					}
					status = CLCommandQueue[0].finish();
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrint("Queue finish failed after image copy\n");
						return -1;
					}
				}
				if (inputScalars.BPType == 4)
					if (!inputScalars.largeDim)
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS };
					else
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };
				else if (inputScalars.BPType == 5) {
					if (inputScalars.pitch)
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };
					else
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5 };
				}
				else
					global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("local[0] = %u\n", local[0]);
					mexPrintBase("local[1] = %u\n", local[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
					mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("local.dimensions() = %u\n", local.dimensions());
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
					mexPrintBase("memSize = %u\n", memSize);
					if (inputScalars.BPType == 4) {
						mexPrintBase("dL = %f\n", inputScalars.dL);
						mexPrintBase("kerroin4[ii] = %f\n", w_vec.kerroin4[ii]);
					}
					mexEval();
				}
				if (inputScalars.offset)
					kernelBP.setArg(kernelIndBPSubIter++, d_T[osa_iter]);
				if (inputScalars.BPType == 5 || inputScalars.BPType == 4) {
					getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_N[ii]));
					status = kernelBP.setArg(kernelIndBPSubIter++, b[ii]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, d[ii]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (inputScalars.BPType == 5) {
						status = kernelBP.setArg(kernelIndBPSubIter++, inputScalars.d_Scale[ii]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
						status = kernelBP.setArg(kernelIndBPSubIter++, inputScalars.dSizeBP);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
					else {
						status = kernelBP.setArg(kernelIndBPSubIter++, w_vec.kerroin4[ii]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
				if (inputScalars.BPType == 4) {
					if (inputScalars.useBuffers)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_output);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_inputImage);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (inputScalars.CT && inputScalars.DSC > 0.f) {
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_angle));
						getErrorString(kernelBP.setArg(kernelIndBPSubIter++, inputScalars.DSC));
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os[uu]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (compSens)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_xFull[0]);
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0)
							status = kernelBP.setArg(kernelIndBPSubIter++, d_x[0]);
						else
							status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (compSens)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_zFull[0]);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ[ee]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
				else {
					if (compSens)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_xFull[0]);
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0)
							status = kernelBP.setArg(kernelIndBPSubIter++, d_x[0]);
						else
							status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (compSens)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_zFull[0]);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, d_inputImage);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os[uu]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ[ee]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if (inputScalars.meanBP) {
						status = kernelBP.setArg(kernelIndBPSubIter++, d_meanBP);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return -1;
						}
					}
				}
			}
			else {
				if ((inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
					global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[osa_iter]) };
				else if (inputScalars.listmode > 0 && compSens)
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) };
				else {
					erotus[0] = length[osa_iter] % local_size[0];

					if (erotus[0] > 0)
						erotus[0] = (local_size[0] - erotus[0]);
					global = { static_cast<cl::size_type>(length[osa_iter] + erotus[0]), 1, 1 };
				}

				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("local[0] = %u\n", local[0]);
					mexPrintBase("local[1] = %u\n", local[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					if (compSens) {
						mexPrintBase("erotusSens[0] = %u\n", erotusSens[0]);
						mexPrintBase("erotusSens[1] = %u\n", erotusSens[1]);
					}
					else {
						mexPrintBase("erotus[0] = %u\n", erotus[0]);
						mexPrintBase("erotus[1] = %u\n", erotus[1]);
					}
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("local.dimensions() = %u\n", local.dimensions());
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("osa_iter = %u\n", osa_iter);
					mexPrintBase("memSize = %u\n", memSize);
					mexEval();
				}
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, d_N[ii]));
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, b[ii]));
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, bmax[ii]));
				getErrorString(kernelBP.setArg(kernelIndBPSubIter++, inputScalars.d_Scale4[ii]));
				status = kernelBP.setArg(kernelIndBPSubIter++, d_output);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os[uu]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (compSens) {
					status = kernelBP.setArg(kernelIndBPSubIter++, d_xFull[0]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					status = kernelBP.setArg(kernelIndBPSubIter++, d_zFull[0]);
					kernelBP.setArg(kernelIndBPSubIter++, inputScalars.rings);
					kernelBP.setArg(kernelIndBPSubIter++, inputScalars.det_per_ring);
				}
				else {
					if ((inputScalars.listmode == 0 && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0))
						status = kernelBP.setArg(kernelIndBPSubIter++, d_x[0]);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
					if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || inputScalars.listmode > 0))
						status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_z[inputScalars.osa_iter0]);
				}
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, length[osa_iter]);
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
					kernelBP.setArg(kernelIndBPSubIter++, d_xyindex[osa_iter]);
					kernelBP.setArg(kernelIndBPSubIter++, d_zindex[osa_iter]);
				}
				if (inputScalars.raw) {
					kernelBP.setArg(kernelIndBPSubIter++, d_L[osa_iter]);
				}
				if (inputScalars.normalization_correction)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[osa_iter]);
				if (inputScalars.scatter)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_scat[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ[ee]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			status = kernelBP.setArg(kernelIndBPSubIter++, no_norm);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.maskBP) {
				if (inputScalars.useBuffers)
					status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBPB);
				else
					if (inputScalars.maskBPZ > 1)
						status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP3);
					else
						status = kernelBP.setArg(kernelIndBPSubIter++, d_maskBP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					if (inputScalars.useBuffers)
						status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBPB);
					else
						if (inputScalars.maskBPZ > 1)
							status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBP3);
						else
							status = kernelSensList.setArg(kernelIndBPSubIter++, d_maskBP);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return -1;
					}
				}
			}
			if (inputScalars.CT)
				kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_long>(length[osa_iter]));
			else {
				status = kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_ulong>(m_size));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, osa_iter);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			getErrorString(kernelBP.setArg(kernelIndBPSubIter++, ii));
		}
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, local, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Backprojection kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
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
	}

	/// <summary>
	/// Get the total global memory of the selected device
	/// </summary>
	inline int64_t getGlobalMem() {
		cl_int status = CL_SUCCESS;
		int64_t mem;
		cl_ulong mem_loc;
		mem = CLDeviceID[0].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}

		mem_loc = CLDeviceID[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (DEBUG) {
			mexPrintBase("mem_loc = %u\n", mem_loc);
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
	inline int computeMRP(const scalarStruct& inputScalars, const uint64_t gSize[]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL median kernel computation");
		cl_uint kernelIndMed = 0U;
		uint64_t erotus[2] = { gSize[0] % localPrior[0], gSize[1] % localPrior[1] };
		cl::NDRange global_size(gSize[0] + (localPrior[0] - erotus[0]), gSize[1] + (localPrior[1] - erotus[1]), gSize[2]);
		CLCommandQueue[0].finish();
		if (DEBUG) {
			mexPrintBase("global_size[0] = %d\n", global_size[0]);
			mexPrintBase("global_size[1] = %d\n", global_size[1]);
			mexPrintBase("global_size[2] = %d\n", global_size[2]);
			mexPrintBase("erotus[0] = %d\n", erotus[0]);
			mexPrintBase("erotus[1] = %d\n", erotus[1]);
			mexPrintBase("gSize[0] = %d\n", gSize[0]);
			mexPrintBase("gSize[1] = %d\n", gSize[1]);
			mexEval();
		}
		kernelMed.setArg(kernelIndMed++, d_inputB);
		kernelMed.setArg(kernelIndMed++, d_W);
		kernelMed.setArg(kernelIndMed++, d_N[0]);
		kernelMed.setArg(kernelIndMed++, d_NOrig);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelMed.setArg(kernelIndMed++, d_maskPrior3);
			else
				kernelMed.setArg(kernelIndMed++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelMed.setArg(kernelIndMed++, d_eFOVIndices);
		cl_int status = CLCommandQueue[0].enqueueNDRangeKernel(kernelMed, cl::NullRange, global_size, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Median filter kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Median kernel launched successfully\n");
		}
		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after MRP kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL median kernel computed");
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
	inline int computeNLM(const scalarStruct& inputScalars, Weighting& w_vec, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL NLM gradient computation");
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
		const cl_int3 searchWindow = { static_cast<cl_int>(w_vec.Ndx) , static_cast<cl_int>(w_vec.Ndy) , static_cast<cl_int>(w_vec.Ndz) };
		const cl_int3 patchWindow = { static_cast<cl_int>(w_vec.Nlx) , static_cast<cl_int>(w_vec.Nly) , static_cast<cl_int>(w_vec.Nlz) };
		//const cl_int3 N = { static_cast<cl_int>(inputScalars.Nx), static_cast<cl_int>(inputScalars.Ny), static_cast<cl_int>(inputScalars.Nz) };
		cl_uint kernelIndNLM = 0ULL;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
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
		}
		kernelNLM.setArg(kernelIndNLM++, d_W);
		if (inputScalars.useImages) {
			kernelNLM.setArg(kernelIndNLM++, d_inputI);
		}
		else {
			kernelNLM.setArg(kernelIndNLM++, d_inputB);
		}
		kernelNLM.setArg(kernelIndNLM++, d_gaussianNLM);
		kernelNLM.setArg(kernelIndNLM++, d_N[0]);
		kernelNLM.setArg(kernelIndNLM++, d_NOrig);
		kernelNLM.setArg(kernelIndNLM++, w_vec.h2);
		kernelNLM.setArg(kernelIndNLM++, inputScalars.epps);
		kernelNLM.setArg(kernelIndNLM++, beta);
		if (w_vec.NLRD || w_vec.NLLange || w_vec.NLGGMRF)
			kernelNLM.setArg(kernelIndNLM++, w_vec.RDP_gamma);
		if (w_vec.NLGGMRF) {
			kernelNLM.setArg(kernelIndNLM++, w_vec.GGMRF_p);
			kernelNLM.setArg(kernelIndNLM++, w_vec.GGMRF_q);
			kernelNLM.setArg(kernelIndNLM++, w_vec.GGMRF_c);
		}
		if (w_vec.NLAdaptive)
			kernelNLM.setArg(kernelIndNLM++, w_vec.NLAdaptiveConstant);
		if (w_vec.NLM_anatomical)
			if (inputScalars.useImages)
				kernelNLM.setArg(kernelIndNLM++, d_urefIm);
			else
				kernelNLM.setArg(kernelIndNLM++, d_uref);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelNLM.setArg(kernelIndNLM++, d_maskPrior3);
			else
				kernelNLM.setArg(kernelIndNLM++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelNLM.setArg(kernelIndNLM++, d_eFOVIndices);
		 //Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelNLM, cl::NullRange, globalPrior, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the NLM kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after NLM kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL NLM gradient computed");
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
	inline int computeRDP(const scalarStruct& inputScalars, const float gamma, const float beta, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL RDP gradient computation");
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed before RDP kernel\n");
			return -1;
		}
		cl_uint kernelIndRDP = 0ULL;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
		if (DEBUG) {
			mexPrintBase("inputScalars.epps = %.9f\n", inputScalars.epps);
			mexPrintBase("gamma = %f\n", gamma);
			mexPrintBase("inputScalars.Nx = %d\n", inputScalars.Nx[0]);
			mexPrintBase("inputScalars.Ny = %d\n", inputScalars.Ny[0]);
			mexPrintBase("inputScalars.Nz * inputScalars.nRekos = %d\n", inputScalars.Nz[0] * inputScalars.nRekos);
			mexPrintBase("globalPrior[0] = %d\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %d\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %d\n", globalPrior[2]);
			mexPrintBase("RDPLargeNeighbor = %d\n", RDPLargeNeighbor);
			mexEval();
		}
		status = kernelRDP.setArg(kernelIndRDP++, d_W);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.useImages) {
			status = kernelRDP.setArg(kernelIndRDP++, d_inputI);
		}
		else {
			status = kernelRDP.setArg(kernelIndRDP++, d_inputB);
		}
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		kernelRDP.setArg(kernelIndRDP++, d_N[0]);
		kernelRDP.setArg(kernelIndRDP++, d_NOrig);
		kernelRDP.setArg(kernelIndRDP++, gamma);
		kernelRDP.setArg(kernelIndRDP++, inputScalars.epps);
		kernelRDP.setArg(kernelIndRDP++, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelRDP.setArg(kernelIndRDP++, d_maskPrior3);
			else
				kernelRDP.setArg(kernelIndRDP++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelRDP.setArg(kernelIndRDP++, d_eFOVIndices);
		if (RDPLargeNeighbor) {
			kernelRDP.setArg(kernelIndRDP++, d_weights);
			if (useRDPRef)
				if (inputScalars.useImages)
					kernelRDP.setArg(kernelIndRDP++, d_RDPrefI);
				else
					kernelRDP.setArg(kernelIndRDP++, d_RDPref);
		}
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelRDP, cl::NullRange, globalPrior, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the RDP kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after RDP kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL RDP gradient computed");
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
	inline int computeGGMRF(const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL GGMRF gradient computation");
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndGGMRF = 0ULL;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
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
		status = kernelGGMRF.setArg(kernelIndGGMRF++, d_W);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.useImages) {
			status = kernelGGMRF.setArg(kernelIndGGMRF++, d_inputI);
		}
		else {
			status = kernelGGMRF.setArg(kernelIndGGMRF++, d_inputB);
		}
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to input GGMRF buffer\n");
			return -1;
		}
		kernelGGMRF.setArg(kernelIndGGMRF++, d_weights);
		kernelGGMRF.setArg(kernelIndGGMRF++, d_N[0]);
		kernelGGMRF.setArg(kernelIndGGMRF++, p);
		kernelGGMRF.setArg(kernelIndGGMRF++, q);
		kernelGGMRF.setArg(kernelIndGGMRF++, c);
		kernelGGMRF.setArg(kernelIndGGMRF++, pqc);
		kernelGGMRF.setArg(kernelIndGGMRF++, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelGGMRF.setArg(kernelIndGGMRF++, d_maskPrior3);
			else
				kernelGGMRF.setArg(kernelIndGGMRF++, d_maskPrior);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelGGMRF, cl::NullRange, globalPrior, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the GGMRF kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after GGMRF kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL GGMRF gradient computed");
		return 0;
	}


	inline int ProxHelperQ(const float alpha, const uint64_t gQ) {
		cl_int status = CL_SUCCESS;
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndProxRDP = 0ULL;
		kernelProxq.setArg(kernelIndProxRDP++, d_qX);
		kernelProxq.setArg(kernelIndProxRDP++, alpha);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxq, cl::NullRange, globalQ, cl::NullRange);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal RDP helper kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after proximal RDP helper kernel\n");
			return -1;
		}
		return status;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TV prior
	/// </summary>
	/// <param name="q the input TV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <returns></returns>
	inline int ProxTVHelperQ(const float alpha, const uint64_t gQ) {
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		cl_int status = CL_SUCCESS;
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndCPTV = 0ULL;
		kernelProxTVq.setArg(kernelIndCPTV++, d_qX);
		kernelProxTVq.setArg(kernelIndCPTV++, d_qY);
		kernelProxTVq.setArg(kernelIndCPTV++, d_qZ);
		kernelProxTVq.setArg(kernelIndCPTV++, alpha);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVq, cl::NullRange, globalQ, cl::NullRange);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after kernel\n");
			return -1;
		}
		return status;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TGV prior
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
	inline int ProxTGVHelperQ(const scalarStruct& inputScalars, const float alpha, const uint64_t globalQ) {
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
		status = (CLCommandQueue[0]).finish();
		kernelProxTGVq.setArg(kernelIndCPTV++, d_rX);
		kernelProxTGVq.setArg(kernelIndCPTV++, d_rY);
		if (!inputScalars.TGV2D)
			kernelProxTGVq.setArg(kernelIndCPTV++, d_rZ);
		kernelProxTGVq.setArg(kernelIndCPTV++, d_rXY);
		if (!inputScalars.TGV2D) {
			kernelProxTGVq.setArg(kernelIndCPTV++, d_rXZ);
			kernelProxTGVq.setArg(kernelIndCPTV++, d_rYZ);
		}
		kernelProxTGVq.setArg(kernelIndCPTV++, alpha);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVq, cl::NullRange, globalQ, cl::NullRange);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after kernel\n");
			return -1;
		}
		return status;
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
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
		if (inputScalars.largeDim)
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
		if (DEBUG) {
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
			mexEval();
		}
		status = (CLCommandQueue[0]).finish();
		kernelProxTVDiv.setArg(kernelIndCPTV++, d_N[0]);
		kernelProxTVDiv.setArg(kernelIndCPTV++, d_NPrior);
		kernelProxTVDiv.setArg(kernelIndCPTV++, d_qX);
		kernelProxTVDiv.setArg(kernelIndCPTV++, d_qY);
		kernelProxTVDiv.setArg(kernelIndCPTV++, d_qZ);
		kernelProxTVDiv.setArg(kernelIndCPTV++, vec_opencl.d_rhs_os[0]);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelProxTVDiv.setArg(kernelIndCPTV++, d_maskPrior3);
			else
				kernelProxTVDiv.setArg(kernelIndCPTV++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelProxTVDiv.setArg(kernelIndCPTV++, d_eFOVIndices);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV divergence kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
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
	inline int ProxTVGrad(const scalarStruct& inputScalars, const float sigma2, const size_t vSize) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TV gradient");
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
		if (inputScalars.largeDim)
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
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
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
			mexPrintBase("vSize = %u\n", vSize);
			mexEval();
		}
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_N[0]);
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_NPrior);
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_inputB);
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_qX);
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_qY);
		kernelProxTVGrad.setArg(kernelIndCPTV++, d_qZ);
		kernelProxTVGrad.setArg(kernelIndCPTV++, sigma2);
		if (vSize > 0) {
			kernelProxTVGrad.setArg(kernelIndCPTV++, d_vX);
			kernelProxTVGrad.setArg(kernelIndCPTV++, d_vY);
			if (!inputScalars.TGV2D)
				kernelProxTVGrad.setArg(kernelIndCPTV++, d_vZ);
		}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelProxTVGrad.setArg(kernelIndCPTV++, d_maskPrior3);
			else
				kernelProxTVGrad.setArg(kernelIndCPTV++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelProxTVGrad.setArg(kernelIndCPTV++, d_eFOVIndices);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVGrad, cl::NullRange, globalPriorEFOV, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TV gradient kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
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
	inline int ProxTGVSymmDeriv(const scalarStruct& inputScalars, const float sigma2) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TGV symmetric derivative");
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTGV = 0ULL;
		if (inputScalars.largeDim)
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
			mexEval();
		}
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_N[0]);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_NPrior);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_vX);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_vY);
		if (!inputScalars.TGV2D)
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_vZ);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rX);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rY);
		if (!inputScalars.TGV2D) {
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rZ);
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rXY);
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rXZ);
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rYZ);
		}
		else
			kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_rXY);
		kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, sigma2);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_maskPrior3);
			else
				kernelProxTGVSymmDeriv.setArg(kernelIndCPTGV++, d_maskPrior);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVSymmDeriv, cl::NullRange, globalPriorEFOV, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV symmetric derivative kernel\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
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
	inline int ProxTGVDiv(const scalarStruct& inputScalars,	const float theta, const float tau) {
		if (inputScalars.verbose >= 3) {
			mexPrint("Starting Proximal TGV divergence");
		}
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTGV = 0ULL;
		if (inputScalars.largeDim)
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		status = (CLCommandQueue[0]).finish();
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_N[0]);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_NPrior);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rX);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rY);
		if (!inputScalars.TGV2D) {
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rZ);
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rXY);
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rXZ);
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rYZ);
		}
		else
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_rXY);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_vX);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_vY);
		if (!inputScalars.TGV2D)
			kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_vZ);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_qX);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_qY);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_qZ);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, theta);
		kernelProxTGVDiv.setArg(kernelIndCPTGV++, tau);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_maskPrior3);
			else
				kernelProxTGVDiv.setArg(kernelIndCPTGV++, d_maskPrior);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Proximal TGV divergence kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
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
	inline int elementWiseComp(const bool mult, const uint64_t size[], const bool D2 = false) {
		cl::NDRange gSize = { static_cast<cl::size_type>(size[0]), static_cast<cl::size_type>(size[1]), static_cast<cl::size_type>(size[2])};
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndE = 0ULL;
		status = (CLCommandQueue[0]).finish();
		if (mult) {
			kernelElementMultiply.setArg(kernelIndE++, d_vector);
			kernelElementMultiply.setArg(kernelIndE++, d_input);
			kernelElementMultiply.setArg(kernelIndE++, static_cast<cl_uchar>(D2));
			// Compute the kernel
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementMultiply, cl::NullRange, gSize, cl::NullRange);
		}
		else {
			kernelElementDivision.setArg(kernelIndE++, d_vector);
			kernelElementDivision.setArg(kernelIndE++, d_input);
			// Compute the kernel
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementDivision, cl::NullRange, gSize, cl::NullRange);
		}
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the element-wise kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
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
	inline int hyperGradient(const scalarStruct& inputScalars, const float sigma, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL hyperbolic prior gradient computation");
		cl_int status = CL_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
		status = (CLCommandQueue[0]).finish();
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		cl_uint kernelIndHyper = 0ULL;
		kernelHyper.setArg(kernelIndHyper++, d_W);
		if (inputScalars.useImages) {
			kernelHyper.setArg(kernelIndHyper++, d_inputI);
		}
		else {
			kernelHyper.setArg(kernelIndHyper++, d_inputB);
		}
		kernelHyper.setArg(kernelIndHyper++, d_N[0]);
		kernelHyper.setArg(kernelIndHyper++, d_NOrig);
		kernelHyper.setArg(kernelIndHyper++, sigma);
		kernelHyper.setArg(kernelIndHyper++, inputScalars.epps);
		kernelHyper.setArg(kernelIndHyper++, beta);
		kernelHyper.setArg(kernelIndHyper++, d_weights);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelHyper.setArg(kernelIndHyper++, d_maskPrior3);
			else
				kernelHyper.setArg(kernelIndHyper++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelHyper.setArg(kernelIndHyper++, d_eFOVIndices);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelHyper, cl::NullRange, globalPrior, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the hyperbolic prior gradient kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after hyperbolic prior gradient kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL hyperbolic prior gradient computed");
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
	inline int TVGradient(const scalarStruct& inputScalars, const float sigma, const float smooth, const float beta, const float C = 0.f, const int type = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL TV gradient computation");
		cl_int status = CL_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
		status = (CLCommandQueue[0]).finish();
		//cl::detail::size_t_array region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("smooth = %f\n", smooth);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		cl_uint kernelIndTV = 0ULL;
		kernelTV.setArg(kernelIndTV++, d_W);
		if (inputScalars.useImages) {
			kernelTV.setArg(kernelIndTV++, d_inputI);
		}
		else {
			kernelTV.setArg(kernelIndTV++, d_inputB);
		}
		kernelTV.setArg(kernelIndTV++, d_N[0]);
		kernelTV.setArg(kernelIndTV++, d_NOrig);
		kernelTV.setArg(kernelIndTV++, sigma);
		kernelTV.setArg(kernelIndTV++, smooth);
		kernelTV.setArg(kernelIndTV++, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.maskBPZ > 1)
				kernelTV.setArg(kernelIndTV++, d_maskPrior3);
			else
				kernelTV.setArg(kernelIndTV++, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			kernelTV.setArg(kernelIndTV++, d_eFOVIndices);
		if (type == 2 || type == 3)
			kernelTV.setArg(kernelIndTV++, C);
		if (type > 0)
			kernelTV.setArg(kernelIndTV++, d_refIm);
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelTV, cl::NullRange, globalPrior, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the TV gradient kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after TV gradient kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL TV gradient computed");
		return 0;
	}


	inline int PoissonUpdate(const scalarStruct& inputScalars, const float lambda, const float epps, const float alpha, const int ii = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL Poisson update (PKMA/MBSREM/BSREM) computation");
		cl_int status = CL_SUCCESS;
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndPoisson = 0ULL;
		global = { inputScalars.Nx[ii] + erotusPDHG[0][ii], inputScalars.Ny[ii] + erotusPDHG[1][ii], inputScalars.Nz[ii] };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].s[2]);
			mexPrintBase("lambda = %.8f\n", lambda);
			mexPrintBase("alpha = %f\n", alpha);
			mexEval();
		}
		kernelPoisson.setArg(kernelIndPoisson++, d_im);
		kernelPoisson.setArg(kernelIndPoisson++, d_rhs);
		kernelPoisson.setArg(kernelIndPoisson++, d_N[ii]);
		kernelPoisson.setArg(kernelIndPoisson++, lambda);
		kernelPoisson.setArg(kernelIndPoisson++, epps);
		kernelPoisson.setArg(kernelIndPoisson++, alpha);
		kernelPoisson.setArg(kernelIndPoisson++, static_cast<cl_uchar>(inputScalars.enforcePositivity));
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPoisson, cl::NullRange, global, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the Poisson update kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after Poisson update kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL Poisson update computed");
		return 0;
	}

	inline int PDHGUpdate(const scalarStruct& inputScalars, const float epps, const float theta, const float tau, const int ii = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenCL PDHG update computation");
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndPDHG = 0ULL;
		global = { inputScalars.Nx[ii] + erotusPDHG[0][ii], inputScalars.Ny[ii] + erotusPDHG[1][ii], inputScalars.Nz[ii] };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].s[2]);
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		kernelPDHG.setArg(kernelIndPDHG++, d_im);
		kernelPDHG.setArg(kernelIndPDHG++, d_rhs);
		kernelPDHG.setArg(kernelIndPDHG++, d_U);
		kernelPDHG.setArg(kernelIndPDHG++, d_N[ii]);
		kernelPDHG.setArg(kernelIndPDHG++, epps);
		kernelPDHG.setArg(kernelIndPDHG++, theta);
		kernelPDHG.setArg(kernelIndPDHG++, tau);
		kernelPDHG.setArg(kernelIndPDHG++, static_cast<cl_uchar>(inputScalars.enforcePositivity));
		// Compute the kernel
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPDHG, cl::NullRange, global, localPrior);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to launch the PDHG update kernel\n");
			return -1;
		}

		status = (CLCommandQueue[0]).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrint("Queue finish failed after PDHG update kernel\n");
			return -1;
		}
		if (inputScalars.verbose >= 3)
			mexPrint("OpenCL PDHG update computed");
		return 0;
	}

};