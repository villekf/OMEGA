/*******************************************************************************************************************************************
* Class object for forward and backward projections.
*
* Copyright (C) 2022 Ville-Veikko Wettenhovi
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
/// Class object for forward and backward projections
/// </summary>
class ProjectorClass {
//private:
	// Local size
	size_t local_size[2];
	// Kernel input indices
	cl_uint kernelInd_MRAMLA = 0;
	cl_uint kernelIndFP = 0;
	cl_uint kernelIndBP = 0;
	cl_uint kernelIndFPSubIter = 0;
	cl_uint kernelIndBPSubIter = 0;
	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
	cl_float3 b, d, bmax;
	// Crystal pitch
	cl_float2 dPitch;
	// Image dimensions
	cl_uint3 d_N;
	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusBP[3];
	// Local and global sizes
	cl::NDRange local, global;
	// Image origin
	cl::detail::size_t_array origin = { 0, 0, 0 };
	// Image format
	cl::ImageFormat format;
	cl::ImageFormat formatMask;

	/// <summary>
	/// This function creates the OpenCL programs for the MBSREM/COSEM/etc. prepass, for forward and backward projections and for NLM/MRP/RDP
	/// </summary>
	inline cl_int createProgram(const char* k_path, cl::Context& af_context, cl::Device& af_device_id, cl::Program& programFP, cl::Program& programBP,
		cl::Program& programMBSREM, cl::Program& programAux, const char* header_directory, scalarStruct& inputScalars, const bool find_lors,
		const RecMethods MethodList, const Weighting& w_vec, const size_t local_size[]) {

		cl_int status = CL_SUCCESS;

		//std::string options = header_directory;
		//options += " -cl-single-precision-constant";

		std::string kernelFile = header_directory;
		std::string kernel_path, kernel_pathBP;
		std::string contentFP, contentBP;
		std::string contentAux;
		std::string options = "-cl-single-precision-constant";
		//options += " -cl-fast-relaxed-math";
		std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
		// Load the header text file
		std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
		if ((inputScalars.tube_width > 0.f && (inputScalars.projector_type == 2u || inputScalars.projector_type == 22u))
			|| inputScalars.projector_type == 3u || inputScalars.projector_type == 33u) {
			if (inputScalars.orthXY)
				options += " -DCRYSTXY";
			if (inputScalars.orthZ)
				options += " -DCRYSTZ";
			std::ifstream sourceHeader1(kernelFile + "general_orth_opencl_functions.h");
			std::string contentHeader1((std::istreambuf_iterator<char>(sourceHeader1)), std::istreambuf_iterator<char>());
			std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
			std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
			contentHeader += contentHeader1 + contentHeader3;
		}

		kernel_path = kernelFile;
		kernel_pathBP = kernelFile;
		//kernel_path += ".cl";
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u 
			|| inputScalars.projector_type == 33u || inputScalars.projector_type == 14) {
			if (!inputScalars.precompute && (inputScalars.n_rays * inputScalars.n_rays3D) > 1)
				kernel_path += "multidevice_siddon_no_precomp.cl";
			else
				kernel_path += "multidevice_kernel.cl";
		}
		else if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41)
			kernel_path += "projectorType4.cl";
		else if (inputScalars.projector_type == 5)
			kernel_path += "projectorType5.cl";
		std::ifstream sourceFile(kernel_path.c_str());
		std::string contentFFP((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
		contentFP = contentHeader + contentFFP;
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 11 || inputScalars.projector_type == 41 || 
			inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			if (!inputScalars.precompute && (inputScalars.n_rays * inputScalars.n_rays3D) > 1)
				kernel_pathBP += "multidevice_siddon_no_precomp.cl";
			else
				kernel_pathBP += "multidevice_kernel.cl";
		}
		else if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14)
			kernel_pathBP += "projectorType4.cl";
		else if (inputScalars.projector_type == 5)
			kernel_pathBP += "projectorType5.cl";
		std::ifstream sourceFileBP(kernel_pathBP.c_str());
		std::string contentFBP((std::istreambuf_iterator<char>(sourceFileBP)), std::istreambuf_iterator<char>());
		contentBP = contentHeader + contentFBP;
		
		// Load the source text file
		// Set all preprocessor definitions
		const bool siddonVal = (inputScalars.projector_type == 1u || inputScalars.projector_type == 11u 
			|| inputScalars.projector_type == 41u || inputScalars.projector_type == 14u) ? true : false;
		if (inputScalars.projector_type == 3u || inputScalars.projector_type == 33u)
			options += " -DVOL";
		if (inputScalars.precompute)
			options += " -DPRECOMPUTE";
		if (inputScalars.raw == 1)
			options += " -DRAW";
		if (inputScalars.maskFP)
			options += " -DMASKFP";
		else if (inputScalars.maskBP)
			options += " -DMASKBP";
		if (inputScalars.projector_type == 2u || inputScalars.projector_type == 22u || inputScalars.projector_type == 3u || inputScalars.projector_type == 33u)
			options += " -DORTH";
		if (inputScalars.attenuation_correction == 1u)
			options += " -DATN";
		if (inputScalars.normalization_correction == 1u)
			options += " -DNORM";
		if (inputScalars.scatter == 1u)
			options += " -DSCATTER";
		if (inputScalars.randoms_correction == 1u)
			options += " -DRANDOMS";
		if (inputScalars.nLayers > 1U)
			options += " -DNLAYERS";
		if (inputScalars.TOF && siddonVal) {
			options += " -DTOF";
		}
		if (inputScalars.CT)
			options += " -DCT";
		else if (inputScalars.PET)
			options += " -DPET";
		else if (inputScalars.SPECT)
			options += " -DSPECT";
		options += (" -DNBINS=" + std::to_string(inputScalars.nBins));
		if (inputScalars.listmode == 1)
			options += " -DLISTMODE";
		else if (inputScalars.listmode == 2)
			options += " -DLISTMODE2";
		if (siddonVal && !inputScalars.precompute && (inputScalars.n_rays * inputScalars.n_rays3D) > 1) {
			options += (" -DN_RAYS=" + std::to_string(inputScalars.n_rays * inputScalars.n_rays3D));
			options += (" -DN_RAYS2D=" + std::to_string(inputScalars.n_rays));
			options += (" -DN_RAYS3D=" + std::to_string(inputScalars.n_rays3D));
		}
		if (find_lors)
			options += " -DFIND_LORS";
		if (inputScalars.PITCH)
			options += " -DPITCH";
		if (((inputScalars.subsets > 1 && inputScalars.subsetType < 8)) && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
			options += " -DSUBSETS";
		//if (projector_type == 1u && use_psf && (precompute || (n_rays * n_rays3D) == 1)) {
		//	options += " -DORTH";
		//	options += " -DCRYSTZ";
		//	uint32_t limit = static_cast<uint32_t>(std::floor(cr_pz / dx));
		//	options += (" -DPSF_LIMIT=" + std::to_string(limit));
		//	options += (" -DX=" + std::to_string(dx));
		//	options += (" -DSIGMA=" + std::to_string(cr_pz / 2.355f));
		//}
		if (local_size[1] > 0ULL) {
			options += (" -DLOCAL_SIZE=" + std::to_string(local_size[0]));
			options += (" -DLOCAL_SIZE2=" + std::to_string(local_size[1]));
		}
		else
			options += (" -DLOCAL_SIZE=" + std::to_string(local_size[0]));
		if (DEBUG) {
			mexPrintf("path = %s\n", kernel_path.c_str());
			mexPrintf("pathBP = %s\n", kernel_pathBP.c_str());
			//mexPrintf("%s\n", options.c_str());
			mexPrintf("file = %s\n", kernelFile.c_str());
			mexEvalString("pause(.0001);");
		}
		// Build subset-based program
		if (inputScalars.projector_type != 4 && inputScalars.projector_type != 5 && inputScalars.projector_type != 6) {
			std::string os_options = options;
			os_options += " -DAF";
			if ((inputScalars.projector_type == 2 || inputScalars.projector_type == 3u || inputScalars.TOF) && inputScalars.dec > 0)
				os_options += (" -DDEC=" + std::to_string(inputScalars.dec));
			os_options += (" -DN_REKOS=" + std::to_string(inputScalars.nRekos));
			if (inputScalars.nRekos == 1)
				os_options += " -DNREKOS1";
			else if (inputScalars.nRekos == 2)
				os_options += " -DNREKOS2";
			os_options += " -DSIDDON";
			std::string os_optionsFP = os_options;
			os_optionsFP += " -DFP";
			if (inputScalars.projector_type < 4) {
				os_optionsFP += " -DBP";
				if (MethodList.MRAMLA || MethodList.MBSREM)
					os_optionsFP += " -DMRAMLA";
				if (MethodList.COSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0u || MethodList.ECOSEM)
					os_optionsFP += " -DCOSEM";
			}

			//status = buildProgram(verbose, k_path, af_context, af_device_id, program_os, atomic_64bit, os_options);
			if (inputScalars.projector_type == 11 || inputScalars.projector_type == 14 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33)
				status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
			else if (inputScalars.projector_type != 41)
				status = buildProgram(inputScalars.verbose, contentBP, af_context, af_device_id, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
			if (inputScalars.projector_type == 11 || inputScalars.projector_type == 41 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33) {
				//os_options = options;
				//os_options += " -DSIDDON";
				os_options += " -DBP";
				status = buildProgram(inputScalars.verbose, contentBP, af_context, af_device_id, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			}
		}
		// Build the prepass phase program
		//if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && (w_vec.MBSREM_prepass ||
		//	MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && inputScalars.projector_type != 4 &&
		//	inputScalars.projector_type != 5 && inputScalars.projector_type != 6) {
		//	options += " -DAF";
		//	options += " -DMBSREM";
		//	if (inputScalars.nRekos == 1)
		//		options += " -DNREKOS1";
		//	else if (inputScalars.nRekos == 2)
		//		options += " -DNREKOS2";
		//	if (MethodList.MRAMLA || MethodList.MBSREM)
		//		options += " -DMRAMLA";
		//	//status = buildProgram(verbose, k_path, af_context, af_device_id, program_mbsrem, atomic_64bit, options);
		//	status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programMBSREM, inputScalars.atomic_64bit, inputScalars.atomic_32bit, options);
		//}
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14 || inputScalars.projector_type == 41) {
			std::string os_options = options;
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41)
				os_options += " -DFP";
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14)
				os_options += " -DBP";
			if (inputScalars.projector_type == 41)
				os_options += " -DPTYPE41";
			os_options += " -DAF";
			os_options += " -DPTYPE4";
			os_options += (" -DNVOXELS=" + std::to_string(NVOXELS));
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41)
				status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			else
				status = buildProgram(inputScalars.verbose, contentBP, af_context, af_device_id, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
		}
		else if (inputScalars.projector_type == 5) {
			std::string os_options = options;
			os_options += " -DPROJ5";
			if (inputScalars.meanFP)
				os_options += " -DMEANDISTANCEFP";
			else if (inputScalars.meanBP)
				os_options += " -DMEANDISTANCEBP";
			os_options += " -DFP";
			os_options += " -DBP";
			os_options += " -DAF";
			status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
		}

		//int nRekos = static_cast<int>(MethodList.LSQR) + static_cast<int>(MethodList.CP);
		//if (nRekos > 0 && inputScalars.projector_type < 4) {
		//	std::string os_options = options;
		//	if ((inputScalars.projector_type == 2 || inputScalars.projector_type == 3u || inputScalars.TOF) && inputScalars.dec > 0)
		//		os_options += (" -DDEC=" + std::to_string(inputScalars.dec));
		//	os_options += (" -DN_REKOS=" + std::to_string(nRekos));
		//	os_options += " -DAF";
		//	if (nRekos == 1)
		//		os_options += " -DNREKOS1";
		//	else
		//		os_options += " -DNREKOS2";
		//	std::string os_optionsFP = os_options;
		//	os_optionsFP += " -DFP";

		//	//status = buildProgram(verbose, k_path, af_context, af_device_id, program_os, atomic_64bit, os_options);
		//	status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programOSFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
		//	std::string os_optionsBP = os_options;
		//	os_optionsBP += " -DBP";
		//	status = buildProgram(inputScalars.verbose, contentFP, af_context, af_device_id, programOSBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsBP);
		//}
		if (MethodList.NLM || MethodList.MRP) {
			std::string auxKernelPath = kernelFile + "auxKernels.cl";
			std::ifstream sourceFileAux(auxKernelPath.c_str());
			std::string contentAux((std::istreambuf_iterator<char>(sourceFileAux)), std::istreambuf_iterator<char>());
			options = "-cl-single-precision-constant";
			options += " -DAF";
			if (MethodList.MRP) {
				options += " -DMEDIAN";
				options += (" -DSEARCH_WINDOW_X=" + std::to_string(w_vec.Ndx));
				options += (" -DSEARCH_WINDOW_Y=" + std::to_string(w_vec.Ndy));
				options += (" -DSEARCH_WINDOW_Z=" + std::to_string(w_vec.Ndz));
			}
			if (MethodList.NLM) {
				options += " -DNLM_";
			}
			status = buildProgram(inputScalars.verbose, contentAux, af_context, af_device_id, programAux, inputScalars.atomic_64bit, inputScalars.atomic_32bit, options);
		}
		return status;
	}

	/// <summary>
	/// Builds the input OpenCL program
	/// </summary>
	/// <param name="verbose the level of verbosity"></param>
	/// <param name="contentFP program code"></param>
	/// <param name="af_context OpenCL context"></param>
	/// <param name="af_device_id OpenCL device ID"></param>
	/// <param name="program the program where to store the built program"></param>
	/// <param name="atomic_64bit are 64-bit (int64) atomics used"></param>
	/// <param name="atomic_32bit are 32-bit (int) atomics used"></param>
	/// <param name="options preprocessor values for the build"></param>
	/// <returns></returns>
	inline cl_int buildProgram(const bool verbose, std::string contentFP, cl::Context& af_context, cl::Device& af_device_id, cl::Program& program,
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
		else
			options += " -DCAST=float";
		if (DEBUG)
			mexPrintf("%s\n", options.c_str());
		if (atomic_64bit) {
			cl::string apu = af_device_id.getInfo<CL_DEVICE_EXTENSIONS>();
			cl::string apu2 = "cl_khr_int64_base_atomics";
			size_t var = apu.find(apu2);
			if (var < 0) {
				options.erase(pituus, options.size() + 1);
				options += " -DCAST=float";
				status = -1;
			}
			else {
				//std::string kernel_path_atom;

				//kernel_path_atom = k_path;
				//kernel_path_atom += ".cl";
				//// Load the source text file
				//std::ifstream sourceFile_atom(kernel_path_atom.c_str());
				//std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
				std::vector<std::string> testi;
				testi.push_back(contentFP);
				cl::Program::Sources source(testi);
				program = cl::Program(af_context, source);
				//try {
				status = program.build(options.c_str());
				if (status == CL_SUCCESS) {
					mexPrintf("OpenCL program (64-bit atomics) built\n");
				}
				//}
				//catch (cl::Error& e) {
					//mexPrintf("%s\n", e.what());
				else {
					mexPrintf("Failed to build 64-bit atomics program.\n");
					if (DEBUG) {
						std::vector<cl::Device> dev;
						af_context.getInfo(CL_CONTEXT_DEVICES, &dev);
						for (int ll = 0; ll < dev.size(); ll++) {
							cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
							if (status != CL_BUILD_ERROR)
								continue;
							std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
							std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
							mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
						}
						//return -1;
					}
					options.erase(pituus, options.size() + 1);
					options += " -DCAST=float";
					//status = -1;
				}
				//}
			}
		}
		else
			status = -1;
		// If not, use 32-bit atomic add (float)
		if (status != CL_SUCCESS) {
			status = CL_SUCCESS;
			atomic_64bit = false;

			//std::string kernel_path;

			//kernel_path = k_path;
			//kernel_path += ".cl";
			//std::fstream sourceFile(kernel_path.c_str());
			//std::string contentFP((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
			std::vector<std::string> testi;
			testi.push_back(contentFP);
			cl::Program::Sources source(testi);
			program = cl::Program(af_context, source);
			//try {
			status = program.build(options.c_str());
			if (status == CL_SUCCESS) {
				mexPrintf("OpenCL program built\n");
			}
			//}
			//catch (cl::Error& e) {
			else {
				mexPrintf("Failed to build OpenCL program.\n");
				std::vector<cl::Device> dev;
				af_context.getInfo(CL_CONTEXT_DEVICES, &dev);
				for (int ll = 0; ll < dev.size(); ll++) {
					cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
					if (status != CL_BUILD_ERROR)
						continue;
					std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
					std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
					mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
				}
			}
			//mexPrintf("%s\n", e.what());
			//status = -1;
		//}
		}
		return status;
	}

	/// <summary>
	/// Creates the necessary OpenCL kernels from the input programs
	/// </summary>
	/// <param name="kernelFP forward projection kernel"></param>
	/// <param name="kernelBP backprojection kernel"></param>
	/// <param name="kernelMBSREM prepass kernel"></param>
	/// <param name="kernelNLM NLM kernel"></param>
	/// <param name="kernelMed MRP kernel"></param>
	/// <param name="kernelRDP RDP kernel"></param>
	/// <param name="programFP program containing forward projection"></param>
	/// <param name="programBP program containing backprojection"></param>
	/// <param name="programMBSREM program containing the prepass phase"></param>
	/// <param name="programAux program containing NLM/MRP/RDP"></param>
	/// <param name="MethodList reconstruction algorithms selected"></param>
	/// <param name="w_vec parameters needed"></param>
	/// <param name="inputScalars scalars needed"></param>
	/// <returns></returns>
	inline cl_int createKernels(cl::Kernel& kernelFP, cl::Kernel& kernelBP, cl::Kernel& kernelMBSREM, cl::Kernel& kernelNLM, cl::Kernel& kernelMed,
		cl::Kernel& kernelRDP, const cl::Program& programFP, const cl::Program& programBP, const cl::Program& programMBSREM,
		const cl::Program& programAux, const RecMethods& MethodList, const Weighting& w_vec, const scalarStruct& inputScalars) {
		cl_int status = CL_SUCCESS;
		// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14 || inputScalars.projector_type == 41) {
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41) {
				kernelFP = cl::Kernel(programFP, "projectorType4Forward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Failed to create projection type 4 FP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose == 2) {
					mexPrintf("OpenCL kernel for projection type 4 FP successfully created\n");
					mexEvalString("pause(.0001);");
				}
			}
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14) {
				if (inputScalars.projector_type == 4)
					kernelBP = cl::Kernel(programFP, "projectorType4Backward", &status);
				else if (inputScalars.projector_type == 14)
					kernelBP = cl::Kernel(programBP, "projectorType4Backward", &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Failed to create projection type 4 BP kernel\n");
					return -1;
				}
				else if (DEBUG || inputScalars.verbose == 2) {
					mexPrintf("OpenCL kernel for projection type 4 BP successfully created\n");
					mexEvalString("pause(.0001);");
				}
			}
		}
		if (inputScalars.projector_type == 5) {
			kernelFP = cl::Kernel(programFP, "projectorType5Forward", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create projection type 5 FP kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose == 2) {
				mexPrintf("OpenCL kernel for projection type 5 FP successfully created\n");
				mexEvalString("pause(.0001);");
			}
			kernelBP = cl::Kernel(programFP, "projectorType5Backward", &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create projection type 5 BP kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose == 2) {
				mexPrintf("OpenCL kernel for projection type 5 BP successfully created\n");
				mexEvalString("pause(.0001);");
			}
		}
		if (inputScalars.projector_type != 4 && inputScalars.projector_type != 5 && inputScalars.projector_type != 6) {


			if (inputScalars.projector_type == 1 && (inputScalars.n_rays * inputScalars.n_rays3D) > 1) {
				kernelBP = cl::Kernel(programFP, "proj1SiddonMultiRay", &status);
			}
			else if (inputScalars.projector_type == 14 || inputScalars.projector_type == 11)
				kernelFP = cl::Kernel(programFP, "proj123SiddonSingleRay", &status);
			if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11)
				kernelBP = cl::Kernel(programBP, "proj123SiddonSingleRay", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create OS-methods kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose == 2) {
				mexPrintf("OpenCL kernel successfully created\n");
				mexEvalString("pause(.0001);");
			}
		}

		// Kernel for the prepass phase needed for MRAMLA, MBSREM, RBI, COSEM, ACOSEM and ECOSEM
		//if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && (w_vec.MBSREM_prepass ||
		//	MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && inputScalars.projector_type != 4 &&
		//	inputScalars.projector_type != 5 && inputScalars.projector_type != 6) {

		//	// Create the prepass kernel
		//	if (inputScalars.projector_type == 1 && (inputScalars.n_rays * inputScalars.n_rays3D) > 1)
		//		kernelMBSREM = cl::Kernel(programMBSREM, "proj1SiddonMultiRay", &status);
		//	else
		//		kernelMBSREM = cl::Kernel(programMBSREM, "proj123SiddonSingleRay", &status);

		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		mexPrintf("Failed to create prepass kernel\n");
		//		return -1;
		//	}
		//	else if (DEBUG || inputScalars.verbose == 2) {
		//		mexPrintf("Prepass kernel successfully created\n");
		//		mexEvalString("pause(.0001);");
		//	}
		//}
		if (MethodList.NLM) {
			kernelNLM = cl::Kernel(programAux, "NLM", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create NLM kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose == 2) {
				mexPrintf("NLM kernel successfully created\n");
				mexEvalString("pause(.0001);");
			}
		}
		if (MethodList.MRP) {
			kernelMed = cl::Kernel(programAux, "medianFilter3D", &status);

			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Failed to create Median kernel\n");
				return -1;
			}
			else if (DEBUG || inputScalars.verbose == 2) {
				mexPrintf("Median kernel successfully created\n");
				mexEvalString("pause(.0001);");
			}
		}
		return status;
	}
public:
	cl::Context af_context;
	cl::Device af_device_id;
	cl::CommandQueue af_queue;
	OpenCL_im_vectors vec_opencl;
	cl::Kernel kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP;
	cl::Buffer d_xcenter, d_ycenter, d_zcenter, d_reko_type, d_V, d_TOFCenter, d_Summ, d_output, d_meanBP, d_meanFP;
	//cl::Buffer d_indices, d_values, d_rowInd;
	cl::Image2D d_maskFP, d_maskBP;
	cl::Image3D d_inputImage, d_imageX, d_imageY, d_atten, d_uref;

	std::vector<cl::Buffer> d_lor;
	std::vector<cl::Buffer> d_L;
	std::vector<cl::Buffer> d_zindex;
	std::vector<cl::Buffer> d_xyindex;
	std::vector<cl::Buffer> d_Sino;
	std::vector<cl::Buffer> d_sc_ra;
	std::vector<cl::Buffer> d_norm;
	std::vector<cl::Buffer> d_scat;
	std::vector<cl::Buffer> d_x;
	std::vector<cl::Buffer> d_z;
	cl_uchar no_norm = 0;
	// Create the projector object
	inline int addProjector(scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, const char* k_path,
		const char* header_directory, const bool find_lors) {
		// Set-up the local group size
		local_size[0] = 64ULL;
		local_size[1] = 1ULL;
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33 || inputScalars.projector_type == 41)
			local_size[0] = 128ULL;
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
			local_size[0] = 8ULL;
			local_size[1] = 8ULL;
		}
		cl_int status = CL_SUCCESS;

		// Create the OpenCL context and command queue and assign the device
		af_context = afcl::getContext(true);
		std::vector<cl::Device> devices = af_context.getInfo<CL_CONTEXT_DEVICES>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		af_device_id = devices[0];
		af_queue = afcl::getQueue(true);
		// For NVIDIA cards, 32 local size seems more optimal with 1D kernelFP
		std::string deviceName = af_device_id.getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		if (NV.compare(deviceName) == 0 && (inputScalars.projector_type == 1 || inputScalars.projector_type == 11) && local_size[1] == 1ULL)
			local_size[0] = 32ULL;

		cl::Program programFP, programBP, programMBSREM, programAux;

		status = createProgram(k_path, af_context, af_device_id, programFP, programBP, programMBSREM, programAux, header_directory,
			inputScalars, find_lors, MethodList, w_vec, local_size);
		if (status != CL_SUCCESS) {
			std::cerr << "Error while creating program" << std::endl;
			return -1;
		}
		else if (DEBUG || inputScalars.verbose == 2) {
			mexPrintf("OpenCL programs successfully created\n");
			mexEvalString("pause(.0001);");
		}

		status = createKernels(kernelFP, kernelBP, kernelMBSREM, kernelNLM, kernelMed, kernelRDP, programFP, programBP, programMBSREM,
			programAux, MethodList, w_vec, inputScalars);
		if (status != CL_SUCCESS) {
			mexPrintf("Failed to create kernels\n");
			return -1;
		}
		else if (DEBUG || inputScalars.verbose == 2) {
			mexPrintf("OpenCL kernels successfully created\n");
			mexEvalString("pause(.0001);");
		}
		format.image_channel_order = CL_R;
		format.image_channel_data_type = CL_FLOAT;
		formatMask.image_channel_order = CL_R;
		formatMask.image_channel_data_type = CL_UNSIGNED_INT8;

		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			erotus[0] = w_vec.size_x % local_size[0];
			erotus[1] = w_vec.size_y % local_size[1];
			if (erotus[1] > 0)
				erotus[1] = (local_size[1] - erotus[1]);
			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
		}


		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type == 14 || MethodList.NLM) {
			erotusBP[0] = inputScalars.Nx % local_size[0];
			erotusBP[1] = inputScalars.Ny % local_size[1];
			if (erotusBP[0] > 0)
				erotusBP[0] = (local_size[0] - erotusBP[0]);
			if (erotusBP[1] > 0)
				erotusBP[1] = (local_size[1] - erotusBP[1]);
		}
		local = { local_size[0] , local_size[1] };
		//kernelIndBP = 0;
		//kernelIndFP = 0;
		//kernelIndFPSubIter = 0;
		//kernelIndBPSubIter = 0;
		//origin = { 0, 0, 0 };
		return 0;
	}

	inline cl_int createAndWriteBuffers(const std::vector<int64_t>& length, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* lor1, const uint16_t* L, const float* Sin, const int64_t* pituus, const float* atten, 
		const float* norm, const float* scat, const float* V, const float* x_center, const float* y_center, const float* z_center, 
		const float* sc_ra, const uint8_t* reko_type, const scalarStruct& inputScalars, const float* TOFCenter, const Weighting& w_vec, 
		const RecMethods& MethodList) {
		cl_int status = CL_SUCCESS;
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(w_vec.size_x) * static_cast<size_t>(w_vec.size_y);
		// Create the necessary buffers
		// Detector coordinates
		d_V = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_V, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (!inputScalars.CT && inputScalars.listmode == 0) {
			d_x[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		//d_y = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_of_x, NULL, &status);
		//if (status != CL_SUCCESS) {
		//	getErrorString(status);
		//	return -1;
		//}
		d_xcenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_center_x, NULL, &status);;
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		d_ycenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_center_y, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		d_zcenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_center_z, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		// Attenuation data for image-based attenuation
		cl::size_type imX = inputScalars.Nx;
		cl::size_type imY = inputScalars.Ny;
		cl::size_type imZ = inputScalars.Nz;
		if (inputScalars.attenuation_correction) {
			d_atten = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (w_vec.NLM_anatomical && MethodList.NLM) {
			d_uref = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.maskFP || inputScalars.maskBP) {
			if (inputScalars.maskFP) {
				imX = w_vec.size_x;
				imY = w_vec.size_y;
				d_maskFP = cl::Image2D(af_context, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
			}
			if (inputScalars.maskBP) {
				imX = inputScalars.Nx;
				imY = inputScalars.Ny;
				if (DEBUG) {
					mexPrintf("imX = %u\n", imX);
					mexPrintf("imY = %u\n", imY);
					mexEvalString("pause(.0001);");
				}
				d_maskBP = cl::Image2D(af_context, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		// TOF bin centers
		d_TOFCenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nBins, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			d_reko_type = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.nRekos, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
			if (DEBUG) {
				mexPrintf("length[kk] = %u\n", length[kk]);
				mexPrintf("kk = %u\n", kk);
				mexEvalString("pause(.0001);");
			}
			if (inputScalars.CT && inputScalars.listmode != 1) {
				if (inputScalars.PITCH)
					d_z[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 6, NULL, &status);
				else
					d_z[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else {
				if (inputScalars.CT || inputScalars.PET)
					if (inputScalars.nLayers > 1)
						d_z[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 3, NULL, &status);
					else
						d_z[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2, NULL, &status);
				else if (kk == inputScalars.osa_iter0)
					d_z[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.CT || inputScalars.listmode > 0) {
				d_x[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 6, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			// How many voxels does each LOR traverse
			if (inputScalars.precompute)
				d_lor[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * vecSize, NULL, &status);
			//else
			//	d_lor[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.size_norm > 1) {
				d_norm[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
			}
			//else {
			//	d_norm[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			//}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.size_scat > 1) {
				d_scat[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
			}
			//else {
			//	d_scat[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			//}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			// Measurement data
			if (inputScalars.TOF && inputScalars.projector_type < 4) {
				if (inputScalars.loadTOF) {
					if ((inputScalars.PET) && inputScalars.listmode == 0)
						d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * w_vec.size_x * w_vec.size_y * length[kk] * inputScalars.nBins, NULL, &status);
					else
						d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * inputScalars.nBins, NULL, &status);
				}
				else {
					if (kk == 0) {
						if ((inputScalars.PET) && inputScalars.listmode == 0)
							d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * w_vec.size_x * w_vec.size_y * length[kk] * inputScalars.nBins, NULL, &status);
						else
							d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * inputScalars.nBins, NULL, &status);
					}
				}
			}
			else if (inputScalars.projector_type < 4) {
				d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.randoms_correction == 1u)
				d_sc_ra[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize, NULL, &status);
			//else
			//	d_sc_ra[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
			if (inputScalars.raw && inputScalars.listmode != 1) {
				//d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
				//d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
				d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && inputScalars.subsets > 1)) {
				d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				//d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
			}
			//else {
			//	d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//	d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//	d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//}
		}

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Buffer creation failed\n");
			mexEvalString("pause(.0001);");
			return -1;
		}
		else if (DEBUG) {
			mexPrintf("Buffer creation succeeded\n");
			mexEvalString("pause(.0001);");
		}


		// assign values to the buffers
		status = af_queue.enqueueWriteBuffer(d_V, CL_FALSE, 0, sizeof(float) * inputScalars.size_V, V);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (!inputScalars.CT && inputScalars.listmode == 0) {
			status = af_queue.enqueueWriteBuffer(d_x[0], CL_FALSE, 0, sizeof(float) * inputScalars.size_of_x, x);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		//status = af_queue.enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(float) * size_of_x, y);
		//if (status != CL_SUCCESS) {
		//	getErrorString(status);
		//	return -1;
		//}
		//if (CT) {
		//	status = af_queue.enqueueWriteBuffer(d_angles, CL_FALSE, 0, sizeof(float) * TotSinos, angles);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return -1;
		//	}
		//}
		status = af_queue.enqueueWriteBuffer(d_xcenter, CL_FALSE, 0, sizeof(float) * inputScalars.size_center_x, x_center);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		status = af_queue.enqueueWriteBuffer(d_ycenter, CL_FALSE, 0, sizeof(float) * inputScalars.size_center_y, y_center);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		status = af_queue.enqueueWriteBuffer(d_zcenter, CL_FALSE, 0, sizeof(float) * inputScalars.size_center_z, z_center);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.attenuation_correction) {
			cl::detail::size_t_array region = { 0, 0, 0 };
			region[0] = inputScalars.Nx;
			region[1] = inputScalars.Ny;
			region[2] = inputScalars.Nz;
			status = af_queue.enqueueWriteImage(d_atten, CL_FALSE, origin, region, 0, 0, atten);
			//status = commandQueues[i].enqueueWriteBuffer(d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (w_vec.NLM_anatomical && MethodList.NLM) {
			cl::detail::size_t_array region = { 0, 0, 0 };
			region[0] = inputScalars.Nx;
			region[1] = inputScalars.Ny;
			region[2] = inputScalars.Nz;
			status = af_queue.enqueueWriteImage(d_uref, CL_FALSE, origin, region, 0, 0, w_vec.NLM_ref);
		}
		if (inputScalars.TOF) {
			status = af_queue.enqueueWriteBuffer(d_TOFCenter, CL_FALSE, 0, sizeof(float) * inputScalars.nBins, TOFCenter);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		//status = af_queue.enqueueWriteBuffer(d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos);
		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (inputScalars.maskFP || inputScalars.maskBP) {
			cl::detail::size_t_array region = { 1, 1, 1 };
			if (inputScalars.maskFP) {
				region[0] = w_vec.size_x;
				region[1] = w_vec.size_y;
				status = af_queue.enqueueWriteImage(d_maskFP, CL_FALSE, origin, region, 0, 0, w_vec.maskFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskBP) {
				//uint8_t* testi = new uint8_t[inputScalars.Nx * inputScalars.Ny];
				//for (int tt = 0; tt < inputScalars.Ny * inputScalars.Nx; tt++)
				//	testi[tt] = static_cast<uint8_t>(1);
				region[0] = inputScalars.Nx;
				region[1] = inputScalars.Ny;
				if (DEBUG) {
					mexPrintf("region[0] = %u\n", region[0]);
					mexPrintf("region[1] = %u\n", region[1]);
					mexPrintf("region[2] = %u\n", region[2]);
					mexEvalString("pause(.0001);");
				}
				status = af_queue.enqueueWriteImage(d_maskBP, CL_FALSE, origin, region, 0, 0, w_vec.maskBP);
				//status = af_queue.enqueueWriteImage(d_maskBP, CL_TRUE, origin, region, 0, 0, testi);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				//delete[] testi;
			}
		}
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			status = af_queue.enqueueWriteBuffer(d_reko_type, CL_FALSE, 0, sizeof(uint8_t) * inputScalars.nRekos, reko_type);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
			if (inputScalars.CT && inputScalars.listmode != 1) {
				if (inputScalars.PITCH)
					status = af_queue.enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 6, &z_det[pituus[kk] * 6]);
				else
					status = af_queue.enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 2, &z_det[pituus[kk] * 2]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else {
				if (inputScalars.PET && inputScalars.listmode == 0)
					if (inputScalars.nLayers > 1)
						status = af_queue.enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 3, &z_det[pituus[kk] * 3]);
					else
						status = af_queue.enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 2, &z_det[pituus[kk] * 2]);
				else if (kk == inputScalars.osa_iter0)
					status = af_queue.enqueueWriteBuffer(d_z[kk], CL_FALSE, 0, sizeof(float) * inputScalars.size_z, z_det);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.CT || inputScalars.listmode > 0) {
				status = af_queue.enqueueWriteBuffer(d_x[kk], CL_FALSE, 0, sizeof(float) * length[kk] * 6, &x[pituus[kk] * 6]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.raw && inputScalars.listmode != 1) {
				//status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
				//status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t), z_index);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
				status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && inputScalars.subsets > 1)) {
				status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &z_index[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				//status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L);
				//if (status != CL_SUCCESS) {
				//	getErrorString(status);
				//	return -1;
				//}
			}
			//else {
			//	status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//	status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t), z_index);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//	status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return -1;
			//	}
			//}
			if (inputScalars.precompute)
				status = af_queue.enqueueWriteBuffer(d_lor[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * vecSize, &lor1[pituus[kk] * vecSize]);
			//else
			//	status = af_queue.enqueueWriteBuffer(d_lor[kk], CL_FALSE, 0, sizeof(uint16_t), lor1);
			status = af_queue.finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.size_norm > 1ULL) {
				status = af_queue.enqueueWriteBuffer(d_norm[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &norm[pituus[kk] * vecSize]);
			}
			//else {
			//	status = af_queue.enqueueWriteBuffer(d_norm[kk], CL_FALSE, 0, sizeof(float) * inputScalars.size_norm, norm);
			//}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.size_scat > 1ULL) {
				status = af_queue.enqueueWriteBuffer(d_scat[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &scat[pituus[kk] * vecSize]);
			}
			//else {
			//	status = af_queue.enqueueWriteBuffer(d_scat[kk], CL_FALSE, 0, sizeof(float) * inputScalars.size_scat, scat);
			//}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (DEBUG) {
				mexPrintf("length[kk] = %d\n", length[kk]);
				if (inputScalars.PITCH && kk > 0) {
					mexPrintf("z_det[pituus[kk] * 6] = %f\n", z_det[pituus[kk] * 6 - 1]);
					mexPrintf("pituus[kk] * 6 = %d\n", pituus[kk] * 6 - 1);
				}
				mexPrintf("pituus[kk] = %d\n", pituus[kk]);
				mexPrintf("erotus = %d\n", pituus[kk + 1] - pituus[kk]);
				mexPrintf("kk = %d\n", kk);
				mexEvalString("pause(.0001);");
			}
			if (inputScalars.TOF && inputScalars.projector_type < 4) {
				if (inputScalars.loadTOF) {
					for (int64_t to = std::int64_t{ 0 }; to < inputScalars.nBins; to++)
						if ((inputScalars.PET) && inputScalars.listmode == 0)
							status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to * w_vec.size_x * w_vec.size_y,
								sizeof(float) * length[kk], &Sin[pituus[kk] * w_vec.size_x * w_vec.size_y + inputScalars.koko * to]);
						else
							status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sin[pituus[kk] + inputScalars.koko * to]);
				}
				else {
					if (kk == inputScalars.osa_iter0) {
						for (int64_t to = std::int64_t{ 0 }; to < inputScalars.nBins; to++)
							if ((inputScalars.PET) && inputScalars.listmode == 0)
								status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sin[pituus[kk] + inputScalars.koko * to]);
							else
								status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to * w_vec.size_x * w_vec.size_y,
									sizeof(float) * length[kk], &Sin[pituus[kk] * w_vec.size_x * w_vec.size_y + inputScalars.koko * to]);
					}
				}
			}
			else if (inputScalars.listmode != 2 && inputScalars.projector_type < 4) {
				if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
					status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, 0, sizeof(float) * length[kk] * w_vec.size_x * w_vec.size_y,
						&Sin[pituus[kk] * w_vec.size_x * w_vec.size_y]);
				else
					status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, 0, sizeof(float) * length[kk], &Sin[pituus[kk]]);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.randoms_correction)
				status = af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_FALSE, 0, sizeof(float) * length[kk] * vecSize, &sc_ra[pituus[kk] * vecSize]);
			//else
			//	status = af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_FALSE, 0, sizeof(float), apu);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			status = af_queue.finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Buffer write failed\n");
			return -1;
		}
		else if (DEBUG) {
			mexPrintf("Buffer write succeeded\n");
			mexEvalString("pause(.0001);");
		}
		return 0;
	}

	inline int createBuffers(scalarStruct& inputScalars, Weighting& w_vec, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* lor1, const uint16_t* L, const int64_t* pituus, const float* atten, const float* norm, const float* scat, 
		const float* V, const float* x_center, const float* y_center, const float* z_center, const float* sc_ra, const float* TOFCenter,
		const std::vector<int64_t>& length, const float* Sino, const uint8_t* reko_type, const RecMethods& MethodList) {
		cl_int status = CL_SUCCESS;
		if (inputScalars.precompute)
			d_lor.resize(inputScalars.subsets);
		if (inputScalars.raw)
			d_L.resize(inputScalars.subsets);
		if (inputScalars.subsetType < 8 && inputScalars.subsets > 1) {
			d_xyindex.resize(inputScalars.subsets);
			d_zindex.resize(inputScalars.subsets);
		}
		if (inputScalars.randoms_correction)
			d_sc_ra.resize(inputScalars.subsets);
		if (inputScalars.normalization_correction)
			d_norm.resize(inputScalars.subsets);
		if (inputScalars.scatter)
			d_scat.resize(inputScalars.subsets);
		d_x.resize(inputScalars.subsets);
		d_z.resize(inputScalars.subsets);
		d_Sino.resize(inputScalars.TOFsubsets);

		cl::size_type imX = inputScalars.Nx;
		cl::size_type imY = inputScalars.Ny;
		cl::size_type imZ = inputScalars.Nz * inputScalars.nRekos2;
		if (inputScalars.projector_type == 5) {
			vec_opencl.d_image_os_int = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, imY + 1, imZ + 1, imX, 0, 0, NULL, &status);
			imX++;
			cl::size_type aY = imY;
			imY = imZ + 1;
			imZ = aY;
		}
		vec_opencl.d_image_os = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
		if (status != CL_SUCCESS) {
			mexPrintf("Failed to create input images\n");
			return -1;
		}

		status = createAndWriteBuffers(length, x, z_det, xy_index, z_index, lor1, L, Sino, pituus, atten, norm, scat, V, x_center, y_center, 
			z_center, sc_ra, reko_type, inputScalars, TOFCenter, w_vec, MethodList);
		if (status != CL_SUCCESS) {
			return -1;
		}
		return 0;
	}

	inline int initializeKernel(scalarStruct& inputScalars, Weighting& w_vec) {
		cl_int status = CL_SUCCESS;
		b = { inputScalars.bx, inputScalars.by, inputScalars.bz };
		d = { inputScalars.dx, inputScalars.dy, inputScalars.dz };
		d_N = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz };
		bmax = { static_cast<float>(inputScalars.Nx) * inputScalars.dx + inputScalars.bx,
			static_cast<float>(inputScalars.Ny) * inputScalars.dy + inputScalars.by,
			static_cast<float>(inputScalars.Nz) * inputScalars.dz + inputScalars.bz };
		dPitch = { w_vec.dPitchX, w_vec.dPitchY };

		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type == 14 || inputScalars.projector_type == 41) {
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type == 41) {
				kernelFP.setArg(kernelIndFP++, d_N);
				kernelFP.setArg(kernelIndFP++, b);
				kernelFP.setArg(kernelIndFP++, w_vec.size_x);
				kernelFP.setArg(kernelIndFP++, w_vec.size_y);
				kernelFP.setArg(kernelIndFP++, dPitch);
			}

			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type == 14) {
				kernelBP.setArg(kernelIndBP++, d_N);
				kernelBP.setArg(kernelIndBP++, b);
				kernelBP.setArg(kernelIndBP++, w_vec.size_x);
				kernelBP.setArg(kernelIndBP++, w_vec.size_y);
				kernelBP.setArg(kernelIndBP++, dPitch);
			}
		}
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14 || inputScalars.projector_type == 41) {
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41) {
				kernelFP.setArg(kernelIndFP++, bmax);
				kernelFP.setArg(kernelIndFP++, inputScalars.dL);
				kernelFP.setArg(kernelIndFP++, inputScalars.d_Scale);
				status = kernelFP.setArg(kernelIndFP++, inputScalars.global_factor);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}

			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14) {
				kernelBP.setArg(kernelIndBP++, d);
				kernelBP.setArg(kernelIndBP++, w_vec.kerroin4);
			}
		}
		else if (inputScalars.projector_type == 5) {
			kernelFP.setArg(kernelIndFP++, d);
			kernelFP.setArg(kernelIndFP++, inputScalars.d_Scale);
			kernelFP.setArg(kernelIndFP++, inputScalars.dSize);

			kernelBP.setArg(kernelIndBP++, d);
			kernelBP.setArg(kernelIndBP++, inputScalars.d_Scale);
			kernelBP.setArg(kernelIndBP++, inputScalars.dSizeBP);
		}
		if (inputScalars.projector_type == 14 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22 || inputScalars.projector_type == 33) {
			// Set the kernelFP parameters that do not change
			status = kernelFP.setArg(kernelIndFP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			kernelFP.setArg(kernelIndFP++, inputScalars.epps);
			kernelFP.setArg(kernelIndFP++, inputScalars.im_dim);
			kernelFP.setArg(kernelIndFP++, d_N);
			kernelFP.setArg(kernelIndFP++, d);
			kernelFP.setArg(kernelIndFP++, b);
			kernelFP.setArg(kernelIndFP++, bmax);
			kernelFP.setArg(kernelIndFP++, w_vec.size_x);
			kernelFP.setArg(kernelIndFP++, inputScalars.det_per_ring);
			kernelFP.setArg(kernelIndFP++, inputScalars.Nxy);
			kernelFP.setArg(kernelIndFP++, inputScalars.sigma_x);
			kernelFP.setArg(kernelIndFP++, dPitch);
			if (inputScalars.projector_type == 22 || inputScalars.projector_type == 33) {
				kernelFP.setArg(kernelIndFP++, inputScalars.tube_width);
				kernelFP.setArg(kernelIndFP++, inputScalars.bmin);
				kernelFP.setArg(kernelIndFP++, inputScalars.bmax);
				kernelFP.setArg(kernelIndFP++, inputScalars.Vmax);
			}
		}
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			status = kernelBP.setArg(kernelIndBP++, inputScalars.global_factor);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			kernelBP.setArg(kernelIndBP++, inputScalars.epps);
			kernelBP.setArg(kernelIndBP++, inputScalars.im_dim);
			kernelBP.setArg(kernelIndBP++, d_N);
			kernelBP.setArg(kernelIndBP++, d);
			kernelBP.setArg(kernelIndBP++, b);
			kernelBP.setArg(kernelIndBP++, bmax);
			kernelBP.setArg(kernelIndBP++, w_vec.size_x);
			kernelBP.setArg(kernelIndBP++, inputScalars.det_per_ring);
			kernelBP.setArg(kernelIndBP++, inputScalars.Nxy);
			kernelBP.setArg(kernelIndBP++, inputScalars.sigma_x);
			kernelBP.setArg(kernelIndBP++, dPitch);
			if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
				kernelBP.setArg(kernelIndBP++, inputScalars.tube_width);
				kernelBP.setArg(kernelIndBP++, inputScalars.bmin);
				kernelBP.setArg(kernelIndBP++, inputScalars.bmax);
				kernelBP.setArg(kernelIndBP++, inputScalars.Vmax);
			}
		}
		if (inputScalars.maskFP || inputScalars.maskBP) {
			if (inputScalars.maskFP) {
				status = kernelFP.setArg(kernelIndFP++, d_maskFP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.maskBP) {
				if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type == 14)
					status = kernelBP.setArg(kernelIndBP++, d_maskBP);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		if (inputScalars.projector_type == 14 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			if (inputScalars.TOF)
				kernelFP.setArg(kernelIndFP++, d_TOFCenter);
			if (inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
				kernelFP.setArg(kernelIndFP++, d_xcenter);
				kernelFP.setArg(kernelIndFP++, d_ycenter);
				kernelFP.setArg(kernelIndFP++, d_zcenter);
				kernelFP.setArg(kernelIndFP++, d_V);
			}
			status = kernelFP.setArg(kernelIndFP++, d_reko_type);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
		}
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			if (inputScalars.TOF)
				kernelBP.setArg(kernelIndBP++, d_TOFCenter);
			if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
				kernelBP.setArg(kernelIndBP++, d_xcenter);
				kernelBP.setArg(kernelIndBP++, d_ycenter);
				kernelBP.setArg(kernelIndBP++, d_zcenter);
				kernelBP.setArg(kernelIndBP++, d_V);
			}
			kernelBP.setArg(kernelIndBP++, d_reko_type);
		}
		if (DEBUG) {
			mexPrintf("kernelIndFP = %u\n", kernelIndFP);
			mexPrintf("kernelIndBP = %u\n", kernelIndBP);
			mexEvalString("pause(.0001);");
			//mexEvalString("pause(2);");
		}
		return 0;
	}

	inline int loadDynamicData(scalarStruct& inputScalars, const std::vector<int64_t>& length, const float* Sino, const float* randomsData,
		const float* scat, const int64_t* pituus) {

		cl_int status = CL_SUCCESS;
		for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
			if (inputScalars.TOF) {
				if (!inputScalars.loadTOF && kk == 0) {
					d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * inputScalars.nBins, NULL, &status);
					for (int64_t to = 0LL; to < inputScalars.nBins; to++)
						status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sino[pituus[kk] + inputScalars.koko * to]);
				}
				else if (inputScalars.loadTOF) {
					for (int64_t to = 0LL; to < inputScalars.nBins; to++)
						status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sino[pituus[kk] + inputScalars.koko * to]);
				}
			}
			else
				status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[kk], &Sino[pituus[kk]]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.randoms_correction) {
				status = af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_TRUE, 0, sizeof(float) * length[kk], &randomsData[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.scatter == 1u) {
				status = af_queue.enqueueWriteBuffer(d_scat[kk], CL_TRUE, 0, sizeof(float) * length[kk], &scat[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		return 0;
	}

	inline int setDynamicKernelData(scalarStruct& inputScalars, Weighting& w_vec) {
		cl_int status = CL_SUCCESS;
		if (inputScalars.projector_type == 14 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u || inputScalars.projector_type == 41) {
			if (inputScalars.projector_type == 14 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
				status = kernelFP.setArg(kernelIndFP++, w_vec.epsilon_mramla);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelFP.setArg(kernelIndFP++, w_vec.size_y);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.attenuation_correction && !inputScalars.CT) {
				status = kernelFP.setArg(kernelIndFP++, d_atten);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
		}
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
				kernelBP.setArg(kernelIndBP++, w_vec.epsilon_mramla);
				status = kernelBP.setArg(kernelIndBP++, w_vec.size_y);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				if (inputScalars.attenuation_correction && !inputScalars.CT)
					kernelBP.setArg(kernelIndBP++, d_atten);
		}
		return 0;
	}

	inline int loadTOFData(scalarStruct& inputScalars, const float* Sino, const int64_t length, const int64_t pituus) {
		cl_int status = CL_SUCCESS;
		d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length * inputScalars.nBins, NULL, &status);
		for (int64_t to = 0LL; to < inputScalars.nBins; to++)
			status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length * to, sizeof(float) * length, &Sino[pituus + inputScalars.koko * to]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		return 0;
	}

	inline void transferSensitivityImage(af::array& apuSum) {
		af::sync();
		d_Summ = cl::Buffer(*apuSum.device<cl_mem>(), true);
	}

	inline int update_opencl_inputs(AF_im_vectors& vec, const scalarStruct& inputScalars) {

		const cl::detail::size_t_array origin = { 0, 0, 0 };
		cl::detail::size_t_array region = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz * inputScalars.nRekos };
		if (inputScalars.projector_type == 5) {
			af::array intIm;
			if (inputScalars.meanFP) {
				af::array im = af::reorder(af::moddims(vec.im_os, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, inputScalars.nRekos), 1, 2, 0, 3);
				vec.meanFP = af::constant(0.f, inputScalars.Nx + inputScalars.Ny);
				vec.meanFP(af::seq(0, inputScalars.Nx - 1)) = af::flat(af::mean(af::mean(im, 0), 1));
				im -= af::tile(vec.meanFP(af::seq(0, inputScalars.Nx - 1)), im.dims(0), im.dims(1), 1);
				intIm = af::constant(0.f, im.dims(0) + 1, im.dims(1) + 1, im.dims(2));
				intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
				if (DEBUG) {
					mexPrintf("intIm.dims(0) = %u\n", intIm.dims(0));
					mexPrintf("intIm.dims(1) = %u\n", intIm.dims(1));
					mexPrintf("intIm.dims(2) = %u\n", intIm.dims(2));
					mexEvalString("pause(.0001);");
				}
				region = { static_cast<cl_ulong>(intIm.dims(0)), static_cast<cl_ulong>(intIm.dims(1)), static_cast<cl_ulong>(intIm.dims(2)) };
				af::sync();
				cl_int status = af_queue.enqueueCopyBufferToImage(cl::Buffer(*(af::flat(intIm)).device<cl_mem>(), true), vec_opencl.d_image_os_int, 0, origin, region);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Integral image 1 copy failed\n");
					mexEvalString("pause(.0001);");
					return -1;
				}
				im = af::reorder(af::moddims(vec.im_os, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz, inputScalars.nRekos), 0, 2, 1, 3);
				vec.meanFP(af::seq(inputScalars.Nx, inputScalars.Nx + inputScalars.Ny)) = af::flat(af::mean(af::mean(im, 0), 1));
				im -= af::tile(vec.meanFP(af::seq(inputScalars.Nx, inputScalars.Nx + inputScalars.Ny)), im.dims(0), im.dims(1), 1);
				af_queue.finish();
				intIm.unlock();
				intIm = af::constant(0.f, im.dims(0) + 1, im.dims(1) + 1, im.dims(2));
				intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
				intIm = af::sat(im);
				region = { static_cast<cl_ulong>(intIm.dims(0)), static_cast<cl_ulong>(intIm.dims(1)), static_cast<cl_ulong>(intIm.dims(2)) };
				af::sync();
				status = af_queue.enqueueCopyBufferToImage(cl::Buffer(*(af::flat(intIm)).device<cl_mem>(), true), vec_opencl.d_image_os, 0, origin, region);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Integral image 2 copy failed\n");
					mexEvalString("pause(.0001);");
					return -1;
				}
				af_queue.finish();
				intIm.unlock();
			}
			else {
				intIm = af::constant(0.f, inputScalars.Ny + 1, inputScalars.Nz + 1, inputScalars.Nx);
				intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os, inputScalars.Nx, 
					inputScalars.Ny, inputScalars.Nz * inputScalars.nRekos), 1, 2, 0));
				if (DEBUG) {
					mexPrintf("intIm.dims(0) = %u\n", intIm.dims(0));
					mexPrintf("intIm.dims(1) = %u\n", intIm.dims(1));
					mexPrintf("intIm.dims(2) = %u\n", intIm.dims(2));
					mexEvalString("pause(.0001);");
				}
				region = { inputScalars.Ny + 1, inputScalars.Nz + 1, inputScalars.Nx };
				intIm = af::flat(intIm);
				//vec.im_os(af::seq(0, (inputScalars.Ny + 1) * (inputScalars.Nz + 1) * 2)) = intIm(af::seq(0, (inputScalars.Ny + 1) * (inputScalars.Nz + 1) * 2));
				af::sync();
				cl::Buffer d_apu = cl::Buffer(*intIm.device<cl_mem>(), true);
				af_queue.finish();
				cl_int status = af_queue.enqueueCopyBufferToImage(d_apu, vec_opencl.d_image_os_int, 0, origin, region);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Integral image 1 copy failed\n");
					mexEvalString("pause(.0001);");
					return -1;
				}
				af_queue.finish();
				intIm.unlock();
				intIm = af::constant(0.f, inputScalars.Nx + 1, inputScalars.Nz + 1, inputScalars.Ny);
				intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os, inputScalars.Nx, 
					inputScalars.Ny, inputScalars.Nz, inputScalars.nRekos), 0, 2, 1, 3));
				intIm = af::flat(intIm);
				region = { inputScalars.Nx + 1, inputScalars.Nz + 1, inputScalars.Ny };
				af::sync();
				cl::Buffer d_apu2 = cl::Buffer(*intIm.device<cl_mem>(), true);
				af_queue.finish();
				status = af_queue.enqueueCopyBufferToImage(d_apu2, vec_opencl.d_image_os, 0, origin, region);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Integral image 2 copy failed\n");
					mexEvalString("pause(.0001);");
					return -1;
				}
				af_queue.finish();
				intIm.unlock();
			}
			//vec.im_os = af::reorder(af::moddims(vec.im_os, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz * n_rekos), 1, 2, 0);
		}
		else {
			af::sync();
			//if (inputScalars.use_psf)
			//	vec_opencl.d_im_os = cl::Buffer(*vec.im_os_blurred.device<cl_mem>(), true);
			//else
			//	vec_opencl.d_im_os = cl::Buffer(*vec.im_os.device<cl_mem>(), true);
			cl_mem* im;
			if (inputScalars.use_psf)
				im = vec.im_os_blurred.device<cl_mem>();
			else
				im = vec.im_os.device<cl_mem>();
			cl::Buffer d_im_os = cl::Buffer(*im, true);
			af_queue.finish();
			cl_int status = af_queue.enqueueCopyBufferToImage(d_im_os, vec_opencl.d_image_os, 0, origin, region);
			af_queue.finish();
			if (inputScalars.use_psf)
				vec.im_os_blurred.unlock();
			else
				vec.im_os.unlock();
			delete im;
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Image copy failed\n");
				mexEvalString("pause(.0001);");
				return -1;
			}
		}
		//af::sync();
		//af_queue.finish();
		//vec_opencl.d_rhs_os = cl::Buffer(*vec.rhs_os.device<cl_mem>(), true);
		af_queue.finish();
		return 0;
	}

	inline int forwardProjection(scalarStruct& inputScalars, Weighting& w_vec, af::array& outputFP, const uint32_t osa_iter,
		const std::vector<int64_t>& length, const uint64_t st, const uint64_t m_size) {
		kernelIndFPSubIter = kernelIndFP;
		cl_int status = CL_SUCCESS;
		d_output = cl::Buffer(*outputFP.device<cl_mem>(), true);
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { w_vec.size_x + erotus[0], w_vec.size_y + erotus[1], static_cast<size_t>(length[osa_iter]) };
		else
			global = { length[osa_iter] + erotus[0], 1, 1 };

		if (DEBUG) {
			mexPrintf("global[0] = %u\n", global[0]);
			mexPrintf("local[0] = %u\n", local[0]);
			mexPrintf("local[1] = %u\n", local[1]);
			mexPrintf("global[1] = %u\n", global[1]);
			mexPrintf("global[2] = %u\n", global[2]);
			mexPrintf("erotus[0] = %u\n", erotus[0]);
			mexPrintf("erotus[1] = %u\n", erotus[1]);
			mexPrintf("global.dimensions() = %u\n", global.dimensions());
			mexPrintf("local.dimensions() = %u\n", local.dimensions());
			mexPrintf("kernelIndFPSubIter = %u\n", kernelIndFPSubIter);
			mexPrintf("kernelIndFP = %u\n", kernelIndFP);
			mexPrintf("m_size = %u\n", m_size);
			mexPrintf("size_x = %u\n", w_vec.size_x);
			mexPrintf("size_y = %u\n", w_vec.size_y);
			mexPrintf("length[osa_iter] = %u\n", length[osa_iter]);
			mexPrintf("listmode = %u\n", inputScalars.listmode);
			mexPrintf("maskBP = %u\n", inputScalars.maskBP);
			mexPrintf("no_norm = %u\n", no_norm);
			mexPrintf("NVOXELS = %u\n", NVOXELS);
			mexEvalString("pause(.0001);");
			//mexEvalString("pause(2);");
		}

		//af::sync();
		af_queue.finish();
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 41) {
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
			if (inputScalars.listmode == 0 && !inputScalars.CT)
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
			status = kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
			if (inputScalars.subsetType < 8 && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				kernelFP.setArg(kernelIndFPSubIter++, d_xyindex[osa_iter]);
				kernelFP.setArg(kernelIndFPSubIter++, d_zindex[osa_iter]);
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
		}
		else if (inputScalars.projector_type == 5) {
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

			}
			kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
		}
		else if (inputScalars.projector_type == 11 || inputScalars.projector_type == 14 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0) {
				status = kernelFP.setArg(kernelIndFPSubIter++, length[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			if (inputScalars.listmode == 0 && !inputScalars.CT)
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[0]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_x[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.CT || inputScalars.PET) && inputScalars.listmode == 0)
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[osa_iter]);
			else
				status = kernelFP.setArg(kernelIndFPSubIter++, d_z[0]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.normalization_correction)
				kernelFP.setArg(kernelIndFPSubIter++, d_norm[osa_iter]);
			if (inputScalars.scatter)
				kernelFP.setArg(kernelIndFPSubIter++, d_scat[osa_iter]);
			status = kernelFP.setArg(kernelIndFPSubIter++, d_Summ);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_uchar>(1));
			if (inputScalars.precompute)
				kernelFP.setArg(kernelIndFPSubIter++, d_lor[osa_iter]);
			if (inputScalars.subsetType < 8 && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				kernelFP.setArg(kernelIndFPSubIter++, d_xyindex[osa_iter]);
				kernelFP.setArg(kernelIndFPSubIter++, d_zindex[osa_iter]);
			}
			if (inputScalars.raw)
				kernelFP.setArg(kernelIndFPSubIter++, d_L[osa_iter]);
			if (inputScalars.TOF && !inputScalars.loadTOF)
				kernelFP.setArg(kernelIndFPSubIter++, d_Sino[0]);
			else
				kernelFP.setArg(kernelIndFPSubIter++, d_Sino[osa_iter]);
			//if (inputScalars.randoms_correction)
			//	kernelFP.setArg(kernelIndFPSubIter++, d_sc_ra[osa_iter]);
			kernelFP.setArg(kernelIndFPSubIter++, vec_opencl.d_image_os);
			//if (projector_type == 7) {
			//	kernelFP.setArg(kernelIndFPSubIter++, d_indices);
			//	kernelFP.setArg(kernelIndFPSubIter++, d_rowInd);
			//	kernelFP.setArg(kernelIndFPSubIter++, d_values);
			//	kernelFP.setArg(kernelIndFPSubIter++, dec);
			//}
			kernelFP.setArg(kernelIndFPSubIter++, d_output);
			kernelFP.setArg(kernelIndFPSubIter++, no_norm);
			kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_ulong>(m_size));
			kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_ulong>(st));
		}
		status = af_queue.enqueueNDRangeKernel(kernelFP, cl::NDRange(), global, local, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG) {
			mexPrintf("Forward projection kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		return 0;
	}

	inline int backwardProjection(scalarStruct& inputScalars, Weighting& w_vec, af::array& outputFP, const uint32_t osa_iter,
		const std::vector<int64_t>& length, const uint64_t st, const uint64_t m_size, af::array& meanBP) {
		cl_int status = CL_SUCCESS;
		kernelIndBPSubIter = kernelIndBP;

		status = af_queue.finish();
		if (inputScalars.projector_type >= 4)
			d_output = cl::Buffer(*outputFP.device<cl_mem>(), true);
		if (inputScalars.projector_type < 4 || inputScalars.projector_type == 41 || inputScalars.projector_type == 11 || inputScalars.projector_type == 22u || inputScalars.projector_type == 33u) {
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				global = { w_vec.size_x + erotus[0], w_vec.size_y + erotus[1], static_cast<size_t>(length[osa_iter]) };
			else {
				erotus[0] = length[osa_iter] % local_size[0];

				if (erotus[0] > 0)
					erotus[0] = (local_size[0] - erotus[0]);
				global = { length[osa_iter] + erotus[0], 1, 1 };
			}

			if (DEBUG) {
				mexPrintf("global[0] = %u\n", global[0]);
				mexPrintf("local[0] = %u\n", local[0]);
				mexPrintf("local[1] = %u\n", local[1]);
				mexPrintf("global[1] = %u\n", global[1]);
				mexPrintf("global[2] = %u\n", global[2]);
				mexPrintf("erotus[0] = %u\n", erotus[0]);
				mexPrintf("erotus[1] = %u\n", erotus[1]);
				mexPrintf("global.dimensions() = %u\n", global.dimensions());
				mexPrintf("local.dimensions() = %u\n", local.dimensions());
				mexPrintf("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
				mexPrintf("m_size = %u\n", m_size);
				mexPrintf("size_x = %u\n", w_vec.size_x);
				mexPrintf("size_y = %u\n", w_vec.size_y);
				mexPrintf("length[osa_iter] = %u\n", length[osa_iter]);
				mexPrintf("st = %u\n", st);
				mexPrintf("listmode = %u\n", inputScalars.listmode);
				mexPrintf("im_dim = %u\n", inputScalars.im_dim);
				mexPrintf("no_norm = %u\n", no_norm);
				mexEvalString("pause(.0001);");
				//mexEvalString("pause(2);");
			}

			// Set kernelBP arguments
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0)
				status = kernelBP.setArg(kernelIndBPSubIter++, length[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.listmode == 0 && !inputScalars.CT)
				status = kernelBP.setArg(kernelIndBPSubIter++, d_x[0]);
			else
				status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if ((inputScalars.CT || inputScalars.PET) && inputScalars.listmode == 0)
				status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
			else
				status = kernelBP.setArg(kernelIndBPSubIter++, d_z[0]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			if (inputScalars.normalization_correction)
				status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[osa_iter]);
			if (inputScalars.scatter)
				status = kernelBP.setArg(kernelIndBPSubIter++, d_scat[osa_iter]);
			status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return -1;
			}
			kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_uchar>(2));
			if (inputScalars.precompute)
				kernelBP.setArg(kernelIndBPSubIter++, d_lor[osa_iter]);
			if (inputScalars.subsetType < 8 && inputScalars.subsets > 1 && inputScalars.listmode == 0) {
				kernelBP.setArg(kernelIndBPSubIter++, d_xyindex[osa_iter]);
				kernelBP.setArg(kernelIndBPSubIter++, d_zindex[osa_iter]);
			}
			if (inputScalars.raw)
				kernelBP.setArg(kernelIndBPSubIter++, d_L[osa_iter]);
			if (inputScalars.TOF && !inputScalars.loadTOF)
				kernelBP.setArg(kernelIndBPSubIter++, d_Sino[0]);
			else
				kernelBP.setArg(kernelIndBPSubIter++, d_Sino[osa_iter]);
			if (inputScalars.randoms_correction && inputScalars.projector_type < 4)
				kernelBP.setArg(kernelIndBPSubIter++, d_sc_ra[osa_iter]);
			if (inputScalars.projector_type < 4)
				kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_image_os);
			else
				kernelBP.setArg(kernelIndBPSubIter++, d_output);
			//if (projector_type == 7) {
			//	kernelBP.setArg(kernelIndBPSubIter++, d_indices);
			//	kernelBP.setArg(kernelIndBPSubIter++, d_rowInd);
			//	kernelBP.setArg(kernelIndBPSubIter++, d_values);
			//	kernelBP.setArg(kernelIndBPSubIter++, dec);
			//}
			kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os);
			kernelBP.setArg(kernelIndBPSubIter++, no_norm);
			kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_ulong>(m_size));
			kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_ulong>(st));
		}
		else {

			if (inputScalars.meanBP && inputScalars.projector_type == 5)
				d_meanBP = cl::Buffer(*meanBP.device<cl_mem>(), true); 

			cl::size_type imX = w_vec.size_x;
			cl::size_type imY = w_vec.size_y;
			cl::size_type imZ = length[osa_iter];
			if (inputScalars.projector_type == 5) {
				imX++;
				imY++;
			}
			cl::detail::size_t_array region = { imX, imY, imZ };
			d_inputImage = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Image creation failed\n");
				mexEvalString("pause(.0001);");
				return -1;
			}

			status = af_queue.enqueueCopyBufferToImage(d_output, d_inputImage, 0, origin, region);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Image copy failed\n");
				mexEvalString("pause(.0001);");
				return -1;
			}
			status = af_queue.finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrintf("Queue finish failed after image copy\n");
				mexEvalString("pause(.0001);");
			}
			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14)
				global = { inputScalars.Nx + erotusBP[0], inputScalars.Ny + erotusBP[1], (inputScalars.Nz + NVOXELS - 1) / NVOXELS };
			else
				global = { inputScalars.Nx + erotusBP[0], inputScalars.Ny + erotusBP[1], inputScalars.Nz };

			if (DEBUG) {
				mexPrintf("global[0] = %u\n", global[0]);
				mexPrintf("local[0] = %u\n", local[0]);
				mexPrintf("local[1] = %u\n", local[1]);
				mexPrintf("global[1] = %u\n", global[1]);
				mexPrintf("global[2] = %u\n", global[2]);
				mexPrintf("erotusBP[0] = %u\n", erotusBP[0]);
				mexPrintf("erotusBP[1] = %u\n", erotusBP[1]);
				mexPrintf("global.dimensions() = %u\n", global.dimensions());
				mexPrintf("local.dimensions() = %u\n", local.dimensions());
				mexPrintf("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
				mexPrintf("m_size = %u\n", m_size);
				mexPrintf("size_x = %u\n", w_vec.size_x);
				mexPrintf("size_y = %u\n", w_vec.size_y);
				mexPrintf("length[osa_iter] = %u\n", length[osa_iter]);
				mexPrintf("st = %u\n", st);
				mexPrintf("listmode = %u\n", inputScalars.listmode);
				mexPrintf("im_dim = %u\n", inputScalars.im_dim);
				mexPrintf("no_norm = %u\n", no_norm);
				mexEvalString("pause(.0001);");
				//mexEvalString("pause(2);");
			}

			if (inputScalars.projector_type == 4 || inputScalars.projector_type == 14) {
				status = kernelBP.setArg(kernelIndBPSubIter++, d_inputImage);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, d_z[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
			}
			else {
				status = kernelBP.setArg(kernelIndBPSubIter++, d_x[osa_iter]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
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
				status = kernelBP.setArg(kernelIndBPSubIter++, vec_opencl.d_rhs_os);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return -1;
				}
				status = kernelBP.setArg(kernelIndBPSubIter++, d_Summ);
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
			kernelBP.setArg(kernelIndBPSubIter++, no_norm);
			kernelBP.setArg(kernelIndBPSubIter++, static_cast<cl_long>(length[osa_iter]));
		}
		status = af_queue.enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, local, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		else if (DEBUG) {
			mexPrintf("Backprojection kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		return 0;
	}

	inline int64_t getGlobalMem() {
		cl_int status = CL_SUCCESS;
		int64_t mem;
		cl_ulong mem_loc;
		mem = af_device_id.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}

		mem_loc = af_device_id.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return -1;
		}
		if (DEBUG) {
			mexPrintf("mem_loc = %u\n", mem_loc);
			mexEvalString("pause(.0001);");
		}
		return mem;
	}

	inline void computeMRP(af::array& padd, af::array& grad, const scalarStruct& inputScalars) {
		uint32_t kernelIndMed = 0U;
		cl::NDRange global_size(padd.dims(0), padd.dims(1), padd.dims(2));
		//cl::NDRange local_size(64U);
		cl::Buffer d_grad = cl::Buffer(*grad.device<cl_mem>(), true);
		cl::Buffer d_padd = cl::Buffer(*padd.device<cl_mem>(), true);
		if (DEBUG) {
			mexPrintf("padd = %f\n", af::sum<float>(padd));
		}
		af_queue.finish();
		kernelMed.setArg(kernelIndMed++, d_padd);
		kernelMed.setArg(kernelIndMed++, d_grad);
		kernelMed.setArg(kernelIndMed++, inputScalars.Nx);
		kernelMed.setArg(kernelIndMed++, inputScalars.Ny);
		kernelMed.setArg(kernelIndMed++, inputScalars.Nz);
		cl_int status = af_queue.enqueueNDRangeKernel(kernelMed, cl::NullRange, global_size, cl::NullRange);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to launch the Median filter kernel\n");
			mexEvalString("pause(.0001);");
		}
		else if (DEBUG) {
			mexPrintf("Median kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
		status = af_queue.finish();
	}

	inline void computeNLM(af::array& grad, const af::array im, const scalarStruct& inputScalars, Weighting& w_vec, const int type) {
		cl_int status = CL_SUCCESS;
		cl::detail::size_t_array region = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz * inputScalars.nRekos };
		const cl_int3 searchWindow = { static_cast<cl_int>(w_vec.Ndx) , static_cast<cl_int>(w_vec.Ndy) , static_cast<cl_int>(w_vec.Ndz) };
		const cl_int3 patchWindow = { static_cast<cl_int>(w_vec.Nlx) , static_cast<cl_int>(w_vec.Nly) , static_cast<cl_int>(w_vec.Nlz) };
		const cl_uint3 N = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz };
		size_t kernelIndNLM = 0ULL;
		cl::NDRange global = { inputScalars.Nx + erotusBP[0], inputScalars.Ny + erotusBP[1], inputScalars.Nz };
		cl::Buffer d_W = cl::Buffer(*grad.device<cl_mem>(), true);
		cl::Image3D d_input = cl::Image3D(af_context, CL_MEM_READ_ONLY, format, inputScalars.Nx, inputScalars.Ny, inputScalars.Nz * inputScalars.nRekos, 0, 0, NULL, &status);
		status = af_queue.enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), d_input, 0, origin, region);
		cl::Buffer d_gaussianNLM = cl::Buffer(*w_vec.gaussianNLM.device<cl_mem>(), true);
		kernelNLM.setArg(kernelIndNLM++, d_W);
		kernelNLM.setArg(kernelIndNLM++, d_input);
		kernelNLM.setArg(kernelIndNLM++, d_gaussianNLM);
		kernelNLM.setArg(kernelIndNLM++, searchWindow);
		kernelNLM.setArg(kernelIndNLM++, patchWindow);
		kernelNLM.setArg(kernelIndNLM++, N);
		kernelNLM.setArg(kernelIndNLM++, w_vec.h2);
		kernelNLM.setArg(kernelIndNLM++, inputScalars.epps);
		kernelNLM.setArg(kernelIndNLM++, type);
		if (w_vec.NLM_anatomical)
			kernelNLM.setArg(kernelIndNLM++, d_uref);
		// Compute the kernel
		status = (af_queue).enqueueNDRangeKernel(kernelNLM, cl::NullRange, global, local);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to launch the NLM kernel\n");
			mexEvalString("pause(.0001);");
		}

		status = (af_queue).finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Queue finish failed after kernel\n");
			mexEvalString("pause(.0001);");
		}
	}

	inline void transferRHS(AF_im_vectors& vec) {
		af::sync();
		vec_opencl.d_rhs_os = cl::Buffer(*vec.rhs_os.device<cl_mem>(), true);
		af_queue.finish();
	}

};

