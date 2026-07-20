/*******************************************************************************************************************************************
* Class object for forward and backward projections. Combined OpenCL and CUDA version.
*
* Select the backend with the CUDA preprocessor definition: define CUDA/HIP for the CUDA/HIP backend, otherwise the OpenCL backend is used.
*
* Copyright (C) 2022-2026 Ville-Veikko Wettenhovi, Niilo Saarlemo
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

// ======== OpenCL / CUDA compatibility macros (generated) ========
#if defined(CUDA) || defined(HIP)
#define STATUS_t CUresult
#define SUCCESS_VALUE CUDA_SUCCESS
#define INT3_t int3
#define UINT2_t uint2
#define DEVBUF_t CUdeviceptr
#else
#define STATUS_t cl_int
#define SUCCESS_VALUE CL_SUCCESS
#define INT3_t cl_int3
#define UINT2_t cl_uint2
#define DEVBUF_t cl::Buffer
#endif
#if defined(CUDA) || defined(HIP)
#define GET_KERNEL(KVAR, PROG, NAME) status = cuModuleGetFunction(&KVAR, PROG, NAME)
#define KCHECK(MSG) CUDA_CHECK(status, MSG, status)
#else
#define GET_KERNEL(KVAR, PROG, NAME) KVAR = cl::Kernel(PROG, NAME, &status)
#define KCHECK(MSG) OCL_CHECK(status, MSG, -1)
#endif
#define CREATE_KERNEL(KVAR, PROG, NAME, MSG) GET_KERNEL(KVAR, PROG, NAME); KCHECK(MSG)
// Unified status check (uses the backend success code); replaces paired CUDA_CHECK / OCL_CHECK.
#define CHECK(STATUS, MSG, RETURN) do { if ((STATUS) != SUCCESS_VALUE) { getErrorString(STATUS); mexPrint(MSG); return RETURN; } } while(0)
// Allocate a device buffer. FLAGS is the OpenCL cl_mem_flags (ignored by CUDA, which has no
// context/flags). Host pointer is always NULL here (data is written separately).
#if defined(CUDA) || defined(HIP)
#define ALLOC_BUFFER(BUF, FLAGS, SIZE) status = cuMemAlloc(&BUF, SIZE)
#else
#define ALLOC_BUFFER(BUF, FLAGS, SIZE) BUF = cl::Buffer(CLContext, FLAGS, SIZE, NULL, &status)
#endif
// Upload SIZE bytes from host SRC into device buffer BUF (non-blocking on OpenCL).
#if defined(CUDA) || defined(HIP)
#define WRITE_BUFFER(BUF, SIZE, SRC) status = cuMemcpyHtoD(BUF, SRC, SIZE)
#else
#define WRITE_BUFFER(BUF, SIZE, SRC) status = CLCommandQueue[0].enqueueWriteBuffer(BUF, CL_FALSE, 0, SIZE, SRC)
#endif
// Append one kernel argument VAR. CUDA pushes its address into the argument vector VEC; OpenCL
// sets it on KERNEL at the running index IDX, reporting any error inline via getErrorString.
#if defined(CUDA) || defined(HIP)
#define KARG(VEC, KERNEL, IDX, VAR) VEC.emplace_back(reinterpret_cast<void*>(&VAR))
#else
#define KARG(VEC, KERNEL, IDX, VAR) getErrorString(KERNEL.setArg(IDX++, VAR))
#endif
#if defined(CUDA) || defined(HIP)
#define ADD_OPT(OPTS, FLAG) OPTS.push_back(FLAG)
#else
#define ADD_OPT(OPTS, FLAG) OPTS += " " FLAG
#endif
// Build option selecting the forward projection. hipRTC force-includes hiprtc_runtime.h, which
// uses FP as a type name, so a command-line -DFP breaks that header. HIP passes -DOMEGA_FP
// instead; general_opencl_functions.h maps it back to FP after the hipRTC prelude.
#if defined(HIP)
#define FP_FLAG "-DOMEGA_FP"
#else
#define FP_FLAG "-DFP"
#endif
// Append an integer-valued build option "FLAG=VALUE". CUDA stores the formatted std::string
// directly in the (std::string) options vector; OpenCL appends it to the options string.
#if defined(CUDA) || defined(HIP)
#define ADD_OPT_INT(OPTS, FLAG, VALUE) OPTS.push_back(FLAG "=" + std::to_string(VALUE))
#else
#define ADD_OPT_INT(OPTS, FLAG, VALUE) OPTS += (" " FLAG "=" + std::to_string(VALUE))
#endif
#if defined(CUDA) || defined(HIP)
#define BACKEND_STR "CUDA"
#else
#define BACKEND_STR "OpenCL"
#endif
// ================================================================

/// <summary>
/// Class object for forward and backward projections. Combined OpenCL and CUDA version
/// </summary>
class ProjectorClass {
	//private:
		// Local size
	size_t local_size[3];
	size_t local_sizePrior[3];
#if defined(CUDA) || defined(HIP)
	// Contents of general_opencl_functions.h. For CUDA/HIP this is handed to NVRTC/HIPRTC as an actual
	// header, so the kernel sources #include it rather than having it concatenated in front of them (as is
	// still done for OpenCL).
	std::string nvrtcKernelHeader;
	// Contents of opencl_functions_orth3D.h, supplied as a distinct header. Non-empty (and passed to
	// NVRTC/HIPRTC) only for the forward/backward projector programs that use it; empty for the auxiliary
	// kernels so it is not included in that build.
	std::string nvrtcOrthHeader;
#endif // END CUDA
	// Kernel input indices
#if defined(CUDA) || defined(HIP)
	unsigned int kernelInd_MRAMLA = 0;
	unsigned int kernelIndFP = 0;
	unsigned int kernelIndBP = 0;
	unsigned int kernelIndFPSubIter = 0;
	unsigned int kernelIndBPSubIter = 0;
	// Crystal pitch
	float2 dPitch;
#else
	cl_uint kernelInd_MRAMLA = 0;
	cl_uint kernelIndFP = 0;
	cl_uint kernelIndBP = 0;
	cl_uint kernelIndFPSubIter = 0;
	cl_uint kernelIndBPSubIter = 0;
	cl_uint kernelIndSens = 0;
#endif // END CUDA
	// Total FOV boundary including multi-resolution
#if defined(CUDA) || defined(HIP)
	float3 totalFOVmin, totalFOVmax;
#else
	cl_float3 totalFOVmin, totalFOVmax;
	// Crystal pitch
	cl_float2 dPitch;
#endif // END CUDA
	// Image dimensions
	INT3_t d_NOrig, d_NPrior;
	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusPrior[3];
	size_t erotusPriorEFOV[3];
	size_t erotusSens[3];
#if defined(CUDA) || defined(HIP)
	size_t kSize = 0ULL;
#endif // END CUDA
	// Local and global sizes
#if defined(CUDA) || defined(HIP)
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
		bool geom5 = false;
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
		int g5Steps = 0;
		int tSteps = 1;
		int attenSize = 0;
	};
	CUDAMemAlloc memAlloc;
	bool useBuffers = true;
#else
	cl::NDRange local, global, localPrior, globalPrior, globalPriorEFOV;
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
	template <typename K, typename T>
	inline K make_vec3(T a, T b, T c) {
		K apu;
		apu.x = a;
		apu.y = b;
		apu.z = c;
		return apu;
#else
	bool constantBuffer = false;

	// Get the OpenCL context for the current platform
	cl_int clGetPlatformsContext(const uint32_t platform, cl::Context & context, std::vector<cl::CommandQueue>&commandQueues, const std::vector<uint32_t>&usedDevices, std::vector<cl::Device>&devices) {
		cl_int status = CL_SUCCESS;

		// Get the number of platforms 
		std::vector<cl::Platform> platforms;
		status = cl::Platform::get(&platforms);
		OCL_CHECK(status, "\n", status);
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
		OCL_CHECK(status, "\n", status);
		// Get device IDs
		std::vector<cl::Device> devices2;
		status = context.getInfo(CL_CONTEXT_DEVICES, &devices2);
		OCL_CHECK(status, "\n", status);
		devices.push_back(devices2[usedDevices[0]]);
		if (DEBUG) {
			mexPrintBase("devices.size() = %u\n", devices.size());
			mexPrintBase("devices[0] = %u\n", devices[0]);
			mexEval();
		}

		// Create the command queues
		// Enable out of order execution (devices can compute kernels at the same time)
		for (size_t i = 0; i < devices.size(); i++) {
			//commandQueues.push_back(cl::CommandQueue(context, devices[i], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status));
			commandQueues.push_back(cl::CommandQueue(context, devices[i], 0, &status));
			OCL_CHECK(status, "\n", status);
		}
		if (DEBUG) {
			mexPrintBase("commandQueues.size() = %u\n", commandQueues.size());
			mexEval();
		}

		for (cl_uint i = 0; i < commandQueues.size(); i++) {
			commandQueues[i].finish();
		}

		return status;
#endif // END CUDA
	}

	/// <summary>
#if defined(CUDA) || defined(HIP)
	/// This function creates the CUDA programs for the forward and backward projections and for NLM/MRP/RDP/TV/hyperbolic
#else
	/// This function creates the OpenCL programs for the forward and backward projections and for NLM/MRP/RDP/TV
#endif // END CUDA
	/// </summary>
#if !defined(CUDA) && !defined(HIP)
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
#endif // END CUDA
	/// <param name="programFP the program to store forward projection program"></param>
	/// <param name="programBP the program to store backprojection program"></param>
	/// <param name="programAux the program to store auxliary (such as priors) programs"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="w_vec specifies some of the special options used"></param>
	/// <param name="local_size the local size"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline nvrtcResult createProgram(CUmodule & programFP, CUmodule & programBP,
		CUmodule & programAux,
#else
	inline cl_int createProgram(cl::Context & CLContext, cl::Device & CLDeviceID, cl::Program & programFP, cl::Program & programBP,
		cl::Program & programAux, cl::Program & programSens,
#endif // END CUDA
		const char* header_directory, scalarStruct & inputScalars, const RecMethods MethodList,
		const Weighting & w_vec, const size_t local_size[], const int type = -1) {

#if defined(CUDA) || defined(HIP)
		int compMajor = 0, compMinor = 0;
		cuDeviceGetAttribute(&compMajor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, CUDeviceID[0]);
		cuDeviceGetAttribute(&compMinor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, CUDeviceID[0]);
#else
		cl_int status = CL_SUCCESS;
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
		nvrtcResult status = NVRTC_SUCCESS;
#else
		std::string deviceName = CLDeviceID.getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		std::string AMD("Advanced Micro Devices, Inc.");
#endif // END CUDA
		std::string kernelFile = header_directory;
		std::string kernel_path, kernel_pathBP;
		std::string contentFP, contentBP;
		std::string contentAux;
#if defined(CUDA) || defined(HIP)
		std::vector<std::string> options;
		int uu = 0;

#if defined(CUDA)
		std::string buffer0 = "--gpu-architecture=compute_" + std::to_string(compMajor) + std::to_string(compMinor);
		options.push_back(buffer0);
		options.push_back("-DCUDA");
#elif defined(HIP)
		// Compile the HIP kernels for the architecture of the current GPU (e.g. gfx90a).
		hipDeviceProp_t hipProp;
		hipGetDeviceProperties(&hipProp, CUDeviceID[0]);
		std::string buffer0 = std::string("--gpu-architecture=") + hipProp.gcnArchName;
		options.push_back(buffer0);
		options.push_back("-DHIP");
#endif
#else
		std::string options = "-cl-single-precision-constant";
		cl::string apu = CLDeviceID.getInfo<CL_DEVICE_EXTENSIONS>();
		cl::string apu2 = "cl_ext_float_atomics";
		ADD_OPT(options, "-DOPENCL");
#endif // END CUDA
		if (inputScalars.useMAD) {
#if defined(CUDA)
			options.push_back("--use_fast_math");
#elif defined(HIP)
			// HIPRTC is Clang-based and does not accept NVRTC's --use_fast_math; use the Clang flag instead.
			options.push_back("-ffast-math");
#else
			ADD_OPT(options, "-cl-fast-relaxed-math");
#endif // END CUDA
			ADD_OPT(options, "-DUSEMAD");
		}
		if ((inputScalars.useImages && inputScalars.FPType != 4 && inputScalars.FPType != 5 && inputScalars.BPType != 5) || (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 5)) {
			ADD_OPT(options, "-DUSEIMAGES");
			inputScalars.useBuffers = false;
#if defined(CUDA) || defined(HIP)
			useBuffers = false;
#endif // END CUDA
		}
		std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
		// Load the header text file
		std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
		// Orthogonal/volume of intersection header, kept as a distinct header (see below). It is only used
		// by the forward/backward projector programs (projectorType123), never by the auxiliary kernels.
		std::string contentOrth;
		bool useOrth = false;
		// Load orthogonal/volume of intersection headers if applicable
		if (inputScalars.FPType == 2 || inputScalars.BPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 3) {
			if (inputScalars.orthXY)
				ADD_OPT(options, "-DCRYSTXY");
			if (inputScalars.orthZ)
				ADD_OPT(options, "-DCRYSTZ");
			std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
			std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
			contentOrth = contentHeader3;
			useOrth = true;
#if !defined(CUDA) && !defined(HIP)
		}
		if (NV.compare(deviceName) == 0)
			ADD_OPT(options, "-DNVIDIA");
		else if (AMD.compare(deviceName) == 0) {
			cl_bool is_integrated = CLDeviceID.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>(&status);
			if (status == CL_SUCCESS && !is_integrated)
				ADD_OPT(options, "-DAMD");
		}
		else if (apu.find(apu2) != std::string::npos) {
			ADD_OPT(options, "-DINTEL");
#endif // END CUDA
		}

#if defined(CUDA) || defined(HIP)
		// Hand the headers to NVRTC/HIPRTC as actual header files; the kernel sources only #include them.
		nvrtcKernelHeader = contentHeader;
		nvrtcOrthHeader = useOrth ? contentOrth : std::string();
		std::string headerPrefix = "#include \"general_opencl_functions.h\"\n";
		if (useOrth)
			headerPrefix += "#include \"opencl_functions_orth3D.h\"\n";
		// The auxiliary kernels never use the orthogonal/volume header, so they only include the general one.
		const std::string headerPrefixAux = "#include \"general_opencl_functions.h\"\n";
#else
		// OpenCL: keep concatenating the header contents in front of each kernel source.
		const std::string headerPrefix = useOrth ? (contentHeader + contentOrth) : contentHeader;
		// The auxiliary kernels only get the general header (orth3D is excluded).
		const std::string& headerPrefixAux = contentHeader;
#endif // END CUDA

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
			contentFP = headerPrefix + contentFFP;
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
			contentBP = headerPrefix + contentFBP;
		}
		if (inputScalars.useParallelBeam)
			ADD_OPT(options, "-DPARALLEL");
		else if (inputScalars.useHelical)
			ADD_OPT(options, "-DHELICAL");

		// Load the source text file
		// Set all preprocessor definitions
		const bool siddonVal = (inputScalars.FPType == 1 || inputScalars.BPType == 1 || inputScalars.FPType == 4 || inputScalars.BPType == 4) ? true : false;
#if !defined(CUDA) && !defined(HIP)
		if (constantBuffer || (inputScalars.listmode > 0 && !inputScalars.indexBased))
			ADD_OPT(options, "-DUSEGLOBAL");
#endif // END CUDA
		if (inputScalars.raw == 1)
			ADD_OPT(options, "-DRAW");
		if (inputScalars.maskFP) {
			ADD_OPT(options, "-DMASKFP");
			if (inputScalars.maskFPZ > 1)
				ADD_OPT(options, "-DMASKFP3D");
		}
		if (inputScalars.maskBP) {
			ADD_OPT(options, "-DMASKBP");
			if (inputScalars.maskBPZ > 1)
				ADD_OPT(options, "-DMASKBP3D");
		}
		if (inputScalars.useTotLength)// && !inputScalars.SPECT)
			ADD_OPT(options, "-DTOTLENGTH");
		if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights)
			ADD_OPT(options, "-DFDK");
		if (inputScalars.offset)
			ADD_OPT(options, "-DOFFSET");
		if (inputScalars.attenuation_correction == 1u && inputScalars.CTAttenuation)
			ADD_OPT(options, "-DATN");
		else if (inputScalars.attenuation_correction == 1u && !inputScalars.CTAttenuation)
			ADD_OPT(options, "-DATNM");
		if (inputScalars.normalization_correction == 1u)
			ADD_OPT(options, "-DNORM");
		if (inputScalars.scatter == 1u)
			ADD_OPT(options, "-DSCATTER");
		if (inputScalars.randoms_correction == 1u)
			ADD_OPT(options, "-DRANDOMS");
		if (inputScalars.nLayers > 1U) {
			if (inputScalars.listmode > 0 && inputScalars.indexBased)
				ADD_OPT_INT(options, "-DNLAYERS", inputScalars.nLayers);
			else
				ADD_OPT_INT(options, "-DNLAYERS", inputScalars.nProjections / (inputScalars.nLayers * inputScalars.nLayers));
		}
		if (inputScalars.TOF) {
			ADD_OPT(options, "-DTOF");
		}
		if (inputScalars.CT)
			ADD_OPT(options, "-DCT");
		else if (inputScalars.PET && inputScalars.listmode == 0)
			ADD_OPT(options, "-DPET");
		else if (inputScalars.SPECT) {
			ADD_OPT(options, "-DSPECT");
		}

		ADD_OPT_INT(options, "-DNBINS", inputScalars.nBins);
		if (inputScalars.listmode == 1)
			ADD_OPT(options, "-DLISTMODE");
		else if (inputScalars.listmode == 2)
			ADD_OPT(options, "-DLISTMODE2");
		if (inputScalars.listmode > 0 && inputScalars.indexBased)
			ADD_OPT(options, "-DINDEXBASED");
		if ((siddonVal && ((inputScalars.n_rays * inputScalars.n_rays3D) > 1)) || inputScalars.SPECT) {
			ADD_OPT_INT(options, "-DN_RAYS", inputScalars.n_rays * inputScalars.n_rays3D);
			ADD_OPT_INT(options, "-DN_RAYS2D", inputScalars.n_rays);
			ADD_OPT_INT(options, "-DN_RAYS3D", inputScalars.n_rays3D);
		}
		if (inputScalars.pitch)
			ADD_OPT(options, "-DPITCH");
		if (((inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
			ADD_OPT(options, "-DSUBSETS");
		if (local_size[1] > 0ULL) {
			ADD_OPT_INT(options, "-DLOCAL_SIZE", local_size[0]);
			ADD_OPT_INT(options, "-DLOCAL_SIZE2", local_size[1]);
		}
		else {
			ADD_OPT_INT(options, "-DLOCAL_SIZE", local_size[0]);
		}
		if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
			ADD_OPT_INT(options, "-DSTYPE", inputScalars.subsetType);
			ADD_OPT_INT(options, "-DNSUBSETS", inputScalars.subsets);
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
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> os_options = options;
#else
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DSIDDON");
			ADD_OPT(os_options, "-DATOMICF");
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> os_optionsFP = os_options;
#else
			std::string os_optionsFP = os_options;
#endif // END CUDA
			ADD_OPT(os_optionsFP, FP_FLAG);
			if (inputScalars.FPType == 3)
				ADD_OPT(os_optionsFP, "-DVOL");
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3)
				ADD_OPT(os_optionsFP, "-DORTH");
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build FP 1-3 program\n");
				}
#if defined(CUDA) || defined(HIP)
				status = buildProgram(inputScalars.verbose, contentFP, programFP, os_optionsFP);
				if (status == NVRTC_SUCCESS && DEBUG) {
#else
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_optionsFP);
				if (status == CL_SUCCESS && DEBUG) {
#endif // END CUDA
					mexPrint("FP 1-3 program built\n");
				}
#if defined(CUDA) || defined(HIP)
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.FPMod = true;
#endif // END CUDA
				}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build BP 1-3 program\n");
				}
				ADD_OPT(os_options, "-DBP");
				if (inputScalars.BPType == 3)
					ADD_OPT(os_options, "-DVOL");
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
					ADD_OPT(os_options, "-DORTH");
#if defined(CUDA) || defined(HIP)
				status = buildProgram(inputScalars.verbose, contentBP, programBP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
#else
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
				if (status == CL_SUCCESS && DEBUG) {
#endif // END CUDA
					mexPrint("BP 1-3 program built\n");
				}
#if defined(CUDA) || defined(HIP)
				else if (status != NVRTC_SUCCESS)
					return status;
				if (status == NVRTC_SUCCESS)
					memAlloc.BPMod = true;
#endif // END CUDA
				}
			}
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> os_options = options;
#else
			std::string os_options = options;
#endif // END CUDA
			if (inputScalars.FPType == 4)
				ADD_OPT(os_options, FP_FLAG);
			ADD_OPT(os_options, "-DPTYPE4");
			// Non-positive interpolation length selects the Joseph-type projector, where the ray is sampled once per voxel plane along
			// its dominant axis (the input dL is then unused)
			if (inputScalars.dL <= 0.f)
				ADD_OPT(os_options, "-DJOSEPH");
			if (!inputScalars.largeDim) {
				ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
				if (inputScalars.useHelical) {
					ADD_OPT_INT(os_options, "-DNVOXELSHELICAL", NVOXELSHELICAL);
				}
			}
#if defined(CUDA) || defined(HIP)
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
#else
			if (inputScalars.FPType == 4)
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
#endif // END CUDA
			if (inputScalars.BPType == 4) {
				os_options = options;
				if (inputScalars.CT) {
					ADD_OPT(os_options, "-DBP4");
					ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
				}
				else {
					ADD_OPT(os_options, "-DPTYPE4");
					ADD_OPT(os_options, "-DATOMICF");
					if (inputScalars.dL <= 0.f)
						ADD_OPT(os_options, "-DJOSEPH");
				}
				ADD_OPT(os_options, "-DBP");
#if defined(CUDA) || defined(HIP)
				status = buildProgram(inputScalars.verbose, contentBP, programBP, os_options);
				if (status == NVRTC_SUCCESS && DEBUG) {
					mexPrint("BP 4 program built\n");
#else
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
#endif // END CUDA
				}
#if defined(CUDA) || defined(HIP)
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
#else
			else if (inputScalars.CT && inputScalars.BPType == 4 && inputScalars.FPType != 4)
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
#endif // END CUDA
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> os_options = options;
#else
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DPROJ5");
			if (inputScalars.meanFP)
				ADD_OPT(os_options, "-DMEANDISTANCEFP");
			if (inputScalars.meanBP)
				ADD_OPT(os_options, "-DMEANDISTANCEBP");
			if (inputScalars.FPType == 5)
				ADD_OPT(os_options, FP_FLAG);
			if (inputScalars.BPType == 5)
				ADD_OPT(os_options, "-DBP");
			// BDD uses precomputed values, but you can also use it with the original on the fly calculations
			// If GEOM5 is not defined, the kernel uses the old format, otherwise it assumes precomputed values
			if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0)
				ADD_OPT(os_options, "-DGEOM5");
			if (inputScalars.pitch) {
				ADD_OPT_INT(os_options, "-DNVOXELS5", 1);
			}
			else {
				ADD_OPT_INT(os_options, "-DNVOXELS5", NVOXELS5);
			}
			ADD_OPT_INT(os_options, "-DNVOXELSFP", NVOXELSFP);
#if defined(CUDA) || defined(HIP)
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
#else
			if (inputScalars.FPType == 5)
				status = buildProgram(inputScalars.verbose, contentFP, CLContext, CLDeviceID, programFP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
			else
				status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programBP, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
#endif // END CUDA
		}
		if (inputScalars.computeSensImag && inputScalars.listmode > 0) {
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> os_options = options;
#else
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DBP");
			ADD_OPT(os_options, "-DATOMICF");
			ADD_OPT(os_options, "-DSENS");
			if (inputScalars.BPType == 3)
				ADD_OPT(os_options, "-DVOL");
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
				ADD_OPT(os_options, "-DORTH");
			if (inputScalars.BPType == 4) {
				ADD_OPT(os_options, "-DPTYPE4");
				ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
				if (inputScalars.dL <= 0.f)
					ADD_OPT(os_options, "-DJOSEPH");
			}
			else
				ADD_OPT(os_options, "-DSIDDON");
#if defined(CUDA) || defined(HIP)
			status = buildProgram(inputScalars.verbose, contentBP, programSens, os_options);
			if (status == NVRTC_SUCCESS)
				memAlloc.SensMod = true;
#else
			status = buildProgram(inputScalars.verbose, contentBP, CLContext, CLDeviceID, programSens, inputScalars.atomic_64bit, inputScalars.atomic_32bit, os_options);
#endif // END CUDA
		}
		// Build prior programs
		if (MethodList.NLM || MethodList.MRP || MethodList.RDP || w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]
			|| MethodList.TV || MethodList.APLS || MethodList.hyperbolic || MethodList.ProxTV || MethodList.ProxTGV || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM ||
			MethodList.CPType || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF || inputScalars.projector_type == 6 || type == 0) {
			if (DEBUG) {
				mexPrint("Building aux programs\n");
			}
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> optionsAux;
			int uu = 0;
			optionsAux.push_back(buffer0);
#if defined(HIP)
			optionsAux.push_back("-DHIP");
#else
			optionsAux.push_back("-DCUDA");
#endif
			if (inputScalars.useMAD) {
#if defined(HIP)
				// HIPRTC is Clang-based and does not accept NVRTC's --use_fast_math; use the Clang flag instead.
				optionsAux.push_back("-ffast-math");
#else
				optionsAux.push_back("--use_fast_math");
#endif
				optionsAux.push_back("-DUSEMAD");
			}
#endif // END CUDA
			std::string auxKernelPath = kernelFile + "auxKernels.cl";
			std::ifstream sourceFileAux(auxKernelPath.c_str());
			std::string contentAAux((std::istreambuf_iterator<char>(sourceFileAux)), std::istreambuf_iterator<char>());
			// The auxiliary kernels only use the general header, not the orthogonal/volume header.
			contentAux = headerPrefixAux + contentAAux;
#if defined(CUDA) || defined(HIP)
			// Ensure the orthogonal/volume header is not supplied to (nor included in) the aux program.
			nvrtcOrthHeader.clear();
			if (inputScalars.use64BitIndices) {
				optionsAux.push_back("-DLTYPE=long long");
				optionsAux.push_back("-DLTYPE3=long3");
			}
#else
			std::string optionsAux;
			optionsAux = "-cl-single-precision-constant";
			ADD_OPT(optionsAux, "-DOPENCL");
			if (inputScalars.useMAD) {
				ADD_OPT(optionsAux, "-cl-fast-relaxed-math");
				ADD_OPT(optionsAux, "-DUSEMAD");
			}
			if (inputScalars.use64BitIndices) {
				ADD_OPT(optionsAux, "-DLTYPE=long");
				ADD_OPT(optionsAux, "-DLTYPE3=long3");
			}
#endif // END CUDA
			if (inputScalars.largeDim)
				ADD_OPT(optionsAux, "-DLARGEDIM");
			if (inputScalars.useExtendedFOV)
				ADD_OPT(optionsAux, "-DEFOV");
			if (inputScalars.useImages)
				ADD_OPT(optionsAux, "-DUSEIMAGES");
			if (type == 2) {
				if (inputScalars.use_psf)
					ADD_OPT(optionsAux, "-DPSF");
			}
			else if (type == 0) {
				if (inputScalars.CT)
					ADD_OPT(optionsAux, "-DCT");
				if (inputScalars.randoms_correction)
					ADD_OPT(optionsAux, "-DRANDOMS");
				if (inputScalars.use_psf)
					ADD_OPT(optionsAux, "-DPSF");
			}
			else
				ADD_OPT(optionsAux, "-DAF");
			if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
				ADD_OPT(optionsAux, "-DMASKPRIOR");
				if (inputScalars.maskBPZ > 1)
					ADD_OPT(optionsAux, "-DMASKBP3D");
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution)
				ADD_OPT(optionsAux, "-DEFOVZ");
			if (MethodList.MRP) {
				ADD_OPT(optionsAux, "-DMEDIAN");
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_X", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_Y", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_Z", w_vec.Ndz);
			}
			if (MethodList.NLM) {
				ADD_OPT(optionsAux, "-DNLM_");
				if (w_vec.NLM_MRP) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 2);
				}
				else if (w_vec.NLTV) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 1);
				}
				else if (w_vec.NLRD) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 3);
				}
				else if (w_vec.NLLange) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 4);
				}
				else if (w_vec.NLLangeFiltered) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 5);
				}
				else if (w_vec.NLGGMRF) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 6);
				}
				else {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 0);
				}
				if (w_vec.NLAdaptive)
					ADD_OPT(optionsAux, "-DNLMADAPTIVE");
				if (w_vec.NLM_anatomical)
					ADD_OPT(optionsAux, "-DNLMREF");
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
				ADD_OPT_INT(optionsAux, "-DPWINDOWX", w_vec.Nlx);
				ADD_OPT_INT(optionsAux, "-DPWINDOWY", w_vec.Nly);
				ADD_OPT_INT(optionsAux, "-DPWINDOWZ", w_vec.Nlz);
			}
			if (MethodList.GGMRF) {
				ADD_OPT(optionsAux, "-DGGMRF");
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
			}
			if (MethodList.hyperbolic) {
				ADD_OPT(optionsAux, "-DHYPER");
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
			}
			if (MethodList.RDP) {
				ADD_OPT(optionsAux, "-DRDP");
				if (w_vec.RDPLargeNeighbor) {
					ADD_OPT(optionsAux, "-DRDPCORNERS");
					ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
					ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
					ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
				}
				if (w_vec.RDP_anatomical)
					ADD_OPT(optionsAux, "-DRDPREF");
			}
			if (MethodList.ProxRDP && w_vec.RDPLargeNeighbor)
				ADD_OPT(optionsAux, "-DRDPCORNERS");
			if (MethodList.TV && !w_vec.data.TV_use_anatomical) {
				ADD_OPT(optionsAux, "-DTVGRAD");
				if (w_vec.data.TVtype == 6)
					ADD_OPT(optionsAux, "-DTVW1");
				else if (w_vec.data.TVtype == 4)
					ADD_OPT(optionsAux, "-DSATV");
				else if (w_vec.data.TVtype == 2)
					ADD_OPT(optionsAux, "-DJPTV");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			else if ((MethodList.TV && w_vec.data.TV_use_anatomical) || MethodList.APLS) {
				ADD_OPT(optionsAux, "-DTVGRAD");
				if (w_vec.data.TVtype == 1)
					ADD_OPT(optionsAux, "-DANATOMICAL1");
				else if (w_vec.data.TVtype == 2)
					ADD_OPT(optionsAux, "-DANATOMICAL2");
				else if (w_vec.data.TVtype == 5 || MethodList.APLS)
					ADD_OPT(optionsAux, "-DANATOMICAL3");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			if (MethodList.ProxTV) {
				ADD_OPT(optionsAux, "-DPROXTV");
				if (w_vec.UseL2Ball)
					ADD_OPT(optionsAux, "-DL2");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			if (MethodList.ProxTGV || MethodList.TGV) {
				ADD_OPT(optionsAux, "-DPROXTV");
				ADD_OPT(optionsAux, "-DPROXTGV");
				if (w_vec.UseL2Ball)
					ADD_OPT(optionsAux, "-DL2");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
				if (!inputScalars.TGV2D)
					ADD_OPT(optionsAux, "-DTGVZ");
			}
			if (MethodList.ProxRDP)
				ADD_OPT(optionsAux, "-DPROXRDP");
			if (local_sizePrior[1] > 0ULL) {
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE", local_sizePrior[0]);
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE2", local_sizePrior[1]);
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE3", local_sizePrior[2]);
			}
			else {
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE", local_sizePrior[0]);
			}
			if (MethodList.PKMA)
				ADD_OPT(optionsAux, "-DPKMA");
			else if (MethodList.MBSREM || MethodList.MRAMLA)
				ADD_OPT(optionsAux, "-DMBSREM");
			else if (MethodList.BSREM || MethodList.RAMLA)
				ADD_OPT(optionsAux, "-DBSREM");
			else if (MethodList.CPType) {
				ADD_OPT(optionsAux, "-DPDHG");
				if (inputScalars.subsets > 1)
					ADD_OPT(optionsAux, "-DSUBSETS");
			}
			if (inputScalars.projector_type == 6)
				ADD_OPT(optionsAux, "-DROTATE");
#if defined(CUDA) || defined(HIP)
			status = buildProgram(inputScalars.verbose, contentAux, programAux, optionsAux);
			if (status == NVRTC_SUCCESS && DEBUG) {
				mexPrint("Aux program built\n");
			}
			else if (status != NVRTC_SUCCESS)
				return status;
			if (status == NVRTC_SUCCESS)
				memAlloc.auxMod = true;
#else
			status = buildProgram(inputScalars.verbose, contentAux, CLContext, CLDeviceID, programAux, inputScalars.atomic_64bit, inputScalars.atomic_32bit, optionsAux);
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	/// Builds the input CUDA program
#else
	/// Builds the input OpenCL program
#endif // END CUDA
	/// </summary>
	/// <param name="verbose the level of verbosity"></param>
	/// <param name="contentFP program code"></param>
#if !defined(CUDA) && !defined(HIP)
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
#endif // END CUDA
	/// <param name="program the program where to store the built program"></param>
#if !defined(CUDA) && !defined(HIP)
	/// <param name="atomic_64bit are 64-bit (int64) atomics used"></param>
	/// <param name="atomic_32bit are 32-bit (int) atomics used"></param>
#endif // END CUDA
	/// <param name="options preprocessor values for the build"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline nvrtcResult buildProgram(const int8_t verbose, std::string & content, CUmodule & module, std::vector<std::string>&options) {
		nvrtcResult status = NVRTC_SUCCESS;
		CUresult status2 = CUDA_SUCCESS;
		nvrtcProgram program;
		if (DEBUG || verbose >= 3) {
			for (int ll = 0; ll < options.size(); ll++)
				mexPrintBase("%s ", options[ll].c_str());
			mexPrintBase("%s\n", "");
#else
	inline cl_int buildProgram(const int8_t verbose, std::string contentFP, cl::Context & CLContext, cl::Device & CLDeviceID, cl::Program & program,
		bool& atomic_64bit, const bool atomic_32bit, std::string options) {
		cl_int status = CL_SUCCESS;
		size_t pituus;
		if (atomic_64bit) {
			pituus = options.length();
			options += " -DCAST=long";
			options += " -DATOMIC";
			ADD_OPT_INT(options, "-DTH", TH);
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		// Supply the headers as actual headers (resolved by the #include directives in the source) instead
		// of concatenating them into the source string. opencl_functions_orth3D.h is only passed when it is
		// used (i.e. for the projector programs; it is empty for the auxiliary kernels).
		std::vector<const char*> headerSources, headerNames;
		headerSources.push_back(nvrtcKernelHeader.c_str());
		headerNames.push_back("general_opencl_functions.h");
		if (!nvrtcOrthHeader.empty()) {
			headerSources.push_back(nvrtcOrthHeader.c_str());
			headerNames.push_back("opencl_functions_orth3D.h");
		}
		status = nvrtcCreateProgram(&program, content.c_str(), "32bit", static_cast<int>(headerNames.size()), headerSources.data(), headerNames.data());
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
#else
		else if (atomic_32bit) {
			options += " -DCAST=int";
			options += " -DATOMIC32";
			ADD_OPT_INT(options, "-DTH", TH32);
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		// Build the program
		std::vector<const char*> optionsC;
		optionsC.reserve(options.size());
		for (const std::string& opt : options)
			optionsC.push_back(opt.c_str());
		status = nvrtcCompileProgram(program, optionsC.size(), optionsC.data());
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
#else
		else {
			options += " -DCAST=float";
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		else if (verbose > 1)
			mexPrint("CUDA program built\n");
			size_t ptxSize;
			status = nvrtcGetPTXSize(program, &ptxSize);
			if (status != NVRTC_SUCCESS) {
				std::cerr << nvrtcGetErrorString(status) << std::endl;
				return status;
#else
		if (DEBUG || verbose >= 3)
			mexPrintBase("%s\n", options.c_str());
		if (atomic_64bit) {
			cl::string apu = CLDeviceID.getInfo<CL_DEVICE_EXTENSIONS>();
			cl::string apu2 = "cl_khr_int64_base_atomics";
			size_t var = apu.find(apu2);
			if (var < 0) {
				options.erase(pituus, options.size() + 1);
				options += " -DCAST=float";
				status = -1;
#endif // END CUDA
			}
#if defined(CUDA) || defined(HIP)
			char* ptx = new char[ptxSize];
			status = nvrtcGetPTX(program, ptx);
			if (status != NVRTC_SUCCESS) {
				std::cerr << nvrtcGetErrorString(status) << std::endl;
				return status;
#else
			else {
				std::vector<std::string> testi;
				testi.push_back(contentFP);
				cl::Program::Sources source(testi);
				program = cl::Program(CLContext, source);
				status = program.build({ CLDeviceID }, options.c_str());
				if (status == CL_SUCCESS && (DEBUG || verbose >= 3)) {
					mexPrint("OpenCL program (64-bit atomics) built\n");
#endif // END CUDA
				}
#if defined(CUDA) || defined(HIP)
				status2 = cuModuleLoadDataEx(&module, ptx, 0, 0, 0);
				CUDA_CHECK(status2, "\n", NVRTC_ERROR_BUILTIN_OPERATION_FAILURE);
#else
				else if (status != CL_SUCCESS) {
					mexPrint("Failed to build 64-bit atomics program.\n");
#endif // END CUDA
					if (DEBUG) {
#if defined(CUDA) || defined(HIP)
						mexPrintBase("ptxSize = %u\n", ptxSize);
#else
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
#endif // END CUDA
						}
#if defined(CUDA) || defined(HIP)
						// Destroy the program.
						status = nvrtcDestroyProgram(&program);
						if (status != NVRTC_SUCCESS) {
							std::cerr << nvrtcGetErrorString(status) << std::endl;
							return status;
#endif // END CUDA
						}
#if defined(CUDA) || defined(HIP)
						delete[] ptx;
#else
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
			status = program.build({ CLDeviceID }, options.c_str());
			if (status == CL_SUCCESS && (DEBUG || verbose >= 3)) {
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
#endif // END CUDA
		return status;
			}

		/// <summary>
#if defined(CUDA) || defined(HIP)
	/// Creates the necessary CUDA kernels from the input programs
#else
	/// Creates the necessary OpenCL kernels from the input programs
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
		inline CUresult createKernels(CUfunction & kernelFP, CUfunction & kernelBP, CUfunction & kernelNLM, CUfunction & kernelMed,
			CUfunction & kernelRDP, CUfunction & kernelGGMRF, const CUmodule & programFP, const CUmodule & programBP, const CUmodule & programAux,
			const RecMethods & MethodList, const Weighting & w_vec, const scalarStruct & inputScalars, const int type = -1) {
			CUresult status = CUDA_SUCCESS;
#else
		inline cl_int createKernels(cl::Kernel & kernelFP, cl::Kernel & kernelBP, cl::Kernel & kernelNLM, cl::Kernel & kernelMed,
			cl::Kernel & kernelRDP, cl::Kernel & kernelGGMRF, const cl::Program & programFP, const cl::Program & programBP, const cl::Program & programAux,
			const cl::Program & programSens, const RecMethods & MethodList, const Weighting & w_vec, const scalarStruct & inputScalars, const int type = -1) {
			cl_int status = CL_SUCCESS;
#endif // END CUDA
			// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
			if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
				if (inputScalars.FPType == 4) {
					CREATE_KERNEL(kernelFP, programFP, "projectorType4Forward", "Failed to create projector type 4 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 4 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 4) {
					if (inputScalars.FPType == 4 && inputScalars.CT)
						GET_KERNEL(kernelBP, programBP, "projectorType4Backward");
					else if (!inputScalars.CT)
						GET_KERNEL(kernelBP, programBP, "projectorType4Forward");
					else
						CREATE_KERNEL(kernelBP, programBP, "projectorType4Backward", "Failed to create projector type 4 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 4 BP successfully created\n");
					}
				}
			}
			if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
				if (inputScalars.FPType == 5) {
					CREATE_KERNEL(kernelFP, programFP, "projectorType5Forward", "Failed to create projector type 5 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 5 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 5) {
					if (inputScalars.FPType == 5)
						GET_KERNEL(kernelBP, programFP, "projectorType5Backward");
					else
						GET_KERNEL(kernelBP, programBP, "projectorType5Backward");
					KCHECK("Failed to create projector type 5 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 5 BP successfully created\n");
					}
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
					GET_KERNEL(kernelFP, programFP, "projectorType123");
					KCHECK("Failed to create projector type 1-3 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 1-3 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					GET_KERNEL(kernelBP, programBP, "projectorType123");
					KCHECK("Failed to create projector type 1-3 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 1-3 BP successfully created\n");
					}
				}
			}
		if (MethodList.NLM) {
			CREATE_KERNEL(kernelNLM, programAux, "NLM", "Failed to create NLM kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("NLM kernel successfully created\n");
			}
		}
		if (MethodList.MRP) {
			CREATE_KERNEL(kernelMed, programAux, "medianFilter3D", "Failed to create Median kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Median kernel successfully created\n");
			}
		}
		if (MethodList.RDP) {
			CREATE_KERNEL(kernelRDP, programAux, "RDPKernel", "Failed to create RDP kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("RDP kernel successfully created\n");
			}
		}
		if (MethodList.GGMRF) {
			CREATE_KERNEL(kernelGGMRF, programAux, "GGMRFKernel", "Failed to create GGMRF kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("GGMRF kernel successfully created\n");
			}
		}
		if (MethodList.TV || MethodList.APLS) {
			CREATE_KERNEL(kernelTV, programAux, "TVKernel", "Failed to create TV kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("TV kernel successfully created\n");
			}
		}
		if (MethodList.hyperbolic) {
			CREATE_KERNEL(kernelHyper, programAux, "hyperbolicKernel", "Failed to create hyperbolic prior kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Hyperbolic prior kernel successfully created\n");
			}
		}
		if (MethodList.PKMA || MethodList.BSREM || MethodList.MBSREM || MethodList.MRAMLA || MethodList.RAMLA) {
			CREATE_KERNEL(kernelPoisson, programAux, "PoissonUpdate", "Failed to create Poisson Update kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Poisson Update kernel successfully created\n");
			}
		}
		if (MethodList.CPType) {
			CREATE_KERNEL(kernelPDHG, programAux, "PDHGUpdate", "Failed to create PDHG Update kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("PDHG Update kernel successfully created\n");
			}
		}
		if (MethodList.ProxTV) {
			GET_KERNEL(kernelProxTVq, programAux, "ProxTVq");
			GET_KERNEL(kernelProxTVDiv, programAux, "ProxTVDivergence");
			GET_KERNEL(kernelProxTVGrad, programAux, "ProxTVGradient");
			CHECK(status, "Failed to create proximal TV kernel\n", status);
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal TV kernel successfully created\n");
			}
		}
		if (MethodList.ProxRDP) {
			GET_KERNEL(kernelProxq, programAux, "Proxq");
			GET_KERNEL(kernelProxRDP, programAux, "ProxRDP");
			GET_KERNEL(kernelProxTrans, programAux, "ProxTrans");
			KCHECK("Failed to create proximal RDP kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal RDP kernel successfully created\n");
			}
		}
		if (MethodList.ProxNLM) {
			GET_KERNEL(kernelProxq, programAux, "Proxq");
			GET_KERNEL(kernelProxNLM, programAux, "ProxNLM");
			GET_KERNEL(kernelProxTrans, programAux, "ProxTrans");
			KCHECK("Failed to create proximal NLM kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal NLM kernel successfully created\n");
			}
		}
		if (MethodList.ProxTGV) {
			GET_KERNEL(kernelProxTVq, programAux, "ProxTVq");
			GET_KERNEL(kernelProxTGVq, programAux, "ProxTGVq");
			GET_KERNEL(kernelProxTVDiv, programAux, "ProxTVDivergence");
			GET_KERNEL(kernelProxTVGrad, programAux, "ProxTVGradient");
			GET_KERNEL(kernelProxTGVDiv, programAux, "ProxTGVDivergence");
			GET_KERNEL(kernelProxTGVSymmDeriv, programAux, "ProxTGVSymmDeriv");
			CHECK(status, "Failed to create proximal TGV kernel\n", status);
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal TGV kernel successfully created\n");
			}
		}
		if (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) {
			CREATE_KERNEL(kernelElementMultiply, programAux, "vectorElementMultiply", "Failed to create element-wise multiplication kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Element-wise kernels successfully created\n");
			}
			CREATE_KERNEL(kernelElementDivision, programAux, "vectorElementDivision", "Failed to create element-wise division kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Element-wise kernels successfully created\n");
			}
		}
#if !defined(CUDA) && !defined(HIP)
		if (type == 0) {
			kernelsumma = cl::Kernel(programAux, "summa", &status);
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			kernelEstimate = cl::Kernel(programAux, "computeEstimate", &status);
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			kernelForward = cl::Kernel(programAux, "forward", &status);
			if (inputScalars.use_psf) {
				kernelPSFf = cl::Kernel(programAux, "Convolution3D_f", &status);
				kernelPSF = cl::Kernel(programAux, "Convolution3D", &status);
			}
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			if (DEBUG || inputScalars.verbose >= 3) {
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
#endif // END CUDA
		if (inputScalars.computeSensImag && inputScalars.listmode > 0) {
			if (inputScalars.BPType == 4)
				GET_KERNEL(kernelSensList, programSens, "projectorType4Forward");
			else
				CREATE_KERNEL(kernelSensList, programSens, "projectorType123", "Failed to create sensitivity image kernels\n");
		}
		if (inputScalars.projector_type == 6) {
			CREATE_KERNEL(kernelRotate, programAux, "rotate", "Failed to create bilinear rotation kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Bilinear rotation kernel successfully created\n");
			}
		}
		return status;
	}
public:
#if defined(CUDA) || defined(HIP)
	std::vector<CUdevice> CUDeviceID;
	std::vector<CUstream> CLCommandQueue;
	CUfunction kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, kernelElementMultiply, kernelElementDivision,
#else
	cl::Context CLContext;
	std::vector<cl::Device> CLDeviceID;
	std::vector<cl::CommandQueue> CLCommandQueue;
	OpenCL_im_vectors vec_opencl;
	cl::Kernel kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, kernelElementMultiply, kernelElementDivision,
#endif // END CUDA
		kernelTV, kernelProxTGVSymmDeriv, kernelProxTGVDiv, kernelProxTGVq, kernelPoisson, kernelPDHG, kernelProxRDP, kernelProxq, kernelProxTrans, kernelProxNLM, kernelGGMRF,
		kernelsumma, kernelEstimate, kernelPSF, kernelPSFf, kernelDiv, kernelMult, kernelForward, kernelSensList, kernelApu, kernelHyper, kernelRotate;
	// Device buffers / pointers common to both backends (plain, non-pointer)
	DEVBUF_t d_xcenter, d_ycenter, d_zcenter, d_V, d_TOFCenter, d_eFOVIndices, d_weights, d_angle, d_g, d_uref, d_maskBPB, d_rayShiftsDetector, d_rayShiftsSource, d_maskPriorB;
	std::vector<DEVBUF_t> d_attenB;
#if defined(CUDA) || defined(HIP)
	CUmodule programFP, programBP, programAux, programSens;
	CUdeviceptr* d_output, * d_meanBP, * d_meanFP, * d_inputB, * d_W, * d_gaussianNLM;
	CUdeviceptr* d_qX, * d_qY, * d_qZ;
	CUdeviceptr* d_rX, * d_rY, * d_rXY, * d_rZ, * d_rXZ, * d_rYZ;
	CUdeviceptr* d_vX, * d_vY, * d_vZ;
	CUdeviceptr* d_vector, * d_input;
	CUdeviceptr* d_im, * d_rhs, * d_U, * d_refIm, * d_RDPref;
	CUtexObject d_maskFP, d_maskBP, d_maskPrior;
	CUtexObject d_inputImage, d_imageX, d_imageY, d_urefIm, d_inputI, d_RDPrefI;
	CUarray atArray, uRefArray, maskArrayBP, maskArrayPrior, BPArray, FPArray, integArrayXY, imArray;
	std::vector<CUtexObject> d_attenIm;
	std::vector<void*> FPArgs, BPArgs, SensArgs;
	CUDA_im_vectors vec_opencl;
#else
	cl::Buffer d_output, d_meanBP, d_meanFP, d_inputB, d_W, d_gaussianNLM;
	cl::Buffer d_qX, d_qY, d_qZ;
	cl::Buffer d_rX, d_rY, d_rXY, d_rZ, d_rXZ, d_rYZ;
	cl::Buffer d_vX, d_vY, d_vZ;
	cl::Buffer d_vector, d_input;
	cl::Buffer d_im, d_rhs, d_U, d_refIm, d_RDPref;
	cl::Image2D d_maskFP, d_maskBP, d_maskPrior;
	cl::Image3D d_maskBP3, d_maskPrior3;
	cl::Image3D d_inputImage, d_urefIm, d_inputI, d_RDPrefI;
	std::vector<cl::Image3D> d_attenIm;
	cl::Buffer d_outputCT;
	size_t memSize = 0ULL;
#endif // END CUDA
	std::chrono::steady_clock::time_point tStartLocal, tStartGlobal, tStartAll;
	std::chrono::steady_clock::time_point tEndLocal, tEndGlobal, tEndAll;
	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
#if defined(CUDA) || defined(HIP)
	std::vector<float3> b, d, bmax;
	// Image dimensions
	std::vector<int3> d_N;
#else
	std::vector<cl_float3> b, d, bmax;
	std::vector<cl_int3> d_N;
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
	std::vector<CUtexObject> d_maskFP3;
	std::vector<CUarray> maskArrayFP;
	std::vector<CUdeviceptr*> d_Summ;
	CUDA_TEXTURE_DESC texDesc;
	CUDA_ARRAY_DESCRIPTOR arr2DDesc;
	CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
	CUDA_RESOURCE_DESC resDesc;
	CUDA_RESOURCE_VIEW_DESC viewDesc;
	unsigned char no_norm = 0;
	int proj6 = 1;
	size_t memSize = 0ULL;
#else
	std::vector<cl::Image3D> d_maskFP3;
	std::vector<cl::Buffer> d_Summ;
	std::vector<cl::Buffer> d_LFull, d_zindexFull, d_xyindexFull;
	std::vector<cl::Buffer> d_meas;
	std::vector<cl::Buffer> d_rand;
	std::vector<cl::Buffer> d_imTemp;
	std::vector<cl::Buffer> d_imFinal;
	// Image origin
	cl::detail::size_t_array origin = { { 0, 0, 0 } };
	cl::detail::size_t_array region = { { 0, 0, 0 } };
	// Image format
	cl::ImageFormat format;
	cl::ImageFormat formatMask;
#endif // END CUDA
	// Vector device buffers common to both backends
	std::vector<DEVBUF_t> d_maskFPB;
	std::vector<DEVBUF_t> d_normFull, d_scatFull, d_xFull, d_zFull;
	std::vector<DEVBUF_t> d_L;
	std::vector<DEVBUF_t> d_zindex, d_xyindex, d_norm, d_atten, d_T;
	// Precomputed per-projection geometry for the BDD backprojection (16 floats per projection, one buffer per subset)
	std::vector<DEVBUF_t> d_geomProj5;
	// Host-side storage for the above geometry
	std::vector<std::vector<float>> geomProj5Host;
	std::vector<std::vector<DEVBUF_t>> d_scat, d_x, d_z, d_trIndex, d_axIndex, d_TOFIndex;
	std::vector<std::vector<size_t>> erotusBP, erotusPDHG;
#if defined(CUDA) || defined(HIP)
	// This is used to define the additional queues/streams in the multi-resolution case
	// Note that these are only used in the multi-resolution case
	std::vector<CUstream> sideQueues;
	CUevent evMain = nullptr;
	std::vector<CUevent> evSide;
	// Persistent per-volume FP input arrays/textures
	// Previously these were always recreated
	std::vector<CUarray> FPArrayCache;
	std::vector<CUtexObject> FPTexCache;
	std::vector<CUarray> intArrayCacheXY;
	std::vector<CUtexObject> intTexCacheXY;
	std::vector<CUarray> intArrayCacheYZ;
	std::vector<CUtexObject> intTexCacheYZ;
#else
	std::vector<cl::CommandQueue> sideQueues;
	// Persistent per-volume FP input images
	std::vector<cl::Image3D> d_imageCache;
	std::vector<cl::Image3D> d_intImageCacheXY;
	std::vector<cl::Image3D> d_intImageCacheYZ;
#endif
	// Cached dimensions of the above (3 entries per volume; third entry 0 = not yet created)
	std::vector<size_t> imageCacheDims, intImageCacheDimsXY, intImageCacheDimsYZ;
	// Cached dimensions of the persistent BP input image/texture (third entry 0 = not yet created)
	size_t BPImageDims[3] = { 0, 0, 0 };
	// ==== End Mod additions ====
#if defined(CUDA) || defined(HIP)
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
		if (memAlloc.attenM) {
			for (int kk = 0; kk < memAlloc.aSteps; kk++) {
				getErrorString(cuMemFree(d_atten[kk]));
			}
		}
		if (memAlloc.V)
			getErrorString(cuMemFree(d_V));
		if (memAlloc.atten && !useBuffers) {
			getErrorString(cuArrayDestroy(atArray));
		}
		for (int tt = 0; tt < memAlloc.tSteps; tt++) {
			if (memAlloc.atten && memAlloc.attenSize < tt) {
				if (useBuffers) {
					getErrorString(cuMemFree(d_attenB[tt]));
				}
				else {
					getErrorString(cuTexObjectDestroy(d_attenIm[tt]));
				}
			}
			if (memAlloc.xSteps >= 0) {
				for (int kk = 0; kk <= memAlloc.xSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_x[tt][kk]));
				}
			}
			if (memAlloc.zType == 0) {
				getErrorString(cuMemFree(d_z[tt][memAlloc.zSteps]));
			}
			else if (memAlloc.zType == 1) {
				for (int kk = 0; kk <= memAlloc.zSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_z[tt][kk]));
				}
			}
			if (memAlloc.extra) {
				for (int kk = 0; kk < memAlloc.eSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_scat[tt][kk]));
				}
			}
			if (memAlloc.indexBased) {
				// When the data is loaded one subset at a time, only a single buffer at [0][0] exists
				int nIndex = memAlloc.iSteps / memAlloc.tSteps;
				if (nIndex == 0)
					nIndex = tt == 0 ? memAlloc.iSteps : 0;
				for (int kk = 0; kk < nIndex; kk++) {
					getErrorString(cuMemFree(d_trIndex[tt][kk]));
					getErrorString(cuMemFree(d_axIndex[tt][kk]));
				}
			}
			if (memAlloc.TOFIndex) {
				// When the data is loaded one subset at a time, only a single buffer at [0][0] exists
				int nTOF = memAlloc.TOFSteps / memAlloc.tSteps;
				if (nTOF == 0)
					nTOF = tt == 0 ? memAlloc.TOFSteps : 0;
				for (int kk = 0; kk < nTOF; kk++) {
					getErrorString(cuMemFree(d_TOFIndex[tt][kk]));
				}
			}
		}
		if (memAlloc.offsetT) {
			for (int kk = 0; kk < memAlloc.oSteps; kk++) {
				getErrorString(cuMemFree(d_T[kk]));
			}
		}
		if (memAlloc.geom5) {
			for (int kk = 0; kk < memAlloc.g5Steps; kk++) {
				getErrorString(cuMemFree(d_geomProj5[kk]));
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
		for (size_t kk = 0; kk < FPTexCache.size(); kk++) {
			if (imageCacheDims[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(FPTexCache[kk]));
				getErrorString(cuArrayDestroy(FPArrayCache[kk]));
			}
		}
		for (size_t kk = 0; kk < intTexCacheXY.size(); kk++) {
			if (intImageCacheDimsXY[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(intTexCacheXY[kk]));
				getErrorString(cuArrayDestroy(intArrayCacheXY[kk]));
			}
		}
		for (size_t kk = 0; kk < intTexCacheYZ.size(); kk++) {
			if (intImageCacheDimsYZ[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(intTexCacheYZ[kk]));
				getErrorString(cuArrayDestroy(intArrayCacheYZ[kk]));
			}
		}
		if (BPImageDims[2] != 0) {
			getErrorString(cuTexObjectDestroy(d_inputImage));
			getErrorString(cuArrayDestroy(BPArray));
		}
		for (size_t kk = 0; kk < sideQueues.size(); kk++)
			getErrorString(cuStreamDestroy(sideQueues[kk]));
		for (size_t kk = 0; kk < evSide.size(); kk++)
			getErrorString(cuEventDestroy(evSide[kk]));
		if (evMain != nullptr)
			getErrorString(cuEventDestroy(evMain));
	}
#else
	cl_uchar no_norm = 0;
	int proj6 = 1;
	~ProjectorClass() {}
#endif // END CUDA

	/// <summary>
	/// Create the additional queues/streams for the multi-resolution case,
	//// The main queue/stream is always the ArrayFire queue/stream, the other volumes use their own ones
	/// </summary>
	/// <param name="n number of additional multi-resolution volumes (inputScalars.nMultiVolumes)"></param>
	inline int initSideQueues(const int n) {
		if (n <= 0 || static_cast<int>(sideQueues.size()) >= n)
			return 0;
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		if (evMain == nullptr) {
			status = cuEventCreate(&evMain, CU_EVENT_DISABLE_TIMING);
			CUDA_CHECK(status, "Failed to create main stream event\n", -1);
		}
		for (int kk = static_cast<int>(sideQueues.size()); kk < n; kk++) {
			CUstream s;
			// Create the stream
			status = cuStreamCreate(&s, CU_STREAM_NON_BLOCKING);
			CUDA_CHECK(status, "Failed to create side stream\n", -1);
			sideQueues.push_back(s);
			CUevent e;
			// Create the event to make sure all relevant computations are complete after the BP
			status = cuEventCreate(&e, CU_EVENT_DISABLE_TIMING);
			CUDA_CHECK(status, "Failed to create side stream event\n", -1);
			evSide.push_back(e);
		}
#else
		cl_int status = CL_SUCCESS;
		for (int kk = static_cast<int>(sideQueues.size()); kk < n; kk++) {
			sideQueues.push_back(cl::CommandQueue(CLContext, CLDeviceID[0], 0, &status));
			OCL_CHECK(status, "Failed to create side command queue\n", -1);
		}
#endif
		return 0;
	}

	/// <summary>
	/// This function forces the main queue/stream to wait for all the side queues/streams
	// This makes sure that the computations will only proceed after all the BPs are done
	/// </summary>
	inline int joinSideQueues() {
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		for (size_t kk = 0; kk < sideQueues.size(); kk++) {
			status = cuEventRecord(evSide[kk], sideQueues[kk]);
			CUDA_CHECK(status, "Failed to record side stream event\n", -1);
			status = cuStreamWaitEvent(CLCommandQueue[0], evSide[kk], 0);
			CUDA_CHECK(status, "Failed to make main stream wait for side stream\n", -1);
		}
#else
		cl_int status = CL_SUCCESS;
		if (sideQueues.size() > 0) {
			std::vector<cl::Event> waitList(sideQueues.size());
			for (size_t kk = 0; kk < sideQueues.size(); kk++) {
				status = sideQueues[kk].enqueueMarkerWithWaitList(nullptr, &waitList[kk]);
				OCL_CHECK(status, "Failed to enqueue side queue marker\n", -1);
				status = sideQueues[kk].flush();
				OCL_CHECK(status, "Failed to flush side queue\n", -1);
			}
			status = CLCommandQueue[0].enqueueBarrierWithWaitList(&waitList);
			OCL_CHECK(status, "Failed to enqueue main queue barrier\n", -1);
		}
#endif
		return 0;
	}
	
	/// <summary>
	/// This function creates the projector class object
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <returns></returns>
	inline int addProjector(scalarStruct & inputScalars, Weighting & w_vec, const RecMethods & MethodList, const char* header_directory, const int type = -1) {
		// Set-up the local group size
#if defined(CUDA) || defined(HIP)
		local_size[0] = 32ULL;
#else
		local_size[0] = 64ULL;
#endif // END CUDA
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
		// Override the above defaults with user-supplied values (per dimension). A negative input value
		// keeps the corresponding default computed above.
		for (int ls = 0; ls < 3; ls++)
			if (inputScalars.localSize[ls] > 0)
				local_size[ls] = static_cast<size_t>(inputScalars.localSize[ls]);
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
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		nvrtcResult status2 = NVRTC_SUCCESS;
#else
		cl_int status = CL_SUCCESS;
#endif // END CUDA
		proj6 = 0;

#if defined(CUDA) || defined(HIP)
		// Create the CUDA/HIP context and stream and assign the device
		int af_id = af::getDevice();
		CUDeviceID.push_back(afcu::getNativeId(af_id));
		CLCommandQueue.push_back(afcu::getStream(CUDeviceID[0]));

		status2 = createProgram(programFP, programBP, programAux, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		if (status2 != NVRTC_SUCCESS) {
			std::cerr << "Error while creating program" << std::endl;
			return -1;
		}
#else
		// Create the OpenCL context and command queue and assign the device
#ifdef AF
		CLContext = afcl::getContext(true);
		CLCommandQueue.push_back(cl::CommandQueue(afcl::getQueue(true), true));
		// Use AF command queue to query for the current OpenCL device
		CLDeviceID.push_back(CLCommandQueue[0].getInfo<CL_QUEUE_DEVICE>(&status));
		OCL_CHECK(status, "\n", -1);
#else
		status = clGetPlatformsContext(inputScalars.platform, CLContext, CLCommandQueue, inputScalars.usedDevices, CLDeviceID);
#endif
		// For NVIDIA cards, 32 local size seems more optimal with 1D kernelFP (unless the user gave an explicit value)
		std::string deviceName = CLDeviceID[0].getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		if (inputScalars.localSize[0] <= 0 && NV.compare(deviceName) == 0 && (inputScalars.projector_type == 1 || inputScalars.projector_type == 11) && local_size[1] == 1ULL)
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
#endif // END CUDA
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("CUDA programs successfully created\n");
		}
#else
		cl_ulong constantBufferSize = CLDeviceID[0].getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>(&status);

		if ((inputScalars.size_of_x + inputScalars.size_z) * sizeof(float) >= constantBufferSize)
			constantBuffer = true;
		if (DEBUG) {
			mexPrintBase("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE = %u\n", constantBufferSize);
			mexPrintBase("(inputScalars.size_of_x + inputScalars.size_z) * sizeof(float) = %u\n", (inputScalars.size_of_x + inputScalars.size_z) * sizeof(float));
			mexPrintBase("inputScalars.size_of_x = %u\n", inputScalars.size_of_x);
			mexPrintBase("inputScalars.size_z = %u\n", inputScalars.size_z);
			mexEval();
		}
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, MethodList, w_vec, inputScalars, type);
		CUDA_CHECK(status, "Failed to create kernels\n", -1);
#else
		cl::Program programFP, programBP, programAux, programSens;

		status = createProgram(CLContext, CLDeviceID[0], programFP, programBP, programAux, programSens, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		OCL_CHECK(status, "Error while creating program\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			mexPrint("CUDA kernels successfully created\n");
#else
			mexPrint("OpenCL programs successfully created\n");
#endif // END CUDA
		}
#if !defined(CUDA) && !defined(HIP)
		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, programSens, MethodList, w_vec, inputScalars, type);
		OCL_CHECK(status, "Failed to create kernels\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("OpenCL kernels successfully created\n");
		}
		format.image_channel_order = CL_A;
		format.image_channel_data_type = CL_FLOAT;
		formatMask.image_channel_order = CL_A;
		formatMask.image_channel_data_type = CL_UNSIGNED_INT8;
#endif // END CUDA

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
#if defined(CUDA) || defined(HIP)
			globalPriorEFOV[0] = (inputScalars.NxPrior + erotusPriorEFOV[0]) / local_sizePrior[0];
			globalPriorEFOV[1] = (inputScalars.NyPrior + erotusPriorEFOV[1]) / local_sizePrior[1];
			globalPriorEFOV[2] = (inputScalars.NzPrior + erotusPriorEFOV[2]) / local_sizePrior[2];
#else
			globalPriorEFOV = { inputScalars.NxPrior + erotusPriorEFOV[0], inputScalars.NyPrior + erotusPriorEFOV[1], inputScalars.NzPrior + erotusPriorEFOV[2] };
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
		local[0] = local_size[0];
		local[1] = local_size[1];
		local[2] = 1;
		localPrior[0] = local_sizePrior[0];
		localPrior[1] = local_sizePrior[1];
		localPrior[2] = local_sizePrior[2];
#else
		local = { local_size[0] , local_size[1] };
		localPrior = { local_sizePrior[0] , local_sizePrior[1], local_sizePrior[2] };
#endif // END CUDA
		erotusPrior[0] = inputScalars.Nx[0] % local_sizePrior[0];
		erotusPrior[1] = inputScalars.Ny[0] % local_sizePrior[1];
		erotusPrior[2] = inputScalars.Nz[0] % local_sizePrior[2];
		if (erotusPrior[0] > 0)
			erotusPrior[0] = (local_sizePrior[0] - erotusPrior[0]);
		if (erotusPrior[1] > 0)
			erotusPrior[1] = (local_sizePrior[1] - erotusPrior[1]);
		if (erotusPrior[2] > 0)
			erotusPrior[2] = (local_sizePrior[1] - erotusPrior[2]);
#if defined(CUDA) || defined(HIP)
		globalPrior[0] = (inputScalars.Nx[0] + erotusPrior[0]) / localPrior[0];
		globalPrior[1] = (inputScalars.Ny[0] + erotusPrior[1]) / localPrior[1];
		globalPrior[2] = (inputScalars.Nz[0] + erotusPrior[2]) / localPrior[2];
		d_NOrig = make_vec3<int3>(static_cast<int>(inputScalars.NxOrig), static_cast<int>(inputScalars.NyOrig), static_cast<int>(inputScalars.NzOrig));
		d_NPrior = make_vec3<int3>(static_cast<int>(inputScalars.NxPrior), static_cast<int>(inputScalars.NyPrior), static_cast<int>(inputScalars.NzPrior));
#else
		globalPrior = { inputScalars.Nx[0] + erotusPrior[0], inputScalars.Ny[0] + erotusPrior[1], inputScalars.Nz[0] + erotusPrior[2] };
		d_NOrig = { static_cast<cl_int>(inputScalars.NxOrig), static_cast<cl_int>(inputScalars.NyOrig), static_cast<cl_int>(inputScalars.NzOrig) };
		d_NPrior = { static_cast<cl_int>(inputScalars.NxPrior), static_cast<cl_int>(inputScalars.NyPrior), static_cast<cl_int>(inputScalars.NzPrior) };
#endif // END CUDA
		dPitch = { w_vec.dPitchX, w_vec.dPitchY };
		if (inputScalars.SPECT) {
			totalFOVmin = { inputScalars.totalFOVxmin, inputScalars.totalFOVymin, inputScalars.totalFOVzmin };
			totalFOVmax = { inputScalars.totalFOVxmax, inputScalars.totalFOVymax, inputScalars.totalFOVzmax };
		}
		b.resize(inputScalars.nMultiVolumes + 1);
		d.resize(inputScalars.nMultiVolumes + 1);
		d_N.resize(inputScalars.nMultiVolumes + 1);
		bmax.resize(inputScalars.nMultiVolumes + 1);
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
#if defined(CUDA) || defined(HIP)
			b[ii] = make_vec3<float3>(inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii]);
			d[ii] = make_vec3<float3>(inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii]);
			d_N[ii] = make_vec3<int3>(static_cast<int>(inputScalars.Nx[ii]), static_cast<int>(inputScalars.Ny[ii]), static_cast<int>(inputScalars.Nz[ii]));
			bmax[ii] = make_vec3<float3>(static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
#else
			b[ii] = { inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii] };
			d[ii] = { inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii] };
			d_N[ii] = { static_cast<cl_int>(inputScalars.Nx[ii]), static_cast<cl_int>(inputScalars.Ny[ii]), static_cast<cl_int>(inputScalars.Nz[ii]) };
			bmax[ii] = { static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
#endif // END CUDA
				static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii],
#if defined(CUDA) || defined(HIP)
				static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii]);
#else
				static_cast<float>(inputScalars.Nz[ii])* inputScalars.dz[ii] + inputScalars.bz[ii] };
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
		if (DEBUG)
			mexPrint("Luuppi valmis\n");
		d_Summ.resize(1);
		d_Summ[0] = nullptr;
#else
		if (d_Summ.size() < 1)
			d_Summ.resize(1);
		region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
#endif // END CUDA
		return 0;
		}

	/// <summary>
#if defined(CUDA) || defined(HIP)
	/// This function first creates the necessary CUDA buffers and then writes the data to them
#else
	/// This function first creates the necessary OpenCL buffers and then writes the data to them
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int createAndWriteBuffers(const std::vector<int64_t>&length, const float* x, const float* z_det, const uint32_t * xy_index,
#else
	inline cl_int createAndWriteBuffers(const std::vector<int64_t>&length, const float* x, const float* z_det, const uint32_t * xy_index,
#endif // END CUDA
		const uint16_t * z_index, const uint16_t * L, const int64_t * pituus, const float* atten, const float* norm, const float* extraCorr,
		const scalarStruct & inputScalars, const Weighting & w_vec, const RecMethods & MethodList) {
		STATUS_t status = SUCCESS_VALUE;
		size_t vecSize = 1;
#if !defined(CUDA) && !defined(HIP)
		cl::size_type imX = inputScalars.Nx[0];
		cl::size_type imY = inputScalars.Ny[0];
		cl::size_type imZ = inputScalars.Nz[0];
#endif // END CUDA
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
#if !defined(CUDA) && !defined(HIP)
		// NLM anatomical reference image
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (inputScalars.useImages)
				d_urefIm = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			else
				d_uref = cl::Buffer(CLContext, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0], NULL, &status);
			OCL_CHECK(status, "\n", -1);
		}
		// The input image for numerous regularization
		// We define the image here and later input the current estimate to this image
		if (MethodList.NLM || MethodList.RDP || MethodList.TV || MethodList.GGMRF || MethodList.APLS || MethodList.hyperbolic || inputScalars.projector_type == 6) {
			if (inputScalars.useImages && !inputScalars.largeDim) {
				d_inputI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
				OCL_CHECK(status, "Failed to create prior image\n", -1);
			}
		}
		// RDP reference image
		if (MethodList.RDP && w_vec.RDPLargeNeighbor && w_vec.RDP_anatomical) {
			if (inputScalars.useImages) {
				d_RDPrefI = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, region[0], region[1], region[2], 0, 0, NULL, &status);
				OCL_CHECK(status, "Failed to create RDP reference image\n", -1);
			}
		}
#endif // END CUDA
		// Create the necessary buffers
		// Distance-based weighting for GGMRF, RDP and hyperbolic prior
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			ALLOC_BUFFER(d_weights, CL_MEM_READ_ONLY, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1);
			CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
			memAlloc.GGMRF = true;
		}
		memAlloc.tSteps = inputScalars.Nt;
		// NLM anatomical reference image
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
			CUDA_CHECK(status, "\n", -1);
			if (inputScalars.useImages)
				memAlloc.NLMRef = 1;
			else
				memAlloc.NLMRef = 2;
#endif // END CUDA
		}
		if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
			if (inputScalars.useBuffers) {
				ALLOC_BUFFER(d_maskPriorB, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ);
			}
			else {
#if defined(CUDA) || defined(HIP)
				std::memset(&texDesc, 0, sizeof(texDesc));
				std::memset(&resDesc, 0, sizeof(resDesc));
				std::memset(&viewDesc, 0, sizeof(viewDesc));
				// Mask images
				if (inputScalars.maskBPZ > 1) {
					std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
					arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
					arr3DDesc.NumChannels = 1;
					arr3DDesc.Height = inputScalars.Nx[0];
					arr3DDesc.Width = inputScalars.Ny[0];
					arr3DDesc.Depth = inputScalars.Nz[0];
					status = cuArray3DCreate(&maskArrayPrior, &arr3DDesc);
					CUDA_CHECK(status, "\n", -1);
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
					CUDA_CHECK(status, "\n", -1);
				}
				else {
					std::memset(&arr2DDesc, 0, sizeof(arr2DDesc));
					arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
					arr2DDesc.NumChannels = 1;
					arr2DDesc.Height = inputScalars.Nx[0];
					arr2DDesc.Width = inputScalars.Ny[0];
					status = cuArrayCreate(&maskArrayPrior, &arr2DDesc);
					CUDA_CHECK(status, "\n", -1);
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
					CUDA_CHECK(status, "\n", -1);
				}
				resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
				resDesc.res.array.hArray = maskArrayPrior;
				texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
				texDesc.flags = CU_TRSF_READ_AS_INTEGER;
				viewDesc.height = inputScalars.Nx[0];
				viewDesc.width = inputScalars.Ny[0];
				if (inputScalars.maskBPZ > 1)
					viewDesc.depth = inputScalars.maskBPZ;
				viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
				status = cuTexObjectCreate(&d_maskPrior, &resDesc, &texDesc, &viewDesc);
				memAlloc.priorMask = true;
#else
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
#endif // END CUDA
			}
			CHECK(status, "\n", -1);
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				ALLOC_BUFFER(d_V, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_V);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.V = true;
#endif // END CUDA
			}
			// Detector coordinates
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				ALLOC_BUFFER(d_x[0][0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.xSteps++;
#endif // END CUDA
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				ALLOC_BUFFER(d_xFull[0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.xFull = true;
#endif // END CUDA
			}
			// Mask images
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						d_maskFPB.resize(inputScalars.subsetsUsed);
						for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
							ALLOC_BUFFER(d_maskFPB[kk], CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk]);
					}
					else {
#if defined(CUDA) || defined(HIP)
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
								CUDA_CHECK(status, "\n", -1);
								cpy3d.Depth = length[kk];
								cpy3d.srcHost = &w_vec.maskFP[pituus[kk] * vecSize];
								cpy3d.dstArray = maskArrayFP[kk];
								status = cuMemcpy3D(&cpy3d);
								CUDA_CHECK(status, "\n", -1);
							}
						}
						else {
							maskArrayFP.resize(1);
							arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr2DDesc.NumChannels = 1;
							arr2DDesc.Height = inputScalars.nRowsD;
							arr2DDesc.Width = inputScalars.nColsD;
							status = cuArrayCreate(&maskArrayFP[0], &arr2DDesc);
							CUDA_CHECK(status, "\n", -1);
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
							CUDA_CHECK(status, "\n", -1);
						}
						resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						texDesc.flags = CU_TRSF_READ_AS_INTEGER;
						viewDesc.height = inputScalars.nRowsD;
						viewDesc.width = inputScalars.nColsD;
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
						CUDA_CHECK(status, "\n", -1);
						memAlloc.maskFP = true;
#else
						imX = inputScalars.nRowsD;
						imY = inputScalars.nColsD;
						imZ = inputScalars.maskFPZ;
						if (imZ > 1) {
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
								d_maskFP3.emplace_back(cl::Image3D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, length[kk], 0, 0, NULL, &status));
						}
						else
							d_maskFP = cl::Image2D(CLContext, CL_MEM_READ_ONLY, formatMask, imX, imY, 0, NULL, &status);
#endif
					}
#if defined(CUDA) || defined(HIP)
					memAlloc.maskFP = true;
#endif
				}
				if (inputScalars.maskBP) {
					if (inputScalars.useBuffers)
						ALLOC_BUFFER(d_maskBPB, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ);
					else {
#if defined(CUDA) || defined(HIP)
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
							CUDA_CHECK(status, "\n", -1);
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
							CUDA_CHECK(status, "\n", -1);
						}
						else {
							arr2DDesc.Format = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
							arr2DDesc.NumChannels = 1;
							arr2DDesc.Height = inputScalars.Nx[0];
							arr2DDesc.Width = inputScalars.Ny[0];
							status = cuArrayCreate(&maskArrayBP, &arr2DDesc);
							CUDA_CHECK(status, "\n", -1);
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
							CUDA_CHECK(status, "\n", -1);
						}
						resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
						resDesc.res.array.hArray = maskArrayBP;
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						texDesc.flags = CU_TRSF_READ_AS_INTEGER;
						viewDesc.height = inputScalars.Nx[0];
						viewDesc.width = inputScalars.Ny[0];
						viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
						if (inputScalars.maskBPZ > 1) {
							viewDesc.depth = inputScalars.maskBPZ;
							texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						}
						if (inputScalars.BPType == 4 && !inputScalars.CT) {
							texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
						}
						status = cuTexObjectCreate(&d_maskBP, &resDesc, &texDesc, &viewDesc);
						CUDA_CHECK(status, "\n", -1);
#else
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
#endif
					}
#if defined(CUDA) || defined(HIP)
					memAlloc.maskBP = true;
#endif
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				ALLOC_BUFFER(d_zFull[0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.zFull = true;
#endif
			}
			if (inputScalars.SPECT) {
				ALLOC_BUFFER(d_rayShiftsDetector, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections);
				CHECK(status, "\n", -1);
				ALLOC_BUFFER(d_rayShiftsSource, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.rayShifts = true;
#endif // END CUDA
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution) {
				ALLOC_BUFFER(d_eFOVIndices, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nz[0]);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.eFOV = true;
#endif // END CUDA
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				ALLOC_BUFFER(d_angle, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.angle = true;
#endif // END CUDA
			}
			// TOF bin centers
			if (inputScalars.TOF) {
				ALLOC_BUFFER(d_TOFCenter, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nBins);
				CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
				memAlloc.TOF = true;
#endif // END CUDA
			}
			for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
				if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
					if (inputScalars.size_atten > inputScalars.im_dim[0] || timestep == 0) {
						if (inputScalars.useBuffers)
							ALLOC_BUFFER(d_attenB[timestep], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0]);
						else {
#if defined(CUDA) || defined(HIP)
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
							CUDA_CHECK(status, "\n", -1);
							CUDA_MEMCPY3D cpy3d;
							std::memset(&cpy3d, 0, sizeof(cpy3d));
							cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
							cpy3d.srcHost = &atten[inputScalars.im_dim[0] * timestep];
							cpy3d.srcPitch = inputScalars.Ny[0] * sizeof(float);
							cpy3d.srcHeight = inputScalars.Nx[0];
							cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
							cpy3d.dstArray = atArray;
							cpy3d.WidthInBytes = inputScalars.Ny[0] * sizeof(float);
							cpy3d.Height = inputScalars.Nx[0];
							cpy3d.Depth = inputScalars.Nz[0];
							status = cuMemcpy3D(&cpy3d);
							CUDA_CHECK(status, "\n", -1);
							resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
							resDesc.res.array.hArray = atArray;
							texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
								texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
								texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
							}
						else
							texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
						viewDesc.height = inputScalars.Nx[0];
						viewDesc.width = inputScalars.Ny[0];
						viewDesc.depth = inputScalars.Nz[0];
						viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
						status = cuTexObjectCreate(&d_attenIm[timestep], &resDesc, &texDesc, &viewDesc);
						CUDA_CHECK(status, "\n", -1);
#else
						imZ = inputScalars.Nz[0];
						d_attenIm[timestep] = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
						OCL_CHECK(status, "\n", -1);
#endif // END CUDA
						}
#if defined(CUDA) || defined(HIP)
					memAlloc.atten = true;
					memAlloc.attenSize++;
#endif // END CUDA
					}
				}
				for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
					if (inputScalars.CT || inputScalars.SPECT) {
						ALLOC_BUFFER(d_x[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 6);
#if defined(CUDA) || defined(HIP)
						CUDA_CHECK(status, "\n", -1);
						memAlloc.xSteps++;
#endif // END CUDA
					}
					// First condition: load all data at once. Second condition: load one subset at a time (only 1 buffer required, loadCoord reloads it for each subset/timestep)
					else if (inputScalars.listmode > 0 && !inputScalars.indexBased && (inputScalars.loadTOF || (kk == 0 && timestep == 0))) {
						ALLOC_BUFFER(d_x[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk + timestep * inputScalars.subsets] * 6);
						CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
						memAlloc.xSteps++;
#endif // END CUDA
						if (DEBUG) {
							mexPrintBase("length[kk + timestep * inputScalars.subsets] * 6 = %u\n", length[kk + timestep * inputScalars.subsets] * 6);
							mexEval();
						}
					}
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode != 1) {
						size_t coef = 2;
						if (inputScalars.useHelical)
							coef = 1;
						else if (inputScalars.pitch)
							coef = 6;
						ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * coef);
						CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
						memAlloc.zType = 1;
						memAlloc.zSteps++;
#endif // END CUDA
					}
					else {
						if (inputScalars.PET && inputScalars.listmode == 0) {
							if (inputScalars.nLayers > 1)
								ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 3);
							else
								ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2);
#if defined(CUDA) || defined(HIP)
							memAlloc.zType = 1;
							memAlloc.zSteps++;
#endif // END CUDA
						}
						else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
							ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z);
#if defined(CUDA) || defined(HIP)
							memAlloc.zType = 0;
							memAlloc.zSteps = kk;
#endif // END CUDA
						}
						else {
							ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z);
#if defined(CUDA) || defined(HIP)
							memAlloc.zType = 0;
							memAlloc.zSteps = kk;
#endif // END CUDA
						}
						CHECK(status, "\n", -1);
					}
					if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) { // Scatter correction buffer
						ALLOC_BUFFER(d_scat[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize);
						CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
						memAlloc.extra = true;
						memAlloc.eSteps++;
#endif // END CUDA
					}
					if (inputScalars.listmode > 0 && inputScalars.indexBased) {
						if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF && timestep == 0)) {
							ALLOC_BUFFER(d_trIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2);
							CHECK(status, "\n", -1);
							ALLOC_BUFFER(d_axIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2);
							CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
							memAlloc.indexBased = true;
							memAlloc.iSteps++;
#endif // END CUDA
						}
					}
					if (inputScalars.listmode > 0 && inputScalars.TOF) {
						if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF && timestep == 0)) {
							ALLOC_BUFFER(d_TOFIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint8_t) * length[kk + timestep * inputScalars.subsets]);
							CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
							memAlloc.TOFIndex = true;
							memAlloc.TOFSteps++;
#endif // END CUDA
						}
					}
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				// Redundancy weighting
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					ALLOC_BUFFER(d_T[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk]);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.offsetT = true;
					memAlloc.oSteps++;
#endif // END CUDA
				}
				if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
					ALLOC_BUFFER(d_geomProj5[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 16);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.geom5 = true;
					memAlloc.g5Steps++;
#endif // END CUDA
				}
				// Normalization weighting (OpenCL condition also requires size_norm > 1; used for both backends)
				if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
					ALLOC_BUFFER(d_norm[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.norm = true;
					memAlloc.nSteps++;
#endif // END CUDA
				}
				// Measurement-based attenuation correction
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					ALLOC_BUFFER(d_atten[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * vecSize);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.attenM = true;
					memAlloc.aSteps++;
#endif // END CUDA
				}
				// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
				// Note that raw data format is not used at the moment
				if (inputScalars.raw && inputScalars.listmode != 1) {
					ALLOC_BUFFER(d_L[kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.raw = true;
					memAlloc.lSteps++;
#endif // END CUDA
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					ALLOC_BUFFER(d_xyindex[kk], CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk]);
					CHECK(status, "\n", -1);
					ALLOC_BUFFER(d_zindex[kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk]);
					CHECK(status, "\n", -1);
#if defined(CUDA) || defined(HIP)
					memAlloc.subInd = true;
					memAlloc.iSteps++;
#endif // END CUDA
				}
			}
		}
		CHECK(status, "Buffer creation failed\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer creation succeeded\n");
		}

		// assign values to the buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			WRITE_BUFFER(d_weights, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, w_vec.weights);
			memSize += (sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1)) / 1048576ULL;
			CHECK(status, "\n", -1);
		}
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (!inputScalars.useImages)
				WRITE_BUFFER(d_uref, sizeof(float) * inputScalars.im_dim[0], w_vec.NLM_ref);
#if !defined(CUDA) && !defined(HIP)
			else {
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
			CHECK(status, "\n", -1);
#endif // END CUDA
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				WRITE_BUFFER(d_V, sizeof(float) * inputScalars.size_V, inputScalars.V);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * inputScalars.size_V) / 1048576ULL;
			}
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				WRITE_BUFFER(d_x[0][0], sizeof(float) * inputScalars.size_of_x, x);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * inputScalars.size_of_x) / 1048576ULL;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				WRITE_BUFFER(d_xFull[0], sizeof(float) * inputScalars.size_of_x, x);
				CHECK(status, "\n", -1);
				WRITE_BUFFER(d_zFull[0], sizeof(float) * inputScalars.size_z, z_det);
				CHECK(status, "\n", -1);
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution) {
				WRITE_BUFFER(d_eFOVIndices, sizeof(uint8_t) * inputScalars.Nz[0], w_vec.eFOVIndices);
				CHECK(status, "\n", -1);
				memSize += (sizeof(uint8_t) * inputScalars.Nz[0]) / 1048576ULL;
			}
			if (inputScalars.maskFP || inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
				if (inputScalars.useBuffers) {
					if (inputScalars.maskFP) {
						if (inputScalars.maskFPZ > 1)
							for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
								WRITE_BUFFER(d_maskFPB[kk], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[kk], &w_vec.maskFP[pituus[kk] * vecSize]);
						else
							WRITE_BUFFER(d_maskFPB[0], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD, w_vec.maskFP);
						CHECK(status, "\n", -1);
					}
					else if (inputScalars.maskBP) {
						WRITE_BUFFER(d_maskBPB, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ, w_vec.maskBP);
						CHECK(status, "\n", -1);
					}
					if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
						WRITE_BUFFER(d_maskPriorB, sizeof(uint8_t) * inputScalars.Nx[0] * inputScalars.Ny[0] * inputScalars.maskBPZ, w_vec.maskPrior);
						CHECK(status, "\n", -1);
					}
				}
#if !defined(CUDA) && !defined(HIP)
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
						OCL_CHECK(status, "\n", -1);
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
						OCL_CHECK(status, "\n", -1);
						memSize += (sizeof(bool) * inputScalars.Nx[0] * inputScalars.Ny[0]) / 1048576ULL;
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
						OCL_CHECK(status, "\n", -1);
						memSize += (sizeof(bool) * inputScalars.Nx[0] * inputScalars.Ny[0]) / 1048576ULL;
					}
				}
#endif // END CUDA
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				WRITE_BUFFER(d_angle, sizeof(float) * inputScalars.nProjections, w_vec.angles);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * inputScalars.nProjections) / 1048576ULL;
			}
			if (inputScalars.TOF) {
				WRITE_BUFFER(d_TOFCenter, sizeof(float) * inputScalars.nBins, inputScalars.TOFCenter);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * inputScalars.nBins) / 1048576ULL;
			}
			if (inputScalars.SPECT) {
				WRITE_BUFFER(d_rayShiftsDetector, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections, w_vec.rayShiftsDetector);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections) / 1048576ULL;
				WRITE_BUFFER(d_rayShiftsSource, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections, w_vec.rayShiftsSource);
				CHECK(status, "\n", -1);
				memSize += (sizeof(float) * 2 * inputScalars.n_rays * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nProjections) / 1048576ULL;
			}

			if (DEBUG) {
				mexPrint("Timestep phase\n");
			}
			for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
				if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
					if (inputScalars.useBuffers) {
						if (inputScalars.size_atten > inputScalars.im_dim[0]) {
							WRITE_BUFFER(d_attenB[timestep], sizeof(float) * inputScalars.im_dim[0], &atten[inputScalars.im_dim[0] * timestep]);
							CHECK(status, "\n", -1);
						}
						else if (timestep == 0)
							WRITE_BUFFER(d_attenB[timestep], sizeof(float) * inputScalars.im_dim[0], atten);
					}
#if !defined(CUDA) && !defined(HIP)
					else {
						cl::detail::size_t_array region = { { 0, 0, 0 } };
						region[0] = inputScalars.Nx[0];
						region[1] = inputScalars.Ny[0];
						region[2] = inputScalars.Nz[0];
						if (inputScalars.size_atten > inputScalars.im_dim[0]) {
							status = CLCommandQueue[0].enqueueWriteImage(d_attenIm[timestep], CL_FALSE, origin, region, 0, 0, &atten[inputScalars.im_dim[0] * timestep]);
						}
						else if (timestep == 0) {
							status = CLCommandQueue[0].enqueueWriteImage(d_attenIm[timestep], CL_FALSE, origin, region, 0, 0, atten);
						}
						OCL_CHECK(status, "\n", -1);
					}
#endif // END CUDA
					memSize += (sizeof(float) * inputScalars.im_dim[0]) / 1048576ULL;
				}
				for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
						size_t kerroin = 2;
						if (inputScalars.pitch)
							kerroin = 6;
						else if (inputScalars.useHelical)
							kerroin = 1;
						WRITE_BUFFER(d_z[timestep][kk], sizeof(float) * length[kk] * kerroin, &z_det[pituus[kk] * kerroin]);
						CHECK(status, "\n", -1);
						memSize += (sizeof(float) * length[kk] * kerroin) / 1048576ULL;
					}
					else {
						if (inputScalars.PET && inputScalars.listmode == 0) {
							int64_t kerroin = 2;
							if (inputScalars.nLayers > 1)
								kerroin = 3;
							WRITE_BUFFER(d_z[timestep][kk], sizeof(float) * length[kk] * kerroin, &z_det[pituus[kk] * kerroin]);
							memSize += (sizeof(float) * length[kk] * kerroin) / 1048576ULL;
						}
						else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased || inputScalars.listmode > 0)) {
							WRITE_BUFFER(d_z[timestep][kk], sizeof(float) * inputScalars.size_z, z_det);
							memSize += (sizeof(float) * inputScalars.size_z) / 1048576ULL;
						}
						CHECK(status, "\n", -1);
					}
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
						WRITE_BUFFER(d_x[timestep][kk], sizeof(float) * length[kk] * 6, &x[pituus[kk] * 6]);
						CHECK(status, "\n", -1);
						memSize += (sizeof(float) * length[kk] * 6) / 1048576ULL;
					}
					else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
						if (inputScalars.loadTOF || (kk == 0 && timestep == 0)) {
							WRITE_BUFFER(d_x[timestep][kk], sizeof(float) * length[kk + timestep * inputScalars.subsets] * 6, &w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6]);
							CHECK(status, "\n", -1);
							if (DEBUG) {
								mexPrintBase("length[kk + timestep * inputScalars.subsets] * 6 = %u\n", length[kk + timestep * inputScalars.subsets] * 6);
								mexPrintBase("pituus[kk + timestep * inputScalars.subsets] * 6 = %u\n", pituus[kk + timestep * inputScalars.subsets] * 6);
								mexPrintBase("w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6] = %f\n", w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6]);
								mexEval();
							}
							memSize += (sizeof(float) * length[kk] * 6) / 1048576ULL;
						}
					}
					if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U) { // Load scatter data
						WRITE_BUFFER(d_scat[timestep][kk], sizeof(float) * length[kk] * vecSize, &extraCorr[pituus[kk] * vecSize + inputScalars.kokoNonTOF * timestep]);
						CHECK(status, "\n", -1);
						memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
					}
					if (inputScalars.listmode > 0 && inputScalars.indexBased) {
						if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF && timestep == 0)) { // First condition: load all data at once. Second condition: load one subset at a time (only 1 buffer required for each timestep).
							WRITE_BUFFER(d_trIndex[timestep][kk], sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2, &w_vec.trIndex[pituus[kk + timestep * inputScalars.subsets] * 2]);
							CHECK(status, "\n", -1);
							memSize += (sizeof(uint16_t) * length[kk] * 2) / 1048576ULL;
							WRITE_BUFFER(d_axIndex[timestep][kk], sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2, &w_vec.axIndex[pituus[kk + timestep * inputScalars.subsets] * 2]);
							CHECK(status, "\n", -1);
							memSize += (sizeof(uint16_t) * length[kk] * 2) / 1048576ULL;
						}
					}
					if (inputScalars.listmode > 0 && inputScalars.TOF) {
						if (inputScalars.loadTOF || (kk == 0 && !inputScalars.loadTOF && timestep == 0)) {
							WRITE_BUFFER(d_TOFIndex[timestep][kk], sizeof(uint8_t) * length[kk + timestep * inputScalars.subsets], &w_vec.TOFIndices[pituus[kk + timestep * inputScalars.subsets]]);
							CHECK(status, "\n", -1);
							memSize += (sizeof(uint8_t) * length[kk]) / 1048576ULL;
						}
					}
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					WRITE_BUFFER(d_T[kk], sizeof(float) * length[kk], &inputScalars.T[pituus[kk]]);
					CHECK(status, "\n", -1);
				}
				// Per projection: s (3), d3 (3), normX (3), normY (3), crossP (3), upperPart (1), i.e. 16 floats.
				// Previously computed in the kernel itself
				if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
					const float indX = static_cast<float>(inputScalars.nRowsD) / 2.f;
					const float indY = static_cast<float>(inputScalars.nColsD) / 2.f;
					const size_t uvStride = inputScalars.pitch ? 6 : 2;
					const float* xs = &x[pituus[kk] * 6];
					const float* uv = &z_det[pituus[kk] * uvStride];
					geomProj5Host[kk].resize(static_cast<size_t>(length[kk]) * 16);
					for (int64_t pp = 0; pp < length[kk]; pp++) {
						const float sX = xs[pp * 6], sY = xs[pp * 6 + 1], sZ = xs[pp * 6 + 2];
						const float dX = xs[pp * 6 + 3], dY = xs[pp * 6 + 4], dZ = xs[pp * 6 + 5];
						float aXx, aXy, aXz, aYx, aYy, aYz;
						if (inputScalars.pitch) {
							aXx = uv[pp * 6] * indX; aXy = uv[pp * 6 + 1] * indX; aXz = uv[pp * 6 + 2] * indX;
							aYx = uv[pp * 6 + 3] * indY; aYy = uv[pp * 6 + 4] * indY; aYz = uv[pp * 6 + 5] * indY;
						}
						else {
							aXx = uv[pp * 2] * indX; aXy = uv[pp * 2 + 1] * indX; aXz = 0.f;
							aYx = 0.f; aYy = 0.f; aYz = indY * w_vec.dPitchY;
						}
						// d3 = d - apuX - apuY
						const float d3x = dX - aXx - aYx, d3y = dY - aXy - aYy, d3z = dZ - aXz - aYz;
						// d2 = apuX - apuY
						const float d2x = aXx - aYx, d2y = aXy - aYy, d2z = aXz - aYz;
						const float nXl = std::sqrt(aXx * aXx + aXy * aXy + aXz * aXz);
						const float nYl = std::sqrt(aYx * aYx + aYy * aYy + aYz * aYz);
						// d3 - d = -apuX - apuY
						const float ddx = d3x - dX, ddy = d3y - dY, ddz = d3z - dZ;
						// crossP = cross(d2, d3 - d)
						const float cx = d2y * ddz - d2z * ddy;
						const float cy = d2z * ddx - d2x * ddz;
						const float cz = d2x * ddy - d2y * ddx;
						// upperPart = dot(crossP, s - d)
						const float up = cx * (sX - dX) + cy * (sY - dY) + cz * (sZ - dZ);
						float* g = &geomProj5Host[kk][pp * 16];
						g[0] = sX; g[1] = sY; g[2] = sZ;
						g[3] = d3x; g[4] = d3y; g[5] = d3z;
						g[6] = aXx / nXl; g[7] = aXy / nXl; g[8] = aXz / nXl;
						g[9] = aYx / nYl; g[10] = aYy / nYl; g[11] = aYz / nYl;
						g[12] = cx; g[13] = cy; g[14] = cz;
						g[15] = up;
					}
					WRITE_BUFFER(d_geomProj5[kk], sizeof(float) * length[kk] * 16, geomProj5Host[kk].data());
					CHECK(status, "\n", -1);
					memSize += (sizeof(float) * length[kk] * 16) / 1048576ULL;
				}
				if (inputScalars.raw && inputScalars.listmode != 1) {
					WRITE_BUFFER(d_L[kk], sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2]);
					CHECK(status, "\n", -1);
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					WRITE_BUFFER(d_zindex[kk], sizeof(uint16_t) * length[kk], &z_index[pituus[kk]]);
					CHECK(status, "\n", -1);
					WRITE_BUFFER(d_xyindex[kk], sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]]);
					CHECK(status, "\n", -1);
					memSize += (sizeof(uint32_t) * length[kk] + sizeof(uint16_t) * length[kk]) / 1048576ULL;
				}
				if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction) {
					WRITE_BUFFER(d_norm[kk], sizeof(float) * length[kk] * vecSize, &norm[pituus[kk] * vecSize]);
					CHECK(status, "\n", -1);
					memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
				}
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					WRITE_BUFFER(d_atten[kk], sizeof(float) * length[kk] * vecSize, &atten[pituus[kk] * vecSize]);
					CHECK(status, "\n", -1);
					memSize += (sizeof(float) * length[kk] * vecSize) / 1048576ULL;
				}
			}
#if defined(CUDA) || defined(HIP)
			status = cuCtxSynchronize();
#else
			status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		}
		CHECK(status, "Buffer write failed\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
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
	inline int createBuffers(scalarStruct & inputScalars, Weighting & w_vec, const float* x, const float* z_det, const uint32_t * xy_index,
		const uint16_t * z_index, const uint16_t * L, const int64_t * pituus, const float* atten, const float* norm, const float* extraCorr,
		const std::vector<int64_t>&length, const RecMethods & MethodList, const int type = 0) {
#if defined(CUDA) || defined(HIP)
		int status = 0;
#else
		cl_int status = CL_SUCCESS;
#endif // END CUDA
		if (inputScalars.raw)
			d_L.resize(inputScalars.subsetsUsed);
		if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1) {
			d_xyindex.resize(inputScalars.subsetsUsed);
			d_zindex.resize(inputScalars.subsetsUsed);
		}
		if (inputScalars.normalization_correction)
			d_norm.resize(inputScalars.subsetsUsed);
		if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
			d_atten.resize(inputScalars.subsetsUsed);
		if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
			d_attenB.resize(inputScalars.Nt);
			d_attenIm.resize(inputScalars.Nt);
		}
		if (inputScalars.projector_type != 6) {
			d_scat.resize(inputScalars.Nt);
			d_x.resize(inputScalars.Nt);
			d_z.resize(inputScalars.Nt);
			d_trIndex.resize(inputScalars.Nt);
			d_axIndex.resize(inputScalars.Nt);
			d_TOFIndex.resize(inputScalars.Nt);
			for (int tt = 0; tt < inputScalars.Nt; tt++) {
				d_scat[tt].resize(inputScalars.subsetsUsed);
				d_x[tt].resize(inputScalars.subsetsUsed);
				d_z[tt].resize(inputScalars.subsetsUsed);
				d_trIndex[tt].resize(inputScalars.subsetsUsed);
				d_axIndex[tt].resize(inputScalars.subsetsUsed);
				d_TOFIndex[tt].resize(inputScalars.subsetsUsed);
			}
		}
		if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
			d_T.resize(inputScalars.subsetsUsed);
		if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
			d_geomProj5.resize(inputScalars.subsetsUsed);
			geomProj5Host.resize(inputScalars.subsetsUsed);
		}

		status = createAndWriteBuffers(length, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, inputScalars, w_vec, MethodList);
#if defined(CUDA) || defined(HIP)
		if (status != 0) {
#else
		if (status != CL_SUCCESS) {
#endif // END CUDA
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
	inline int initializeKernel(scalarStruct & inputScalars, Weighting & w_vec) {
#if defined(CUDA) || defined(HIP)
		int status = 0;
#else
		cl_int status = CL_SUCCESS;
#endif // END CUDA

		if (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.FPType == 7) {
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nRowsD);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nColsD);
			KARG(FPArgs, kernelFP, kernelIndFP, dPitch);
			if (inputScalars.useHelical) {
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.helicalRadius);
			}
		}

		if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || inputScalars.BPType == 7) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nRowsD);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nColsD);
			KARG(BPArgs, kernelBP, kernelIndBP, dPitch);
			if (inputScalars.useHelical)
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.helicalRadius);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nRowsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nColsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, dPitch);
			}
		}
		if (inputScalars.FPType == 4) {
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.dL);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.global_factor);
		}
		if (inputScalars.BPType == 4 && !inputScalars.CT) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.dL);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.global_factor);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.dL);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.global_factor);
			}
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			// Set the kernelFP parameters that do not change
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.global_factor);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.epps);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nRowsD);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.det_per_ring);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.sigma_x);

			if (inputScalars.SPECT) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_rayShiftsDetector);
				KARG(FPArgs, kernelFP, kernelIndFP, d_rayShiftsSource);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffA);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffB);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffC);
				KARG(FPArgs, kernelFP, kernelIndFP, totalFOVmin);
				KARG(FPArgs, kernelFP, kernelIndFP, totalFOVmax);
			}

			KARG(FPArgs, kernelFP, kernelIndFP, dPitch);
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (inputScalars.FPType == 2) {
					KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.tube_width);
				}
				else {
					KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.cylRadiusProj3);
				}
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.bmin);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.bmax);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.Vmax);
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.global_factor);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.epps);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nRowsD);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.det_per_ring);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.sigma_x);

			if (inputScalars.SPECT) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_rayShiftsDetector);
				KARG(BPArgs, kernelBP, kernelIndBP, d_rayShiftsSource);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffA);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffB);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffC);
				KARG(BPArgs, kernelBP, kernelIndBP, totalFOVmin);
				KARG(BPArgs, kernelBP, kernelIndBP, totalFOVmax);
			}

			KARG(BPArgs, kernelBP, kernelIndBP, dPitch);
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (inputScalars.BPType == 2) {
					KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.tube_width);
				}
				else {
					KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.cylRadiusProj3);
				}
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.bmin);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.bmax);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.Vmax);
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.global_factor);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.epps);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nRowsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.det_per_ring);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.sigma_x);
				KARG(SensArgs, kernelSensList, kernelIndSens, dPitch);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					if (inputScalars.BPType == 2) {
						KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.tube_width);
					}
					else
						KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.cylRadiusProj3);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.bmin);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.bmax);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.Vmax);
				}

			}
		}
		if (DEBUG) {
			mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
			mexPrintBase("inputScalars.helicalRadius = %f\n", inputScalars.helicalRadius);
			mexEval();
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if (inputScalars.TOF) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_TOFCenter);
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_V);
			}
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nColsD);
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if (inputScalars.TOF)
				KARG(BPArgs, kernelBP, kernelIndBP, d_TOFCenter);
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_V);
			}
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nColsD);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.TOF)
					KARG(SensArgs, kernelSensList, kernelIndSens, d_TOFCenter);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					KARG(SensArgs, kernelSensList, kernelIndSens, d_V);
				}
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nColsD);
			}
		}
		if ((inputScalars.BPType == 4 || inputScalars.FPType == 4) && !inputScalars.CT && inputScalars.TOF) {
			if (inputScalars.FPType == 4) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_TOFCenter);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.sigma_x);
			}
			if (inputScalars.BPType == 4) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_TOFCenter);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.sigma_x);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					KARG(SensArgs, kernelSensList, kernelIndSens, d_TOFCenter);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.sigma_x);
				}
			}
		}
		if (DEBUG) {
#if defined(CUDA) || defined(HIP)
			mexPrintBase("kernelIndFP = %u\n", FPArgs.size());
			mexPrintBase("kernelIndBP = %u\n", BPArgs.size());
#else
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("kernelIndBP = %u\n", kernelIndBP);
#endif // END CUDA
			mexEval();
		}
		return 0;
	}

	template <typename T>
	inline int loadCoord(uint32_t currentSubset, scalarStruct & inputScalars, const int64_t length, const T * listCoord, const T * listCoordAx = nullptr, const uint8_t * TOFIndices = nullptr) {
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.indexBased) {
#if defined(CUDA) || defined(HIP)
			// CUDA/HIP must explicitly release the previous allocation; OpenCL's cl::Buffer reassignment does this itself.
			getErrorString(cuMemFree(d_trIndex[0][0]));
			getErrorString(cuMemFree(d_axIndex[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_trIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2);
			CHECK(status, "\n", -1);
			ALLOC_BUFFER(d_axIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_trIndex[0][0], sizeof(uint16_t) * length * 2, listCoord);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_axIndex[0][0], sizeof(uint16_t) * length * 2, listCoordAx);
			CHECK(status, "\n", -1);
		}
		else {
#if defined(CUDA) || defined(HIP)
			getErrorString(cuMemFree(d_x[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_x[0][0], CL_MEM_READ_ONLY, sizeof(float) * length * 6);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_x[0][0], sizeof(float) * length * 6, listCoord);
			CHECK(status, "\n", -1);
		}
		if (inputScalars.TOF) {
#if defined(CUDA) || defined(HIP)
			getErrorString(cuMemFree(d_TOFIndex[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_TOFIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint8_t) * length);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_TOFIndex[0][0], sizeof(uint8_t) * length, TOFIndices);
			CHECK(status, "\n", -1);
		}
		return 0;
#if !defined(CUDA) && !defined(HIP)
	}

	int computeConvolutionF(const scalarStruct & inputScalars, const int ii = 0) {
		int status = CL_SUCCESS;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		cl_uint kernelInd = 0U;
		KARG(kArgs, kernelPSFf, kernelInd, d_imFinal[ii]);
		KARG(kArgs, kernelPSFf, kernelInd, d_imTemp[ii]);
		KARG(kArgs, kernelPSFf, kernelInd, d_g);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_x);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_y);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSFf, cl::NDRange(), globalC, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");

		return status;
	}

	template <typename T>
	int computeConvolution(const scalarStruct & inputScalars, cl::Buffer & input, const int ii = 0, const T & cType = 0) {
		int status = CL_SUCCESS;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		cl::Buffer d_BPApu = cl::Buffer(CLContext, CL_MEM_READ_WRITE, sizeof(T) * inputScalars.im_dim[ii], NULL, &status);
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		cl_uint kernelInd = 0U;
		KARG(kArgs, kernelPSF, kernelInd, input);
		KARG(kArgs, kernelPSF, kernelInd, d_BPApu);
		KARG(kArgs, kernelPSF, kernelInd, d_g);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_x);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_y);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSF, cl::NDRange(), globalC, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");
		status = CLCommandQueue[0].enqueueCopyBuffer(d_BPApu, input, 0, 0, sizeof(T) * inputScalars.im_dim[ii]);
		OCL_CHECK(status, "\n", -1);
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);

		return status;
	}

	int computeForward(const scalarStruct & inputScalars, const std::vector<int64_t>&length, const uint32_t osa_iter, const uint32_t timestep = 0) {
		int status = CL_SUCCESS;
		cl::NDRange localF = { 64, 1, 1 };
		int indD = osa_iter + timestep * inputScalars.subsets;
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { inputScalars.nRowsD * inputScalars.nColsD * static_cast<size_t>(length[indD]) * inputScalars.nBins, 1, 1 };
		else
			if (inputScalars.listmode == 0)
				global = { static_cast<cl::size_type>(length[indD]) * inputScalars.nBins, 1, 1 };
			else
				global = { static_cast<cl::size_type>(length[indD]), 1, 1 };

		size_t erotusF = global[0] % localF[0];
		if (erotusF > 0)
			erotusF = localF[0] - erotusF;
		global = { static_cast<cl::size_type>(global[0] + erotusF), 1, 1 };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("global.dimensions() = %u\n", global.dimensions());
			mexPrintBase("length[indD] = %u\n", length[indD]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexEval();
		}
		cl_uint kernelInd = 0U;
		KARG(kArgs, kernelForward, kernelInd, d_output);
		KARG(kArgs, kernelForward, kernelInd, d_meas[osa_iter]);
		if (inputScalars.CT)
			KARG(kArgs, kernelForward, kernelInd, d_outputCT);
		if (inputScalars.randoms_correction)
			KARG(kArgs, kernelForward, kernelInd, d_rand[osa_iter]);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelForward, cl::NDRange(), global, localF, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward step kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		if (inputScalars.verbose >= 3)
			mexPrint("Forward step computed");

		return status;
	}

	int computeEstimate(const scalarStruct & inputScalars, const int ii = 0, const int uu = 0, const int timestep = 0) {
		int status = CL_SUCCESS;
		global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		cl_uint kernelInd = 0U;
		KARG(kArgs, kernelEstimate, kernelInd, d_Summ[uu]);
		KARG(kArgs, kernelEstimate, kernelInd, vec_opencl.d_rhs_os[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, d_imFinal[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, inputScalars.epps);
		KARG(kArgs, kernelEstimate, kernelInd, d_N[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, no_norm);
		if (inputScalars.CT)
			KARG(kArgs, kernelEstimate, kernelInd, inputScalars.flat);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelEstimate, cl::NDRange(), global, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Estimate step kernel launched successfully\n");
		}
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
		if (inputScalars.verbose >= 3)
			mexPrint("Estimate step computed");

		return status;
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int forwardProjection(scalarStruct & inputScalars, Weighting & w_vec, uint32_t osa_iter, uint32_t timestep, const std::vector<int64_t>&length1, uint64_t m_size, int ii = 0, const int uu = 0) {
#else
	inline int forwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const uint32_t timestep, const std::vector<int64_t>&length, uint64_t m_size, const int32_t ii = 0, const int uu = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
#if defined(CUDA) || defined(HIP)
		std::vector<int64_t> length = length1;
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kTemp = FPArgs;
		if (inputScalars.FPType == 5) {
			global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
			global[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1]) / local[1];
			global[2] = length[osa_iter + timestep * inputScalars.subsets];
		}
		else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
			global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
			global[2] = length[osa_iter + timestep * inputScalars.subsets];
		}
#else
		kernelIndFPSubIter = kernelIndFP;
		cl_int status = CL_SUCCESS;
		if (inputScalars.FPType == 5)
			global = { inputScalars.nRowsD + erotus[0], (inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1], static_cast<size_t>(length[osa_iter + timestep * inputScalars.subsets]) };
		else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[osa_iter + timestep * inputScalars.subsets]) };
#endif // END CUDA
		else {
			erotus[0] = length[osa_iter + timestep * inputScalars.subsets] % local_size[0];

			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
#if defined(CUDA) || defined(HIP)
			global[0] = (length[osa_iter + timestep * inputScalars.subsets] + erotus[0]) / local[0];
			global[1] = 1;
			global[2] = 1;
#else
			global = { static_cast<cl::size_type>(length[osa_iter + timestep * inputScalars.subsets] + erotus[0]), 1, 1 };
		}
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}

		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("local[0] = %u\n", local[0]);
			mexPrintBase("local[1] = %u\n", local[1]);
			mexPrintBase("local[2] = %u\n", local[2]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d[ii].s0 = %f\n", d[ii].x);
			mexPrintBase("d[ii].s1 = %f\n", d[ii].y);
			mexPrintBase("d[ii].s2 = %f\n", d[ii].z);
			mexPrintBase("b[ii].s0 = %f\n", b[ii].x);
			mexPrintBase("b[ii].s1 = %f\n", b[ii].y);
			mexPrintBase("b[ii].s2 = %f\n", b[ii].z);
			mexPrintBase("bmax[ii].s0 = %f\n", bmax[ii].x);
			mexPrintBase("bmax[ii].s1 = %f\n", bmax[ii].y);
			mexPrintBase("bmax[ii].s2 = %f\n", bmax[ii].z);
#else
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
#endif // END CUDA
			mexPrintBase("kernelIndFPSubIter = %u\n", kernelIndFPSubIter);
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
			mexPrintBase("size_y = %u\n", inputScalars.nColsD);
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("FPType = %u\n", inputScalars.FPType);
			if (inputScalars.FPType == 4) {
				mexPrintBase("dL = %f\n", inputScalars.dL);
#if defined(CUDA) || defined(HIP)
				mexPrintBase("d_Scale4[ii].s[0] = %f\n", inputScalars.d_Scale4[ii].x);
				mexPrintBase("d_Scale4[ii].s[1] = %f\n", inputScalars.d_Scale4[ii].y);
				mexPrintBase("d_Scale4[ii].s[2] = %f\n", inputScalars.d_Scale4[ii].z);
#else
				mexPrintBase("d_Scale4[ii].s[0] = %f\n", inputScalars.d_Scale4[ii].s[0]);
				mexPrintBase("d_Scale4[ii].s[1] = %f\n", inputScalars.d_Scale4[ii].s[1]);
				mexPrintBase("d_Scale4[ii].s[2] = %f\n", inputScalars.d_Scale4[ii].s[2]);
#endif // END CUDA
			}
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter + timestep * inputScalars.subsets]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexPrintBase("maskBP = %u\n", inputScalars.maskBP);
			mexPrintBase("maskFP = %u\n", inputScalars.maskFP);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("ii = %u\n", ii);
			mexPrintBase("NVOXELS = %u\n", NVOXELS);
			mexPrintBase("NVOXELS5 = %u\n", NVOXELS5);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexPrintBase("timestep = %u\n", timestep);
			mexPrintBase("memSize = %u\n", memSize);
			mexPrintBase("subsetType = %u\n", inputScalars.subsetType);
			mexEval();
		}
#if defined(CUDA) || defined(HIP)
		CUevent tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
		}
#ifndef AF
		// Synchronization is not necessary when using AF, but only when using the custom reconstructions
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "\n", -1);
#endif
#else
#ifndef AF
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
#endif
#endif // END CUDA
		if (!inputScalars.CT && (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 || inputScalars.FPType == 4)) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_atten[osa_iter]);
			}
			else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.size_atten > inputScalars.im_dim[0]) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenB[timestep]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenIm[timestep]);
				}
				else {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenB[0]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenIm[0]);
				}
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.FPType == 4) {
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_N[ii]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, b[ii]);
			if (inputScalars.FPType == 5) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.dSize[ii]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d[ii]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.d_Scale[ii]);
			}
			else {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, bmax[ii]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.d_Scale4[ii]);
			}
		}
		if (inputScalars.FPType == 4) {
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			}
			if ((inputScalars.CT || inputScalars.PET)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			}
			else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[0][0]);
			}
			else {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			}
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP3[osa_iter]);
					}
					else {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
					}
			}
			KARG(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.raw) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_L[osa_iter]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.det_per_ring);
			}
			if (inputScalars.normalization_correction)
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[osa_iter]);
			if (inputScalars.scatter) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_scat[timestep][osa_iter]);
			}
			KARG(kTemp, kernelFP, kernelIndFPSubIter, no_norm);
#if !defined(CUDA) && !defined(HIP)
			m_size = static_cast<cl_ulong>(m_size);
#endif // END CUDA
			KARG(kTemp, kernelFP, kernelIndFPSubIter, m_size);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, osa_iter);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, ii);
		}
		else if (inputScalars.FPType == 5) {
			if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os_int);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			if (inputScalars.meanFP) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_meanFP);
			}
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP3[osa_iter]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
			}
			if (inputScalars.normalization_correction) {
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[osa_iter]);
			}
			KARG(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
		}
		else if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
			if (inputScalars.maskFP) {
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1)
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP3[osa_iter]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
			}
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			}
			else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			if (inputScalars.normalization_correction) {
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[osa_iter]);
			}
			if (inputScalars.scatter)
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_scat[timestep][osa_iter]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_Summ[uu]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_N[ii]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d[ii]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, b[ii]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, bmax[ii]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.raw) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_L[osa_iter]);
			}
			if (inputScalars.useBuffers) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_im);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, no_norm);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, m_size);
			//status = kernelFP.setArg(kernelIndFPSubIter++, static_cast<cl_ulong>(m_size));
			KARG(kTemp, kernelFP, kernelIndFPSubIter, osa_iter);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, ii);
		}
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelFP, global[0], global[1], global[2], local[0], local[1], local[2], 0, CLCommandQueue[0], kTemp.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch forward projection kernel\n", -1);
#else
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelFP, cl::NDRange(), global, local, NULL);
		OCL_CHECK(status, "\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward projection kernel launched successfully\n");
		}
#if defined(CUDA) || defined(HIP)
#ifndef AF
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "\n", -1);
#endif
#else
#ifndef AF
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
#endif
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("Forward projection completed in %f seconds\n", milliseconds);
			// Previously, the events weren't correctly destroyed
			cuEventDestroy(tStart);
			cuEventDestroy(tEnd);
#else
			CLCommandQueue[0].finish();
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("Forward projection completed in %f seconds\n", tDiff);
#endif // END CUDA
		}
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
#if defined(CUDA) || defined(HIP)
	inline int backwardProjection(scalarStruct & inputScalars, Weighting & w_vec, uint32_t osa_iter, uint32_t timestep, std::vector<int64_t>&length, uint64_t m_size, const bool compSens = false, int ii = 0, const int uu = 0,
		const int queueIdx = 0, const bool newInput = true) {
		if (inputScalars.verbose >= 3 || DEBUG)
#else
	inline int backwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const uint32_t timestep, const std::vector<int64_t>&length, const uint64_t m_size, const bool compSens = false, const int32_t ii = 0, const int uu = 0,
		int ee = -1, const int queueIdx = 0, const bool newInput = true) {
		if (inputScalars.verbose >= 3)
#endif // END CUDA
			mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kTemp = BPArgs;
#else
		if (ee < 0)
			ee = uu;
		cl_int status = CL_SUCCESS;
		kernelIndBPSubIter = kernelIndBP;
#endif // END CUDA
		if (inputScalars.listmode > 0 && compSens) {
			kernelApu = kernelBP;
			kernelBP = kernelSensList;
#if defined(CUDA) || defined(HIP)
			kTemp = SensArgs;
#else
			kernelIndBPSubIter = kernelIndSens;
#endif // END CUDA
		}
		int indD = osa_iter + timestep * inputScalars.subsets;
#if defined(CUDA) || defined(HIP)
		CUevent tStart, tEnd;
#else
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}

		if (!inputScalars.CT && (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.BPType == 4)) {
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_atten[osa_iter]);
			}
			else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.size_atten > inputScalars.im_dim[0]) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenB[timestep]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenIm[timestep]);
				}
				else {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenB[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenIm[0]);
				}
			}
		}

		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
#if defined(CUDA) || defined(HIP)
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
				global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
				global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
				global[2] = length[indD];
			}
			else if (inputScalars.listmode > 0 && compSens) {
				global[0] = static_cast<size_t>(inputScalars.det_per_ring + erotusSens[0]) / local[0];
				global[1] = (static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1]) / local[1];
				global[2] = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.nLayers);
			}
#else
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
				global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[indD]) };
			else if (inputScalars.listmode > 0 && compSens)
				if (inputScalars.nLayers > 1)
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1],
					static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.nLayers) };
				else
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) };
#endif // END CUDA
			else {
				erotus[0] = length[indD] % local_size[0];

				if (erotus[0] > 0)
					erotus[0] = (local_size[0] - erotus[0]);
#if defined(CUDA) || defined(HIP)
				global[0] = (length[indD] + erotus[0]) / local[0];
				global[1] = 1;
				global[2] = 1;
#else
				global = { static_cast<cl::size_type>(length[indD] + erotus[0]), 1, 1 };
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
				mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#else
				mexPrintBase("global.dimensions() = %u\n", global.dimensions());
				mexPrintBase("local.dimensions() = %u\n", local.dimensions());
				mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
#endif // END CUDA
				mexPrintBase("m_size = %u\n", m_size);
				mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
				mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
				mexPrintBase("length[indD] = %u\n", length[indD]);
				mexPrintBase("listmode = %u\n", inputScalars.listmode);
				mexPrintBase("rings = %u\n", inputScalars.rings);
				mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("no_norm = %u\n", no_norm);
				mexPrintBase("ii = %u\n", ii);
				mexPrintBase("memSize = %u\n", memSize);
				mexPrintBase("compSens = %u\n", compSens);
				mexPrintBase("osa_iter = %u\n", osa_iter);
				mexPrintBase("timestep = %u\n", timestep);
				mexEval();
			}

			// Set kernelBP arguments
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1)
							subset = osa_iter;
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFPB[subset]);
					}
					else
						if (inputScalars.maskFPZ > 1) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP3[osa_iter]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP);
				}
				if (inputScalars.maskBP) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB);
					}
					else
#if !defined(CUDA) && !defined(HIP)
						if (inputScalars.maskBPZ > 1) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3);
						}
						else
#endif // END CUDA
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP);
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0)
				KARG(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
			if (compSens) {
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.rings);
			}
			else {
				if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
				}
				else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
				if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT)) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
				}
				else if (inputScalars.indexBased || inputScalars.listmode > 0) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[0][0]);
				}
				else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			}
			if (compSens) {
				if (inputScalars.normalization_correction)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_normFull[0]);
				if (inputScalars.scatter)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scatFull[0]);
			}
			else {
				if (inputScalars.normalization_correction)
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[osa_iter]);
				if (inputScalars.scatter)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scat[timestep][osa_iter]);
			}
#if defined(CUDA) || defined(HIP)
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[uu]);
#else
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
#endif // END CUDA
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d[ii]);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, bmax[ii]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.raw)
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_L[osa_iter]);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, no_norm);
#if defined(CUDA) || defined(HIP)
			KARG(kTemp, kernelBP, kernelIndBPSubIter, m_size);
#else
			KARG(kTemp, kernelBP, kernelIndBPSubIter, static_cast<cl_ulong>(m_size));
#endif // END CUDA
			KARG(kTemp, kernelBP, kernelIndBPSubIter, osa_iter);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, ii);
			}
		else {
			if (inputScalars.CT) {

				// Projector type 7 always reads the projection data from a buffer, no image/texture is needed
				if (!inputScalars.useBuffers && inputScalars.BPType != 7) {
#if defined(CUDA) || defined(HIP)
					std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
					arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
					arr3DDesc.NumChannels = 1;
					arr3DDesc.Height = inputScalars.nColsD;
					arr3DDesc.Width = inputScalars.nRowsD;
					arr3DDesc.Depth = length[indD];
#else
					cl::size_type imX = inputScalars.nRowsD;
					cl::size_type imY = inputScalars.nColsD;
					cl::size_type imZ = length[indD];
#endif // END CUDA
					if (inputScalars.BPType == 5) {
#if defined(CUDA) || defined(HIP)
						arr3DDesc.Height++;
						arr3DDesc.Width++;
#else
						imX++;
						imY++;
#endif // END CUDA
					}
#if defined(CUDA) || defined(HIP)
					if (DEBUG) {
						mexPrintBase("arr3DDesc.NumChannels= %u\n", arr3DDesc.NumChannels);
						mexPrintBase("arr3DDesc.Height = %u\n", arr3DDesc.Height);
						mexPrintBase("arr3DDesc.Width = %u\n", arr3DDesc.Width);
						mexPrintBase("arr3DDesc.Depth = %u\n", arr3DDesc.Depth);
						mexEval();
						//mexEvalString("pause(2);");
#else
					cl::detail::size_t_array region = { imX, imY, imZ };
					// Make the inputs persistent such that they are only created when required
					if (newInput) {
						if (BPImageDims[0] != imX || BPImageDims[1] != imY || BPImageDims[2] != imZ) {
							d_inputImage = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
							OCL_CHECK(status, "Image creation failed\n", -1);
							BPImageDims[0] = imX;
							BPImageDims[1] = imY;
							BPImageDims[2] = imZ;
						}
						status = CLCommandQueue[0].enqueueCopyBufferToImage(d_output, d_inputImage, 0, origin, region);
						OCL_CHECK(status, "Image copy failed\n", -1);
					}
#endif // END CUDA
					}
#if defined(CUDA) || defined(HIP)
				// Create only when necessary
				if (newInput) {
					if (BPImageDims[0] != arr3DDesc.Width || BPImageDims[1] != arr3DDesc.Height || BPImageDims[2] != arr3DDesc.Depth) {
						if (BPImageDims[2] != 0) {
							getErrorString(cuTexObjectDestroy(d_inputImage));
							getErrorString(cuArrayDestroy(BPArray));
						}
						status = cuArray3DCreate(&BPArray, &arr3DDesc);
						CUDA_CHECK(status, "Array creation failed\n", -1);
						if (DEBUG)
							mexPrint("Array creation succeeded\n");
						CUDA_RESOURCE_DESC resDescIm;
						std::memset(&resDescIm, 0, sizeof(resDescIm));
						std::memset(&texDesc, 0, sizeof(texDesc));
						std::memset(&viewDesc, 0, sizeof(viewDesc));
						viewDesc.height = inputScalars.nColsD;
						viewDesc.width = inputScalars.nRowsD;
						viewDesc.depth = length[indD];
						viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
						resDescIm.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
						resDescIm.res.array.hArray = BPArray;
						texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
						texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
						if (inputScalars.BPType == 4) {
							texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
						}
						else {
							texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
							viewDesc.height++;
							viewDesc.width++;
						}
						status = cuTexObjectCreate(&d_inputImage, &resDescIm, &texDesc, &viewDesc);
						CUDA_CHECK(status, "Image creation failed\n", -1);
						BPImageDims[0] = arr3DDesc.Width;
						BPImageDims[1] = arr3DDesc.Height;
						BPImageDims[2] = arr3DDesc.Depth;
					}
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
					cpy3d.Depth = length[indD];
					if (inputScalars.BPType == 5) {
						cpy3d.srcPitch += sizeof(float);
						cpy3d.srcHeight++;
						cpy3d.WidthInBytes += sizeof(float);
						cpy3d.Height++;
					}
					status = cuMemcpy3DAsync(&cpy3d, CLCommandQueue[0]);
					CUDA_CHECK(status, "Array mem copy failed\n", -1);
					if (DEBUG)
						mexPrint("Array mem copy succeeded\n");
				}
				}
			if (inputScalars.BPType == 4) {
				global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
				global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
#else
				if (inputScalars.BPType == 4)
#endif // END CUDA
					if (!inputScalars.largeDim)
						if (!inputScalars.useHelical)
#if defined(CUDA) || defined(HIP)
							global[2] = (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS;
#else
							global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS };
#endif // END CUDA
						else
#if defined(CUDA) || defined(HIP)
							global[2] = (inputScalars.Nz[ii] + NVOXELSHELICAL - 1) / NVOXELSHELICAL;
#else
							global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELSHELICAL - 1) / NVOXELSHELICAL };
#endif // END CUDA
					else
#if defined(CUDA) || defined(HIP)
						global[2] = inputScalars.Nz[ii];
			}
#else
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };
#endif // END CUDA
				else if (inputScalars.BPType == 5) {
#if defined(CUDA) || defined(HIP)
					if (inputScalars.pitch) {
						global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / local[0];
						global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / local[1];
						global[2] = inputScalars.Nz[ii];
#else
					if (inputScalars.pitch)
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };
					else
						global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5 };
#endif // END CUDA
					}
#if defined(CUDA) || defined(HIP)
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
#else
				else
					global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };
#endif // END CUDA

				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("local[0] = %u\n", local[0]);
					mexPrintBase("local[1] = %u\n", local[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
					mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
#if defined(CUDA) || defined(HIP)
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#else
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("local.dimensions() = %u\n", local.dimensions());
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
#endif // END CUDA
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[indD] = %u\n", length[indD]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("memSize = %u\n", memSize);
					if (inputScalars.BPType == 4) {
						mexPrintBase("dL = %f\n", inputScalars.dL);
						mexPrintBase("kerroin4[ii] = %f\n", w_vec.kerroin4[ii]);
					}
					mexEval();
				}
				if (inputScalars.offset && inputScalars.BPType != 7)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_T[osa_iter]);
				if (inputScalars.BPType == 5 || inputScalars.BPType == 4 || inputScalars.BPType == 7) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d[ii]);
					if (inputScalars.BPType == 5) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.d_Scale[ii]);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.dSizeBP);
					}
					else if (inputScalars.BPType == 4) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, w_vec.kerroin4[ii]);
					}
				}
				if (inputScalars.BPType == 4) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputImage);
					if (inputScalars.CT && inputScalars.DSC > 0.f) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_angle);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.DSC);
					}
					//mexPrint("2!!!!\n");
					KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
					//mexPrint("3!!!!\n");
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					}
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
#if defined(CUDA) || defined(HIP)
					//mexPrint("4!!!!\n");
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[uu]);
					//mexPrint("5!!!!\n");
#else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
#endif // END CUDA
				}
				else {
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					}
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
					// Only when the kernel was built with -DGEOM5
					if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_geomProj5[osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputImage);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
#if defined(CUDA) || defined(HIP)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[uu]);
#else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
#endif // END CUDA
					if (inputScalars.meanBP) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_meanBP);
					}
				}
				if (inputScalars.normalization_correction)
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[osa_iter]);
			}
			else {
#if defined(CUDA) || defined(HIP)
				if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
					global[0] = (inputScalars.nRowsD + erotus[0]) / local[0];
					global[1] = (inputScalars.nColsD + erotus[1]) / local[1];
					global[2] = length[indD];
				}
				else if (inputScalars.listmode > 0 && compSens) {
					global[0] = static_cast<size_t>(inputScalars.det_per_ring + erotusSens[0]) / local[0];
					global[1] = (static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1]) / local[1];
					global[2] = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings);
				}
#else
				if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
					global = { inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], static_cast<size_t>(length[indD]) };
				else if (inputScalars.listmode > 0 && compSens)
					global = { static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) };
#endif // END CUDA
				else {
					erotus[0] = length[indD] % local_size[0];

					if (erotus[0] > 0)
						erotus[0] = (local_size[0] - erotus[0]);
#if defined(CUDA) || defined(HIP)
					global[0] = (length[indD] + erotus[0]) / local[0];
					global[1] = 1;
					global[2] = 1;
#else
					global = { static_cast<cl::size_type>(length[indD] + erotus[0]), 1, 1 };
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#else
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("local.dimensions() = %u\n", local.dimensions());
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
#endif // END CUDA
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[indD] = %u\n", length[indD]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("osa_iter = %u\n", osa_iter);
					mexPrintBase("memSize = %u\n", memSize);
					mexEval();
				}
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, bmax[ii]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.d_Scale4[ii]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
				if (compSens) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.rings);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.det_per_ring);
				}
				else {
					if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT || (inputScalars.listmode > 0 && !inputScalars.indexBased))) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
					}
					else if (inputScalars.indexBased || inputScalars.listmode > 0) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[0][0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][inputScalars.osa_iter0]);
				}
				if (inputScalars.maskFP || inputScalars.maskBP) {
					if (inputScalars.maskFP) {
						if (inputScalars.useBuffers) {
							int subset = 0;
							if (inputScalars.maskFPZ > 1)
								subset = osa_iter;
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFPB[subset]);
						}
						else
							if (inputScalars.maskFPZ > 1) {
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP3[osa_iter]);
							}
							else
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP);
					}
					if (inputScalars.maskBP) {
						if (inputScalars.useBuffers) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB);
						}
						else
#if !defined(CUDA) && !defined(HIP)
							if (inputScalars.maskBPZ > 1) {
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3);
							}
							else
#endif // END CUDA
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP);
					}
				}
				KARG(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsetsUsed > 1 && inputScalars.listmode == 0) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xyindex[osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zindex[osa_iter]);
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
					if (!inputScalars.loadTOF) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[0][0]);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[0][0]);
					}
					else {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[timestep][osa_iter]);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[timestep][osa_iter]);
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.TOF) {
					if (!inputScalars.loadTOF) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[0][0]);
					}
					else {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[timestep][osa_iter]);
					}
				}
				if (inputScalars.raw) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_L[osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.det_per_ring);
				}
				if (inputScalars.normalization_correction)
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[osa_iter]);
				if (inputScalars.scatter)
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scat[timestep][osa_iter]);
#if defined(CUDA) || defined(HIP)
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[uu]);
#else
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
#endif // END CUDA
				}
			KARG(kTemp, kernelBP, kernelIndBPSubIter, no_norm);
			if (inputScalars.CT && inputScalars.maskBP && (inputScalars.BPType == 4 || inputScalars.BPType == 5 || inputScalars.BPType == 7)) {
				if (inputScalars.useBuffers) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB);
				}
				else
#if !defined(CUDA) && !defined(HIP)
					if (inputScalars.maskBPZ > 1) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3);
					}
					else
#endif // END CUDA
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP);
			}
			if (inputScalars.CT) {
#if defined(CUDA) || defined(HIP)
				KARG(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
#else
				KARG(kTemp, kernelBP, kernelIndBPSubIter, static_cast<cl_long>(length[indD]));
#endif // END CUDA
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kTemp, kernelBP, kernelIndBPSubIter, m_size);
#else
				KARG(kTemp, kernelBP, kernelIndBPSubIter, static_cast<cl_ulong>(m_size));
#endif // END CUDA
				KARG(kTemp, kernelBP, kernelIndBPSubIter, osa_iter);
			}
			KARG(kTemp, kernelBP, kernelIndBPSubIter, ii);
			}
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size()) {
			status = cuEventRecord(evMain, CLCommandQueue[0]);
			CUDA_CHECK(status, "Failed to record main stream event\n", -1);
			status = cuStreamWaitEvent(sideQueues[queueIdx - 1], evMain, 0);
			CUDA_CHECK(status, "Failed to make side stream wait for the main stream\n", -1);
			status = cuLaunchKernel(kernelBP, global[0], global[1], global[2], local[0], local[1], local[2], 0, sideQueues[queueIdx - 1], kTemp.data(), 0);
		}
		else
			status = cuLaunchKernel(kernelBP, global[0], global[1], global[2], local[0], local[1], local[2], 0, CLCommandQueue[0], kTemp.data(), 0);
		CUDA_CHECK(status, "Failed to launch backprojection kernel\n", -1);
#else
		if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size()) {
			cl::Event evMain;
			status = CLCommandQueue[0].enqueueMarkerWithWaitList(nullptr, &evMain);
			OCL_CHECK(status, "Failed to enqueue main queue marker\n", -1);
			std::vector<cl::Event> waitList = { evMain };
			status = sideQueues[queueIdx - 1].enqueueBarrierWithWaitList(&waitList);
			OCL_CHECK(status, "Failed to enqueue side queue barrier\n", -1);
			status = sideQueues[queueIdx - 1].enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, local, NULL);
			OCL_CHECK(status, "\n", -1);
			status = sideQueues[queueIdx - 1].flush();
			OCL_CHECK(status, "Failed to flush side queue\n", -1);
		}
		else {
			status = CLCommandQueue[0].enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, local, NULL);
			OCL_CHECK(status, "\n", -1);
		}
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Backprojection kernel launched successfully\n");
		}
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size() ? sideQueues[queueIdx - 1] : CLCommandQueue[0]);
		
#ifndef AF
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Synchronization failed after backprojection\n", -1);
#endif
#else
#ifndef AF
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
#endif
#endif // END CUDA
		if (inputScalars.listmode > 0 && compSens) {
			kernelBP = kernelApu;
		}
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("Backprojection completed in %f seconds\n", milliseconds);
			cuEventDestroy(tStart);
			cuEventDestroy(tEnd);
#else
			if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size())
				sideQueues[queueIdx - 1].finish();
			else
				CLCommandQueue[0].finish();
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("Backprojection completed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		return 0;
		}

	/// <summary>
	/// Release buffers needed only by the initial computation of the sensitivity image including all measurements (e.g. image-based preconditioners 2-3)
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <returns></returns>
	inline void releaseBuffer(const scalarStruct & inputScalars) {
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
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		size_t mem;
		size_t memF;
		unsigned long long mem_loc;
		status = cuMemGetInfo(&memF, &mem);
		CUDA_CHECK(status, "\n", -1);
		int apu = 0;
		cuDeviceGetAttribute(&apu, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, CUDeviceID[0]);
#else
		cl_int status = CL_SUCCESS;
		int64_t mem;
		cl_ulong mem_loc;
		mem = CLDeviceID[0].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status);
		OCL_CHECK(status, "\n", -1);

		mem_loc = CLDeviceID[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status);
#endif // END CUDA
		CHECK(status, "\n", -1);
		if (DEBUG) {
#if defined(CUDA) || defined(HIP)
			mexPrintBase("mem_loc = %u\n", apu);
			mexPrintBase("memFree = %u\n", memF);
#else
			mexPrintBase("mem_loc = %u\n", mem_loc);
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int computeMRP(const scalarStruct & inputScalars, const uint64_t global_size[]) {
		std::vector<void*> kArgs;
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
	inline int computeMRP(const scalarStruct & inputScalars, const uint64_t gSize[]) {
		cl_uint kernelIndMed = 0U;
		uint64_t erotus[2] = { gSize[0] % localPrior[0], gSize[1] % localPrior[1] };
		cl::NDRange global_size(gSize[0] + (localPrior[0] - erotus[0]), gSize[1] + (localPrior[1] - erotus[1]), gSize[2]);
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint("Starting " BACKEND_STR " median kernel computation");
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		unsigned int gSize[3];
		unsigned int erotus[2];
		erotus[0] = localPrior[0] - (global_size[0] % localPrior[0]);
		erotus[1] = localPrior[1] - (global_size[1] % localPrior[1]);
		gSize[0] = (global_size[0] + erotus[0]) / localPrior[0];
		gSize[1] = (global_size[1] + erotus[1]) / localPrior[1];
		gSize[2] = global_size[2];
		status = cuCtxSynchronize();
#else
		CLCommandQueue[0].finish();
#endif // END CUDA
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
		KARG(kArgs, kernelMed, kernelIndMed, d_inputB);
		KARG(kArgs, kernelMed, kernelIndMed, d_W);
		KARG(kArgs, kernelMed, kernelIndMed, d_N[0]);
		KARG(kArgs, kernelMed, kernelIndMed, d_NOrig);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelMed, kernelIndMed, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelMed, kernelIndMed, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelMed, kernelIndMed, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelMed, kernelIndMed, d_eFOVIndices);
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelMed, gSize[0], gSize[1], gSize[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
#else
		cl_int status = CLCommandQueue[0].enqueueNDRangeKernel(kernelMed, cl::NullRange, global_size, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Median filter kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Median kernel launched successfully\n");
		}
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Synchronize failed after MRP kernel\n", -1);
#else
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after MRP kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA MRP kernel completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL MRP kernel computed in %f seconds\n", tDiff);
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int computeNLM(const scalarStruct & inputScalars, Weighting & w_vec, float beta, const int kk = 0) {
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
		}
#else
	inline int computeNLM(const scalarStruct & inputScalars, Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA NLM gradient computation");
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
#else
			mexPrint("Starting OpenCL NLM gradient computation");
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			tStart = std::chrono::steady_clock::now();
		}
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
		const cl_int3 searchWindow = { static_cast<cl_int>(w_vec.Ndx) , static_cast<cl_int>(w_vec.Ndy) , static_cast<cl_int>(w_vec.Ndz) };
		const cl_int3 patchWindow = { static_cast<cl_int>(w_vec.Nlx) , static_cast<cl_int>(w_vec.Nly) , static_cast<cl_int>(w_vec.Nlz) };
		cl_uint kernelIndNLM = 0ULL;
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
#if defined(CUDA) || defined(HIP)
			globalPrior[2] = Nz;
			NzOrig = d_N[0].z;
			d_N[0].z = Nz;
#else
			globalPrior = { globalPrior[0], globalPrior[1], Nz };
			NzOrig = d_N[0].s2;
			d_N[0].s2 = Nz;
#endif // END CUDA
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Nlz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Nlz + w_vec.Ndz, Nz - w_vec.Nlz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Nlz + w_vec.Ndz, Nz };
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
		//const int3 searchWindow = { static_cast<int>(w_vec.Ndx) , static_cast<int>(w_vec.Ndy) , static_cast<int>(w_vec.Ndz) };
		//const int3 patchWindow = { static_cast<int>(w_vec.Nlx) , static_cast<int>(w_vec.Nly) , static_cast<int>(w_vec.Nlz) };
#endif // END CUDA
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
			mexPrintBase("Nz = %u\n", Nz);
			mexPrintBase("kk = %u\n", kk);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("nOffset.x = %u\n", nOffset.x);
			mexPrintBase("nOffset.y = %u\n", nOffset.y);
			mexPrintBase("d_N[0].z = %u\n", d_N[0].z);
#else
			mexPrintBase("nOffset.x = %u\n", nOffset.s0);
			mexPrintBase("nOffset.y = %u\n", nOffset.s1);
			mexPrintBase("d_N[0].z = %u\n", d_N[0].s2);
#endif // END CUDA
			mexPrintBase("w_vec.RDP_gamma = %f\n", w_vec.RDP_gamma);
			mexPrintBase("useImages = %d\n", inputScalars.useImages);
			mexEval();
		}
		KARG(kArgs, kernelNLM, kernelIndNLM, d_W);
		if (inputScalars.useImages) {
			KARG(kArgs, kernelNLM, kernelIndNLM, d_inputI);
		}
		else {
			KARG(kArgs, kernelNLM, kernelIndNLM, d_inputB);
		}
		KARG(kArgs, kernelNLM, kernelIndNLM, d_gaussianNLM);
		KARG(kArgs, kernelNLM, kernelIndNLM, d_N[0]);
		KARG(kArgs, kernelNLM, kernelIndNLM, d_NOrig);
		KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.h2);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelNLM, kernelIndNLM, apu);
#else
		KARG(kArgs, kernelNLM, kernelIndNLM, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelNLM, kernelIndNLM, beta);
		if (w_vec.NLRD || w_vec.NLLange || w_vec.NLGGMRF)
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.RDP_gamma);
		if (w_vec.NLGGMRF) {
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_p);
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_q);
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_c);
		}
		if (w_vec.NLAdaptive)
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.NLAdaptiveConstant);
		if (w_vec.NLM_anatomical)
			if (inputScalars.useImages) {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_urefIm);
			}
			else {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_uref);
			}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelNLM, kernelIndNLM, d_eFOVIndices);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelNLM, kernelIndNLM, nOffset);
		//Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelNLM, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the NLM kernel\n", status);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after NLM kernel\n", status);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelNLM, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the NLM kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after NLM kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA NLM gradient completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL NLM gradient computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			d_N[0].z = NzOrig;
#else
			d_N[0].s2 = NzOrig;
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int computeRDP(const scalarStruct & inputScalars, float gamma, const Weighting & w_vec, float beta, const int kk = 0, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
#else
	inline int computeRDP(const scalarStruct & inputScalars, const float gamma, const Weighting & w_vec, const float beta, const int kk = 0, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA RDP gradient computation");
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL RDP gradient computation");
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if !defined(CUDA) && !defined(HIP)
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed before RDP kernel\n", -1);
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
#if defined(CUDA) || defined(HIP)
			globalPrior[2] = Nz;
			NzOrig = d_N[0].z;
			d_N[0].z = Nz;
#else
			globalPrior = { globalPrior[0], globalPrior[1], Nz };
			NzOrig = d_N[0].s2;
			d_N[0].s2 = Nz;
#endif // END CUDA
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, NzOrig };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { (Nz - NzOrig) / 2, (Nz + NzOrig) / 2 };
		else if (inputScalars.largeDim)
			nOffset = { Nz - NzOrig, Nz };
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		cl_uint kernelIndRDP = 0ULL;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
#endif // END CUDA
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
		KARG(kArgs, kernelRDP, kernelIndRDP, d_W);
		if (inputScalars.useImages) {
			KARG(kArgs, kernelRDP, kernelIndRDP, d_inputI);
		}
		else {
			KARG(kArgs, kernelRDP, kernelIndRDP, d_inputB);
		}
		KARG(kArgs, kernelRDP, kernelIndRDP, d_N[0]);
		KARG(kArgs, kernelRDP, kernelIndRDP, d_NOrig);
		KARG(kArgs, kernelRDP, kernelIndRDP, gamma);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelRDP, kernelIndRDP, apu);
#else
		KARG(kArgs, kernelRDP, kernelIndRDP, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelRDP, kernelIndRDP, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelRDP, kernelIndRDP, d_eFOVIndices);
		if (RDPLargeNeighbor) {
			KARG(kArgs, kernelRDP, kernelIndRDP, d_weights);
			if (useRDPRef)
				if (inputScalars.useImages) {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_RDPrefI);
				}
				else {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_RDPref);
				}
		}
		if (inputScalars.largeDim)
			KARG(kArgs, kernelRDP, kernelIndRDP, nOffset);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelRDP, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the RDP kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after RDP kernel\n", -1);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelRDP, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the RDP kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after RDP kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA RDP gradient completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL RDP gradient computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			d_N[0].z = NzOrig;
#else
			d_N[0].s2 = NzOrig;
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int computeGGMRF(const scalarStruct & inputScalars, float p, float q, float c, float pqc, const Weighting & w_vec, float beta, const int kk = 0) {
#else
	inline int computeGGMRF(const scalarStruct & inputScalars, const float p, const float q, const float c, const float pqc, const Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA GGMRF gradient computation");
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL GGMRF gradient computation");
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if !defined(CUDA) && !defined(HIP)
		CLCommandQueue[0].finish();
		cl_int status = CL_SUCCESS;
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
#if defined(CUDA) || defined(HIP)
			globalPrior[2] = Nz;
			NzOrig = d_N[0].z;
			d_N[0].z = Nz;
#else
			globalPrior = { globalPrior[0], globalPrior[1], Nz };
			NzOrig = d_N[0].s2;
			d_N[0].s2 = Nz;
#endif // END CUDA
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz };
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
#else
		cl_uint kernelIndGGMRF = 0ULL;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
#endif // END CUDA
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
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_W);
		if (inputScalars.useImages) {
			KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_inputI);
		}
		else {
			KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_inputB);
		}
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_weights);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_N[0]);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, p);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, q);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, c);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, pqc);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPrior);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelGGMRF, kernelIndGGMRF, nOffset);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelGGMRF, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the GGMRF kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after GGMRF kernel\n", -1);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelGGMRF, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the GGMRF kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after GGMRF kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA GGMRF gradient completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL GGMRF gradient computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			d_N[0].z = NzOrig;
#else
			d_N[0].s2 = NzOrig;
#endif // END CUDA
		return 0;
	}


#if defined(CUDA) || defined(HIP)
	inline int ProxHelperQ(float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
		unsigned int kernelIndProxRDP = 0ULL;
#else
	inline int ProxHelperQ(const float alpha, const uint64_t gQ) {
		cl_int status = CL_SUCCESS;
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndProxRDP = 0ULL;
#endif // END CUDA
		KARG(kArgs, kernelProxq, kernelIndProxRDP, d_qX);
		KARG(kArgs, kernelProxq, kernelIndProxRDP, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal RDP helper kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after proximal RDP helper kernel\n", -1);
		return 0;
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal RDP helper kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after Proximal RDP helper kernel\n", -1);
		return status;
#endif // END CUDA
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TV prior
	/// </summary>
	/// <param name="q the input TV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTVHelperQ(float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
#else
	inline int ProxTVHelperQ(const float alpha, const uint64_t gQ) {
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		cl_int status = CL_SUCCESS;
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndCPTV = 0ULL;
#endif // END CUDA
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TV kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after proximal TV kernel\n", -1);
		return 0;
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal TV kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after kernel\n", -1);
		return status;
#endif // END CUDA
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TGV prior
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVHelperQ(const scalarStruct & inputScalars, float alpha, const uint64_t globalQ) {
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		status = cuCtxSynchronize();
#else
	inline int ProxTGVHelperQ(const scalarStruct & inputScalars, const float alpha, const uint64_t globalQ) {
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rX);
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rZ);
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rXY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rXZ);
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rYZ);
		}
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TGV kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after TGV kernel\n", -1);
		return 0;
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal TGV kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after kernel\n", -1);
		return status;
#endif // END CUDA
	}

	/// <summary>
	/// Divergence of the TV prior
	/// </summary>
	/// <param name="im the input array from where the divergence is computed"></param>
	/// <param name="input the backprojection, to which the divergence is added"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTVDiv(const scalarStruct & inputScalars) {
#else
	inline int ProxTVDiv(const scalarStruct & inputScalars, uint32_t timestep = 0) { // TODO: remove default argument
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TV divergence");
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
#else
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
#endif // END CUDA
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
#else
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
#endif // END CUDA
			mexEval();
		}
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed before divergence kernel\n", -1);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_N[0]);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_NPrior);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, vec_opencl.d_rhs_os[0]);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_eFOVIndices);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TV divergence kernel\n", -1);

		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		OCL_CHECK(status, "Failed to launch the Proximal TV divergence kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after divergence kernel\n", -1);
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
#if defined(CUDA) || defined(HIP)
	inline int ProxTVGrad(const scalarStruct & inputScalars, float sigma2, const size_t vSize) {
#else
	inline int ProxTVGrad(const scalarStruct & inputScalars, const float sigma2, const size_t vSize) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TV gradient");
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
#else
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTV = 0ULL;
#endif // END CUDA
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
#else
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
#endif // END CUDA
			mexPrintBase("vSize = %u\n", vSize);
			mexEval();
		}
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_N[0]);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_NPrior);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_inputB);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, sigma2);
		if (vSize > 0) {
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vX);
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vY);
			if (!inputScalars.TGV2D)
				KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vZ);
		}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_eFOVIndices);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVGrad, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVGrad, cl::NullRange, globalPriorEFOV, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Proximal TV gradient kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after gradient kernel\n", -1);
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
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVSymmDeriv(const scalarStruct & inputScalars, float sigma2) {
#else
	inline int ProxTGVSymmDeriv(const scalarStruct & inputScalars, const float sigma2) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint("Starting Proximal TGV symmetric derivative");
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
#else
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTGV = 0ULL;
#endif // END CUDA
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
		unsigned int kernelIndCPTGV = 0ULL;
#else
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
#endif // END CUDA
			mexEval();
		}
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_N[0]);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_NPrior);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vX);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vZ);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rX);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rZ);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXY);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXZ);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rYZ);
		}
		else
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXY);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, sigma2);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPrior);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVSymmDeriv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVSymmDeriv, cl::NullRange, globalPriorEFOV, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Proximal TGV symmetric derivative kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			mexPrint("Proximal TGV symmetric derivative kernel launched successfully\n");
#else
			mexPrint("Proximal TV gradient kernel launched successfully\n");
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after symmetric derivative kernel\n", -1);
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
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVDiv(const scalarStruct & inputScalars, float theta, float tau) {
#else
	inline int ProxTGVDiv(const scalarStruct & inputScalars, const float theta, const float tau) {
#endif // END CUDA
		if (inputScalars.verbose >= 3) {
			mexPrint("Starting Proximal TGV divergence");
		}
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
#else
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndCPTGV = 0ULL;
#endif // END CUDA
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			globalPriorEFOV[2] = inputScalars.Nz[0];
		std::vector<void*> kArgs;
#else
			globalPriorEFOV = { globalPriorEFOV[0], globalPriorEFOV[1], inputScalars.Nz[0] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[0].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[0].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[0].s[2]);
#endif // END CUDA
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_N[0]);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_NPrior);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rZ);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXY);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXZ);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rYZ);
		}
		else
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXY);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vZ);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qY);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qZ);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, theta);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, tau);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPrior);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TGV divergence kernel\n", -1);

		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		OCL_CHECK(status, "Failed to launch the Proximal TGV divergence kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after divergence kernel\n", -1);
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
#if defined(CUDA) || defined(HIP)
	inline int elementWiseComp(const bool mult, const uint64_t size[], bool D2 = false) {
		const unsigned int gSize[3] = { static_cast<unsigned int>(size[0]), static_cast<unsigned int>(size[1]), static_cast<unsigned int>(size[2]) };
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
		CUDA_CHECK(status, "Failed to synchronize before element-wise kernel\n", -1);
#else
	inline int elementWiseComp(const bool mult, const uint64_t size[], const bool D2 = false) {
		cl::NDRange gSize = { static_cast<cl::size_type>(size[0]), static_cast<cl::size_type>(size[1]), static_cast<cl::size_type>(size[2]) };
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndE = 0ULL;
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		if (mult) {
			KARG(kArgs, kernelElementMultiply, kernelIndE, d_vector);
			KARG(kArgs, kernelElementMultiply, kernelIndE, d_input);
#if defined(CUDA) || defined(HIP)
			KARG(kArgs, kernelElementMultiply, kernelIndE, D);
#else
			KARG(kArgs, kernelElementMultiply, kernelIndE, static_cast<cl_uchar>(D2));
#endif // END CUDA
			// Compute the kernel
#if defined(CUDA) || defined(HIP)
			status = cuLaunchKernel(kernelElementMultiply, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
#else
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementMultiply, cl::NullRange, gSize, cl::NullRange);
#endif // END CUDA
		}
		else {
			KARG(kArgs, kernelElementDivision, kernelIndE, d_vector);
			KARG(kArgs, kernelElementDivision, kernelIndE, d_input);
			// Compute the kernel
#if defined(CUDA) || defined(HIP)
			status = cuLaunchKernel(kernelElementDivision, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
#else
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementDivision, cl::NullRange, gSize, cl::NullRange);
#endif // END CUDA
		}
		CHECK(status, "Failed to launch the element-wise kernel\n", -1);
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after element-wise kernel\n", -1);
		return 0;
	}

	/// <summary>
	/// The gradient of hyperbolic prior
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma adjustable weighting parameter"></param>
	/// <param name="beta regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int hyperGradient(const scalarStruct & inputScalars, float sigma, const Weighting & w_vec, float beta, const int kk = 0) {
#else
	inline int hyperGradient(const scalarStruct & inputScalars, const float sigma, const Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA hyperbolic prior gradient computation");
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL hyperbolic prior gradient computation");
		cl_int status = CL_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
#if defined(CUDA) || defined(HIP)
			globalPrior[2] = Nz;
			NzOrig = d_N[0].z;
			d_N[0].z = Nz;
#else
			globalPrior = { globalPrior[0], globalPrior[1], Nz };
			NzOrig = d_N[0].s2;
			d_N[0].s2 = Nz;
#endif // END CUDA
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz };
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
		std::vector<void*> kArgs;
#else
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndHyper = 0ULL;
#endif // END CUDA
		KARG(kArgs, kernelHyper, kernelIndHyper, d_W);
		if (inputScalars.useImages) {
			KARG(kArgs, kernelHyper, kernelIndHyper, d_inputI);
		}
		else {
			KARG(kArgs, kernelHyper, kernelIndHyper, d_inputB);
		}
#if defined(CUDA) || defined(HIP)
		float smooth = inputScalars.epps;
#endif // END CUDA
		KARG(kArgs, kernelHyper, kernelIndHyper, d_N[0]);
		KARG(kArgs, kernelHyper, kernelIndHyper, d_NOrig);
		KARG(kArgs, kernelHyper, kernelIndHyper, sigma);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelHyper, kernelIndHyper, smooth);
#else
		KARG(kArgs, kernelHyper, kernelIndHyper, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelHyper, kernelIndHyper, beta);
		KARG(kArgs, kernelHyper, kernelIndHyper, d_weights);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelHyper, kernelIndHyper, d_eFOVIndices);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelHyper, kernelIndHyper, nOffset);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelHyper, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the hyperbolic prior gradient kernel\n", -1);

		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after hyperbolic prior gradient kernel\n", -1);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelHyper, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the hyperbolic prior gradient kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after hyperbolic prior gradient kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA hyperbolic prior gradient completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL hyperbolic prior gradient computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			d_N[0].z = NzOrig;
#else
			d_N[0].s2 = NzOrig;
#endif // END CUDA
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
#if defined(CUDA) || defined(HIP)
	inline int TVGradient(const scalarStruct & inputScalars, float sigma, float smooth, const Weighting & w_vec, float beta, const int kk = 0, float C = 0.f, const int type = 0) {
#else
	inline int TVGradient(const scalarStruct & inputScalars, const float sigma, const float smooth, const Weighting & w_vec, const float beta, const int kk = 0, const float C = 0.f, const int type = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA TV gradient computation");
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL TV gradient computation");
		cl_int status = CL_SUCCESS;
		if (inputScalars.largeDim)
			globalPrior = { globalPrior[0], globalPrior[1], inputScalars.Nz[0] };
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).finish();
		//cl::detail::size_t_array region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("smooth = %f\n", smooth);
			mexPrintBase("beta = %f\n", beta);
			if (type == 2 || type == 3)
				mexPrintBase("C = %f\n", C);
			mexEval();
		}
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
#if defined(CUDA) || defined(HIP)
			globalPrior[2] = Nz;
			NzOrig = d_N[0].z;
			d_N[0].z = Nz;
#else
			globalPrior = { globalPrior[0], globalPrior[1], Nz };
			NzOrig = d_N[0].s2;
			d_N[0].s2 = Nz;
#endif // END CUDA
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - 1 };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { 1, Nz - 1 };
		else if (inputScalars.largeDim)
			nOffset = { 1, Nz };
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#else
		cl_uint kernelIndTV = 0ULL;
#endif // END CUDA
		KARG(kArgs, kernelTV, kernelIndTV, d_W);
		if (inputScalars.useImages) {
			KARG(kArgs, kernelTV, kernelIndTV, d_inputI);
		}
		else {
			KARG(kArgs, kernelTV, kernelIndTV, d_inputB);
		}
		KARG(kArgs, kernelTV, kernelIndTV, d_N[0]);
		KARG(kArgs, kernelTV, kernelIndTV, d_NOrig);
		KARG(kArgs, kernelTV, kernelIndTV, sigma);
		KARG(kArgs, kernelTV, kernelIndTV, smooth);
		KARG(kArgs, kernelTV, kernelIndTV, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelTV, kernelIndTV, d_maskPriorB);
			}
			else
#if !defined(CUDA) && !defined(HIP)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelTV, kernelIndTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelTV, kernelIndTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelTV, kernelIndTV, d_eFOVIndices);
		if (type == 2 || type == 3)
			KARG(kArgs, kernelTV, kernelIndTV, C);
		if (type > 0)
			KARG(kArgs, kernelTV, kernelIndTV, d_refIm);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelTV, kernelIndTV, nOffset);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelTV, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the TV gradient kernel\n", -1);
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after TV gradient kernel\n", -1);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelTV, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the TV gradient kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after TV gradient kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA TV gradient completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL TV gradient computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		if (inputScalars.largeDim)
#if defined(CUDA) || defined(HIP)
			d_N[0].z = NzOrig;
#else
			d_N[0].s2 = NzOrig;
#endif // END CUDA
		return 0;
	}


#if defined(CUDA) || defined(HIP)
	inline int PoissonUpdate(const scalarStruct & inputScalars, float lambda, float epps, float alpha, const int ii = 0) {
#else
	inline int PoissonUpdate(const scalarStruct & inputScalars, const float lambda, const float epps, const float alpha, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA Poisson update (PKMA/MBSREM/BSREM) computation");
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL Poisson update (PKMA/MBSREM/BSREM) computation");
		cl_int status = CL_SUCCESS;
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
		status = cuCtxSynchronize();
		global[0] = (inputScalars.Nx[ii] + erotusPDHG[0][ii]) / localPrior[0];
		global[1] = (inputScalars.Ny[ii] + erotusPDHG[1][ii]) / localPrior[1];
		global[2] = inputScalars.Nz[ii];
		bool apu = inputScalars.enforcePositivity;
#else
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndPoisson = 0ULL;
		global = { inputScalars.Nx[ii] + erotusPDHG[0][ii], inputScalars.Ny[ii] + erotusPDHG[1][ii], inputScalars.Nz[ii] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
			mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
			mexPrintBase("erotusPDHG[0] = %u\n", erotusPDHG[0][ii]);
			mexPrintBase("erotusPDHG[1] = %u\n", erotusPDHG[1][ii]);
			mexPrintBase("localPrior[0] = %u\n", localPrior[0]);
			mexPrintBase("localPrior[1] = %u\n", localPrior[1]);
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].s[2]);
#endif // END CUDA
			mexPrintBase("lambda = %.8f\n", lambda);
			mexPrintBase("alpha = %f\n", alpha);
			mexEval();
		}
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_im);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_rhs);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_N[ii]);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, lambda);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, epps);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, alpha);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelPoisson, kernelIndPoisson, apu);
#else
		KARG(kArgs, kernelPoisson, kernelIndPoisson, static_cast<cl_uchar>(inputScalars.enforcePositivity));
#endif // END CUDA
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelPoisson, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the Poisson update kernel\n", -1);

		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPoisson, cl::NullRange, global, localPrior);
		OCL_CHECK(status, "Failed to launch the Poisson update kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after Poisson update kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA Poisson update completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL Poisson update computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		return 0;
	}

#if defined(CUDA) || defined(HIP)
	inline int PDHGUpdate(const scalarStruct & inputScalars, float epps, float theta, float tau, const int ii = 0) {
#else
	inline int PDHGUpdate(const scalarStruct & inputScalars, const float epps, const float theta, const float tau, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA PDHG update computation");
		CUresult status = CUDA_SUCCESS;
		CUevent tStart, tEnd;
#else
			mexPrint("Starting OpenCL PDHG update computation");
		cl_int status = CL_SUCCESS;
		std::chrono::steady_clock::time_point tStart;
		std::chrono::steady_clock::time_point tEnd;
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventCreate(&tStart, CU_EVENT_DEFAULT);
			cuEventCreate(&tEnd, CU_EVENT_DEFAULT);
#else
			tStart = std::chrono::steady_clock::now();
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
		global[0] = (inputScalars.Nx[ii] + erotusPDHG[0][ii]) / localPrior[0];
		global[1] = (inputScalars.Ny[ii] + erotusPDHG[1][ii]) / localPrior[1];
		global[2] = inputScalars.Nz[ii];
		bool apu = inputScalars.enforcePositivity;
#else
		status = (CLCommandQueue[0]).finish();
		cl_uint kernelIndPDHG = 0ULL;
		global = { inputScalars.Nx[ii] + erotusPDHG[0][ii], inputScalars.Ny[ii] + erotusPDHG[1][ii], inputScalars.Nz[ii] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].s[2]);
#endif // END CUDA
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_im);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_rhs);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_U);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_N[ii]);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, epps);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, theta);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, tau);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelPDHG, kernelIndPDHG, apu);
#else
		KARG(kArgs, kernelPDHG, kernelIndPDHG, static_cast<cl_uchar>(inputScalars.enforcePositivity));
#endif // END CUDA
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tStart, CLCommandQueue[0]);
		status = cuLaunchKernel(kernelPDHG, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, CLCommandQueue[0]);
		CUDA_CHECK(status, "Failed to launch the PDHG update kernel\n", -1);

		status = cuCtxSynchronize();
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPDHG, cl::NullRange, global, localPrior);
		OCL_CHECK(status, "Failed to launch the PDHG update kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
#endif // END CUDA
		CHECK(status, "Queue finish failed after PDHG update kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			cuEventSynchronize(tEnd);
			float milliseconds = 0;
			cuEventElapsedTime(&milliseconds, tStart, tEnd);
			milliseconds /= 1000.f;
			mexPrintBase("CUDA PDHG update completed in %f seconds\n", milliseconds);
#else
			tEnd = std::chrono::steady_clock::now();
			const std::chrono::duration<double> tDiff = tEnd - tStart;
			mexPrintBase("OpenCL PDHG update computed in %f seconds\n", tDiff);
#endif // END CUDA
		}
		return 0;
	}

#if defined(CUDA) || defined(HIP)
	inline int rotateCustom(const scalarStruct & inputScalars, float cosa, float sina, const int ii = 0) {
#else
	inline int rotateCustom(const scalarStruct & inputScalars, const float cosa, const float sina, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3)
#if defined(CUDA) || defined(HIP)
			mexPrint("Starting CUDA bilinear image rotation computation");
		CUresult status = CUDA_SUCCESS;
		std::vector<void*> kArgs;
		global[0] = (inputScalars.Nx[ii] + erotusPrior[0]) / localPrior[0];
		global[1] = (inputScalars.Ny[ii] + erotusPrior[1]) / localPrior[1];
		global[2] = inputScalars.Nz[ii];
#else
			mexPrint("Starting OpenCL bilinear image rotation computation");
		cl_int status = CL_SUCCESS;
		cl_uint kernelIndRot = 0ULL;
		global = { inputScalars.Nx[ii] + erotusPrior[0], inputScalars.Ny[ii] + erotusPrior[1], inputScalars.Nz[ii] };
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
#if defined(CUDA) || defined(HIP)
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].x);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].y);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].z);
#else
			mexPrintBase("d_N.s[0] = %u\n", d_N[ii].s[0]);
			mexPrintBase("d_N.s[1] = %u\n", d_N[ii].s[1]);
			mexPrintBase("d_N.s[2] = %u\n", d_N[ii].s[2]);
#endif // END CUDA
			mexEval();
		}

		KARG(kArgs, kernelRotate, kernelIndRot, d_rhs);
		if (!inputScalars.useBuffers) {
			KARG(kArgs, kernelRotate, kernelIndRot, d_inputI);
		}
		else {
			KARG(kArgs, kernelRotate, kernelIndRot, d_im);
		}
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].x);
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].y);
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].z);
#else
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].s[0]);
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].s[1]);
		KARG(kArgs, kernelRotate, kernelIndRot, d_N[ii].s[2]);
#endif // END CUDA
		KARG(kArgs, kernelRotate, kernelIndRot, cosa);
		KARG(kArgs, kernelRotate, kernelIndRot, sina);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelRotate, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the bilinear image rotation kernel\n", -1);
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Queue finish failed after bilinear image rotation kernel\n", -1);
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
#else
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelRotate, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the bilinear image rotation kernel\n", -1);
		status = (CLCommandQueue[0]).finish();
		OCL_CHECK(status, "Queue finish failed after bilinear image rotation kernel\n", -1);
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint(BACKEND_STR " bilinear image rotation computed");
		return 0;
		}
#if defined(CUDA) || defined(HIP)
	inline int transferTex(const scalarStruct & inputScalars, CUdeviceptr * input, const bool RDP = false, const uint32_t Nz = 1) {

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
		arr3DDesc.Depth = Nz;
		status = cuArray3DCreate(&imArray, &arr3DDesc);
		CUDA_CHECK(status, "Failed to create image array\n", -1);
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
		cpy3d.Depth = Nz;
		status = cuMemcpy3D(&cpy3d);
		CUDA_CHECK(status, "Failed to copy image array\n", -1);
		resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
		resDesc.res.array.hArray = imArray;
		texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
		viewDesc.height = inputScalars.Nx[0];
		viewDesc.width = inputScalars.Ny[0];
		viewDesc.depth = Nz;
		viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
		if (RDP)
			status = cuTexObjectCreate(&d_RDPrefI, &resDesc, &texDesc, &viewDesc);
		else
			status = cuTexObjectCreate(&d_inputI, &resDesc, &texDesc, &viewDesc);
		CUDA_CHECK(status, "Image copy failed\n", -1);
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Synchronization failed\n", -1);
		if (DEBUG)
			mexPrint("Synchronization completed\n");
		return 0;
	}

#endif // END CUDA
	};


