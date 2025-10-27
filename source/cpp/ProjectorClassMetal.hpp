#pragma once
#include "structs.h"
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include "Metal.hpp"
#include "kernelParams.hpp"
#include "PreprocessorUtils.cpp"
#include <simd/simd.h>

class ProjectorClass {
	// Local size
	size_t local_size[3];
	size_t local_sizePrior[3];

	// Crystal pitch
	simd::float2 dPitch;

	// Image dimensions
	simd::int3 d_NOrig, d_NPrior;

	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusPrior[3];
	size_t erotusPriorEFOV[3];
	size_t erotusSens[3];

	// Local and global sizes
	simd::int3 local, global, localPrior, globalPrior, globalPriorEFOV; // or use MTL::Size

public:
	NS::SharedPtr<MTL::Device> mtlDevice;
    NS::SharedPtr<MTL::CommandQueue> mtlCommandQueue;
	NS::SharedPtr<MTL::CommandBuffer> commandBufferFP, commandBufferBP;
	NS::SharedPtr<MTL::ComputeCommandEncoder> kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, kernelElementMultiply, kernelElementDivision, 
		kernelTV, kernelProxTGVSymmDeriv, kernelProxTGVDiv, kernelProxTGVq, kernelPoisson, kernelPDHG, kernelProxRDP, kernelProxq, kernelProxTrans, kernelProxNLM, kernelGGMRF,
		kernelsumma, kernelEstimate, kernelPSF, kernelPSFf, kernelDiv, kernelMult, kernelForward, kernelSensList, kernelApu, kernelHyper, kernelRotate;
	NS::SharedPtr<MTL::Buffer> d_xcenter, d_ycenter, d_zcenter, d_V, d_TOFCenter, d_output, d_meanBP, d_meanFP, d_eFOVIndices, d_weights, d_inputB, d_W, d_gaussianNLM;
	NS::SharedPtr<MTL::Buffer> d_angle;
	std::chrono::steady_clock::time_point tStartLocal, tStartGlobal, tStartAll;
	std::chrono::steady_clock::time_point tEndLocal, tEndGlobal, tEndAll;
	METAL_im_vectors vec_opencl;

	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
	std::vector<simd::float3> b, d, bmax;
	std::vector<simd::uint3> d_N;

	std::vector<NS::SharedPtr<MTL::Buffer>> d_maskFP;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_LFull, d_zindexFull, d_xyindexFull, d_normFull, d_scatFull, d_xFull, d_zFull;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_L;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_Summ;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_meas;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_rand;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_imTemp;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_imFinal;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_zindex;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_xyindex;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_trIndex;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_axIndex;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_TOFIndex;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_norm;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_scat;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_x;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_z;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_atten;
	std::vector<NS::SharedPtr<MTL::Buffer>> d_T;
	NS::SharedPtr<MTL::Buffer> d_rayShiftsSource, d_rayShiftsDetector, d_attenB, d_maskBP;

	// Image origin
	MTL::Origin origin = MTL::Origin(0, 0, 0);
	MTL::Size region = MTL::Size(0, 0, 0);

	std::vector<std::vector<size_t>> erotusBP, erotusPDHG;
	u_char no_norm = 0;
	int proj6 = 1;

	// Here init device and create programs (=libraries in Metal)
	int createProgram(
		NS::SharedPtr<MTL::Library>& libFP,
		NS::SharedPtr<MTL::Library>& libBP,
		NS::SharedPtr<MTL::Library>& libAux,
		NS::SharedPtr<MTL::Library>& libSens,
		const char* header_directory,
		scalarStruct& inputScalars,
		const RecMethods MethodList,
		const Weighting& w_vec,
		const size_t local_size[],
		const int type = -1
	);

	// Create kernels (=ComputeCommandEncoder in Metal)
	int createKernels(
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelFP,
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelBP,
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelNLM,
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelMed,
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelRDP,
		NS::SharedPtr<MTL::ComputeCommandEncoder>& kernelGGMRF,
		const NS::SharedPtr<MTL::Library>& libFP,
		const NS::SharedPtr<MTL::Library>& libBP,
		const NS::SharedPtr<MTL::Library>& libAux,
		const NS::SharedPtr<MTL::Library>& libSens,
		const RecMethods& MethodList,
		const Weighting& w_vec,
		const scalarStruct& inputScalars,
		const int type = -1
	);

    int addProjector(
        scalarStruct& inputScalars,
        Weighting& w_vec,
        const RecMethods& MethodList,
        const char* header_directory,
        const int type = -1
    );

    int createBuffers(
		scalarStruct& inputScalars,
		Weighting& w_vec,
		const float* x,
		const float* z_det,
		const uint32_t* xy_index,
		const uint16_t* z_index,
		const uint16_t* L,
		const int64_t* pituus,
		const float* atten,
		const float* norm,
		const float* extraCorr, 
		const std::vector<int64_t>& length,
		const RecMethods& MethodList,
		const int type = 0
	);

	// Inputs constant values to the kernels, i.e. values that do not change in each time step or iteration
    int initializeKernel(
		scalarStruct& inputScalars, // various scalar parameters defining the build parameters and what features to use
		Weighting& w_vec // specifies some of the special options/parameters used
	);

	// Sets kernel parameters that do not change per iteration but only per time step
	int setDynamicKernelData(
		scalarStruct& inputScalars, // various scalar parameters defining the build parameters and what features to use
		Weighting& w_vec // specifies some of the special options/parameters used
	);

	// Compute the forward projection for the selected projector type
    int forwardProjection(
		const scalarStruct& inputScalars, // various scalar parameters defining the build parameters and what features to use
		Weighting& w_vec, // specifies some of the special options/parameters used
		const uint32_t osa_iter, // current subset (sub-iteration)
		const std::vector<int64_t>& length, // the number of measurements/projection/sinograms per subset
		const uint64_t m_size, // for projector types 1-3, the total number of LORs
		const int32_t ii = 0,
		const int uu = 0
	);

	// Compute the backprojection for the selected projector type
    int backwardProjection(
		const scalarStruct& inputScalars, // various scalar parameters defining the build parameters and what features to use
		Weighting& w_vec, // specifies some of the special options/parameters used
		const uint32_t osa_iter, // current subset (sub-iteration)
		const std::vector<int64_t>& length, // the number of measurements/projection/sinograms per subset
		const uint64_t m_size, // for projector types 1-3, the total number of LORs
		const bool compSens = false, // if true, computes the sensitivity image as well
		const int32_t ii = 0,
		const int uu = 0,
		int ee = -1
	);
};